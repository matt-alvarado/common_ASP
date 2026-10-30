!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! OrgPhaseReactions.h                                                       !!
!!                                                                           !!
!! First-order particle-phase organic reactions (oligomerisation, nitrate    !!
!! hydrolysis, photolysis of particle-phase organics ...).                   !!
!!                                                                           !!
!! Input deck: OrgPhaseReactions.in (optional; with no file, or no reaction  !!
!! lines, nothing happens).  One reaction per line, ';' separated:           !!
!!                                                                           !!
!!   Reactant ; Rate ; Products                                              !!
!!                                                                           !!
!!   Reactant  an organic-phase species (in OrgPhaseChems.in).               !!
!!   Rate      either a number [1/s], or                                     !!
!!             PHOTO <n> <scale>   rate = scale * (current rate [1/s] of    !!
!!             reaction n of GasChemicalMechanism.in, counted in the order  !!
!!             the reaction lines appear), e.g. a photolysis frequency.     !!
!!   Products  'coef1 Name1 + coef2 Name2 ...'.  Coefficients are moles per !!
!!             mole of reactant.  A product that is an organic-phase        !!
!!             species is added to the particle; any other product must be  !!
!!             a gas-phase species and is released to the gas phase (HNO3,  !!
!!             ORA1 ...).  Empty product list: the reactant is only lost.   !!
!!                                                                           !!
!!   Comment lines start with '!'.                                           !!
!!                                                                           !!
!! Each particle is advanced with the exact solution of  dA/dt = -k A  over  !!
!! the step, A(t+dt) = A(t) exp(-k dt); the lost moles reappear as products. !!
!! Species that are present in both the hydrophobic and the hydrophilic     !!
!! organic lists (OrgPhaseChems.in and HydrophilicOrgChems.in) react in      !!
!! both pools, the products going to the pool of the reactant when they    !!
!! exist in it.                                                              !!
!!                                                                           !!
!! Called by StepCondensationAll() after the organic condensation steps.    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE StepOrgPhaseReactions (TimestepSize)

	USE InfrastructuralCode, ONLY : ERROR, Transcript, GetFileHandle, Int2Str
	USE ModelParameters,     ONLY : InputDeckSubDir, Avogadro
	USE GridPointFields,     ONLY : GetTemp, AddToGridCellChemBurden, GasPhaseChemicalRates
	USE Aerosols,            ONLY : Particles, Particle
	USE Chemistry,           ONLY : FindChem, HowManyEvolveGasChems

	IMPLICIT NONE

	REAL*8, INTENT(IN) :: TimestepSize		!! Seconds

	INTEGER, PARAMETER :: MaxRxn = 200, MaxProd = 8

	!! Reaction table, read once
	LOGICAL, SAVE :: Initialized = .FALSE.
	INTEGER, SAVE :: NumRxn = 0
	INTEGER, SAVE :: Reactant(MaxRxn), RateKind(MaxRxn), PhotoRxn(MaxRxn), NumProd(MaxRxn)
	INTEGER, SAVE :: ReactantAq(MaxRxn)
	REAL*8,  SAVE :: RateValue(MaxRxn)		!! 1/s, or the scale for PHOTO
	REAL*8,  SAVE :: ProdCoef(MaxRxn, MaxProd)
	INTEGER, SAVE :: ProdOrg(MaxRxn, MaxProd), ProdAq(MaxRxn, MaxProd), ProdGas(MaxRxn, MaxProd)

	INTEGER :: FH, Status, I, J, K, P, Pos, NTok
	CHARACTER(len=1024) :: Line, Fld(3), Rest, Item
	CHARACTER(len=64)   :: Nm
	LOGICAL :: FileExists
	REAL*8  :: Rate, Frac, Lost, Temperature
	TYPE(Particle), POINTER :: Cur

	IF (.NOT. Initialized) THEN
		Initialized = .TRUE.
		NumRxn = 0
		INQUIRE (FILE=TRIM(InputDeckSubDir)//'OrgPhaseReactions.in', EXIST=FileExists)
		IF (FileExists) THEN
			FH = GetFileHandle()
			OPEN (UNIT=FH, FILE=TRIM(InputDeckSubDir)//'OrgPhaseReactions.in', STATUS='OLD')
			DO
				READ (FH, "(a)", IOSTAT=Status) Line
				IF (Status .NE. 0) EXIT
				Line = ADJUSTL(Line)
				IF (LEN_TRIM(Line) .EQ. 0 .OR. Line(1:1) .EQ. '!') CYCLE
				Pos = INDEX(Line, '!')
				IF (Pos .GT. 0) Line = Line(1:Pos-1)
				!! split into three fields at ';'
				Rest = Line
				DO K = 1, 3
					Pos = INDEX(Rest, ';')
					IF (Pos .EQ. 0) THEN
						Fld(K) = ADJUSTL(Rest)
						Rest = ''
					ELSE
						Fld(K) = ADJUSTL(Rest(1:Pos-1))
						Rest = Rest(Pos+1:)
					END IF
				END DO
				IF (LEN_TRIM(Fld(1)) .EQ. 0) CYCLE
				NumRxn = NumRxn + 1
				IF (NumRxn .GT. MaxRxn) CALL ERROR("Too many reactions in OrgPhaseReactions.in (max 200)")
				Reactant(NumRxn)    = FindChem(TRIM(Fld(1)), 2)
				ReactantAq(NumRxn)  = FindChem(TRIM(Fld(1)), 3, .TRUE.)
				!! rate
				IF (Fld(2)(1:5) .EQ. 'PHOTO') THEN
					RateKind(NumRxn) = 2
					READ (Fld(2)(6:), *, IOSTAT=Status) PhotoRxn(NumRxn), RateValue(NumRxn)
					IF (Status .NE. 0) CALL ERROR("Bad PHOTO rate in OrgPhaseReactions.in: "//TRIM(Line))
				ELSE
					RateKind(NumRxn) = 1
					READ (Fld(2), *, IOSTAT=Status) RateValue(NumRxn)
					IF (Status .NE. 0) CALL ERROR("Bad rate in OrgPhaseReactions.in: "//TRIM(Line))
				END IF
				!! products
				NumProd(NumRxn) = 0
				Rest = ADJUSTL(Fld(3))
				DO WHILE (LEN_TRIM(Rest) .GT. 0)
					Pos = INDEX(Rest, '+')
					IF (Pos .EQ. 0) THEN
						Item = ADJUSTL(Rest)
						Rest = ''
					ELSE
						Item = ADJUSTL(Rest(1:Pos-1))
						Rest = ADJUSTL(Rest(Pos+1:))
					END IF
					IF (LEN_TRIM(Item) .EQ. 0) CYCLE
					NumProd(NumRxn) = NumProd(NumRxn) + 1
					IF (NumProd(NumRxn) .GT. MaxProd) CALL ERROR("Too many products in OrgPhaseReactions.in: "//TRIM(Line))
					P = NumProd(NumRxn)
					READ (Item, *, IOSTAT=Status) ProdCoef(NumRxn, P), Nm
					IF (Status .NE. 0) CALL ERROR("Bad product in OrgPhaseReactions.in: "//TRIM(Line))
					ProdOrg(NumRxn, P) = FindChem(TRIM(Nm), 2, .TRUE.)
					ProdAq(NumRxn, P)  = FindChem(TRIM(Nm), 3, .TRUE.)
					ProdGas(NumRxn, P) = FindChem(TRIM(Nm), 0, .TRUE.)
					IF (ProdOrg(NumRxn, P) .LE. 0 .AND. ProdGas(NumRxn, P) .LE. 0) &
						CALL ERROR("Product "//TRIM(Nm)//" of OrgPhaseReactions.in is neither an organic-phase nor a gas-phase species")
				END DO
			END DO
			CLOSE (FH)
			CALL Transcript("OrgPhaseReactions.in: "//TRIM(Int2Str(NumRxn))//" particle-phase organic reactions read")
		END IF
	END IF

	IF (NumRxn .EQ. 0) RETURN

	Temperature = GetTemp()

	DO I = 1, NumRxn
		IF (RateKind(I) .EQ. 2) THEN
			Rate = RateValue(I) * GasPhaseChemicalRates(PhotoRxn(I))
		ELSE
			Rate = RateValue(I)
		END IF
		IF (Rate .LE. 0.0D0) CYCLE
		Frac = 1.0D0 - EXP(-Rate * TimestepSize)

		Cur => Particles%First
		DO WHILE (ASSOCIATED(Cur))
			IF (Cur%NumberOfParticles .GT. 1.0D-40 .AND. .NOT. Cur%Dry) THEN
				!! hydrophobic pool
				IF (Reactant(I) .GT. 0) THEN
					Lost = Cur%OrgChems(Reactant(I)) * Frac
					IF (Lost .GT. 0.0D0) THEN
						Cur%OrgChems(Reactant(I)) = Cur%OrgChems(Reactant(I)) - Lost
						DO P = 1, NumProd(I)
							IF (ProdOrg(I, P) .GT. 0) THEN
								Cur%OrgChems(ProdOrg(I, P)) = Cur%OrgChems(ProdOrg(I, P)) + ProdCoef(I, P) * Lost
							ELSE IF (ProdGas(I, P) .LE. HowManyEvolveGasChems) THEN
								CALL AddToGridCellChemBurden (ProdGas(I, P), ProdCoef(I, P) * Lost * Avogadro * Cur%NumberOfParticles)
							END IF
						END DO
					END IF
				END IF
				!! hydrophilic pool
				IF (ReactantAq(I) .GT. 0) THEN
					Lost = Cur%AqOrgChems(ReactantAq(I)) * Frac
					IF (Lost .GT. 0.0D0) THEN
						Cur%AqOrgChems(ReactantAq(I)) = Cur%AqOrgChems(ReactantAq(I)) - Lost
						DO P = 1, NumProd(I)
							IF (ProdAq(I, P) .GT. 0) THEN
								Cur%AqOrgChems(ProdAq(I, P)) = Cur%AqOrgChems(ProdAq(I, P)) + ProdCoef(I, P) * Lost
							ELSE IF (ProdOrg(I, P) .GT. 0) THEN
								Cur%OrgChems(ProdOrg(I, P)) = Cur%OrgChems(ProdOrg(I, P)) + ProdCoef(I, P) * Lost
							ELSE IF (ProdGas(I, P) .LE. HowManyEvolveGasChems) THEN
								CALL AddToGridCellChemBurden (ProdGas(I, P), ProdCoef(I, P) * Lost * Avogadro * Cur%NumberOfParticles)
							END IF
						END DO
					END IF
				END IF
			END IF
			Cur => Cur%Next
		END DO
	END DO

	!! activity coefficients follow the new composition
	Cur => Particles%First
	DO WHILE (ASSOCIATED(Cur))
		IF (Cur%NumberOfParticles .GT. 1.0D-40 .AND. .NOT. Cur%Dry) CALL UpdateHydrophobicUNIFAC (Cur, Temperature)
		Cur => Cur%Next
	END DO

END SUBROUTINE StepOrgPhaseReactions
