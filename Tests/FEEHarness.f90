!! ASP (c), 2004-2026, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! FEEHarness.f90
!!
!! Test harness for the electrolyte equilibrium solver (FindElectrolyteEquilibrium).
!! It is not part of the model: it is compiled only by Tests/Makefile, for two test
!! programs (see Tests/README.md):
!!
!!   asp_feedump  the ASP model (DevelopmentDriver.f90) with FEEStateHook pointed at
!!                FEEDumpHook, which writes sampled solver states:
!!                  ASP_FEE_DUMP=<file>      write the states to <file> (off if unset)
!!                  ASP_FEE_DUMP_EVERY=<n>   every n-th call (default 100)
!!                  ASP_FEE_DUMP_ALL=1       also calls that update the activity coefficients
!!                Only calls after the first equilibration are sampled.
!!
!!   FEEReplay    re-solves every dumped state with the solver the deck selects, times
!!                the solver calls and compares each result with the dumped one.
!!
!! One record holds everything FindElectrolyteEquilibriumCore reads (record kind 1, with
!! grid temperature, pressure and water vapour) or writes (kind 2).

MODULE FEEHarness

	USE Aerosols, ONLY : Particle, Particles, FEEStateHook, FindElectrolyteEquilibriumCore,	&
	                     EquilibrateAqueousDissociationReaction, EquilibriumConstantsRatio,	&
	                     DeliquesenceRH, EqCoeff, EqViolation, InitSolidCircuits, NCirc, CircN, CircRxn,	&
	                     CircCoef, CircuitApplies, StableSideAvailable, CompletesUnstableSet,	&
	                     KusikMeissner,							&
	                     FEECountOuter, FEECountEADR, FEECountECREval, FEECountNewtonFail,	&
	                     FEECountNewtonSing, FEECountNewtonMaxIt, FEECountNewtonDrop,		&
	                     FEECountNewtonCycles, FEECountPairRefuse, FEECountResolved,		&
	                     FEECountDry, FEECountPreAssoc

	IMPLICIT NONE

	PRIVATE
	PUBLIC :: FEEHarnessAttach, ReplayFEE

	!Dump settings (read from the environment at the first solver call)
	LOGICAL, SAVE   :: DumpSettingsRead = .FALSE., DumpOn = .FALSE., DumpAll = .FALSE., DumpThisCall = .FALSE.
	INTEGER, SAVE   :: DumpEvery = 100, DumpUnit = -1
	INTEGER*8, SAVE :: DumpCount = 0

	!ViolationDiag tallies, per reaction
	INTEGER, SAVE :: VDViol(200) = 0, VDMoved(200) = 0, VDFixed(200) = 0, VDStuck(200) = 0
	REAL*8, SAVE  :: VDSize(200) = 0.

CONTAINS

!Point the model's test hook at the dump routine.  Called from C (FEEDumpAttach.c) before the
!ASP main program starts, so it only sets the pointer; the settings are read at the first call.
SUBROUTINE FEEHarnessAttach () BIND(C, NAME="fee_harness_attach")
	IMPLICIT NONE
	FEEStateHook => FEEDumpHook
END SUBROUTINE FEEHarnessAttach

SUBROUTINE FEEDumpHook (P, UpdateThermo, FirstEquilibration, ReturnType, Before)
	IMPLICIT NONE
	TYPE(PARTICLE),POINTER :: P
	LOGICAL :: UpdateThermo, FirstEquilibration, Before
	INTEGER :: ReturnType
	CHARACTER(LEN=512) :: V
	INTEGER :: L, ST

	IF (.NOT. DumpSettingsRead) THEN
		DumpSettingsRead = .TRUE.
		CALL GET_ENVIRONMENT_VARIABLE ("ASP_FEE_DUMP_ALL", V, L, ST)
		DumpAll = (ST .EQ. 0 .AND. L .GT. 0)
		CALL GET_ENVIRONMENT_VARIABLE ("ASP_FEE_DUMP_EVERY", V, L, ST)
		IF (ST .EQ. 0 .AND. L .GT. 0) READ (V, *) DumpEvery
		CALL GET_ENVIRONMENT_VARIABLE ("ASP_FEE_DUMP", V, L, ST)
		IF (ST .EQ. 0 .AND. L .GT. 0) THEN
			DumpOn = .TRUE.
			OPEN (NEWUNIT=DumpUnit, FILE=TRIM(V), FORM="UNFORMATTED", ACCESS="STREAM", STATUS="REPLACE")
		END IF
	END IF
	IF (.NOT. DumpOn) RETURN

	IF (Before) THEN
		DumpThisCall = .FALSE.
		IF ((.NOT.UpdateThermo .OR. DumpAll) .AND. .NOT.FirstEquilibration) THEN
			DumpCount = DumpCount + 1
			DumpThisCall = (MOD(DumpCount, INT(DumpEvery,8)) .EQ. 0)
		END IF
	END IF
	IF (DumpThisCall) CALL WriteFEEState (P, UpdateThermo, FirstEquilibration, ReturnType, Before)
END SUBROUTINE FEEDumpHook

SUBROUTINE WriteFEEState (P, UpdateThermo, FirstEquilibration, ReturnType, Pre)
	USE GridPointFields, ONLY : GetTemp, GetPress, GetGridCellChemBurden
	IMPLICIT NONE
	TYPE(PARTICLE),POINTER :: P
	LOGICAL :: UpdateThermo, FirstEquilibration, Pre
	INTEGER :: ReturnType, IU, IF1
	IU = DumpUnit
	IF (Pre) THEN
		IF1 = 0 ; IF (UpdateThermo) IF1 = IF1 + 1 ; IF (FirstEquilibration) IF1 = IF1 + 2
		WRITE (IU) 1, IF1, ReturnType, GetTemp(), GetPress(), GetGridCellChemBurden(1)
	ELSE
		WRITE (IU) 2, 0, ReturnType
	END IF
	WRITE (IU) P%NumberOfParticles, P%Dry, P%Temperature, P%EmbryoRadius, P%InsolubleRadius, &
	           P%EffectiveRadius, P%IonicStr, P%WaterActivity, P%InorgSolnDensity, P%SolutionDensity, &
	           P%InsolubleDensity, P%ParticleDensity, P%SurfaceTension
	WRITE (IU) SIZE(P%AqChems), SIZE(P%GammaMixed), SIZE(P%OrgChems), SIZE(P%AqOrgChems)
	WRITE (IU) P%AqChems, P%GammaMixed, P%OrgChems, P%AqOrgChems
END SUBROUTINE WriteFEEState

SUBROUTINE ReadFEEState (IU, P, Kind, Flags, RT, T, Pr, H2O, IOS)
	IMPLICIT NONE
	INTEGER :: IU, Kind, Flags, RT, IOS, NA, NG, NO, NQ
	TYPE(PARTICLE),POINTER :: P
	REAL*8 :: T, Pr, H2O
	READ (IU, IOSTAT=IOS) Kind, Flags, RT
	IF (IOS .NE. 0) RETURN
	IF (Kind .EQ. 1) READ (IU, IOSTAT=IOS) T, Pr, H2O
	READ (IU, IOSTAT=IOS) P%NumberOfParticles, P%Dry, P%Temperature, P%EmbryoRadius, P%InsolubleRadius, &
	           P%EffectiveRadius, P%IonicStr, P%WaterActivity, P%InorgSolnDensity, P%SolutionDensity, &
	           P%InsolubleDensity, P%ParticleDensity, P%SurfaceTension
	READ (IU, IOSTAT=IOS) NA, NG, NO, NQ
	IF (NA .NE. SIZE(P%AqChems) .OR. NG .NE. SIZE(P%GammaMixed) .OR. NO .NE. SIZE(P%OrgChems) &
	    .OR. NQ .NE. SIZE(P%AqOrgChems)) THEN
		IOS = -99 ; RETURN
	END IF
	READ (IU, IOSTAT=IOS) P%AqChems, P%GammaMixed, P%OrgChems, P%AqOrgChems
END SUBROUTINE ReadFEEState

!Violation V of reaction R's equilibrium condition, as in EqViolation (Thermodynamics.h);
!.FALSE. if EqViolation skips R.
LOGICAL FUNCTION ReactionViolation (R, InParticle, RH, V)
	USE Chemistry, ONLY : AqEquilibriaList
	USE ModelParameters, ONLY : WaterEquilibriumIndex
	IMPLICIT NONE
	INTEGER :: R
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8  :: RH, V, Ratio
	LOGICAL :: Solid
	ReactionViolation = .FALSE. ; V = 0.
	IF (R .EQ. WaterEquilibriumIndex) RETURN	!water: the solvers leave it once OH- or H+ is negligible
	IF (DABS(AqEquilibriaList(R,9) + 1.) .LT. 1.0D-6) RETURN
	Solid = (DABS(AqEquilibriaList(R,6) - 1.) .LT. 1.0D-5)
	IF (Solid) THEN
		IF (RH .GE. DeliquesenceRH(R, InParticle)) RETURN
	END IF
	Ratio = EquilibriumConstantsRatio (R, InParticle)
	IF (.NOT. (Ratio .GT. 0. .AND. Ratio .LT. 1.0D11)) RETURN
	ReactionViolation = .TRUE.
	V = ABS(1. - Ratio)
	IF (AqEquilibriaList(R,1) .GT. 0.) THEN
		IF (InParticle%AqChems(INT(AqEquilibriaList(R,1))) .LE. 1.0D-28) V = MAX(0.D0, Ratio - 1.)
	END IF
	!an absent solid that would complete the unstable side of a solid circuit must not precipitate
	!(with consistent activities it would be undersaturated); not a violation
	IF (Solid .AND. V .GT. 0. .AND. Ratio .GT. 1.) THEN
		IF (InParticle%AqChems(INT(AqEquilibriaList(R,1))) .LE. 1.0D-28) THEN
			IF (CompletesUnstableSet(R, InParticle)) V = 0.
		END IF
	END IF
END FUNCTION ReactionViolation

!For every reaction that violates the equilibrium conditions by more than 2 x tolerance in a
!result, apply the single-reaction solver to a copy (Tmp) and record whether it moves the state
!and whether that fixes the reaction.
SUBROUTINE ViolationDiag (InParticle, Tmp)
	USE Chemistry, ONLY : HowManyAqEqReactions
	USE GridPointFields, ONLY : GetRelativeHumidity
	USE ModelParameters, ONLY : AqThermoEquilibriumError
	IMPLICIT NONE
	TYPE(PARTICLE),POINTER :: InParticle, Tmp
	INTEGER :: R, Ret
	REAL*8  :: RH, V, V2
	RH = GetRelativeHumidity ()
	DO R = 1, HowManyAqEqReactions
		IF (.NOT. ReactionViolation (R, InParticle, RH, V)) CYCLE
		IF (V .LE. 2.*AqThermoEquilibriumError) CYCLE
		VDViol(R) = VDViol(R) + 1 ; VDSize(R) = VDSize(R) + V
		Tmp%AqChems = InParticle%AqChems ; Tmp%GammaMixed = InParticle%GammaMixed
		Tmp%NumberOfParticles = InParticle%NumberOfParticles ; Tmp%Dry = InParticle%Dry
		Tmp%Temperature = InParticle%Temperature ; Tmp%WaterActivity = InParticle%WaterActivity
		Tmp%EmbryoRadius = InParticle%EmbryoRadius ; Tmp%SolutionDensity = InParticle%SolutionDensity
		Tmp%SurfaceTension = InParticle%SurfaceTension ; Tmp%IonicStr = InParticle%IonicStr
		CALL EquilibrateAqueousDissociationReaction (R, Tmp, OptErrorTolerance = AqThermoEquilibriumError, ReturnType = Ret)
		IF (ANY(Tmp%AqChems .NE. InParticle%AqChems)) THEN
			VDMoved(R) = VDMoved(R) + 1
			V2 = ABS(1. - EquilibriumConstantsRatio(R, Tmp))
			IF (V2 .LE. AqThermoEquilibriumError) VDFixed(R) = VDFixed(R) + 1
		ELSE
			VDStuck(R) = VDStuck(R) + 1
		END IF
	END DO
END SUBROUTINE ViolationDiag

!Number of circuits whose unstable side is entirely present in the particle
INTEGER FUNCTION CountUnstableSets (InParticle)
	USE Chemistry, ONLY : AqEquilibriaList
	IMPLICIT NONE
	TYPE(PARTICLE),POINTER :: InParticle
	INTEGER :: L, I
	REAL*8  :: S
	LOGICAL :: AllPresent
	IF (.NOT. ALLOCATED(CircN)) CALL InitSolidCircuits (SIZE(InParticle%AqChems))
	CountUnstableSets = 0
	DO L = 1, NCirc
		S = 0.
		DO I = 1, CircN(L)
			S = S + CircCoef(I,L) * LOG(EqCoeff(CircRxn(I,L), InParticle, 1))
		END DO
		IF (ABS(S) .LT. 1.0D-8) CYCLE
		IF (.NOT. CircuitApplies(L, InParticle)) CYCLE
		IF (.NOT. StableSideAvailable(L, S, InParticle)) CYCLE
		AllPresent = .TRUE.
		DO I = 1, CircN(L)
			IF (CircCoef(I,L)*S .LE. 0.) CYCLE
			IF (InParticle%AqChems(INT(AqEquilibriaList(CircRxn(I,L),1))) .LE. 1.0D-28) AllPresent = .FALSE.
		END DO
		IF (AllPresent) CountUnstableSets = CountUnstableSets + 1
	END DO
END FUNCTION CountUnstableSets

!Replay the dump DumpFile through FindElectrolyteEquilibriumCore with the solver selected now,
!time only the solver calls, compare each result with the dumped result and write FEEReplay.txt.
!ResultFile (if not blank) receives each replayed result; VDiag adds the per-reaction
!violation diagnosis; NRep repeats the replay (for timing).  The model must be initialised
!(chemistry, particles) first: Particles%First is used as the work particle.
SUBROUTINE ReplayFEE (DumpFile, ResultFile, VDiag, NRep)
	USE GridPointFields, ONLY : GetPress, SetTempField, SetPressField, ReplaceGridCellChemBurden
	USE ModelParameters, ONLY : PressureScale, AqThermoEquilibriumError, AqInnerSolver, AqDevCache
	IMPLICIT NONE
	CHARACTER(LEN=*) :: DumpFile, ResultFile
	LOGICAL :: VDiag
	INTEGER :: NRep
	INTEGER :: IU, Kind, Flags, RT, RTref, IOS, N, NSame, NRT, OU
	INTEGER :: IRep, RU, KB, KHist(5)
	INTEGER*8 :: K0
	INTEGER :: NUnstRef, NUnstOut, NUnstIn
	INTEGER :: NVR, NVO, NVRWet, NVOWet, NDryRef, NDryOut
	INTEGER :: IR2
	REAL*8 :: VR, VO, VRSum, VOSum
	TYPE(PARTICLE),POINTER :: P, Ref, Pre
	REAL*8, ALLOCATABLE :: SolveG(:)
	LOGICAL :: UT
	REAL*8 :: T, Pr, H2O, Dum1, Dum2, Dum3, C0, C1, TSolve, D, DMax, DSum, Scale
	REAL*8, ALLOCATABLE :: Out(:)

	P => Particles%First
	ALLOCATE (Ref)
	ALLOCATE (Ref%AqChems(SIZE(P%AqChems)), Ref%GammaMixed(SIZE(P%GammaMixed)), &
	          Ref%OrgChems(SIZE(P%OrgChems)), Ref%AqOrgChems(SIZE(P%AqOrgChems)))
	ALLOCATE (Out(SIZE(P%AqChems)))
	ALLOCATE (Pre)
	ALLOCATE (Pre%AqChems(SIZE(P%AqChems)), Pre%GammaMixed(SIZE(P%GammaMixed)), &
	          Pre%OrgChems(SIZE(P%OrgChems)), Pre%AqOrgChems(SIZE(P%AqOrgChems)), SolveG(SIZE(P%GammaMixed)))
	OPEN (NEWUNIT=IU, FILE=TRIM(DumpFile), FORM="UNFORMATTED", ACCESS="STREAM", STATUS="OLD")
	N = 0 ; NSame = 0 ; NRT = 0 ; TSolve = 0. ; DMax = 0. ; DSum = 0.
	FEECountEADR = 0 ; FEECountOuter = 0 ; FEECountECREval = 0 ; KHist = 0 ; FEECountNewtonFail = 0
	FEECountNewtonSing = 0 ; FEECountNewtonMaxIt = 0 ; FEECountNewtonDrop = 0 ; FEECountNewtonCycles = 0
	FEECountPairRefuse = 0 ; NUnstRef = 0 ; NUnstOut = 0 ; NUnstIn = 0 ; FEECountResolved = 0 ; FEECountDry = 0 ; FEECountPreAssoc = 0
	NVR = 0 ; NVO = 0 ; VRSum = 0. ; VOSum = 0. ; NVRWet = 0 ; NVOWet = 0 ; NDryRef = 0 ; NDryOut = 0
	IF (LEN_TRIM(ResultFile) .GT. 0) OPEN (NEWUNIT=RU, FILE=TRIM(ResultFile), FORM="UNFORMATTED", &
	                                       ACCESS="STREAM", STATUS="REPLACE")
	DO IRep = 1, NRep
	REWIND (IU)
	DO
		CALL ReadFEEState (IU, P, Kind, Flags, RT, T, Pr, H2O, IOS)
		IF (IOS .NE. 0) EXIT
		IF (Kind .NE. 1) STOP "ReplayFEE: record out of order"
		CALL SetTempField (T)
		IF (GetPress() .NE. Pr) CALL SetPressField (Pr / PressureScale)
		IF (GetPress() .NE. Pr) STOP "ReplayFEE: cannot restore pressure exactly"
		CALL ReplaceGridCellChemBurden (1, H2O)
		IF (IRep .EQ. 1 .AND. CountUnstableSets(P) .GT. 0) NUnstIn = NUnstIn + 1
		!Calls that update the activity coefficients solve with KusikMeissner(input) and recompute them
		!at the end; the equilibrium checks below use the coefficients the solve used.
		UT = (MOD(Flags,2) .EQ. 1)
		IF (UT .AND. IRep .EQ. 1) THEN
			Pre%AqChems = P%AqChems ; Pre%GammaMixed = P%GammaMixed ; Pre%OrgChems = P%OrgChems
			Pre%AqOrgChems = P%AqOrgChems ; Pre%NumberOfParticles = P%NumberOfParticles ; Pre%Dry = P%Dry
			Pre%Temperature = P%Temperature ; Pre%WaterActivity = P%WaterActivity ; Pre%IonicStr = P%IonicStr
			Pre%EmbryoRadius = P%EmbryoRadius ; Pre%SolutionDensity = P%SolutionDensity
			Pre%SurfaceTension = P%SurfaceTension
			CALL KusikMeissner (Pre)
			SolveG = Pre%GammaMixed
		END IF
		K0 = FEECountOuter
		CALL CPU_TIME (C0)
		CALL FindElectrolyteEquilibriumCore (P, MOD(Flags,2) .EQ. 1, Flags .GE. 2, RT)
		CALL CPU_TIME (C1)
		TSolve = TSolve + (C1 - C0)
		KB = 1 ; IF (FEECountOuter-K0 .GT. 10) KB = 2 ; IF (FEECountOuter-K0 .GT. 30) KB = 3
		IF (FEECountOuter-K0 .GT. 60) KB = 4 ; IF (FEECountOuter-K0 .GT. 100) KB = 5
		KHist(KB) = KHist(KB) + 1
		Out = P%AqChems
		IF (LEN_TRIM(ResultFile) .GT. 0 .AND. IRep .EQ. 1) WRITE (RU) Out, RT
		CALL ReadFEEState (IU, Ref, Kind, Flags, RTref, Dum1, Dum2, Dum3, IOS)
		IF (IOS .NE. 0 .OR. Kind .NE. 2) STOP "ReplayFEE: missing result record"
		N = N + 1
		IF (ALL(TRANSFER(Out,0_8,SIZE(Out)) .EQ. TRANSFER(Ref%AqChems,0_8,SIZE(Out)))) NSame = NSame + 1
		IF (RT .EQ. RTref) NRT = NRT + 1
		IF (IRep .EQ. 1) THEN
			Ref%Temperature = P%Temperature
			IF (UT) THEN
				Ref%GammaMixed = SolveG ; P%GammaMixed = SolveG
			END IF
			IF (CountUnstableSets(Ref) .GT. 0) NUnstRef = NUnstRef + 1
			IF (CountUnstableSets(P) .GT. 0) NUnstOut = NUnstOut + 1
			VR = EqViolation(Ref) ; VO = EqViolation(P)
			VRSum = VRSum + VR ; VOSum = VOSum + VO
			IF (VR .GT. 2.*AqThermoEquilibriumError) NVR = NVR + 1
			IF (VO .GT. 2.*AqThermoEquilibriumError) NVO = NVO + 1
			!equilibrium conditions are only meaningful with a solution phase; count wet results
			IF (P%AqChems(1) .LT. 1.0D-3 * SUM(P%AqChems(2:))) THEN
				NDryOut = NDryOut + 1
			ELSE IF (VO .GT. 2.*AqThermoEquilibriumError) THEN
				NVOWet = NVOWet + 1
			END IF
			IF (Ref%AqChems(1) .LT. 1.0D-3 * SUM(Ref%AqChems(2:))) THEN
				NDryRef = NDryRef + 1
			ELSE IF (VR .GT. 2.*AqThermoEquilibriumError) THEN
				NVRWet = NVRWet + 1
			END IF
			IF (VDiag) CALL ViolationDiag (P, Ref)
		END IF
		!difference of each amount, relative to the particle's total dissolved moles
		Scale = SUM(ABS(Ref%AqChems(2:))) + 1.0D-300
		D = MAXVAL(ABS(Out(2:) - Ref%AqChems(2:))) / Scale
		DMax = MAX(DMax, D) ; DSum = DSum + D
	END DO
	END DO
	IF (LEN_TRIM(ResultFile) .GT. 0) CLOSE (RU)
	CLOSE (IU)
	OPEN (NEWUNIT=OU, FILE="FEEReplay.txt", STATUS="REPLACE")
	WRITE (OU, '(A,I10)')     "records          ", N
	WRITE (OU, '(A,I10)')     "bitwise same     ", NSame
	WRITE (OU, '(A,I10)')     "same ReturnType  ", NRT
	WRITE (OU, '(A,ES12.4)')  "solver CPU s     ", TSolve
	WRITE (OU, '(A,ES12.4)')  "max  |d|/total   ", DMax
	WRITE (OU, '(A,ES12.4)')  "mean |d|/total   ", DSum / MAX(N,1)
	WRITE (OU, '(A,L2)')      "newton inner     ", AqInnerSolver .EQ. 1
	WRITE (OU, '(A,L2)')      "deviation cache  ", AqDevCache .EQ. 1
	WRITE (OU, '(A,I10)')     "repeats          ", NRep
	WRITE (OU, '(A,F10.2)')   "outer iter/call  ", DBLE(FEECountOuter)/MAX(N,1)
	WRITE (OU, '(A,F10.2)')   "1-rxn solves/call", DBLE(FEECountEADR)/MAX(N,1)
	WRITE (OU, '(A,F10.2)')   "ratio evals/call ", DBLE(FEECountECREval)/MAX(N,1)
	WRITE (OU, '(A,5I8)')     "outer iter <=10,30,60,100,>100", KHist
	WRITE (OU, '(A,4I10)')    "coupled: fallbacks, singular, max iter, cycles", FEECountNewtonFail, FEECountNewtonSing, &
	                          FEECountNewtonMaxIt, FEECountNewtonCycles
	WRITE (OU, '(A,I10)')     "coupled: DRH-ordered drops accepted", FEECountNewtonDrop
	WRITE (OU, '(A,I10)')     "coupled: dry particles put into solids", FEECountDry
	WRITE (OU, '(A,I10)')     "coupled: dissolved electrolytes set to equilibrium first", FEECountPreAssoc
	WRITE (OU, '(A,I10)')     "precipitations refused (pair rule)", FEECountPairRefuse
	WRITE (OU, '(A,I10)')     "calls with an unstable set resolved", FEECountResolved
	WRITE (OU, '(A,I6)')      "solid circuits", NCirc
	WRITE (OU, '(A,3I10)')    "records with an unstable solid set: input, reference, replay", NUnstIn, NUnstRef, NUnstOut
	WRITE (OU, '(A,2ES12.4)') "mean equilibrium violation: reference, replay", VRSum/MAX(N,1), VOSum/MAX(N,1)
	WRITE (OU, '(A,2I10)')    "records violating by > 2 x tolerance: reference, replay", NVR, NVO
	WRITE (OU, '(A,2I10)')    "dry results (water < 1e-3 of electrolyte): reference, replay", NDryRef, NDryOut
	WRITE (OU, '(A,2I10)')    "wet results violating by > 2 x tolerance: reference, replay", NVRWet, NVOWet
	IF (VDiag) THEN
		WRITE (OU, '(A)') "violations in the replay result, per reaction: rxn, count, mean size, solver moves it, fixes it, stuck"
		DO IR2 = 1, 200
			IF (VDViol(IR2) .GT. 0) WRITE (OU, '(I4,I8,ES11.3,3I8)') IR2, VDViol(IR2), VDSize(IR2)/VDViol(IR2), &
				VDMoved(IR2), VDFixed(IR2), VDStuck(IR2)
		END DO
	END IF
	CLOSE (OU)
END SUBROUTINE ReplayFEE

END MODULE FEEHarness
