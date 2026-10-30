!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ParticleAcidity.h                                                         !!
!!                                                                           !!
!! FillAcidState: for every particle (in particle-list order, the order of   !!
!! HeteroRadii / HeteroNumberConc) the acid-catalysed epoxide ring-opening   !!
!! rate  kp = sum_i k_i [acid_i] [nucleophile_i]  [1/s]  and the fraction of !!
!! it that goes to organosulfate.  Used by rate flag 18 through the hook     !!
!! Chemistry:AcidStateHook, which SetOrganicDissolution registers.           !!
!!                                                                           !!
!! Concentrations are molarities of the particle's aqueous solution          !!
!! (mol/L, same basis as AerosolPH): H+, HSO4-, SO4--, NO3-, water.          !!
!! Rate constants [M-1 s-1 per M of nucleophile] from CMAQ AERO6/7           !!
!! (Eddingsaas et al. 2010; Piletic et al. 2013):                            !!
!!                   H+ acid        HSO4- acid                               !!
!!   water          9.0e-4         1.31e-5                                   !!
!!   sulfate        8.83e-3 (IEPOX), 2.0e-4 (IMAE)   2.92e-6                 !!
!!   nitrate        2.0e-4         2.92e-6          (IEPOX kind 1 only)      !!
!! kinds: 1 IEPOX (AERO6 / CRACMM), 2 IMAE, 3 IEPOX (AERO7i, no nitrate).    !!
!! Sulfate is SO4-- only (CMAQ subtracts HSO4- from total sulfate).          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FillAcidState (NumberBins, KPart, FracOS)

	USE ModelParameters, ONLY : ProtonIndex
	USE Chemistry,       ONLY : FindChem
	USE Aerosols,        ONLY : Particles, Particle, AqueousSolutionMass

	IMPLICIT NONE

	INTEGER, INTENT(in)  :: NumberBins
	REAL*8,  INTENT(out) :: KPart(3,NumberBins), FracOS(3,NumberBins)

	REAL*8, PARAMETER :: KH_WATER = 9.0D-4, KH_SO4 = 8.83D-3, KH_NUC = 2.0D-4
	REAL*8, PARAMETER :: KHSO4_WATER = 1.31D-5, KHSO4_NUC = 2.92D-6

	LOGICAL, SAVE :: First = .TRUE.
	INTEGER, SAVE :: IHSO4 = 0, ISO4 = 0, INO3 = 0
	TYPE(Particle), POINTER :: Cur
	INTEGER :: I, K
	REAL*8  :: Mass, Fac, H, HS, S, N, W, os, other

	IF (First) THEN
		First = .FALSE.
		IHSO4 = FindChem("HSO4-", 1, .TRUE.)
		ISO4  = FindChem("SO4--", 1, .TRUE.)
		INO3  = FindChem("NO3-",  1, .TRUE.)
	END IF

	KPart  = 0.0D0
	FracOS = 0.0D0

	Cur => Particles%First
	I = 0
	DO WHILE (ASSOCIATED(Cur) .AND. I .LT. NumberBins)
		I = I + 1
		IF (Cur%NumberOfParticles .GT. 1.0D-40 .AND. .NOT. Cur%Dry) THEN
			Mass = AqueousSolutionMass(Cur)
			IF (Mass .GT. 0.0D0) THEN
				Fac = Cur%SolutionDensity*1000.0D0/Mass		!mol/L per (mol per particle)
				H = 0.0D0; HS = 0.0D0; S = 0.0D0; N = 0.0D0
				IF (ProtonIndex .GT. 0) H  = MAX(Cur%AqChems(ProtonIndex), 0.0D0)*Fac
				IF (IHSO4 .GT. 0)       HS = MAX(Cur%AqChems(IHSO4), 0.0D0)*Fac
				IF (ISO4  .GT. 0)       S  = MAX(Cur%AqChems(ISO4),  0.0D0)*Fac
				IF (INO3  .GT. 0)       N  = MAX(Cur%AqChems(INO3),  0.0D0)*Fac
				W = MAX(Cur%AqChems(1), 0.0D0)*Fac

				!kind 1: IEPOX, water + sulfate + nitrate
				os    = KH_SO4*H*S + KHSO4_NUC*HS*S
				other = KH_WATER*H*W + KH_NUC*H*N + KHSO4_WATER*HS*W + KHSO4_NUC*HS*N
				KPart(1,I) = os + other
				IF (KPart(1,I) .GT. 0.0D0) FracOS(1,I) = os/KPart(1,I)
				!kind 3: IEPOX, water + sulfate
				other = KH_WATER*H*W + KHSO4_WATER*HS*W
				KPart(3,I) = os + other
				IF (KPart(3,I) .GT. 0.0D0) FracOS(3,I) = os/KPart(3,I)
				!kind 2: IMAE
				os    = KH_NUC*H*S + KHSO4_NUC*HS*S
				other = KH_WATER*H*W + KHSO4_WATER*HS*W
				KPart(2,I) = os + other
				IF (KPart(2,I) .GT. 0.0D0) FracOS(2,I) = os/KPart(2,I)
			END IF
		END IF
		Cur => Cur%Next
	END DO

END SUBROUTINE FillAcidState
