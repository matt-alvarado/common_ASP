!! ASP (c), 2004-2013, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Coagulation.f90
!! This file contains the coagulation numerical integration routine 
!! for the internally mixed, moving-center size distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 08/26  2006   Matt Alvarado     Fixed bug in CollisionEfficiency          !!
!! 07/26  2007   Matt Alvarado     Added RecalculateRadius to                !!
!!                                  monodisperse coag	                     !!
!! 02/17  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 24670d0: CoagTest H1: precompute   !!
!!                                 single-particle kernel quantities once    !!
!!                                 per section                               !!
!! 10/28  2026   Claude Sonnet 5.5 Commit c20e470: CoagTest H4: collect      !!
!!                                 species-independent coagulation terms     !!
!!                                 once per step                             !!
!! 10/28  2026   Claude Sonnet 5.5 Commit cf81a74: CoagTest H6: skip the     !!
!!                                 duplicate RecalculateRadius pass when     !!
!!                                 Coag calls RegridAerosol                  !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 89d60a6: CoagTest H5: build        !!
!!                                 VolumeFractions with one bracket scan per !!
!!                                 (I,J)                                     !!
!! 11/02  2026   Claude Opus 5.5   Patch 37: regrid before the semi-      !!
!!                                 implicit coagulation step when the     !!
!!                                 occupied sections are out of size order!!
!!                                 (fixes the 'Summed Total Volume        !!
!!                                 concentration is lost' stop)           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES				    !!
!! 1. ModelParameters			    !!
!! 2. InfrastructuralCode		    !!
!! 3. GridPointFields			    !!
!! 4. Chemistry				    !!
!! 5. Aerosols				    !!
!! 6. OutputRoutines			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES						    !!
!! 1. CoagulationKernels.h					    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. REAL*8 FUNCTION CollisionEfficiency (Particle1, Particle2)	!!
!! 2. SUBROUTINE StepSectionalCoagulationJacobson (TimeStep)            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE Coagulation

	PRIVATE

	REAL*8, PARAMETER :: CE = 1.0 !Coagulation Efficiency
					!(fraction of collisions where
					!particles stick together)
	
	PUBLIC :: StepSectionalCoagulationJacobson
	PUBLIC :: CoagRegridCount
	INTEGER*8, SAVE :: CoagRegridCount = 0	!Patch 37: regrids done before coagulation
			  
CONTAINS

INCLUDE "CoagulationKernels.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Collision Efficiency 
!! using formula from Jacobson, "Fundamentals of Atmospheric Modeling"
!! , 2nd ed., p. 510 (2005)
REAL*8 FUNCTION CollisionEfficiency (Particle1, Particle2)

	USE Aerosols,	ONLY : Particle, TerminalVelocity, ReynoldsNumber

	USE ModelParameters, ONLY : g
	USE InfrastructuralCode, ONLY : Error
	IMPLICIT NONE

	!! External Variables
	TYPE(Particle), POINTER :: Particle1, Particle2

	!! Internal Variables
	REAL*8  :: Stokes, Reynolds, Ev, Ea

	!From Jacobson, "Fundamentals of Atmospheric Modeling", p. 510
	IF(Particle2%EffectiveRadius .GE. Particle1%EffectiveRadius) THEN
		Stokes = TerminalVelocity(Particle1) &
                         *ABS(TerminalVelocity(Particle2) - TerminalVelocity(Particle1)) &
                         /Particle2%EffectiveRadius/g
		Reynolds = ReynoldsNumber(Particle2)
		Ea = (Stokes**2.0)/((Stokes+0.5)**2.0)
		IF(Stokes .LE. 1.214) THEN
			Ev = 0.
		ELSE
			Ev = (1 + (0.75*LOG(2*Stokes))/(Stokes-1.214) )**(-2.0)
		END IF

		CollisionEfficiency = (60*Ev + Ea*Reynolds)/(60 + Reynolds)

	ELSE
		Stokes = TerminalVelocity(Particle2)*ABS(TerminalVelocity(Particle1) - TerminalVelocity(Particle2))/Particle1%EffectiveRadius/g
		Reynolds = ReynoldsNumber(Particle1)
		Ea = (Stokes**2.0)/((Stokes+0.5)**2.0)
		IF(Stokes .LE. 1.214) THEN
			Ev = 0.
		ELSE
			Ev = (1 + (0.75*LOG(2*Stokes))/(Stokes-1.214) )**(-2.0)
		END IF

		CollisionEfficiency = (60*Ev + Ea*Reynolds)/(60 + Reynolds)
	
	
	END IF

	IF (CollisionEfficiency .GT. 1.0 .OR. CollisionEfficiency .LT. 0.0) THEN
		WRITE(*,*) CollisionEfficiency
		CALL ERROR("CollisionEfficiency in Coagulation.f90 is greater than 1 or less than 0.")
	END IF
END FUNCTION CollisionEfficiency 

SUBROUTINE StepSectionalCoagulationJacobson (TimeStep)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Matt Alvarado created this by simplifying StepCoagulation to only 
!! work for a sectional distribution and to always use mean values.
!! This follows the semi-implicit method of Jacobson - see
!! Jacobson, "Fundamentals of Atmospheric Modeling"
!! , 2nd ed., Section 15.2, (2005)

	USE InfrastructuralCode, ONLY : INT2STR, ERROR

	USE GridPointFields,     ONLY : GetDynamicViscosityOfAir,	&
					GetTemp, GetAirDensity

	USE Aerosols,            ONLY : Particles, Particle,		&
				SortAerosolAtGridPointForCoagulation,	&
					RecalculateRadius,		&
					RegridAerosol

	USE ModelParameters,	  ONLY : ThermoBulkMode,		  &
					 Pi, &
					 ConstantKernel, &
					 DoCoagulationFlag, &
					 Monodisperse, &
					 CoagKernel

	USE Chemistry,	          ONLY : HowManyAqChems,		&
					 HowManyAqCations,	&
				         HowManyAqAnions, &
					 AqMolecularMass, &
					 HowManyOrgChems, &
					 OrgMolecularMass, &
					 OrgDensity, &
					 HowManyAqOrgChems, &
					 AqOrgMolecularMass, &
					 AqOrgDensity, &
					 SolidSaltDensity, &
					 AqPhaseChemicalNames
	USE OutputRoutines

	IMPLICIT NONE

	!! External Variables
	REAL*8, INTENT(IN)  :: TimeStep

	!! Internal Variables
	INTEGER	:: HowManyLinks, I, J, K, L, M, Q, Count, JJ

	REAL*8	:: AA, BB, II, LL, Kfac, &
		   BigV, LowerSum, UpperSum, InnerSum, TotalAqOrgVolConc, & 
                   TotalOrgVolConc, TotalAqVolConc, &
		   InVolConc, OutVolConc, SolnSum, AvgSolutionDensity, & 
                   OutSumVolumeConc, InSumVolumeConc, StoreDensity, Monokernel

	REAL*8, ALLOCATABLE	:: CoagArray(:,:), SingleParticleVolume(:), &
                                   VolumeConc(:), &
				   VolumeFractions(:,:,:), NumConc(:), &
                                   AqVolConc(:,:), &
				   OrgVolConc(:,:), AqOrgVolConc(:,:), &
                                   InitNumConc(:), SumVolConc(:)

	!! Single-particle kernel quantities, computed once per section per step
	REAL*8, ALLOCATABLE	:: PreR(:), PreD(:), PreMass(:), PreDel(:), PreRe(:), PreVt(:)
	LOGICAL, ALLOCATABLE	:: PreSect(:)
	REAL*8	:: PreDynVisc, PreTemp, PreAirDensity
	!! Species-independent coagulation terms (collected once per step)
	REAL*8, ALLOCATABLE	:: LowerK(:), EntW(:)
	INTEGER, ALLOCATABLE	:: KStart(:), EntI(:), EntJ(:)
	LOGICAL, ALLOCATABLE	:: EntLast(:)
	INTEGER	:: NEnt, NGroupStart, N
	REAL*8	:: SPVK
	INTEGER :: NRegrid	!Patch 37
	LOGICAL :: OutOfOrder
	REAL*8  :: PrevRadius

	TYPE(Particle), POINTER :: Particle1, Particle2, ParticleK, & 
                                   Runner, ReplicaPlaceholder, &
		DifferenceParticle1, DifferenceParticle2, DifferenceRunner , &
			           ParticleI, ParticleJ

	!IF main input deck said to skip coagulation, skip it
	IF(.NOT.DoCoagulationFlag) RETURN
	
	!! If we're in bulk mode, coagulation doesn't do anything
	IF (ThermoBulkMode) RETURN

	!! Sort the gridcell, sections in order and first
	CALL SortAerosolAtGridPointForCoagulation ()

	!! Patch 37: the semi-implicit scheme needs the sections in ascending order of their actual
	!! single-particle volume (the volume partition brackets v_i + v_j between neighbouring sections,
	!! and gains to section K use the already-updated smaller sections).  Condensation can move a
	!! section's particles past its neighbour; the volume check then stops the run ("Summed Total
	!! Volume concentration is lost during step").  If the occupied sections are out of order, regrid
	!! first (RegridAerosol moves particles at most one section per call, so repeat up to 10 times).
	DO NRegrid = 1, 10
		OutOfOrder = .FALSE. ; PrevRadius = -1.0
		Particle1 => Particles%First
		DO WHILE (ASSOCIATED(Particle1))
			IF (Particle1%NumberOfParticles .GT. 0.0D0) THEN
				IF (Particle1%EffectiveRadius .LT. PrevRadius) OutOfOrder = .TRUE.
				PrevRadius = Particle1%EffectiveRadius
			END IF
			Particle1 => Particle1%Next
		END DO
		IF (.NOT. OutOfOrder) EXIT
		CALL RegridAerosol ()
		CALL SortAerosolAtGridPointForCoagulation ()
		CoagRegridCount = CoagRegridCount + 1
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! If there are no aerosol at the gridcell, !!
	!! then there is nothing to do.  Return.    !!
	Particle1 => Particles%First !!!!!!!





5	IF (Particle1%NumberOfParticles .EQ. 0.0D0) THEN
		IF (ASSOCIATED(Particle1%Next)) THEN
			Particle1 => Particle1%Next
			GOTO 5
		END IF
		RETURN
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Everything should now be sorted           !!
	!! Back up to the head of the aerosol list   !!
	Particle1 => Particles%First !!!!!!!!
	Particle2 => Particle1
	Runner    => Particle1 !! In the second pass, this may need an initial position


	!! COUNT the number of LINKS in the gridcell (total number of sections)
	HowManyLinks =  0
15	HowManyLinks = HowManyLinks + 1

	IF(ASSOCIATED(Particle2%Next)) THEN
		Particle2 => Particle2%Next
		GOTO 15
	END IF
	Particle2 => Particle1


	!! BUILD N VECTOR for Single particle volumes
	ALLOCATE (SingleParticleVolume(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate SingleParticleVolume at size "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N VECTOR for Volume Concentrations
	ALLOCATE (VolumeConc(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate VolumeConc at size "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N VECTOR for Number Concentrations
	ALLOCATE (NumConc(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate NumConc at size "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N VECTOR for Initial Number Concentrations
	ALLOCATE (InitNumConc(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate NumConc at size "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N VECTOR for Summed Volume Concentrations
	ALLOCATE (SumVolConc(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate NumConc at size "//TRIM(INT2STR(HowManyLinks)))

	
	!! BUILD N by N ARRAY for COAGULATION KERNELS
	ALLOCATE (CoagArray(HowManyLinks,HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate CoagArray at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N by M ARRAY for Species Volume Concentrations
	ALLOCATE (AqVolConc(HowManyLinks, HowManyAqChems+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate AqVolConc at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyAqChems+HowManyAqCations+HowManyAqAnions)))

	!! BUILD N by M ARRAY for Species Volume Concentrations
	ALLOCATE (OrgVolConc(HowManyLinks, HowManyOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate OrgVolConc at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyOrgChems)))

	!! BUILD N by M ARRAY for Species Volume Concentrations
	ALLOCATE (AqOrgVolConc(HowManyLinks, HowManyAqOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate AqOrgVolConc at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyAqOrgChems)))

	
	!! BUILD N by N by N ARRAY for VolumeFractions
	ALLOCATE (VolumeFractions(HowManyLinks, HowManyLinks, HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepCoagulation(), couldn't allocate VolumeFractions at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyLinks))//" by "//TRIM(INT2STR(HowManyLinks)))

	

	!! If we're in monodisperse mode, do monodisperse coag
	IF (Monodisperse) THEN
		Particle1 => Particles%First
		Particle2 => Particle1
		IF(ConstantKernel) THEN
			MonoKernel = CoagKernel	
		ELSE
			MonoKernel = (CoagKernBrownianDiff   (Particle1, Particle2) &
					                  + CoagKernGravCollection (Particle1, Particle2)) !! Add Brownian and Convective Contributions
		END IF
		
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(1,I) = Particle1%AqChems(I)*Particle1%NumberofParticles !mol/cm3 air
		END DO
		
		DO I = 1, HowManyOrgChems
			OrgVolConc(1,I) = Particle1%OrgChems(I)*Particle1%NumberofParticles !mol/cm3 air
		END DO

		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(1,I) = Particle1%AqOrgChems(I)*Particle1%NumberofParticles !mol/cm3 air
		END DO

		!Explicit coagulation integration
		Particle1%NumberofParticles = Particle1%NumberofParticles - Timestep*MonoKernel*Particle1%NumberofParticles**2
		
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			Particle1%AqChems(I) = AqVolConc(1,I)/Particle1%NumberofParticles !mol/particle
		END DO
		
		DO I = 1, HowManyOrgChems
			Particle1%OrgChems(I) = OrgVolConc(1,I)/Particle1%NumberofParticles !mol/particle
		END DO

		DO I = 1, HowManyAqOrgChems
			Particle1%AqOrgChems(I) = AqOrgVolConc(1,I)/Particle1%NumberofParticles !mol/particle
		END DO
		
		!Does radius and density
		CALL RecalculateRadius(Particle1)		
			
		RETURN		
	END IF


	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! LOOP AND FILL CoagArray with Kernel Values   !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Quantities that depend on a single particle are computed once per
	!! section here, rather than inside the N*N pair loop.
	ALLOCATE (PreR(HowManyLinks), PreD(HowManyLinks), PreMass(HowManyLinks), PreDel(HowManyLinks), &
	          PreRe(HowManyLinks), PreVt(HowManyLinks), PreSect(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate single-particle kernel arrays at size "//TRIM(INT2STR(HowManyLinks)))
	PreSect = .FALSE.
	PreR = 0.; PreD = 0.; PreMass = 0.; PreDel = 0.; PreRe = 0.; PreVt = 0.
	IF (.NOT. ConstantKernel) THEN
		PreDynVisc    = GetDynamicViscosityOfAir ()
		PreTemp       = GetTemp ()
		PreAirDensity = GetAirDensity ()
		Particle1 => Particles%First
		DO I = 1, HowManyLinks
			IF (Particle1%NumberOfParticles .NE. 0.0D0 .AND. Particle1%EffectiveRadius .NE. 0. .AND. Particle1%Sectional) THEN
				PreSect(I) = .TRUE.
				PreR(I)    = Particle1%EffectiveRadius
				CALL CoagPrecomputeParticle (Particle1, PreDynVisc, PreTemp, PreD(I), PreMass(I), PreDel(I), PreRe(I), PreVt(I))
			END IF
			Particle1 => Particle1%Next
		END DO
	END IF
	Particle1 => Particles%First
	Particle2 => Particles%First




	DO I = 1, HowManyLinks
		DO J = 1, HowManyLinks

			!! If either of these is zero we don't have to actually calculate anything
			IF (Particle1%NumberOfParticles .EQ. 0.0D0  .OR.  Particle2%NumberOfParticles .EQ. 0.0D0) THEN
				CoagArray(I,J) = 0.
				GOTO 40
			END IF

			IF (Particle1%EffectiveRadius   .EQ. 0. .OR. Particle2%EffectiveRadius   .EQ. 0. .OR. &
				Particle1%NumberOfParticles .EQ. 0.0D0 .OR. Particle2%NumberOfParticles .EQ. 0.0D0) THEN

				CoagArray(I,J) = 0.

			ELSE

				
				IF(ConstantKernel) THEN
					CoagArray(I,J) = CoagKernel

				ELSE IF (PreSect(I) .AND. PreSect(J)) THEN
					!! Both sectional: use the per-section quantities computed above
					CoagArray(I,J) = CoagKernPairSectional (PreR(I), PreD(I), PreMass(I), PreDel(I), PreRe(I), PreVt(I), &
					                                        PreR(J), PreD(J), PreMass(J), PreDel(J), PreRe(J), PreVt(J), &
					                                        PreDynVisc, PreTemp, PreAirDensity)

				ELSE
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					!! Add together compenents of the Coagulation Kernel (units cm3 part.-1 s-1) !!
					CoagArray(I,J) = (CoagKernBrownianDiff   (Particle1, Particle2) &
					                  + CoagKernGravCollection (Particle1, Particle2)) !! Add Brownian and Convective Contributions
					!write(*,*) 'CRL: CoagKernBrownianDiff   (Particle1, Particle2) =',CoagKernBrownianDiff   (Particle1, Particle2)
					!write(*,*) 'CRL: CoagKernGravCollection (Particle1, Particle2) =',CoagKernGravCollection (Particle1, Particle2)
					!write(*,*) ''

				END IF
			END IF
			
			IF(CoagArray(I,J) .LT. 0.0) WRITE(*,*) "CoagArray", I, J, CoagArray(I,J)
				
			!! Step the second particle (always either smaller than the larger or with a lower distribution id)
40			Particle2 => Particle2%Next
		END DO

		!! Advance the larger particle and reset the smaller
		Particle1 => Particle1%Next
		Particle2 => Particles%First

	END DO
	!write(*,*) 'CRL: After CoagArray =',CoagArray

!	WRITE(*,*) "Coag Check 1"


	!! Get average solution density
	!! Use this so that mass is conserved in transfers between particles
	Particle1 => Particles%First
	SolnSum =  0
	Count = 0
22	IF(Particle1%Numberofparticles .NE. 0.0) THEN
		SolnSum = SolnSum + Particle1%SolutionDensity
		Count = Count + 1
	END IF
	IF(ASSOCIATED(Particle1%Next)) THEN
		Particle1 => Particle1%Next
		GOTO 22
	END IF

	AvgSolutionDensity = SolnSum/Count
	
	!! Calculate single particle volume
	!!and volume concentrations for each component
	Particle1 => Particles%First
	JJ = 1
	InSumVolumeConc = 0.0
45	IF (Particle1%numberofparticles .GT. 0.0D0) THEN
		SingleParticleVolume(JJ) = (4.0*Pi/3.0)*(Particle1%effectiveradius)**3 !cm3
		
		NumConc(JJ) = Particle1%numberofparticles
		
		InitNumConc(JJ) = Particle1%numberofparticles 

		VolumeConc(JJ) = SingleParticleVolume(JJ)*Particle1%numberofparticles !cm3/cm3


		IF(VolumeConc(JJ) .LT. 0.0) WRITE(*,*) "Negative Volume Conc!"
		
		!Convert aqueous concentrations to cm3 solution/cm3 air
		
		!Water
		AqVolConc(JJ,1)  = Particle1%AqChems(1) * Particle1%NumberOfParticles*AqMolecularMass(1)/AvgSolutionDensity !Particle1%SolutionDensity

		!Solid Salts
		DO I = 2, HowManyAqChems
			AqVolConc(JJ,I)  = Particle1%AqChems(I) * Particle1%NumberOfParticles*AqMolecularMass(I)/SolidSaltDensity(I)
			IF(AqVolConc(JJ,I) .LT. 0.0) WRITE(*,*) "Negative AqVolume Conc!"	
			IF(AqVolConc(JJ,I) .GT. 1.0e30) WRITE(*,*) "Error", AqPhaseChemicalNames(I), AqVolConc(JJ,I)
		END DO

		!Inorganic Ions
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(JJ,I)  = Particle1%AqChems(I) * Particle1%NumberOfParticles*AqMolecularMass(I)/AvgSolutionDensity !Particle1%SolutionDensity
			IF(AqVolConc(JJ,I) .LT. 0.0) WRITE(*,*) "Negative AqVolume Conc!"	
		END DO

		!Convert org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyOrgChems
			OrgVolConc(JJ,I)  = Particle1%OrgChems(I) * Particle1%NumberOfParticles*OrgMolecularMass(I)/OrgDensity(I)
			IF(OrgVolConc(JJ,I) .LT. 0.0) WRITE(*,*) "Negative OrgVolume Conc!"	
		END DO
	
		!Convert hydrophilic org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(JJ,I)  = Particle1%AqOrgChems(I) * Particle1%NumberOfParticles*AqOrgMolecularMass(I)/AqOrgDensity(I)
			IF(AqOrgVolConc(JJ,I) .LT. 0.0) THEN
				WRITE(*,*) "Negative AqOrgVolume Conc!", JJ, I, Particle1%NumberOfParticles, Particle1%AqOrgChems(I)
			END IF
		END DO
					

		TotalAqVolConc = 0.0
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			TotalAqVolConc = TotalAqVolConc + AqVolConc(JJ,I)
		END DO

		TotalOrgVolConc = 0.0
		DO I = 1, HowManyOrgChems
			TotalOrgVolConc = TotalOrgVolConc + OrgVolConc(JJ,I)
		END DO

		TotalAqOrgVolConc = 0.0
		DO I = 1, HowManyAqOrgChems
			TotalAqOrgVolConc = TotalAqOrgVolConc + AqOrgVolConc(JJ,I)
		END DO

		InSumVolumeConc = InSumVolumeConc + (TotalAqVolConc+TotalOrgVolConc+TotalAqOrgVolConc)
	
	ELSE
		
		SingleParticleVolume(JJ) = (4.0*Pi/3.0)*(Particle1%edges(1)**3.0)
		
		NumConc(JJ) = 0.0
		InitNumConc(JJ) = 0.0
		
		VolumeConc(JJ) = 0.0
		
		!Convert aqueous concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(JJ,I)  = 0.0
		END DO
		!Convert organic concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyOrgChems
			OrgVolConc(JJ,I)  = 0.0
		END DO
		!Convert aqueous organic concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(JJ,I)  = 0.0
		END DO
	END IF

	!WRITE(*,*) VolumeConc(JJ), NumConc(JJ)

	!! Go to the next particle
	IF (ASSOCIATED(Particle1%Next)) THEN
		Particle1 => Particle1%Next
		JJ = JJ+1
		GOTO 45
	END IF

!	WRITE(*,*) "Coag Check 2"

	InVolConc = 0.0
	DO JJ = 1, HowManyLinks
		InVolConc = InVolConc + VolumeConc(JJ)
	END DO
	!WRITE(*,*) "In: ", InVolConc
	
	!!Calculate the volume fractions 

	!! For a combined volume BigV, VolumeFractions(I,J,K) can only be non-zero where
	!!   (a) SingleParticleVolume(K) <= BigV < SingleParticleVolume(K+1)   (K < N): the formula below,
	!!       and then K+1 gets 1 minus that value (the K+1 test above cannot also be true there), or
	!!   (b) K = N and BigV >= SingleParticleVolume(N): 1.
	!! Every other element is zero.  This is exactly what the original K loop evaluated (same
	!! comparisons, so NaN volumes behave identically); it finds the brackets in one pass over K
	!! instead of testing three conditions at every K.
	VolumeFractions = 0.0
	DO I = 1, HowManyLinks
		DO J = 1, HowManyLinks
			!Calculate volume of the combined particle
			BigV = SingleParticleVolume(I) + SingleParticleVolume(J)

			SPVK = SingleParticleVolume(1)
			DO K = 1, HowManyLinks-1
				IF (BigV .GE. SPVK .AND. BigV .LT. SingleParticleVolume(K+1)) THEN
					VolumeFractions(I,J,K) = (SingleParticleVolume(K)/BigV)*((SingleParticleVolume(K+1)-BigV)/(SingleParticleVolume(K+1)-SingleParticleVolume(K)))
					IF(VolumeFractions(I,J,K) .GT. 1.0) THEN
						WRITE(*,*) "Volume Fraction out of bounds"
						STOP
					END IF
					VolumeFractions(I,J,K+1) = 1-VolumeFractions(I,J,K)
				END IF
				SPVK = SingleParticleVolume(K+1)
			END DO
			IF (BigV .GE. SingleParticleVolume(HowManyLinks)) VolumeFractions(I,J,HowManyLinks) = 1.0
		END DO
	END DO

!	WRITE(*,*) "Coag Check 3"
	!! The weights VolumeFractions(I,J,K)*CoagArray(I,J) and the conditions under which a (I,J,K)
	!! term contributes do not depend on the chemical species, and neither does LowerSum.
	!! They are collected once here (in the same I, J, K order as the original loops) and reused
	!! for the total volume and for every species, instead of being re-derived for each species.
	ALLOCATE (LowerK(HowManyLinks), KStart(HowManyLinks+1), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate LowerK/KStart at size "//TRIM(INT2STR(HowManyLinks)))

	NEnt = 0
	DO K = 1, HowManyLinks
		KStart(K) = NEnt + 1
		LowerSum = 0
		DO J = 1, HowManyLinks
			IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(K,J) .NE. 0.0 .AND. VolumeFractions(K,J,K) .NE. 1.0) THEN
				LowerSum = LowerSum + (1-VolumeFractions(K,J,K))*CoagArray(K,J)*InitNumConc(J)
			END IF
		END DO
		LowerK(K) = LowerSum
		IF (K .GT. 1) THEN
			DO J = 1, K
				IF (InitNumConc(J) .GT. 0.0D0) THEN
					DO I = 1, K-1
						IF (CoagArray(I,J) .NE. 0.0 .AND. VolumeFractions(I,J,K) .NE. 0.0) NEnt = NEnt + 1
					END DO
				END IF
			END DO
		END IF
	END DO
	KStart(HowManyLinks+1) = NEnt + 1

	ALLOCATE (EntI(MAX(NEnt,1)), EntJ(MAX(NEnt,1)), EntW(MAX(NEnt,1)), EntLast(MAX(NEnt,1)), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate coagulation term lists at size "//TRIM(INT2STR(NEnt)))

	NEnt = 0
	DO K = 2, HowManyLinks
		DO J = 1, K
			IF (InitNumConc(J) .GT. 0.0D0) THEN
				NGroupStart = NEnt + 1
				DO I = 1, K-1
					IF (CoagArray(I,J) .NE. 0.0 .AND. VolumeFractions(I,J,K) .NE. 0.0) THEN
						NEnt = NEnt + 1
						EntI(NEnt) = I
						EntJ(NEnt) = J
						EntW(NEnt) = VolumeFractions(I,J,K)*CoagArray(I,J)
						EntLast(NEnt) = .FALSE.
					END IF
				END DO
				IF (NEnt .GE. NGroupStart) EntLast(NEnt) = .TRUE.
			END IF
		END DO
	END DO

	!Calculate change in total volume conc. and use that to get change in number
	DO K = 1, HowManyLinks
		LowerSum = LowerK(K)

		UpperSum=0
		InnerSum=0
		DO N = KStart(K), KStart(K+1)-1
			InnerSum = InnerSum + EntW(N)*VolumeConc(EntI(N))*InitNumConc(EntJ(N))
			IF (EntLast(N)) THEN
				UpperSum = UpperSum + InnerSum
				InnerSum = 0
			END IF
		END DO

		IF (LowerSum .LT. 0.0 ) CALL ERROR("Problem in LowerSum in Coag routine, early")

		VolumeConc(K) = (VolumeConc(K)+TimeStep*UpperSum)/(1+TimeStep*LowerSum)

		IF (VolumeConc(K) .GT. 0.0 .AND. SingleParticleVolume(K) .GT. 0.0) THEN
			!Number concentration estimated from change in volume concentration
			NumConc(K) = VolumeConc(K)/SingleParticleVolume(K)
		END IF

	END DO

!	WRITE(*,*) "Coag Check 4"

	!Calculate change in volume conc. of individual Aq. species
	DO q = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		DO K = 1, HowManyLinks
			LowerSum = LowerK(K)

			UpperSum=0
			InnerSum=0
			DO N = KStart(K), KStart(K+1)-1
				InnerSum = InnerSum + EntW(N)*AqVolConc(EntI(N),q)*InitNumConc(EntJ(N))
				IF (EntLast(N)) THEN
					IF (InnerSum .LT. 0.0 ) CALL ERROR("Problem in InnerSum in Coag routine")
					UpperSum = UpperSum + InnerSum
					InnerSum = 0
				END IF
			END DO

			IF (LowerSum .LT. 0.0 ) CALL ERROR("Problem in LowerSum in Coag routine")
			IF (UpperSum .LT. 0.0 ) CALL ERROR("Problem in UpperSum in Coag routine")

			AqVolConc(K,q) = (AqVolConc(K,q)+TimeStep*UpperSum)/(1+TimeStep*LowerSum)

			IF(AqVolConc(K,Q) .LT. 0.0) CALL ERROR("Negative inorganic concentrations from Coag Routine.")

		END DO
	END DO

	!Calculate change in volume conc. of individual hydrophobic organic species
	DO q = 1, HowManyOrgChems
		DO K = 1, HowManyLinks
			LowerSum = LowerK(K)

			UpperSum=0
			InnerSum=0
			DO N = KStart(K), KStart(K+1)-1
				InnerSum = InnerSum + EntW(N)*OrgVolConc(EntI(N),q)*InitNumConc(EntJ(N))
				IF (EntLast(N)) THEN
					UpperSum = UpperSum + InnerSum
					InnerSum = 0
				END IF
			END DO

			IF (LowerSum .LT. 0.0 ) CALL ERROR("Problem in LowerSum in Coag routine")
			OrgVolConc(K,q) = (OrgVolConc(K,q)+TimeStep*UpperSum)/(1+TimeStep*LowerSum)
			IF(OrgVolConc(K,Q) .LT. 0.0) CALL ERROR("Negative hydrophobic organic concentrations from Coag Routine.")
		END DO
	END DO

	!Calculate change in volume conc. of individual hydrophilic organic species
	DO q = 1, HowManyAqOrgChems
		DO K = 1, HowManyLinks
			LowerSum = LowerK(K)

			UpperSum=0
			InnerSum=0
			DO N = KStart(K), KStart(K+1)-1
				InnerSum = InnerSum + EntW(N)*AqOrgVolConc(EntI(N),q)*InitNumConc(EntJ(N))
				IF (EntLast(N)) THEN
					UpperSum = UpperSum + InnerSum
					InnerSum = 0
				END IF
			END DO

			IF (LowerSum .LT. 0.0 ) CALL ERROR("Problem in LowerSum in Coag routine")
			AqOrgVolConc(K,q) = (AqOrgVolConc(K,q)+TimeStep*UpperSum)/(1+TimeStep*LowerSum)
			IF(AqOrgVolConc(K,Q) .LT. 0.0) CALL ERROR("Negative aqueous organic concentrations from Coag Routine.")

		END DO
	END DO

	DEALLOCATE(LowerK, KStart, EntI, EntJ, EntW, EntLast)

	OutSumVolumeConc = 0.0
	DO JJ = 1, HowManyLinks
		!Volume Conc Check
		TotalAqVolConc = 0.0
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			TotalAqVolConc = TotalAqVolConc + AqVolConc(JJ,I)
			!IF(AqVolConc(JJ,I) .NE. 0.0) WRITE(*,*) JJ, AqPhaseChemicalNames(I), AqVolConc(JJ,I) 
		END DO

		TotalOrgVolConc = 0.0
		DO I = 1, HowManyOrgChems
			TotalOrgVolConc = TotalOrgVolConc + OrgVolConc(JJ,I)
		END DO

		TotalAqOrgVolConc = 0.0
		DO I = 1, HowManyAqOrgChems
			TotalAqOrgVolConc = TotalAqOrgVolConc + AqOrgVolConc(JJ,I)
		END DO

		SumVolConc(JJ) = TotalAqVolConc+TotalOrgVolConc+TotalAqOrgVolConc
		!WRITE(*,*) JJ, SumVolConc(JJ)
		OutSumVolumeConc = OutSumVolumeConc + (TotalAqVolConc+TotalOrgVolConc+TotalAqOrgVolConc)

	END DO

	OutVolConc = 0.0
	DO JJ = 1, HowManyLinks
		OutVolConc = OutVolConc + VolumeConc(JJ)
	END DO
	!WRITE(*,*) "Volume Ratio", OutVolConc/InVolConc
	IF(ABS((OutVolConc/InVolConc)-1.0) .GT. 0.0001) CALL ERROR("Error in StepSectionalCoagulationJacobson(). Total Volume concentration" &
			           //" is lost during step.")
	
	!WRITE(*,*) "Sum Volume Ratio", OutSumVolumeConc/InSumVolumeConc
	IF(ABS((OutSumVolumeConc/InSumVolumeConc)-1.0) .GT. 0.0001) CALL ERROR("Error in StepSectionalCoagulationJacobson(). Summed Total Volume concentration" &
			           //" is lost during step.")


!	WRITE(*,*) "Coag Check 5"

	!Save adjustments back to particle list
	!!Calculate single particle volume
	!!and volume concentrations for each component
	Particle1 => Particles%First
	JJ = 1
55	IF (NumConc(JJ) .GT. 0.0D0) THEN
		!WRITE(*,*) 'In Gt 0 Coauglation.f90 Particle1%AqChem(1) = ',Particle1%AqChems(1),NumConc(JJ),Particle1%NumberOfParticles

		Particle1%numberofparticles = NumConc(JJ)
		!Convert aqueous concentrations back to mol/particle
		IF(Particle1%SolutionDensity .EQ. 0.0) THEN
			Particle1%SolutionDensity = StoreDensity
		END IF
		StoreDensity = Particle1%SolutionDensity
		Particle1%AqChems(1) = AqVolConc(JJ,1)*AvgSolutionDensity/AqMolecularMass(1)/ Particle1%NumberOfParticles	
		!write(*,*) 'AqVolConc(JJ,1) = ',AqVolConc(JJ,1)
		!write(*,*) 'AvgSolutionDensity = ',AvgSolutionDensity 
		!write(*,*) 'AqMolecularMass(1) = ',AqMolecularMass(1) 
		!write(*,*) 'Particle1%NumberOfParticles = ',Particle1%NumberOfParticles 






		DO I = 2, HowManyAqChems
			Particle1%AqChems(I) = AqVolConc(JJ,I)*SolidSaltDensity(I)/AqMolecularMass(I)/ Particle1%NumberOfParticles		
		END DO
		
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			Particle1%AqChems(I) = AqVolConc(JJ,I)*AvgSolutionDensity/AqMolecularMass(I)/ Particle1%NumberOfParticles			
		END DO
		
		!Convert organic concentrations back to mol/particle
		DO I = 1, HowManyOrgChems
			Particle1%OrgChems(I) = OrgVolConc(JJ,I)*OrgDensity(I)/OrgMolecularMass(I)/ Particle1%NumberOfParticles		
		END DO

		!Convert hydrophilic organic concentrations back to mol/particle
		DO I = 1, HowManyAqOrgChems
			Particle1%AqOrgChems(I) = AqOrgVolConc(JJ,I)*AqOrgDensity(I)/AqOrgMolecularMass(I)/ Particle1%NumberOfParticles		
		END DO

		!! This will do RADIUS, EffectiveRadius, & DENSITY
		!write(*,*) 'Before RecalculateRadius(Particle1) Particle1%AqChems(1) = ',Particle1%AqChems(1),NumConc(JJ)
		CALL RecalculateRadius(Particle1) !Calls ParticleDensity!
		!write(*,*) 'After RecalculateRadius(Particle1)Particle1%AqChems(1) = ',Particle1%AqChems(1),NumConc(JJ)

	ELSE
		!WRITE(*,*) 'In Else Coauglation.f90 Particle1%AqChem(1) = ',Particle1%AqChems(1),NumConc(JJ)
		Particle1%numberofparticles = 0.0D0
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			Particle1%AqChems(I) = 0.0
		END DO
		
		DO I = 1, HowManyOrgChems
			Particle1%OrgChems(I) = 0.0
		END DO
		DO I = 1, HowManyAqOrgChems
			Particle1%AqOrgChems(I) = 0.0
		END DO
		!WRITE(*,*) 'In Coauglation.f90 Particle1%AqChem(1) = ',Particle1%AqChems(1),NumConc(JJ)
	END IF



	!! Go to the next particle
	IF (ASSOCIATED(Particle1%Next)) THEN
		Particle1 => Particle1%Next
		JJ = JJ+1
		GOTO 55
	END IF

!	WRITE(*,*) "Coag Check 6"

	!Regrid the aerosol to make sure sizes work out correctly
	!(Also calls aqueous equilibrium, and resets temperature)
	!write(*,*) 'Before RegridAerosol in Coagulation.f90 calling spillbeans...'
	!call spillbeans()
	CALL RegridAerosol(RadiiCurrent = .TRUE.)
	!write(*,*) 'After RegridAerosol in Coagulation.f90 calling spillbeans...'
	!call spillbeans()	
	!! Houseclean memory
	DEALLOCATE(CoagArray)
	DEALLOCATE(PreR, PreD, PreMass, PreDel, PreRe, PreVt, PreSect)
	DEALLOCATE(SingleParticleVolume)
	DEALLOCATE(VolumeFractions)
	DEALLOCATE(VolumeConc)
	DEALLOCATE(NumConc)
	DEALLOCATE(AqVolConc)
	DEALLOCATE(OrgVolConc)
	DEALLOCATE(AqOrgVolConc)
			
!	WRITE(*,*) "Coag Check 7"
	RETURN


   	!Particle1 => particles%first
	!DO WHILE(associated(Particle1))
	!	write(*,*) 'CRL: In StepSectionalCoagulation - current%AqChems = ',Particle1%AqChems
	!	write(*,*) 'CRL: In StepSectionalCoagulation - current%numberofparticles = ',Particle1%numberofparticles
	!	Particle1 => Particle1%next
	!END DO !END STEP PARTICLE DILUTION



END SUBROUTINE StepSectionalCoagulationJacobson



END MODULE Coagulation
