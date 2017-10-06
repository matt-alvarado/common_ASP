!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! StepASP.f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						             !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 06/27  2007   Matt Alvarado     Added InputConcentrations
!!		                   Added OutputConcentrations
!!                                 Added InitializeASP
!! 07/16  2007   Matt Alvarado     Set any aerosol concentration below
!!					1.0e-30 mol/particle to 0.
!! 07/18  2007   Matt Alvarado     Corrected output concentrations for
!!					solid electrolytes
!! 08/31  2007   Matt Alvarado     Added aerosol size to interface
!! 10/09  2007   Matt Alvarado     Changed zero criteria for 
!!                                     aerosol concentrations in 
!!				       InputConcentrations to 
!!                                    1.0e-40 mol/particle
!! 10/12  2007   Matt Alvarado     Removed Scale from photorates: 
!!                                     now in chemtrop.F
!! 02/19  2009   Matt Alvarado     Changed DFPORT to IFPORT 
!!                                     (to use ifort compiler)
!! 08/26  2010   Matt Alvarado     Comment out IFPORT (to use pgf90 compiler)
!! 02/15  2012   Matt Alvarado     Removed Eulerian Coords
!! 05/03  2012   Matt Alvarado     Updated call to AerosolOptProps to include
!!                                  BackScatCoeff
!! 08/16  2012   Matt Alvarado     Updated call to AerosolOptProps to include
!!                                  SubExtCoeff and SubSSA
!! 08/17  2012   Matt Alvarado     Updated call to AerosolOptProps to include
!!                                  SubExtCoeff and SubSSA
!! 11/08  2012   Matt Alvarado     Updated call to AerosolOptProps to include
!!                                  mixing flag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES	                            !!
!! 1. ModelParameters			    !!
!! 2. GridPointFields			    !!
!! 3. Aerosols				    !!
!! 4. Condensation			    !!
!! 5. Time				    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	    !!	
!! 1. SUBROUTINE StepASPOnce ()					    !!
!! 2. SUBROUTINE InputConcentrations
!! 3. SUBROUTINE OutputConcentrations
!! 4. SUBROUTINE InitializeASP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE StepASP

		IMPLICIT NONE

		PRIVATE

		PUBLIC ::	ASPInterface, InitializeASP

	CONTAINS


	SUBROUTINE StepASPOnce (TimeStep)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine recalculates the gas-phase reactions rates,  !!
	!! then steps gas chemistry, condensation, coagulation, and     !!
	!! updates the optical properties of the aerosol.	        !!
	!! All calculations are done in gridbox 1,1,1 MJA 052307	!!
        !! Removed Eulerian Grid Points MJA 021512                      !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


		!USE IFPORT,				 ONLY : RTC, ETIME
		USE ModelParameters
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: TimeStep
		
		!! Internal Variables
		INTEGER :: I, NumBins
		REAL*8 :: TEMP
	
		TYPE(Particle), POINTER :: cur
		
		!Calculate Temperature	
		TEMP = GetTemp()
					
		!Step Gas Chemistry
		CALL StepGasChemistry(TimeStep,1)
			
		!Step Condensation, if desired
		cur => particles%first
		IF(.TRUE.) THEN
			CALL StepCondensationAll (TimeStep, 1)
			!Note above calls RegridAerosol automatically
			!with recalculates radius and density
			!WRITE(*,*) "Cond Okay"
		END IF

		!Step Coagulation, if desired
		IF (.TRUE.) THEN	
			CALL StepSectionalCoagulationJacobson(TimeStep)
			!Note that above already calls RegridAerosol
			!which recalculates radius and density
			!WRITE(*,*) "Coag Okay"
		END IF

		!Step Optical Properties Particles
		CALL SortAerosolAtGridPointForCoagulation () 
		cur => particles%first
		I = 1
		DO WHILE(associated(cur))
    			
			!Force particle temperature to be env. Temperature
			!WARNING: This is a kludge, since the program is 
                        !calculating a huge temperature for the largest 
                        !particles, and I'm not sure why
			cur%Temperature = TEMP
					
			!Recalculate optical parameters
			!WRITE(*,*) "Call Optical, Particle #", I
			CALL ShellRefIndAndRad(cur)
			!WRITE(*,*) "Optical Okay"
				
			I = I + 1
		
			cur => cur%next
				
		END DO !END STEP OPTICAL PROPERTIES
	
		RETURN
	END SUBROUTINE StepASPOnce

	SUBROUTINE InputConcentrations(Temp, Press, Dens, GasConc, &
				AeroNumConc, AeroMassConc, inPhotoRates)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine takes in Temperature, pressure,		!!
	!! gas and aerosol concentrations, and photolysis rates from	!!
	!! a 3D model and updates the appropriate concentrations in ASP !!
	!! grid box 1,1,1. MJA 052307					!!
        !! Removed Eulerian coordinates, MJA 021512                     !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
	  !USE IFPORT,				ONLY : RTC, ETIME
	  USE ModelParameters,	ONLY : HowManyBins
	  USE InfrastructuralCode	
	  USE GridPointFields
	  USE Chemistry
	  USE Aerosols		
	  USE Condensation
	  USE OutputRoutines
	  USE Coagulation

	  IMPLICIT NONE

	  !! External Variables
	  REAL*8, INTENT(IN) :: Temp  !(in K)
	  REAL*8, INTENT(IN) :: Press !(in mbar)
	  REAL*8, INTENT(IN) :: Dens  !(in kg/m3 dry air)
	  REAL*8 :: GasConc(HowManyEvolveGasChems)
	  REAL*8, INTENT(IN) :: AeroNumConc(HowManyBins)
	  REAL*8 :: AeroMassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
	  REAL*8, INTENT(IN) :: inPhotoRates(12)

	  !! Internal Variables
	  INTEGER :: I, J, L, RxnIndex, NumAq, GG
	  REAL*8 :: HeteroNumberConc(HowManyBins)
          REAL*8 :: HeteroRadii(HowManyBins)
	  REAL*8 :: Scale

	  TYPE(Particle), POINTER :: cur
		
	  !Set temperature
	  CALL SetTempField (Temp)

	  !Set Pressure
	  CALL SetPressField (Press)
		
	  !Update chemical concentrations
	  DO I = 1, HowManyEvolveGasChems
             IF(GasConc(I) .LT. 0.) GasConc(I) = 0.
	  END DO 
	  CALL UpdateChemicalConcentrations(GasConc(1:HowManyEvolveGasChems))
		
	  !Update aerosol concentrations
	  CALL SortAerosolAtGridPointForCoagulation () 
	  cur => particles%first !current particle
	  I = 1
		
	  NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
	  DO WHILE(associated(cur))
            cur%numberofparticles = AeroNumConc(I)
	    IF(cur%Numberofparticles .GT. 0.) THEN
	       IF(cur%dry) THEN
		  cur%dry = .FALSE.
	       END IF

	       DO J = 1, NumAq
	          IF(AeroMassConc(I,J) .LT. 0.) AeroMassConc(I,J) = 0.
	          cur%AqChems(J) = AeroMassConc(I,J)*Dens*1.0e-12 &
	                  /(cur%Numberofparticles*AqMolecularMass(J))
		  IF(cur%AqChems(J).LT.1.0e-40) cur%AqChems(J)=0.
	       END DO

	       DO J = 1, HowManyOrgChems
		  IF(AeroMassConc(I,J+NumAq) .LT. 0.) AeroMassConc(I,J+NumAq) = 0.
		  Cur%OrgChems(J) = AeroMassConc(I,J+NumAq)*Dens*1.0e-12 &
                        /(cur%Numberofparticles*OrgMolecularMass(J))
		  IF(Cur%OrgChems(J) .LT. 1.0e-40) Cur%OrgChems(J) = 0.
	       END DO
			
	       DO J = 1, HowManyAqOrgChems
		  IF(AeroMassConc(I,J+NumAq+HowManyOrgChems) .LT. 0.) &
                     AeroMassConc(I,J+NumAq+HowManyOrgChems) = 0.
		  Cur%AqOrgChems(J) = AeroMassConc(I,J+NumAq+HowManyOrgChems) &
                     *Dens*1.0e-12 &
		     /(cur%Numberofparticles*AqOrgMolecularMass(J))
		  IF(Cur%AqOrgChems(J) .LT. 1.0e-40) Cur%AqOrgChems(J) = 0.
	       END DO
	    ELSE !No particles
	       cur%numberofparticles = 0.
	       DO J = 1, NumAq
		 Cur%AqChems(J) = 0.
	       END DO
				
	       DO J = 1, HowManyOrgChems
		 Cur%OrgChems(J) = 0.
	       END DO
				
	       DO J = 1, HowManyAqOrgChems
		 Cur%AqOrgChems(J) = 0.
	       END DO				
	    END IF
		
       	    !Force particle temperature to be env. Temperature
	    cur%Temperature = Temp
				
	    !Equilibrate, recalculate radius and density
	    !WRITE(*,*) "Before Recauculate Radius", cur%effectiveradius
	    CALL FindElectrolyteEquilibrium (cur, UpdateThermo=.TRUE., FirstEquilibration=.TRUE., ReturnType=GG)
	    CALL RecalculateRadius(cur)		
	    !WRITE(*,*) "After Recauculate Radius", cur%effectiveradius
	    I = I + 1
	    cur => cur%next		
          END DO
	  !WRITE(*,*) "After Aero Chems Input"
		
	  !Recalculate all reaction rates
	  CALL RecalculateReactionRates(IFGAS,IFAQ)
		
	  !WRITE(*,*) "Before Photorates"
	  !Reset Photolysis Rates
	  DO L = 1, GasPhaseReactionRateArraySize(1)
						
	    !! Only do photolysis reactions
	    IF (INT(GasPhaseReactionRates(L,1)) .EQ. 1 .AND. &
                INT(GasPhaseReactionRates(L,3)) .EQ. 0.) THEN	
								
		!Find reaction index number
		RxnIndex = INT(GasPhaseReactionRates(L,4))

		SELECT CASE (RxnIndex)
		  CASE(4) !NO2 => NO + O
		    GasPhaseChemicalRates(L) = inPhotoRates(1)
		  CASE(5) !NO3 => NO + O2
		    GasPhaseChemicalRates(L) = inPhotoRates(2)
		  CASE(6) !NO3 => NO2 + O
		    GasPhaseChemicalRates(L) = inPhotoRates(3)		
		  CASE(3) !O3 => O + O2
	            GasPhaseChemicalRates(L) = inPhotoRates(4)
		  CASE(2) !O3 => O1D + O2
		    GasPhaseChemicalRates(L) = inPhotoRates(5)
		  CASE(12) !HONO => Products
		    GasPhaseChemicalRates(L) = inPhotoRates(6)
		  CASE(11) !H2O2 => 2OH
		    GasPhaseChemicalRates(L) = inPhotoRates(7)
		  CASE(15) !HCHO => 2HO2 + CO
		    GasPhaseChemicalRates(L) = inPhotoRates(8)
		  CASE(16) !HCHO => H2 + CO
		    GasPhaseChemicalRates(L) = inPhotoRates(9)
                  CASE(17) !ALD2 => Products
		    GasPhaseChemicalRates(L) = inPhotoRates(10)
		  CASE(23) !KETL => Products
		    GasPhaseChemicalRates(L) = inPhotoRates(11)
	          CASE(22) !MGLY => Products
		    GasPhaseChemicalRates(L) = inPhotoRates(12)
		END SELECT
	     END IF
          END DO
	  !WRITE(*,*) "After Photorates"

	  !Get particle number conc. and radii for hetero. rate calculation
	  cur => particles%first
	  I = 1
	  DO WHILE(associated(cur))
	    HeteroNumberConc(I) = cur%NumberofParticles
            HeteroRadii(I) = cur%EffectiveRadius/100.0 !convert from cm to m
	    !ResetHeteroRates expects radii in units of m
	    cur => cur%next
	    I = I+1
	  END DO
				
	  !Calculate Heterogeneous rates from size dist. info
	  CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HowManyBins, .FALSE.)
	  !WRITE(*,*) "Hetero Okay"

	END SUBROUTINE InputConcentrations

	SUBROUTINE OutputConcentrations(Temp, Press, Dens, GasConc, &
					AeroNumConc, AeroMassConc, &
					ExtCoeff, SingScat, Assym, Radius)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine outputs Temperature, pressure,		!!
	!! gas and aerosol concentrations, and optical properties from	!!
	!! a 3D model and updates the appropriate concentrations in ASP !!
	!! grid box 1,1,1. MJA 052307					!!
        !! Removed Eulerian coordinates, MJA 021512                     !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
		!USE IFPORT,				ONLY : RTC, ETIME
		USE ModelParameters,	ONLY : HowManyBins
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: Temp, Press, Dens
		REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym,BackScatCoeff,SubExtCoeff, SubSSA, &
                                          ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                                          SubExtCoeffSep, SubSSASep
		REAL*8 :: GasConc(HowManyEvolveGasChems)
		REAL*8 :: AeroNumConc(HowManyBins)
		REAL*8 :: AeroMassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
		REAL*8 :: Radius(HowManyBins)
		
		!! Internal Variables
		INTEGER :: I, J, NumAq
		REAL*8 :: SaltConc(HowManyAqChems)
	
		TYPE(Particle), POINTER :: Cur
		
		
		!Output Temperature and Pressure
		Temp = GetTemp()
		Press = GetPress()
	
		!Output Gas chemical concentrations
		DO I = 1, HowManyEvolveGasChems
			GasConc(I) = GridGasChem(I)
		END DO

		!Update aerosol concentrations
		CALL SortAerosolAtGridPointForCoagulation () 
		cur => particles%first
		I = 1
		NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
		DO WHILE(associated(cur))
    		AeroNumConc(I) = cur%numberofparticles
			Radius(I) = cur%effectiveradius
			
			!Water
			AeroMassConc(I,1) = cur%AqChems(1)*(cur%Numberofparticles*AqMolecularMass(1)) &
					                    /(Dens*1.0e-12)
												
			!Do this step so that we only pass neutral electrolytes
			!to dynamics code - helps prevent odd charge balance problems
			!WRITE(*,*) "Before output HypotheticalElectrolyteConcentrations"
			SaltConc = HypotheticalElectrolyteConcentrations (cur)
			!WRITE(*,*) SaltConc(5)+SaltConc(8)+SaltConc(12), cur%AqChems(5)+cur%AqChems(8)+cur%Aqchems(12)
		
			IF(cur%Numberofparticles .GT. 0.) THEN
				DO J = 2, HowManyAqChems
					AeroMassConc(I,J) = (SaltConc(J)+Cur%AqChems(J))*(cur%Numberofparticles*AqMolecularMass(J)) &
					                    /(Dens*1.0e-12)
				END DO
				
				DO J = HowManyAqChems+1, NumAq
					AeroMassConc(I,J) = 0.
				END DO

				DO J = 1, HowManyOrgChems
						AeroMassConc(I,J+NumAq) = Cur%OrgChems(J)*(cur%Numberofparticles*OrgMolecularMass(J)) &
					                    /(Dens*1.0e-12)
				END DO
				
				DO J = 1, HowManyAqOrgChems
						AeroMassConc(I,J+NumAq+HowManyOrgChems) = Cur%AqOrgChems(J)*(cur%Numberofparticles*AqOrgMolecularMass(J)) &
					                    /(Dens*1.0e-12)
				END DO
			ELSE !No particles
				AeroNumConc(I) = 0.
				DO J = 1, HowManyAqChems
					AeroMassConc(I,J) = 0.
				END DO
				
				DO J = 1, HowManyOrgChems
					AeroMassConc(I,J+HowManyAqChems) = 0.
				END DO
				
				DO J = 1, HowManyAqOrgChems
					AeroMassConc(I,J+HowManyAqChems+HowManyOrgChems) = 0.
				END DO				
			END IF
		
			I = I + 1
			cur => cur%next
				
		END DO

		!Output Average Aerosol optical properties, 0 flag assumes core in shell mixing
		CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSSA, 0)
		!WRITE(*,*) ExtCoeff(1), SingScat(1), Assym(1)
		!PAUSE

	END SUBROUTINE OutputConcentrations

	SUBROUTINE InitializeASP()

	!! Include modules
	!USE IFPORT,				 ONLY : RTC, ETIME
	USE ModelParameters
	USE InfrastructuralCode	
	USE GridPointFields
	USE Chemistry
	USE Aerosols		
	USE Condensation
	USE OutputRoutines
	USE Coagulation
	
	implicit none

	!! Define the local variables
	integer :: i, j, k, l,q, r, NumBins
        INTEGER :: time(8)
	REAL*4  :: etimearray(2)
	TYPE(Particle), POINTER :: Cur, Env
	CHARACTER (len = 8)     :: ErrorDate
	CHARACTER (len = 10)    :: ErrorTime

	!! Initialize the random number generator 
        CALL Random_Seed()

		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize ()
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()
	
	!! Send a Welcome Message
	CALL Transcript ("**********************************************************")
	CALL Transcript ("** Welcome to AER's ASP Model developed by Matt Alvarado  **")
	CALL Transcript ("** (mjalvara@mit.edu).                                  **")
	CALL Transcript ("**********************************************************")
	CALL Date_and_Time(ErrorDate, ErrorTime)
	CALL Transcript ("*********************************************")
	CALL Transcript ("** Run Started at                          **")
	CALL Transcript ("** "//ErrorTime(1:2)//":"//ErrorTime(3:4)//":"//ErrorTime(5:6)//" on "//ErrorDate(5:6)//		&
					  "/"//ErrorDate(7:8)//"/"//ErrorDate(1:4)//" "//"                 "//"**")
	CALL Transcript ("*********************************************")
				
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams

		!! INITIALIZE the DISSOLUTION routine
		CALL SetAllDissolution

		!! INITIALIZE the PARTICLES	
		CALL InitializeParticlesSectional
				
		!!WARNING Force initial water equilibrium between gas and aerosol
		!!(Uses open system for water eq.)
		!CALL EquilibrateInternallyAtGridPoint(1,1,1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
		CALL RegridAerosol ()

	END SUBROUTINE InitializeASP

	SUBROUTINE ASPInterface(Timestep, Temp, Press, Dens, GasConc, &
				NumConc, MassConc, outPhotoRates, &
				ExtCoeff, SingScat, Assym, Radius)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine connects ASP to an external 3D Eulerian model!! 
	!! MJA 062707							!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
		!USE IFPORT,				ONLY : RTC, ETIME
		USE	ModelParameters,	ONLY : HowManyBins
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: Timestep !(in s)
		REAL*8 :: Temp  !(in K)
		REAL*8 :: Press !(in mbar)
		REAL*8 :: Dens  !(in kg/m3 dry air)
		REAL*8 :: GasConc(HowManyEvolveGasChems)
		REAL*8 :: NumConc(HowManyBins)
		REAL*8 :: MassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
		REAL*8 :: Radius(HowManyBins)
		REAL*8 :: outPhotoRates(12)
		REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym
		INTEGER :: I
	
		!WRITE(*,*) "Before Conc Input"
		CALL InputConcentrations(Temp, Press, Dens, GasConc, &
									NumConc, MassConc, outPhotoRates)
		!WRITE(*,*) "After Conc Input"
		CALL StepASPOnce (TimeStep)
		!WRITE(*,*) "Before Conc Output"
		CALL OutputConcentrations(Temp, Press, Dens, GasConc, &
									NumConc, MassConc, &
									ExtCoeff, SingScat, Assym, Radius)
		!WRITE(*,*) "After Conc Output"

	END SUBROUTINE ASPInterface
END MODULE StepASP
