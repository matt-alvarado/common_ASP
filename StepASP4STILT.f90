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
!! 05/18  2016   Matt Alvarado     Made changes to allow aerosols to interface 
!!                                  with STILT-Chem
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
!! 5. SUBROUTINE ASPInterface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE StepASP

    IMPLICIT NONE

            ! CMB: aerosol flag is private
    PRIVATE :: isAerosolProc !, parcelIndex, timeIndex

            ! CMB: add last two methods
    PUBLIC :: ASPInterface, InitializeASP, setAerosolProc, getAerosolProc
    integer :: isAerosolProc !, parcelIndex, timeIndex
CONTAINS

    ! CMB add
    integer function getAerosolProc() result(aeroProc)
        aeroProc = isAerosolProc
    end function
    subroutine setAerosolProc(aeroProc)
        integer, intent(in) :: aeroProc
        isAerosolProc = aeroProc
    end subroutine

   
        
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
		
        ! CMB add
        integer :: doAerosols != getAerosolProc()
        doAerosols = getAerosolProc()
        ! end CMB add
		
		
        !Calculate Temperature	
        TEMP = GetTemp()
					
        !Step Gas Chemistry
                !print *, 'Stepping gas chemistry'
                !print *, 'Time step = ', timestep
        !print *, 'GasConc(599) at start of InputConcentrations = ', GasConc(599)
        !print *, 'GridGasChem(599) at start of StepAspOnce = ', GridGasChem(599)  
        !call spillbeans(tout_unit=100)     
        CALL StepGasChemistry(TimeStep,1)
        !call spillbeans(tout_unit=101)
        !print *, 'GridGasChem(599) after StepGasChemistry = ', GridGasChem(599)
		
        if (doAerosols .gt. 0) then
            !Step Condensation, if desired
            cur => particles%first
                    
            IF(.TRUE.) THEN
                !print *, 'calling stepcondensationall from StepASP()'
                !print *, 'GridGasChem(599) before StepCondensationAll = ', GridGasChem(599) 
                !call spillbeans(tout_unit=103)
                CALL StepCondensationAll (TimeStep, 1)
                !call spillbeans(tout_unit=104)
                !print *, 'GridGasChem(599) after StepCondensationAll = ', GridGasChem(599) 
                !Note above calls RegridAerosol automatically
                !with recalculates radius and density
                !WRITE(*,*) "Cond Okay"
            END IF
    
            !Step Coagulation, if desired
            !if (.false.) then   ! cmb change
            IF (.TRUE.) THEN
                 !           print *, 'calling stepsectionalcoagulationjacobson'
                !print *, 'GridGasChem(599) before StepSectionalCoagulationJacobson = ', GridGasChem(599)             
                !call spillbeans(tout_unit=105)
                CALL StepSectionalCoagulationJacobson(TimeStep)
                !call spillbeans(tout_unit=106)
                !print *, 'GridGasChem(599) before StepSectionalCoagulationJacobson = ', GridGasChem(599) 
                !Note that above already calls RegridAerosol
                !which recalculates radius and density
                !WRITE(*,*) "Coag Okay"
            END IF
    
            !Step Optical Properties Particles
                    !print *, 'calling sortaerosolatgridpointforcoagulation from stepasp'
            !print *, 'GridGasChem(599) before SortAerosolAtGridPointForCoagulation = ', GridGasChem(599)        
            !call spillbeans(tout_unit=107)
            CALL SortAerosolAtGridPointForCoagulation ()
            !call spillbeans(tout_unit=108)
            !print *, 'done with sortaerosolatgridpointforcoagulation from stepasp'
            !print *, 'GridGasChem(599) before SortAerosolAtGridPointForCoagulation = ', GridGasChem(599) 
            cur => particles%first
            I = 1
            DO WHILE(associated(cur))
                    
                !Force particle temperature to be env. Temperature
                !WARNING: This is a kludge, since the program is 
                            !calculating a huge temperature for the largest 
                            !particles, and I'm not sure why
                            !print *, 'applying a temperature kludge'
                cur%Temperature = TEMP
                        
                !Recalculate optical parameters
                !WRITE(*,*) "Call Optical, Particle #", I
                ! CMB (AER, Inc): This call significantly slows down
                !                 aerosol processing.  Comment out until
                !                 we figure everything out.
                !CALL ShellRefIndAndRad(cur)
                !WRITE(*,*) "Optical Okay"
                    
                I = I + 1
            
                cur => cur%next
                    
            END DO !END STEP OPTICAL PROPERTIES
        endif
        !print *, 'GridGasChem(599) at end of StepAspOnce = ', GridGasChem(599) 
        RETURN
    END SUBROUTINE StepASPOnce

    SUBROUTINE InputConcentrations(Temp, Press, Dens, GasConc, &
        AeroNumConc, AeroMassConc, &
        !PhotoRates, &
        sza, Transmissivity)
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
        REAL*8 :: GasConc(HowManyEvolveGasChems) !(in molecules/cm3)
        REAL*8, INTENT(IN) :: AeroNumConc(HowManyBins) !(in particles/cm3)
        REAL*8 :: AeroMassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems) !(in ug/m3)
        REAL*8, INTENT(IN) :: Transmissivity
        !REAL*8, INTENT(IN) :: PhotoRates(12)

        !! Internal Variables
        INTEGER :: I, J, L, RxnIndex, NumAq, GG, iii
        REAL*8 :: HeteroNumberConc(HowManyBins)
        REAL*8 :: HeteroRadii(HowManyBins)
        REAL*8 :: Scale

        TYPE(Particle), POINTER :: cur
	  
        ! CMB add
        real*8, intent(in) :: sza
        integer :: doAerosols
        doAerosols = getAerosolProc()
        ! end CMB add
	  	
        !Set temperature
                  !print *, 'setting temperature = ', temp
        CALL SetTempField (Temp)

        !Set Pressure
                  !print *, 'setting pressure = ', press
        CALL SetPressField (Press)
		
        !Update chemical concentrations
                  !print *, 'updating gas concentrations'
        DO I = 1, HowManyEvolveGasChems
            IF(GasConc(I) .LT. 0.) GasConc(I) = 0.
                     !if (gasconc(i) .lt. 1.0e-29) gasconc(i) = 1.0e-30
        END DO 
        !print *, 'GasConc(599) at start of InputConcentrations = ', GasConc(599)
        !print *, 'GridGasChem(599) at start of InputConcentrations = ', GridGasChem(599)
        CALL UpdateChemicalConcentrations(GasConc(1:HowManyEvolveGasChems))
        !print *, 'GasConc after UpdateChemicalConcentrations(599) = ', GasConc(599)
        !print *, 'GridGasChem(599) after UpdateChemicalConcentrations = ', GridGasChem(599)
	  
        if (doAerosols .gt. 0) then
            !call spillbeans(tout_unit=356)
	        	
            !Update aerosol concentrations
                !print *, 'calling sort aerosol at grid point for coag'
            !print *, 'GasConc(599) at start of SortAerosolAtGridPointForCoagulation = ', GasConc(599)
            !print *, 'GridGasChem(599) at start of SortAerosolAtGridPointForCoagulation = ', GridGasChem(599)    
            CALL SortAerosolAtGridPointForCoagulation () 
            !print *, 'GasConc after SortAerosolAtGridPointForCoagulation = ', GasConc(599)
            !print *, 'GridGasChem(599) after SortAerosolAtGridPointForCoagulation = ', GridGasChem(599)
            cur => particles%first !current particle
            I = 1
            
            NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
            DO WHILE(associated(cur))
                cur%numberofparticles = AeroNumConc(I)
    
                                ! CB: for debugging
                                !print *, 'Particle ID = ', cur%particleID
                                !print *, 'I = ', I
                                !print *, 'AeroNumConc(I) = ', AeroNumConc(I)
                                !write(*,*) "cur%numberofparticles = ", cur%numberofparticles
                                !write(*,*) 'Initial NH3 (aq) AeroMassConc(1, 3) = ', AeroMassConc(1,3)
                !Force minimum water
                IF(AeroMassConc(I,1) .LT. 0.001) 	AeroMassConc(I,1) = 0.001						   
                    !write(*,*) '    Corrected AeroMassConc(I, 1) = ', AeroMassConc(I,1)
                    ! end CB
    
                    ! CB: Modify this condition to accept floating-point tolerance
                    !     equality situation
                !IF(cur%Numberofparticles .GT. 0.) THEN
                if(cur%numberofparticles .gt. 1.0E-40) then
                         !print *, '    ******** number of particles > 0 ********'
                    IF(cur%dry) THEN
                        cur%dry = .FALSE.
                    END IF
    
                    DO J = 1, NumAq
                        IF(AeroMassConc(I,J) .LT. 0.) AeroMassConc(I,J) = 0. ! 1.0e-30
                        cur%AqChems(J) = AeroMassConc(I,J)*1.0e-12 &
                            /(cur%Numberofparticles*AqMolecularMass(J))
                        !write(*,*) '    From Array: j = ', j, ', cur%AqChems(j) = ', cur%aqchems(j)
                        IF(cur%AqChems(J).LT.1.0e-40) cur%AqChems(J)=0.  !cur%aqChems(j) = 1.0e-30
                        !write(*,*) '    After Checking For 0: j = ', j, ', cur%AqChems(j) = ', cur%aqchems(j)
                    END DO
                    call flush(6)
                    DO J = 1, HowManyOrgChems
                        IF(AeroMassConc(I,J+NumAq) .LT. 0.) AeroMassConc(I,J+NumAq) = 0.
                        Cur%OrgChems(J) = AeroMassConc(I,J+NumAq)*1.0e-12 &
                            /(cur%Numberofparticles*OrgMolecularMass(J))
                        !write(*,*) '    From Array: j = ', j, ', cur%orgChems(j) = ', cur%orgchems(j)    
                        IF(Cur%OrgChems(J) .LT. 1.0e-40) Cur%OrgChems(J) = 0.
                        !write(*,*) '    After Checking for 0: j = ', j, ', cur%orgChems(j) = ', cur%orgchems(j)
                    END DO
                    call flush(6)
                    DO J = 1, HowManyAqOrgChems
                        IF(AeroMassConc(I,J+NumAq+HowManyOrgChems) .LT. 0.) &
                            AeroMassConc(I,J+NumAq+HowManyOrgChems) = 0.
                        !)
                        Cur%AqOrgChems(J) = AeroMassConc(I,J+NumAq+HowManyOrgChems) &
                            *1.0e-12 &
                            /(cur%Numberofparticles*AqOrgMolecularMass(J))
                        !write(*,*) '    From Array: j = ', j, ', cur%AqOrgChems(j) = ', cur%aqorgchems(j)    
                        IF(Cur%AqOrgChems(J) .LT. 1.0e-40) Cur%AqOrgChems(J) = 0.
                        !write(*,*) '    After Checking For 0: j = ', j, ', cur%AqOrgChems(j) = ', cur%aqorgchems(j)
                    END DO
                    call flush(6)
                ELSE !No particles
                         !print *, '**** number of particles not > 0 ********'
                    cur%numberofparticles = 0.
                    !cur%numberofparticles = 1.0e-30
                    DO J = 1, NumAq
                        Cur%AqChems(J) = 0.
                    !cur%aqchems(j) = 1.0e-30
                    END DO
                    
                    DO J = 1, HowManyOrgChems
                        Cur%OrgChems(J) = 0.
                    !cur%orgchems(j) = 1.0e-30
                    END DO
                    
                    DO J = 1, HowManyAqOrgChems
                        Cur%AqOrgChems(J) = 0.
                    !cur%aqorgchems(j) = 1.0e-30
                    END DO				
                END IF
            
                    !Force particle temperature to be env. Temperature
                cur%Temperature = Temp
                    
                !Equilibrate, recalculate radius and density
                ! CMB (AER, Inc): Add if condition for this block to not do a find
                !                 here if no water.
                ! Nope, causes some weird errors?
                !if (cur%aqchems(j) .gt. 1.0e-34) then
                !WRITE(*,*) "Before FindElectrolyteEquilibrium in InputConcentrations"
                !call flush(6)
                !print *, 'GasConc(599) at start of FindElectrolyteEquilibrium = ', GasConc(599)
                !print *, 'GridGasChem(599) at start of FindElectrolyteEquilibrium = ', GridGasChem(599) 
                !do iii = 1, HowManyBins
                !    print *, 'IVOC1 (Org)[', iii, '] = ', AeroMassConc(iii, 56)
                !    print *, 'IVOC1 (AqOrg)[', iii, '] = ', AeroMassConc(iii, 78)
                !    print *, '    NumConc[', iii, '] = ', AeroNumConc(iii)
                !enddo   
                CALL FindElectrolyteEquilibrium (cur, UpdateThermo=.TRUE., &
                    FirstEquilibration=.TRUE., &
                    ReturnType=GG)
                                                 
                !print *, 'GasConc(599) at end of FindElectrolyteEquilibrium = ', GasConc(599)
                !print *, 'GridGasChem(599) at end of FindElectrolyteEquilibrium = ', GridGasChem(599) 
                !do iii = 1, HowManyBins
                !    print *, 'IVOC1 (Org)[', iii, '] = ', AeroMassConc(iii, 56)
                !    print *, 'IVOC1 (AqOrg)[', iii, '] = ', AeroMassConc(iii, 78)
                !    print *, '    NumConc[', iii, '] = ', AeroNumConc(iii)
                !enddo                                        
                !write(*,*) 'Done with FindElectrolyteEquilibrium in InputConcentrations'
                !call flush(6)
                !endif
                CALL RecalculateRadius(cur)		
!                print *, 'GasConc(599) at end of RecalculateRadius = ', GasConc(599)
!                print *, 'GridGasChem(599) at end of RecalculateRadius = ', GridGasChem(599)   
!                do iii = 1, HowManyBins
!                    print *, 'IVOC1 (Org)[', iii, '] = ', AeroMassConc(iii, 56)
!                    print *, 'IVOC1 (AqOrg)[', iii, '] = ', AeroMassConc(iii, 78)
!                    print *, '    NumConc[', iii, '] = ', AeroNumConc(iii)
!                enddo                                      
                !WRITE(*,*) "After Recauculate Radius", cur%effectiveradius
                I = I + 1
                cur => cur%next		
            END DO
            !call spillbeans(tout_unit=358)
                  !WRITE(*,*) "After Aero Chems Input"
        else
            !write(*,*) 'Not doing aerosols'
        endif
	  	!call spillbeans(tout_unit=765)
	  	!stop
	  	
        !Recalculate all reaction rates
        !print *, 'CMB: Calculating reaction rates: (ifgas, ifaq) = ', ifgas, ifaq
        CALL RecalculateReactionRates(IFGAS,IFAQ)
        !rint *, 'GasConc(599) after RecalculateReactionRates = ', GasConc(599)
        !print *, 'GridGasChem(599) after RecalculateReactionRates = ', GridGasChem(599)
         
        ! CMB: Recalc photo rates (new)
        call RecalculatePhotoRatesSZA(sza, Transmissivity, .false.)
        !print *, 'GasConc(599) = ', GasConc(599)
        !print *, 'GridGasChem(599) = ', GridGasChem(599)
        ! CMB: Debugging...
        !      print *, 'Solar Zenith Angle: ', sza
        !      do l = 1, GasPhaseReactionRateArraySize(1)
        !        print *, 'L = ', L
        !        if (int(GasPhaseReactionRates(L, 1)) .eq. 1 .and. &
        !                int(GasPhaseReactionRates(L,3)) .eq. 0) then
        !            RxnIndex = int(GasPhaseReactionRates(l,4))
        !            
        !            print *, '  Photolysis Reaction'
        !            select case(RxnIndex)
        !            case(2)
        !                print *, '    RxnIndex = 2, GasPhaseChemicalRates(L) = ', &
        !                        GasPhaseChemicalRates(l)
        !            case(6)
        !                print *, '    RxnIndex = 6, GasPhaseChemicalRates(L) = ', &
        !                        GasPhaseChemicalRates(l)
        !            end select
        !        endif   
        !      enddo
        !      stop
        ! end CMB
      		
        !WRITE(*,*) "Before Photorates"
        !Reset Photolysis Rates
        !	  DO L = 1, GasPhaseReactionRateArraySize(1)
        !						
        !	    !! Only do photolysis reactions
        !	    IF (INT(GasPhaseReactionRates(L,1)) .EQ. 1 .AND. &
        !                INT(GasPhaseReactionRates(L,3)) .EQ. 0.) THEN	
        !								
        !		!Find reaction index number
        !		RxnIndex = INT(GasPhaseReactionRates(L,4))
        !
        !		SELECT CASE (RxnIndex)
        !		  CASE(4) !NO2 => NO + O
        !		    GasPhaseChemicalRates(L) = PhotoRates(1)
        !		  CASE(5) !NO3 => NO + O2
        !		    GasPhaseChemicalRates(L) = PhotoRates(2)
        !		  CASE(6) !NO3 => NO2 + O
        !		    GasPhaseChemicalRates(L) = PhotoRates(3)		
        !		  CASE(3) !O3 => O + O2
        !	            GasPhaseChemicalRates(L) = PhotoRates(4)
        !		  CASE(2) !O3 => O1D + O2
        !		    GasPhaseChemicalRates(L) = PhotoRates(5)
        !		  CASE(12) !HONO => Products
        !		    GasPhaseChemicalRates(L) = PhotoRates(6)
        !		  CASE(11) !H2O2 => 2OH
        !		    GasPhaseChemicalRates(L) = PhotoRates(7)
        !		  CASE(15) !HCHO => 2HO2 + CO
        !		    GasPhaseChemicalRates(L) = PhotoRates(8)
        !		  CASE(16) !HCHO => H2 + CO
        !		    GasPhaseChemicalRates(L) = PhotoRates(9)
        !                  CASE(17) !ALD2 => Products
        !		    GasPhaseChemicalRates(L) = PhotoRates(10)
        !		  CASE(23) !KETL => Products
        !		    GasPhaseChemicalRates(L) = PhotoRates(11)
        !	          CASE(22) !MGLY => Products
        !		    GasPhaseChemicalRates(L) = PhotoRates(12)
        !		END SELECT
        !	     END IF
        !          END DO
        !WRITE(*,*) "After Photorates"

        if (doAerosols .gt. 0) then
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
            !call spillbeans(tout_unit=360)    
            !Calculate Heterogeneous rates from size dist. info
            CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HowManyBins, .FALSE.)
            !print *, 'GasConc(599) after RecalculateHeteroRates = ', GasConc(599) 
            !print *, 'GridGasChem(599) after RecalculateHeteroRates = ', GridGasChem(599)
            !do iii = 1, HowManyBins
            !    print *, 'IVOC1 (Org)[', iii, '] = ', AeroMassConc(iii, 56)
            !    print *, 'IVOC1 (AqOrg)[', iii, '] = ', AeroMassConc(iii, 78)
            !    print *, '    NumConc[', iii, '] = ', AeroNumConc(iii)
            !enddo      
         !WRITE(*,*) "Hetero Okay"
            !call spillbeans(tout_unit=362)
        endif
        !print *, 'GasConc(599) after InputConcentrations = ', GasConc(599) 
        !print *, 'GridGasChem(599) after inputConcentrations = ', GridGasChem(599)
    END SUBROUTINE InputConcentrations

    SUBROUTINE OutputConcentrations(Temp, Press, Dens, GasConc, &
        AeroNumConc, AeroMassConc, &
        ExtCoeff, SingScat, Assym, Radius, TermVel)
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
        REAL*8 :: GasConc(HowManyEvolveGasChems)!molecules/cm3
        REAL*8 :: AeroNumConc(HowManyBins)!particles/cm3
        REAL*8 :: AeroMassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems) !ug/m3
        REAL*8 :: Radius(HowManyBins) !cm
        REAL*8 :: TermVel(HowManyBins) !m/s
		
        !! Internal Variables
        INTEGER :: I, J, NumAq
        REAL*8 :: SaltConc(HowManyAqChems)
	
        TYPE(Particle), POINTER :: Cur
		
        ! CMB add
        integer :: doAerosols, iii != getAerosolProc()
        real*8 :: filler_value, ro2_t, rco3_t
        doAerosols = getAerosolProc()
        filler_value = 1.0D-40
        ! end CMB add
		
        !Output Temperature and Pressure
        Temp = GetTemp()
        Press = GetPress()
		
        ! CMB (AER): Bug fix: Apply peroxy getters
        !ro2_t = GetTotalPeroxy()
        !rco3_t = GetTotalAcylPeroxy()
        !write(678,*) 'In OutputConcentrations(): ro2_t = ', ro2_t
        !write(678,*) 'In OutputConcentrationS(): rco3_t = ', rco3_t
        ! end CMB add
        !print *, 'GasConc(599) at start of output concentrations = ', GasConc(599)
        !print *, 'GridGasChem(599) at start of output concentrations = ', GridGasChem(599)
        !Output Gas chemical concentrations
        DO I = 1, HowManyEvolveGasChems
            GasConc(I) = GridGasChem(I)
        END DO
        !print *, 'GasConc(599) = ', GasConc(599)
        !do iii = 1, HowManyBins
        !    print *, 'IVOC1 (Org)[', iii, '] = ', AeroMassConc(iii, 56)
        !    print *, 'IVOC1 (AqOrg)[', iii, '] = ', AeroMassConc(iii, 78)
        !    print *, '    NumConc[', iii, '] = ', AeroNumConc(iii)
        !enddo
        
        if (doAerosols .gt. 0) then
            !Update aerosol concentrations
            !write(*,*) 'Calling SortAerosolAtGridPointForCoagulation() from OutputConcentrations'
            
            CALL SortAerosolAtGridPointForCoagulation () 
            !print *, 'GasConc(599) after SortAerosolAtGridPointForCoagulation = ', GasConc(599)
            !print *, 'GridGasChem(599) after SortAerosolAtGridPointForCoagulation = ', GridGasChem(599)
            !do iii = 1, HowManyBins
            !    print *, 'IVOC1 (Org)[', iii, '] = ', AeroMassConc(iii, 56)
            !    print *, 'IVOC1 (AqOrg)[', iii, '] = ', AeroMassConc(iii, 78)
            !    print *, '    NumConc[', iii, '] = ', AeroNumConc(iii)
            !enddo
            
            cur => particles%first
            I = 1
            NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
            DO WHILE(associated(cur))
                AeroNumConc(I) = cur%numberofparticles
                Radius(I) = cur%effectiveradius
                !write(*,*) "Radius (um): ", Radius(I)*1e4
                TermVel(I) = TerminalVelocity(cur)/100.0 !Convert from cm/s to m/s
                !write(*,*) "Terminal velocity (m/s)", TermVel(I)
                !Water
                AeroMassConc(I,1) = cur%AqChems(1)*(&
                    cur%Numberofparticles*AqMolecularMass(1)) &
                    /(1.0e-12)
                               
                                !Force minimum water
                IF(AeroMassConc(I,1) .LT. 0.001) 	AeroMassConc(I,1) = 0.001						   
                
                !Do this step so that we only pass neutral electrolytes
                !to dynamics code - helps prevent odd charge balance problems
                ! CMB (uncomment write statements for now
                !write(*,*) 'AeroMassConc(i,1) = ', aeromassconc(i,1)
                !WRITE(*,*) "Before output HypotheticalElectrolyteConcentrations"
                SaltConc = HypotheticalElectrolyteConcentrations (cur)
                !WRITE(*,*) SaltConc(2), SaltConc(3), SaltConc(4), SaltConc(5), SaltConc(6)
                !WRITE(*,*) cur%AqChems(2), cur%AqChems(3), cur%AqChems(4), cur%AqChems(5), cur%AqChems(6) 
                !WRITE(*,*) findchem('NH3', 1)
                IF(cur%Numberofparticles .GT. 0.) then !filler_value) then !0.) THEN
                    !Be careful with NH3 to avoid double counting, MJA, 2016-06-14
                    DO J = 2, HowManyAqChems
                        if (J .EQ. findchem('NH3', 1)) then
                            AeroMassConc(I,J) = (Cur%AqChems(J))*(cur%Numberofparticles*AqMolecularMass(J)) &
                                /(1.0e-12)						
                        else 
                            AeroMassConc(I,J) = (SaltConc(J)+Cur%AqChems(J))*(cur%Numberofparticles*AqMolecularMass(J)) &
                                /(1.0e-12)
                        endif
                    END DO
                    
                    DO J = HowManyAqChems+1, NumAq
                        AeroMassConc(I,J) = 0.
                        !aeromassconc(i,j) = 1.0e-30
                    END DO
    
                    DO J = 1, HowManyOrgChems
                        AeroMassConc(I,J+NumAq) = Cur%OrgChems(J)*(cur%Numberofparticles*OrgMolecularMass(J)) &
                            /(1.0e-12)
                    END DO
                    
                    DO J = 1, HowManyAqOrgChems
                        AeroMassConc(I,J+NumAq+HowManyOrgChems) = Cur%AqOrgChems(J)*(cur%Numberofparticles*AqOrgMolecularMass(J)) &
                            /(1.0e-12)
                    END DO
                ELSE !No particles
                    AeroNumConc(I) = 0.
                    !aeronumconc(i) = 1.0e-30
                    DO J = 1, HowManyAqChems
                        AeroMassConc(I,J) = 0.
                        !aeromassconc(i,j) = 1.0e-30
                    END DO
                    
                    DO J = 1, HowManyOrgChems
                        AeroMassConc(I,J+HowManyAqChems) = 0.
                        !AeroMassConc(I,J+HowManyAqChems) = 1.0e-30
                    END DO
                    
                    DO J = 1, HowManyAqOrgChems
                        AeroMassConc(I,J+HowManyAqChems+HowManyOrgChems) = 0.
                        !AeroMassConc(I,J+HowManyAqChems+HowManyOrgChems) = 1.0e-30
                    END DO				
                END IF
            
                I = I + 1
                cur => cur%next
                    
            END DO
    
            !Output Average Aerosol optical properties, 0 flag assumes core in shell mixing
            !print *, 'calling aerosoloptprop within outputconcentrations()'
            call flush(6)
            !CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSSA, 0)
            !WRITE(*,*) ExtCoeff(1), SingScat(1), Assym(1)
            !PAUSE
        endif
        !print *, 'GasConc(599) at end of OutputConcentrations = ', GasConc(599)    
        !print *, 'GridGasChem(599) at end of OutputConcentrations = ', GridGasChem(599)
        !do iii = 1, HowManyBins
        !    print *, 'IVOC1 (Org)[', iii, '] = ', AeroMassConc(iii, 56)
        !    print *, 'IVOC1 (AqOrg)[', iii, '] = ', AeroMassConc(iii, 78)
        !    print *, '    NumConc[', iii, '] = ', AeroNumConc(iii)
        !enddo
    END SUBROUTINE OutputConcentrations

    !	SUBROUTINE InitializeASP(doAerosols, Transmissivity, sza)
    SUBROUTINE InitializeASP(num_weights, mass_weights, doAerosols, Transmissivity, sza)
    
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

        ! CMB add
        integer :: ii !, derp
        real*8, intent(in), optional :: sza, Transmissivity
        !logical, intent(in), optional :: doAerosols
        integer, intent(in), optional :: doAerosols
        !logical :: todoAerosols
        integer :: todoAerosols
            ! end CMB add

        !! Define the local variables
        integer :: i, j, k, l, seedarray(2),q, r, NumBins
        INTEGER :: time(8), sizerand
        REAL*4  :: etimearray(2)
        TYPE(Particle), POINTER :: Cur, Env
        CHARACTER (len = 8)     :: ErrorDate
        CHARACTER (len = 10)    :: ErrorTime
        ! 
        ! 20 September 2016 CMB: Define some variables needed for new computations
        !                        However, we don't know how big these are yet.
        !                        I can either make these allocatable or just large
        !                        static size.  I'm choosing the former for right now.
        !
        integer :: allo_stat, iter, num_aeros
        real*8 :: temp, press, dens, total_sum, row_sum
        real*8, dimension(451) :: ExtCoeff, SingScat, Assym
        real*8, allocatable, intent(out) :: num_weights(:), mass_weights(:)
        REAL*8, allocatable :: GasConc(:) !HowManyEvolveGasChems)
        REAL*8, allocatable :: NumConc(:) !HowManyBins)
        REAL*8, allocatable :: MassConc(:, :) !HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
        REAL*8, allocatable :: Radius(:) !HowManyBins)
        REAL*8, allocatable :: TermVel(:) !HowManyBins)
        ! end CMB
         
        ! CMB add
        todoAerosols = 0
        if (present(doAerosols)) todoAerosols = doAerosols
        !print *, 'todoAerosols = ', todoAerosols
        call setAerosolProc(todoAerosols)
        !derp = getAerosolProc()
        !if (derp .eq. 0) then
        !    print *, 'aerosol checker does not work'
        !    stop
        !endif
        ! end CMB add
    
        !! Initialize the random number generator 
        call DATE_AND_TIME(values=time)     
    
        ! Get the current time (for pgi compiler)
        seedarray(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
    
        !seedarray(1) = rtc() !for ifort compiler
        seedarray(2) = FLOOR(SQRT(REAL(seedarray(1))))

        ! CMB (AER): Modify this to play well with gfortran
        !CALL Random_Seed(put=seedarray)
        sizerand = 2
        CALL Random_Seed(size=sizerand)!, put=seedarray)

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !! SET the DOMAIN SIZE !!
        CALL SetDomainSize ()
	
        !! Call subroutines to set module data
        CALL SetFileHandleCounter(12)
        CALL ReadMainInputDeck()
	
        !! Send a Welcome Message
        CALL Transcript ("**********************************************************")
        CALL Transcript ("** Welcome to the ASP Model developed by Matt Alvarado  **")
        CALL Transcript ("** (mjalvara@mit.edu).                                  **")
        CALL Transcript ("** ASP is an updated, expanded version of MELAM         **")
        CALL Transcript ("** by H. D. Steele (donnan@mit.edu).                    **")
        CALL Transcript ("**********************************************************")
        CALL Date_and_Time(ErrorDate, ErrorTime)
        CALL Transcript ("*********************************************")
        CALL Transcript ("** Run Started at                          **")
    
            ! CB: Modified the below call so that vim doesn't have a syntax
            ! highlighting heart attack
        CALL Transcript ("** "//ErrorTime(1:2)//":"//ErrorTime(3:4)//":"//&
            ErrorTime(5:6)//" on "//ErrorDate(5:6)//&
            "/"//ErrorDate(7:8)//"/"//ErrorDate(1:4)//" "//" "//"**")
        CALL Transcript ("******************************************")
				
        ! INITIALIZE the CHEMISTRY
            !print *, 'CMB: Setting Chemistry Params' 
        CALL SetChemistryParams
        !print *, 'CMB: Done with SetChemistryParams'

        ! CB: This is NEEDED and done only within
        ! DevelopmentDriver.f90, which is the main().  We're
        ! not calling that in STILT-ASP...
        !print *, 'CMB: Reading PHotolysis Rates'
        call ReadPhotolysisRates
        !print *, 'CMB: Recalculating Reaction Rates'
        call RecalculateReactionRates(ifgas, ifaq)
        ! end cb; we may need more...

        ! CMB: call recalculate photo rates sza after call of recalculatereactionrates
        !      with "FirstTime" set to TRUE
        if (present(sza)) call RecalculatePhotoRatesSZA(sza, Transmissivity, .TRUE.)
            ! end CMB add
    
        !! INITIALIZE the DISSOLUTION routine
            !print *, 'CMB: Setting all dissolution...'
        CALL SetAllDissolution

        !! INITIALIZE the PARTICLES
            !print *, 'CMB: initializing the particles'        
        CALL InitializeParticlesSectional
        !
        ! 20 September 2016 CMB: According to Matt Alvarado, the following functions are 
        !                        required to initialize the background.  Took these from
        !                        DevelopmentDriver.f90.  Hopefully I put them in the
        !                        right place... 
        !! Call Boundary distributions
        CALL ReadBoundaryDistributionsOrganic
        CALL PopulateParticlesSectionsRightAway(.TRUE.)
        CALL SortBoundaryPart
    				
        !!WARNING Force initial water equilibrium between gas and aerosol
        !!(Uses open system for water eq.)
        !CALL EquilibrateInternallyAtGridPoint(1,1,1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
            !print *, 'CMB: Regridding aerosol'
        
        CALL RegridAerosol ()
	
        !
        ! 20 September 2016 CMB: Matt says to call OutputConcentrations to populate
        !                        the aerosol    
        num_aeros = getHowManyAqChems() + getHowManyAqCations() + getHowManyAqAnions() + &
            getHowManyOrgChems() + getHowManyAqOrgChems() 
        !print *, 'CB: Total Number of Aerosols: ', num_aeros
                     
        temp = GetTemp()        ! GetTempFromGrid?
        press = GetPress()      ! GetPressFromGrid?
        dens = GetAirDensity()  ! ?
    
        !print *, 'CB: Temp = ', temp
        !print *, 'CB: Pressure = ', press
        !print *, 'CB: Density = ', dens
    
        allocate(GasConc(HowManyEvolveGasChems), stat=allo_stat)
        if (allo_stat .ne. 0) call error("Cannot allocate GasConc()")
    
        allocate(NumConc(HowManyBins), num_weights(HowManyBins), &
            mass_weights(HowManyBins), stat=allo_stat)
        if (allo_stat .ne. 0) call error("Cannot allocate NumConc and/or weights")
    
        allocate(MassConc(HowManyBins, num_aeros), stat=allo_stat)
        if (allo_stat .ne. 0) call error("Cannot allocate MassConc")
    
        allocate(Radius(HowManyBins), TermVel(HowManyBins), stat=allo_stat)
        if (allo_stat .ne. 0) call error("Cannot allocate NumConc and/or weights")
    
        GasConc = 0.0D0
        NumConc = 0.0D0
        MassConc = 0.0D0
        Radius = 0.0D0
        TermVel = 0.0D0
        total_sum = 0.0D0
        row_sum = 0.0D0
    
        call OutputConcentrations(temp, press, dens, GasConc, NumConc, &
            MassConc, ExtCoeff, SingScat, Assym, Radius, &
            TermVel)
    
        ! NumConc is in #particles/cm3.  MassConc is in ug/m3.  The "weights" factor is 
        ! in units of #particles/ug after this conversion?
        !
        ! UPDATE: now including a mass factor too                          
        total_sum = sum(MassConc(:,:))
        !print *, 'CB: Total Mass of all particles = ', total_sum
        if (total_sum .lt. 1e-30) then
            call error("MassConc() sum is < 1e-30 after call to OutputConcentrations")
        else
            num_weights(:) = NumConc(:)*1.0e6/total_sum
            do iter = 1, HowManyBins
                mass_weights(iter) = sum(MassConc(iter, :))/total_sum
            enddo    
        endif
    
        ! TODO: Deallocate these, or return to STILT-ASP?
        deallocate(NumConc)
        deallocate(MassConc)
        deallocate(GasConc)
        deallocate(Radius)
        deallocate(TermVel)    
    ! end CMB add                               
    
    END SUBROUTINE InitializeASP

    SUBROUTINE ASPInterface(Timestep, tTemp, Press, Dens, GasConc, &
        NumConc, MassConc, &
        !PhotoRates, &
        ExtCoeff, SingScat, Assym, Radius, TermVel, sza, Transmissivity, itt, jp, &
        do_aerosols) !, &
        !it, jp)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This subroutine connects ASP to an external 3D Eulerian model!! 
        !! MJA 062707							!!
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
        REAL*8 :: Timestep !(in s)
        REAL*8 :: tTemp  !(in K)
        REAL*8 :: Press !(in mbar)
        REAL*8 :: Dens  !(in kg/m3 dry air)
        REAL*8 :: GasConc(HowManyEvolveGasChems)
        REAL*8 :: NumConc(HowManyBins)
        REAL*8 :: MassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
        REAL*8 :: Radius(HowManyBins)
        REAL*8 :: TermVel(HowManyBins)
        integer :: itt, jp
        !REAL*8 :: PhotoRates(12)
		
        ! CMB add: 
        real*8, intent(in) :: sza, Transmissivity
        integer, intent(in), optional :: do_aerosols
        type(Particle), POINTER :: cur
        real*8 :: total_num_part
        !logical :: want_aerosols
        ! end CMB add

        ! CB: Ifort complains unless these dimensions match
        !REAL*8, DIMENSION(18) :: ExtCoeff, SingScat, Assym
        REAL*8, DIMENSION(451):: ExtCoeff, SingScat, Assym
        INTEGER :: I, J, K

        ! CMB add
        !want_aerosols = .false.
        if (present(do_aerosols)) call setAerosolProc(do_aerosols)
	            ! end CMB add

        !write(345, *) 'Inputs into ASP...'
        !write(345, *) '------------------'
        !write(345, *) 'Temperature: ', tTemp
        !write(345, *) 'Pressure: ', Press
        !write(345, *) 'Density: ', Dens
        !write(345, *) 'is, jp = ', itt, jp
        !write(345, *) 'Timestep: ', timestep
        !write(345, *) 'Solar Zenith Angle: ', sza
        !write(345, *) 'Transmissivity: ', transmissivity
        !write(345, *) 'Process Aerosols? ', getAerosolProc()
        !write(345, *) 'GasPhase (molec/cm3)'
        !write(345, '(A5, 1x, A16, 1x, A25)') 'Index', 'Species Name', 'molec/cm3'
        !do i = 1, HowManyEvolveGasChems
        !write(345, '(I5, 1x, A16, 1x, E25.15)'), I, GasPhaseChemicalNames(i)(1:len(trim(&
        !        GasPhaseChemicalNames(i)))), GasConc(I) 
        !enddo
!        write(345, *) 'AerosolPhase (ug/m3)'
!        write(345, *) 'Aerosol Concentrations'
!        write(345, '(A5,1x,A16,1x,A18,9(1x,A10))') 'Index', 'Species Name', 'Mass (ug/m3) Bin 1', 'Bin 2', &
!                'Bin 3', 'Bin 4', 'Bin 5', 'Bin 6', 'Bin 7', 'Bin 8', 'Bin 9', 'Bin 10'
!                                
!        do i = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems
!            if (i .le. HowManyAqChems) then
!                j = i
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, AqPhaseChemicalNames(j)(1:len(trim(AqPhaseChemicalNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            else if (i .le. HowManyAqChems + HowManyAqCations) then
!                j = i - HowManyAqChems
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, AqCationNames(j)(1:len(trim(AqCationNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            else if (i .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions) then
!                j = i - HowManyAqChems - HowManyAqCations
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, AqAnionNames(j)(1:len(trim(AqAnionNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            else if (i .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions + HowManyOrgChems) then
!                j = i - HowManyAqChems - HowManyAqCations - HowManyAqAnions
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, OrgPhaseChemicalNames(j)(1:len(trim(OrgPhaseChemicalNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            else
!                j = i - HowManyAqChems - HowManyAqCations - HowManyAqAnions - HowManyOrgChems
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, AqOrgPhaseChemicalNames(j)(1:len(trim(AqOrgPhaseChemicalNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            endif           
!        enddo
!        do i = 1, HowManyBins
!            write(345, '(A,I2,A,E20.5)') 'Total Num Conc (#particles/cm3, bin ', i, '): ', &
!                    NumConc(i)
!        enddo
!        write(345, *) 'Total Number of Particles: ', sum(NumConc(:))

!                            write(345, '(A5,1x,A16,1x,A25,1x,A25)') 'Index', 'Species Name', 'Prior to ASP', 'After ASP'
!                            write(345, '(A5,1x,A16,1x,A25,1x,A25)') '-----', '------------', '------------', '---------'
!                            do i = 1, num_asp_aerosols
!                                write(345, '(I5,1x,A16,1x,F25.15,1x,F25.15)') &
!                                        i, AerosolPhaseChemicalNames(i)(1:len(trim(&
!                                        AerosolPhaseChemicalNames(i)))), previous_time_asp_aerosols_massconc(jp,1,i), &
!                                        tempMassConc(1,i)
!                            enddo
!                            write(345,'(A5,1x,A16,F25.15,1x,F25.15)') 'N/A', 'NumConc', &
!                                    previous_time_asp_aerosols_numconc(jp,1), &
!                                    tempNumConc(1)
!                            write(345, '(A,E20.5)') 'Total Num Conc (#particles/cm3): ', &
!                                    previous_time_asp_aerosols_numconc(jp, 1)     
        !write(*,*) "MassConc(1,1) = ", MassConc(1,1)
        !call spillbeans(tout_unit=199)
		!write(*,*) 'Sum of NumConc(:) upon entering ASPInterface: ', sum(NumConc(:))
		!write(*,*) 'In StepASP: itt, jp = ', itt, jp
		call setTimestepIndex(itt)
		call setParcelIndex(jp)
		!write(*,*) 'After Setters: itt, jp = ', getTimeStepIndex(), getParcelIndex()
		
        CALL InputConcentrations(tTemp, Press, Dens, GasConc, &
            NumConc, MassConc, &
            !PhotoRates, 
            sza, Transmissivity)
        !stop
!        call flush(6)
!        
!        write(345, *) 'After InputConcentrations...'
!        write(345, *) '----------------------------'
!        write(345, *) 'GasPhase (molec/cm3)'
!        write(345, '(A5, 1x, A16, 1x, A25)') 'Index', 'Species Name', 'molec/cm3'
!        do i = 1, HowManyEvolveGasChems
!        write(345, '(I5, 1x, A16, 1x, E25.15)'), I, GasPhaseChemicalNames(i)(1:len(trim(&
!                GasPhaseChemicalNames(i)))), GridGasChem(I) 
!        enddo
!        write(345, *) 'AerosolPhase (ug/m3)'
!        write(345, *) 'Aerosol Concentrations'
!        !write(345, '(A5,1x,A16,1x,A18,9(1x,A10))') 'Index', 'Species Name', 'Mass (ug/m3) Bin 1', 'Bin 2', &
!        !        'Bin 3', 'Bin 4', 'Bin 5', 'Bin 6', 'Bin 7', 'Bin 8', 'Bin 9', 'Bin 10'
!        
!        total_num_part = 0.0D0
!        cur => particles%first
!        i = 1
!        do while (associated(cur))
!            write(345, *) '--> Particle/Bin Index I = ', I
!            write(345, *) 'Number of Particles: ', cur%numberofparticles
!            total_num_part = total_num_part + cur%numberofparticles
!            
!            do k = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions + HowManyOrgChems + &
!                      HowManyAqOrgChems
!                if (k .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions) then
!                    j = k
!                    write(345, '(4x, A16, 1x, E10.3)') AqPhaseChemicalNames(j), cur%aqchems(j)
!                else if (k .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions + HowManyOrgChems) then
!                    j = k - HowManyAqChems - HowManyAqCations - HowManyAqAnions
!                    write(345, '(4x, A16, 1x, E10.3)') OrgPhaseChemicalNames(j), cur%orgchems(j)
!                else
!                    j = k - HowManyAqChems - HowManyAqCations - HowManyAqAnions - HowManyOrgChems
!                    write(345, '(4x, A16, 1x, E10.3)') AqOrgPhaseChemicalNames(j), cur%aqorgchems(j)
!                endif          
!            enddo
!            i = i + 1  
!            cur => cur%next         
!        enddo                        
!        write(345, *) 'Total # of Particles: ', total_num_part
!        
        !call GasOutput(300, 301, 302, timestep, .true., .true.)
        !print *, 'GridGasChem(line 918, StepASP.f90) = ', GridGasChem(599)	
        !call spillbeans(tout_unit=203)						
        !WRITE(678,*) "ASPInterface: After InputConcentrations(), calling StepASpOnce()"
		
        CALL StepASPOnce (TimeStep)
        
!        write(345, *) 'After StepAspOnce...'
!        write(345, *) '----------------------------'
!        write(345, *) 'GasPhase (molec/cm3)'
!        write(345, '(A5, 1x, A16, 1x, A25)') 'Index', 'Species Name', 'molec/cm3'
!        do i = 1, HowManyEvolveGasChems
!        write(345, '(I5, 1x, A16, 1x, E25.15)'), I, GasPhaseChemicalNames(i)(1:len(trim(&
!                GasPhaseChemicalNames(i)))), GridGasChem(I) 
!        enddo
!        write(345, *) 'AerosolPhase (ug/m3)'
!        write(345, *) 'Aerosol Concentrations'
!        
!        total_num_part = 0.0D0
!        cur => particles%first
!        i = 1
!        do while (associated(cur))
!            write(345, *) '--> Particle/Bin Index I = ', I
!            write(345, *) 'Number of Particles: ', cur%numberofparticles
!            total_num_part = total_num_part + cur%numberofparticles
!            
!            do k = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions + HowManyOrgChems + &
!                      HowManyAqOrgChems
!                if (k .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions) then
!                    j = k
!                    write(345, '(4x, A16, 1x, E10.3)') AqPhaseChemicalNames(j), cur%aqchems(j)
!                else if (k .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions + HowManyOrgChems) then
!                    j = k - HowManyAqChems - HowManyAqCations - HowManyAqAnions
!                    write(345, '(4x, A16, 1x, E10.3)') OrgPhaseChemicalNames(j), cur%orgchems(j)
!                else
!                    j = k - HowManyAqChems - HowManyAqCations - HowManyAqAnions - HowManyOrgChems
!                    write(345, '(4x, A16, 1x, E10.3)') AqOrgPhaseChemicalNames(j), cur%aqorgchems(j)
!                endif          
!            enddo
!            i = i + 1 
!            cur => cur%next          
!        enddo                        
!        write(345, *) 'Total # of Particles: ', total_num_part
        
        !call spillbeans(tout_unit=203)
        !call GasOutput(304, 205, 206, timestep, .true., .true.)
        !print *, 'GridGasChem(line 925, StepASP.f90) = ', GridGasChem(599)
        !WRITE(678,*) "ASPInterface: After StepASPOnce(); calling OutputConcentrations()" !Before Conc Output; AFTER STEPPING"
        !write(*,*) "After step ASP once: aeromassconc(1,1) = ", massconc(1,1)
        !call flush(6)
        CALL OutputConcentrations(tTemp, Press, Dens, GasConc, &
            NumConc, MassConc, &
            ExtCoeff, SingScat, Assym, Radius, TermVel)
            
!        write(345, *) 'After OutputConcentrations...'
!        write(345, *) '----------------------------'
!        write(345, *) 'GasPhase (molec/cm3)'
!        write(345, '(A5, 1x, A16, 1x, A25)') 'Index', 'Species Name', 'molec/cm3'
!        do i = 1, HowManyEvolveGasChems
!        write(345, '(I5, 1x, A16, 1x, E25.15)'), I, GasPhaseChemicalNames(i)(1:len(trim(&
!                GasPhaseChemicalNames(i)))), GasConc(I) 
!        enddo
!        write(345, *) 'AerosolPhase (ug/m3)'
!        write(345, *) 'Aerosol Concentrations'
!        write(345, '(A5,1x,A16,1x,A18,9(1x,A10))') 'Index', 'Species Name', 'Mass (ug/m3) Bin 1', 'Bin 2', &
!                'Bin 3', 'Bin 4', 'Bin 5', 'Bin 6', 'Bin 7', 'Bin 8', 'Bin 9', 'Bin 10'
!        do i = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems
!            if (i .le. HowManyAqChems) then
!                j = i
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, AqPhaseChemicalNames(j)(1:len(trim(AqPhaseChemicalNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            else if (i .le. HowManyAqChems + HowManyAqCations) then
!                j = i - HowManyAqChems
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, AqCationNames(j)(1:len(trim(AqCationNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            else if (i .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions) then
!                j = i - HowManyAqChems - HowManyAqCations
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, AqAnionNames(j)(1:len(trim(AqAnionNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            else if (i .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions + HowManyOrgChems) then
!                j = i - HowManyAqChems - HowManyAqCations - HowManyAqAnions
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, OrgPhaseChemicalNames(j)(1:len(trim(OrgPhaseChemicalNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            else
!                j = i - HowManyAqChems - HowManyAqCations - HowManyAqAnions - HowManyOrgChems
!                write(345, '(I5, 1x, A16, 1x, E18.5, 9(1x,E10.3))') &
!                        i, AqOrgPhaseChemicalNames(j)(1:len(trim(AqOrgPhaseChemicalNames(j)))), &
!                        MassConc(1, i), &
!                        MassConc(2, i), &
!                        MassConc(3, i), &
!                        MassConc(4, i), &
!                        MassConc(5, i), &
!                        MassConc(6, i), &
!                        MassConc(7, i), &
!                        MassConc(8, i), &
!                        MassConc(9, i), &
!                        MassConc(10, i)
!            endif           
!        enddo
!        do i = 1, HowManyBins
!            write(345, '(A,I2,A,E20.5)') 'Total Num Conc (#particles/cm3, bin ', i, '): ', &
!                    NumConc(i)
!        enddo
!        write(345, *) 'Total Number of Particles: ', sum(NumConc(:))
        
        
        
        
        
        
        
!        write(345, *) 'After OutputConcentrations...'
!        write(345, *) '----------------------------'
!        write(345, *) 'GasPhase (molec/cm3)'
!        write(345, '(A5, 1x, A16, 1x, A25)') 'Index', 'Species Name', 'molec/cm3'
!        do i = 1, HowManyEvolveGasChems
!        write(345, '(I5, 1x, A16, 1x, E25.15)'), I, GasPhaseChemicalNames(i)(1:len(trim(&
!                GasPhaseChemicalNames(i)))), GasConc(I) 
!        enddo
!        write(345, *) 'AerosolPhase (ug/m3)'
!        write(345, *) 'Aerosol Concentrations'
!        write(345, '(A5,1x,A16,1x,A18,9(1x,A10))') 'Index', 'Species Name', 'Mass (ug/m3) Bin 1', 'Bin 2', &
!                'Bin 3', 'Bin 4', 'Bin 5', 'Bin 6', 'Bin 7', 'Bin 8', 'Bin 9', 'Bin 10'
!                                
!       total_num_part = 0.0D0
!        cur => particles%first
!        i = 1
!        do while (associated(cur))
!            write(345, *) '--> Particle/Bin Index I = ', I
!            write(345, *) 'Number of Particles: ', cur%numberofparticles
!            total_num_part = total_num_part + cur%numberofparticles
!            
!            do k = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions + HowManyOrgChems + &
!                      HowManyAqOrgChems
!                if (k .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions) then
!                    j = k
!                    write(345, '(4x, A16, 1x, E10.3)') AqPhaseChemicalNames(j), cur%aqchems(j)
!                else if (k .le. HowManyAqChems + HowManyAqCations + HowManyAqAnions + HowManyOrgChems) then
!                    j = k - HowManyAqChems - HowManyAqCations - HowManyAqAnions
!                    write(345, '(4x, A16, 1x, E10.3)') OrgPhaseChemicalNames(j), cur%orgchems(j)
!                else
!                    j = k - HowManyAqChems - HowManyAqCations - HowManyAqAnions - HowManyOrgChems
!                    write(345, '(4x, A16, 1x, E10.3)') AqOrgPhaseChemicalNames(j), cur%aqorgchems(j)
!                endif          
!            enddo
!            i = i + 1  
!            cur => cur%next         
!        enddo                        
!        write(345, *) 'Total # of Particles: ', total_num_part
        
        !write(*,*) 'Sum of NumConc(:) upon leaving ASPInterface: ', sum(NumConc(:))
        !call spillbeans(tout_unit=202)
        !print *, 'GridGasChem(line 933, StepASP.f90) = ', GridGasChem(599)
            !call GasOutput(207, 208, 209, timestep, .true., .true.)
	    !WRITE(*,*) "After Conc Output"
	    !write(*,*) 

    END SUBROUTINE ASPInterface
END MODULE StepASP
