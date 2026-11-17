!! ASP (c), 2004-2013, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ReadInputDeck.f90
!! This file holds parameters and constants germaine to the 
!! aerosol definitions and internal processes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY			  			             !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 11/06  2012   Matt Alvarado     Updated main input deck to give model more!!
!!                                 flexibility without recompiling.          !!
!! 07/03  2013   Matt Alvarado     Updated main input deck to allow interp   !!
!!                                 between two arbitrary photolysis points   !!
!! 07/03  2013   Matt Alvarado     Updated main input deck to allow input    !!
!!                                   of CO emission rate to Eulerian Box     !!
!!                                   routines.                               !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 189da4e: N11 input-deck switch: !!
!!                                 GDTableMode (1 = table, 0 = full      !!
!!                                 DGAUS8 integral) and GDTableCells,    !!
!!                                 optional last lines of ASPInputDeck.in!!
!! 11/02  2026   Claude Opus 5.5   Patch 28: optional ASPInputDeck.in line!!
!!                                 after GDTableCells sets AqInnerSolver  !!
!! 11/02  2026   Claude Opus 5.5   Patch 30: optional ASPInputDeck.in line!!
!!                                 after AqInnerSolver sets AqOuterSolver !!
!! 11/17  2026   Claude Opus 5.5   Optional "Name = value" settings at the!!
!!                                 end of the deck (ReadOptionalDeckSet-  !!
!!                                 tings); settings echoed to transcript  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES		                    !!
!! 1. Chemistry	       	                    !!
!! 2. InfrastructuralCode		    !!
!! 3. GridPointFields			    !!				
!! 4. ModelParameters			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	    !!	
!! 1. SUBROUTINE ReadMainInputDeck()				    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!!
!! CMB: Modified the conditionals so that ints were no longer compared
!!      against bools
!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! READ INPUT DECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Read the main input deck and set the appropriate !!
	!! parameters throughout the model.		    !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ReadMainInputDeck ()

		USE Chemistry,	         ONLY : SetPhaseIfs

		USE InfrastructuralCode, ONLY : GetFileHandle,		&
						ReturnFileHandle,	&
						GetLine,		&
						GetToken,		&
						IsInteger,		&
						IsReal,			&
						STR2INT,		&
						STR2REAL,		&
						StripToken,		&
						EOF,			&
						ERROR,WARN,INT2STR,	&
						Transcript

		USE GridPointFields,     ONLY : InitializeGridPointFields, &
						SetWaterVaporField 
		USE ModelParameters

		IMPLICIT NONE

		INTEGER   :: I,FH,FHT, Q
		LOGICAL   :: IFGAS,IFAQ,IFSOLID
		CHARACTER(len=512) :: Token
		CHARACTER(len=512)::CurrLine,StoreLine
		REAL*8    :: T,P,II, RH, BRH, Num, Wind, length, width
                REAL*8  :: lights_on_time, lights_off_time
                INTEGER :: lightswitch_flag

		FHT = GetFileHandle()
		FH  = GetFileHandle()
	        OPEN(UNIT=FH, FILE=TRIM(MainDeckDir)//'ASPInputDeck.in', STATUS='OLD')

		!!! -- Specify Input Deck Location
		CALL SetInputDeckSubDir(TRIM(MainDeckDir))
		
		!!! -- Specify Output Deck Location
		CALL SetOutputDeckSubDir(TRIM(MainOutputDir))
                CALL SetTranscriptFH(FHT)
		OPEN(UNIT=FHT, FILE=TRIM(MainOutputDir)//"Transcript.TXT")
		


                !!SECTION 1. PROCESS CONTROL SETTINGS
                !!
                !!Here we decide which processes (gas-phase chemistry, condensation, 
                !!coagulation, etc.) are modeled and how they are calculated.

		!!! -- Specify Which Phases have kinetic chemistry
		CurrLine = GetLine(FH)

			!!! Gas Phase
                !print *, 'CMB: Gas Phase'
			CALL GetToken(CurrLine,";",Token)

                        !print *, 'Token = "', Token(1:len(trim(token))), '"'
                        ! CMB (AER): Modifed to pass a certain portion
                        ! of the token
                        ! UPDATE: No, I just changed the IsInteger
                        ! calls, they should have "nots" when errors are
                        ! sought
			if (.not. isInteger(token)) &
				CALL ERROR("Must Supply an Integer in ASPInputDeck.in to indicate "//&
						   "whether Gas Phase Chemistry Should be Included")
			I = STR2INT(Token)
                !print *, 'within routine - I = ', I, ' - within routine'
			IF (I < 0 .OR. I > 1) &
				CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
						   "whether Gas Phase Chemistry Should be Included.  Received "//StripToken(Token))
			IF (I .EQ. 0) THEN
				IFGAS = .FALSE.
				CALL SetGasChemistryFlag (.FALSE.)
			ELSE
				IFGAS = .TRUE. 
				CALL SetGasChemistryFlag (.TRUE.)
			END IF

			!!! Aqueous Phase
                            !print *, 'CMB: Aqueous Phase;'
			CALL GetToken(CurrLine,";",Token)

			if (.not. isInteger(token)) &
				CALL ERROR("Must Supply an Integer in ASPInputDeck.in to indicate"//&
						   "whether Aqueous Phase Chemistry Should be Included")
			I = STR2INT(Token)
			IF (I < 0 .OR. I > 1) &
				CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
						   "whether Aqueous Phase Chemistry Should be Included.  Received "//StripToken(Token))
			IF (I .EQ. 0) THEN
				IFAQ = .FALSE.
				CALL SetAqChemistryFlag (.FALSE.)
			ELSE
				IFAQ = .TRUE.
				CALL SetAqChemistryFlag (.TRUE.)
			END IF
                !print *, 'CMB: SetPhaseIfs'
		CALL SetPhaseIfs(IFGAS,IFAQ)
		
		!!! Poly or Mono Disperse?
                !print *, 'CMB: Poly or Mono disperse?'
		Token = GetLine(FH)

		if (.not. isInteger(token)) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Aeresol are Monodisperse.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Aeresol are Monodisperse.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetMonodisperse(.FALSE.)
		ELSE
			CALL SetMonodisperse(.TRUE.)
		END IF

		!!! Condensation?
                !print *, 'CMB: Condensation?'
		Token = GetLine(FH)

                !print *, 'Token = "', token(1:len(trim(token))), '"'
                !print *, 'ISInteger (token) = ', isInteger(token)
		IF (.not. IsInteger(Token)) then ! .EQ. 0) then !&
                    !print *, 'is integer is not 0'
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Condensation is considered.  Received "//StripToken(Token))
                endif
		I = STR2INT(Token)
                !print *, 'I = ', I
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Condensation is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetCondensationFlag(.FALSE.)
		ELSE
			CALL SetCondensationFlag(.TRUE.)
		END IF

		!!! Inorganic Dissolution?
                !print *, 'Inorganic dissolution?'
		Token = GetLine(FH)

		if (.not. isInteger(token)) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Dissolution is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Dissolution is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetDissolutionFlag(.FALSE.)
		ELSE
			CALL SetDissolutionFlag(.TRUE.)
		END IF

		!!! Hydrophobic Organic Dissolution?
                !print *, 'hydrophobic organic dissolution?'
		Token = GetLine(FH)

		if (.not. isInteger(token)) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Hydrophobic Organic Dissolution is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Hydrophobic Organic Dissolution is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetHydrophobicOrgDissolutionFlag(.FALSE.)
		ELSE
			CALL SetHydrophobicOrgDissolutionFlag(.TRUE.)
		END IF

		!!! Hydrophilic Organic Dissolution?
                !print *, 'CMB: hydrophilic organic dissolution?'
		Token = GetLine(FH)

		if (.not. isInteger(token)) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Hydrophilic Organic Dissolution is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Hydrophilic Organic Dissolution is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetHydrophilicOrgDissolutionFlag(.FALSE.)
		ELSE
			CALL SetHydrophilicOrgDissolutionFlag(.TRUE.)
		END IF

		!! Is dissolution equilibrium or flux?
		Token = GetLine(FH)
                !print *, 'IS dissolution equilibrium or flux?'
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether disolution is an equilibrium 0 or flux 1 process. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether disolution is an equilibrium (0) or flux (1) process. Received "//StripToken(Token))		
		CALL SetDissolutionEquilibriumOrFlux (I)
		
		!! Is Ammonia dissolution equilibrium or flux?
                !print *, 'CMB: ammonium dissolution equilib or flux?'
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether disolution is an equilibrium 0 or flux 1 process. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether ammonia disolution is an equilibrium (0) or flux (1) process. Received "//StripToken(Token))	
		IF (I .EQ. 0) THEN
			CALL SetAmmoniaFlux(.FALSE.)
		ELSE
			CALL SetAmmoniaFlux(.TRUE.)
		END IF
		
		
		!!Metastable Inorganic Aerosol?
                !print *, 'CMB: MEtastable inorganic aerosol?'
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether inorganic aerosol is metastable (1) or solid salts are allowed (0). Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether inorganic aerosol is metastable (1) or solid salts are allowed (0). Received "//StripToken(Token))	
		IF (I .EQ. 0) THEN
			CALL SetMetastable(.FALSE.)
		ELSE
			CALL SetMetastable(.TRUE.)
		END IF

		!!Ignore corrections to diffusivity? (Also sets curvature corrections to 1.0!)
		Token = GetLine(FH)
                !print *, 'CMB: Ignore corrections to diffusivity?'
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether diffusivity corrections are ignored (1) or included (0). Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether diffusivity corrections are ignored (1) or included (0).. Received "//StripToken(Token))	
		IF (I .EQ. 0) THEN
			CALL SetIgnoreDiffuseCorrect(.FALSE.)
		ELSE
			CALL SetIgnoreDiffuseCorrect(.TRUE.)
		END IF
		
		!!! Coagulation?
		Token = GetLine(FH)
                !print *, 'CMB: Coagulation?'
		if (.not. isInteger(token)) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Coagulation is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Coagulation is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetCoagulationFlag(.FALSE.)
		ELSE
			CALL SetCoagulationFlag(.TRUE.)
		END IF

		!!Constant Coagulation Kernel?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether coagulationkernel is constant (1) or calculated (0) process. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether coagulation kernel is constant (1) or calculated (0) process. Received "//StripToken(Token))	
		IF (I .EQ. 0) THEN
			CALL SetConstantKernel(.FALSE.)
		ELSE
			CALL SetConstantKernel(.TRUE.)
		END IF

                !!**************************************************************************
                !!SECTION 2. RUN CONTROL SETTINGS
                !!
                !!Here we decide between a Langrangian parcel run, Eulerian box run, smog
                !!chamber run, or other special run formats.

	        !!Lagrangian Parcel Model?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if the Lagrangian Parcel Model is to be used. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if the Lagrangian Parcel Model is to be used. Received "//StripToken(Token))	             
                CALL SetLagrangianParcelModel(I)

	        !!Eulerian Box Model?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if the Eulerian Box Model is to be used. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if the Eulerian Box Model is to be used. Received "//StripToken(Token))	             
                CALL SetEulerianBoxModel(I)

	        !!Smog Chamber Model?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if the Smog Chamber Model is to be used. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if the Smog Chamber Model is to be used. Received "//StripToken(Token))	             
                CALL SetSmogChamberModel(I)

	        !!Optical Properties Model?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if the Optical Properties Model is to be used. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if the Optical Properties Model is to be used. Received "//StripToken(Token))	             
                CALL SetOpticalPropertiesModel(I)

	        !!Thermo Test?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if a Thermodynamics Test is to be run. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if a Thermodynamics Test is to be run. Received "//StripToken(Token))	             
                CALL SetThermoTest(I)

	        !!Coag Test?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if a Coagulation Test is to be run. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if a Coagulation Test is to be run. Received "//StripToken(Token))	             
                CALL SetCoagTest(I)


	        !!Cond Test?
		Token = GetLine(FH)
                write(*,*) Token
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate "//&
			"if a Condensation Test is to be run. Received "//StripToken(Token))	
		I = STR2INT(Token)
                write(*,*) I
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate "//&
			"if a Condensation Test is to be run. Received "//StripToken(Token))	             
                CALL SetCondTest(I)
 
                !!Subroutine Model?
		Token = GetLine(FH)
                write(*,*) Token
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate "//&
			"if ASP is to be run as a subroutine. Received "//StripToken(Token))	
		I = STR2INT(Token)
                write(*,*) I
		IF (I .LT. 0 .OR. I .GT. 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate "//&
			"if ASP is to be run as a subroutine. Received "//StripToken(Token))	             
                CALL SetSubroutineModel(I)

                !!*************************************************************************
                !!SECTION 3. RUN PARAMETERS
                !!
                !!Here we set temperature, pressure, RH, and other parameters needed
                !!for the run types given above
                !
                ! Modified by CMB (AER): Had to put "not" in front of
                ! each in order to reflect changes in STILT-ASP

		!! -- Get TEMPERATURE
		Token = GetLine(FH)

		!! Assume it is an input deck if not a real
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load Temperature from an input "//&
					   "deck.  As much as I'd like to oblige, that isn't implemented yet.")
		ELSE
			T = STR2REAL(Token)
		END IF

		!! -- Get PRESSURE
		Token = GetLine(FH)

		!! Assume it is an input deck if not a real
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load Pressure from an input "//&
					   "deck.  As much as I'd like to oblige, that isn't implemented yet.")
		ELSE
			P = STR2REAL(Token)
		END IF

		!! -- Get INITIAL RELATIVE HUMIDITY
		Token = GetLine(FH)

		!! Assume it is an input deck if not a real
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load RH from an input "//&
					   "deck.  As much as I'd like to oblige, that isn't implemented yet.")
		ELSE
			II = STR2REAL(Token)
			RH = II
		END IF

		!! -- Get BACKGROUND RELATIVE HUMIDITY
		Token = GetLine(FH)

		!! Assume it is an input deck if not a real
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load RH from an input "//&
					   "deck.  As much as I'd like to oblige, that isn't implemented yet.")
		ELSE
			II = STR2REAL(Token)
			BRH = II
		END IF

		!! Set the grid structure, t and p
		CALL InitializeGridPointFields (T,P)
		CALL SetWaterVaporField (RH)
		CALL SetBackgroundRH(BRH)

                !!Get PHOTO FLAG
                Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			" how photolysis is to be calculated. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			" how photolysis is to be calculated. Received "//StripToken(Token))

	        CALL SetPhotoFlag(I)   

                !!Get Start time (hours in UTC)
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added "not"
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real start time.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetStartTime(II)
		END IF

                !!Get End time for phtolysis interpolation (hours in UTC)
                !!Only used if PhotoFlag = 1
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real end time for photolysis interpolation.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetPhotoEndTime(II)
		END IF

	        !!Get Run Length (in min)
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real run length.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetRunLength(II)
		END IF

	        !!Get Chem Step (in s)
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real chemistry time step.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetChemStep(II)
		END IF

	        !!Get Mix Step (in s)
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real mixing time step.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetMixStep(II)
		END IF

	        !!Get Cond Step (in s)
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real condensation time step.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetCondStep(II)
		END IF

	        !!Get Coag Step (in s)
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real coagulation time step.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetCoagStep(II)
		END IF

	        !!Get Inversion height for deposition calcs
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real inversion height.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetInversionHeight(II)
		END IF

	        !!Include Aerosols?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if aerosols are to be included. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if aerosols are to be included. Received "//StripToken(Token))	             
                CALL SetIncludeAerosols(I)

	        !!Get Hetero propoerties for bulk or gas only runs
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real number conc.")
		ELSE
			II = STR2REAL(Token)
                        Num = II
		END IF

                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real aerosol radii.")
		ELSE
			II = STR2REAL(Token)
                        CALL SetHetero(Num, II)
		END IF

	        !!Get Lagranigan Mixing param (s)
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real Lagrangian mixing scale.")
		ELSE
			II = STR2REAL(Token)
                        Call SetLagrangianMix(II)
		END IF

	        !!Get Eulerian Box Parameters
                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real wind speed.")
		ELSE
			II = STR2REAL(Token)
                        wind = II
		END IF

                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real box length.")
		ELSE
			II = STR2REAL(Token)
                        length = II
		END IF

                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real box width.")
		ELSE
			II = STR2REAL(Token)
                        width = II

		END IF

                Token = GetLine(FH)

		!! If not a real, give error
                ! CMB (AER): Added not
		if (.not. isReal(token)) then
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real CO Emission Rate.")
		ELSE
			II = STR2REAL(Token)
                        
                        Call SetEulerianBox(wind,length,width, II)
		END IF

                !!Add in lightswitch flag for chamber esperiments
                Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
	                CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if ASP is to be run for a lights on-off chamber experiment. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
			"if ASP is to be run a lights on-off chamber experiment. Received "//StripToken(Token))

                lightswitch_flag = I

                Token = GetLine(FH)
		!! If not a real, give error
		IF (.NOT.IsReal(Token)) THEN
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real time to turn on the lights.")
		ELSE
			II = STR2REAL(Token)
                        lights_on_time = II

		END IF

                Token = GetLine(FH)
		!! If not a real, give error
		IF (.NOT.IsReal(Token)) THEN
			CALL ERROR("In ASPInputDeck.in, You tried to load a non-real time to turn off the lights.")
		ELSE
			II = STR2REAL(Token)
                        lights_off_time = II

		END IF

                CALL LightSwitchParams(lightswitch_flag, lights_on_time, lights_off_time)

		!! Optional solver and evaluation settings at the end of the deck
		CALL ReadOptionalDeckSettings ()
		CALL EchoDeckSettings ()

                !print *, 'CMB: END MAIN INPUT DECK'

	CONTAINS

	!! Optional settings after the last required line of ASPInputDeck.in, read until the end
	!! of the file.  Each is a line
	!!     Name = value        (any order; names are not case sensitive)
	!! Older decks may instead give the first four as positional lines, in this order:
	!!     GDTableMode, GDTableCells, AqInnerSolver, AqOuterSolver
	!! Positional lines must come before any "Name = value" line.  Settings not given keep
	!! their defaults (ModelParameters.f90):
	!!   GDTableMode   1    1 = tabulated Gibbs-Duhem integral, 0 = full DGAUS8 integral
	!!   GDTableCells  256  table cells per reaction (16 to 65536)
	!!   AqInnerSolver 1    single-reaction solver: 1 = Newton, 0 = Jacobson bisection
	!!   AqOuterSolver 1    electrolyte equilibrium: 1 = coupled Newton solve, 0 = worst-first loop
	!!   AqPairRule    -1   solid-circuit precipitation rule: 1 = on, 0 = off,
	!!                      -1 = on with the coupled solve (AqOuterSolver = 1) only
	!!   AqExitFix     0    worst-first loop: 0 = early exit on the last reaction's deviation,
	!!                      1 = no early exit, 2 = also check every reaction at the end
	!!   AqDevCache    1    worst-first loop deviation cache: 1 = on, 0 = off (same results)
	!!   GasEvalLinked 0    gas ODE/Jacobian evaluation: 0 = flat arrays, 1 = linked lists (same results)
	!!   ThermoScanUp  0    ThermoTest RH scan: 0 = down 0.90 -> 0.20, 1 = up 0.20 -> 0.90
	SUBROUTINE ReadOptionalDeckSettings ()
		IMPLICIT NONE
		INTEGER :: NPositional, EqPos
		LOGICAL :: Named
		CHARACTER(len=32), PARAMETER :: PositionalName(4) = (/ "GDTableMode  ", "GDTableCells ", &
		                                                       "AqInnerSolver", "AqOuterSolver" /)

		NPositional = 0
		Named = .FALSE.
		DO
			CurrLine = GetLine(FH)
			IF (CurrLine .EQ. EOF) EXIT
			EqPos = INDEX(CurrLine, "=")
			IF (EqPos .GT. 0) THEN
				Named = .TRUE.
				CALL SetDeckSetting (CurrLine(1:EqPos-1), CurrLine(EqPos+1:))
			ELSE
				IF (Named) CALL ERROR ("In ASPInputDeck.in, positional setting "//TRIM(CurrLine)// &
					" after a 'Name = value' line; use 'Name = value' for all settings after the first one")
				NPositional = NPositional + 1
				IF (NPositional .GT. 4) CALL ERROR ("In ASPInputDeck.in, more than four positional settings "// &
					"at the end of the deck (GDTableMode, GDTableCells, AqInnerSolver, AqOuterSolver); "// &
					"use 'Name = value' lines for the others.  Received "//TRIM(CurrLine))
				CALL SetDeckSetting (PositionalName(NPositional), CurrLine)
			END IF
		END DO
	END SUBROUTINE ReadOptionalDeckSettings

	!! Check and store one optional setting
	SUBROUTINE SetDeckSetting (InName, InValue)
		IMPLICIT NONE
		CHARACTER*(*) :: InName, InValue
		CHARACTER(len=64)  :: Name, UName
		CHARACTER(len=512) :: Value
		INTEGER :: K, V, VMin, VMax

		Name  = StripToken(InName)
		Value = InValue
		DO K = 1, LEN(Value)		!a deck with DOS line ends leaves a carriage return
			IF (Value(K:K) .EQ. CHAR(13)) Value(K:K) = " "
		END DO
		Value = StripToken(Value)
		UName = Name
		DO K = 1, LEN_TRIM(UName)
			IF (UName(K:K) .GE. "a" .AND. UName(K:K) .LE. "z") UName(K:K) = ACHAR(IACHAR(UName(K:K)) - 32)
		END DO

		SELECT CASE (TRIM(UName))
		CASE ("GDTABLEMODE")   ; VMin = 0  ; VMax = 1
		CASE ("GDTABLECELLS")  ; VMin = 16 ; VMax = 65536
		CASE ("AQINNERSOLVER") ; VMin = 0  ; VMax = 1
		CASE ("AQOUTERSOLVER") ; VMin = 0  ; VMax = 1
		CASE ("AQPAIRRULE")    ; VMin = -1 ; VMax = 1
		CASE ("AQEXITFIX")     ; VMin = 0  ; VMax = 2
		CASE ("AQDEVCACHE")    ; VMin = 0  ; VMax = 1
		CASE ("GASEVALLINKED") ; VMin = 0  ; VMax = 1
		CASE ("THERMOSCANUP")  ; VMin = 0  ; VMax = 1
		CASE DEFAULT
			CALL ERROR ("In ASPInputDeck.in, unknown setting "//TRIM(Name)//".  Known settings: GDTableMode, "// &
				"GDTableCells, AqInnerSolver, AqOuterSolver, AqPairRule, AqExitFix, AqDevCache, "// &
				"GasEvalLinked, ThermoScanUp")
		END SELECT
		IF (.NOT. IsInteger(Value)) CALL ERROR ("In ASPInputDeck.in, "//TRIM(Name)// &
			" must be an integer. Received "//TRIM(Value))
		V = STR2INT(Value)
		IF (V .LT. VMin .OR. V .GT. VMax) CALL ERROR ("In ASPInputDeck.in, "//TRIM(Name)//" must be between "// &
			TRIM(INT2STR(VMin))//" and "//TRIM(INT2STR(VMax))//". Received "//TRIM(Value))

		SELECT CASE (TRIM(UName))
		CASE ("GDTABLEMODE")   ; GDTableMode   = V
		CASE ("GDTABLECELLS")  ; GDTableCells  = V
		CASE ("AQINNERSOLVER") ; AqInnerSolver = V
		CASE ("AQOUTERSOLVER") ; AqOuterSolver = V
		CASE ("AQPAIRRULE")    ; AqPairRule    = V
		CASE ("AQEXITFIX")     ; AqExitFix     = V
		CASE ("AQDEVCACHE")    ; AqDevCache    = V
		CASE ("GASEVALLINKED") ; GasEvalLinked = V
		CASE ("THERMOSCANUP")  ; ThermoScanUp  = V
		END SELECT
	END SUBROUTINE SetDeckSetting

	!! Write the solver and evaluation settings in use to the transcript and standard output
	SUBROUTINE EchoDeckSettings ()
		IMPLICIT NONE
		CALL Transcript ("Solver and evaluation settings (ASPInputDeck.in; defaults where not given):")
		CALL Transcript ("  GDTableMode   = "//TRIM(INT2STR(GDTableMode))//"   (1 = tabulated Gibbs-Duhem integral, 0 = DGAUS8)")
		CALL Transcript ("  GDTableCells  = "//TRIM(INT2STR(GDTableCells)))
		CALL Transcript ("  AqInnerSolver = "//TRIM(INT2STR(AqInnerSolver))//"   (1 = Newton, 0 = Jacobson bisection)")
		CALL Transcript ("  AqOuterSolver = "//TRIM(INT2STR(AqOuterSolver))//"   (1 = coupled Newton solve, 0 = worst-first loop)")
		CALL Transcript ("  AqPairRule    = "//TRIM(INT2STR(AqPairRule))//"  (1 = on, 0 = off, -1 = on with the coupled solve only)")
		CALL Transcript ("  AqExitFix     = "//TRIM(INT2STR(AqExitFix))//"   (worst-first loop exit)")
		CALL Transcript ("  AqDevCache    = "//TRIM(INT2STR(AqDevCache))//"   (worst-first loop deviation cache)")
		CALL Transcript ("  GasEvalLinked = "//TRIM(INT2STR(GasEvalLinked))//"   (0 = flat arrays, 1 = linked lists)")
		CALL Transcript ("  ThermoScanUp  = "//TRIM(INT2STR(ThermoScanUp))//"   (ThermoTest RH scan: 0 = down, 1 = up)")
	END SUBROUTINE EchoDeckSettings

	END SUBROUTINE ReadMainInputDeck
