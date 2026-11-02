!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Thermodynamics.h
!! This file contains many functions relating to the aqueous phase
!! thermodynamics of the inorganic ions and salts, including the associated
!! water content. The internal inorganic equilibrium routine is also contained
!! here.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						             !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History	           	     !!
!! 07/19  2007   Matt Alvarado     Added check to 
!!                                  EquilibrateAqueousDissociationReaction 
!!				    to return if (NH4)3H(SO4)2 not in eq 
!!                                  but SO4-- concentration is approaching 0.
!!				   Added check to 
!!                                  EquilibrateAqueousDissociationReaction 
!!				    to return if H2O not in eq but OH-
!!			            concentration is approaching 0 after 
!!                                  25,000 iterations.	                     !!
!! 08/01  2007    Matt Alvarado     Removed low water check from 
!!				FindElectrolyteEquilibrium, since it sometimes 
!!				gave unphysical answers.		     !!
!! 08/02  2007    Matt Alvarado     Removing low water check causes problems. 
!!                                  Put it back in with a water activity check 
!!                                  that seems to do the right thing.  
!! 09/26  2007    Matt Alvarado     Added check to 
!!                                  EquilibrateAqueousDissociationReaction 
!!				to return if (NH4)2SO4 not in eq but SO4--
!!				concentration is approaching 0.
!! 10/05  2007    Matt Alvarado     Set FindElectrolyte equilibrium to just 
!!                                  exit after max iterations
!!				    rather than error.
!!				    Set EquilibrateAqueousDissociationReaction 
!!                                  to just exit after max iterations rather 
!!                                  than error.
!! 10/09  2007    Matt Alvarado     Set MeanActivity to reset any negative 
!!                                  molalities to 0.
!! 10/14  2010    Matt Alvarado     Removed OPTIONAL type declarations
!! 02/16  2012    Matt Alvarado     Removed Eulerian grids, making ASP       !!
!!                                 a one-box model or subroutine.            !!
!! 11/08  2012    Matt Alvarado     Removed some activity coefficient updates!! 
!!                                  in FindElectrolyteEquilibrium to speed   !!
!!                                  it up                                    !!
!! 01/29  2013    Matt Alvarado     Removed EquilibriumWaterContent, replaced!!
!!                                     with EquilibriumWaterContentAmount    !!
!! 10/28  2026   Claude Sonnet 5.5 Commit bc80201: N1-N5: remove redundant   !!
!!                                 setup and equilibrium work (StripToken,   !!
!!                                 MakeJacobian, MeanActivity, solver loop)  !!
!! 10/28  2026   Claude Sonnet 5.5 Commit a56950d: N6-N11: integrand caches, !!
!!                                 EqCoeff cache, tabulated Gibbs-Duhem      !!
!!                                 integral (N11 changes results ~2e-4)      !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 8766c40: Patch 17: guard          !!
!!                                 tabulated Gibbs-Duhem integral           !!
!!                                 against NaN or negative ionic            !!
!!                                 strength (DGAUS8 fallback)               !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 46b4665: Patch 18 (CondTest       !!
!!                                 C1): per-reaction chemical-name flags    !!
!!                                 replace string compares in the           !!
!!                                 dissociation solver                      !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 86d3a78: Patch 19 (CondTest       !!
!!                                 C2): memoise repeated x**y (PowMemo)     !!
!!                                 in MeanActivity and                      !!
!!                                 EquilibriumConstantsRatio                !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 189da4e: N11 input-deck switch: !!
!!                                 GDTableMode (1 = table, 0 = full      !!
!!                                 DGAUS8 integral) and GDTableCells,    !!
!!                                 optional last lines of ASPInputDeck.in!!
!! 10/28  2026   Claude Sonnet 5.5 Commit f08134b: Patch 24: reuse the      !!
!!                                 EquilibriumConstantsRatio value at     !!
!!                                 the entry of EquilibrateAqueous        !!
!!                                 DissociationReaction                   !!
!! 10/29  2026   Claude Sonnet 5.5 Commit 9e33a67: Patch 25: hoist the   !!
!!                                 per-reaction AqEquilibriaList reads    !!
!!                                 and index arithmetic out of the        !!
!!                                 EquilibrateAqueousDissociationReaction !!
!!                                 iteration loop                         !!
!! 11/02  2026   Claude Opus 5.5   Patch 26: by-value memo of             !!
!!                                 EquilibriumConstantsRatio (bit-for-bit !!
!!                                 key; byte-identical output)            !!
!! 11/02  2026   Claude Opus 5.5   Patch 27: FindElectrolyteEquilibrium   !!
!!                                 test harness: ASP_FEE_DUMP writes      !!
!!                                 sampled solver states, ASP_FEE_REPLAY  !!
!!                                 re-solves them; the solver is now      !!
!!                                 FindElectrolyteEquilibriumCore behind a!!
!!                                 wrapper                                !!
!! 11/02  2026   Claude Opus 5.5   Patch 28: Newton single-reaction solve !!
!!                                 in                                     !!
!!                                 EquilibrateAqueousDissociationReaction !!
!!                                 (safeguarded, inside the mass-flux     !!
!!                                 bracket); selected by AqInnerSolver,   !!
!!                                 default Newton                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		     !!
!! 1. FUNCTION PrintAqEqReaction (ReactIndex)                                !!
!! 2. SUBROUTINE PrintParticleContents                                       !!
!!     (InParticle, FH, Initialize, ToScreen, CommentString)                 !!
!! 3. SUBROUTINE DumpParticleContentsAtError (InParticle, OutFileName, &     !!
!!              InRelativeHumidity, InRadius, InpH, InDoWaterActivity, &     !!
!!		InDensity, InDoSurfaceTension, InIonicStrength)              !!
!! 4. FUNCTION EqCoeff (WhichEquilibrium, InParticle , Phase)                !!
!! 5. FUNCTION DeliquesenceRH (WhichEquilibrium, InParticle)                 !!
!! 6. SUBROUTINE FindAqElectrolyteEquilibriumForGridPoint                    !!
!!      (UpdateThermo)                                                       !!
!! 7. SUBROUTINE FindElectrolyteEquilibrium (InParticle, UpdateThermo, &     !!
!!       FirstEquilibration, ReturnType)                                     !!
!! 8. RECURSIVE SUBROUTINE EquilibrateAqueousDissociationReaction            !!
!!       (WhichEquilibriumRxn, InParticle, OptErrorTolerance, ReturnType)    !!
!! 9. SUBROUTINE ReformElectrolyte(WhichRxn, InParticle)                     !!
!!10. FUNCTION CheckEquilibriumConvergenceForAqDissociationRxn               !!
!!      (WhichEquilibriumRxn, InParticle, ErrorTolerance)                    !!
!!11. FUNCTION EquilibriumConstantsRatio (WhichEquilibriumRxn, InParticle)   !!
!!12. SUBROUTINE CalculateIonicStrength (InParticle)                         !!
!!13. FUNCTION MeanActivity (WhichEquilibrium, InParticle)                   !!
!!14. FUNCTION ChargeFraction (CationChg, AnionChg)                          !!
!!15. FUNCTION CalcQ (EqRxn, InParticle)                                     !!
!!16. SUBROUTINE KusikMeissner (InParticle)                                  !!
!!17. FUNCTION CalcLogGammaPure (K, InParticle)                              !!
!!18. SUBROUTINE UpdateWaterActivity (InParticle)                            !!
!!19. FUNCTION WaterResidual (InParticle, SpecifiedRH)                       !!
!!20. SUBROUTINE EquilibriumWaterContent (InParticle, OpenSystem,      &     !!
!!      ReturnType, ForceEquilibrium, EquilibToThisRH)                       !!
!!21. SUBROUTINE EquilibriumWaterContentAmount                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSSC!!
!! Print an equilibrium reaction   !!
!! (for some debugging / reporting !!
!! purpose)						   !!
 CHARACTER (len = 5000) FUNCTION PrintAqEqReaction (ReactIndex)

	USE Chemistry, ONLY : AqPhaseChemicalNames,     &
				AqCationNames,		&
				AqAnionNames,		&
				HowManyAqEqReactions,   &
				AqEquilibriaList

	USE InfrastructuralCode, ONLY : ERROR,REAL2STR,INT2STR

	IMPLICIT NONE

	INTEGER :: ReactIndex

	IF (ReactIndex .GT. HowManyAqEqReactions) &
	CALL ERROR ("Gave PrintAqEqReaction() in Thermodynamics.h an unacceptable equilibrium index.")

	IF (AqEquilibriaList(ReactIndex, 1) .GT. 0) THEN
		PrintAqEqReaction = TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(ReactIndex,1))))
	ELSE 
		PrintAqEqReaction = ""
	END IF

	IF ((AqEquilibriaList(ReactIndex, 1) .GT. 0) .AND. (AqEquilibriaList(ReactIndex, 8) .GT. 0)) &
		PrintAqEqReaction = TRIM(PrintAqEqReaction) // " & "

	IF (AqEquilibriaList(ReactIndex, 8) .GT. 0) THEN
		PrintAqEqReaction = TRIM(PrintAqEqReaction) // " " //TRIM(REAL2STR(AqEquilibriaList(ReactIndex,8),1)) // " H2O + "
	END IF


	PrintAqEqReaction =	TRIM(PrintAqEqReaction)//										&
						" <=> " //														&
						TRIM(REAL2STR(AqEquilibriaList(ReactIndex,4),1)) // " " //		&
						TRIM(AqCationNames(INT(AqEquilibriaList(ReactIndex,2)))) //			&
						" & "//															&
						TRIM(REAL2STR(AqEquilibriaList(ReactIndex,5),1)) // " " //		&
						TRIM(AqAnionNames(INT(AqEquilibriaList(ReactIndex,3))))

	RETURN
END FUNCTION PrintAqEqReaction 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Print the concentrations of the particle to the screen and perhaps a file !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PrintParticleContents (InParticle, FH, Initialize, ToScreen, &
                                  CommentString)

	USE GridPointFields,     ONLY : GetRelativeHumidity
	USE InfrastructuralCode, ONLY : INT2STR,			&
					REAL2STR,			&
					StripToken,			&
				        ERROR,WARN,			&
					Transcript									

	USE Chemistry, ONLY :		HowManyAqChems,		&
					HowManyAqCations,	&
					HowManyAqAnions,	&
					AqPhaseChemicalNames,&
					AqCationNames,		&
					AqAnionNames,		&
					HowManyAqEqReactions,&
					AqEquilibriaList

	USE ModelParameters, ONLY : moles, micron, cm, grams, dyncm, Avogadro

	IMPLICIT NONE

	!! External Variable
	Type(Particle), POINTER :: InParticle
	CHARACTER*(*) :: CommentString
	INTEGER :: FH	! file pointer
	LOGICAL :: Initialize ! if true then do title line, else do values
	LOGICAL :: ToScreen ! if true then dump each result to the screen as well...

	!! Internal Variables
	INTEGER :: I
	INTEGER, PARAMETER :: DecimalPlaces = 8
	LOGICAL :: RelativeHumidity = .TRUE.
	LOGICAL :: Radius = .TRUE.
	LOGICAL :: pH = .TRUE.
	LOGICAL :: WaterActivity = .TRUE.
	LOGICAL :: Density = .TRUE.
	LOGICAL :: DoSurfaceTension = .TRUE.
	LOGICAL :: DoActivities = .TRUE.
	LOGICAL :: DoIonicStrength = .TRUE.
	LOGICAL :: RawWaterMoles = .TRUE.

	!! This should only be called after one time step!!
	!! Initial concentrations for some chemcials that will eventually
	!! be included may be 0, which screws things up.


	!! RH
	IF (RelativeHumidity) THEN
		IF (Initialize)      WRITE (FH, '(a)', advance ='no') "RH, "
		IF (.NOT.Initialize) WRITE (FH, '(a)', advance ='no') & 
                    TRIM(REAL2STR(GetRelativeHumidity(),DecimalPlaces))//", "
		IF (.NOT.Initialize .AND. ToScreen) &
                    PRINT *, "RH                 : ", &
                    TRIM(REAL2STR(GetRelativeHumidity(),DecimalPlaces))
	ENDIF

	!! Radius
	IF (Radius) THEN
		IF (Initialize)      WRITE (FH, '(a)', advance ='no') "Effective Radius, "
		IF (.NOT.Initialize) WRITE (FH, '(a)', advance ='no') &
                    TRIM(REAL2STR(InParticle%Effectiveradius/micron,DecimalPlaces))//", "
		IF (.NOT.Initialize .AND. ToScreen) &
                    PRINT *, "Effective Radius             : ",TRIM(REAL2STR(InParticle%Effectiveradius/micron,DecimalPlaces))
	ENDIF

	!! pH
	IF (pH) THEN
		IF (Initialize)      WRITE (FH, '(a)', advance ='no') "pH, "
		IF (.NOT.Initialize) WRITE (FH, '(a)', advance ='no') &
                    TRIM(REAL2STR(AerosolPH (InParticle),DecimalPlaces))//", "
		IF (.NOT.Initialize .AND. ToScreen) &
                    PRINT *, "pH                 : ",TRIM(REAL2STR(AerosolPH (InParticle),DecimalPlaces))
	ENDIF

	!! pH
	IF (DoIonicStrength) THEN
		IF (Initialize)      WRITE (FH, '(a)', advance ='no') "IonicStr, "
		IF (.NOT.Initialize) WRITE (FH, '(a)', advance ='no') &
                    TRIM(REAL2STR(InParticle%IonicStr,DecimalPlaces))//", "
		IF (.NOT.Initialize .AND. ToScreen) &
                    PRINT *, "Ionic Strength     : ",TRIM(REAL2STR(InParticle%IonicStr,DecimalPlaces))
	ENDIF

	!! Density
	IF (Density) THEN
		IF (Initialize)      WRITE (FH, '(a)', advance ='no') "Density, "
		IF (.NOT.Initialize) WRITE (FH, '(a)', advance ='no') &
                    TRIM(REAL2STR(ParticleDensity (InParticle),DecimalPlaces))//", "
		IF (.NOT.Initialize .AND. ToScreen) &
		PRINT *, "Density             : ",&
                    TRIM(REAL2STR(ParticleDensity (InParticle)*cm*cm*cm/grams,DecimalPlaces))
	ENDIF

	!! Surface Tension
	IF (DoSurfaceTension) THEN
		IF (Initialize)      WRITE (FH, '(a)', advance ='no') "Surface Tension, "
		IF (.NOT.Initialize) WRITE (FH, '(a)', advance ='no') TRIM(REAL2STR(InParticle%SurfaceTension,DecimalPlaces))//", "
		IF (.NOT.Initialize .AND. ToScreen) &
		PRINT *, "Surface Tension     : ",TRIM(REAL2STR(InParticle%SurfaceTension/dyncm,DecimalPlaces))
	ENDIF

	!! Water
	IF(RawWaterMoles) THEN
		IF (Initialize)		  WRITE (FH, '(a)', advance ='no') "Raw H2O, "
		IF (.NOT. Initialize) WRITE (FH, '(a)', advance ='no') TRIM(REAL2STR(InParticle%AqChems(1)*Avogadro,DecimalPlaces))//", "
		IF (.NOT.Initialize .AND. ToScreen) PRINT *, "Raw H2O             : ",TRIM(REAL2STR(InParticle%AqChems(1),DecimalPlaces))
	END IF

	!! Loop Over Chemicals for Molalities
	DO I = 2, HowManyAqChems
		IF (InParticle%AqChems(I) .GT. 0.) THEN
			IF (Initialize)		  WRITE (FH, '(a)', advance ='no') "m("//TRIM(AqPhaseChemicalNames(I))//"), "
			IF (.NOT. Initialize) WRITE (FH, '(a)', advance ='no') TRIM(REAL2STR(Molality(I,InParticle),DecimalPlaces))//", "
			IF (.NOT.Initialize .AND. ToScreen) &
			PRINT *, "m("//AqPhaseChemicalNames(I)//") : ",TRIM(REAL2STR(Molality(I,InParticle),DecimalPlaces))
		ENDIF
	END DO

	DO I = 1, HowManyAqCations
		IF (InParticle%AqChems(I+HowManyAqChems) .GT. 0.) THEN
			IF (Initialize)		  WRITE (FH, '(a)', advance ='no') "m("//TRIM(AqCationNames(I))//"), "
			IF (.NOT. Initialize) WRITE (FH, '(a)', advance ='no') TRIM(REAL2STR(CationMolality(I,InParticle),DecimalPlaces))//", "
			IF (.NOT.Initialize .AND. ToScreen) &
			PRINT *, "m("//AqCationNames(I)//") : ",TRIM(REAL2STR(CationMolality(I,InParticle),DecimalPlaces))
		ENDIF
	END DO

	DO I = 1, HowManyAqAnions
		IF (InParticle%AqChems(I+HowManyAqChems+HowManyAqCations) .GT. 0.) THEN
			IF (Initialize)		  WRITE (FH, '(a)', advance ='no') "m("//TRIM(AqAnionNames(I))//"), "
			IF (.NOT. Initialize) WRITE (FH, '(a)', advance ='no') TRIM(REAL2STR(AnionMolality(I,InParticle),DecimalPlaces))//", "
			IF (.NOT.Initialize .AND. ToScreen) &
			PRINT *, "m("//AqAnionNames(I)//") : ",TRIM(REAL2STR(AnionMolality(I,InParticle),DecimalPlaces))
		ENDIF
	END DO

	!! Loop Over Eq Reactions for Activities
	IF (DoActivities) THEN
	DO I = 1, HowManyAqEqReactions

		IF (Initialize)	THEN
			IF (AqEquilibriaList(I,1) .GT. 0.) THEN
			  WRITE (FH, '(a)', advance ='no') "Gmix("//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))//"), "
			ELSE
			  WRITE (FH, '(a)', advance ='no') "Gmix(H2O), "
		END IF ; END IF
		IF (.NOT. Initialize) WRITE (FH, '(a)', advance ='no') TRIM(REAL2STR(InParticle%GammaMixed(I),DecimalPlaces))//", "
		IF (.NOT.Initialize .AND. ToScreen)	THEN
			IF (AqEquilibriaList(I,1) .GT. 0.) THEN
				PRINT *, "GammaMixed("//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))//") : ",	&
						 TRIM(REAL2STR(InParticle%GammaMixed(I),DecimalPlaces))
			ELSE
				PRINT *, "GammaMixed(H2O) : ",TRIM(REAL2STR(InParticle%GammaMixed(I),DecimalPlaces))
		END IF ; END IF
	END DO
	END IF

	!! Loop Over Eq Reactions for Activities
	!IF (PRESENT(CommentString)) THEN
		WRITE (FH, '(a)', advance ='no') TRIM(CommentString)//", "
		IF (ToScreen) PRINT *, CommentString
	!END IF

	!! Terminate the line
	WRITE (FH, '(a)') ""

	RETURN

END SUBROUTINE PrintParticleContents


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Print the contents of a problematic particle to a file before an error    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DumpParticleContentsAtError (InParticle, OutFileName,    & 
             InRelativeHumidity, InRadius, InpH, InDoWaterActivity, &
	     InDensity, InDoSurfaceTension, InIonicStrength)

	USE GridPointFields,     ONLY : GetRelativeHumidity

	USE InfrastructuralCode, ONLY : INT2STR,			&
					REAL2STR,			&
					StripToken,			&
					Transcript,			&
					WARN,				&
					GetFileHandle,		        &
					ReturnFileHandle

	USE Chemistry,           ONLY : HowManyAqChems, HowManyAqCations,&
                                        HowManyAqAnions,                &
                                        AqPhaseChemicalNames,           &
				        AqCationNames, AqAnionNames,    &
                                        HowManyAqEqReactions,           &
                                        AqEquilibriaList

	USE ModelParameters,     ONLY : moles, micron, cm, grams, dyncm,&
                                        OutPutDeckSubDir

	IMPLICIT NONE

	!! External Variable
	Type(Particle),POINTER    :: InParticle
	CHARACTER*(*)     :: OutFileName
	LOGICAL :: InRelativeHumidity 
	LOGICAL :: InRadius           
	LOGICAL :: InpH               
	LOGICAL :: InDensity          
	LOGICAL :: InDoSurfaceTension 
	LOGICAL :: InIonicStrength
	LOGICAL :: InDoWaterActivity

	!! Internal Variables
	  INTEGER :: I, GG
	INTEGER :: FH	! file pointer
	CHARACTER (len = 8)     :: ErrorDate
	CHARACTER (len = 10)    :: ErrorTime
	INTEGER, PARAMETER      :: DecimalPlaces = 7
	LOGICAL			:: RelativeHumidity 
	LOGICAL			:: Radius           
	LOGICAL			:: pH               
	LOGICAL			:: Density          
	LOGICAL			:: DoSurfaceTension 
	LOGICAL			:: IonicStrength    
	LOGICAL			:: DoWaterActivity 

	FH = GetFileHandle ()
	OPEN(UNIT=FH, FILE=TRIM(OutputDeckSubDir)//TRIM(OutFileName))

	CALL WARN ("Printing details of problematic aerosol to "//TRIM(OutFileName))

	CALL Date_and_Time(ErrorDate, ErrorTime)
	WRITE (FH, '(a)') "*********************************************"
	WRITE (FH, '(a)') "** This file  contains a profile  of the   **"
	WRITE (FH, '(a)') "** particle that lead MELAM to err at      **"
	WRITE (FH, '(a)') "** "//ErrorTime(1:2)//":"//ErrorDate(3:4)//":"//ErrorDate(5:6)//" on "//ErrorDate(5:6)//		&
					  "/"//ErrorDate(7:8)//"/"//ErrorDate(1:4)//" "//"                 "//"**"
	WRITE (FH, '(a)') "*********************************************"

	IF (InParticle%Dry)	  WRITE (FH, '(a)') "Particle is Dry "
	IF (RelativeHumidity) WRITE (FH, '(a)') "RH              = "//&
            TRIM(REAL2STR(GetRelativeHumidity (),DecimalPlaces))
	IF (Radius)           WRITE (FH, '(a)') &
            "Effecitve Radius          = "//&
            TRIM(REAL2STR(InParticle%Effectiveradius/micron,DecimalPlaces))
	IF (pH)               WRITE (FH, '(a)') "pH              = "//&
            TRIM(REAL2STR(AerosolPH (InParticle),DecimalPlaces))
	IF (Density)		  &
	WRITE (FH, '(a)') "Solution Density= "//&
            TRIM(REAL2STR(ParticleDensity (InParticle)*cm*cm*cm/grams,DecimalPlaces))
	IF (Density)		  &
	WRITE (FH, '(a)') "Particle Density= "//&
            TRIM(REAL2STR(InParticle%ParticleDensity*cm*cm*cm/grams,DecimalPlaces))
	IF (DoSurfaceTension) WRITE (FH, '(a)') "Surface Tension = "//&
            TRIM(REAL2STR(InParticle%surfaceTension/dyncm,DecimalPlaces))
	IF (IonicStrength)    WRITE (FH, '(a)') "Ionic Strength  = "//&
            TRIM(REAL2STR(InParticle%IonicStr,DecimalPlaces))
        WRITE (FH, '(a)') "# of Particles  = "//TRIM(REAL2STR(InParticle%NumberOfParticles,DecimalPlaces))
        WRITE (FH, '(a)') "Distribution #  = "//TRIM(INT2STR(InParticle%ParticleDistribution))
	IF (InParticle%Sectional) &
	WRITE (FH, '(a)') "Sectional Particle with Edges "//TRIM(REAL2STR(InParticle%Edges(1),DecimalPlaces))// &
			  " and "//TRIM(REAL2STR(InParticle%Edges(2),DecimalPlaces))
	IF (.NOT.InParticle%Sectional) WRITE (FH, '(a)') "Lagrangian Particle with ID "//TRIM(INT2STR(InParticle%ParticleID))

	IF (DoWaterActivity)  THEN
		CALL UpdateWaterActivity(InParticle)
		WRITE (FH, '(a)') "Water Activity  = "//TRIM(REAL2STR(InParticle%WaterActivity,DecimalPlaces))
	END IF

	WRITE (FH, '(a)') ""
	WRITE (FH, '(a)') "** Moles of Water **"
	WRITE (FH, '(a)') "H20 Content: "//TRIM(REAL2STR(InParticle%AqChems(1)/moles))

	WRITE (FH, '(a)') ""
	WRITE (FH, '(a)') "** Molalities **"
	DO I = 2, HowManyAqChems
		IF (InParticle%AqChems(I) .GT. 0.) &
			WRITE (FH, '(a)') "m("//TRIM(AqPhaseChemicalNames(I))//")  ="//TRIM(REAL2STR(Molality(I,InParticle),DecimalPlaces))
		IF (InParticle%AqChems(I) .LT. 0.) &
			WRITE (FH, '(a)') "m("//TRIM(AqPhaseChemicalNames(I))//")  = (negative concentration)"

	END DO

	DO I = 1, HowManyAqCations
		IF (InParticle%AqChems(I+HowManyAqChems) .GT. 0.) &
			WRITE (FH, '(a)') "m("//TRIM(AqCationNames(I))//") ="//TRIM(REAL2STR(CationMolality(I,InParticle),DecimalPlaces))
		IF (InParticle%AqChems(I+HowManyAqChems) .LT. 0.) &
			WRITE (FH, '(a)') "m("//TRIM(AqCationNames(I))//") = (negative concentration)"
	END DO

	DO I = 1, HowManyAqAnions
		IF (InParticle%AqChems(I+HowManyAqChems+HowManyAqCations) .GT. 0.) &
			WRITE (FH, '(a)') "m("//TRIM(AqAnionNames(I))//") ="//TRIM(REAL2STR(AnionMolality(I,InParticle),DecimalPlaces))
		IF (InParticle%AqChems(I+HowManyAqChems+HowManyAqCations) .LT. 0.) &
			WRITE (FH, '(a)') "m("//TRIM(AqAnionNames(I))//") = (negative concentration)"
	END DO

	WRITE (FH, '(a)') ""
	WRITE (FH, '(a)') "** Activity Coefficients **"
	DO I = 1, HowManyAqEqReactions
		IF (AqEquilibriaList(I,1) .GT. 0) THEN
			WRITE (FH, '(a)') "GammaMixed("//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))//") = "//  &
							  TRIM(REAL2STR(InParticle%GammaMixed(I),DecimalPlaces))
		ELSE
			WRITE (FH, '(a)') "GammaMixed(H2O) = "//TRIM(REAL2STR(InParticle%GammaMixed(I),DecimalPlaces))
		END IF
	END DO

	!! Terminate the line
	WRITE (FH, '(a)') ""

	CLOSE(FH)
	CALL ReturnFileHandle(FH)

	RETURN

END SUBROUTINE DumpParticleContentsAtError

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! An Equilibrium Coefficient is a function of temperature and !!
!! three user defined constants.			       !!
!!							       !!
!! Returns a positive value unless the equilibrium reaction is !!
!! undefined, in which case it returns zero (indicating that   !!
!! no dissociation occurs).				       !!
!!							       !!
!! The Phase Flag is either Aqueous-Phase (1) or Gas-Aerosol   !!
!! Transfer (2).					       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION EqCoeff (WhichEquilibrium, InParticle , Phase)

	USE Chemistry, ONLY : AqEquilibriaList
	USE ModelParameters, ONLY : EqCoeffTref 

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle
	INTEGER :: WhichEquilibrium, Phase

	!! Internal Variables
	REAL*8  :: TempRatio

	!! Shouldn't ever trip this, but may somewhere
	!! These equilibria are infinitely dissociating
	!! (like sulfuric acid might be assumed to be)
	IF (AqEquilibriaList(WhichEquilibrium, 9) .EQ. -1.) THEN
		EqCoeff = -1.
		RETURN
	END IF

	!Patch N7: the value depends only on the reaction and the particle temperature
	!(EqCoeffTref is a PARAMETER), so remember it per reaction and reuse while T is unchanged.
	IF (Phase .EQ. 1) THEN
		IF (.NOT. ALLOCATED(EqCoeffCacheT)) THEN
			ALLOCATE (EqCoeffCacheT(SIZE(AqEquilibriaList,1)), EqCoeffCacheV(SIZE(AqEquilibriaList,1)))
			EqCoeffCacheT = -1.0D300
			EqCoeffCacheV = 0.0D0
		END IF
		IF (InParticle%Temperature .EQ. EqCoeffCacheT(WhichEquilibrium)) THEN
			EqCoeff = EqCoeffCacheV(WhichEquilibrium)
			RETURN
		END IF
	END IF

	TempRatio = EqCoeffTref / InParticle%Temperature

	!! It is an aqueous equilibrium reaction
	IF (Phase .EQ. 1) THEN
	EqCoeff = AqEquilibriaList(WhichEquilibrium, 9)	 &
	  * EXP(AqEquilibriaList(WhichEquilibrium, 10) * (TempRatio - 1.) + &
	  AqEquilibriaList(WhichEquilibrium, 11) * &
          (LOG(TempRatio) - TempRatio + 1.))
		EqCoeffCacheT(WhichEquilibrium) = InParticle%Temperature
		EqCoeffCacheV(WhichEquilibrium) = EqCoeff
	END IF

	RETURN
END FUNCTION EqCoeff 

REAL*8 FUNCTION DeliquesenceRH (WhichEquilibrium, InParticle)

	USE Chemistry, ONLY : AqEquilibriaList

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle
	INTEGER :: WhichEquilibrium

	!! Internal Variables
	REAL*8  :: InvTemperature, InvRefTemp


	InvTemperature = 1.0/InParticle%Temperature
	InvRefTemp = 1.0/298.0 

	DeliquesenceRH = AqEquilibriaList(WhichEquilibrium, 20)*exp(AqEquilibriaList(WhichEquilibrium, 12) &
						*(InvTemperature-InvRefTemp))
	!WRITE(*,*) "DRH :", WhichEquilibrium, DeliquesenceRH

	RETURN
END FUNCTION DeliquesenceRH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Loop over all of the particles at a grid point and !! 
!! find the electrolyte equilibria for each.          !!
SUBROUTINE FindAqElectrolyteEquilibriumForGridPoint (UpdateThermo)

	IMPLICIT NONE

	INTEGER :: GG
	LOGICAL :: UpdateThermo

	TYPE(PARTICLE),POINTER :: Current

	Current => Particles%First

        CALL FindElectrolyteEquilibrium (Current, UpdateThermo, FirstEquilibration=.FALSE., ReturnType = GG)

	DO WHILE (Associated(Current%Next))
		Current => Current%Next
		CALL FindElectrolyteEquilibrium (Current, UpdateThermo, &
                     FirstEquilibration=.FALSE., ReturnType = GG)
	END DO


	RETURN
END SUBROUTINE FindAqElectrolyteEquilibriumForGridPoint 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! --- FIND ELECTROLYTE EQUILIBRIUM ---					    !!
!!									    !!
!! Iterates over various partitionings between ions and electrolytes until  !!
!! it finds an equiliribum solution.					    !!
!! That equilibrium is defined using Kusik-Meissner activities, user-	    !!
!! defined equilibrium coefficients, and an interative solver.		    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FindElectrolyteEquilibriumCore (InParticle, UpdateThermo, &
                                       FirstEquilibration, ReturnType)

	USE Chemistry,	ONLY : HowManyAqEqReactions,			&
				HowManyAqCations,			&
				HowManyAqAnions,			&
				HowManyAqChems,				&
				AqCationCharge,				&
				AqAnionCharge,				&
				AqEquilibriaList 

	USE InfrastructuralCode, ONLY : INT2STR,REAL2STR,		&
					TRANSCRIPT,			&
					ERROR, WARN

	USE ModelParameters,     ONLY : AqThermoEquilibriumError,	&
					AqThermoNumbAllowedIterations,	&
					AerosolWaterEquilibriumRHThreshold,&
					WaterEquilibriumIndex,		&
					SmallestAerosolPossible,        &
					MinimumWater

	USE Time,		ONLY : BeginningTime,			&
				       CurrentTime 
									
	USE GridPointFields,     ONLY : GetRelativeHumidity

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle
	LOGICAL      :: FirstEquilibration
	INTEGER      :: ReturnType
	LOGICAL                :: UpdateThermo


        
	!! Internal Variables
	INTEGER :: I, J, K, WorstEqIndex, DemoteAfterThisIteration
	REAL*8  :: II, JJ, KK, LL, CC, RH, WorstEqValue, AnTot, CatTot, ErrorTol

    ! CMB (AER): Gfortran complains at being able to perform int/bool/double
    !            conversions implicitly; let's provide more context.
    logical :: majorReaction(howManyAqEqReactions), &
               initiallyEquilibrated(HowManyAqEqReactions), &
               localFirstEquilibration, &
               alreadyEquilibrated
    integer :: numbInitIterations, localNumbIterations
!	LOGICAL :: MajorReaction(HowManyAqEqReactions), InitiallyEquilibrated(HowManyAqEqReactions),	&
!			   LocalFirstEquilibration, LocalNumbIterations, NumbInitIterations, AlreadyEquilibrated

	!! If equilibration of a particluar reaction affects
	!! the ionic concentration (cationic or anionic separately) 
	!! by less than this fraction, then it is thought to be 
	!! of secondary importance and is ignored until the end,
	!! at which point it is equilibrated once.
	REAL*8, PARAMETER :: CritMajMin = AqThermoEquilibriumError / 100.

        ! cb add
        real*8, parameter :: zero_thresh = 1.1E-40	! TODO
	
        !! Use Scaffolding Code?
	LOGICAL :: Scaffolding
        
        Scaffolding = .FALSE. ! cb
        !Scaffolding = .true. ! cb in for debugging
        cattot = 0.0D0
	IF(Scaffolding) WRITE(*,*) "Entering FindElectrolyteEquilibrium"
	
	!! If this is an empty section, we have no truck with it
        ! Modified by CB: This doesn't work correctly since it's performing
        !                 equality of a float/double with an int.  Let's try
        !                 this modification... 
		!
		! Yep, have to keep this in...otherwise I
	!IF (InParticle%NumberOfParticles .EQ. 0) then
        if (inparticle%numberofparticles .le. zero_thresh) then
            !print *, 'CB: Leaving FindElectrolyteEquilibrium since no particles'
            RETURN
        else
            !print *, '		CB: numparticles = ', inparticle%numberofparticles
        endif



	IF(Scaffolding) WRITE(*,*) "Me 1"

	! CMB (AER, Inc): Revise to provide more error description
	if (inparticle%aqchems(1) .le. zero_thresh) then
		call warn("FindElectrolyteEquilibrium() was called for a particle that "// &
				"contains no water.  InParticle%AqChems(1) = "// &
				trim(real2str(InParticle%AqChems(1)))// &
				".  Zero thresh is currently 1.0e-40...")
		return 
	endif					
				
!	IF (InParticle%AqChems(1) .LE. zero_thresh) & !1.0E-17) &   ! CMB mod for zero equality
!			CALL ERROR("FindElectrolyteEquilibrium() was called for a particle that "// &
!					   "contains no water.  This isn't ok...")

	!! Check optional argument
	! CMB: ?
	!IF (.NOT.PRESENT(FirstEquilibration)) THEN
	!	LocalFirstEquilibration = .FALSE.
	!ELSE
		LocalFirstEquilibration = FirstEquilibration
	!END IF
	
	!write(*,*) 'LocalFirstEquilibration: ', LocalFirstEquilibration
	!Skip if water content is very small and this isn't an initialization step
	IF(.NOT.(LocalFirstEquilibration)) THEN
	CC = CurvatureCorrection (InParticle)
	RH = GetRelativeHumidity ()
        ! CMB (AER): Comment this out
	!WRITE(*,*) "Water Activity: ", InParticle%WaterActivity
	IF(InParticle%AqChems(1) .LT. MinimumWater .AND. .NOT.(LocalFirstEquilibration) &
		.AND. InParticle%WaterActivity*CC .LT. RH) THEN
		ReturnType = 2
	
		RETURN
	END IF
	END IF


	IF(Scaffolding) WRITE(*,*) "Before Kusik", LocalFirstEquilibration
	!! Update the Mixed Activity Coefficients before get started unless not doing that
	IF (UpdateThermo) CALL KusikMeissner (InParticle)
	IF(Scaffolding) WRITE(*,*) "After Kusik"

	IF (.NOT.LocalFirstEquilibration) THEN
		LocalNumbIterations = AqThermoNumbAllowedIterations
		ErrorTol = AqThermoEquilibriumError/2.
		IF (InParticle%EmbryoRadius .GE. SmallestAerosolPossible) THEN
			NumbInitIterations      = 3.
		ELSE
			NumbInitIterations      = 3.
		END IF
	ELSE
		LocalNumbIterations     = AqThermoNumbAllowedIterations * 4.
		NumbInitIterations      = 30.
		ErrorTol = AqThermoEquilibriumError/10.
	END IF

	DemoteAfterThisIteration = 10 !CEILING (LocalNumbIterations / 8.)

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Equilibrate each reaction before we selectively iterate       !!	
	!! This, at least, will fix all of the infinitely dissociating   !!
	!! reactions  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF(Scaffolding) WRITE(*,*) "Before First loop"
	AlreadyEquilibrated = .TRUE.
	DO I = 1, HowManyAqEqReactions

		!WRITE(*,*) J
		CALL EquilibrateAqueousDissociationReaction (I, InParticle, OptErrorTolerance = ErrorTol, ReturnType=J)
	        IF(Scaffolding) WRITE(*,*) I,J
		IF (J .EQ. 1) THEN
			AlreadyEquilibrated      = .FALSE.
			InitiallyEquilibrated(I) = .FALSE.
	 	ELSE
			InitiallyEquilibrated(I) = .TRUE.
		END IF
		!MJA 11-08-2012
                !IF (UpdateThermo) CALL KusikMeissner (InParticle)
	END DO

	IF (Scaffolding) WRITE(*,*) "First Loop Okay.", AlreadyEquilibrated

	!ReturnType ususally not present
	!IF (PRESENT(ReturnType)) THEN
		IF (AlreadyEquilibrated .OR. InParticle%EmbryoRadius .LT. SmallestAerosolPossible) THEN
			ReturnType = 2
		ELSE
			ReturnType = 1
		END IF
	!END IF

	IF (AlreadyEquilibrated) RETURN


	IF (Scaffolding) WRITE(*,*) "Before Second Loop."
	
	!Here, iterate about 30 times, all reactions
	DO K = 1, NumbInitIterations-1
		DO I = 1, HowManyAqEqReactions
			!IF (Scaffolding) WRITE(*,*) I,K
			CALL EquilibrateAqueousDissociationReaction (I, InParticle, OptErrorTolerance = ErrorTol, ReturnType = J)
			IF (Scaffolding) WRITE(*,*) "Dissociation Okay"
                        !MJA 11-08-2012
			!IF (UpdateThermo) CALL KusikMeissner (InParticle)
			IF (Scaffolding) WRITE(*,*) "KusikMeissner Okay"
		END DO
	END DO

	IF (Scaffolding) WRITE(*,*) "Second Loop Okay."

	!! Count iterations
	K = 0

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! We begin by assuming that each of the     !!
	!! reactions is essential to the equilibrium.!!
	!! BUT, if it turns out that one of them is  !!
	!! diddling around the margin, then we'll    !!
	!! demote it and equilibrate it at the end   !!
25	DO I = 1, HowManyAqEqReactions  !!!!!!!!!!!!!!!
		MajorReaction(I) = .TRUE.
	END DO

	!! Loop until equilibrated
	DO

	K = K+1  !! tally the iterations
	IF(SCAFFOLDING) WRITE(*,*) "NumIter", K
	
	
	!! See if any reactions are out of equilibrium
	!! and keep track of the reaction farthest from equilibrium
	WorstEqIndex = 0  ;  WorstEqValue = AqThermoEquilibriumError
	DO I = 1, HowManyAqEqReactions

		!! Skip the minor reactions in this process
		IF (.NOT. MajorReaction(I)) THEN
                        ! CMB (AER, Inc): Comment this out since it 
                        !                 comes out all the time
			!WRITE(*,*) "Minor Reaction", I
			CYCLE
		END IF

		!! Get the deviation from Equilibrium
		! CMB: use double precision
		!II = ABS(1.-EquilibriumConstantsRatio (I, InParticle))
		II = dabs(1. - EquilibriumConstantsRatio(I, InParticle))

		!! If it is the worst the loop has seen, keep track of it
		IF (II .GT. WorstEqValue) THEN
			WorstEqIndex = I
			WorstEqValue = II
		END IF
	END DO

        ! CMB (AER, Inc): this comes up often; let's comment it out
	!WRITE(*,*) "Worst: ", WorstEqIndex, WorstEqValue

	!If Worst equilibrium is below error tolerance, we're done
	IF (II .LT. AqThermoEquilibriumError .AND. K .GT. 5*DemoteAfterThisIteration) THEN
		EXIT
	END IF

	IF(SCAFFOLDING) WRITE(*,*) "Check 1"

	!! Just exit if eq. not reached (MJA, 100507)
	!! If this is the first equilibration, it may be much harder 
	!! to equilibrate so give all the time it needs
	IF (K .GT. LocalNumbIterations) THEN
		EXIT
	END IF
	
	IF(SCAFFOLDING) WRITE(*,*) "Check 2"
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! If done iterating the major reactions, then			!!
	!! equilibrate each of the minor ones once and return   !!
	IF (WorstEqIndex .EQ. 0) THEN  !!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! Equilibrate each of the minor reactions once
		!! Ordering here matters.  If any of the reactions are
		!! conflicting with the water equilibration reaction, it is
		!! very useful to do that one first (since some, like 
		!! Ammonia<=>Ammonium, are often written as compound reactions
		!! with the water reaction)
		IF (WaterEquilibriumIndex .NE. 0 .AND. .NOT.MajorReaction(WaterEquilibriumIndex)) &
			CALL EquilibrateAqueousDissociationReaction(WaterEquilibriumIndex, InParticle, &
                                OptErrorTolerance =AqThermoEquilibriumError, ReturnType=J)

		DO I = 1, HowManyAqEqReactions
			IF (.NOT. MajorReaction(I) .AND. I.NE.WaterEquilibriumIndex) &
                            CALL EquilibrateAqueousDissociationReaction (I,InParticle, &
                                    OptErrorTolerance =AqThermoEquilibriumError, ReturnType=J)
		END DO

		!! Make sure this doesn't mess up the major reactions
		DO I = 1, HowManyAqEqReactions

		IF (MajorReaction(I)) THEN
			IF (ABS(1.-EquilibriumConstantsRatio (I, InParticle)) .GT. AqThermoEquilibriumError) &
				GOTO 25 !Go to top of main loop and try again
		END IF
		END DO

		EXIT
	END IF

	IF(SCAFFOLDING) WRITE(*,*) "Check 3"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Capture total molality of ions and return if in equilibrium. !!
	!! This total is used for the entire equilibration, which       !!
	!! assumes that it is close enough to equilibrium to run this   !!
	!! sort of a filter. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (K .EQ. DemoteAfterThisIteration) THEN
	CatTot = 0.  ;  AnTot = 0. 
	DO I = 1, HowManyAqEqReactions

			CALL EquilibrateAqueousDissociationReaction (I, InParticle, OptErrorTolerance =AqThermoEquilibriumError, ReturnType=J)

		II = 0. ; JJ = 0.
		DO J = 1, HowManyAqCations
			II = II + InParticle%AqChems(J+HowManyAqChems)
		END DO

		DO J = 1, HowManyAqAnions
			JJ = JJ + InParticle%AqChems(J+HowManyAqChems+HowManyAqCations)
		END DO

		!! Get the Total Ionic Molality
		IF (II .GT. CatTot) CatTot = II
		IF (JJ .GT. AnTot)  AnTot  = JJ

	END DO
	END IF

	!IF(SCAFFOLDING) WRITE(*,*) "Check 4, Furthest form EQ: ", WorstEqIndex, WorstEqValue 

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! EQUILIBRATE THE WORST REACTION !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Save the current molalities so we can 
	!! see if this turned out to be a minor
	!! reaction in the overall system
	II = InParticle%AqChems(INT(AqEquilibriaList(WorstEqIndex,2))+HowManyAqChems)
	JJ = InParticle%AqChems(INT(AqEquilibriaList(WorstEqIndex,3))+HowManyAqChems+HowManyAqCations)
	IF(AqEquilibriaList(WorstEqIndex,22) .NE. 0.) &
            LL = InParticle%AqChems(INT(AqEquilibriaList(WorstEqIndex,22))+HowManyAqChems+HowManyAqCations)
	!! Equilibrate that worst reaction
	CALL EquilibrateAqueousDissociationReaction (WorstEqIndex, InParticle, &
                OptErrorTolerance =AqThermoEquilibriumError, ReturnType=J)

	!! Now evaluate whether this equilibration changed anything enough
	!! to continue working this as part of the major equilibrium problem.
	!! IF NOT -- Then put it in reserve and equilibrate it once again at the end
	IF(AqEquilibriaList(WorstEqIndex,22) .EQ. 0.) THEN
		!Normal Case
		IF ( MAX(ABS((II - InParticle%AqChems(INT(AqEquilibriaList(WorstEqIndex,2))+HowManyAqChems)) / CatTot),	&
			 ABS((JJ - InParticle%AqChems(INT(AqEquilibriaList(WorstEqIndex,3))+&
                                     HowManyAqChems+HowManyAqCations)) / AnTot))  &
			 .LT. CritMajMin  .AND.	K .GE. DemoteAfterThisIteration) &
			MajorReaction(WorstEqIndex) = .FALSE.
	ELSE
		!Levitocite Case
		IF ( MAX(ABS((II - InParticle%AqChems(INT(AqEquilibriaList(WorstEqIndex,2))+&
                                        HowManyAqChems)) / CatTot),					&
			 ABS((JJ - InParticle%AqChems(INT(AqEquilibriaList(WorstEqIndex,3))+&
                                     HowManyAqChems+HowManyAqCations)) / AnTot), &
			 ABS((LL - InParticle%AqChems(INT(AqEquilibriaList(WorstEqIndex,22))+&
                                     HowManyAqChems+HowManyAqCations)) / AnTot))  &
			 .LT. CritMajMin  .AND.	K .GE. DemoteAfterThisIteration) &
			MajorReaction(WorstEqIndex) = .FALSE.
	END IF

	IF(SCAFFOLDING) WRITE(*,*) "Check 5"

	END DO  ! end the main loop

	!! If only minor reactions were initially out of equilibrium,
	!! that's as good as the entire system being in equilibrium.
	!IF (PRESENT(ReturnType)) THEN
		ReturnType = 2
		DO J = 1, HowManyAqEqReactions
			IF (MajorReaction(J) .AND. .NOT.InitiallyEquilibrated(J)) THEN
                                ! CMB (AER, Inc): Commented the write statement out
				!WRITE(*,*) J, InitiallyEquilibrated(J)
				ReturnType = 1
			END IF
		END DO
	!END IF

	IF (UpdateThermo) THEN
		!! Update the Mixed Activity Coefficients after each equilibration
		CALL KusikMeissner (InParticle)
		CALL RecalculateRadius(InParticle)
	END IF


	!! Report total number of iterations
	IF (Scaffolding) CALL TRANSCRIPT ("")
	IF (Scaffolding) CALL TRANSCRIPT ("Called "//&
                TRIM(INT2STR(J))//" equilibrations of individual reactions before "// &
		"equilibrium for the system was achieved")

	RETURN
END SUBROUTINE FindElectrolyteEquilibriumCore

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! --- EQUILIBRATE ONE AQUOEOUS DISSOCIATION REACTION ---      !!
!!			  				       !!
!! Iterate a particular reaction until it reaches equilibrium. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Patch C2: x**y with a memo of the last (x,y) per (slot, reaction).  The result is exactly what
!pow would return for the same operands, so output is bit-identical; it only skips repeated
!identical calls (activity coefficients and water activity are constant across solver iterations).
!Only strictly positive finite bases are memoised.

!Patch 27: thin wrapper.  With no ASP_FEE_DUMP set it only calls the unchanged solver.
SUBROUTINE FindElectrolyteEquilibrium (InParticle, UpdateThermo, &
                                       FirstEquilibration, ReturnType)
	IMPLICIT NONE
	TYPE(PARTICLE),POINTER :: InParticle
	LOGICAL :: UpdateThermo, FirstEquilibration
	INTEGER :: ReturnType
	LOGICAL :: Dump

	IF (.NOT. AqSwitchesRead) CALL ReadAqSwitches ()
	Dump = .FALSE.
	IF (FEEDumpOn .AND. .NOT.UpdateThermo .AND. .NOT.FirstEquilibration) THEN
		FEEDumpCount = FEEDumpCount + 1
		Dump = (MOD(FEEDumpCount, INT(FEEDumpEvery,8)) .EQ. 0)
	END IF
	IF (Dump) CALL WriteFEEState (InParticle, UpdateThermo, FirstEquilibration, ReturnType, .TRUE.)
	CALL FindElectrolyteEquilibriumCore (InParticle, UpdateThermo, FirstEquilibration, ReturnType)
	IF (Dump) CALL WriteFEEState (InParticle, UpdateThermo, FirstEquilibration, ReturnType, .FALSE.)
END SUBROUTINE FindElectrolyteEquilibrium

SUBROUTINE ReadAqSwitches ()
	USE ModelParameters, ONLY : AqInnerSolver
	IMPLICIT NONE
	CHARACTER(LEN=512) :: V
	INTEGER :: L, ST
	AqSwitchesRead = .TRUE.
	AqNewtonInner = (AqInnerSolver .EQ. 1)
	CALL GET_ENVIRONMENT_VARIABLE ("ASP_AQ_INNER", V, L, ST)
	IF (ST .EQ. 0 .AND. L .GT. 0) THEN
		IF (TRIM(V) .EQ. "newton")   AqNewtonInner = .TRUE.
		IF (TRIM(V) .EQ. "jacobson") AqNewtonInner = .FALSE.
	END IF
	CALL GET_ENVIRONMENT_VARIABLE ("ASP_FEE_DUMP_EVERY", V, L, ST)
	IF (ST .EQ. 0 .AND. L .GT. 0) READ (V, *) FEEDumpEvery
	CALL GET_ENVIRONMENT_VARIABLE ("ASP_FEE_DUMP", V, L, ST)
	IF (ST .EQ. 0 .AND. L .GT. 0) THEN
		FEEDumpOn = .TRUE.
		OPEN (NEWUNIT=FEEDumpUnit, FILE=TRIM(V), FORM="UNFORMATTED", ACCESS="STREAM", STATUS="REPLACE")
	END IF
END SUBROUTINE ReadAqSwitches

!Patch 27: one record = everything FindElectrolyteEquilibriumCore reads (Pre) or writes (Post).
SUBROUTINE WriteFEEState (P, UpdateThermo, FirstEquilibration, ReturnType, Pre)
	USE GridPointFields, ONLY : GetTemp, GetPress, GetGridCellChemBurden
	IMPLICIT NONE
	TYPE(PARTICLE),POINTER :: P
	LOGICAL :: UpdateThermo, FirstEquilibration, Pre
	INTEGER :: ReturnType, IU, IF1
	IU = FEEDumpUnit
	IF (Pre) THEN
		IU = FEEDumpUnit
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

!Patch 27: replay a dump through FindElectrolyteEquilibriumCore with the solver selected now,
!time only the solver calls, and compare each result with the dumped result.
SUBROUTINE MaybeReplayFEE ()
	USE GridPointFields, ONLY : GetTemp, GetPress, SetTempField, SetPressField, &
	                            GetGridCellChemBurden, ReplaceGridCellChemBurden
	USE ModelParameters, ONLY : PressureScale
	IMPLICIT NONE
	CHARACTER(LEN=512) :: V
	INTEGER :: L, ST, IU, Kind, Flags, RT, RTref, IOS, N, NSame, NRT, I, OU
	TYPE(PARTICLE),POINTER :: P, Ref
	REAL*8 :: T, Pr, H2O, Dum1, Dum2, Dum3, C0, C1, TSolve, D, DMax, DSum, Scale
	REAL*8, ALLOCATABLE :: Out(:)

	CALL GET_ENVIRONMENT_VARIABLE ("ASP_FEE_REPLAY", V, L, ST)
	IF (ST .NE. 0 .OR. L .EQ. 0) RETURN
	IF (.NOT. AqSwitchesRead) CALL ReadAqSwitches ()
	P => Particles%First
	ALLOCATE (Ref)
	ALLOCATE (Ref%AqChems(SIZE(P%AqChems)), Ref%GammaMixed(SIZE(P%GammaMixed)), &
	          Ref%OrgChems(SIZE(P%OrgChems)), Ref%AqOrgChems(SIZE(P%AqOrgChems)))
	ALLOCATE (Out(SIZE(P%AqChems)))
	OPEN (NEWUNIT=IU, FILE=TRIM(V), FORM="UNFORMATTED", ACCESS="STREAM", STATUS="OLD")
	N = 0 ; NSame = 0 ; NRT = 0 ; TSolve = 0. ; DMax = 0. ; DSum = 0.
	DO
		CALL ReadFEEState (IU, P, Kind, Flags, RT, T, Pr, H2O, IOS)
		IF (IOS .NE. 0) EXIT
		IF (Kind .NE. 1) STOP "ReplayFEE: record out of order"
		CALL SetTempField (T)
		IF (GetPress() .NE. Pr) CALL SetPressField (Pr / PressureScale)
		IF (GetPress() .NE. Pr) STOP "ReplayFEE: cannot restore pressure exactly"
		CALL ReplaceGridCellChemBurden (1, H2O)
		CALL CPU_TIME (C0)
		CALL FindElectrolyteEquilibriumCore (P, MOD(Flags,2) .EQ. 1, Flags .GE. 2, RT)
		CALL CPU_TIME (C1)
		TSolve = TSolve + (C1 - C0)
		Out = P%AqChems
		CALL ReadFEEState (IU, Ref, Kind, Flags, RTref, Dum1, Dum2, Dum3, IOS)
		IF (IOS .NE. 0 .OR. Kind .NE. 2) STOP "ReplayFEE: missing result record"
		N = N + 1
		IF (ALL(TRANSFER(Out,0_8,SIZE(Out)) .EQ. TRANSFER(Ref%AqChems,0_8,SIZE(Out)))) NSame = NSame + 1
		IF (RT .EQ. RTref) NRT = NRT + 1
		!difference of each amount, relative to the particle's total dissolved moles
		Scale = SUM(ABS(Ref%AqChems(2:))) + 1.0D-300
		D = MAXVAL(ABS(Out(2:) - Ref%AqChems(2:))) / Scale
		DMax = MAX(DMax, D) ; DSum = DSum + D
	END DO
	CLOSE (IU)
	OPEN (NEWUNIT=OU, FILE="FEEReplay.txt", STATUS="REPLACE")
	WRITE (OU, '(A,I10)')     "records          ", N
	WRITE (OU, '(A,I10)')     "bitwise same     ", NSame
	WRITE (OU, '(A,I10)')     "same ReturnType  ", NRT
	WRITE (OU, '(A,ES12.4)')  "solver CPU s     ", TSolve
	WRITE (OU, '(A,ES12.4)')  "max  |d|/total   ", DMax
	WRITE (OU, '(A,ES12.4)')  "mean |d|/total   ", DSum / MAX(N,1)
	WRITE (OU, '(A,L2)')      "newton inner     ", AqNewtonInner
	CLOSE (OU)
	STOP "ReplayFEE done"
END SUBROUTINE MaybeReplayFEE


REAL*8 FUNCTION PowMemo (Slot, R, Base, Expo)
	USE Chemistry, ONLY : AqEquilibriaList
	IMPLICIT NONE
	INTEGER :: Slot, R
	REAL*8  :: Base, Expo
	IF (.NOT. ALLOCATED(PowMemoBase)) THEN
		ALLOCATE (PowMemoBase(3,SIZE(AqEquilibriaList,1)), PowMemoExp(3,SIZE(AqEquilibriaList,1)), &
		          PowMemoVal(3,SIZE(AqEquilibriaList,1)))
		PowMemoBase = -1.0D0 ; PowMemoExp = 0.0D0 ; PowMemoVal = 0.0D0
	END IF
	IF (Base .GT. 0.0D0 .AND. Base .LE. HUGE(Base)) THEN
		IF (Base .EQ. PowMemoBase(Slot,R) .AND. Expo .EQ. PowMemoExp(Slot,R)) THEN
			PowMemo = PowMemoVal(Slot,R)
			RETURN
		END IF
		PowMemo = Base ** Expo
		PowMemoBase(Slot,R) = Base ; PowMemoExp(Slot,R) = Expo ; PowMemoVal(Slot,R) = PowMemo
	ELSE
		PowMemo = Base ** Expo
	END IF
END FUNCTION PowMemo

!Patch C1: chemical-name tests of EquilibrateAqueousDissociationReaction depend only on the reaction.
SUBROUTINE InitEqRxnNameFlags ()
	USE Chemistry, ONLY : AqEquilibriaList, AqPhaseChemicalNames, AqAnionNames
	IMPLICIT NONE
	INTEGER :: R, I1, I3
	ALLOCATE (EqRxnAnionIsCO3(SIZE(AqEquilibriaList,1)), EqRxnIsLevitocite(SIZE(AqEquilibriaList,1)), &
	          EqRxnIsAmmSulfate(SIZE(AqEquilibriaList,1)))
	EqRxnAnionIsCO3 = .FALSE. ; EqRxnIsLevitocite = .FALSE. ; EqRxnIsAmmSulfate = .FALSE.
	DO R = 1, SIZE(AqEquilibriaList,1)
		I1 = INT(AqEquilibriaList(R,1)) ; I3 = INT(AqEquilibriaList(R,3))
		IF (I1 .GE. 1 .AND. I1 .LE. SIZE(AqPhaseChemicalNames)) THEN
			EqRxnIsLevitocite(R) = (AqPhaseChemicalNames(I1) .EQ. "(NH4)3H(SO4)2")
			EqRxnIsAmmSulfate(R) = (AqPhaseChemicalNames(I1) .EQ. "(NH4)2SO4")
		END IF
		IF (I3 .GE. 1 .AND. I3 .LE. SIZE(AqAnionNames)) EqRxnAnionIsCO3(R) = (AqAnionNames(I3) .EQ. "CO3--")
	END DO
END SUBROUTINE InitEqRxnNameFlags

RECURSIVE SUBROUTINE EquilibrateAqueousDissociationReaction &
        (WhichEquilibriumRxn, InParticle, OptErrorTolerance, ReturnType)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyAqChems, &
                                    HowManyAqCations, &
                                    AqPhaseChemicalNames,&
                                    HowManyAqAnions, AqAnionNames

	USE ModelParameters, ONLY : AqThermoEquilibriumError,		&
				AqThermoNumbAllowedIterations,	        &
				ProtonIndex,				&
				HydroxyIndex,				&
				WaterEquilibriumIndex

	USE GridPointFields, ONLY : GetRelativeHumidity

	USE InfrastructuralCode, ONLY : INT2STR, Real2Str, Transcript,  &
                                        ERROR, WARN, IsNaN

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichEquilibriumRxn
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: OptErrorTolerance
	INTEGER :: ReturnType  
        !! 1 if had to equilibrate, 2 if already in equilibrium

	!! Local Variables
	INTEGER :: I, J
	LOGICAL :: HaveEqRatio, Converged, StateChanged	!Patch N5, Patch 24
	REAL*8  :: R1, R2, R3, R4, R5, R8, R9, R22, R23, R24, R25, HalfIter	!Patch 25
	INTEGER :: IEl, ICatI, IAnI, IAn2I, IHydI	!Patch 25
	LOGICAL :: NotInfDiss, LCO3, LLev, LAmm	!Patch 25
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, &
                   ActivityCoefficients, ErrorTolerance, DRH, RH
	LOGICAL :: WaterRxn
	!Patch 28 (Newton single-reaction solve)
	REAL*8  :: B0Cat, B0An, B0An2, B0Hyd, B0El, B0W, X, XLo, XHi, XN, DG
	INTEGER :: NTry
	LOGICAL :: SolidRxn, Bad

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .false. 

	!Patch C1: fill the per-reaction chemical-name flags once
	IF (.NOT. ALLOCATED(EqRxnAnionIsCO3)) CALL InitEqRxnNameFlags ()
	IF (.NOT. AqSwitchesRead) CALL ReadAqSwitches ()	!Patch 28

        !! If this is an empty section, we have no truck with it
	!IF (InParticle%NumberOfParticles .EQ. 0 .OR. InParticle%Dry) RETURN
	if (inparticle%numberofparticles .lt. 1e-12 .or. inparticle%dry) then
		!print *, 'empty section, no truck'
		return
	endif
	
        ReturnType = 2
        II = 0.0D0 
	DRH = DeliquesenceRH (WhichEquilibriumRxn, InParticle)

	RH = GetRelativeHumidity ()

!!DRH and INFINITELY DISSOCIATING CHECK!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! If this is a reaction that isn't supposed to equilibrate                !!
!! (infinitely dissociating or a solid that has a DRH below the RH.)  Then !!
!! push into the dissociated form.                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! CMB (AER, Inc): More floating-point equality checks...

! 5 Oct 2016 CMB (AER): Fixing my own introduced bugs...
!    if ((aqequilibrialist(whichequilibriumrxn, 9) .le. -0.99999 .and. &
!    		aqequilibrialist(whichequilibriumrxn,9) .gt. -1.000001) .or. &
!    		((aqequilibrialist(whichequilibriumrxn, 6) .gt. 0.99999999 .and. &
!    				aqequilibrialist(whichequilibriumrxn, 6) .lt. 1.00000001) .and. &
!    				RH .ge. DRH)) then
	if ((dabs(AqEquilibriaList(WhichEquilibriumRxn, 9) + 1.) .lt. 0.000001) .or. &
			((dabs(AqEquilibriaList(WhichEquilibriumRxn, 6) - 1.) .lt. 0.000001) .and. &
					RH .ge. DRH)) then					
!	IF (AqEquilibriaList(WhichEquilibriumRxn,9) .EQ. -1. .OR. &
!		(AqEquilibriaList(WhichEquilibriumRxn,6) .EQ. 1.0 .AND. RH .GE. DRH)) THEN 

! 5 Oct 2016 CMB (AER): Fixed my own bug
!		if (aqequilibrialist(whichequilibriumrxn, 6) .gt. 0.99999999 .and. &
!				aqequilibrialist(whichequilibriumrxn, 6) .lt. 1.00000001) return
!		if (inparticle%aqchems(int(aqequilibrialist(whichequilibriumrxn,1))) .le. 1.0e-40) &
!			return
		if (dabs(AqEquilibriaList(WhichEquilibriumRxn, 1) + 1.) .lt. 0.000001) return
		if (dabs(InParticle%AqChems(int(AqEquilibriaList(WhichEquilibriumRxn, 1)))) .le. 1.0e-40) &
			return
!		IF (AqEquilibriaList(WhichEquilibriumRxn,1) .EQ. -1.) RETURN
!		IF (InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) .EQ. 0.) RETURN

                ReturnType = 1

		!! The change amount, two cases: equations with and without water
		! CMB (AER, Inc): Floating-point equality fix
!		IF (AqEquilibriaList(WhichEquilibriumRxn,8) .EQ. 0) THEN
!			II = InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1)))
!		ELSE
!			II = MIN(InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))),	& ! water limited, or
!			         InParticle%AqChems(1)/AqEquilibriaList(WhichEquilibriumRxn,8))	  ! electrolyte limited?
!		END IF
		IF (dabs(AqEquilibriaList(WhichEquilibriumRxn,8)) .lt. 0.00001) THEN
			II = InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1)))
		ELSE
			II = MIN(InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))),	& ! water limited, or
					 InParticle%AqChems(1)/AqEquilibriaList(WhichEquilibriumRxn,8))	  ! electrolyte limited?
		END IF
				
		!If it is already dissociated, it was already in equilibrium
                ! CMB (AER, Inc): Floating-point equality doesn't work right 
                !                 here
		!IF(II .EQ. 0.0) THEN
                if (dabs(ii) .le. 1.0E-40) then
			!IF (PRESENT(ReturnType)) ReturnType = 2
                        ReturnType = 2
			RETURN
		END IF
		
		!! Add to ions:
		!! CATION:
		InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2)))	=&
		InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2))) +&
		II * AqEquilibriaList(WhichEquilibriumRxn,4)

		!! ANION:
		InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) =&
		InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) +&
		II * AqEquilibriaList(WhichEquilibriumRxn,5)

		!! Second Anion for levitocite:
		! cmb
		if (dabs(AqEquilibriaList(WhichEquilibriumRxn,22)) .gt. 0.00001) then
		!IF (AqEquilibriaList(WhichEquilibriumRxn,22) .NE. 0.) THEN
			InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,22))) =&
			InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,22))) +&
			II * AqEquilibriaList(WhichEquilibriumRxn,23)	
		END IF

		!! Hydrates: Water equilibrium
		! cmb
		if (dabs(AqEquilibriaList(WhichEquilibriumRxn, 24)) .gt. 0.00001) then
		!IF (AqEquilibriaList(WhichEquilibriumRxn,24) .NE. 0.) THEN
			InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,24))) =&
			InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,24))) +&
			II * AqEquilibriaList(WhichEquilibriumRxn,25)	
		END IF
		
		!! And Adjust WATER CONTENT appropriately (for rxns that ionize water)
		! cmb
		if (inparticle%aqchems(1) .lt. 0.0) then
			call warn("InParticle%AqChems(1) < 0.0 in EquilibrateAqueousDissociationReaction")
		endif ! end cb
		InParticle%AqChems(1) = InParticle%AqChems(1) - II * AqEquilibriaList(WhichEquilibriumRxn,8)
		! CMB (AER, Inc): Put in 0 check
		!if (inparticle%aqchems(1) .lt. 0) inparticle%aqchems(1) = 0.0
		if (inparticle%aqchems(1) .lt. 0.0) then
			call warn("InParticle%AqChems(1) < 0.0 in EquilibrateAqueousDissociationReaction")
		endif ! end cb
				
		!! -- Then the Electrolyte
		IF (AqEquilibriaList(WhichEquilibriumRxn,1) .GT. 0.) &
		InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1)))	= &
					InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) - II
		! CMB (AER, Inc): Put in 0 check
		!if (inparticle%aqchems(1) .lt. 0) inparticle%aqchems(1) = 0.0
		
		!! Ask and warn if it is very dry
                ! CMB (AER, Inc): Add the RH print out, revise to make it a 
                !                 compound if statement
		IF (InParticle%AqChems(1) .EQ. 0 .OR. &
                        InParticle%AqChems(1) .LE. 1.0e-40) then
		    CALL WARN("An aerosol lost all of its water while equilibrating reaction <<"// &
                            TRIM(PrintAqEqReaction(WhichEquilibriumRxn))// &
                            ">>, which means that, in the best case, "// &
			    "the initial or current environmental humidity is extremely low.  RH = "// &
                            TRIM(REAL2STR(RH))//", II = "//TRIM(REAL2STR(II)))
                endif

		RETURN
	END IF

	!!Matt Alvarado added this to check for solid electrolyte reactions that are already in
	!!Equilibrium, as the solid conc = 0.0 and the solution is subsaturated
	!Patch 24: this EqRatio is reused below (the particle is not modified until the neutralisation step)
	StateChanged = .FALSE.
	EqRatio = EquilibriumConstantsRatio (WhichEquilibriumRxn, InParticle)
	! CMB: mod for floating-point inequality
        if (dabs(AqEquilibriaList(WhichEquilibriumRxn,6) - 1.0) &
                .lt. 0.00001) then
            if (EqRatio .lt. 1.0) then
                if (InParticle%AqChems(nint(aqequilibrialist(&
                        whichequilibriumrxn,1))) .lt. 1e-100) then
                    return
                endif
            endif
        endif
	!if (dabs(AqEquilibriaList(WhichEquilibriumRxn,6) - 1.0) .lt. 0.00001 &
	!!IF (AqEquilibriaList(WhichEquilibriumRxn,6) .EQ. 1.0 & !A solid electrolyte
	!    .AND. EqRatio .LT. 1.0 & !The solution is sub-saturated
	!	.AND. InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) .LT. 1e-100) THEN !There is no solid left (less than a molecule per particle)
	!	RETURN
	!END IF
	
	!!Matt Alvarado added this to skip carbonate equilibrium 
	!!when CO3-- concentration approaches 0.
	! cmb
	IF (dabs(AqEquilibriaList(WhichEquilibriumRxn,8) - 1.0) .lt. 0.00001 & !A hydrate reaction
	!IF (AqEquilibriaList(WhichEquilibriumRxn,8) .EQ. 1.0 & !A hydrate reaction
		.AND. EqRxnAnionIsCO3(WhichEquilibriumRxn) & !CO2 equilibrium only
	    .AND. EqRatio .GT. 1.0 & !The solution is super-saturated
		.AND. InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) &
                    .LT. 1e-30) THEN !There is no anion left
		!WRITE(*,*) EqRatio, WhichEquilibriumRxn, AqAnionNames(AqEquilibriaList(WhichEquilibriumRxn,3))
		RETURN
	END IF

	!For levitocite equilibrium as SO4-- approaches 0 (acidic aerosol)
	IF(AqEquilibriaList(WhichEquilibriumRxn,22) .GT. 0. & !levitocite
		.AND. EqRatio .GT. 1.0 & !The solution is supersaturated
		.AND. InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,22))) &
                    .LT. 1e-28) THEN !Not enough SO4--
		IF(EqRxnIsLevitocite(WhichEquilibriumRxn)) THEN
			RETURN
		END IF
	END IF

	!For (NH4)2SO4 equilibrium as SO4-- approaches 0 (acidic aerosol)
	! cmb
	IF (dabs(AqEquilibriaList(WhichEquilibriumRxn,6) - 1.0) .lt. 0.00001) THEN !A solid electrolyte
	!IF (AqEquilibriaList(WhichEquilibriumRxn,6) .EQ. 1.0) THEN !A solid electrolyte
		IF(EqRatio .GT. 1.0 & !The solution is supersaturated
			.AND. InParticle%AqChems(HowManyAqChems+HowManyAqCations+&
                            INT(AqEquilibriaList(WhichEquilibriumRxn,3))) .LT. 1e-28) THEN !Not enough SO4--
			IF(EqRxnIsAmmSulfate(WhichEquilibriumRxn)) THEN
				RETURN
			END IF
		END IF
	END IF
		
	!For water equilibrium as OH- approaches 0. (very high pH)
	IF(AqEquilibriaList(WhichEquilibriumRxn,1) .LT. 0.0 & !water (index is -1.0)
		.AND. EqRatio .GT. 1.0 & !The solution is supersaturated
		.AND. InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(&
                            WhichEquilibriumRxn,3))) .LT. 1e-40) THEN !Not enough OH-
		!write(*,*) 'Not enough OH-'
		RETURN
	END IF
		
	!! Set the error tolerance according to the 
	!! subroutine's optional input
	!IF (PRESENT(OptErrorTolerance)) THEN
		ErrorTolerance = OptErrorTolerance
	!ELSE
	!	ErrorTolerance = AqThermoEquilibriumError
	!END IF

	!! If all of the concentrations are zero, then equilibrium has effectively been established.
	! CMB (AER, Inc): Adjust for floating-point inequality
	!IF  (AqEquilibriaList(WhichEquilibriumRxn,1) .GT. 0) THEN
	!	IF ((Molality (INT(AqEquilibriaList(WhichEquilibriumRxn,1)), InParticle) .EQ. 0) .AND.	 &
    !		(CationMolality (INT(AqEquilibriaList(WhichEquilibriumRxn,2)), InParticle) .EQ. 0.  .OR.     &
	!	 AnionMolality  (INT(AqEquilibriaList(WhichEquilibriumRxn,3)), InParticle) .EQ. 0.) )RETURN
	!END IF
	IF  (AqEquilibriaList(WhichEquilibriumRxn,1) .GT. 0) THEN
		IF ((dabs(Molality(INT(AqEquilibriaList(WhichEquilibriumRxn,1)), InParticle)) .le. 1.0e-40) .AND.	 &
				(dabs(CationMolality(INT(AqEquilibriaList(WhichEquilibriumRxn,2)), InParticle)) .le. 1.0e-40  .OR.     &
				 dabs(AnionMolality(INT(AqEquilibriaList(WhichEquilibriumRxn,3)), InParticle)) &
				 .le. 1.0e-40)) then
			 !write(*,*) 'All Molalities are <= 1.0e-40; equilibrium has been established'
			!call flush(6) 
			RETURN
	    else
	    	!write(*,*) 'Molality = ', Molality(int(AqEquilibriaList(WhichEquilibriumRxn,1)),&
	    	!		InParticle)
			!write(*,*) 'CationMolality = ', CationMolality(INT(AqEquilibriaList(&
			!		WhichEquilibriumRxn,2)), InParticle)
			!write(*,*) 'AnionMolality = ', AnionMolality(INT(AqEquilibriaList(&
			!		WhichEquilibriumRxn,2)), InParticle)		
	    endif
	END IF
		
	IF(SCAFFOLDING) WRITE(*,*) "Before Eq Ratio"
	!! If equilibrium has been established, then it has been established...
	!Patch 24: EqRatio from the entry check above is still valid (nothing modifies the particle in between)
	IF(SCAFFOLDING) WRITE(*,*) "After Eq Ratio"

	! CMB: mod for double precision
	if (dabs(eqratio - 1) .lt. errortolerance) return
	!IF (ABS(EqRatio-1) .LT. ErrorTolerance) RETURN

	!IF (PRESENT(ReturnType)) ReturnType = 1
        ReturnType = 1

	!! NEUTRALIZE.
	!! If the reaction produces either a proton or a hydroxy ion and there is a water 
	!! reaction (as is proper), and there is an excess concentration of the other ion,
	!! Then push enough to the right to allow the water reaction to neutralize the
	!! excess protons.
	IF (WaterEquilibriumIndex .GT. 0 .AND. WhichEquilibriumRxn .NE. WaterEquilibriumIndex) THEN
		IF (AqEquilibriaList(WhichEquilibriumRxn,2) .EQ. ProtonIndex-HowManyAqChems .AND. & !If the cation is a proton
			InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) &
                            .GT. InParticle%AqChems(HydroxyIndex)   .AND. &
			InParticle%AqChems(HydroxyIndex) .GT. InParticle%AqChems(ProtonIndex) .AND. &
		    InParticle%AqChems(HydroxyIndex) .GT. 0.) THEN			


			StateChanged = .TRUE.	!Patch 24
			!! Make enough H+ to neutralize the OH-
			dX = MIN (InParticle%AqChems(HydroxyIndex) / AqEquilibriaList(WhichEquilibriumRxn,4),&
					  InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))))
		    !print *, 'dX1 (making enough H+) = ', dX
			call flush(6)
			InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2)))	&
					 = InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2)))&
					   + dX * AqEquilibriaList(WhichEquilibriumRxn,4)
	    	!print *, 'first one = ', InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2)))
			InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) = &
                            InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) &
					   + dX * AqEquilibriaList(WhichEquilibriumRxn,5)
	    	!print *, 'second one = ', InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3)))
			!! -- Then the Electrolyte
			InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) = &
					InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) - dX
			!print *, 'third one = ', InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1)))
			!! -- Then the Water Content
			!write(*,*) 'electrolyte after dX = ', inparticle%aqchems(int(aqequilibrialist(&
	     	!		whichequilibriumrxn,1)))
			InParticle%AqChems(1) = InParticle%AqChems(1) - dX * AqEquilibriaList(WhichEquilibriumRxn,8)
			!write(*,*) 'water content = ', inparticle%aqchems(1)
			
			!! -- Then call the water equilibrium			
			call flush(6)
			CALL EquilibrateAqueousDissociationReaction (WaterEquilibriumIndex, InParticle, &
                                OptErrorTolerance =AqThermoEquilibriumError, ReturnType=J)

		END IF

		IF (AqEquilibriaList(WhichEquilibriumRxn,3) .EQ. HydroxyIndex-HowManyAqChems-HowManyAqCations  .AND.	&
			InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) .GT. &
                            InParticle%AqChems(ProtonIndex) .AND. &
			InParticle%AqChems(ProtonIndex) .GT. InParticle%AqChems(HydroxyIndex)   .AND. &
		    InParticle%AqChems(ProtonIndex) .GT. 0.) THEN

			StateChanged = .TRUE.	!Patch 24
			!! Make enough H+ to neutralize the OH-
			dX = MIN (InParticle%AqChems(ProtonIndex) / AqEquilibriaList(WhichEquilibriumRxn,5),	&
					  InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))))
			!print *, 'dX2 (making enough H+) = ', dX
			InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2)))	&
					 = InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2))) &
					   + dX * AqEquilibriaList(WhichEquilibriumRxn,4)

			InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) = &
                            InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) &
					   + dX * AqEquilibriaList(WhichEquilibriumRxn,5)

			!! -- Then the Electrolyte
			InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) = &
					InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) - dX

			!! -- Then the Water Content
			InParticle%AqChems(1) = InParticle%AqChems(1) - dX * AqEquilibriaList(WhichEquilibriumRxn,8)

			!! -- Then call the water equilibrium
			CALL EquilibrateAqueousDissociationReaction (WaterEquilibriumIndex, &
                                InParticle, OptErrorTolerance =AqThermoEquilibriumError, ReturnType=J)

		END IF
	END IF

	IF(SCAFFOLDING) WRITE(*,*) "After Neutralize"

	!! I counts the number of iterations
	I = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine generally follows the Mass Flux Iteration method reviewed in !!
!! Jacobson, 1999 (the textbook) and reviewed earlier in Jacobson et al. 1996!!
!! and Villars 1959                                                          !!
!!									     !!
!! Jacobson, M.Z., A. Tabazadeh, and R.P. Turco, Simulating equilibrium      !!
!! within aerosols and nonequilibrium between gases and aerosols,            !!
!! Journal of Geophysical Research, 101 (D4), 9079-9091, 1996.		     !!
!!									     !!
!! Villars, D.S., A method of successive approximations for computing        !!
!! combustion equilibria on a high speed digital computer,                   !!
!! Journal of Physical Chemistry, 63, 521-5, 1959.			     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 1: Find the most aberrant ratio of concentration to                  !!
!!         stoicheometric coefficient                                        !!
!! (cf., Jacobson 1999, p.497)						     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (AqEquilibriaList(WhichEquilibriumRxn,22) .GT. 0) THEN
		!Levitocite
		IF(.NOT. EqRxnIsLevitocite(WhichEquilibriumRxn)) THEN
			CALL ERROR("Non-levitocite complex salt detected in EquilibrateAqueousDissociationReaction."//&
                               " This is not allowed.")
		END IF
		Qnumer = MIN( InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2)))/	&
				  AqEquilibriaList(WhichEquilibriumRxn,4),														&
				  InParticle%AqChems(HowManyAqChems+HowManyAqCations+&
                                      INT(AqEquilibriaList(WhichEquilibriumRxn,3)))/&
				  AqEquilibriaList(WhichEquilibriumRxn,5),														&
				  InParticle%AqChems(HowManyAqChems+HowManyAqCations+&
                                      INT(AqEquilibriaList(WhichEquilibriumRxn,22)))/&
				  AqEquilibriaList(WhichEquilibriumRxn,23))

	ELSE
		Qnumer = MIN( InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2)))&
				  /AqEquilibriaList(WhichEquilibriumRxn,4),														&
				  InParticle%AqChems(HowManyAqChems+HowManyAqCations+&
                                      INT(AqEquilibriaList(WhichEquilibriumRxn,3)))&
				  /AqEquilibriaList(WhichEquilibriumRxn,5))
	END IF
	
	!! If there is a non-water electrolyte, then water or the electrolyte may limit the reaction
	IF (AqEquilibriaList(WhichEquilibriumRxn,1) .GT. 0) THEN
		IF (AqEquilibriaList(WhichEquilibriumRxn,8) .GT. 0) THEN
			Qdenom = MIN(InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))),	&
				InParticle%AqChems(1)/AqEquilibriaList(WhichEquilibriumRxn,8))
		ELSE
			Qdenom = InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1)))
		END IF
		WaterRxn = .FALSE.
	ELSE
		Qdenom   = InParticle%AqChems(1)*EqCoeff(WhichEquilibriumRxn, InParticle, 1)
		WaterRxn = .TRUE.
	END IF


	!! STEP 2: Calculate the Mass Step Sizes by which to Correct Concentrations
	Z  = (Qdenom + Qnumer) / 2.	
	dX = Qdenom - Z				!! The Mass Flux Factor

	!! Loop over Corrections until convergence if there is anything to do
	! CMB: floating-point equality check
	if (.not. (dabs(qnumer) .le. 1.0e-40 .and. dabs(qdenom) .le. 1.0e-40)) then
	!IF (.NOT. (Qnumer .EQ. 0 .AND. Qdenom .EQ. 0)) THEN
	!Patch N5: the convergence test evaluates EquilibriumConstantsRatio on the particle state;
	!the previous iteration ended by evaluating the same function on the same (unmodified)
	!state.  Reuse that EqRatio instead of recomputing it (not valid for the infinitely
	!dissociating branch, which the test handles differently).
	!Patch 25: the AqEquilibriaList entries and index arithmetic used in the iteration loop are constants of the reaction;
	!read/convert them once instead of on every iteration (same values, same operations).
	R1 = AqEquilibriaList(WhichEquilibriumRxn,1)
	R2 = AqEquilibriaList(WhichEquilibriumRxn,2)
	R3 = AqEquilibriaList(WhichEquilibriumRxn,3)
	R4 = AqEquilibriaList(WhichEquilibriumRxn,4)
	R5 = AqEquilibriaList(WhichEquilibriumRxn,5)
	R8 = AqEquilibriaList(WhichEquilibriumRxn,8)
	R9 = AqEquilibriaList(WhichEquilibriumRxn,9)
	R22 = AqEquilibriaList(WhichEquilibriumRxn,22)
	R23 = AqEquilibriaList(WhichEquilibriumRxn,23)
	R24 = AqEquilibriaList(WhichEquilibriumRxn,24)
	R25 = AqEquilibriaList(WhichEquilibriumRxn,25)
	IEl   = INT(R1)
	ICatI = HowManyAqChems+INT(R2)
	IAnI  = HowManyAqChems+HowManyAqCations+INT(R3)
	IAn2I = HowManyAqChems+HowManyAqCations+INT(R22)
	IHydI = INT(R24)
	NotInfDiss = (DABS(R9 - 1) .GT. 1.0e-40)
	LCO3 = EqRxnAnionIsCO3(WhichEquilibriumRxn)
	LLev = EqRxnIsLevitocite(WhichEquilibriumRxn)
	LAmm = EqRxnIsAmmSulfate(WhichEquilibriumRxn)
	HalfIter = AqThermoNumbAllowedIterations/2.
	!Patch 24: if neutralisation did not change the particle, the entry EqRatio is also the first iteration's ratio
	HaveEqRatio = (.NOT. StateChanged) .AND. (NotInfDiss)
	IF (.NOT. (AqNewtonInner .AND. NotInfDiss)) THEN
	DO
		IF (HaveEqRatio) THEN
			Converged = (ErrorTolerance .GE. ABS(EqRatio - 1.))
		ELSE
			Converged = CheckEquilibriumConvergenceForAqDissociationRxn (WhichEquilibriumRxn, InParticle, ErrorTolerance)
		END IF
		IF (Converged) EXIT



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 3: Adjust Molalities: !!
	!!							  !!
	!! -- Ions First			  !!
	InParticle%AqChems(ICatI)					&
			 = InParticle%AqChems(ICatI)		&
			   + dX * R4

	InParticle%AqChems(IAnI)				&
			 = InParticle%AqChems(IAnI)	&
			   + dX * R5

	!WRITE(*,*) InParticle%AqChems(HowManyAqChems+HowManyAqCations+AqEquilibriaList(WhichEquilibriumRxn,3)), AqAnionNames(AqEquilibriaList(WhichEquilibriumRxn,3))
	
	IF (R22 .GT. 0) THEN
		!Levitocite
		InParticle%AqChems(IAn2I)				&
				= InParticle%AqChems(IAn2I)	&
				+ dX * R23
	END IF

	!! Hydrates: Water equilibrium
	! CMB (floating-point equality check)
	if (dabs(R24) .gt. 1.0e-40) then
	!IF (AqEquilibriaList(WhichEquilibriumRxn,24) .NE. 0.) THEN
		InParticle%AqChems(IHydI) =		&
				InParticle%AqChems(IHydI)		&
				+ dX * R25	
	END IF

	
	
	IF(SCAFFOLDING) WRITE(*,*) "Before ion check"
	!! Because we have to guess the dX size for some reactions (like for water), we could
	!! occasionally get negative concentrations here.  We will not allow this.  Should only
	!! happen with water reactions or the like, in which case these will go negative but the
	!! "electrolyte" concentration will not.
10	IF (InParticle%AqChems(ICatI) .LT. 0. .OR.		&
		InParticle%AqChems(IAnI) .LT. 0.) THEN

		dX = dX / 2.

		InParticle%AqChems(ICatI)					&
					 = InParticle%AqChems(ICatI)		&
					   - dX * R4

		InParticle%AqChems(IAnI)				&
				 = InParticle%AqChems(IAnI)	&
				   - dX * R5

		GOTO 10

	END IF
	IF(SCAFFOLDING) WRITE(*,*) "After ion check"

	!! -- Then the Electrolyte
	IF (R1 .GT. 0) THEN
		InParticle%AqChems(IEl) = InParticle%AqChems(IEl) - dX
		
		IF (InParticle%AqChems(IEl) .LT. 0.0) THEN
		    WRITE(*,*) "Negative Electrolyte Concentration"
			WRITE(*,*) "Electrolyte Conc.: ", InParticle%AqChems(IEl), " Electrolyte ", AqPhaseChemicalNames(IEl)
			WRITE(*,*) "Cation Conc.: ", InParticle%AqChems(ICatI)
			WRITE(*,*) "Anion Conc.: ", InParticle%AqChems(HowManyAqChems+INT(R3))
			WRITE(*,*) "dX = ", dX
			WRITE(*,*) "EqRatio = ", EqRatio
			WRITE(*,*) "Qdenom = ", Qdenom
			WRITE(*,*) "Qnumer = ", Qnumer
			WRITE(*,*) "Iteration Number: ", I
			write(*,*) 'AqEquilibriaList(WhichEquilibriumRxn, 1) = ', &
				AqEquilibriaList(WhichEquilibriumRxn, 1)
			
			! CMB: let's fake this too
			! update: maybe not?
			!inparticle%aqchems(int(aqequilibrialist(WhichEquilibriumRxn,1))) = 0.0D0
			!InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2))) = 0.0D0
			!InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) = 0.0D0
			STOP
						
		END IF
	
	END IF
	
	!! -- Then the Water Content
	InParticle%AqChems(1) = InParticle%AqChems(1) - dX * R8

	!! Do an error check and suggest a problem and maybe a hunch
	IF (InParticle%AqChems(1) .LE. 0 ) THEN

		CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle--NegativeWaterContent.txt", 	&
						InRelativeHumidity=.TRUE., InPh=.TRUE., InDoSurfaceTension=.FALSE., &
                                                 InDensity=.FALSE.,InDoWaterActivity=.FALSE.,	&
						InRadius = .FALSE., InIonicStrength = .FALSE.)

		IF(AerosolPH(InParticle) .EQ. 14) THEN
			CALL ERROR ("We achieved a negative water content and the pH of the particle is 14 while equilibrating an "// &
						"electrolyte in EquilibrateAqueousDissociationReaction(), which probably means that there was "// &
						"no way for ionic protons to form in the particle dumped to ")
		ELSE
			CALL ERROR ("We achieved a negative water content while equilibrating an electrolyte in "// &
						"EquilibrateAqueousDissociationReaction(), which probably means that there was "// &
						"no way for ionic protons to form in the particle dumped to ")
	END IF ; END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Z = 0.5*Z

	!! The recalculation is based on whether the equilibrium coefficient 
	!! is greater than or less than one	
	EqRatio = EquilibriumConstantsRatio (WhichEquilibriumRxn, InParticle)
	
	! CMB add
	! 7 October modification: set ratio to 1, don't stop'
	if (EqRatio .NE. EqRatio) then
		write(*, *) 'EqRatio is NaN for EquilibriumRxn#', WhichEquilibriumRxn
		write(*,*) 'Setting to 1 and returning'
		call flush(6)
		EqRatio = 1.0D0
		return
		!exit
		!stop
	endif
	
	IF (EqRatio .GT. 1.) dX = -1.*Z
	IF (EqRatio .LT. 1.) dX = Z
	

	!! Test to see if it's a water reaction that's
	!! no longer converging towards the solution.
	!! If it is (and this should only happen after large
	!! RH jumps, such as in validation runs, then we
	!! increase Z, which was a guess initially and so is 
	!! subject to error.  (other electrolytes have firm 
	!! constraints).  (the I .GT. 1 thing is so it has a
	!! previous value to compare to.)
	!!
	!! the minimum iteration number is because this can 
	!! occassionally happen early in the equilibration process
        if ( (dabs(II - EqRatio) .lt. ErrorTolerance/10.0)) then
            if (I .gt. 25) then
                if (WhichEquilibriumRxn .eq. 1) then
                    Z = Z/sqrt(errortolerance)
                else
                    !write(*,*) 'Bumping Reaction #', WhichEquilibriumRxn
                endif
            endif
        endif
	!if ( (dabs(II - EqRatio) .lt. ErrorTolerance/10.0) .and. &
	!		(I .gt. 25.) .and. &
	!		(WhichEquilibriumRxn .eq. 1) ) then
	!IF (ABS(II-EqRatio) .LT. ErrorTolerance/10. .AND. I .GT. 25. .AND. WhichEquilibriumRxn .EQ. 1)  THEN
	!	IF (WhichEquilibriumRxn .NE. 1) WRITE(*,*) "Bumping Reaction #", WhichEquilibriumRxn		
	!	Z = Z/SQRT(ErrorTolerance)
	!END IF

	!!Check to see if Electrolyte almost completely dissociated
	IF (R1 .GT. 0) THEN
		!For neutral Electrolytes
		IF (R1 .LE. HowManyAqChems .AND. InParticle%AqChems(IEl) .LT. 1.0e-28) THEN !Much Less than a molecule per particle
			!WRITE(*,*) "Solid Electrolyte ", trim(AqPhaseChemicalNames(AqEquilibriaList(WhichEquilibriumRxn,1))), " is completely dissociated"
			!InParticle%AqChems(AqEquilibriaList(WhichEquilibriumRxn,1)) = 0.0
			RETURN
		END IF
		
		!For ionic Electrolytes that may also approach 0 (ex. HSO4-)
		IF (R1 .GT. HowManyAqChems .AND. InParticle%AqChems(IEl) .LT. 1.0e-28) THEN !Much Less than a molecule per particle
			!WRITE(*,*) "Non-Solid Electrolyte ", trim(AqPhaseChemicalNames(AqEquilibriaList(WhichEquilibriumRxn,1))), " is completely dissociated"
			!InParticle%AqChems(AqEquilibriaList(WhichEquilibriumRxn,1)) = 0.0
			RETURN
		END IF

		!For carbonate equilibrium as CO3-- approaches 0
		IF(R8 .EQ. 1.0 & !A hydrate reaction
		    .AND. LCO3 & !CO2 equilibrium only
			.AND. InParticle%AqChems(INT(R3)) .LT. 1e-28) THEN
			RETURN
		END IF
		
		!For levitocite equilibrium as SO4-- approaches 0
		IF(LLev & !levitocite
		    .AND. InParticle%AqChems(IAn2I) .LT. 1e-28) THEN
			RETURN
		END IF

		!For (NH4)2SO4 equilibrium as SO4-- approaches 0
		IF(LAmm & 
		    .AND. InParticle%AqChems(IAnI) .LT. 1e-28) THEN
			RETURN
		END IF

	END IF

	!For water equilibrium as OH- approaches 0. (very high pH)
	IF(R1 .LT. 0. & !water has index -1.0
		.AND. InParticle%AqChems(IAnI) .LT. 1e-40 & !Not enough OH-
		.AND. I .GT. HalfIter) THEN !Made a good-faith effort
		RETURN
	END IF

	!! Store this EqRatio for testing against the next one, in the next loop
	II = EqRatio
	HaveEqRatio = (NotInfDiss)	!Patch N5

	!! Count the number of attempts before convergence
	I = I+1
	
	!! Just exit if no convergence after nax iterations (MJA, 100507)
	! CMB (AER): Tried return instead of exit?
	IF (I .GT. AqThermoNumbAllowedIterations) THEN
		!write(*,*) 'reached max number of iterations'
		!call flush(6)
		EXIT
	END IF

	END DO
	ELSE
	!Patch 28: Newton single-reaction solve (AqInnerSolver = 1, the default).  The mass-flux iteration moves the
	!extent X by a halving step Z (bisection).  Here X is found by Newton's method on
	!G(X) = ln(EqRatio), kept inside the same bracket [-Qnumer, Qdenom], with a bisection step
	!whenever Newton would leave it.  Every state is base + stoichiometry * X, so mass and charge
	!are conserved exactly, as in the mass-flux iteration.  The activity coefficients and the
	!water activity are held fixed, as above.  Convergence test and early returns are the same.
	SolidRxn = (DABS(AqEquilibriaList(WhichEquilibriumRxn,6) - 1.0) .LT. 0.00001)
	B0Cat = InParticle%AqChems(ICatI) ; B0An = InParticle%AqChems(IAnI)
	B0An2 = 0. ; IF (R22 .GT. 0) B0An2 = InParticle%AqChems(IAn2I)
	B0Hyd = 0. ; IF (DABS(R24) .GT. 1.0e-40) B0Hyd = InParticle%AqChems(IHydI)
	B0El  = 0. ; IF (R1 .GT. 0) B0El = InParticle%AqChems(IEl)
	B0W   = InParticle%AqChems(1)
	XLo = -Qnumer ; XHi = Qdenom ; X = 0.
	IF (.NOT. HaveEqRatio) EqRatio = EquilibriumConstantsRatio (WhichEquilibriumRxn, InParticle)
	I = 1
	DO
		IF (ErrorTolerance .GE. ABS(EqRatio - 1.)) EXIT
		IF (EqRatio .GT. 1.) THEN
			XHi = X
		ELSE
			XLo = X
		END IF
		!dG/dX: each species enters ln Q as stoich*ln(n/W); W is the water amount
		DG = R4*R4/MAX(InParticle%AqChems(ICatI),1.0D-300) + R5*R5/MAX(InParticle%AqChems(IAnI),1.0D-300) &
		     + (R4 + R5)*R8/MAX(InParticle%AqChems(1),1.0D-300)
		IF (R22 .GT. 0) DG = DG + R23*R23/MAX(InParticle%AqChems(IAn2I),1.0D-300) &
		                        + R23*R8/MAX(InParticle%AqChems(1),1.0D-300)
		IF (R1 .GT. 0 .AND. .NOT.SolidRxn) DG = DG + 1./MAX(InParticle%AqChems(IEl),1.0D-300) &
		                                         - R8/MAX(InParticle%AqChems(1),1.0D-300)
		IF (EqRatio .GT. 0. .AND. EqRatio .LT. 1.0D300 .AND. DG .GT. 0.) THEN
			XN = X - LOG(EqRatio)/DG
		ELSE
			XN = 0.5*(XLo + XHi)
		END IF
		!A solid that is still undersaturated with all of it dissolved: dissolve it all
		IF (SolidRxn .AND. R1 .GT. 0 .AND. EqRatio .LT. 1. .AND. XN .GE. XHi .AND. XHi .EQ. B0El) THEN
			XN = XHi
		ELSE IF (.NOT. (XN .GT. XLo .AND. XN .LT. XHi)) THEN
			XN = 0.5*(XLo + XHi)
		END IF
		DO NTry = 1, 60
			InParticle%AqChems(ICatI) = B0Cat + XN*R4
			InParticle%AqChems(IAnI)  = B0An  + XN*R5
			IF (R22 .GT. 0) InParticle%AqChems(IAn2I) = B0An2 + XN*R23
			IF (DABS(R24) .GT. 1.0e-40) InParticle%AqChems(IHydI) = B0Hyd + XN*R25
			IF (R1 .GT. 0) InParticle%AqChems(IEl) = B0El - XN
			InParticle%AqChems(1) = B0W - XN*R8
			Bad = (InParticle%AqChems(ICatI) .LT. 0. .OR. InParticle%AqChems(IAnI) .LT. 0.)
			IF (R22 .GT. 0) Bad = Bad .OR. (InParticle%AqChems(IAn2I) .LT. 0.)
			IF (Bad) THEN
				XLo = XN ; XN = 0.5*(XLo + XHi) ; CYCLE
			END IF
			Bad = .FALSE.
			IF (R1 .GT. 0) Bad = (InParticle%AqChems(IEl) .LT. 0.)
			Bad = Bad .OR. (InParticle%AqChems(1) .LE. 0.)
			IF (Bad) THEN
				XHi = XN ; XN = 0.5*(XLo + XHi) ; CYCLE
			END IF
			EXIT
		END DO
		X = XN
		EqRatio = EquilibriumConstantsRatio (WhichEquilibriumRxn, InParticle)
		if (EqRatio .NE. EqRatio) then
			write(*, *) 'EqRatio is NaN for EquilibriumRxn#', WhichEquilibriumRxn
			write(*,*) 'Setting to 1 and returning'
			call flush(6)
			EqRatio = 1.0D0
			return
		endif
		IF (R1 .GT. 0) THEN
			IF (InParticle%AqChems(IEl) .LT. 1.0e-28) RETURN !Much less than a molecule per particle
			IF (R8 .EQ. 1.0 .AND. LCO3 .AND. InParticle%AqChems(INT(R3)) .LT. 1e-28) RETURN
			IF (LLev .AND. InParticle%AqChems(IAn2I) .LT. 1e-28) RETURN
			IF (LAmm .AND. InParticle%AqChems(IAnI) .LT. 1e-28) RETURN
		END IF
		IF (R1 .LT. 0. .AND. InParticle%AqChems(IAnI) .LT. 1e-40 .AND. I .GT. HalfIter) RETURN
		I = I+1
		IF (I .GT. AqThermoNumbAllowedIterations) EXIT
	END DO
	END IF ; END IF

	!! Report number of iterations if transcribing
	IF (SCAFFOLDING) CALL TRANSCRIPT("Number of Iterations to Converge Reaction "//  &
									 TRIM(INT2STR(WhichEquilibriumRxn))//": "//TRIM(INT2STR(I)))

	RETURN
END SUBROUTINE EquilibrateAqueousDissociationReaction

!! When doing dissolution reactions, it is necessary to reform the electrolyte initially 
!! during the first step of the equilibration process.  This lets the electrolyte escape
!! the aerosol phase if it is going to without many, many iterations (as would be necessary
!! to pull the ionized form to the electrolytic form and then out of the aerosol).  This
!! is called from EquilibrateGridPoint() in CondensationRelatedFunctions.h
SUBROUTINE ReformElectrolyte(WhichRxn, InParticle)

	USE Chemistry,       ONLY : AqEquilibriaList,				&
								HowManyAqEqReactions,			&
								HowManyAqChems,					&
								HowManyAqCations

	USE ModelParameters, ONLY : ProtonIndex,					&
								HydroxyIndex,					&
								WaterEquilibriumIndex

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn
	TYPE(PARTICLE),POINTER :: InParticle

	INTEGER :: RxnType
	REAL*8  :: dX

	RxnType = -1

	!! If this is an empty section, we have no truck with it
	IF (InParticle%NumberOfParticles .EQ. 0) RETURN

	!! Five Possibilities:
	!! 0. the reaction number is out of bounds or
	!!    there is no water reaction
	IF (WhichRxn .LE. 0 .OR. WhichRxn .GT. HowManyAqEqReactions) THEN
		RETURN

	!! there is no water reaction
	ELSE IF (WaterEquilibriumIndex .LE. 0) THEN

		!! 1. nothing should be done to the water equation
		RxnType = 1

	ELSE IF (AqEquilibriaList(WhichRxn,2) .EQ. ProtonIndex - HowManyAqChems) THEN
		
		!! 2. a proton is involved and is the limiting reagent
		RxnType = 2

		!! 1. a proton is involved but it is not the limiting reactant in reformation
		IF (InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichRxn,3))) / AqEquilibriaList(WhichRxn,5) .LT. &
			InParticle%AqChems(ProtonIndex) / AqEquilibriaList(WhichRxn,4)) &
			RxnType = 1

	!! the reaction incorporates a hydroxy ion
	ELSE IF (AqEquilibriaList(WhichRxn,3) .EQ. HydroxyIndex - HowManyAqChems - HowManyAqCations) THEN

		!! 3. a hyrdoxy ion is involved and is the limiting reagent
		RxnType = 3

		!! 1. a hyroxy ion is involved but it is not the limiting reactant in reformation
		IF (InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichRxn,2))) / AqEquilibriaList(WhichRxn,4) .LT. &
			InParticle%AqChems(HydroxyIndex) / AqEquilibriaList(WhichRxn,5)) &
			RxnType = 1

	END IF

	!! Now Adjust the water reaction
	SELECT CASE (RxnType)

	CASE (2)

		!! DISSOCIATE ENOUGH WATER SO HAVE ENOUGH H+ to PROCEED
		dX = MAX(0., InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichRxn,3)))/		&
				     AqEquilibriaList(WhichRxn,5)*AqEquilibriaList(WhichRxn,4)-InParticle%AqChems(ProtonIndex))

		InParticle%AqChems(1)            = InParticle%AqChems(1) - dX
		InParticle%AqChems(ProtonIndex)  = InParticle%AqChems(ProtonIndex) + dX
		InParticle%AqChems(HydroxyIndex) = InParticle%AqChems(HydroxyIndex) + dX

	CASE (3)

		!! DISSOCIATE ENOUGH WATER SO HAVE ENOUGH OH- to PROCEED
		dX = MAX(0., InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichRxn,2)))/AqEquilibriaList(WhichRxn,4)*	&
					 AqEquilibriaList(WhichRxn,5)-InParticle%AqChems(HydroxyIndex))

		InParticle%AqChems(1)            = InParticle%AqChems(1) - dX
		InParticle%AqChems(ProtonIndex)  = InParticle%AqChems(ProtonIndex) + dX
		InParticle%AqChems(HydroxyIndex) = InParticle%AqChems(HydroxyIndex) + dX



	END SELECT
	
	!! Now adjust the reaction itself
	dX = MIN(InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichRxn,2)))/AqEquilibriaList(WhichRxn,4), &
		     InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichRxn,3)))/AqEquilibriaList(WhichRxn,5))

	!! Water
	InParticle%AqChems(1) = InParticle%AqChems(1) + dX * AqEquilibriaList(WhichRxn,8)

	!! Electrolyte
	InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,1))) = InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,1))) + dX

	!! Cation
	InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,2))+HowManyAqChems) =												&
													InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,2))+HowManyAqChems) &
													- dX * AqEquilibriaList(WhichRxn,4)

	!! Anion
	InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,3))+HowManyAqChems+HowManyAqCations) =										&
									InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,3))+HowManyAqChems+HowManyAqCations)		&
									- dX * AqEquilibriaList(WhichRxn,5)
	
	!Levitocite case - second anion
	IF(AqEquilibriaList(WhichRxn,22) .NE. 0) THEN
		InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,22))+HowManyAqChems+HowManyAqCations) =										&
									InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,22))+HowManyAqChems+HowManyAqCations)		&
									- dX * AqEquilibriaList(WhichRxn,23)
	END IF

	!! Hydrates: Water equilibrium
	IF (AqEquilibriaList(WhichRxn,24) .NE. 0.) THEN
		InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,24))) =		&
			InParticle%AqChems(INT(AqEquilibriaList(WhichRxn,24)))		&
			- dX * AqEquilibriaList(WhichRxn,25)	
	END IF

	RETURN
END SUBROUTINE ReformElectrolyte


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns .TRUE. if it thinks a given dissociation reaction has !!
!! converged, .FALSE. otherwise									 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION CheckEquilibriumConvergenceForAqDissociationRxn (WhichEquilibriumRxn, InParticle, ErrorTolerance)

	USE Chemistry,       ONLY : AqEquilibriaList ! AqIonGrid, 

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichEquilibriumRxn
	REAL*8  :: ErrorTolerance
	TYPE(PARTICLE),POINTER :: InParticle

	!! These are infinitely dissociating reactions
	! CMB: floating-point inequality check
	if (dabs(aqequilibrialist(WhichEquilibriumRxn, 9) - 1) .le. 1.0e-40) then
	!IF (AqEquilibriaList(WhichEquilibriumRxn,9) .EQ. -1.) THEN
		! CMB: floating-point check
		IF (dabs(InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1)))) .le. 1.0e-40) THEN
		!IF (InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) .EQ. 0) THEN
			CheckEquilibriumConvergenceForAqDissociationRxn = .TRUE.
		ELSE
			CheckEquilibriumConvergenceForAqDissociationRxn = .FALSE.
		END IF
		RETURN
	END IF

	!! Convergence has been acceptably reached if the ratio of ideal and actual equilibrium ratio
	!! are within the error value "AqThermoEquilibriumError" found in ModelParameters
	CheckEquilibriumConvergenceForAqDissociationRxn = .FALSE.
	IF(ErrorTolerance .GE. abs(EquilibriumConstantsRatio (WhichEquilibriumRxn, InParticle) - 1.))	&
		CheckEquilibriumConvergenceForAqDissociationRxn = .TRUE.
	
	RETURN
END FUNCTION CheckEquilibriumConvergenceForAqDissociationRxn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Check the Ratio between the Existing and Ideal Equilibrium   !!
!! Constants.  This is used by several other routines.		    !!
!!															    !!
!! When the equations are read, the program tries to figure out !!
!! whether an uncharged electrolyte is dissociating (in which   !!
!! case the denominator activity coefficient is simply 1), and  !!
!! if not it stores an equivalent activity coefficient ratio.   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Patch 26: EquilibriumConstantsRatio is a pure function of the values listed in its memo key below
!(the electrolyte, cation, anion and levitocite-anion amounts, the water amount, the number
!concentration, the Dry flag, the temperature, the water activity and the two activity
!coefficients it can read).  FindElectrolyteEquilibrium re-evaluates every major reaction on
!every iteration although only the reaction just equilibrated changed the particle, so the
!same arguments recur.  The key is compared bit for bit (TRANSFER to INTEGER*8), so a hit
!returns exactly the value a fresh evaluation would.  NOT thread-safe (like patches N7/N11/C2).
REAL*8 FUNCTION EquilibriumConstantsRatio (WhichEquilibriumRxn, InParticle)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyAqChems, HowManyAqCations

	IMPLICIT NONE

	INTEGER :: WhichEquilibriumRxn
	TYPE(PARTICLE),POINTER :: InParticle

	INTEGER*8 :: Key(ECRMemoKeyLen)
	INTEGER   :: R, NA, Idx, G2

	R  = WhichEquilibriumRxn
	NA = SIZE(InParticle%AqChems)
	IF (.NOT. ALLOCATED(ECRMemoKey)) THEN
		ALLOCATE (ECRMemoKey(ECRMemoKeyLen, SIZE(AqEquilibriaList,1)), ECRMemoVal(SIZE(AqEquilibriaList,1)), &
		          ECRMemoOK(SIZE(AqEquilibriaList,1)))
		ECRMemoOK = .FALSE.
	END IF

	Key = 0_8
	Idx = INT(AqEquilibriaList(R,1))
	IF (Idx .GE. 1 .AND. Idx .LE. NA) Key(1) = TRANSFER(InParticle%AqChems(Idx), 0_8)
	Idx = HowManyAqChems + INT(AqEquilibriaList(R,2))
	IF (Idx .GE. 1 .AND. Idx .LE. NA) Key(2) = TRANSFER(InParticle%AqChems(Idx), 0_8)
	Idx = HowManyAqChems + HowManyAqCations + INT(AqEquilibriaList(R,3))
	IF (Idx .GE. 1 .AND. Idx .LE. NA) Key(3) = TRANSFER(InParticle%AqChems(Idx), 0_8)
	Idx = HowManyAqChems + HowManyAqCations + INT(AqEquilibriaList(R,22))
	IF (AqEquilibriaList(R,22) .NE. 0. .AND. Idx .GE. 1 .AND. Idx .LE. NA) &
		Key(4) = TRANSFER(InParticle%AqChems(Idx), 0_8)
	Key(5) = TRANSFER(InParticle%AqChems(1), 0_8)
	Key(6) = TRANSFER(InParticle%NumberOfParticles, 0_8)
	IF (InParticle%Dry) Key(7) = 1_8
	Key(8) = TRANSFER(InParticle%Temperature, 0_8)
	Key(9) = TRANSFER(InParticle%WaterActivity, 0_8)
	Key(10) = TRANSFER(InParticle%GammaMixed(R), 0_8)
	G2 = FLOOR(AqEquilibriaList(R,16))
	IF (G2 .GE. 1 .AND. G2 .LE. SIZE(InParticle%GammaMixed)) Key(11) = TRANSFER(InParticle%GammaMixed(G2), 0_8)

	IF (ECRMemoOK(R)) THEN
		IF (ALL(Key .EQ. ECRMemoKey(:,R))) THEN
			EquilibriumConstantsRatio = ECRMemoVal(R)
			RETURN
		END IF
	END IF

	EquilibriumConstantsRatio = EquilibriumConstantsRatioEval (WhichEquilibriumRxn, InParticle)
	ECRMemoKey(:,R) = Key
	ECRMemoVal(R)   = EquilibriumConstantsRatio
	ECRMemoOK(R)    = .TRUE.

END FUNCTION EquilibriumConstantsRatio

!Patch 26: the unchanged evaluation, called by the memo above on a miss.
REAL*8 FUNCTION EquilibriumConstantsRatioEval (WhichEquilibriumRxn, InParticle)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyAqChems, HowManyAqCations

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichEquilibriumRxn
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8	:: MM, Eq, AwPow, Dum1

	!! If this should be an infinitely dissociating reaction, then it either
	!! is in equilibrium or it is infinitely far away.
	! CMB: Mod for floating-point inequality
	if (dabs(aqequilibrialist(WhichEquilibriumRxn,9) - 1) .le. 1.0e-40) then
	!IF (AqEquilibriaList(WhichEquilibriumRxn,9) .EQ. -1.) THEN
		if (dabs(inparticle%aqchems(int(aqequilibrialist(WhichEquilibriumRxn, 1))))&
				.le. 1.0e-40) then
		!IF (InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) .EQ. 0) THEN
			EquilibriumConstantsRatioEval = 1.
		ELSE
			EquilibriumConstantsRatioEval = 0.000001
		END IF
		RETURN
	END IF

	!! If there is none of the electrolyte and none of one of the constituent ions, 
	!! then equilibrium has been established because everything is zero, essentially
	IF (AqEquilibriaList(WhichEquilibriumRxn,1) .GT. 0) THEN
		IF (InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibriumRxn,1))) .EQ. 0)  THEN
			IF((InParticle%AqChems(HowManyAqChems+INT(AqEquilibriaList(WhichEquilibriumRxn,2))) .EQ. 0.) .OR.	 &
			   (InParticle%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(WhichEquilibriumRxn,3))) .EQ. 0.))  THEN
					EquilibriumConstantsRatioEval = 1.
					RETURN
	END IF ; END IF; END IF


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Prepare the needed quantities !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Get the molality of the electrolyte and the equilibrium at this temperature
	Eq = EqCoeff (WhichEquilibriumRxn, InParticle, 1)  !! The "1" indicates a aqueous dissociation reaction

	!! If there is only water, then get we don't need the molality
	IF (AqEquilibriaList(WhichEquilibriumRxn,8) .GT. 0 .AND. AqEquilibriaList(WhichEquilibriumRxn,1) .EQ. -1) THEN
		MM = 1.
	!If it is a solid electrolyte, set MM (molality of undissociated species) = 1
	ELSE IF (AqEquilibriaList(WhichEquilibriumRxn,6) .EQ. 1.0) THEN
		MM = 1.	
	ELSE	
		MM = Molality (INT(AqEquilibriaList(WhichEquilibriumRxn,1)), InParticle) 
	END IF


	!! If there is a water molecule in the equilibrium, then we
	!! need to calculate the water activity.
	IF (AqEquilibriaList(WhichEquilibriumRxn,8) .GT. 0  ) THEN
		!Patch N8: pow(x,1.0) is exactly x
		IF (AqEquilibriaList(WhichEquilibriumRxn,8) .EQ. 1.0) THEN
			AwPow = InParticle%WaterActivity
		ELSE
			AwPow = PowMemo(2, WhichEquilibriumRxn, InParticle%WaterActivity, AqEquilibriaList(WhichEquilibriumRxn,8))
		END IF
	ELSE IF (AqEquilibriaList(WhichEquilibriumRxn,24) .NE. 0 ) THEN
		!Hydrates
		IF (AqEquilibriaList(WhichEquilibriumRxn,25) .EQ. 1.0) THEN
			AwPow = InParticle%WaterActivity
		ELSE
			AwPow = PowMemo(3, WhichEquilibriumRxn, InParticle%WaterActivity, AqEquilibriaList(WhichEquilibriumRxn,25))
		END IF
	ELSE	
		AwPow = 1.
	END IF

	!! If the solution is set to blow up, then establish a large artificial value.
	IF (MM .eq. 0 .or. Eq .eq. 0) THEN
		EquilibriumConstantsRatioEval = 1.e12
		RETURN
	END IF


	!! This triggers for the bisulfate-like reactions that require 
	!! activity coefficient ratios.
	IF ((AqEquilibriaList(WhichEquilibriumRxn, 17) .GT. 0) .AND. &
		(AqEquilibriaList(WhichEquilibriumRxn, 18) .EQ. 0) .AND. &
		(AqEquilibriaList(WhichEquilibriumRxn, 22) .EQ. 0) .AND. &
		(AqEquilibriaList(WhichEquilibriumRxn, 24) .EQ. 0))	THEN

	  EquilibriumConstantsRatioEval  =																								&
		CationMolality(INT(AqEquilibriaList(WhichEquilibriumRxn,2)),InParticle) ** AqEquilibriaList(WhichEquilibriumRxn, 4) *	&
	    AnionMolality (INT(AqEquilibriaList(WhichEquilibriumRxn,3)),InParticle) ** AqEquilibriaList(WhichEquilibriumRxn, 5) *	&
		InParticle%GammaMixed(FLOOR(AqEquilibriaList(WhichEquilibriumRxn, 16)))													&
		** ANINT(10000. * (AqEquilibriaList(WhichEquilibriumRxn, 16) - FLOOR(AqEquilibriaList(WhichEquilibriumRxn, 16)))) /		&
		InParticle%GammaMixed(WhichEquilibriumRxn)																				&
		** ANINT(10000. * (AqEquilibriaList(WhichEquilibriumRxn, 17) - FLOOR(AqEquilibriaList(WhichEquilibriumRxn, 17)))) /		&
		AwPow / MM / Eq

	ELSE IF (AqEquilibriaList(WhichEquilibriumRxn,6) .EQ. 1.0 .AND. &
			 (AqEquilibriaList(WhichEquilibriumRxn, 24) .EQ. 0)) THEN
		!Solid Formation Reaction (non-hydrate)
		EquilibriumConstantsRatioEval = MeanActivity (WhichEquilibriumRxn, InParticle) / Eq 
	
	ELSE IF(AqEquilibriaList(WhichEquilibriumRxn,6) .EQ. 1.0 .AND. &
			 (AqEquilibriaList(WhichEquilibriumRxn, 24) .NE. 0)) THEN
		!Solid formation reaction (hydrate)
		EquilibriumConstantsRatioEval = MeanActivity (WhichEquilibriumRxn, InParticle)*AwPow / Eq 

	ELSE
		!Canonical case
		Dum1 = MeanActivity (WhichEquilibriumRxn, InParticle) / MM / Eq / AwPow
		IF (Dum1 .LT. 0.0) THEN
			WRITE(*,*) "Negative EquilibriumConstantsRatio", WhichEquilibriumRxn
			WRITE(*,*) "Mean Activity", MeanActivity (WhichEquilibriumRxn, InParticle)
			WRITE(*,*) "Eq", Eq, "MM", MM
			WRITE(*,*) "Water Activity", AwPow
			STOP
		END IF
		!Patch N3: Dum1 above is this same expression; do not call MeanActivity again.
		EquilibriumConstantsRatioEval = Dum1
	END IF

	RETURN
END FUNCTION EquilibriumConstantsRatioEval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Fill the Ion Fraction Vectors for a given particle. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CalculateIonicStrength (InParticle)

	USE Chemistry, ONLY : HowManyAqCations,			&
	                      HowManyAqAnions,			&
			      AqCationCharge,			&
			      AqAnionCharge

	IMPLICIT NONE

	!! External Variables
	TYPE (Particle) :: InParticle

	!! Internal Variables
	INTEGER :: I
	REAL*8  :: II

	!! If this is an empty section, we have no truck with it
	! CMB: floating point inequality
	if (inparticle%numberofparticles .le. 1.0e-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0 .OR. InParticle%Dry) RETURN

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Ionic Strength Fractions and Ionic Strengths !!
	!! Are functions of the Cationic and Anionic	!!
	!! Strengths.  Calculate those first.			!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	II = 0.
	DO I = 1, HowManyAqCations
		II = II + CationMolality(I, InParticle) * AqCationCharge(I) * AqCationCharge(I)
	END DO

	DO I = 1, HowManyAqAnions
		II = II + AnionMolality(I, InParticle) * AqAnionCharge(I) * AqAnionCharge(I)
	END DO

	!! And then load the total ionic strength
	InParticle%IonicStr = II / 2.
	
	RETURN
END SUBROUTINE CalculateIonicStrength 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Mean Activity of an Electrolyte !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION MeanActivity (WhichEquilibrium, InParticle)

	USE Chemistry, ONLY : AqEquilibriaList, AqCationNames, AqAnionNames, &
						  HowManyAqChems, HowManyAqCations

	USE InfrastructuralCode


	IMPLICIT NONE

	!! External Variables
	INTEGER					:: WhichEquilibrium
	TYPE(PARTICLE),POINTER	:: InParticle
	REAL*8 :: CatM, AnM, LevM	!Patch N4: molalities evaluated once, not up to 4 times
	REAL*8 :: E4, E5, PC, PA	!Patch N8

	!! If this is an empty section, we have no truck with it
        ! CMB (AER, Inc): Floating-point equality issue
        if (inparticle%numberofparticles .le. 1.0E-40) then
	!IF (InParticle%NumberOfParticles .EQ. 0) THEN
		MeanActivity = 0.
		RETURN
	END IF

	!Check for negative Molality
        ! CMB (AER, Inc): In order to debug some problems, we are going to uncomment the "Call DumpParticleContentsAtError" 
        !                 and "Stop" lines.  Need to fill in some arguments first...
	CatM = CationMolality(INT(AqEquilibriaList(WhichEquilibrium,2)),InParticle)
	IF (CatM .LT. 0.0) THEN
		!WRITE(*,*) "Negative Cation Molality, Reaction # ", TRIM(INT2STR(WhichEquilibrium)), " ", TRIM(AqCationNames(INT(AqEquilibriaList(WhichEquilibrium,2)))), " ", TRIM(REAL2STR(CationMolality(INT(AqEquilibriaList(WhichEquilibrium,2)),InParticle)))
		!WRITE(*,*) "Resetting to 0."
		InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibrium,2)) + HowManyAqChems) = 0.0
		
		CALL DumpParticleContentsAtError (InParticle, "NegMolality.txt", .false., &
				.false., .false., .false., .false., .false., .false.)
		STOP
	ELSE
	  AnM = AnionMolality(INT(AqEquilibriaList(WhichEquilibrium,3)),InParticle)
	END IF
	IF (AnM .LT. 0.0) THEN
		!WRITE(*,*) "Negative Anion Molality, Reaction # ", TRIM(INT2STR(WhichEquilibrium)), " ", TRIM(AqAnionNames(INT(AqEquilibriaList(WhichEquilibrium,3)))), " ", TRIM(REAL2STR(AnionMolality(INT(AqEquilibriaList(WhichEquilibrium,3)),InParticle)))
		!WRITE(*,*) "Resetting to 0."
		InParticle%AqChems(INT(AqEquilibriaList(WhichEquilibrium,3)) + HowManyAqChems + HowManyAqCations) = 0.0

		CALL DumpParticleContentsAtError (InParticle, "NegMolality.txt", .false., &
				.false., .false., .false., .false., .false., .false.)
		STOP
	END IF


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Activity Coefficient is defined as:					!!
	!!									!!
	!! a_i,j = m_i ^ nu_i * m_j ^ nu_j * Gamma_i,j ^ (nu_i + nu_j)		!!
	!!									!!
	!! where m is molality, nu is a stoicheometric function, and Gamma is	!!
	!! the mean activity coefficient					!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Patch N8: pow(x,1.0) is exactly x; skip the pow call for unit stoichiometry
	E4 = AqEquilibriaList(WhichEquilibrium, 4)
	E5 = AqEquilibriaList(WhichEquilibrium, 5)
	IF (E4 .EQ. 1.0) THEN
		PC = CatM
	ELSE
		PC = CatM**E4
	END IF
	IF (E5 .EQ. 1.0) THEN
		PA = AnM
	ELSE
		PA = AnM**E5
	END IF

	IF (AqEquilibriaList(WhichEquilibrium,22) .EQ. 0.) THEN
		MeanActivity = PC * &
				   PA * &
				   PowMemo(1, WhichEquilibrium, InParticle%GammaMixed(WhichEquilibrium), AqEquilibriaList(WhichEquilibrium, 4) &
				   + AqEquilibriaList(WhichEquilibrium, 5))
	END IF

	!Correction for Levitocite
	IF (AqEquilibriaList(WhichEquilibrium,22) .NE. 0.) THEN
		MeanActivity = PC * &
				   PA * &
				   AnionMolality(INT(AqEquilibriaList(WhichEquilibrium,22)),InParticle) **AqEquilibriaList(WhichEquilibrium, 23) * &
				   InParticle%GammaMixed(WhichEquilibrium)**(AqEquilibriaList(WhichEquilibrium, 4) &
				   + AqEquilibriaList(WhichEquilibrium, 5)+ AqEquilibriaList(WhichEquilibrium, 23))
	END IF


	
	RETURN
END FUNCTION MeanActivity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Charge Fraction (Z Fraction) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ChargeFraction (CationChg, AnionChg)

	USE Chemistry, ONLY : AqCationCharge, AqAnionCharge

	IMPLICIT NONE

	INTEGER :: CationChg, AnionChg

	ChargeFraction = (CationChg+AnionChg)*(CationChg+AnionChg)/2./CationChg/AnionChg

	RETURN
END FUNCTION ChargeFraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Q Parameter for the K-M Scheme !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CalcQ (EqRxn, InParticle)

	USE Chemistry,		 ONLY : AqCationCharge,		&
				        AqAnionCharge,		&
					HowManyAqCations,	&
					HowManyAqAnions,	&
					HowManyAqChems,		&
					KMRef

	USE ModelParameters, ONLY : ThermoTref

	IMPLICIT NONE

	INTEGER :: EqRxn
	TYPE (Particle), INTENT(IN) :: InParticle	!Patch N10

	IF (KMRef(EqRxn,6) .LE. HowManyAqChems) THEN
		CalcQ = KMRef(EqRxn,6) * (1. + KMRef(EqRxn,7) *(InParticle%Temperature-ThermoTref) /	    &
				abs(AqCationCharge(INT(KMRef(EqRxn,2)))*AqAnionCharge(INT(KMRef(EqRxn,3)))))
	ELSE
		IF(KMRef(EqRxn,6) .LE. HowManyAqChems+HowManyAqCations) THEN

			CalcQ = KMRef(EqRxn,6) * (1. + KMRef(EqRxn,7) * (InParticle%Temperature-ThermoTref) /		  &
					abs((AqCationCharge(INT(KMRef(EqRxn,2)))-AqCationCharge(INT(KMRef(EqRxn,1))-HowManyAqChems)) &
					*AqAnionCharge(INT(KMRef(EqRxn,3)))))
		ELSE

			CalcQ = KMRef(EqRxn,6) *	(1. + KMRef(EqRxn,7) * (InParticle%Temperature-ThermoTref) /	&
					abs(AqCationCharge(INT(KMRef(EqRxn,2)))														&
					*(AqAnionCharge(INT(KMRef(EqRxn,3)))-AqAnionCharge(INT(KMRef(EqRxn,1))-HowManyAqChems-HowManyAqCations))))
		END IF
	END IF

	RETURN
END FUNCTION CalcQ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Composite Activity Coefficient Using the !!
!! Kusik & Meissner approximation.			  !!
!!							  !!
!! This routine owes a good deal to one developed by	  !!
!! Tim Resch (1995)	in his thesis model at MIT.	  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE KusikMeissner (InParticle)

	USE Chemistry, ONLY : AqEquilibriaList,		&
			AqCationCharge,			&
			AqAnionCharge,			&
			HowManyAqCations,		&
			HowManyAqAnions,		&
			HowManyAqChems,			&
			HowManyAqEqReactions,		&
			HowManyKMParams,		&
			KMRef,                          &  
                        AqCationNames,AqAnionNames

	USE InfrastructuralCode, ONLY :  ERROR,         &
                                      int2str,real2str

	USE ModelParameters, ONLY : MaxIonicStrKM

	IMPLICIT NONE

	!! External Variables
	TYPE (Particle) :: InParticle

	!! Internal Variables
	INTEGER :: I, J, K, allocation_error,												&
			   CationicIndex(HowManyAqEqReactions),AnionicIndex(HowManyAqEqReactions),	&
			   CatCharge(HowManyAqEqReactions),AnCharge(HowManyAqEqReactions),			&
			   KMCatCharge(HowManyKMParams),KMAnCharge(HowManyKMParams)

    REAL*8  :: CationicStrFr(HowManyAqEqReactions),AnionicStrFr(HowManyAqEqReactions),	&
			   SUM1, SUM2, Ia, Ic, LGP1, LGP2, LGP3, II, JJ, FauxIonicStr, TempIonicStr

    REAL*8, ALLOCATABLE  :: LogGammaPure(:)

	!! If this is an empty section, we have no truck with it
	! cmb
	!IF (InParticle%NumberOfParticles .EQ. 0 .OR. InParticle%Dry) RETURN
	IF (InParticle%NumberOfParticles .le. 1.0e-40 .OR. InParticle%Dry) RETURN

	!! Update the Ionic Strength Values for the particle
	CALL CalculateIonicStrength (InParticle)
	!write(*,*) 'IonicStrength: ', inparticle%ionicstr

	!! If the ionic strength is out of range, then
	!! limit it for this calculation
	IF (InParticle%IonicStr .GT. MaxIonicStrKM) THEN
		TempIonicStr        = InParticle%IonicStr
		InParticle%IonicStr = MaxIonicStrKM 
	ELSE
		TempIonicStr = 0.
	END IF

!! Water reaction is one unless have another way to determine (like KOH, etc)
!! 

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Presort the chemical charges, !!
	!! and the ion indices,          !!
	!! so can account for bisulfate- !!
	!! like charges appropriately.   !!
	!! Just use these for inflating  !!
	!! activity coefficients		 !!
	DO I = 1, HowManyAqEqReactions !!!!

		!! If this is positive, we are to take a post-mixing result
		IF (AqEquilibriaList(I,19) .GT. 0) then
			!write(*,*) 'AqEquilibriaList(I, 19) = 0; ', i, AqEquilibriaList(i,19)
			CYCLE
		endif

	    CatCharge(I) = INT(AqCationCharge(INT(AqEquilibriaList(I,2))))
		AnCharge(I)  = INT(-1 * AqAnionCharge(INT(AqEquilibriaList(I,3))))

		CationicIndex(I) = INT(AqEquilibriaList(I,2))
		AnionicIndex(I)  = INT(AqEquilibriaList(I,3))

		!! Check to see if the electrolyte is an ion
		IF (AqEquilibriaList(I,1) .GT. HowManyAqChems) THEN
			IF (AqEquilibriaList(I,1) .LE. HowManyAqChems + HowManyAqCations) THEN
				CationicIndex(I) = INT(AqEquilibriaList(I,1))-HowManyAqChems
				CatCharge(I) = CatCharge(I) - AqCationCharge(CationicIndex(I))
			ELSE
				AnionicIndex(I)  = INT(AqEquilibriaList(I,1))-HowManyAqChems-HowManyAqCations
				AnCharge(I)  = AnCharge(I)  + AqAnionCharge(AnionicIndex(I))
			END IF
		END IF
	END DO

	!! Do another charge index for the KM Equations, rather than the
	!! equilibria for the model
	DO I = 1, HowManyKMParams 

	    KMCatCharge(I) = INT(AqCationCharge(INT(KMRef(I,2))))
		KMAnCharge(I)  = INT(-1 * AqAnionCharge(INT(KMRef(I,3))))

		!! Check to see if the electrolyte is an ion
		IF (KMRef(I,1) .GT. HowManyAqChems) THEN
			IF (KMRef(I,1) .LE. HowManyAqChems + HowManyAqCations) THEN
				KMCatCharge(I) = KMCatCharge(I) - AqCationCharge(INT(KMRef(I,1))-HowManyAqChems)
			ELSE
				KMAnCharge(I)  = KMAnCharge(I)  + AqAnionCharge(INT(KMRef(I,1))-HowManyAqChems-HowManyAqCations)
			END IF
		END IF
	END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The K-M logG values for the pure substances are calculated with the	     !!
!! Meissner equation that resembles:					     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! logGamma = (log[1+B(1+0.1I)^q-B] + -0.5107*sqrt(I)/(1.0+c*sqrt(I))	     !!
!!  where b = 0.75-0.065*q						     !!
!!		  c = 1.0 + 0.055*q*exp(-0.023 I^3)			     !!
!!		  q = qr * (1. + qt*(T-T0)/zx zy)			     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! And stored temporarily in an array					     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ALLOCATE (LogGammaPure(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of LogGammaPure Failed in KusikMeissner()")

	!! Initialize this in case we don't specify it (needed in the summations)
	DO I = 1, HowManyAqEqReactions
		LogGammaPure(I) = 0.
	END DO

	!! then we shift to consider the reactions only, which may or may not
	!! include all possible pairings, as they need to be user defined and
	!! underspecifying results in a warning, not an error.  User is warned
	!! if (s)he doesn't specify as many reactions as there may be.
	!! This has implications for effective mixing of the activities as well
	!! as the species that are 
	DO K = 1, HowManyAqEqReactions

		!! If this is positive, we are to take a post-mixing result
		IF (AqEquilibriaList(K,19) .GT. 0) CYCLE

		!! Only calculate an activity via composite when AEL(K,18) = 0
		!! 17 needn't equal zero in that case, since it could be positive
		!! if we have a partially-dissociating reaction, which is dereferenced
		!! during equilibration
		! CMB (AER): floating point inequality
		if (dabs(aqequilibrialist(K, 18)) .le. 1.0e-40) then
		!IF (AqEquilibriaList(K,18) .EQ. 0.) THEN

			!! If only 16 is filled, then this is a simple case and we calculate
			!! the activity based on that indexed value.  If 16 and 17 but not 18
			!! are filled, then we have an ion dissociating into other ions (such
			!! as bisulfate) in which case the activity index for the actual reaction
			!! is in 17, as it will eventually be used in the denominator of the 
			!! equilibration calculation.
			! CMB (AER): floating point inequality
			IF (dabs(AqEquilibriaList(K,17)) .le. 1.0e-40) THEN
			!IF (AqEquilibriaList(K,17) .EQ. 0) THEN
				LogGammaPure(K) = CalcLogGammaPure(FLOOR(AqEquilibriaList(K,16)), InParticle)
			ELSE
				LogGammaPure(K) = CalcLogGammaPure(FLOOR(AqEquilibriaList(K,17)), InParticle)
			END IF
		
		ELSE

			!! If it is a composite activity coefficient, then we treat the composite
			!! as if it were the actual activity coefficient, for mixing purposes and
			!! the rest.

			I = FLOOR(AqEquilibriaList(K,16))
			JJ = 10000.*(AqEquilibriaList(K,16)-FLOOR(AqEquilibriaList(K,16)))
			LGP1 = CalcLogGammaPure(I, InParticle)*(KMCatCharge(I)+KMAnCharge(I))*JJ

			I = FLOOR(AqEquilibriaList(K,18))
			JJ = 10000.*(AqEquilibriaList(K,18)-FLOOR(AqEquilibriaList(K,18)))
			LGP2 = CalcLogGammaPure(I, InParticle)*(KMCatCharge(I)+KMAnCharge(I))*JJ

			I = FLOOR(AqEquilibriaList(K,17))
			JJ = 10000.*(AqEquilibriaList(K,17)-FLOOR(AqEquilibriaList(K,17)))
			LGP3 = CalcLogGammaPure(I, InParticle)*(KMCatCharge(I)+KMAnCharge(I))*JJ

			LogGammaPure(K) = (LGP1+LGP2-LGP3)*(1/(AqEquilibriaList(K,4)+AqEquilibriaList(K,5)))  &
							  /(CatCharge(K)+AnCharge(K))

		END IF
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Then we find the results for a mixed particle	!!
	!! using the mixing rule of Patwardhan and Kumar:	!!
	!!													!!
	!! ------------------------------------------------	!!
	!!               z1      xj (z1+zj)^2				!!
	!! logG(mix) = -----*Sum[------------ logG(1j),{j}]	!!
	!!             z1+z2        2 z1 zj					!!
	!!													!!
	!!               z2      xi (zi+z2)^2				!!
	!!           + -----*Sum[------------ logG(i2),{i}]	!!
	!!             z1+z2        2 zi z2					!!
	!! ------------------------------------------------	!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Loop over all of the electrolytes and determine the mixed coefficients
	DO K = 1, HowManyAqEqReactions

	  !! If this is positive, we are to take a post-mixing result
	  IF (AqEquilibriaList(K,19) .GT. 0) CYCLE

	  !! Reset FauxIonicStr.  For cases where all of the pairs of ions are
	  !! specified as reactions, this is overkill computationally.  But it 
	  !! is necessary for the cases in which it is not.
	  FauxIonicStr = 0.

	  !! No need to calculate the mixture if the ionic strength is zero! (Which it oughtn't ever be.)
	  ! CMB (floating point inequality check)
	  IF (dabs(InParticle%IonicStr) .le. 1.d-40) THEN
		  InParticle%GammaMixed(K) = 1.
	  ELSE

	
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !! Use the mixed activity coefficient scheme from Kusik and Meissner, 1978 !!
	  !! modified slightly to deal with partial dissociation.                    !!

	  SUM1 = 0. 
	  DO J = 1, HowManyAqEqReactions

		!! If this is positive, we are to take a post-mixing result
		IF (AqEquilibriaList(J,19) .GT. 0) CYCLE

		IF (CationicIndex(K) .EQ. CationicIndex(J))	THEN
			SUM1 = SUM1 + ChargeFraction(CatCharge(K),AnCharge(J))*LogGammaPure(J)*				&
				   AnionMolality (AnionicIndex(J), InParticle) * AnCharge(J) * AnCharge(J) / 2.
	
			FauxIonicStr = FauxIonicStr + AnCharge(J)*AnCharge(J)*AnionMolality(AnionicIndex(J),InParticle)
		END IF
	  END DO

	  SUM2 = 0.
	  DO J = 1, HowManyAqEqReactions

		!! If this is positive, we are to take a post-mixing result
		IF (AqEquilibriaList(J,19) .GT. 0) CYCLE

		IF (AnionicIndex(K) .EQ. AnionicIndex(J)) THEN
			SUM2 = SUM2 + ChargeFraction(CatCharge(J),AnCharge(K))*	LogGammaPure(J)*				&
				   CationMolality (CationicIndex(J), InParticle) * CatCharge(J) * CatCharge(J) / 2.

			FauxIonicStr = FauxIonicStr + CatCharge(J)*CatCharge(J)*CationMolality(CationicIndex(J),InParticle)
		END IF
	  END DO

	  FauxIonicStr = FauxIonicStr / 2.


	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !! The 10.** and the leading Zi*Zj come from inverting from !!
	  !! reduced activity coefficient to normal for storage       !!
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !! Changed this to FauxIonicStr from real one calculated at !!
	  !! beginning of the equilibration in July, 2003 !!!!!!!!!!!!!!
	  IF (FauxIonicStr .EQ. 0) THEN
		  InParticle%GammaMixed(K) = 0.
	  ELSE
		  InParticle%GammaMixed(K) = 10.**(CatCharge(K)*AnCharge(K)*					&
									 (CatCharge(K)*SUM1+AnCharge(K)*SUM2)				&
									 /(CatCharge(K) + AnCharge(K))/ FauxIonicStr)		

	  END IF
	  END IF
	END DO

	!! Loop over all of the electrolytes and transfer the mixed coefficients if necessary
	DO K = 1, HowManyAqEqReactions
	  IF (AqEquilibriaList(K,19) .GT. 0) THEN
		IF (AqEquilibriaList(K,22) .GT. 0) THEN 
			!Levitocite	
			InParticle%GammaMixed(K) = ( (InParticle%GammaMixed(INT(AqEquilibriaList(K,16)))**3.0) * InParticle%GammaMixed(INT(AqEquilibriaList(K,17))) ) **(1.0/4.0)	
		ELSE
			!Bisulfate
			InParticle%GammaMixed(K) = InParticle%GammaMixed(INT(AqEquilibriaList(K,19)))
		END IF
	 END IF
	END DO


	DEALLOCATE(LogGammaPure) !! This array was only to be used in this calculation routine

	!! Replace the ionic strength
	IF (TempIonicStr .GT. 0) THEN
		InParticle%IonicStr = TempIonicStr 
	END IF

	RETURN
END SUBROUTINE KusikMeissner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Log Gamma Pure for a particular !!
!! reaction indexed in KMRef.		     !!
!! This is simply the Kusik-Meissner         !!
!! formulation				     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL *8 FUNCTION CalcLogGammaPure (K, InParticle)
     
        USE InfrastructuralCode, ONLY : IsNaN

	IMPLICIT NONE

	!! External Variables
	INTEGER :: K

	!! Internal Variables
	REAL *8 :: B, C, Q, SqrtI
	TYPE (Particle) :: InParticle

	!! If this is an empty section, we have no truck with it
        ! CB: Another floating-point modification adjustment...
	!IF (InParticle%NumberOfParticles .EQ. 0 .OR. InParticle%Dry) RETURN
        if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return

	IF (K .EQ. 0) THEN
		CalcLogGammaPure = 0.
		RETURN
	END IF

	! CMB: expand upon this a bit
	IF(IsNAN(InParticle%IonicStr)) then
		WRITE(*,*) "Problem in CalcLogGammaPure!"
        write(*,*) "Ionic Strength is NaN for this particle."
		write(*,*) "Water Concentration: ", inparticle%aqchems(1)
	endif
	
	!! This is the temperature corrected Kusik-Meissner Coefficient
	Q = CalcQ(K, InParticle)
	B = 0.75-0.065*Q
	C = 1.0 + 0.055*Q*EXP(-0.023*InParticle%IonicStr**3.)
	SqrtI = SQRT(InParticle%IonicStr)

	!! Give the log of the reduced activity coefficient for the pure electrolyte
	CalcLogGammaPure = LOG10(1.+B*(1.+0.1*InParticle%IonicStr)**Q-B) - 0.5107*SqrtI/(1.0+SqrtI*C)

	!! Rounding Errors are normal when Ionic Strength is zero and this value goes to zero as well.
	IF (ABS(CalcLogGammaPure) .LE. 10.**(-12)) CalcLogGammaPure = 0.

	RETURN 
END FUNCTION CalcLogGammaPure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Water Activity of a Particular Aerosol Droplet !!
!! of mixed electrolyte composition.				!!
!!								!!
!! The PerturbWater input adds that amount of water (a multiple !!
!! to the  particle, and is used to get a sense of the derivate !!
!! of the residual water activity.				!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Patch N11.  The integral from 0 to I of EvalIdLogGamma depends only on the   !!
!! reaction's (Q, Stoich, Charges), i.e. on reaction and temperature.  Instead  !!
!! of an adaptive DGAUS8 integration on every call, tabulate the cumulative     !!
!! integral in s = sqrt(I) (where the integrand is analytic), using 8-point     !!
!! Gauss-Legendre per cell, and evaluate by cubic Hermite interpolation using   !!
!! the exact derivative.  The table covers s in [0, SHi]; SHi grows on demand   !!
!! (up to sqrt(MaxIonicStrKM)).  Not thread-safe.                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GibbsDuhemIntegralTabulated (K, Q, Stoich, Charges, IonicStrength, Result, IERR)

	USE ModelParameters, ONLY : WaterContentPrecision, MaxIonicStrKM, GDTableCells
	USE Chemistry,       ONLY : HowManyAqEqReactions

	IMPLICIT NONE

	INTEGER :: K, IERR
	REAL*8  :: Q(3), Stoich(4), Charges(4), IonicStrength, Result

	EXTERNAL EvalIdLogGamma
	REAL*8   :: EvalIdLogGamma

	REAL*8, PARAMETER :: GLX(4) = (/ 0.1834346424956498D0, 0.5255324099163290D0, &
	                                 0.7966664774136267D0, 0.9602898564975363D0 /)
	REAL*8, PARAMETER :: GLW(4) = (/ 0.3626837833783620D0, 0.3137066458778873D0, &
	                                 0.2223810344533745D0, 0.1012285362903763D0 /)

	REAL*8  :: Key(11), S, SHi, SMax, H, T, T2, T3, X, Sum, Gm, FJ, FJ1, GJ, GJ1
	INTEGER :: J, L, N

	N = GDTableCells
	IERR = 0
	SMax = SQRT(MaxIonicStrKM)

	IF (.NOT. ALLOCATED(GDTabOK)) THEN
		ALLOCATE (GDTabOK(HowManyAqEqReactions), GDTabBad(HowManyAqEqReactions), GDTabSHi(HowManyAqEqReactions))
		ALLOCATE (GDTabF(0:N,HowManyAqEqReactions), GDTabG(0:N,HowManyAqEqReactions), GDTabKey(11,HowManyAqEqReactions))
		GDTabOK = .FALSE.
		GDTabBad = .FALSE.
		GDTabSHi = 0.0D0
	END IF

	S = SQRT(IonicStrength)
	Key(1:3) = Q ; Key(4:7) = Stoich ; Key(8:11) = Charges

	!! A parameter change (temperature) invalidates the table
	IF (GDTabOK(K) .OR. GDTabBad(K)) THEN
		IF (ANY(Key .NE. GDTabKey(:,K))) THEN
			GDTabOK(K) = .FALSE.
			GDTabBad(K) = .FALSE.
		END IF
	END IF

	IF (.NOT. GDTabBad(K)) THEN
		!! (Re)build if missing or too short
		IF (.NOT. GDTabOK(K) .OR. S .GT. GDTabSHi(K)) THEN
			IF (GDTabOK(K)) THEN
				SHi = MIN(SMax, MAX(2.0D0*GDTabSHi(K), 1.25D0*S))
			ELSE
				SHi = MIN(SMax, MAX(1.0D0, 1.25D0*S))
			END IF
			H = SHi / N
			GDTabF(0,K) = 0.0D0
			GDTabG(0,K) = 0.0D0
			DO J = 1, N
				Sum = 0.0D0
				DO L = 1, 4
					X = (J-1)*H + 0.5D0*H*(1.0D0 - GLX(L))
					Sum = Sum + GLW(L) * 2.0D0*X*EvalIdLogGamma(X*X, Q, Stoich, Charges)
					X = (J-1)*H + 0.5D0*H*(1.0D0 + GLX(L))
					Sum = Sum + GLW(L) * 2.0D0*X*EvalIdLogGamma(X*X, Q, Stoich, Charges)
				END DO
				GDTabF(J,K) = GDTabF(J-1,K) + 0.5D0*H*Sum
				X = J*H
				GDTabG(J,K) = 2.0D0*X*EvalIdLogGamma(X*X, Q, Stoich, Charges)
			END DO
			GDTabKey(:,K) = Key
			GDTabSHi(K) = SHi
			GDTabOK(K) = .TRUE.
			!! A non-finite table (integrand invalid somewhere in range) is never used
			IF (ANY(GDTabF(:,K) .NE. GDTabF(:,K)) .OR. ANY(ABS(GDTabF(:,K)) .GT. 1.0D100)) THEN
				GDTabOK(K) = .FALSE.
				GDTabBad(K) = .TRUE.
			END IF
		END IF
	END IF

	!! A NaN or negative ionic strength (e.g. zero-water particle) must not index the table: DGAUS8 handles it as the original did
	IF (GDTabOK(K) .AND. .NOT. GDTabBad(K) .AND. S .GE. 0.0D0) THEN
		H = GDTabSHi(K) / N
		J = MIN(INT(S/H), N-1)
		T = (S - J*H)/H
		T2 = T*T ; T3 = T2*T
		FJ = GDTabF(J,K) ; FJ1 = GDTabF(J+1,K) ; GJ = GDTabG(J,K) ; GJ1 = GDTabG(J+1,K)
		Result = (2.0D0*T3 - 3.0D0*T2 + 1.0D0)*FJ + (T3 - 2.0D0*T2 + T)*H*GJ &
		       + (-2.0D0*T3 + 3.0D0*T2)*FJ1 + (T3 - T2)*H*GJ1
	ELSE
		CALL DGAUS8 (EvalIdLogGamma, 0., IonicStrength, Q, Stoich, Charges, &
		             WaterContentPrecision, Result, IERR)
	END IF

	RETURN
END SUBROUTINE GibbsDuhemIntegralTabulated

SUBROUTINE UpdateWaterActivity (InParticle) 

	USE Chemistry,		ONLY : HowManyAqEqReactions,     &
					HowManyAqChems,		 &
					HowManyAqCations,	 &
					HowManyAqAnions,	 &
					AqCationCharge,		 &
					AqAnionCharge,		 &
					AqEquilibriaList,	 &
					KMRef

	USE ModelParameters,		ONLY : WaterContentPrecision, &
						MaxIonicStrKm, GDTableMode
	USE InfrastructuralCode,	ONLY : WARN, ERROR

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER   :: InParticle

	!! Internal Variables
	INTEGER :: K, IERR, CatIndex, AnIndex
	REAL*8  :: LogAw, Q(3),	Stoich(4), Charges(4),	&
			   InorganicIntegral,					&	! The Value of the Integrated Value of the Inorganic Portion
			   II, JJ,								&
			   PatKumWeights(HowManyAqEqReactions,3), &
			   IonicStrength

	EXTERNAL EvalIdLogGamma
	REAL*8 :: EvalIdLogGamma

	!! If this is an empty section, we have no truck with it
	IF (InParticle%NumberOfParticles .EQ. 0 .OR. InParticle%Dry) RETURN

	!! Update the ionic strength
	CALL CalculateIonicStrength(InParticle)

	!! This is the counter for the logorhythm of the Water Activity, to be added to within the loop
	LogAw = 0.

	!! Get the weightings for the Patwardhan and Kumar mixing rule
	PatKumWeights = GetPatwardhandAndKumarWeightings (InParticle)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Water activity depends on each of the allowed ionic interactions. !!
	DO K = 1, HowManyAqEqReactions

		IF(AqEquilibriaList(K,22) .NE. 0) THEN
			!Skip Levitocite
			CYCLE
		END IF
		
		
		!! The integrated form of the Gibbs Duhem equation depends on q.  Compute
		!! the temperature-corrected value now.
		IF (AqEquilibriaList(K,18) .EQ. 0) THEN	
			IF (AqEquilibriaList(K,17) .EQ. 0) THEN			! a simple reaction
				Q(1) = CalcQ(FLOOR(AqEquilibriaList(K,16)), InParticle)
				Q(2) = 0 ; Q(3) = 0
				Stoich(1)  = 0 ; Stoich(2)  = 0 ; Stoich(3)  = 0 ; Stoich(4)  = 0
				Charges(1) = 0 ; Charges(2) = 0 ; Charges(3) = 0 ; Charges(4) = 0
			ELSE											! a dissociating ion
				Q(1) = CalcQ(FLOOR(AqEquilibriaList(K,17)), InParticle)
				Q(2) = 0 ; Q(3) = 0
				Stoich(1)  = 0 ; Stoich(2)  = 0 ; Stoich(3)  = 0 ; Stoich(4)  = 0
				Charges(1) = 0 ; Charges(2) = 0 ; Charges(3) = 0 ; Charges(4) = 0
			END IF
		ELSE												! it is a composite activity coefficient

			Q(1) = CalcQ(FLOOR(AqEquilibriaList(K,16)), InParticle)
			Q(2) = CalcQ(FLOOR(AqEquilibriaList(K,18)), InParticle)
			Q(3) = CalcQ(FLOOR(AqEquilibriaList(K,17)), InParticle)

			!! The definitions of all of these 

			Stoich(1) = ANINT((AqEquilibriaList(K,16) - FLOOR(AqEquilibriaList(K,16))) * 100000000.)/10000.
			Stoich(2) = ANINT((AqEquilibriaList(K,18) - FLOOR(AqEquilibriaList(K,18))) * 100000000.)/10000.
			Stoich(3) = ANINT((AqEquilibriaList(K,17) - FLOOR(AqEquilibriaList(K,17))) * 100000000.)/10000.
			Stoich(4) = AqEquilibriaList(K,4) + AqEquilibriaList(K,5)

			Charges(1) = AqCationCharge(INT(KMRef(FLOOR(AqEquilibriaList(K,16)),2))) * &
										abs(AqAnionCharge(INT(KMRef(FLOOR(AqEquilibriaList(K,16)),3))))
			Charges(2) = AqCationCharge(INT(KMRef(FLOOR(AqEquilibriaList(K,18)),2))) * &
										abs(AqAnionCharge(INT(KMRef(FLOOR(AqEquilibriaList(K,18)),3))))
			Charges(3) = AqCationCharge(INT(KMRef(FLOOR(AqEquilibriaList(K,17)),2))) * &
										abs(AqAnionCharge(INT(KMRef(FLOOR(AqEquilibriaList(K,17)),3))))
			Charges(4) = AqCationCharge(FLOOR(AqEquilibriaList(K,2))) * &
										abs(AqAnionCharge(FLOOR(AqEquilibriaList(K,3))))
			!IF(K .EQ. 30) THEN
				!WRITE(*,*) Stoich(1), Stoich(2), Stoich(3), Stoich(4)
				!WRITE(*,*) Charges(1), Charges(2), Charges(3), Charges(4)
				!WRITE(*,*) PatKumWeights(K,1), PatKumWeights(K,2), PatKumWeights(K,3) 
			!END IF
		END IF

		!!If the electrolyte is not present, cycle (Matt Alvarado, 6/30/06)
		IF (PatKumWeights(K,1) .LE. 0.0) CYCLE
		
		!WRITE(*,*) K, PatKumWeights(K,1)
		
		!Don't let ionic strength grow too large
		IF(InParticle%IonicStr .GT. MaxIonicStrKM) THEN
			IonicStrength = MaxIonicStrKM
		ELSE
			IonicStrength = InParticle%IonicStr
		END IF
		
		!! This routine is housed in DGAUS8.f.  It's an adaptive time step integrator, modified to pass dummy arguments to FUN.
		!Patch N11: tabulated Gibbs-Duhem integral (falls back to DGAUS8 if the table is unusable)
	    IF (GDTableMode .EQ. 1) THEN
	        CALL GibbsDuhemIntegralTabulated (K, Q, Stoich, Charges, IonicStrength, InorganicIntegral, IERR)
	    ELSE
	        !! Input-deck switch: full adaptive integration on every call (original behaviour)
	        CALL DGAUS8 (EvalIdLogGamma, 0., IonicStrength, Q, Stoich, Charges, &
	                     WaterContentPrecision, InorganicIntegral, IERR)
	    END IF


		!! Consider the Flagged Result of the Integration
		!! 1  is normal return
		!! -1 is a normal return of 0 given because the bounds of the finite integral are so close together.
		!! 2  is an "abnormal reult" (aka failure)
		IF (IERR .EQ. 2) THEN
			CALL DumpParticleContentsAtError(InParticle, "FailedWaterActivityIntegratingParticle.txt", 	&
						InRelativeHumidity=.TRUE., InPh=.TRUE., InDoSurfaceTension=.TRUE., &
                                                 InDensity=.TRUE.,InDoWaterActivity=.TRUE.,	&
						InRadius = .TRUE., InIonicStrength = .TRUE.)
			CALL ERROR("Integration of the Gibbs-Duhem Inorganic Integral Failed in UpdateWaterActivity(), a "//	&
					   "LagrangianAerosols.f90 / Thermodynamics.h Subroutine.")
		END IF

		!! 0.0360376 is Mw / 500, 2.3025... is ln(10.)...
		LogAw = LogAw - 0.0360376*PatKumWeights(K,1)*(InParticle%IonicStr/abs(2.302585093*				&
				AqCationCharge(INT(PatKumWeights(K,2)))*AqAnionCharge(INT(PatKumWeights(K,3)))) + InorganicIntegral)

		!!!
		!!!
		!!!
		!!!
		!!! THIS DOESN"T INCLUDE NON-ELECTROLYTIC VALUES AT ALL!  THERE SHOULD BE ANOTHER INTEGRAL
		!!!
		!!!
		!!!
		!!!

	END DO

	InParticle%WaterActivity = 10.**(LogAw)

	RETURN
END SUBROUTINE UpdateWaterActivity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Residual Between the Ambient RH and the Local	!!
!! Water Activity corrected for curvature effects				!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION WaterResidual (InParticle, SpecifiedRH)

	USE GridPointFields,     ONLY : GetRelativeHumidity
	USE infrastructuralcode, ONLY : ERROR

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8       :: SpecifiedRH

	!! Internal Variables
	REAL*8 :: RH

	CALL UpdateWaterActivity(InParticle)

!	IF (PRESENT(SpecifiedRH) .AND. SpecifiedRH .NE. -1.) THEN
	IF (SpecifiedRH .NE. -1.) THEN
		RH = SpecifiedRH
	ELSE
		RH = GetRelativeHumidity ()
	END IF

	!! The curvature effects are proscribed by a simple equation from 
	!! Pruppacher and Klett:
	!!  RH = a_w * exp(2 Mw surftens_s/a / RT density r)
	WaterResidual = InParticle%WaterActivity*CurvatureCorrection(InParticle) &
					- RH

	RETURN
END FUNCTION WaterResidual

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Use Bisection to Match the ambient relative humidity                 !!
!! to the local water activity.					        !!
!!								        !!
!! ForceWaterEquilibrium will make the routine use the equilibrium      !!
!! approach no matter whether the RH is above the threshold that        !!
!! would tell the routine to use the integration approach.              !!
!! MATT'S NOTE: This differs from EquilibriumWaterContent in two ways:  !!
!! 1. This returns the total inorganic water per particle as a variable,!! 
!!  rather than adjusting the water content of the particle directly.   !!
!! 2. This never adjusts the input RH; As such, it simply callculates   !!
!!	   what water content would be in equilibrium with the input    !!
!!     relative humidity.					        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibriumWaterContentAmount (InParticle, EquilibToThisRH, &
                            ReturnType, ForceEquilibrium, WaterContent)

	USE Chemistry,ONLY : HowManyAqCations, HowManyAqAnions, HowManyAqChems

	USE ModelParameters,     ONLY : AqThermoNumbAllowedIterations,	&
					WaterContentPrecision,		&
					AerosolWaterEquilibriumRHThreshold,  &
					OpenSystemForCondensation, Avogadro, &
					SmallestAerosolPossible,            &
                                        MinimumWater

	USE InfrastructuralCode, ONLY : Error, IsNaN

	USE GridPointFields,	 ONLY : GetRelativeHumidity, GetGasBurden, &
					ReplaceGridCellChemBurden

	USE Time,		ONLY : CurrentTime, BeginningTime

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8	     :: EquilibToThisRH
	LOGICAL      :: ForceEquilibrium
	INTEGER      :: ReturnType
	REAL*8	     :: WaterContent
	

	!! Internal Variables 
	REAL*8  :: InitialWaterContent, Residual, Change, ChangeFac, II, RH
	INTEGER :: I, N, Direction
	LOGICAL :: InternalOpenSystem, InternalForceEquilibrium

	!! If this is an empty section, we have no truck with it
	IF (InParticle%NumberOfParticles .EQ. 0) RETURN

	!IF (PRESENT(ReturnType)) 
        ReturnType = 2

	!Set RH level to input value
	RH = EquilibToThisRH
	
	!Assume an open system
	InternalOpenSystem = .TRUE.
	

	InternalForceEquilibrium = ForceEquilibrium

	!Store Initial Water Content - we'll restore to this later
	!NOTE (Matt Alvarado, 6/06): What I've done here is allowed 
	!the particle water content
	!to follow the interation, to fit with how the functions
	!for water activity and size are written. I then
	!restore their initial values at the end, and return the 
	!calculated equilibrium inorganic water content. 
	InitialWaterContent = InParticle%AqChems(1)
	
	!Initialize water content
	WaterContent = InParticle%AqChems(1)

	!! Decide on a direction and loop until the correct order of magnitude
	InParticle%AqChems(1) = WaterContent
        
	CALL RecalculateRadius (InParticle)
	Residual = WaterResidual (InParticle, RH)


	!! If the particle is already in equilibrium, exit
	IF (ABS(Residual) .LE. WaterContentPrecision) RETURN

	!IF (PRESENT(ReturnType)) 
        ReturnType = 1

	IF (IsNaN(Residual)) THEN
		CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle.txt",InDoSurfaceTension=.FALSE., &
                      InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
		      InDensity=.TRUE., InIonicStrength=.TRUE.)
		CALL ERROR("In the initial call to WaterResidual, EquilibriumWaterContentAmount()"// &
                           " in Thermodynamics.h, the water residual "// &
		           "calculation returned NaN.  Something is quite wrong.")
	END IF

	IF (Residual .GT. 0) THEN
		Direction = -1
		ChangeFac = 0.5
	ELSE
		Direction = 1
		ChangeFac = 2.
	END IF

	!! Overshoot the correct value by one iteration (of doubling or halving)
        !WRITE(*,*) "Before First Loop"
	N = 1
	DO WHILE (Residual * Direction .LE. 0.) 

		WaterContent = WaterContent * ChangeFac 

		!Calculate Residual
		InParticle%AqChems(1) = WaterContent
                !WRITE(*,*) " Before Rad in first loop:  ",InParticle%AqChems(1) 
		CALL RecalculateRadius (InParticle)
		Residual = WaterResidual (InParticle, RH)

                !If particle is shrinking below the minimum water, stop loop 
                !and return as if equilibrated
                IF(Residual .GT. 0.0 .AND. InParticle%AqChems(1) .LT. MinimumWater) THEN
                     ReturnType = 2
                     RETURN
                ENDIF
		
		IF (IsNaN(Residual)) THEN
			CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle.txt",InDoSurfaceTension=.FALSE., &
                      InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
		      InDensity=.TRUE., InIonicStrength=.TRUE.)
			CALL ERROR("In the initial search loop of EquilibriumWaterContentAmount() in Thermodynamics.h, the water residual "// &
					   "calculation returned NaN.  Something is quite wrong.")
		END IF

		IF (N .GT. 10*AqThermoNumbAllowedIterations) THEN
			CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle.txt",InDoSurfaceTension=.FALSE., &
                      InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
		      InDensity=.TRUE., InIonicStrength=.TRUE.)
			CALL ERROR ("Took too many iterations to find the right order of magnitude for water content in "// &
			            "first loop of EquilibriumWaterContentAmount().  Increase AqThermoNumbAllowedIterations in ModelParameters.f90 "// &
				    "if you think nothing is wrong.  Perhaps this is called for too high of a relative humidity "// &
				    "(won't equilibrate if RH > S').")
		END IF

		N = N + 1
                !WRITE(*,*) "First Loop: ", N

	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Move on to a secondary loop that refines our search !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Define Change initially as something that will retreat a full step   !!
!! in less than infinite steps of dividing in half and walking (as      !!
!! simply dividing the last step by two would)				!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (ChangeFac .EQ. 0.5) THEN 
		Change = WaterContent
	ELSE 
		Change = WaterContent/2.
	END IF

        !WRITE(*,*) "Before Second Loop"

	!! Converge on the correct water content
	N = 1
	DO
		!! Reduce the Step Size
		Change = Change / 2.

		!! Recalculate the Residual and Loop Again
		InParticle%AqChems(1) = WaterContent
		CALL RecalculateRadius (InParticle)
		Residual = WaterResidual (InParticle, RH)
		!WRITE(*,*) Residual, Change, InParticle%AqChems(1)

                !If particle is shrinking below the minimum water, stop loop 
                !and return as if equilibrated
                IF(Residual .GT. 0.0 .AND. InParticle%AqChems(1) .LT. MinimumWater) THEN
                     ReturnType = 2
                     RETURN
                ENDIF

		!! If the appropriate precision has been reached, then return the solution
		IF (ABS(Residual) .LE. WaterContentPrecision/10) EXIT

		IF (IsNaN(Residual)) THEN
			CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle.txt",InDoSurfaceTension=.FALSE., &
                      InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
		      InDensity=.TRUE., InIonicStrength=.TRUE.)
			CALL ERROR("In the main search loop of EquilibriumWaterContentAmount() in Thermodynamics.h, the water residual "// &
					   "calculation returned NaN.  Something is quite wrong.")
		END IF

		!! Step in the correct direction
		IF (Residual .GT. 0.) THEN
			WaterContent = WaterContent - Change
		ELSE 
			WaterContent = WaterContent + Change
		END IF


		!! If the particle is smaller than a threshold radius, this equilibration becomes
		!! increasingly difficult and decreasingly important.  So just give the old 
		!! college try and then return.
		IF ((N .GT. AqThermoNumbAllowedIterations/10.) .AND. (InParticle%EmbryoRadius .LT. SmallestAerosolPossible)) RETURN

		!!Just return if change approaches 0 and you are fairly close to eq.
		IF(Change .EQ. 0.0 .AND. ABS(Residual) .LE. WaterContentPrecision*10) EXIT
		
		!! Keep count of the number of loops and exit if too many
		IF (N .GT. AqThermoNumbAllowedIterations) THEN
		    ReturnType = 2
                    EXIT
		END IF

		N = N + 1
        !WRITE(*,*) "EquilibriumWaterContentAmount iterated ", N, " times."
	END DO
	
	!Restore particle to initial state
	InParticle%AqChems(1) = InitialWaterContent
        !WRITE(*,*) "Last Rad Call"
	CALL RecalculateRadius (InParticle)
	RETURN
END SUBROUTINE EquilibriumWaterContentAmount
