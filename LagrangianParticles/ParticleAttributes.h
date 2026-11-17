!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ParticleAttributes.h
!! Contains functions for calculating particle properties 
!! (density, radius, mass
!! curvature correction, Knudsen number, Reynolds number, etc.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 07/24  2006   Matt Alvarado     Removed "ISNAN" to fit pgf90		     !!
!! 09/07  2006   Matt Alvarado   Changed RecalculateRadius to give lower edge!!
!!				      when number conc. = 0.		     !!
!! 11/01  2006   Matt Alvarado     Created ShellRefIndAndRad
!! 05/18  2007   Matt Alvarado     Revised ShellRefIndAndRad
!!				   Created AerosolOptProp
!! 07/13  2007   Matt Alvarado     Added line to  ShellRefIndAndRad
!!	 			      to keep core rad below shell rad
!! 07/17  2007   Matt Alvarado     Added line to SurfaceTension
!!				      to keep value above that of pure Ethanol
!! 07/18  2007   Matt Alvarado     Removed calls to ParticleDensity from
!!				CurvatureCorrection and OrgCurvatureCorrection
!!			 Just use values in particle data structure instead
!!				   In ParticleDensity, if ion loop fails, force
!!				      density of ionic solution to 1.0 g/cm3
!! 07/26  2007   Matt Alvarado     Fixed AbsCoreRad calculation in 
!!                                 ShellRefIndAndRad
!! 09/07  2007   Matt Alvarado     Fixed shell refractive index calcultion
!! 09/17  2007   Matt Alvarado     Removed write statement from ParticleDensity
!! 10/01  2007   Matt Alvarado     Set ParticleDensity to skip bins 
!!                                 with 0 particles
!! 10/03  2007   Matt Alvarado     Set Denomsum limit in SurfaceTension to 
!!                                 1.0e40
!! 10/04  2007   Matt Alvarado     Protected ShellRefIndAndRad from near zero 
!!                                  number concentrations
!! 10/11  2007   Matt Alvarado     Set HypotheticalElectrolyteConcentrations 
!!                                 to warn of charge
!!                                 imbalance, but not stop.
!! 10/15  2007   Matt Alvarado     Fixed Assymetry parameter in 
!!                                 ShellRefIndAndRad
!! 08/30  2010   Matt Alvarado     Changed ShellRefIndandRad, 
!!                                 AerosolOpticalDepth, and AerosolOptProp 
!!                                 to allow 62 photolysis wavelengths
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/27  2012   Matt Alvarado     Changed ShellRefIndandRad, 
!!                                  AerosolOpticalDepth, and AerosolOptProp 
!!                                  to allow 451 bins between 250 nm and 700 nm
!!                                 ShellRefIndandRad now used OPAC wavelength
!!                                  dependent complex refractive indices
!! 02/29  2012   Matt Alvarado     Changed ShellRefIndandRad to use different
!!                                  wavelength dependent refractive indices
!!                                  for different electrolytes based on flags
!!                                  in AqPhaseChems.in. Also added refractive
!!                                  index for dust.
!! 05/03  2012   Matt Alvarado     Changed AerosolOptProp to calculate
!!                                  backscattering coefficient.
!! 05/08  2012   Matt Alvarado     Added special routine for dry BC only particles
!!                                  to particle density.
!! 05/03  2012   Matt Alvarado     Fixed calculation of hemispherical
!!                                  backscattering coefficient in AerosolOptProp
!  08/17  2012   Matt Alvarado     Added capability for external mixtures to ShellRefIndAndRad
!! 11/08  2012   Matt Alvarado     Added Maxwell-Garnett mixing rule to ShellRefIndAndRad
!! 01/24  2013   Matt Alvarado     Updated OrgCurvatureCorrection to give 1.0
!!                                  if particle number = 0.0
!! 07/10  2013   Matt Alvarado     Made it so only RecalculateRadius resets 
!!                                  the particle%surfacetension by calling 
!!                                  SurfaceTension, and CurvatureCorrection
!!                                  uses the particle%surfacetension
!! 10/28  2026   Claude Sonnet 5.5 Commit 8cb0991: Phase 1.1: hoist OPAC     !!
!!                                 tables to module PARAMETERs               !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 3717b42: Phase 1.2: hoist per-     !!
!!                                 category refractive-index lookups         !!
!! 10/28  2026   Claude Sonnet 5.5 Commit 7737090: Phase 1.3: MIEV0 for      !!
!!                                 homogeneous-sphere Mie calls (adds        !!
!!                                 MieSupport.f; link MIEV0.f +              !!
!!                                 MieSupport.f)                             !!
!! 11/17  2026   Claude Opus 5.5   Fast-TUV Mie table as PARAMETER MieLUT;  !!
!!                                 module OPAC tables; volume and shell/core !!
!!                                 index code shared with ShellRefIndAndRad  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE ResetAllOriginPressureLevels (OriginPressureLevel)
!! 2. FUNCTION ParticleKnudsen (InParticle)
!! 3. FUNCTION ParticleKnudsenForEnergy (InParticle)
!! 4. FUNCTION CunninghamSlipCorrectionFactor (InParticle)
!! 5. FUNCTION AerosolPH (InParticle)
!! 6. FUNCTION ParticleDensity (InParticle) !Sets radius and effective radius
!! 7. FUNCTION ParticleMass (InParticle)
!! 8. FUNCTION GetPatwardhandAndKumarWeightings (InParticle)
!! 9. FUNCTION HypotheticalElectrolyteConcentrations (InParticle)
!!10. FUNCTION TerminalVelocity (InParticle)
!!11. SUBROUTINE RecalculateRadius (InParticle) !Calls particle density
!!12. FUNCTION CurvatureCorrection (InParticle)
!!13. FUNCTION SurfaceTension (InParticle)
!!14. FUNCTION ReynoldsNumber (InParticle)
!!15. REAL*8 FUNCTION OrgCurvatureCorrection (InParticle)
!!16. REAL*8 FUNCTION MolalityAqCarbon(InParticle)
!!17. REAL*8 FUNCTION AqueousSolutionMass
!!18. SUBROUTINE ShellRefIndAndRad(InParticle)
!!19. REAL*8 FUNCTION AerosolOpticalDepth()
!!20. SUBROUTINE AerosolOptProp
!!21. SUBROUTINE OpticalVolumeConcentrations (shared by 18 and 22)
!!    SUBROUTINE ShellAndCoreIndices (shared by 18 and 22)
!!22. SUBROUTINE ShellRefIndAndRadFASTTUV (lookup table MieLUT)
!!23. SUBROUTINE AerosolOptPropFASTTUV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! Modified by CMB (AER): Gfortran complains when "Particle" is used both
!                        as a structure type declarator and a local variable
!                        name; I have replaced the variable names with
!                        "p_Particle"

SUBROUTINE ResetAllOriginPressureLevels (OriginPressureLevel)

	IMPLICIT NONE

	Type(Particle),POINTER :: p_Particle
	REAL    :: OriginPressureLevel

	p_Particle => Particles%First

	IF (ASSOCIATED(p_Particle)) p_Particle%OriginPressureLevel = OriginPressureLevel
		
	DO WHILE (ASSOCIATED(p_Particle%Next))
		p_Particle => p_Particle%Next
		IF (ASSOCIATED(p_Particle)) p_Particle%OriginPressureLevel = OriginPressureLevel
	END DO

	RETURN
END SUBROUTINE ResetAllOriginPressureLevels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Knudsen Number of the Particle	!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleKnudsen (InParticle)

	USE GridPointFields, ONLY : MeanFreePathOfAir

	IMPLICIT NONE

	Type(Particle),POINTER :: InParticle

	ParticleKnudsen = MeanFreePathOfAir ()/ InParticle%EffectiveRadius
	
	!write(*,*) 'CRL: ParticleKnudsen = ',ParticleKnudsen
	!write(*,*) 'CRL: MeanFreePathOfAir = ',MeanFreePathOfAir ()
	!write(*,*) 'CRL: InParticle%EffectiveRadius= ',InParticle%EffectiveRadius
 	!write(*,*) ''

END FUNCTION ParticleKnudsen


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Knudsen Number Of Energy of the Particle !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleKnudsenForEnergy (InParticle)

	USE GridPointFields, ONLY : MolecularThermalDiffusivity, & 
                                    GetThermalVelocityOfAir

    IMPLICIT NONE

	Type(Particle),POINTER :: InParticle

	ParticleKnudsenForEnergy = 3 *MolecularThermalDiffusivity() &
			/ GetThermalVelocityOfAir()   &
			/ InParticle%EffectiveRadius

	RETURN
END FUNCTION ParticleKnudsenForEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the CUNNINGHAM SLIP CORRECTION FACTOR, which !!
!! has a number of possible forms.  This function selects !!
!! the appropriate one and returns the value.			  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CunninghamSlipCorrectionFactor (InParticle)

    IMPLICIT NONE

	Type(Particle),POINTER :: InParticle
	REAL*8		   :: Knudsen

	Knudsen = ParticleKnudsen (InParticle)

	!! From Kasten 1968 (via Jacobson, 1997)
	CunninghamSlipCorrectionFactor = 1 + Knudsen * (1.249 + 0.42 * dexp(-0.87 / Knudsen))

	!write(*,*) 'CRL: CunninghamSlipCorrectionFactor = ',CunninghamSlipCorrectionFactor
	!write(*,*) 'CRL: Knudsen = ',Knudsen
	
	RETURN
END FUNCTION CunninghamSlipCorrectionFactor



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns the pH of an aerosol. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION AerosolPH (InParticle)

	USE ModelParameters, ONLY : ProtonIndex

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle

	!!Internal Variables
	REAL*8 :: ProtonMolarity


	IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) THEN
		AerosolPH = 14.

		RETURN
	END IF

	IF (InParticle%AqChems(ProtonIndex) .GT. 0) THEN
		!Molarity in mol/L solution
		ProtonMolarity = InParticle%AqChems(ProtonIndex)*InParticle%SolutionDensity*1000/AqueousSolutionMass(InParticle)
		
		AerosolPH = -1.*LOG10(ProtonMolarity)

	ELSE
		AerosolPH = 14.

	END IF


END FUNCTION AerosolPH



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the PARTICLE DENSITY for a lagrangian particle !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleDensity (InParticle)
!MJA 05-008-2012 Added special routine for dry BC only particles

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, &
                                        moles, pi, avogadro,    &
                                        SmallestAerosolPossible, & 
									micron ! cmb add

	USE Chemistry,		ONLY : IFSOLID, AqMolecularMass,       &
                                       HowManyAqChems,HowManyAqCations,&
				       HowManyAqAnions,                &
                                       HowManyAqEqReactions,           & 
                                       AqEquilibriaList,	       &
				       AqAnionCharge, AqCationCharge,  &
				       HowManyOrgChems,                &
                                       OrgMolecularMass, OrgDensity,   &
                                       SolidSaltDensity,               &
                                       HowManyAqOrgChems,              &
                                       AqOrgMolecularMass,             &
                                       AqOrgDensity

	USE InfrastructuralCode, ONLY : ERROR,REAL2STR, WARN, IsNaN, int2str

	! cmb add
	!use outputroutines, only : spillbeans
	! 
	
	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I,K, SkipFlag, C
	REAL*8  ::Volume,Mass,ionMass,PSIij, TempElectrolytes(HowManyAqChems),&
	BinaryDensity, NumerSum, DenomSum, TempMolality, &
  	PatKumWeights(HowManyAqEqReactions,3), XX, YY, Residual, dR, &
	SolutionDensity, InsolMass, InsolVolume, InsolDensity, InorgSolnDensity

	REAL*8, PARAMETER :: ResidualThreshold = 0.0025 !**3.

        LOGICAL :: ForceSolnDens

	real*8 :: sum_electrolytes

        !Use this flag to force the aqueous solution density to be 1.0 g/cm3 and skip a lot of code
        ForceSolnDens = .TRUE.
        
        ! CB: Modify this equality check to have a designated error 
        !     tolerance        
        !IF (InParticle%NumberOfParticles .LE. 0. ) THEN
        if (inparticle%numberofparticles .le. 1.0E-40) then
			ParticleDensity = 1.0
			InParticle%EmbryoRadius    = 0.
			InParticle%EffectiveRadius = 0.
			InParticle%InsolubleRadius = 0.
			InParticle%ParticleDensity = 1.
			InParticle%SolutionDensity = 1.
			InParticle%InsolubleDensity = 1.
			InParticle%InorgSolnDensity = 1.
			RETURN
		Else
			! cb add
			!call system("/bin/rm -f fort.2222")
			!call spillbeans(tout_unit=2222)
			!write(*,*) 'Nonzero # of particles (', inparticle%numberofparticles, ')'
			!write(*,*) ''
			!write(*, *)"Number ",trim(INT2STR(I)),",  ",inparticle%numberofparticles,inparticle%sectional
			!write(*, *)"Effective Radius: ",trim(real2str(inparticle%EffectiveRadius/micron)),",  (",trim(real2str(inparticle%Edges(1)/micron))," to ",trim(real2str(inparticle%Edges(2)/micron)),")"
			!write(*, *)"Embryo Radius: ",trim(real2str(inparticle%EmbryoRadius/micron))
			!write(*, *)"Insoluble Radius: ",trim(real2str(inparticle%InsolubleRadius/micron))
			!write(*, *)"Absorbing Core Radius: ",trim(real2str(inparticle%AbsCoreRad/micron))
			!write(*, *)"Shell Real Ref. Ind. 550 nm",trim(real2str(inparticle%ShellRealRefracInd(1)))
			!write(*, *)"Shell Imag. Ref. Ind. 550 nm",trim(real2str(inparticle%ShellImagRefracInd(1)))
			!write(*, *)"Extinction Efficiency: 550 nm",trim(real2str(inparticle%ExtEff(1)))
			!write(*, *)"Scattering Efficiency: 550 nm",trim(real2str(inparticle%ScaEff(1)))
			!write(*, *)"BackScattering Efficiency: 550 nm",trim(real2str(inparticle%BackScaEff(1)))
			!write(*, *)"Assymetry Parameter: 550 nm",trim(real2str(inparticle%AssymParam(1)))
		
			!write(*, *)"Temperature ",trim(real2str(inparticle%Temperature))
			!write(*, *)"Ionic Strength ",trim(real2str(inparticle%IonicStr))
			!write(*, *)"Water Activity ",trim(real2str(inparticle%WaterActivity))
			!write(*, *)"SolutionDensity ",trim(real2str(inparticle%SolutionDensity))
			!!write(*, *)"InsolubleDensity ",trim(real2str(inparticle%InsolubleDensity))
			!write(*, *)"ParticleDensity ",trim(real2str(inparticle%ParticleDensity))
			!write(*, *)"Ph ", trim(real2str(AerosolPh(inparticle)))
			!write(*, *)"Curvature Correction ", trim(real2str(CurvatureCorrection(inparticle)))
			!write(*, *)"OrgCurvature Correction ", trim(real2str(OrgCurvatureCorrection(inparticle)))	
			!write(*, *)"Solution Surface Tension", trim(real2str(inparticle%SurfaceTension))
			
			!write(*, *)"Temperature ", trim(real2str(inparticle%temperature))
			
			if (inparticle%numberofparticles .GT. 0.0) then
				!write(*,*) '--> calculating HypotheticalElectrolyteConcentrations <--'
				!		SaltConc = HypotheticalElectrolyteConcentrations (inparticle)
			endif   
			!write(*,*) ''
			!write(*,*) 'AqChems:'
			!write(*,*) '--------'
			!do i = 1, size(inparticle%aqchems)
		!		write(*,*) 'i = ', i, ', inparticle%aqchems(i) = ', inparticle%aqchems(i)
		!	enddo
		!	write(*,*) ''
			!write(*,*) 'OrgChems:'
			!write(*,*) '---------'
			!do i = 1, size(inparticle%orgchems)
		!		write(*,*) 'i = ', i, ', inparticle%orgchems(i) = ', inparticle%orgchems(i)
		!	enddo
		!	write(*,*) ''
		!	write(*,*) 'AqOrgChems:'
		!	write(*,*) '-----------'
		!	do i = 1, size(inparticle%aqorgchems)
		!		write(*,*) 'i = ', i, ', inparticle%aqorgchems(i) = ', inparticle%aqorgchems(i)
		!	enddo
			! end cb add
		endif

	!MJA 05-08-2012 Special routine for dry (BC only) particles
        IF (InParticle%Dry) THEN
		ParticleDensity = OrgDensity(1) !BC density
		InsolVolume = InParticle%OrgChems(1)*OrgMolecularMass(1)/OrgDensity(1) 
                InParticle%EmbryoRadius    = 0.0
		InParticle%InsolubleRadius = (0.75 * InsolVolume / pi)**0.33333333333
                InParticle%EffectiveRadius = InParticle%InsolubleRadius

		InParticle%ParticleDensity = OrgDensity(1)
		InParticle%SolutionDensity = OrgDensity(1)
		InParticle%InsolubleDensity = OrgDensity(1)
		InParticle%InorgSolnDensity = OrgDensity(1)
		RETURN
	END IF

        Mass = 0 ; Volume = 0; ionMass = 0; NumerSum = 0; DenomSum = 0
        IF (ForceSolnDens) THEN
            InorgSolnDensity = 1.0
            DO I = 1, HowManyAqChems+HowManyAqAnions+HowManyAqCations
                Mass = Mass    + AqMolecularMass(I)*InParticle%AqChems(I)
            ENDDO
            Volume = Mass*InorgSolnDensity
        ELSE        


           !! Get electrolyte equivalents of solution containing dissociated ions,
           !! which assumes a two-parameter fit proscribed in Resch, 1995 
           !! (assuming that the density of water does not change with Temperature).
           !! This is a SOLUTION density, which includes the incorporated water.
           !! So we scale up by the mass of the water as well!
           TempElectrolytes = HypotheticalElectrolyteConcentrations (InParticle)
           ! cmb in: 
           !write(*,*) 'CMB: Output from HypotheticalElectrolyteConcentrations...'
           call flush(6)
           ! end cmb
           do i = 1, size(TempElectrolytes)
              !write(*,*) 'i = ', i, ', TempElectrolytes = ', TempElectrolytes(i)
              sum_electrolytes = sum_electrolytes + TempElectrolytes(i)
           enddo
           ! end cmb in
	
           !For some reason, This routine returns NaN if there isn't a write
           !statement here. MJA, 070507
           WRITE(1234,*) TempElectrolytes(1)
	   ! cb add
           close(1234)

           !! Get the weightings for the Patwardhan and Kumar mixing rule
           !WRITE(*,*)'Just Before PatKumWeights spilling beans'
           !call spillbeans()
           PatKumWeights = GetPatwardhandAndKumarWeightings (InParticle)

           !! Loop over each reaction and consider the appropriate equation
           DO I = 1, HowManyAqEqReactions
              !write(*,*) 'i (of HowManyAqEqReactions) = ', i

              !! Some quasi-electrolytes (such as HSO4-) are not quantifyable in this
              !! density parameterization and should have been pushed to upper level
              !! electrolytes anyways (such as H2SO4).  So we ignore those electrolytes!
              ! CMB (AER, Inc): floating-point inequality check again?
              IF (AqEquilibriaList(I,1) .LE. HowManyAqChems .AND. AqEquilibriaList(I,1) .ge. 1.0e-40) THEN ! GT 0 disallows water
                 !IF (AqEquilibriaList(I,1) .LE. HowManyAqChems .AND. AqEquilibriaList(I,1) .GT. 0.) THEN ! GT 0 disallows water

                 !! This is the equivalent molality of the associated electrolyte
                 !write(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))),I
                 !Write(*,*) "Temp Electrolytes: ", TempElectrolytes(AqEquilibriaList(I,1))
                 !Write(*,*) "InParticle%AqChems(1) ", InParticle%AqChems(1)
                 !Write(*,*) "moles ", moles
                 !Write(*,*) "grams ", grams
                 !Write(*,*) "AqMolecularMass(1) ", AqMolecularMass(1)			! cmb: changed int to nint?  not sure if necessary...
                 TempMolality = TempElectrolytes(nINT(AqEquilibriaList(I,1))) / InParticle%AqChems(1) / &
                      (moles / grams * AqMolecularMass(1)) * 1000.

                 !Mass    = Mass    + AqMolecularMass(INT(AqEquilibriaList(I,1)))*TempElectrolytes(INT(AqEquilibriaList(I,1)))
                 !ionMass = ionMass + AqMolecularMass(INT(AqEquilibriaList(I,1)))*TempElectrolytes(INT(AqEquilibriaList(I,1)))
                 Mass    = Mass    + AqMolecularMass(nINT(AqEquilibriaList(I,1)))*TempElectrolytes(nINT(AqEquilibriaList(I,1)))
                 ionMass = ionMass + AqMolecularMass(nINT(AqEquilibriaList(I,1)))*TempElectrolytes(nINT(AqEquilibriaList(I,1)))

                 IF(IsNaN (TempMolality)) THEN
                    WRITE(*,*) "Problem in Temp Molality"
                    WRITE(*,*) InParticle%NumberOfParticles, I, AqEquilibriaList(I,13)
                    WRITE(*,*) TempElectrolytes(INT(AqEquilibriaList(I,1))), TempMolality, InParticle%AqChems(1)
                 END IF

                 !WRITE(*,*) TempMolality, I
                 !! Binary density is an exponential function of the two parameters input into the system
                 IF(TempMolality .LT. 100.0) THEN

                    !WRITE(*,*) "Before Exp, Rxn: ", I
                    !WRITE(*,*) "TempMolality: ", TempMolality
                    BinaryDensity = AqEquilibriaList(I,12) + EXP(-1.*AqEquilibriaList(I,13)*TempMolality)*		&
                         (1.*grams/cm/cm/cm-AqEquilibriaList(I,12))

                    !This prevents the binary density of NH4OH from going below 0.9 g/cm3
                    IF(AqEquilibriaList(I,13) .LT. 0. .AND. BinaryDensity .LT. 0.9) BinaryDensity = 0.9

                 ELSE IF(TempMolality .GE. 100.0 .AND. AqEquilibriaList(I,13) .GT. 0.) THEN
                    !If the molality gets very large, and exponent goes to zero
                    !(This keeps EXP from giving an overflow error)
                    BinaryDensity = AqEquilibriaList(I,12)
                 ELSE !For NH4OH at high molality, keep it above 0.9 g/cm3
                    BinaryDensity = 0.9	
                 END IF

                 !! Employ Mixing Rule of Patwardhan and Kumar (AIChE J. 39, 711-714) to Mix these ions appropriately
                 PSIij =  1000. * PatKumWeights(I,1) + AqMolecularMass(INT(AqEquilibriaList(I,1))) * TempMolality
                 !write(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(I,1)))
                 !write(*,*) 'PatKumWeights(I,1) = ', PatKumWeights(I,1)
                 !write(*,*) 'AqMolecularMass(INT(AqEquilibriaList(I,1))) = ',AqMolecularMass(INT(AqEquilibriaList(I,1)))	
                 !write(*,*) 'TempMolality = ',TempMolality

                 IF(IsNaN(PatKumWeights(I,1))) THEN
                    CALL WARN("Problem in PatKumWeights! Get NaN!")
                 END IF
                 !! Composite Density is the Numerator Sum divided by the Denomenator Sum, although we will track it as Volume
                 NumerSum = NumerSum + PSIij
                 DenomSum = DenomSum + PSIij / BinaryDensity
              END IF
           END DO

           !WRITE(*,*) "After Ion loop"
           !PAUSE
           !! Now use the composite density to add volume to the system
           Volume = Volume + (ionMass + AqMolecularMass(1)*InParticle%AqChems(1)) * DenomSum / NumerSum
           Mass   = Mass   + AqMolecularMass(1)*InParticle%AqChems(1)
           !write(*,*) 'Volume = ', Volume
           !write(*,*) 'ionMass = ',ionMass	
           !write(*,*) 'AqMolecularMass(1) = ',AqMolecularMass(1)
           !write(*,*) 'InParticle%AqChems(1) = ',InParticle%AqChems(1)
           !write(*,*) 'DenomSum = ',DenomSum
           !write(*,*) 'NumerSum = ',NumerSum

           !WRITE(*,*) "Ion solution Volume: ", Volume, DenomSum, NumerSum

           !If for some reason the above doesn't work, set inorganic solution density to 1.0
           IF(ISNaN(Volume)) THEN
              CALL WARN("Setting ion solution density to 1.0 g/cm3")
              Volume = Mass 
              !Volume = 1.0 
           END IF

           !This only includes water and inorganic ions (needed for optical property calculation)
           InorgSolnDensity = Mass/Volume
           !WRITE(*,*) "Inorganic Solution Density: ", InorgSolnDensity
        ENDIF


	!Now add in contributions from Hydrophilic Organics (assumed constant density)
	!! This assumes a fixed density for these aqueous phase chemicals
	!WRITE(*,*) "Before Aq Org loop Volume = ", Volume
	DO I = 1, HowManyAqOrgChems
		!write(*,*)
		Mass   = Mass   + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
		Volume = Volume + InParticle%AqOrgChems(I)*AqOrgMolecularMass(I)/AqOrgDensity(I) 
		!write(*,*) 'AqOrgPhaseChemicalNames(I) = ',AqOrgPhaseChemicalNames(I)
		!write(*,*) 'Volume =', Volume
		!write(*,*) 'Mass = ',Mass	
		!write(*,*) 'AqOrgMolecularMass(I) = ',AqOrgMolecularMass(I)
		!write(*,*) 'InParticle%AqOrgChems(I) = ',InParticle%AqOrgChems(I) 
		!write(*,*) 'AqOrgDensity(I)  = ',AqOrgDensity(I) 
		!write(*,*) 'DenomSum = ',DenomSum
		!write(*,*) 'NumerSum = ',NumerSum

	END DO
	
	!WRITE(*,*) "After Aq Org loop", Volume
! 28 September 2016 CMB: expanding this statement to let the user know why there is a problem
! 7 October 2016 CMB: just set them to 1 again...ugh...
!	IF(ISNaN(Volume)) STOP
!	if (IsNaN(Volume)) then
!		call error("In ParticleAttributes.h; volume is NaN after Aq Org loop, stop statement called")
!	endif
	if (isNaN(volume)) then
		call warn("Setting aq org solution density to 1.0 g/cm3")
		call flush(6)
		volume = mass
	endif

	!! This is purely the solution density (water, ions, non-solid electrolytes (e.g. NH3), aq organics)
	SolutionDensity =  Mass/Volume

	!!Calculate radius of the solution embryo
	InParticle%EmbryoRadius = (0.75 * Volume / pi)**0.33333333333
	!WRITE(*,*) "Embryo Radius: ", InParticle%EmbryoRadius
	
	!! Consider the contributions of Solid salts (except ions and water)
	InsolVolume = 0.
	InsolMass = 0.
	DO I = 2, HowManyAqChems
		SkipFlag = 0
		!Skip non-solid electrolytes (NH3, CO2, etc).
		DO C = 1, HowManyAqEqReactions
			! CMB
			if (dabs(aqequilibrialist(c,1) - I) .lt. 0.000001 .and. &
					dabs(AqEquilibriaList(c,6)) .lt. 1.0e-12) then
			!IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				SkipFlag = 1
			END IF
		END DO
			
		IF(SkipFlag .EQ. 0) THEN
			InsolMass   = InsolMass   + AqMolecularMass(I)*InParticle%AqChems(I)
			InsolVolume = InsolVolume + InParticle%AqChems(I)*AqMolecularMass(I)/SolidSaltDensity(I)

			Mass   = Mass   + AqMolecularMass(I)*InParticle%AqChems(I)
			Volume = Volume + InParticle%AqChems(I)*AqMolecularMass(I)/SolidSaltDensity(I)
		END IF
	END DO
	
	!WRITE(*,*) "After Solid Salts loop", Volume
        ! CMB: Revise to add more info to exit
		! update again: set to 1 if NaN
!	IF(ISNaN(Volume)) STOP
        if (isnan(volume)) then
            call warn("volume is NaN after solid salts loop")
		call flush(6)
        	volume = mass	! cmb add
		endif

	 ! cmb add
 !778 continue
 	 ! end cmb add
	 
	!Now add in controbutions from hydrophobic organics	
	DO I = 1, HowManyOrgChems
		!write(*,*) 'In loop'
		!write(*,*) 'InsolMass = ',InsolMass
		!write(*,*) 'OrgMolecularMass(I) = ',OrgMolecularMass(I)	

		IF(ISNaN(InParticle%OrgChems(I))) THEN
			InParticle%OrgChems(I) = 0.
		END IF

		!write(*,*) 'InParticle%OrgChems(I) = ',InParticle%OrgChems(I),OrgPhaseChemicalNames(I)
		!write(*,*) 'InsolVolume = ',InsolVolume
		!write(*,*) 'OrgDensity(I) = ',OrgDensity(I)
	
		InsolMass   = InsolMass   + OrgMolecularMass(I)*InParticle%OrgChems(I)
		InsolVolume = InsolVolume + InParticle%OrgChems(I)*OrgMolecularMass(I)/OrgDensity(I) 
	
		!WRITE(*,*) I, InsolVolume
		Mass   = Mass   + OrgMolecularMass(I)*InParticle%OrgChems(I)

		Volume = Volume + InParticle%OrgChems(I)*OrgMolecularMass(I)/OrgDensity(I) 
	END DO
		
!	WRITE(*,*) "After Hydrophobic Org loop", Volume
        ! CMB: Revise to add more info to exit
		! update: make NaN into 1 and don't exit'
!	IF(ISNaN(Volume)) STOP
        if (isnan(volume)) then
            !call error("In ParticleAttributes.h; volume is NaN after hydrophobic org loop")
		    call warn("volume is NaN after hydrophobic org loop")
        	volume = mass
		endif
	
	!!This is the average density of the non-aqeous part (solid salts, BC, and hydrophobic organics)
	InsolDensity = InsolMass/InsolVolume
	InParticle%InsolubleRadius = (0.75 * InsolVolume / pi)**0.33333333333

	!! This is the average density of the particle
	ParticleDensity = Mass / Volume
	!IF(Mass .LT. 0.0 ) WRITE(*,*) "Mass (<0), Volume = ", Mass, Volume

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! RESET EFFECTIVE RADIUS, too. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	InParticle%EffectiveRadius = (0.75 * Volume / pi)**0.33333333333


	InParticle%InorgSolnDensity = InorgSolnDensity 
	InParticle%SolutionDensity = SolutionDensity
	InParticle%InsolubleDensity = InsolDensity
	InParticle%ParticleDensity = ParticleDensity

	!WRITE(*,*) "End Particle Density", Volume, Mass
	


	RETURN
END FUNCTION ParticleDensity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Particle Mass !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleMass (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro
	USE Chemistry,		 ONLY : AqMolecularMass, HowManyAqChems,	&
					HowManyAqAnions, HowManyAqCations,	&
					AqCationMass, AqAnionMass, &
					HowManyOrgChems, OrgMolecularMass, &
					HowManyAqOrgChems, AqOrgMolecularMass 

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I

	ParticleMass = 0. 

	!! Consider the contributions of Aqueous Chemicals (except ions)
	!! This assumes a fixed density for these aqueous phase chemicals
	DO I = 1, HowManyAqChems
		ParticleMass   = ParticleMass + AqMolecularMass(I)*InParticle%AqChems(I)
	END DO

	DO I = 1, HowManyAqCations
		ParticleMass   = ParticleMass + AqCationMass(I)*InParticle%AqChems(HowManyAqChems+I)
	END DO

	DO I = 1, HowManyAqAnions
		ParticleMass   = ParticleMass + AqAnionMass(I)*InParticle%AqChems(HowManyAqChems+HowManyAqCations+I)
	END DO

	DO I = 1, HowManyOrgChems
		ParticleMass   = ParticleMass + OrgMolecularMass(I)*InParticle%OrgChems(I)
	END DO

	DO I = 1, HowManyAqOrgChems
		ParticleMass   = ParticleMass + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
	END DO

	RETURN
END FUNCTION ParticleMass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get a vector of weightings for the Patwardhan and Kumar	   !!
!! y_ij factors corrected for partial dissociation, incomplete !!
!! reaction sets, etc.										   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns 1). y_ij  2). Cation Index 3). Anion Index		   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Patwardhan and Kumar (AIChE J. 39, 711-714)				   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GetPatwardhandAndKumarWeightings (InParticle)

	USE Chemistry, ONLY : HowManyAqEqReactions,	&
			      HowManyAqChems,		&
			      HowManyAqCations,		&
			      HowManyAqAnions,		&
			      AqEquilibriaList,		&
			      AqCationCharge,		&
			      AqAnionCharge,            &
			      AqPhaseChemicalNames

	USE InfrastructuralCode, ONLY : WARN, IsNaN

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: GetPatwardhandAndKumarWeightings(HowManyAqEqReactions,3)
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: K, CatIndex, AnIndex
	REAL*8  :: II, JJ,														  &
			   CationicStrFr(HowManyAqCations),AnionicStrFr(HowManyAqAnions), &
			   CationChgFrac(HowManyAqCations),AnionChgFrac(HowManyAqAnions)

	!! Update the ionic strength
	!WRITE(*,*) "CRL: Before CalculateIonic Strength... "
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%IonicStr
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%AqChems
	!CALL CalculateIonicStrength(InParticle)
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%IonicStr
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%AqChems
	IF(IsNAN(InParticle%IonicStr)) CALL WARN("Ionic Strength Error in PKweightings!")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make Charge Fraction Vectors !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Calculate Cationic First:
	II = 0.
	DO K = 1, HowManyAqCations
		CationChgFrac(K) = CationMolality(K, InParticle) * AqCationCharge(K)
		II               = II + CationChgFrac(K)

	END DO
	DO K = 1, HowManyAqCations
		CationChgFrac(K) = CationChgFrac(K) /II
		IF(IsNAN(CationChgFrac(K))) CALL WARN("CationChgFrac error in PKweightings!")
		IF(IsNAN(CationChgFrac(K))) WRITE(*,*)  CationMolality(K, InParticle),AqCationCharge(K),II
	END DO

	!! And Anionic Next:
	II = 0.
	DO K = 1, HowManyAqAnions
		AnionChgFrac(K) = AnionMolality(K, InParticle) * AqAnionCharge(K)
		II              = II + AnionChgFrac(K)
	END DO
	DO K = 1, HowManyAqAnions
		AnionChgFrac(K) = AnionChgFrac(K) / II
		IF(IsNAN(AnionChgFrac(K))) CALL WARN("AnionChgFrac error in PKweightings!")
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make Ionic Strength Fraction Vectors !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (Inparticle%IonicStr .EQ. 0.) THEN
		DO K = 1, HowManyAqCations
			CationicStrFr(K) = 0.
		END DO
		DO K = 1, HowManyAqAnions
			AnionicStrFr(K)  = 0.
		END DO
	ELSE
		DO K = 1, HowManyAqCations
			CationicStrFr(K) = CationMolality(K, InParticle) * AqCationCharge(K) * AqCationCharge(K) / 2. / InParticle%IonicStr
			IF(IsNAN(CationicStrFr(K))) CALL WARN("CationicStrFr error in PKweightings!")

		END DO
		DO K = 1, HowManyAqAnions
			AnionicStrFr(K)  = AnionMolality(K, InParticle)  * AqAnionCharge(K) * AqAnionCharge(K) /2. / InParticle%IonicStr
			IF(IsNAN(AnionicStrFr(K))) CALL WARN("AnionicStrFr error in PKweightings!")
		END DO
	END IF

	
	!! This will be the sum of all the y_ij Patwardhan and Kumar mixing rule terms
	JJ = 0.

	!! Loop over all reactions
	DO K = 1, HowManyAqEqReactions
		

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! For partial dissociation, consider the reactions as if the partially !!
		!! dissociated form exchanges directly with the fully dissociated form, !!
		!! which is not how it would have been written.							!!
		IF(AqEquilibriaList(K,1) .LE. HowManyAqChems) THEN ! electrolyte uncharged
			CatIndex = AqEquilibriaList(K,2)
			AnIndex  = AqEquilibriaList(K,3)
		ELSE											   ! electrolyte charged

			!! Is the partially dissociated "electrolyte" a cation or an anion?
			IF (AqEquilibriaList(K,1) .LE. HowManyAqChems+HowManyAqCations) THEN  ! It's a Cation
				CatIndex = AqEquilibriaList(K,1)-HowManyAqChems
				AnIndex  = AqEquilibriaList(K,3)
			ELSE																  ! It's an Anion
				CatIndex = AqEquilibriaList(K,2)
				AnIndex  = AqEquilibriaList(K,1)-HowManyAqChems-HowManyAqCations
			END IF
		END IF
			
		!! Store the indices used to calculate this	
		GetPatwardhandAndKumarWeightings(K,2) = REAL(CatIndex)
		GetPatwardhandAndKumarWeightings(K,3) = REAL(AnIndex)

		
		!Exclude Levitocite reaction from this calculation
		! cmb
		if (dabs(AqEquilibriaList(K, 22)) .gt. 1e-40) then
		!IF(AqEquilibriaList(K,22) .NE. 0) THEN
			IF (AqPhaseChemicalNames(INT(AqEquilibriaList(K,1))) .EQ. "(NH4)3H(SO4)2") THEN
				GetPatwardhandAndKumarWeightings(K,1) = 0.
			END IF
		
		!Normal Calculation
		ELSE
			!! Calculate the Appropriate Weightings for the Patwardhan and Kumar (AIChE J. 39, 711-714)
			!! Mixed Activity Coefficient:
			!!
			!!   ln(Aw_mix) = Sum[(i,j), (Isf[i]*Cf[j]+Isf[j]*Cf[i])*ln(Aw[ij])]
			GetPatwardhandAndKumarWeightings(K,1) = AnionicStrFr(AnIndex)*CationChgFrac(CatIndex) + &
												CationicStrFr(CatIndex)*AnionChgFrac(AnIndex)
		
		END IF
		
		JJ = JJ + GetPatwardhandAndKumarWeightings(K,1)
	END DO

	!! Normalize the y_ij estimates to equal one
	DO K = 1, HowManyAqEqReactions
		GetPatwardhandAndKumarWeightings(K,1) = GetPatwardhandAndKumarWeightings(K,1) / JJ 
	END DO

	RETURN
	!write(*,*) 'GetPatwardhandAndKumarWeightings =',GetPatwardhandAndKumarWeightings 
END FUNCTION GetPatwardhandAndKumarWeightings


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Estimate Dissociated Electrolyte Concentrations !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For some routines, it is necessary to find the hypothetical
!! electrolyte concentrations that a set of dissociated ions 
!! would make up.  This function provides a vector of such electrolyte
!! abundances, in the form of a vector matching the AqChems vector
!! in the particle structure.
!!
!! The algorithm is taken from Jacobson's Fundamentals of Atmospheric
!! Modeling, 1999 page 494
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION HypotheticalElectrolyteConcentrations (InParticle)

	USE Chemistry, ONLY : HowManyAqChems, HowManyAqCations, &
                              HowManyAqAnions, HowManyAqEqReactions, &
                              AqEquilibriaList, & 
                              AqPhaseChemicalNames, AqAnionNames, &
                              AqCationNames

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterationsForHypLoop, &
                                    moles,WaterEquilibriumIndex

	! CMB: add warn
	USE InfrastructuralCode, ONLY : ERROR, REAL2STR, TRANSCRIPT, warn

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: HypotheticalElectrolyteConcentrations(HowManyAqChems)
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	REAL*8  :: TempChem(HowManyAqChems+HowManyAqCations+HowManyAqAnions), DeltaMoles, Sum, ContentSum, ExcludeRxnInit
	INTEGER :: I, J, K, C, ExcludeRxn 

	!! counters and settings
	K = 0
	ExcludeRxn = 0

	!! First set all electrolytes to zero, unless they are non-solid
	!write(*,*) '(HypotheticalElectrolyteConcentrations) Particle ID: ', inparticle%particleid
10	ContentSum=0
	DO I = 1, HowManyAqChems
		TempChem(I) = 0.
		!write(*,*) 'inparticle%aqchems(', i, ') = ', inparticle%aqchems(i)
		if (inparticle%aqchems(i) .lt. 0) then
			call warn("InParticle%AqChems("//trim(real2str(real(i)))//") < 0 ("//&
					trim(real2str(inparticle%aqchems(i)))//"in HypoElectConc; ID="//trim(real2str(&
					real(inparticle%particleid))))
		endif
		!Check for non-solid electrolyte (NH3, CO2, etc), and add in.
		DO C = 1, HowManyAqEqReactions
			if ((abs(AqEquilibriaList(C, 1) - I) .lt. 0.000001) .and.(&
					abs(AqEquilibriaList(c, 6)) .lt. 0.000001)) then
			!IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				!if (inparticle%aqchems(i) .lt. 0) then
				!	call warn("InParticle%AqChems("//trim(real2str(real(i)))//") < 0 ("//&
				!			trim(real2str(inparticle%aqchems(i)))//"in HypoElectConc")
				!endif
				!write(*,*) 'assigning to temp chem for C index using value = ', c, i, inparticle%aqchems(i)
				TempChem(I) = InParticle%AqChems(I)
			END IF
		END DO			
		 
		ContentSum = ContentSum + InParticle%AqChems(I)
	END DO

	!WRITE(*,*) "Check 1 ContentSum: ", ContentSum
	!call flush(6)

	!! Then take the Ion Concentrations as is.
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		TempChem(I) = InParticle%AqChems(I)
		ContentSum  = ContentSum + InParticle%AqChems(I)
	END DO

	!WRITE(*,*) "Check 2 ContentSum: ", ContentSum

	!WRITE(*,*) "Check 3"

	!! Loop over the reactions and adjust the equations
	DO J = 1, 2
	DO I = 1, HowManyAqEqReactions
		  !WRITE(*,*) I, ExcludeRxn, AqEquilibriaList(I,1), HowManyAqChems
		IF (I .EQ. ExcludeRxn) THEN
			CYCLE

		!! Some reactions may have ions serving as electrolytes
		!! We must identify these, and then push them the other direction
		ELSE IF (AqEquilibriaList(I,1) .GT. HowManyAqChems) THEN

			!WRITE(*,*) TempChem (AqEquilibriaList(I,1))
			!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), TempChem (HowManyAqChems+AqEquilibriaList(I,2)), TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,4)
			!WRITE(*,*) AqAnionNames(AqEquilibriaList(I,3)), TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)), TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,5)
		
			TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) =&
                            TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) + &
                            TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,4)
			TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) = &
			    TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) + &
			    TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,5)
			TempChem (INT(AqEquilibriaList(I,1))) = TempChem(INT(AqEquilibriaList(I,1))) - &
                            TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,5)

			!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), TempChem (HowManyAqChems+AqEquilibriaList(I,2))
			!WRITE(*,*) AqAnionNames(AqEquilibriaList(I,3)), TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3))
			!WRITE(*,*) TempChem (AqEquilibriaList(I,1))


		!! Other Reactions are only water dissociations
		ELSE IF (AqEquilibriaList(I,1) .LE. 0) THEN

			!! The adjustment parameter is limited by the concentration of each ion divided
			!! by the appropriate stoicheometric coefficient.
			DeltaMoles  = MIN(TempChem(HowManyAqChems+INT(AqEquilibriaList(I,2)))/AqEquilibriaList(I,4),&
                                TempChem(HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3)))/&
                                    AqEquilibriaList(I,5))

			!! Now adjust all of the concentrations appropriately
			IF (DeltaMoles .GT. 0) THEN
				!! Don't adjust back towards water at all, just get rid of it
				TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) = &
					TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) - &
                                        DeltaMoles*AqEquilibriaList(I,4)
				TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) = &
					TempChem (HowManyAqChems+HowManyAqCations + &
                                                INT(AqEquilibriaList(I,3))) - DeltaMoles * AqEquilibriaList(I,5)
			END IF

		!! And still others are canonical
		ELSE
			! CMB: more floating-point equality fixes
			!if (abs(AqEquilibriaList(I,22)) .le. 1.0e-40) then
			IF (AqEquilibriaList(I,22) .EQ. 0) THEN
			
				!! The adjustment parameter is limited by the concentration of each ion divided
				!! by the appropriate stoicheometric coefficient.
				DeltaMoles  = MIN(TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) /&
                                        AqEquilibriaList(I,4),	& 
                                        TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) /&
                                        AqEquilibriaList(I,5))
                                !WRITE(*,*) "Delta Moles: ", DeltaMoles
				!! Now adjust all of the concentrations appropriately
				IF (DeltaMoles .GT. 0) THEN
					TempChem (INT(AqEquilibriaList(I,1))) = TempChem(INT(AqEquilibriaList(I,1))) +&
                                            DeltaMoles
                                        TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) = &
                                            TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) - &
                                            DeltaMoles * AqEquilibriaList(I,4)
					TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) =&
                                            TempChem (HowManyAqChems+HowManyAqCations + &
                                                    INT(AqEquilibriaList(I,3))) - DeltaMoles * AqEquilibriaList(I,5)
				END IF
		

			ELSE
				!Levitocite case: Skip this
				!write(*,*) 'line 830 in particleattributes.h'
				CYCLE						
			END IF
		
		END IF
	END DO
	END DO

	!WRITE(*,*) "Check 4"

	!! Check to make sure have accounted for all of the ions
	SUM = 0.
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		!IF(IsNan(TempChem(I))) THEN
		!	WRITE(*,*) "Error in HypElec!"
		!	STOP
		!END IF
		SUM = SUM + TempChem(I)
	END DO

	!! This is an error
	!! I changed zero limit to 1.0e-3 from 1.0e-12 to eliminate some problems (Matt Alvarado, 5/10/2006)
	!! Note that this does not error so long as the amount of ions left is less than 1/1000 of the total aqueous moles!
	!! For dry particles, just runs!
	IF (SUM .GT. ContentSum*1.e-3) THEN
	!! It should not be caught in an infinite loop!

		!WRITE(*,*) "Sum, ContentSum: ", SUM, ContentSum*1.0e-3

		!! Show us what was wrong if the calculation doesn't work right.  (Probably the ContentSum*1.e-12 is not right below.)
		CALL TRANSCRIPT("")
		CALL TRANSCRIPT("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		CALL TRANSCRIPT("+++ There are problems in Hypothetical Electrolyte Concentrations +++")
		CALL TRANSCRIPT("")

		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions

		IF (TempChem(I) .gt. 0) then
		IF (I .LE. HowManyAqChems) CALL TRANSCRIPT(TRIM(AqPhaseChemicalNames(I))//&
                        TRIM(REAL2STR(TempChem(I)/moles)))
		IF (I .GT. HowManyAqChems .AND. I .LE. HowManyAqChems+HowManyAqCations) &
		    CALL TRANSCRIPT(TRIM(AqCationNames(I-HowManyAqChems))//TRIM(REAL2STR(TempChem(I)/moles)))
		IF (I .GT. HowManyAqChems+HowManyAqCations .AND. I .LE. HowManyAqChems+HowManyAqCations+HowManyAqAnions) &
		    CALL TRANSCRIPT(TRIM(AqAnionNames(I-HowManyAqChems-HowManyAqCations))//&
                            TRIM(REAL2STR(TempChem(I)/moles)))
		END IF

		END DO

		CALL DumpParticleContentsAtError(InParticle, "FailedHypotheticalElectrolyteConcentrationsParticle.txt",	&
						 InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., &
                                                 InDoWaterActivity=.FALSE.,InDensity=.FALSE.,InDoSurfaceTension=.FALSE., &
                                                 InIonicStrength=.TRUE.)
		!Removed this so that routines warns of charge imbalance 
		!in transcript, but does not stop
		! CMB: replace with a warn
		! never mind, it gives a seg fault?
		!CALL warn ("HypotheticalElectrolyteConcentrations() is stuck in what appears to be an infinite loop.")
		!CALL ERROR ("HypotheticalElectrolyteConcentrations() is stuck in what appears to be an infinite loop.")

	END IF

	!! Transfer the good stuff into a shorter vector for output
	DO I = 1, HowManyAqChems
		!WRITE(*,*) I, TempChem(I)
		HypotheticalElectrolyteConcentrations(I) = TempChem(I)
		!WRITE(*,*) 'HypotheticalElectrolyteConcentrations =',TempChem(I), I,AqPhaseChemicalNames
	END DO

	RETURN
END FUNCTION HypotheticalElectrolyteConcentrations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the TERMINAL FALL VELOCITY for particles !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION TerminalVelocity (InParticle)

	USE GridPointFields, ONLY : GetAirDensity,GetDynamicViscosityOfAir

	IMPLICIT NONE

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Internal Vairables
	REAL*8 :: Dum

	TerminalVelocity = ReynoldsNumber (InParticle) * 0.5             &
					   * GetDynamicViscosityOfAir()  &
					   / GetAirDensity()		 &
					   / InParticle%EffectiveRadius

	RETURN
END FUNCTION TerminalVelocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Recalculate the Radius of the Particle !!
SUBROUTINE RecalculateRadius (InParticle)

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, IFAQ, OrgPhaseChemicalNames
	!USE outputroutines, only : spillbeans
	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Type local variables
	REAL*8  :: II

	! CMB (AER, Inc): Floating-point inequality check
	if (dabs(inparticle%numberofparticles) .le. 1.0e-40) then
	!IF (InParticle%NumberOfParticles .EQ. 0. ) THEN
		IF (InParticle%Edges(1) .GT. 0.) THEN
			InParticle%EmbryoRadius    = InParticle%Edges(1)
			InParticle%InsolubleRadius = 0.
			InParticle%EffectiveRadius = InParticle%Edges(1)
		ELSE
			InParticle%EmbryoRadius    = 0.01*InParticle%Edges(2)
			InParticle%InsolubleRadius = 0.
			InParticle%EffectiveRadius = 0.01*InParticle%Edges(2)
		END IF

		InParticle%ParticleDensity = 1.
		InParticle%SolutionDensity = 1.
		InParticle%InsolubleDensity = 1.
		InParticle%InorgSolnDensity = 1.
		RETURN
	END IF

	!! Calling the Density Routine Resets the Radius, so just do that
	!Write(*,*) "In RecalculateRadius InParticle%AqChems(1) ", InParticle%AqChems(1)
	!Write(*,*) "In RecalculateRadius InParticle%OrgChems(14) ", InParticle%OrgChems(14),OrgPhaseChemicalNames(14)
	!Write(*,*) "Just before calling ParticleDensity spilling beans"
	!call spillbeans()
	II = ParticleDensity (InParticle)
	!Write(*,*) "End ParticleDensity() "
        !! Reset the Surface Tension as well
        II = SurfaceTension (InParticle)	
	!Write(*,*) "End SurfaceTension() "
	!write(*,*) 'Exiting RecalculateRadius()'
	RETURN
END SUBROUTINE RecalculateRadius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The correction to pressure from curvature, used for water activity and in condensation. !!
!! It is also known as the "Kelvin Effect"						   !!
!!											   !!
!! The optional argument PerturbWater is a multiple by which to change the water content   !!
!! so we can get a sense of the derivative of the water residual.			   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CurvatureCorrection (InParticle) 

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : Rstar, micron, ThermoBulkMode, IgnoreCorrectionstoDiffusivity
	USE GridPointFields, ONLY : SurfaceTensionofWater

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER  :: InParticle
	REAL*8 :: TotalMoles
	INTEGER :: I

	!!Type Internal variables
	REAL*8 :: TotalDens
	!WRITE(*,*) "Curvature Started", SurfaceTension(InParticle), ParticleDensity(InParticle)
	
	!Ignore this for bulk mode, no particles, or if ordered to.
	IF (ThermoBulkMode .OR. InParticle%NumberofParticles .LE. 0.0D0 .OR. IgnoreCorrectionstoDiffusivity) THEN
	   CurvatureCorrection = 1.
	   !write(*,*) 'since InParticle%NumberofParticles should be le zero, CurvatureCorrection = 1',CurvatureCorrection, InParticle%NumberofParticles 	
	ELSE
	   !! Do this in stages because ParticleDensity forces a radius + solution density update...
	   !WRITE(*,*) "Before Dens"
	   TotalDens  = InParticle%SolutionDensity
	   !WRITE(*,*) "After Dens" 
		
	   !WRITE(*,*) "Curvature Correction: ", InParticle%NumberofParticles
	   !write(*,*) 'InParticle%NumberofParticles should be gt zero', InParticle%NumberofParticles
	   CurvatureCorrection = 2.*AqMolecularMass(1)/(InParticle%Temperature*Rstar*InParticle%SolutionDensity)
	   !WRITE(*,*) CurvatureCorrection, InParticle%EmbryoRadius, SurfaceTension(InParticle)
		
	   CurvatureCorrection = EXP(CurvatureCorrection*InParticle%SurfaceTension/InParticle%EmbryoRadius)
	   !WRITE(*,*) "Curvature Correction: ", CurvatureCorrection
	END IF

	!WRITE(*,*) "Curvature Okay", CurvatureCorrection 
	RETURN
	
END FUNCTION CurvatureCorrection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Surface Tension depends on Temperature and Electrolytic Properties  !!
!! of the Aerosol's constituent species.			       !!
!!								       !!
!! Follow Li and Lu (2001): Chem. Eng. Sci., 56: 2879-2888, who use    !!
!! a Gibbs dividing surface, define surface excess, etc.  It is a 2-   !!
!! parameter-per-electrolyte system and the data is taken in when the  !!
!! electrolytes are defined.					       !!
!!								       !!
!! Many other systems seem to use the sig = sig_w + Beta I formulation !!
!! from Pruppacher and Klett (1997), but that is significantly less    !!
!! effective than this technique.				       !!
!!                                                                     !!
!! Units are dyn / cm                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION SurfaceTension (InParticle)

	USE GridPointFields, ONLY : SurfaceTensionOfWater, &
                                    GetRelativeHumidity

	USE Chemistry,       ONLY : AqEquilibriaList,      &
                                    HowManyAqEqReactions
	
        USE ModelParameters, ONLY : UseDonnanSurfaceTension, &
                                    Rstar,                   &
                                    HowManyAqChems,          &
                                    HowManyAqAnions,         &
                                    HowManyAqCations   
	
        USE InfrastructuralCode, ONLY : Error, REAL2STR
	IMPLICIT NONE

	!! External Variables
	TYPE (Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I
	REAL*8  :: II, Sum1, DenomSum, TotalIonMolality, Test

	! CMB: equality check
	if (inparticle%numberofparticles .le. 1.0e-40) return
	!IF (InParticle%NumberOfParticles .EQ. 0.0) RETURN

	!! CALL FindElectrolyteEquilibrium (InParticle)
	CALL CalculateIonicStrength(InParticle)
        CALL KusikMeissner(InParticle)

	IF(UseDonnanSurfaceTension) THEN !If true, Donnan; If False, Jacobson formulas
		!!Donnan's Parameterization
		!! 1. + SUM (K_i * Activity_i)
		DenomSum = 1.
		DO I = 1, HowManyAqEqReactions
			DenomSum = DenomSum + MeanActivity(I, InParticle) * AqEquilibriaList(I,15)
		END DO

!		IF (Denomsum .GT. 1.0e40) THEN
!			WRITE(*,*) DenomSum, InParticle%numberofParticles
!			CALL DumpParticleContentsAtError(InParticle, "BadActivitiesInSurfaceTensionParticle.txt", &
!                                                 InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., &
!                                                 InDoWaterActivity=.TRUE.,InDensity=.TRUE.,InDoSurfaceTension=.FALSE., &
!                                                 InIonicStrength=.TRUE.)
!			CALL ERROR ("In SurfaceTension(), the mean activity coefficients appear to have grown without bounds.  "// &
!					"This is likely because the water is leaching out of the system towards zero.  Perhaps the "// &
!					"relative humidity for the particle is too small?  I think it is "//trim(real2str(100*         &
!					GetRelativeHumidity()))//  &
!					"%.  Another possibility is that you are in the water equilibration routine and it is dumping "// &
!					"water from the particle to try to reduce the water activity, and it got too dry.")
!		END IF

		!! Loop over each reaction and consider the appropriate equation
		Sum1 = 0.
		DO I = 1, HowManyAqEqReactions

			!! Ignore dissociation of 
			!IF ((AqEquilibriaList(I,1) .GT. HowManyAqChems) .OR. (AqEquilibriaList(I,1) .LT. 0)) CYCLE

			!! AqIonGrid(I,14) is Surface Excess and AqIonGrid(I,15) is Equilibrium Adsorption Coefficient
			!WRITE(*,*) "Mean Act.: ", MeanActivity(I, InParticle)
			!WRITE(*,*) "Ka: ", AqEquilibriaList(I,15)
			!WRITE(*,*) "DenomSum: ", DenomSum
			Sum1 = Sum1 + AqEquilibriaList(I,14) * LOG (1. - AqEquilibriaList(I,15)	&  
				* MeanActivity(I, InParticle) / DenomSum) 

		END DO

		!! Surface Tesion of the Solution is the Surface Tension of Water + R * T * Sum1
		SurfaceTension = SurfaceTensionOfWater() &
	                 + Rstar*InParticle%Temperature*Sum1
	
		!Now, add in effect of aqueous organics using formula in
		! Jacobson, "Fund. of Atm. Modeling", 2nd ed., p. 534
                !Based on Fachini et al., Nature, 401, 257-259, 1999.
 		SurfaceTension = SurfaceTension - 0.0187*InParticle%Temperature*LOG(1+628.14*MolalityAqCarbon(InParticle))
	
	ELSE !Use formulas from Jacobson, "Fund. of Atm. Modeling", 2nd ed., p. 534		
		SurfaceTension = SurfaceTensionOfWater() 
		!Ions
		TotalIonMolality = 0.
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			TotalIonMolality = TotalIonMolality + Molality(I, InParticle)
		END DO
                !Based on Pruppacher and Klett, 1997.
		SurfaceTension = SurfaceTension + 1.7*TotalIonMolality
		!Organics
                !Based on Fachini et al., Nature, 401, 257-259, 1999.
		SurfaceTension = SurfaceTension - 0.0187*InParticle%Temperature*LOG(1+628.14*MolalityAqCarbon(InParticle))

	END IF

	!!Force a minimum Surface tension
	!!Use the surface tension of Benzene at 298 K 
	!!(Also used in OrgCurvatureCorrection as the surface tension)		
	!!See Seinfeld and Pandis (1998), Table 10.2, p. 558
	IF(SurfaceTension .LT. 28.21) SurfaceTension = 28.21

	!! Units are dyn / cm
	InParticle%SurfaceTension = SurfaceTension

	RETURN
END FUNCTION SurfaceTension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the PARTICLE REYNOLDS NUMBER !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Uses a three-regime (slip, continuum, and transition) !!
!! from Beard's parameterization. (see Jacobson, 1997)   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tested against Fig. 16.3 in Jacobson, 1999 and it     !!
!! matches well.					 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ReynoldsNumber (InParticle)

	USE GridPointFields, ONLY : GetAirDensity,GetDynamicViscosityOfAir
	USE ModelParameters, ONLY : g, grams, cm
	IMPLICIT NONE

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Type local variables
	REAL*8 :: Dv, Cs, AirDensity, X, Y

	Dv         = GetDynamicViscosityOfAir()
	Cs         = CunninghamSlipCorrectionFactor(InParticle)
	AirDensity = GetAirDensity()

	!! Estimate the Reynolds Number and then Calculate the Real Value, 
	!! this is the estimated terminal velocity
	ReynoldsNumber = 2*InParticle%EffectiveRadius*InParticle%EffectiveRadius*(InParticle%ParticleDensity-AirDensity)*g*Cs/(9*DV)

	!! The Estimated Reynolds Number
	ReynoldsNumber = 2.*InParticle%EffectiveRadius*AirDensity*ReynoldsNumber/DV 

	!! The temporary is the Reynolds Number unless Re > 0.01
	IF (ReynoldsNumber > 0.01 .and. ReynoldsNumber < 300) THEN

		X = dLog((InParticle%ParticleDensity-AirDensity)*AirDensity*InParticle%EffectiveRadius**3.*32.*g/(3.*DV*DV))

		!! The Values here taken from Prupp and Klett, p.417 and originate from Beard's paper (1976)
		Y = -.318657d1+(0.992696d0)*X-(0.153193d-2)*X**2.-(0.987059d-3)*X**3.	&
			-(0.578878d-3)*X**4.+(0.855176d-4)*X**5.-(0.327815d-5)*X**6.

		ReynoldsNumber = Cs*dExp(Y)

	ENDIF

	IF (ReynoldsNumber .GE. 300.) THEN

		!! Use X to signify Np**1/6
		X = (InParticle%SurfaceTension**3.*AirDensity*AirDensity/(DV**4. * (InParticle%ParticleDensity-AirDensity) * g))**(1./6.)

		Y = dLog(16./3.*X*InParticle%EffectiveRadius*InParticle%EffectiveRadius*	&
				 (InParticle%ParticleDensity-AirDensity)*g/InParticle%SurfaceTension)

		ReynoldsNumber = Cs*X*dExp(-5.00015+5.23778*Y-2.04914*Y**2.+0.475294*Y**3.-0.0542819*Y**4.+0.00238449*Y**5.)

	END IF

	RETURN
END FUNCTION ReynoldsNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The correction to pressure from curvature of the organic phase.			   !!
!! It is also known as the "Kelvin Effect"						   !!
!! This uses approximate values for molecular mass and surface tension                     !!
!! of the organic phase. Since the curvature correction is mainly for                      !!
!! water above RH = 98%, this should be of minor importance.				   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION OrgCurvatureCorrection (InParticle) 

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : Rstar, micron, ThermoBulkMode, IgnoreCorrectionstoDiffusivity
	USE GridPointFields, ONLY : SurfaceTensionofWater

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER  :: InParticle
	REAL*8 :: TotalMoles
	INTEGER :: I

	!!Type Internal variables
	REAL*8 :: TotalDens, SurfTens, Radius
	
	!Ignore this for bulk mode or no insoluble species
	IF (ThermoBulkMode .OR. IgnoreCorrectionstoDiffusivity .OR. InParticle%NumberofParticles .EQ. 0.0) THEN
		OrgCurvatureCorrection = 1.

	ELSE
		TotalDens = InParticle%InsolubleDensity	
		!!Assume a molecular mass of 200 g/mol
		OrgCurvatureCorrection = 2.*200.0/(InParticle%Temperature*Rstar*InParticle%InsolubleDensity)

		!Keep radius above 1 nm for calculation
		IF(InParticle%InsolubleRadius .LE. 0.001*micron) THEN
			Radius = 0.001*micron
		ELSE
			Radius = InParticle%InsolubleRadius
		END IF
		
		!!Assume surface Tension of Benzene, Seinfeld and Pandis, p. 522
		SurfTens = 28.21
		OrgCurvatureCorrection = EXP(OrgCurvatureCorrection*SurfTens/Radius)

	END IF

	!WRITE(*,*) "OrgCurvature Okay"
	RETURN
	
END FUNCTION OrgCurvatureCorrection

!Molality of carbon (mol/kg H2O) in aqueous phase for surface tension calculation (Matt Alvarado, 6/2006)
REAL*8 FUNCTION MolalityAqCarbon(InParticle)

	USE Chemistry,				ONLY :	AqOrgNumCarbon
	USE ModelParameters
	
	IMPLICIT NONE

	REAL*8  :: Total
	INTEGER :: I
	TYPE(PARTICLE),POINTER :: InParticle

	
	Total = 0.0
	
	DO I = 1, HowManyAqOrgChems
		Total = Total + InParticle%AqOrgChems(I)*AqOrgNumCarbon(I)
	END DO
	
	MolalityAqCarbon = Total / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

	RETURN

END FUNCTION MolalityAqCarbon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Mass of aqueous phase (g/particle) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION AqueousSolutionMass (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro
	USE Chemistry,			 ONLY : AqMolecularMass, HowManyAqChems,	&
									HowManyAqAnions, HowManyAqCations,	&
									AqCationMass, AqAnionMass, &
									HowManyAqOrgChems, AqOrgMolecularMass, &
									HowManyAqEqReactions, AqEquilibriaList

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I, C

	AqueousSolutionMass = 0. 

	!Water
	AqueousSolutionMass   = AqueousSolutionMass + AqMolecularMass(1)*InParticle%AqChems(1)

	!Non-solid Electrolytes
	DO I = 2, HowManyAqChems
		!Check for non-solid electrolyte (NH3, CO2, etc), and add in.
		DO C = 1, HowManyAqEqReactions
			IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				AqueousSolutionMass = AqueousSolutionMass + AqMolecularMass(I)*InParticle%AqChems(I)
			END IF
		END DO					 
	END DO

	DO I = 1, HowManyAqCations
		AqueousSolutionMass   = AqueousSolutionMass + AqCationMass(I)*InParticle%AqChems(HowManyAqChems+I)
	END DO

	DO I = 1, HowManyAqAnions
		AqueousSolutionMass   = AqueousSolutionMass + AqAnionMass(I)*InParticle%AqChems(HowManyAqChems+HowManyAqCations+I)
	END DO

	DO I = 1, HowManyAqOrgChems
		AqueousSolutionMass   = AqueousSolutionMass + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
	END DO
	!write(*,*) 'AqueousSolutionMass   = ',AqueousSolutionMass 
	RETURN
END FUNCTION AqueousSolutionMass

!! Patch 1.3: homogeneous-sphere Mie efficiencies via MIEV0 (Wiscombe) instead of DMiLay.
!! DMiLay with RindCo = RindSh is a homogeneous sphere of radius Rshell; MIEV0 computes that case
!! directly (no core/shell layering).  Outputs map as QEXT, QSCA, GQSC (= g*Qsca), and
!! QBS = |SBACK|**2 / (Pi*x**2) (DMiLay's definition of its backscatter output).
!! MIEV0.f is single-precision source; like the rest of ASP it relies on default REAL = 8 bytes.
SUBROUTINE MieHomogeneousSphere (Radius, WVNO, Ind, QExt, QSca, QBack, GQSca)
	IMPLICIT NONE
	REAL*8  :: Radius, WVNO, QExt, QSca, QBack, GQSca
	COMPLEX :: Ind
	REAL*8, PARAMETER :: PIc = 3.14159265358979324D0
	REAL*8  :: XX, MIMCUT, SPIKE, XMU(1), PMOM(0:1,1)
	COMPLEX :: SFORW, SBACK, S1(1), S2(1), TFORW(2), TBACK(2)
	LOGICAL :: PRNT(2), PERFCT, ANYANG
	INTEGER :: NUMANG, NMOM, IPOLZN, MOMDIM
	EXTERNAL MIEV0

	XX     = Radius*WVNO
	MIMCUT = 0.0D0
	PRNT(1) = .FALSE. ; PRNT(2) = .FALSE.
	XMU(1) = 0.0D0 ; PMOM = 0.0D0
	PERFCT = .FALSE. ; ANYANG = .FALSE. ; NUMANG = 0 ; NMOM = 0 ; IPOLZN = 0 ; MOMDIM = 1
	!! (all MIEV0 inputs are passed as variables: its first-call self-test overwrites them temporarily)
	CALL MIEV0 (XX, Ind, PERFCT, MIMCUT, ANYANG, NUMANG, XMU, NMOM, IPOLZN, MOMDIM, PRNT, &
	            QExt, QSca, GQSca, PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK, SPIKE)
	QBack = (REAL(SBACK)**2 + AIMAG(SBACK)**2) / (PIc*XX*XX)
	RETURN
END SUBROUTINE MieHomogeneousSphere

!! This subroutine calculates the volume-average refractive indices for
!! the reflective shell and then calls DMILAY for the efficienies and 
!! assymetry parameter
!! Volume concentrations (cm3 per cm3 of air) used for the optical properties: water, solid
!! salts and ions (AqVolConc), insoluble organics (OrgVolConc; index 1 is the BC core) and
!! hydrophilic organics (AqOrgVolConc), the aqueous solution (water and ions) and the shell
!! (the particle volume less the BC core).  Shared by ShellRefIndAndRad and
!! ShellRefIndAndRadFASTTUV.
SUBROUTINE OpticalVolumeConcentrations (InParticle, AqVolConc, OrgVolConc, AqOrgVolConc, &
                                        SolutionVolConc, ShellVolumeConc)
	USE ModelParameters
	USE Chemistry, ONLY : AqMolecularMass, SolidSaltDensity, &
				OrgMolecularMass, OrgDensity, &
				AqOrgMolecularMass, AqOrgDensity
	IMPLICIT NONE

	TYPE(Particle),POINTER :: InParticle
	REAL*8 :: AqVolConc(:), OrgVolConc(:), AqOrgVolConc(:), SolutionVolConc, ShellVolumeConc
	INTEGER :: I
	REAL*8	:: SingleParticleVolume, VolumeConc

	SingleParticleVolume = (4.0*Pi/3.0)*(InParticle%effectiveradius)**3 !cm3
	VolumeConc = SingleParticleVolume*InParticle%numberofparticles !cm3/cm3

	!Water
	AqVolConc(1)  = InParticle%AqChems(1) * &
                        InParticle%NumberOfParticles*AqMolecularMass(1)/ &
                        InParticle%SolutionDensity
	SolutionVolConc = AqVolConc(1)

	!Solid Salts
	DO I = 2, HowManyAqChems
		AqVolConc(I)  = InParticle%AqChems(I) * &
                                InParticle%NumberOfParticles*AqMolecularMass(I)/ &
                                SolidSaltDensity(I)
	END DO

	!Inorganic Ions
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		AqVolConc(I)  = InParticle%AqChems(I) * &
                                InParticle%NumberOfParticles*AqMolecularMass(I)/ &
                                InParticle%SolutionDensity
		SolutionVolConc = SolutionVolConc + AqVolConc(I)
	END DO

	!Insoluble organics
	DO I = 1, HowManyOrgChems
		OrgVolConc(I)  = InParticle%OrgChems(I) * &
                                 InParticle%NumberOfParticles*OrgMolecularMass(I)/ &
                                 OrgDensity(I)
	END DO

	!Hydrophilic organics
	DO I = 1, HowManyAqOrgChems
		AqOrgVolConc(I)  = InParticle%AqOrgChems(I) * &
                                   InParticle%NumberOfParticles*AqOrgMolecularMass(I) &
                                   /AqOrgDensity(I)
	END DO

	!Shell volume conc (Total - BC)
	ShellVolumeConc = VolumeConc - OrgVolConc(1)
END SUBROUTINE OpticalVolumeConcentrations

!! Shell and core refractive indices at one wavelength (nm), interpolated in the OPAC tables
!! (LagrangianAerosols.f90) between entries w_ind_low and w_ind_low+1 with weights a_w, b_w.
!! The core is BC (OrgChems(1)), the shell everything else.
!!  - No shell (all BC, or a dry particle): HasShell = .FALSE.; shell index = core index; the
!!    effective radius is set to the BC radius and the core radius to 0.99 of it.
!!  - Otherwise the shell index (also stored in ShellReal/ImagRefracInd(ib)) is the volume-
!!    average dielectric constant (SumEr, SumEi; Jacobson, Fund. of Atm. Modeling, 2nd ed.,
!!    p. 310) of the solution (Tang et al., JGR 102, 23269 and Tang, JGR 102, 1883, 1997),
!!    solid salts (by RefIndFlag) and organics.  Without BC the core index is the shell index
!!    and the core radius 0.01 of the effective radius; with BC the core radius is that of the
!!    BC volume, at most 0.9999 of the effective radius.
!! The organic terms are added as one term per class when OrgVolConcSum and AqOrgVolConcSum
!! are given (ShellRefIndAndRad), species by species otherwise (ShellRefIndAndRadFASTTUV);
!! the two orders round differently, so each routine keeps its own.
SUBROUTINE ShellAndCoreIndices (InParticle, ib, Wavelength, w_ind_low, AqVolConc, OrgVolConc, &
                                AqOrgVolConc, SolutionVolConc, ShellVolumeConc, a_w, b_w, &
                                SumEr, SumEi, SolarShellInd, SolarCoreInd, HasShell, &
                                OrgVolConcSum, AqOrgVolConcSum)
	USE ModelParameters
	USE Chemistry, ONLY : AqMolecularMass, &
                                CationIonicRefraction, &
                                AnionIonicRefraction,  &
                                RefIndFlag
	IMPLICIT NONE

	TYPE(Particle),POINTER :: InParticle
	INTEGER :: ib, w_ind_low
	REAL*8  :: Wavelength, AqVolConc(:), OrgVolConc(:), AqOrgVolConc(:), SolutionVolConc, &
	           ShellVolumeConc, a_w, b_w, SumEr, SumEi
	COMPLEX :: SolarShellInd, SolarCoreInd
	LOGICAL :: HasShell
	REAL*8, OPTIONAL :: OrgVolConcSum, AqOrgVolConcSum

	INTEGER :: I, CatIdx
	REAL*8  :: n_real, n_imag, V_molal, R_molal, Er, Ei, Sum, SolnMolWeight, MolalRefraction, &
	           MolVol, xxx, Factor, SolutionRefracIndex(2)
	REAL*8  :: MolFrac(1+HowManyAqCations+HowManyAqAnions)
	REAL*8  :: CatN_real(0:3), CatN_imag(0:3), CatEr(0:3), CatEi(0:3)

	!Calculate weightings of wavelength-dependent refractive index points
	a_w = (WAVE(w_ind_low+1)-Wavelength) &
	       /(WAVE(w_ind_low+1)-WAVE(w_ind_low))
	b_w = 1-a_w

	!Calculate Core Refractive Index
	n_real = a_w*SOOT_REAL(w_ind_low)+b_w*SOOT_REAL(w_ind_low+1)
	n_imag = -1.0*(a_w*SOOT_IMAG(w_ind_low)+b_w*SOOT_IMAG(w_ind_low+1))
	SolarCoreInd = CMPLX(n_real, -1.0*n_imag)

	!All core (no shell)
	HasShell = .NOT. (ShellVolumeConc .LE. 0 .OR. InParticle%Dry)
	IF (.NOT. HasShell) THEN
		!Shell and core refractive indices equal.
		SolarShellInd = SolarCoreInd
		!Calculate Core radius
		InParticle%AbsCoreRad = OrgVolConc(1)/InParticle%NumberofParticles
		InParticle%AbsCoreRad = (0.75*InParticle%AbsCoreRad/Pi)**(1.0/3.0) !In cm
		!Set Shell radius to above, and core a little less than shell
		InParticle%EffectiveRadius = InParticle%AbsCoreRad
		InParticle%AbsCoreRad = 0.99*InParticle%AbsCoreRad
		RETURN
	END IF

	!Shell refractive index: volume average dielectric constant mixing rule
	SumEr = 0.
	SumEi = 0.

	!Water refractive index and molal Refraction
	n_real = a_w*WATER_REAL(w_ind_low)+b_w*WATER_REAL(w_ind_low+1)
	n_imag = -1.0*(a_w*WATER_IMAG(w_ind_low)+b_w*WATER_IMAG(w_ind_low+1))
	V_molal = 18.0674 !Tang, 1997: density = 0.9971 g/cm3, MW = 18.015 g/mol, V = MW/dens
	R_molal = V_molal*(n_real**2 - 1.0)/(n_real**2 + 2)

	!Average refractive index for the ionic solution (Tang et al., 1997)
	!Step 1: Calculate mole fractions in solution
	Sum = InParticle%AqChems(1)
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		Sum = Sum + InParticle%AqChems(I)
	END DO
	MolFrac(1) = InParticle%AqChems(1)/Sum
	DO I = 2, 1+HowManyAqCations+HowManyAqAnions
		MolFrac(I) = InParticle%AqChems(HowManyAqChems+I-1)/Sum
	END DO

	!Step 2: Calculate average molecular weight
	SolnMolWeight = MolFrac(1)*AqMolecularMass(1)
	DO I = 2, 1+HowManyAqCations+HowManyAqAnions
		SolnMolWeight = SolnMolWeight + MolFrac(I)*AqMolecularMass(HowManyAqChems+I-1)
	END DO

	!Step 3: Calculate average molal refraction
	MolalRefraction = MolFrac(1)*R_molal !See Tang et. al, 1997
	DO I = 2, 1+HowManyAqCations
		MolalRefraction = MolalRefraction + MolFrac(I)*CationIonicRefraction(I-1)
	END DO
	DO I = 2+HowManyAqCations, 1+HowManyAqCations+HowManyAqAnions
		MolalRefraction = MolalRefraction + MolFrac(I)*AnionIonicRefraction(I-HowManyAqCations-1)
	END DO

	!Step 4: Calculate solution molal volume
	MolVol = SolnMolWeight/InParticle%InorgSolnDensity

	!Step 5: Calculate Refractive Index
	xxx = MolalRefraction/MolVol
	SolutionRefracIndex(1) = SQRT((2.*xxx+1.)/(1.-xxx))
	SolutionRefracIndex(2) = n_imag !Imaginary index for water

	!Step 6: Calculate Er and Ei
	Er = SolutionRefracIndex(1)**2 - SolutionRefracIndex(2)**2
	Ei = 2*SolutionRefracIndex(1)*SolutionRefracIndex(2)
	SumEr = SumEr + (SolutionVolConc/ShellVolumeConc)*Er
	SumEi = SumEi + (SolutionVolConc/ShellVolumeConc)*Ei

	!Solid salts: one refractive index per RefIndFlag category (0 to 3).  Category 0 mixes
	!WATER at w_ind_low with SULF at w_ind_low+1, as in the original code.
	CatN_real(0) = a_w*WATER_REAL(w_ind_low)+b_w*SULF_REAL(w_ind_low+1)
	CatN_imag(0) = -1.0*(a_w*WATER_IMAG(w_ind_low)+b_w*SULF_IMAG(w_ind_low+1))
	CatN_real(1) = a_w*SULF_REAL(w_ind_low)+b_w*SULF_REAL(w_ind_low+1)
	CatN_imag(1) = -1.0*(a_w*SULF_IMAG(w_ind_low)+b_w*SULF_IMAG(w_ind_low+1))
	CatN_real(2) = a_w*SEASALT_REAL(w_ind_low)+b_w*SEASALT_REAL(w_ind_low+1)
	CatN_imag(2) = -1.0*(a_w*SEASALT_IMAG(w_ind_low)+b_w*SEASALT_IMAG(w_ind_low+1))
	CatN_real(3) = a_w*DUST_REAL(w_ind_low)+b_w*DUST_REAL(w_ind_low+1)
	CatN_imag(3) = -1.0*(a_w*DUST_IMAG(w_ind_low)+b_w*DUST_IMAG(w_ind_low+1))
	DO CatIdx = 0, 3
		CatEr(CatIdx) = CatN_real(CatIdx)**2 - CatN_imag(CatIdx)**2
		CatEi(CatIdx) = 2*CatN_real(CatIdx)*CatN_imag(CatIdx)
	END DO
	DO I = 2, HowManyAqChems
		SumEr = SumEr + (AqVolConc(I)/ShellVolumeConc)*CatEr(RefIndFlag(I))
		SumEi = SumEi + (AqVolConc(I)/ShellVolumeConc)*CatEi(RefIndFlag(I))
	END DO

	!Insoluble organics (skipping the BC core) and hydrophilic organics, both as WASO
	n_real = a_w*WASO_REAL(w_ind_low)+b_w*WASO_REAL(w_ind_low+1)
	n_imag = -1.0*(a_w*WASO_IMAG(w_ind_low)+b_w*WASO_IMAG(w_ind_low+1))
	Er = n_real**2 - n_imag**2
	Ei = 2*n_real*n_imag
	IF (PRESENT(OrgVolConcSum)) THEN
		SumEr = SumEr + (OrgVolConcSum/ShellVolumeConc)*Er
		SumEi = SumEi + (OrgVolConcSum/ShellVolumeConc)*Ei
		SumEr = SumEr + (AqOrgVolConcSum/ShellVolumeConc)*Er
		SumEi = SumEi + (AqOrgVolConcSum/ShellVolumeConc)*Ei
	ELSE
		DO I = 2, HowManyOrgChems
			SumEr = SumEr + (OrgVolConc(I)/ShellVolumeConc)*Er
			SumEi = SumEi + (OrgVolConc(I)/ShellVolumeConc)*Ei
		END DO
		DO I = 1, HowManyAqOrgChems
			SumEr = SumEr + (AqOrgVolConc(I)/ShellVolumeConc)*Er
			SumEi = SumEi + (AqOrgVolConc(I)/ShellVolumeConc)*Ei
		END DO
	END IF

	!Shell refractive index for the core-in-shell mixing rule
	Factor = SQRT(SumEr**2 + SumEi**2)
	InParticle%ShellRealRefracInd(ib) = SQRT((Factor + SumEr)/2.0)
	InParticle%ShellImagRefracInd(ib) = SQRT((Factor - SumEr)/2.0)
	SolarShellInd = CMPLX(InParticle%ShellRealRefracInd(ib),-1.0*InParticle%ShellImagRefracInd(ib))

	IF (OrgVolConc(1) .LE. 0.0) THEN !No core
		!Set core and shell refractive indices equal
		SolarCoreInd = SolarShellInd
		!Set small core radius
		InParticle%AbsCoreRad = 0.01*InParticle%EffectiveRadius
	ELSE !Both Shell and Core
		!Calculate Core radius
		InParticle%AbsCoreRad = OrgVolConc(1)/InParticle%NumberofParticles
		InParticle%AbsCoreRad = (0.75*InParticle%AbsCoreRad/Pi)**(1.0/3.0) !In cm
		!This shouldn't happen, but be safe
		IF(InParticle%AbsCoreRad .GT. 0.9999*InParticle%EffectiveRadius) &
			InParticle%AbsCoreRad = 0.9999*InParticle%EffectiveRadius
	ENDIF
END SUBROUTINE ShellAndCoreIndices

SUBROUTINE ShellRefIndAndRad(InParticle)
! 08/30/2010 Changed for the Photolyis Wavebands MJA
! 02/27/2012 Changed to wavelength-dependent refractive indices
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added Maxwell Garnett mixing rule
	USE ModelParameters
	USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, ERROR
	
	IMPLICIT NONE
	
	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!!External Subroutines 
        EXTERNAL  DMiLay

	!! Internal Variables
	INTEGER, PARAMETER  :: MaxAng = 100
 	INTEGER	:: I, NumAng, ib, w_ind_low
	REAL*8	:: Factor, SumEr, SumEi, Wavelength, ShellVolumeConc, SolutionVolConc
        REAL*8  :: ErBC, EiBC, ErVA, EiVA, bcvolrat, ErMG, EiMG
	LOGICAL :: HasShell

	REAL*8, ALLOCATABLE	:: AqVolConc(:), &
				OrgVolConc(:), AqOrgVolConc(:)

        !Sum of organic and hydrophilic-organic volume concentrations used by
        !the organic shell terms (ShellAndCoreIndices). OrgVolConc and
        !AqOrgVolConc do not depend on wavelength, so these sums are computed
        !once per particle (see below) rather than once per species per band.
        REAL*8  :: OrgVolConcSum, AqOrgVolConcSum

        !Wavelength-dependent complex refractive indices (WAVE, SOOT_REAL,
        !SOOT_IMAG, WATER_REAL, WATER_IMAG, WASO_REAL, WASO_IMAG,
        !SEASALT_REAL, SEASALT_IMAG, SULF_REAL, SULF_IMAG, DUST_REAL,
        !DUST_IMAG) are module-level PARAMETER arrays declared in
        !LagrangianAerosols.f90 (OPAC v3.1a data; see that file for the full
        !citation).
        REAL*8 :: a_w, b_w, n_real, n_imag

	COMPLEX :: SolarShellInd, SolarCoreInd, VAInd
        COMPLEX :: Eshell, EBC, Factor2, EMG, MGInd
	REAL*8	:: CosAng(MaxAng), M1( MaxAng, 2), &
			M2( MaxAng, 2 ), S21( MaxAng, 2 ), &
		D21( MaxAng, 2 )

	IF (InParticle%numberofparticles .LE. 1.0e-6) THEN    !No particles, use dummy numbers
		InParticle%ExtEff = 1.0
		InParticle%ScaEff = 1.0
		InParticle%BackScaEff = 1.0
		InParticle%AssymParam = 1.0
                RETURN          
	END IF

	ALLOCATE (AqVolConc(HowManyAqChems+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "AqVolConc at size "//	&
	     TRIM(INT2STR(HowManyAqChems+HowManyAqCations+HowManyAqAnions)))
	ALLOCATE (OrgVolConc(HowManyOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate OrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyOrgChems)))
	ALLOCATE (AqOrgVolConc(HowManyAqOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate AqOrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))

	CALL OpticalVolumeConcentrations (InParticle, AqVolConc, OrgVolConc, AqOrgVolConc, &
	                                  SolutionVolConc, ShellVolumeConc)

	!Index 1 of OrgVolConc is BC (the core), so the insoluble organics
	!shell term starts at index 2.
	!(explicit loops, matching the order of the original per-species sums)
	OrgVolConcSum = 0.0
	DO I = 2, HowManyOrgChems
		OrgVolConcSum = OrgVolConcSum + OrgVolConc(I)
	END DO
	AqOrgVolConcSum = 0.0
	DO I = 1, HowManyAqOrgChems
		AqOrgVolConcSum = AqOrgVolConcSum + AqOrgVolConc(I)
	END DO

        !Calculate shell and core refractive index for each wavelength bin
        DO ib = 1, 451 !For photolysis between 250 and 700 nm, every 1 nm

           Wavelength = (250+(ib-1)*1) !in nm

           !Index of point on array WAVE closest to, but less than,
           !the wavelength for interpolation
           w_ind_low = INT(ib/50)+1

           CALL ShellAndCoreIndices (InParticle, ib, Wavelength, w_ind_low, AqVolConc, OrgVolConc, &
                                     AqOrgVolConc, SolutionVolConc, ShellVolumeConc, a_w, b_w, &
                                     SumEr, SumEi, SolarShellInd, SolarCoreInd, HasShell, &
                                     OrgVolConcSum, AqOrgVolConcSum)

           IF (HasShell) THEN
            !!For volume average dielectric constant mixing rule
	    n_real = a_w*SOOT_REAL(w_ind_low)+b_w*SOOT_REAL(w_ind_low+1)
            n_imag = -1.0*(a_w*SOOT_IMAG(w_ind_low)+b_w*SOOT_IMAG(w_ind_low+1))
            ErBC = n_real**2 - n_imag**2
            EiBC = 2*n_real*n_imag
            BCVolRat = (InParticle%AbsCoreRad/InParticle%EffectiveRadius)**3.0
            ErVA = ErBC*BCVolRat+SumEr*(1.0-BCVolRat)
            EiVA = EiBC*BCVolRat+SumEi*(1.0-BCVolRat)
            Factor = SQRT(ErVA**2 + EiVA**2)
	    InParticle%VARealRefracInd(ib) = SQRT((Factor + ErVA)/2.0)
	    InParticle%VAImagRefracInd(ib) = SQRT((Factor - ErVA)/2.0)            
            VAInd = CMPLX(InParticle%VARealRefracInd(ib),-1.0*InParticle%VAImagRefracInd(ib))

            !!For Maxwell-Garnett mixing rule
            Eshell =  CMPLX(SumEr,-1.0*SumEi)
            EBC =  CMPLX(ErBC,-1.0*EiBC)
            Factor2 = (EBC-Eshell)/(EBC+2.0*Eshell)
            EMG = Eshell*(1.0+(3.0*BCVolRat*Factor2/(1.0-BCVolRat*Factor2)))
            ErMG = REAL(EMG)
            EiMG = AIMAG(EMG)
            Factor = SQRT(ErMG**2 + EiMG**2)
	    InParticle%MGRealRefracInd(ib) = SQRT((Factor + ErMG)/2.0)
	    InParticle%MGImagRefracInd(ib) = SQRT((Factor - ErMG)/2.0)            
            MGInd = CMPLX(InParticle%MGRealRefracInd(ib),-1.0*InParticle%MGImagRefracInd(ib))
           ENDIF

	  NumAng = 0.

          !DMiLay needs 2*Pi divided by the wavelength
          Wavelength = 2*PI/((Wavelength)*nm)

!    DMiLay Documentation				
!    I N P U T   A R G U M E N T S

!    (Definition:  size parameter = sphere circumference / wavelength )

!      Rshell      radius of shell

!      Rcore       radius of core

!      WVNO        2*pi / wavelength

!      RindSh      COMPLEX refractive index of shell (negative
!                     imaginary part)

!      RindCo      COMPLEX refractive index of core (negative
!                     imaginary part)

!      MU          array of cosines of scattering angles (angles between
!                     directions of incident and scattered radiation).
!                     For angles between 90 and 180 degrees, use the
!                     supplement (180-angle) of the angle instead, so
!                     0.le.MU.le.1 (see comments below on M1,M2,21,D21)

!      NumAng      Number of scattering angles for which computations
!                     are required; should not exceed MaxAng
!                     (NOTE:  NumAng=0 will suppress the calculation
!                      of the scattering matrix quantitities  M1, M2,
!                      S21, D21 and save a lot of computer time)

!      MaxAng      First dimension of M1,M2,21,D21 in calling program



!    O U T P U T   A R G U M E N T S

!      (Definitions for these arguments can be found in H.C. van de
!       Hulst, Light Scattering By Small Particles, Dover Press, New
!       York, 1981 (reprint of 1957 edition); abbreviated VDH below)

!      QEXT     Efficiency factor for extinction (VDH Sec 9.32)
!               (same as corresponding quantity in MIEV)

!      Qsca     Efficiency factor for scattering (VDH Sec 9.32)
!               (same as corresponding quantity in MIEV)

!      GQS!     average(cosine theta) * Qsca (VDH Sec 9.32)
!                  (<cos theta> is usually denoted by g, hence
!                   the name of the variable)
!               (same as corresponding quantity in MIEV)

!      QBS      Backscatter cross section.
!               ( Re(SBACK)**2 + Im(SBACK)**2 ) / (Pi*XSHELL**2)
!               where the corresponding quantity from MIEV is

!               SBACK = 0.5*sum(n=1 to inf)((-1)**(n+1)(2n+1)(an-bn))

!               and an,bn are ACOE,BCOE below.

!      M1(j,k)  Element M1 of scattering matrix F' (VDH Sec 5.14);
!                  M1(j,1) refers to angle with cosine MU(j); 
!                  M1(j,2) refers to supplement of that angle.
!               (Be sure to type REAL in calling program.)

!      M2(j,k)  Element M2 of scattering matrix F' (VDH Sec 5.14);
!                  M2(j,1) refers to angle with cosine MU(j); 
!                  M2(j,2) refers to supplement of that angle.
!               (Be sure to type REAL in calling program.)

!     S21(j,k)  Element S21 of scattering matrix F' (VDH Sec 5.14);
!                  S21(j,1) refers to angle with cosine MU(j); 
!                  S21(j,2) refers to supplement of that angle.

!     D21(j,k)  Element D21 of scattering matrix F' (VDH Sec 5.14);
!                  D21(j,1) refers to angle with cosine MU(j); 
!                  D21(j,2) refers to supplement of that angle.

!          WRITE(*,*) "Call 1"
	  CALL DMiLay( InParticle%AbsCoreRad, &
			InParticle%EffectiveRadius, &
			Wavelength, &
			SolarShellInd, & 
			SolarCoreInd, &
			CosAng, &
			NumAng, &
			InParticle%ExtEff(ib), &
			InParticle%ScaEff(ib), &
			InParticle%BackScaEff(ib), &
			InParticle%AssymParam(ib), & 
                        M1, M2, S21, D21, MAXANG )

	   !IF (ib .eq. 419 .or. ib .eq. 420) THEN
                !WRITE(*,*) "ExtEff: ", InParticle%ExtEff(ib)
                !WRITE(*,*) "ScaEff: ", InParticle%ScaEff(ib)
                !WRITE(*,*) "SSA: ", InParticle%ScaEff(ib)/InParticle%ExtEff(ib)
                !WRITE(*,*) "AssymParam: ", InParticle%AssymParam(ib)
           !ENDIF

	   !DMiLay gives the assymetry parameter times the scattering efficiency
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParam(ib) = InParticle%AssymParam(ib)/InParticle%ScaEff(ib)

          !Now call for a hypothetical particle of
          !same effective size but all BC
!          WRITE(*,*) "Call 2"
	  CALL MieHomogeneousSphere( InParticle%EffectiveRadius, Wavelength, SolarCoreInd, &
			InParticle%ExtEffBC(ib), &
			InParticle%ScaEffBC(ib), &
			InParticle%BackScaEffBC(ib), &
			InParticle%AssymParamBC(ib) )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParamBC(ib) = &
                InParticle%AssymParamBC(ib)/InParticle%ScaEffBC(ib)

          !Now call for a hypothetical particle of
          !same effective size but all shell
!          WRITE(*,*) "Call 3"
	  CALL MieHomogeneousSphere( InParticle%EffectiveRadius, Wavelength, SolarShellInd, &
			InParticle%ExtEffShell(ib), &
			InParticle%ScaEffShell(ib), &
			InParticle%BackScaEffShell(ib), &
			InParticle%AssymParamShell(ib) )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParamShell(ib) = &
                InParticle%AssymParamShell(ib)/InParticle%ScaEffShell(ib)

!          WRITE(*,*) "Call 4"
          !Now call for volume average dielectric constant mixing rule
	  CALL MieHomogeneousSphere( InParticle%EffectiveRadius, Wavelength, VAInd, &
			InParticle%ExtEffVA(ib), &
			InParticle%ScaEffVA(ib), &
			InParticle%BackScaEffVA(ib), &
			InParticle%AssymParamVA(ib) )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParamVA(ib) = &
                InParticle%AssymParamVA(ib)/InParticle%ScaEffVA(ib)

          !Now call for Maxwell-Garnett mixing rule
	  CALL MieHomogeneousSphere( InParticle%EffectiveRadius, Wavelength, MGInd, &
			InParticle%ExtEffMG(ib), &
			InParticle%ScaEffMG(ib), &
			InParticle%BackScaEffMG(ib), &
			InParticle%AssymParamMG(ib) )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParamMG(ib) = &
                InParticle%AssymParamMG(ib)/InParticle%ScaEffMG(ib)
          
	END DO

        DEALLOCATE(AqVolConc,OrgVolConc,AqOrgVolConc, STAT = I)

	RETURN
END SUBROUTINE ShellRefIndAndRad

!Calculates aerosol optical depth for a grid point at 550 nm.
!Can only be used in inversion mode!

! Modified by CMB (AER): Make Gfortran cooperate with type safety
!                        as described above                         
REAL*8 FUNCTION AerosolOpticalDepth ()
! 08/30/2010 Changed for photolysis radiation bins MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE
	
	Type(Particle),POINTER :: p_Particle

	AerosolOpticalDepth = 0.
	
	p_Particle => Particles%First

	DO WHILE (ASSOCIATED(p_Particle))

           AerosolOpticalDepth = AerosolOpticalDepth + &
                p_Particle%NumberofParticles*&
                    (PI*p_Particle%EffectiveRadius**2)*p_Particle%ExtEff(301)

           p_Particle => p_Particle%Next

	END DO
	
	AerosolOpticalDepth = AerosolOpticalDepth*InversionHeight
	
	RETURN
END FUNCTION AerosolOpticalDepth

!Calculates total aerosol extinction coefficient (m**-1) 
!average single scattering albedo, and 
!average assymetry factor 
!for a grid point at 18 radiative bands.
SUBROUTINE AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, &
                          SubExtCoeff, SubSingScat, Flag)!&
                          !ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                          !SubExtCoeffSep, SubSingScatSep)
! 08/30/2010 Changed for the photolysis wavebands MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
! 05/03/2012 Changed to calculate backscattering coefficient, units
!             for extcoeff and backscatcoeff now Mm^-1
! 08/16/2012 Fixed backscattering coefficient, added submicron extinction & SSA
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added flag to decide between mixing rules
!            Flag = 0, core-in-shell (default)
!            Flag = 1, external mixture of BC and rest
!            Flag = 2, volume average dielectric constant mixing rule
!            Flag = 3, Maxwell-Garnett mixing rule
!
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE

        !Input Variables
        INTEGER :: Flag
	!Output Variables
	REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym, BackScatCoeff
	REAL*8, DIMENSION(451) :: SubExtCoeff, SubSingScat
        !Below are for external mixture of BC and "shell"
!	REAL*8, DIMENSION(451) :: ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep
!	REAL*8, DIMENSION(451) :: SubExtCoeffSep, SubSingScatSep

	!Internal variables
	INTEGER :: ib
	REAL*8, DIMENSION(451) :: ScatCoeff, SubScatCoeff
!	REAL*8, DIMENSION(451) :: ScatCoeffSep, SubScatCoeffSep
	REAL*8 :: AssymSum,xxx, hbs_frac, assym_sq, bc_vol_frac

        ! CMB: Gfortran doesn't like the additional definitions of the 
        !      public variables, nor does it like "Particle" being used
        !      both as a type and a variable name
	Type(Particle),POINTER :: p_Particle

        DO ib = 1, 451
		ExtCoeff(ib) = 0.0
		ScatCoeff(ib) = 0.0
		BackScatCoeff(ib) = 0.0
		SubExtCoeff(ib) = 0.0
		SubScatCoeff(ib) = 0.0
                AssymSum = 0.0
!		ExtCoeffSep(ib) = 0.0
!		ScatCoeffSep(ib) = 0.0
!		BackScatCoeffSep(ib) = 0.0
!		SubExtCoeffSep(ib) = 0.0
!		SubScatCoeffSep(ib) = 0.0
!                AssymSumSep = 0.0
		p_Particle => Particles%First

		DO WHILE (ASSOCIATED(p_Particle))
                     IF (p_Particle%EffectiveRadius .gt. 0) THEN 
                       bc_vol_frac = (p_Particle%AbsCoreRad/&
                               p_Particle%EffectiveRadius)**3
                     ELSE
                       bc_vol_frac = 0.0
                     ENDIF
		     xxx = p_Particle%NumberofParticles*&
                        (PI*p_Particle%EffectiveRadius**2)

                     IF (FLAG .EQ. 0) THEN !Core in Shell  
			ExtCoeff(ib) = ExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
			AssymSum = AssymSum + &
                            xxx*p_Particle%ScaEff(ib)*p_Particle%AssymParam(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
                        ENDIF

                     ELSE IF (Flag .EQ. 1) THEN !External mixture
                        ExtCoeff(ib) = ExtCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ExtEffBC(ib)&
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
			AssymSum = AssymSum + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib)*&
                                p_Particle%AssymParamBC(ib) &
                                       +(1.0-bc_vol_frac)*xxx* &
                                       p_Particle%ScaEffShell(ib)*&
                                       p_Particle%AssymParamShell(ib) 
                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN	
                          SubExtCoeff(ib) = SubExtCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ExtEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
                        ENDIF
	              ELSE IF (FLAG .EQ. 2) THEN !Volume-averaged internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffVA(ib)*p_Particle%AssymParamVA(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
                        ENDIF
                      ELSE IF (FLAG .EQ. 3) THEN !Maxwell-Garnett internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffMG(ib)*&
                            p_Particle%AssymParamMG(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
                        ENDIF
                      ENDIF

                      p_Particle => p_Particle%Next
		END DO

		!Formula from p.327, eqn 9.103 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		SingScat(ib) = ScatCoeff(ib)/ExtCoeff(ib)
		SubSingScat(ib) = SubScatCoeff(ib)/SubExtCoeff(ib)                
!		SingScatSep(ib) = ScatCoeffSep(ib)/ExtCoeffSep(ib)
!		SubSingScatSep(ib) = SubScatCoeffSep(ib)/SubExtCoeffSep(ib)   
		!IF(ib .EQ. 419 .OR. ib .EQ. 420) WRITE(*,*) "SingScat: ", SingScat(ib), ScatCoeff(ib), ExtCoeff(ib) 
		
		!Formula from p.329, eqn 9.117 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		Assym(ib) = AssymSum/ScatCoeff(ib)
!		AssymSep(ib) = AssymSumSep/ScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "Assym: ", Assym(1)

		!Convert ExtCoeff and BackScatCoeff from 1/cm to 1/Mm
		ExtCoeff(ib) = 1.0e8*ExtCoeff(ib)
		SubExtCoeff(ib) = 1.0e8*SubExtCoeff(ib)
!		ExtCoeffSep(ib) = 1.0e8*ExtCoeffSep(ib)
!		SubExtCoeffSep(ib) = 1.0e8*SubExtCoeffSep(ib)


		!hbs_frac = hemispherical backscattering fraction
                !See Wiscombe and Grams (1976) J. Atmos. Sci. 33(12), pp. 2440-2451.
                !Integrate Eq. 23 using Henyey-Greenstein phase function from Eq. 16
                assym_sq = assym(ib)**2
                hbs_frac = 0.5*((assym_sq-1)/assym(ib))*(1.0/SQRT(assym_sq+2*assym(ib)+1)-1.0/SQRT(assym_sq+1))
		BackScatCoeff(ib) = hbs_frac*ScatCoeff(ib)
                BackScatCoeff(ib) = 1.0e8*BackScatCoeff(ib)

!                assym_sq = assymSep(ib)**2
!                hbs_frac = 0.5*((assym_sq-1)/assymSep(ib))*(1.0/SQRT(assym_sq+2*assymSep(ib)+1)-1.0/SQRT(assym_sq+1))
!		BackScatCoeffSep(ib) = hbs_frac*ScatCoeffSep(ib)
!                BackScatCoeffSep(ib) = 1.0e8*BackScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "ExtCoeff: ", ExtCoeff(1)

	END DO
	
	RETURN
END SUBROUTINE AerosolOptProp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FastTUV aerosol optics for SAM-ASP-TUV                                     !!
!!                                                                           !!
!! ShellRefIndAndRadFASTTUV and AerosolOptPropFASTTUV compute the particle   !!
!! and population optical properties at the 7 tropospheric Fast-TUV         !!
!! wavelength bins (297.7-736.3 nm) from a lookup table of DMiLay results   !!
!! by wavelength bin, effective-radius bin and core/shell radius ratio,     !!
!! instead of the full-spectrum Mie calculation in ShellRefIndAndRad and    !!
!! AerosolOptProp. They are used by StepASP4SAM (SAM-ASP-TUV), which passes !!
!! the 7-bin extinction, single scattering albedo and asymmetry to TUV.     !!
!! Moved here from the ASP copy bundled in SAM-ASP-TUV, where they replaced !!
!! ShellRefIndAndRad (file ParticleAttributes-Table-3.h there).             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ShellRefIndAndRadFASTTUV(InParticle)
! 08/30/2010 Changed for the Photolyis Wavebands MJA
! 02/27/2012 Changed to wavelength-dependent refractive indices
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added Maxwell Garnett mixing rule
	USE ModelParameters
	USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, ERROR
	
	IMPLICIT NONE
	
	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
 	INTEGER	:: I, ib, sizebin, ratiobin
	REAL*8	:: radratio, SumEr, SumEi, a_w, b_w, ShellVolumeConc, SolutionVolConc
	LOGICAL :: HasShell
	COMPLEX :: SolarShellInd, SolarCoreInd
	REAL*8, ALLOCATABLE	:: AqVolConc(:), OrgVolConc(:), AqOrgVolConc(:)

	!For each Fast-TUV wavelength bin, the OPAC table entry (WAVE) just below its centre
	INTEGER, PARAMETER :: FastTUVWaveLow(NFastTUVWvl) = (/ 1, 2, 2, 3, 4, 8, 10 /)

	!Upper edges of the effective-radius bins (cm); bin 1 starts above 0
	REAL*8, PARAMETER :: MieSizeEdge(10) = (/ 2.5e-06, 4.0e-06, 6.3e-06, 1.0e-05, 1.59e-05, &
	                                          2.52e-05, 4.01e-05, 6.36e-05, 1.01e-04, 1.0e5 /)
	!Upper edges of the core/shell radius-ratio bins; bin 1 starts above 1e-09, and ratios
	!above 0.99 (cores capped at 0.9999 of the radius, dry particles at 0.99) use bin 5
	REAL*8, PARAMETER :: MieRatioEdge(5) = (/ 0.2, 0.4, 0.6, 0.8, 0.99 /)

	!DMiLay results by quantity (ExtEff, ScaEff, BackScaEff, AssymParam * ScaEff), core/shell
	!radius-ratio bin, effective-radius bin and Fast-TUV wavelength bin
	!Wavelength bin 1: per size bin (1-10) and ratio bin (1-5): ExtEff, ScaEff, BackScaEff, g*ScaEff
	REAL*8, PARAMETER :: MieLUTWvl1(4,5,10) = RESHAPE((/ &
		6.338627489052332E-003, 1.241488234603144E-003, 1.433321141549675E-004, 1.728520368530504E-005, &	! size 1, ratio 1
		1.248333996193211E-002, 1.267648706378087E-003, 1.463999009108589E-004, 1.747987158759047E-005, &	! size 1, ratio 2
		3.545044720169458E-002, 1.377484024610864E-003, 1.592507689480405E-004, 1.840337212544841E-005, &	! size 1, ratio 3
		8.452503919162944E-002, 1.680195179957653E-003, 1.944862060809522E-004, 2.162433372144669E-005, &	! size 1, ratio 4
		0.162185165001693, 2.364855797138718E-003, 2.733733945103151E-004, 3.193401515405426E-005, &	! size 1, ratio 5
		7.324842246942065E-002, 5.714947852366278E-002, 5.396783139919124E-003, 5.307740447917354E-003, &	! size 2, ratio 1
		9.549220396408416E-002, 5.848915239976296E-002, 5.540367277664674E-003, 5.367213844836445E-003, &	! size 2, ratio 2
		0.178405027599639, 6.393871254897240E-002, 6.113377856825348E-003, 5.651958365539543E-003, &	! size 2, ratio 3
		0.353384365393497, 7.844087674108738E-002, 7.567518690627478E-003, 6.685904844738792E-003, &	! size 2, ratio 4
		0.625060843079395, 0.109054267629784, 1.035389876492209E-002, 9.995562341826034E-003, &	! size 2, ratio 5
		0.351291783830091, 0.319748227133345, 1.898839944583400E-002, 7.789879551901473E-002, &	! size 3, ratio 1
		0.402493239923126, 0.325522247831505, 1.962465127322233E-002, 7.811687181908000E-002, &	! size 3, ratio 2
		0.584608494330860, 0.347368523258932, 2.174915876805711E-002, 8.009690764578657E-002, &	! size 3, ratio 3
		0.938110307991438, 0.404399010676221, 2.554968452519550E-002, 9.233244460421390E-002, &	! size 3, ratio 4
		1.43868408414879, 0.515062275233131, 2.879153764666507E-002, 0.133777979016359, &	! size 3, ratio 5
		1.37342361010497, 1.30034385267227, 3.763082625039913E-003, 0.800965518269872, &	! size 4, ratio 1
		1.43928900487113, 1.27608250244860, 6.069771997363317E-003, 0.778293600986693, &	! size 4, ratio 2
		1.65044041551035, 1.17975775055226, 9.242745474808526E-003, 0.704765161285258, &	! size 4, ratio 3
		1.99089397695059, 1.06950509863887, 3.094153704386623E-003, 0.641970780289727, &	! size 4, ratio 4
		2.40116615784776, 1.11349947132658, 4.742024510115556E-004, 0.665018755374916, &	! size 4, ratio 5
		3.38575015020377, 3.25077565269196, 1.508446412793669E-002, 2.37766688732998, &	! size 5, ratio 1
		3.27409021012275, 2.98143230924048, 5.598333667586948E-003, 2.10816040360092, &	! size 5, ratio 2
		3.01731868016347, 2.34486286815046, 3.183202528173745E-003, 1.60144535705915, &	! size 5, ratio 3
		2.70047140243071, 1.54161483767586, 1.651629848283973E-003, 1.14005705699940, &	! size 5, ratio 4
		2.76456387593862, 1.33122951526543, 3.571474194449051E-003, 1.02423843119691, &	! size 5, ratio 5
		4.18136701800153, 3.97025755785278, 0.108606663672993, 2.91867188270452, &	! size 6, ratio 1
		4.06394415032188, 3.67756794759130, 9.350836275077198E-002, 2.66177385186489, &	! size 6, ratio 2
		3.61794215147068, 2.85166425306824, 5.061092198825547E-002, 1.97047089335660, &	! size 6, ratio 3
		2.84701689988633, 1.57140918123535, 1.004723256882462E-002, 1.21514996089066, &	! size 6, ratio 4
		2.54452505400219, 1.12185830190530, 2.046874128868622E-003, 0.979975694893448, &	! size 6, ratio 5
		1.87098454155520, 1.57262814971042, 0.179455739020470, 0.826587691449864, &	! size 7, ratio 1
		2.11331735314581, 1.64429534097236, 0.184961775133334, 0.893850657724618, &	! size 7, ratio 2
		2.65893542470647, 1.83796372572409, 0.169899454145463, 1.21605889053948, &	! size 7, ratio 3
		2.83157834003745, 1.63313473394288, 1.550668915559024E-002, 1.31003612229362, &	! size 7, ratio 4
		2.42950597324128, 1.07532366873380, 7.534829502916470E-004, 1.01785913956035, &	! size 7, ratio 5
		2.75273381059395, 2.29287661020004, 0.332542993634415, 1.87727764549639, &	! size 8, ratio 1
		2.78522801824329, 2.17747302639268, 0.384586401022679, 1.76367656922482, &	! size 8, ratio 2
		2.32899233368725, 1.43454170410489, 0.300743383637634, 1.11524168831126, &	! size 8, ratio 3
		2.38808091024147, 1.17829949927142, 0.110810966562384, 1.05724712475111, &	! size 8, ratio 4
		2.41621435374200, 1.14849985328176, 4.253146086490113E-003, 1.07742063648435, &	! size 8, ratio 5
		2.44858113182375, 1.89873736543775, 1.756280136843464E-002, 1.55333735008986, &	! size 9, ratio 1
		2.55761296962624, 1.89816188211727, 6.173942284503940E-002, 1.55863693680222, &	! size 9, ratio 2
		2.11680519282049, 1.24821332222136, 8.819402620649464E-002, 0.960133085076927, &	! size 9, ratio 3
		2.41818605856077, 1.32259918890826, 4.781133523560435E-003, 1.23957418310785, &	! size 9, ratio 4
		2.36876691735927, 1.28176946657756, 9.201960147680920E-003, 1.13906085474838, &	! size 9, ratio 5
		2.13297337686945, 1.35787060904386, 1.548551042971653E-003, 1.21616908455416, &	! size 10, ratio 1
		2.20928800758987, 1.39399897848554, 8.294012329455452E-003, 1.25407274982272, &	! size 10, ratio 2
		2.22156244494289, 1.32399449175411, 6.473783882947021E-003, 1.19914631263325, &	! size 10, ratio 3
		2.18138867980128, 1.19396817796879, 4.716963136782338E-003, 1.11352663602915, &	! size 10, ratio 4
		2.18569477923327, 1.20660859994311, 3.016967617235626E-003, 1.10586310013231 /), &	! size 10, ratio 5
		(/ 4, 5, 10 /))
	!Wavelength bin 2: per size bin (1-10) and ratio bin (1-5): ExtEff, ScaEff, BackScaEff, g*ScaEff
	REAL*8, PARAMETER :: MieLUTWvl2(4,5,10) = RESHAPE((/ &
		5.139038309134193E-003, 1.062224852832295E-003, 1.229431498678993E-004, 1.368792181601226E-005, &	! size 1, ratio 1
		1.103167815585350E-002, 1.085095451234660E-003, 1.256286612884928E-004, 1.384556118672044E-005, &	! size 1, ratio 2
		3.304657208718486E-002, 1.181018035298606E-003, 1.368679159545767E-004, 1.459447195855821E-005, &	! size 1, ratio 3
		8.003271280431039E-002, 1.444726916079719E-003, 1.676181042182596E-004, 1.720509234474033E-005, &	! size 1, ratio 4
		0.154227758362364, 2.039345294924047E-003, 2.363008157370650E-004, 2.554250517508710E-005, &	! size 1, ratio 5
		6.174351265724396E-002, 4.900113880311652E-002, 4.714337877280814E-003, 4.209664958810828E-003, &	! size 2, ratio 1
		8.264246166699793E-002, 5.017934629863428E-002, 4.841424949812438E-003, 4.258892657080572E-003, &	! size 2, ratio 2
		0.160648446320607, 5.497998991483244E-002, 5.350354387040322E-003, 4.493809289003795E-003, &	! size 2, ratio 3
		0.325671135931675, 6.775085956735458E-002, 6.647893983457620E-003, 5.338322497605263E-003, &	! size 2, ratio 4
		0.582615408336245, 9.476760901990416E-002, 9.165054441082652E-003, 8.028235353939209E-003, &	! size 2, ratio 5
		0.305464876517266, 0.280575984311160, 1.789571418393630E-002, 6.261560271941657E-002, &	! size 3, ratio 1
		0.353870970813702, 0.286290024924826, 1.850049093473012E-002, 6.291873757478936E-002, &	! size 3, ratio 2
		0.527508642710369, 0.308021878809883, 2.058334313826214E-002, 6.496036939386565E-002, &	! size 3, ratio 3
		0.868170462402061, 0.363582058736275, 2.456105876110214E-002, 7.564743502150444E-002, &	! size 3, ratio 4
		1.35450097228911, 0.470022020429068, 2.864258940056085E-002, 0.111199968822955, &	! size 3, ratio 5
		1.20301001368349, 1.14643428189651, 3.728795170429543E-003, 0.683668626620534, &	! size 4, ratio 1
		1.27324128141990, 1.13052692814481, 5.953346666972397E-003, 0.666327205880667, &	! size 4, ratio 2
		1.50912825191453, 1.07089188764741, 9.428179277320462E-003, 0.615974866884009, &	! size 4, ratio 3
		1.91045026916680, 1.01633546416832, 4.507098562386480E-003, 0.587384510776133, &	! size 4, ratio 4
		2.36268918025467, 1.08522888394559, 1.196208984032071E-003, 0.626970551154550, &	! size 4, ratio 5
		3.17683312725760, 3.06691675151759, 1.638071099334255E-002, 2.21469171886986, &	! size 5, ratio 1
		3.10316297314846, 2.83578663035769, 6.341459097105847E-003, 1.99640121691916, &	! size 5, ratio 2
		2.92370396649955, 2.26842320042052, 3.443317796962763E-003, 1.55074999195171, &	! size 5, ratio 3
		2.70807959371097, 1.54648151179000, 2.312392736724962E-003, 1.13771473209496, &	! size 5, ratio 4
		2.77208430333744, 1.34549800213284, 6.126951669554335E-003, 1.01984374784182, &	! size 5, ratio 5
		4.23984112710875, 4.05787742959576, 0.146014487633035, 2.95175309258135, &	! size 6, ratio 1
		4.08513798438302, 3.72078037219820, 0.125042913708215, 2.63671204196320, &	! size 6, ratio 2
		3.61430185207521, 2.86519233639460, 7.132297716563619E-002, 1.92704590242080, &	! size 6, ratio 3
		2.89251558294860, 1.58849906905663, 1.408075881283735E-002, 1.22953821185668, &	! size 6, ratio 4
		2.56567682815404, 1.13796710201124, 4.484374058573750E-004, 0.984674160656133, &	! size 6, ratio 5
		2.13135989024255, 1.85793408755914, 0.241193857011960, 0.971547494130008, &	! size 7, ratio 1
		2.40289721652576, 1.94844833740984, 0.269022640075835, 1.07243353795922, &	! size 7, ratio 2
		2.98625319288557, 2.15647180032080, 0.216819365604938, 1.41424419435123, &	! size 7, ratio 3
		2.98413938438816, 1.72378363078381, 2.513497377169763E-003, 1.38572260018926, &	! size 7, ratio 4
		2.43456544414595, 1.07471266757135, 9.341586908463232E-004, 1.01384372016860, &	! size 7, ratio 5
		2.98515002132892, 2.58582574596983, 0.215299572699643, 2.03039850011861, &	! size 8, ratio 1
		2.88156975776157, 2.31574174515928, 0.265991554779424, 1.78625929718757, &	! size 8, ratio 2
		2.33874988836856, 1.46334441697541, 0.215984644393330, 1.03370071092364, &	! size 8, ratio 3
		2.41403552508252, 1.22434799971637, 2.522146742431766E-002, 1.11515601761215, &	! size 8, ratio 4
		2.41343216825036, 1.13517179373856, 3.398463151616537E-003, 1.07000209174110, &	! size 8, ratio 5
		2.57288363479180, 2.10963425218824, 1.204908203279876E-002, 1.73101183532299, &	! size 9, ratio 1
		2.49352858893840, 1.90095836014293, 5.480643006265104E-002, 1.53622154748801, &	! size 9, ratio 2
		2.09685865156925, 1.27034787442852, 8.413073035068880E-002, 0.957746167416338, &	! size 9, ratio 3
		2.38692448995987, 1.29576188661496, 1.134166856340495E-003, 1.20948776303605, &	! size 9, ratio 4
		2.38285108786522, 1.28014803616328, 8.369357332960934E-003, 1.13828009609778, &	! size 9, ratio 5
		2.23764411470239, 1.53567895609872, 3.392344224483384E-002, 1.36350832165984, &	! size 10, ratio 1
		2.09578050154409, 1.33898693533697, 4.483088457007030E-002, 1.16858182418855, &	! size 10, ratio 2
		2.10355282540924, 1.23516698011055, 4.588265698179459E-002, 1.08635277102341, &	! size 10, ratio 3
		2.24895160836331, 1.26160576658518, 4.644847274671248E-003, 1.17875993716182, &	! size 10, ratio 4
		2.18383476886376, 1.20213075138530, 5.250066425087228E-003, 1.08774762885791 /), &	! size 10, ratio 5
		(/ 4, 5, 10 /))
	!Wavelength bin 3: per size bin (1-10) and ratio bin (1-5): ExtEff, ScaEff, BackScaEff, g*ScaEff
	REAL*8, PARAMETER :: MieLUTWvl3(4,5,10) = RESHAPE((/ &
		4.263139559487195E-003, 8.678188365096641E-004, 1.007407481705718E-004, 1.011486858243403E-005, &	! size 1, ratio 1
		9.824738790186843E-003, 8.865980196457613E-004, 1.029492819769158E-004, 1.023204761827078E-005, &	! size 1, ratio 2
		3.059731242702119E-002, 9.653255345575724E-004, 1.121898933991893E-004, 1.078930308842013E-005, &	! size 1, ratio 3
		7.490825388297630E-002, 1.181492663995140E-003, 1.374515268544293E-004, 1.273109799514702E-005, &	! size 1, ratio 4
		0.144812679292337, 1.668311505157806E-003, 1.938561160319200E-004, 1.892483759965921E-005, &	! size 1, ratio 5
		5.057006973339805E-002, 4.010205846162691E-002, 3.943011051862751E-003, 3.115401483151503E-003, &	! size 2, ratio 1
		6.977025758065473E-002, 4.107031545316867E-002, 4.048318605115704E-003, 3.152597047211763E-003, &	! size 2, ratio 2
		0.141512475631317, 4.502610859996364E-002, 4.472038199127779E-003, 3.329597172343329E-003, &	! size 2, ratio 3
		0.293670507177789, 5.557005162808929E-002, 5.560477419200620E-003, 3.960465083758706E-003, &	! size 2, ratio 4
		0.531433841467238, 7.800226375329375E-002, 7.708573944335241E-003, 5.963814556042968E-003, &	! size 2, ratio 5
		0.255624221481006, 0.235302311794726, 1.622437278838499E-002, 4.699747032971522E-002, &	! size 3, ratio 1
		0.299901356490615, 0.240498741906031, 1.676200456901536E-002, 4.731412025801772E-002, &	! size 3, ratio 2
		0.460186402890855, 0.260466841879183, 1.867559148941026E-002, 4.915210649956857E-002, &	! size 3, ratio 3
		0.779050110889278, 0.311184486362578, 2.258814091761381E-002, 5.766194544095191E-002, &	! size 3, ratio 4
		1.23996624514364, 0.408463054261145, 2.734083581546594E-002, 8.561118615865378E-002, &	! size 3, ratio 5
		1.01940542275800, 0.974707004521188, 6.393422006093996E-003, 0.543270423557770, &	! size 4, ratio 1
		1.09150328194008, 0.966509038146904, 8.311286263230834E-003, 0.531535831562413, &	! size 4, ratio 2
		1.33852810361431, 0.938731322831333, 1.162230282168859E-002, 0.502398103621807, &	! size 4, ratio 3
		1.78429194952155, 0.936569432400333, 8.035549488415072E-003, 0.506131319667848, &	! size 4, ratio 4
		2.29222835390614, 1.03395440348948, 3.580206860432335E-003, 0.568648666960359, &	! size 4, ratio 5
		2.82373243492915, 2.73627678098817, 2.762871640954166E-002, 1.90456440796543, &	! size 5, ratio 1
		2.79171757067844, 2.55402760711118, 1.875195187120392E-002, 1.75134734125139, &	! size 5, ratio 2
		2.70058090013906, 2.08043513665998, 1.294640567268690E-002, 1.40223289260488, &	! size 5, ratio 3
		2.67205845603090, 1.52431176163280, 5.801715152035907E-003, 1.10659639703817, &	! size 5, ratio 4
		2.78628380135758, 1.36074130619291, 9.037377844402734E-003, 1.01355098208938, &	! size 5, ratio 5
		4.02361663987596, 3.88486261377668, 8.836415432823291E-002, 2.84931045496053, &	! size 6, ratio 1
		3.81636850766442, 3.49034762729065, 7.591212543258335E-002, 2.45544094723994, &	! size 6, ratio 2
		3.30904943558053, 2.60632644752455, 3.168884184523759E-002, 1.69868864562301, &	! size 6, ratio 3
		2.81698823339211, 1.54303889094819, 8.778628905193386E-003, 1.17769066859975, &	! size 6, ratio 4
		2.60591462359312, 1.16408535385012, 1.110688542706662E-003, 0.995936288016512, &	! size 6, ratio 5
		2.13140593071492, 1.91090021632114, 0.260019140372145, 0.952632522132422, &	! size 7, ratio 1
		2.40423282417051, 2.00506467553308, 0.245418726819286, 1.05361578839653, &	! size 7, ratio 2
		2.88402673625340, 2.08256710393457, 0.134760983236147, 1.29297896211798, &	! size 7, ratio 3
		2.90341053850232, 1.66754570953437, 3.716049871002671E-002, 1.28930458133611, &	! size 7, ratio 4
		2.44319477905224, 1.07440846412043, 1.180109235404892E-004, 1.00793025770467, &	! size 7, ratio 5
		2.81760720579818, 2.51729276797066, 0.197318430580763, 2.04044574055271, &	! size 8, ratio 1
		2.57488480960150, 2.10668718533517, 0.201193345762590, 1.63060098489225, &	! size 8, ratio 2
		2.04767334628316, 1.24793885026032, 0.120656751271596, 0.890123772127963, &	! size 8, ratio 3
		2.45837383604067, 1.28060031799697, 5.244784741941551E-002, 1.14765636485615, &	! size 8, ratio 4
		2.41085234555214, 1.12036950583058, 3.010789121369003E-003, 1.06140943853389, &	! size 8, ratio 5
		2.60877262142503, 2.19230972755405, 3.817992804076797E-002, 1.77911080621722, &	! size 9, ratio 1
		2.33891354601800, 1.78239959331719, 9.912359476817148E-002, 1.37900071570137, &	! size 9, ratio 2
		2.21209387722090, 1.39426224157774, 0.130532300116988, 1.05801687815812, &	! size 9, ratio 3
		2.32551374554615, 1.24141997828906, 1.670943852470345E-003, 1.14203276923661, &	! size 9, ratio 4
		2.39805184963928, 1.27305434602451, 9.189984740622667E-003, 1.13577122111321, &	! size 9, ratio 5
		2.14964998564504, 1.51295303925601, 5.334558397666915E-002, 1.31845874388974, &	! size 10, ratio 1
		2.13670642026765, 1.43234784607461, 3.581665921291712E-002, 1.24015288508466, &	! size 10, ratio 2
		2.07566598312208, 1.23342500902051, 5.009895263700193E-002, 1.07244466563173, &	! size 10, ratio 3
		2.22863863650280, 1.23323597862107, 9.968545595441249E-003, 1.1632950291983, &	! size 10, ratio 4
		2.17256371810077, 1.18063013582568, 7.661975129484521E-003, 1.06154342330829 /), &	! size 10, ratio 5
		(/ 4, 5, 10 /))
	!Wavelength bin 4: per size bin (1-10) and ratio bin (1-5): ExtEff, ScaEff, BackScaEff, g*ScaEff
	REAL*8, PARAMETER :: MieLUTWvl4(4,5,10) = RESHAPE((/ &
		2.753364167085392E-003, 4.725331281362864E-004, 5.525631584254712E-005, 4.071316893356675E-006, &	! size 1, ratio 1
		7.425505208959161E-003, 4.827810529546341E-004, 5.646618826504992E-005, 4.118701937006665E-006, &	! size 1, ratio 2
		2.486664441636313E-002, 5.256939654885124E-004, 6.152496371390335E-005, 4.344494468088785E-006, &	! size 1, ratio 3
		6.203756948746719E-002, 6.431720886111547E-004, 7.532868641386177E-005, 5.130076429292264E-006, &	! size 1, ratio 4
		0.120590912593262, 9.070519208908492E-004, 1.061388237443969E-004, 7.627445733771990E-006, &	! size 1, ratio 5
		2.864233227498378E-002, 2.186962766642239E-002, 2.265966881205732E-003, 1.257057690430274E-003, &	! size 2, ratio 1
		4.366756379127423E-002, 2.239247385477439E-002, 2.324102770485585E-003, 1.272431377629956E-003, &	! size 2, ratio 2
		9.988853168620168E-002, 2.454187372940983E-002, 2.560548197627754E-003, 1.345304400266801E-003, &	! size 2, ratio 3
		0.219654871327424, 3.030065847186051E-002, 3.178346543821798E-003, 1.600874676180828E-003, &	! size 2, ratio 4
		0.408067820596703, 4.272803842922383E-002, 4.445949619160706E-003, 2.407762067481284E-003, &	! size 2, ratio 5
		0.147288160966071, 0.134465686275423, 1.097615308225398E-002, 1.948109916119950E-002, &	! size 3, ratio 1
		0.180354610625155, 0.137771689092346, 1.131710925779624E-002, 1.968055173124623E-002, &	! size 3, ratio 2
		0.302347408115515, 0.150817512007478, 1.260948478125482E-002, 2.068048504157519E-002, &	! size 3, ratio 3
		0.553448056892250, 0.184252437596009, 1.559196483684770E-002, 2.456498493945252E-002, &	! size 3, ratio 4
		0.929974307432744, 0.250380937135547, 2.033968194994756E-002, 3.699096772651982E-002, &	! size 3, ratio 5
		0.652413758309375, 0.626191688703548, 1.778549309279416E-002, 0.253968122802817, &	! size 4, ratio 1
		0.721025076713686, 0.630468248902443, 1.893515198426141E-002, 0.251574796453445, &	! size 4, ratio 2
		0.956324844257156, 0.646146505558791, 2.162672347371767E-002, 0.248984659038393, &	! size 4, ratio 3
		1.39820340055190, 0.705373479886065, 2.214629728784083E-002, 0.277005623997083, &	! size 4, ratio 4
		1.98374226343678, 0.840267636255933, 1.760556340296673E-002, 0.366428447602450, &	! size 4, ratio 5
		2.21481592015604, 2.16240977794318, 3.898035959863895E-002, 1.33502856826726, &	! size 5, ratio 1
		2.24271052456984, 2.06800362941874, 3.629342249738617E-002, 1.27799769373332, &	! size 5, ratio 2
		2.24314254153314, 1.70972523872547, 2.928371976972988E-002, 1.06585367139569, &	! size 5, ratio 3
		2.36865146176109, 1.35898568682125, 1.235527109455989E-002, 0.906115632551481, &	! size 5, ratio 4
		2.73420713362149, 1.33278582460604, 9.967977942811165E-003, 0.925948717186274, &	! size 5, ratio 5
		4.13294149320973, 4.03032111503355, 7.360828078472970E-002, 2.95878817637063, &	! size 6, ratio 1
		3.93657977690474, 3.66770358783291, 6.378862352038903E-002, 2.61411380470694, &	! size 6, ratio 2
		3.45846369760482, 2.76610072707306, 2.661771858894812E-002, 1.88735403979881, &	! size 6, ratio 3
		2.86381723273350, 1.60383375935866, 5.473898240651396E-003, 1.22698810833429, &	! size 6, ratio 4
		2.69304955122330, 1.25033059908016, 6.660793655015043E-003, 1.02464397830354, &	! size 6, ratio 5
		2.95525036302991, 2.81071944848402, 0.144722870003543, 1.79750526281799, &	! size 7, ratio 1
		3.07622075305655, 2.74767979573311, 0.114086144198558, 1.75771823060325, &	! size 7, ratio 2
		3.13443353727213, 2.38962795526458, 6.656313299481739E-002, 1.56192006850129, &	! size 7, ratio 3
		2.81932474143716, 1.59110473436932, 1.192733327901523E-002, 1.21178915781235, &	! size 7, ratio 4
		2.47216837452084, 1.07757752023344, 1.080634786921352E-003, 0.987493204096343, &	! size 7, ratio 5
		2.27724011917898, 2.03534636405276, 0.296937265639801, 1.33451093586322, &	! size 8, ratio 1
		2.06821111017703, 1.64151843375264, 0.270707651280783, 0.953114456017294, &	! size 8, ratio 2
		2.07020163745579, 1.29124608353427, 0.170095675320752, 0.702221463350847, &	! size 8, ratio 3
		2.72568302842573, 1.51534078721248, 4.813809898166456E-002, 1.25743137980625, &	! size 8, ratio 4
		2.41019470324832, 1.09009441391071, 1.305040745713129E-003, 1.04080910045869, &	! size 8, ratio 5
		2.04176676954650, 1.74112647693926, 4.227385751022267E-002, 1.29350216586826, &	! size 9, ratio 1
		2.12983639968120, 1.66122289119653, 6.941325831696683E-002, 1.23517403923357, &	! size 9, ratio 2
		2.65516028103601, 1.88314938012702, 5.469645993818988E-002, 1.54336518258899, &	! size 9, ratio 3
		2.13538718383939, 1.05843352518967, 1.092278492080876E-004, 0.929098956036963, &	! size 9, ratio 4
		2.42230784450767, 1.22969039581765, 7.302126663631858E-003, 1.11814100390391, &	! size 9, ratio 5
		2.23747218460542, 1.72439047618241, 6.409481896871368E-002, 1.47381594101870, &	! size 10, ratio 1
		2.06562574336866, 1.45518814322546, 1.205943583104343E-002, 1.20942166076782, &	! size 10, ratio 2
		2.10833500508912, 1.31014119463384, 7.973096319157684E-002, 1.10213226816747, &	! size 10, ratio 3
		2.11640145406860, 1.11670002015824, 1.780475953408376E-002, 1.03170971466352, &	! size 10, ratio 4
		2.15041186500276, 1.11393480479371, 7.529827215511730E-003, 1.04393547659690 /), &	! size 10, ratio 5
		(/ 4, 5, 10 /))
	!Wavelength bin 5: per size bin (1-10) and ratio bin (1-5): ExtEff, ScaEff, BackScaEff, g*ScaEff
	REAL*8, PARAMETER :: MieLUTWvl5(4,5,10) = RESHAPE((/ &
		2.156470752780501E-003, 2.423510895145302E-004, 2.850558992752902E-005, 1.497584622399152E-006, &	! size 1, ratio 1
		6.011100304721312E-003, 2.475458340153061E-004, 2.912078589763700E-005, 1.514855723801651E-006, &	! size 1, ratio 2
		2.039744966943594E-002, 2.692705667022480E-004, 3.169077807365558E-005, 1.597228459102876E-006, &	! size 1, ratio 3
		5.105526786922314E-002, 3.285750068981799E-004, 3.868990015467943E-005, 1.883101191165021E-006, &	! size 1, ratio 4
		9.935247288747233E-002, 4.615071333409716E-004, 5.430775616233305E-005, 2.788408085006948E-006, &	! size 1, ratio 5
		1.669458937281732E-002, 1.120209690850707E-002, 1.209098321487180E-003, 4.628219868201358E-004, &	! size 2, ratio 1
		2.837001252781295E-002, 1.146191869364921E-002, 1.238543761667780E-003, 4.684342739501339E-004, &	! size 2, ratio 2
		7.205891157290933E-002, 1.253454013758639E-002, 1.359193147791398E-003, 4.950143410363731E-004, &	! size 2, ratio 3
		0.165305049113183, 1.541858216571215E-002, 1.678107315956416E-003, 5.873438238030849E-004, &	! size 2, ratio 4
		0.312501798369903, 2.170880988992300E-002, 2.350547480252358E-003, 8.776614745490395E-004, &	! size 2, ratio 5
		8.048939535738715E-002, 7.043825104650552E-002, 6.468355714037413E-003, 7.268310940692170E-003, &	! size 3, ratio 1
		0.104420575985344, 7.218307088083935E-002, 6.653368078192959E-003, 7.353345246323822E-003, &	! size 3, ratio 2
		0.193572106344506, 7.919845083348445E-002, 7.380076781413513E-003, 7.762654644333494E-003, &	! size 3, ratio 3
		0.381053745280321, 9.748035902239889E-002, 9.168685947886111E-003, 9.244252821987917E-003, &	! size 3, ratio 4
		0.670302709765434, 0.135022061326670, 1.243710257185139E-002, 1.391330405960398E-002, &	! size 3, ratio 5
		0.399820651387827, 0.379883274074818, 2.020361226829106E-002, 0.103873792829933, &	! size 4, ratio 1
		0.455403541898221, 0.386772613189943, 2.098751310955096E-002, 0.104049427405662, &	! size 4, ratio 2
		0.651763858720752, 0.412126744622607, 2.344456978061815E-002, 0.106475504137779, &	! size 4, ratio 3
		1.02870049467777, 0.476681422317556, 2.713394500818145E-002, 0.123147638691622, &	! size 4, ratio 4
		1.55707178531904, 0.598791865162430, 2.878177527669097E-002, 0.177857561447019, &	! size 4, ratio 5
		1.62099816410536, 1.57416424958781, 9.008312064238813E-003, 0.984514280597314, &	! size 5, ratio 1
		1.68427078358998, 1.53976687599924, 1.108278493238152E-002, 0.955625625132237, &	! size 5, ratio 2
		1.86220978318272, 1.38391969898805, 1.271187494500821E-002, 0.846451907574883, &	! size 5, ratio 3
		2.12280361098686, 1.18823762138111, 3.356462294672066E-003, 0.736746348599184, &	! size 5, ratio 4
		2.50059318282653, 1.20340374509707, 1.311716856288350E-003, 0.746498201158639, &	! size 5, ratio 5
		3.54078596830612, 3.45847685303443, 3.198360178172965E-002, 2.51481322206534, &	! size 6, ratio 1
		3.38232496270729, 3.13603488550063, 2.177095831839554E-002, 2.17920072508563, &	! size 6, ratio 2
		3.05185203148685, 2.41454769078163, 1.359885498022188E-002, 1.61285574825469, &	! size 6, ratio 3
		2.71020889915230, 1.56121405249262, 4.120296335153626E-003, 1.15128600049189, &	! size 6, ratio 4
		2.79775087376873, 1.34534074511176, 1.266605735950606E-003, 1.05237794231973, &	! size 6, ratio 5
		3.79839753863319, 3.67542762340395, 6.231530104948533E-002, 2.64329165893099, &	! size 7, ratio 1
		3.70555289052490, 3.39630630241535, 5.870560457402267E-002, 2.39134165736987, &	! size 7, ratio 2
		3.32135560077935, 2.61225061502338, 1.745859699296794E-002, 1.72559001421307, &	! size 7, ratio 3
		2.76597986500317, 1.52787613178405, 2.020718458469559E-003, 1.15470793231692, &	! size 7, ratio 4
		2.52951896303713, 1.10551449919514, 4.493611108713668E-003, 0.977491723508472, &	! size 7, ratio 5
		1.98269153957613, 1.75070724654365, 0.596312156974750, 0.814896613634202, &	! size 8, ratio 1
		2.13994561469170, 1.71505990064196, 0.561875869574628, 0.811718423690809, &	! size 8, ratio 2
		2.71924620632794, 1.91427232338367, 0.450655248663453, 1.11990749586247, &	! size 8, ratio 3
		2.99844439004540, 1.72913526143464, 7.159558927873695E-002, 1.39133971455852, &	! size 8, ratio 4
		2.42621425011645, 1.07919795900881, 1.832738097184638E-004, 1.02579094906008, &	! size 8, ratio 5
		2.51437206079227, 2.23110602644525, 0.210648672235153, 1.75170270019410, &	! size 9, ratio 1
		2.71636878080886, 2.25775825646040, 0.194005985128595, 1.79854133717791, &	! size 9, ratio 2
		2.45775870854036, 1.66858194943098, 7.243187861157524E-002, 1.31182512900495, &	! size 9, ratio 3
		2.25680078236096, 1.10370970849717, 1.730107992007260E-002, 0.978909629046713, &	! size 9, ratio 4
		2.41970148126958, 1.16623188123296, 4.823501579071294E-003, 1.08759811723919, &	! size 9, ratio 5
		2.37604420245805, 1.90518173079963, 6.002949337120361E-002, 1.59334484865810, &	! size 10, ratio 1
		2.13488567918209, 1.55166581678293, 1.588575944612918E-002, 1.24950632183058, &	! size 10, ratio 2
		2.40555145061762, 1.60765353248846, 9.845384462499929E-002, 1.35685469246507, &	! size 10, ratio 3
		2.29395827727994, 1.28530950905199, 5.914862772105029E-003, 1.18513281110847, &	! size 10, ratio 4
		2.19533391582167, 1.14522369107149, 1.248136068711735E-003, 1.11011499130817 /), &	! size 10, ratio 5
		(/ 4, 5, 10 /))
	!Wavelength bin 6: per size bin (1-10) and ratio bin (1-5): ExtEff, ScaEff, BackScaEff, g*ScaEff
	REAL*8, PARAMETER :: MieLUTWvl6(4,5,10) = RESHAPE((/ &
		1.745990730877233E-003, 7.391174987435961E-005, 8.751096344873538E-006, 2.526572623971385E-007, &	! size 1, ratio 1
		4.440939533503097E-003, 7.546527233444613E-005, 8.935722060557604E-006, 2.555310567509572E-007, &	! size 1, ratio 2
		1.449796849247914E-002, 8.192307596275294E-005, 9.702702563658139E-006, 2.692444256782151E-007, &	! size 1, ratio 3
		3.593629234104791E-002, 9.935145909251058E-005, 1.176989932640274E-005, 3.163832137897422E-007, &	! size 1, ratio 4
		6.974674587368791E-002, 1.380002289731614E-004, 1.634264801108306E-005, 4.634947479846343E-007, &	! size 1, ratio 5
		8.003059671339458E-003, 3.403547228659068E-003, 3.844268516783187E-004, 7.809869126870609E-005, &	! size 2, ratio 1
		1.562305355258488E-002, 3.478762198184714E-003, 3.931456551665856E-004, 7.902165898845314E-005, &	! size 2, ratio 2
		4.411449864602950E-002, 3.789242520650595E-003, 4.289828086256152E-004, 8.339861704319302E-005, &	! size 2, ratio 3
		0.104963438764179, 4.620191164186989E-003, 5.240207479496873E-004, 9.839738585336499E-005, &	! size 2, ratio 4
		0.201198483519751, 6.436651325675930E-003, 7.281489100011894E-004, 1.449333341115026E-004, &	! size 2, ratio 5
		2.958896551755912E-002, 2.161440615214868E-002, 2.241389981346530E-003, 1.235340318836642E-003, &	! size 3, ratio 1
		4.362928843777036E-002, 2.212396060954086E-002, 2.298035921877446E-003, 1.250376698653324E-003, &	! size 3, ratio 2
		9.617834969558516E-002, 2.420001626704581E-002, 2.526221633178067E-003, 1.321561677527488E-003, &	! size 3, ratio 3
		0.208229242962370, 2.966670186588255E-002, 3.111995027755152E-003, 1.567100051923124E-003, &	! size 3, ratio 4
		0.384882292364829, 4.126765149541759E-002, 4.294187520066071E-003, 2.325522236346111E-003, &	! size 3, ratio 5
		0.147289932438826, 0.132262363149071, 1.083488882099055E-002, 1.900024659026915E-002, &	! size 4, ratio 1
		0.178214429058716, 0.135513949560788, 1.116901832652692E-002, 1.920097283927704E-002, &	! size 4, ratio 2
		0.292463937989071, 0.148240061353616, 1.242417396618586E-002, 2.019830645151947E-002, &	! size 4, ratio 3
		0.528374841368533, 0.180210496949550, 1.526841985223957E-002, 2.394712939964825E-002, &	! size 4, ratio 4
		0.883780761355177, 0.242377608224781, 1.973428625847923E-002, 3.563398716507265E-002, &	! size 4, ratio 5
		0.655866354883852, 0.625022675951519, 1.773540365132158E-002, 0.253551650833792, &	! size 5, ratio 1
		0.721449871409231, 0.630362557417295, 1.888538677092272E-002, 0.251659400008581, &	! size 5, ratio 2
		0.947940060498512, 0.649338032078007, 2.157892979403054E-002, 0.250707767889194, &	! size 5, ratio 3
		1.37699128271348, 0.710893696073668, 2.209180435872176E-002, 0.280300293979218, &	! size 5, ratio 4
		1.95333915266329, 0.843706535935239, 1.734450453694417E-002, 0.370309187029721, &	! size 5, ratio 5
		2.21155597637799, 2.15005296254870, 3.842774817988816E-002, 1.32837581497928, &	! size 6, ratio 1
		2.24185436275420, 2.06553593504804, 3.562590509315052E-002, 1.27824338108958, &	! size 6, ratio 2
		2.25389419640571, 1.73367031909003, 2.875169477837651E-002, 1.08150467565139, &	! size 6, ratio 3
		2.38277647808416, 1.39719657679278, 1.276567136944368E-002, 0.927166406690145, &	! size 6, ratio 4
		2.75039792073941, 1.36715119187155, 9.849083369577715E-003, 0.947433996589785, &	! size 6, ratio 5
		4.12200574278341, 4.00174704435463, 7.175774120730871E-002, 2.94243861164802, &	! size 7, ratio 1
		3.93526119047506, 3.65772322316503, 6.286334426448967E-002, 2.61403874572427, &	! size 7, ratio 2
		3.47453652819076, 2.78560024289033, 2.579655482451309E-002, 1.90668298076763, &	! size 7, ratio 3
		2.90209025062548, 1.65039032725206, 4.996210538022457E-003, 1.25721943136780, &	! size 7, ratio 4
		2.71994191970323, 1.27893743131628, 6.444624170997618E-003, 1.04453110948951, &	! size 7, ratio 5
		2.95992740299252, 2.79135476125437, 0.140381758860169, 1.79503415913209, &	! size 8, ratio 1
		3.07041978015654, 2.72580740171228, 0.111857495874945, 1.75174584094881, &	! size 8, ratio 2
		3.11849224338763, 2.36906649499724, 6.687453488181348E-002, 1.54733378731707, &	! size 8, ratio 3
		2.80619569449233, 1.57642238421905, 1.038199592688869E-002, 1.19575528096164, &	! size 8, ratio 4
		2.47860056315199, 1.08036129578137, 1.181422890651614E-003, 0.987927800861476, &	! size 8, ratio 5
		2.26311639619112, 1.98787413805281, 0.262094922520487, 1.32193503285332, &	! size 9, ratio 1
		2.06766010618677, 1.61556765233881, 0.241008695412522, 0.962079949615984, &	! size 9, ratio 2
		2.07358117616456, 1.27957179697689, 0.150505123505038, 0.719201438179157, &	! size 9, ratio 3
		2.71080656248684, 1.49868666323326, 4.339057737712454E-002, 1.24784272622554, &	! size 9, ratio 4
		2.41482825444001, 1.09375391218603, 1.237477201533048E-003, 1.04390195376583, &	! size 9, ratio 5
		2.27037858457026, 1.81605110214765, 5.576530762260055E-002, 1.48495089386339, &	! size 10, ratio 1
		2.05247974159767, 1.47654312653192, 0.125019234074214, 1.15154802069341, &	! size 10, ratio 2
		2.46409375074626, 1.65435627485793, 0.170163724610554, 1.38854140585308, &	! size 10, ratio 3
		2.29582319104450, 1.24745592736670, 1.843349596407969E-002, 1.13258828871739, &	! size 10, ratio 4
		2.29684611193584, 1.26036871862716, 5.638748954638020E-003, 1.13886041054344 /), &	! size 10, ratio 5
		(/ 4, 5, 10 /))
	!Wavelength bin 7: per size bin (1-10) and ratio bin (1-5): ExtEff, ScaEff, BackScaEff, g*ScaEff
	REAL*8, PARAMETER :: MieLUTWvl7(4,5,10) = RESHAPE((/ &
		1.822525031820004E-003, 3.301443988157131E-005, 3.919444865237387E-006, 7.547501792487873E-008, &	! size 1, ratio 1
		3.976783558301614E-003, 3.370752387684337E-005, 4.001932571094586E-006, 7.633503523585032E-008, &	! size 1, ratio 2
		1.201545347492798E-002, 3.658147190747974E-005, 4.343829425194149E-006, 8.043689139507312E-008, &	! size 1, ratio 3
		2.915180842891664E-002, 4.430108466274274E-005, 5.261370472705972E-006, 9.447242553291143E-008, &	! size 1, ratio 4
		5.618368287014076E-002, 6.134127290403799E-005, 7.283375779648253E-006, 1.379899582095978E-007, &	! size 1, ratio 5
		6.350211255988941E-003, 1.516929500974525E-003, 1.745174041177409E-004, 2.332714440848649E-005, &	! size 2, ratio 1
		1.227184437690569E-002, 1.549921682674426E-003, 1.783785489322480E-004, 2.360067697563908E-005, &	! size 2, ratio 2
		3.439992174863217E-002, 1.686098531844482E-003, 1.942696076786347E-004, 2.489889273816199E-005, &	! size 2, ratio 3
		8.164294863008409E-002, 2.049861918610712E-003, 2.364616196182227E-004, 2.932714034242508E-005, &	! size 2, ratio 4
		0.156339333036700, 2.845401925807674E-003, 3.276771175075854E-004, 4.299545171349382E-005, &	! size 2, ratio 5
		1.778645879661125E-002, 9.630624272648958E-003, 1.047199790254237E-003, 3.693591132836390E-004, &	! size 3, ratio 1
		2.817176459979360E-002, 9.850994744658644E-003, 1.072254685596641E-003, 3.738316430915966E-004, &	! size 3, ratio 2
		6.703824339568021E-002, 1.075298035076666E-002, 1.174049194434285E-003, 3.949759463885883E-004, &	! size 3, ratio 3
		0.150070514170024, 1.313740414208335E-002, 1.438868893717062E-003, 4.671294722475529E-004, &	! size 3, ratio 4
		0.281435205955893, 1.825633648332662E-002, 1.990015549345189E-003, 6.887746698199976E-004, &	! size 3, ratio 5
		7.511537565154930E-002, 6.037039703762549E-002, 5.661249225289418E-003, 5.763441149550745E-003, &	! size 4, ratio 1
		9.608087097199536E-002, 6.185787194192215E-002, 5.819799652834955E-003, 5.832797077916601E-003, &	! size 4, ratio 2
		0.174315781725701, 6.780390801172184E-002, 6.439851071284188E-003, 6.163905103196722E-003, &	! size 4, ratio 3
		0.339577525959736, 8.305973787255339E-002, 7.951594246883955E-003, 7.325304521773114E-003, &	! size 4, ratio 4
		0.596353723504103, 0.114051583657035, 1.072263296429436E-002, 1.088757311624777E-002, &	! size 4, ratio 5
		0.367213460887450, 0.338129863781664, 1.942458711762533E-002, 8.548555144248844E-002, &	! size 5, ratio 1
		0.416908850227615, 0.345048113260993, 2.015557361575887E-002, 8.586748845454591E-002, &	! size 5, ratio 2
		0.594332181245844, 0.370529389685645, 2.251890121891220E-002, 8.865139920209882E-002, &	! size 5, ratio 3
		0.939935921136212, 0.432391959185205, 2.636799177643430E-002, 0.103202829981282, &	! size 5, ratio 4
		1.43145499993234, 0.545537007405742, 2.893519960973154E-002, 0.149244833889450, &	! size 5, ratio 5
		1.45816059241186, 1.38995743041681, 4.843928180701920E-003, 0.863975288788868, &	! size 6, ratio 1
		1.52390102005966, 1.36972244072900, 6.940145129539475E-003, 0.843783357047611, &	! size 6, ratio 2
		1.73431959302154, 1.27869184206174, 9.603752900582917E-003, 0.773201646713378, &	! size 6, ratio 3
		2.06314074085207, 1.16139320271965, 2.764357017880969E-003, 0.706554457231506, &	! size 6, ratio 4
		2.45645089565132, 1.19339683407815, 6.017506001740586E-004, 0.722065171180426, &	! size 6, ratio 5
		3.45951821027362, 3.33678681914849, 1.893751813256706E-002, 2.44167514622609, &	! size 7, ratio 1
		3.34310639991427, 3.06831050588359, 1.005370487819418E-002, 2.16582495070475, &	! size 7, ratio 2
		3.07849309828382, 2.43162727982638, 6.015725161101985E-003, 1.65243326082859, &	! size 7, ratio 3
		2.75540223704259, 1.61943232037147, 1.834465461441918E-003, 1.18930306573657, &	! size 7, ratio 4
		2.82286734244914, 1.39140112703455, 2.767346080934591E-003, 1.07134563671328, &	! size 7, ratio 5
		4.05511210350117, 3.86942445930530, 8.359854479853722E-002, 2.83651622731413, &	! size 8, ratio 1
		3.94244140337470, 3.58637491300409, 7.429315165712792E-002, 2.58680033029294, &	! size 8, ratio 2
		3.51060204855541, 2.77177769656355, 3.406103442355219E-002, 1.89292821121038, &	! size 8, ratio 3
		2.80140587662272, 1.54935066760967, 5.461204278651640E-003, 1.17934121311270, &	! size 8, ratio 4
		2.54983639259369, 1.12259967128459, 3.408981499861641E-003, 0.979790120132008, &	! size 8, ratio 5
		1.78841120186596, 1.50285152352527, 0.256009325703574, 0.724721746210302, &	! size 9, ratio 1
		2.00646138691121, 1.55006876500885, 0.240575656297192, 0.770042000363660, &	! size 9, ratio 2
		2.54595095593060, 1.73068609919718, 0.198090514450071, 1.07049084651018, &	! size 9, ratio 3
		2.83394281860531, 1.61873736873039, 3.734013425565764E-002, 1.28031725207590, &	! size 9, ratio 4
		2.43717515107114, 1.08346654076758, 5.314656322385932E-004, 1.02602074598043, &	! size 9, ratio 5
		2.45011338214344, 1.95615262558493, 3.891223075565580E-002, 1.59939995051104, &	! size 10, ratio 1
		2.49067116114546, 1.87809892686616, 0.105901825728398, 1.52462700764817, &	! size 10, ratio 2
		2.04981903293732, 1.20916697342055, 0.108478335713496, 0.912513517238510, &	! size 10, ratio 3
		2.40722168261355, 1.30811672639663, 5.825940715847113E-003, 1.22516460991257, &	! size 10, ratio 4
		2.37166743841577, 1.27322805612039, 9.146157931267608E-003, 1.13626397409937 /), &	! size 10, ratio 5
		(/ 4, 5, 10 /))
	REAL*8, PARAMETER :: MieLUT(4,5,10,NFastTUVWvl) = RESHAPE((/ MieLUTWvl1, MieLUTWvl2, MieLUTWvl3, &
	                     MieLUTWvl4, MieLUTWvl5, MieLUTWvl6, MieLUTWvl7 /), (/ 4, 5, 10, NFastTUVWvl /))

	IF (InParticle%numberofparticles .LE. 1.0e-6) THEN    !No particles, use dummy numbers
		InParticle%ExtEff = 1.0
		InParticle%ScaEff = 1.0
		InParticle%BackScaEff = 1.0
		InParticle%AssymParam = 1.0
                RETURN          
	END IF

	ALLOCATE (AqVolConc(HowManyAqChems+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRadFASTTUV(), couldn't allocate "// &
                               "AqVolConc at size "//	&
	     TRIM(INT2STR(HowManyAqChems+HowManyAqCations+HowManyAqAnions)))
	ALLOCATE (OrgVolConc(HowManyOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRadFASTTUV(), couldn't allocate OrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyOrgChems)))
	ALLOCATE (AqOrgVolConc(HowManyAqOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRadFASTTUV(), couldn't allocate AqOrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))

	CALL OpticalVolumeConcentrations (InParticle, AqVolConc, OrgVolConc, AqOrgVolConc, &
	                                  SolutionVolConc, ShellVolumeConc)

	DO ib = 1, NFastTUVWvl

	   !Shell and core indices and radii (the table uses only the radii)
	   CALL ShellAndCoreIndices (InParticle, ib, FastTUVWavelengths(ib), FastTUVWaveLow(ib), &
	                             AqVolConc, OrgVolConc, AqOrgVolConc, SolutionVolConc, &
	                             ShellVolumeConc, a_w, b_w, SumEr, SumEi, SolarShellInd, &
	                             SolarCoreInd, HasShell)

	   !Effective-radius bin
	   IF (.NOT. (InParticle%EffectiveRadius .GT. 0. .AND. &
	              InParticle%EffectiveRadius .LE. MieSizeEdge(10))) THEN
	      WRITE(*,*) "InParticle%EffectiveRadius = ",InParticle%EffectiveRadius
	      call error("Unknown InParticle%EffectiveRadius range. STOP")
	   END IF
	   sizebin = 1
	   DO WHILE (InParticle%EffectiveRadius .GT. MieSizeEdge(sizebin))
	      sizebin = sizebin + 1
	   END DO

	   !Core/shell radius-ratio bin
	   radratio = InParticle%AbsCoreRad/InParticle%EffectiveRadius
	   IF (.NOT. (radratio .GT. 1e-09)) THEN
	      WRITE(*,*) "InParticle%AbsCoreRad = ",InParticle%AbsCoreRad
	      WRITE(*,*) "InParticle%EffectiveRadius = ",InParticle%EffectiveRadius
	      WRITE(*,*) "radratio = ",radratio
	      call error("Unknown radratio range. STOP")
	   END IF
	   ratiobin = 1
	   DO WHILE (ratiobin .LT. 5 .AND. radratio .GT. MieRatioEdge(ratiobin))
	      ratiobin = ratiobin + 1
	   END DO

	   InParticle%ExtEff(ib)     = MieLUT(1,ratiobin,sizebin,ib)
	   InParticle%ScaEff(ib)     = MieLUT(2,ratiobin,sizebin,ib)
	   InParticle%BackScaEff(ib) = MieLUT(3,ratiobin,sizebin,ib)
	   !The table holds the asymmetry parameter times the scattering efficiency (as DMiLay
	   !returns it), so divide by the scattering efficiency here
	   InParticle%AssymParam(ib) = MieLUT(4,ratiobin,sizebin,ib)/InParticle%ScaEff(ib)

	END DO

	DEALLOCATE(AqVolConc, OrgVolConc, AqOrgVolConc, STAT = I)

	RETURN
END SUBROUTINE ShellRefIndAndRadFASTTUV


SUBROUTINE AerosolOptPropFASTTUV (ExtCoeff, SingScat, Assym, BackScatCoeff, &
                          SubExtCoeff, SubSingScat, Flag)!&
                          !ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                          !SubExtCoeffSep, SubSingScatSep)
! 08/30/2010 Changed for the photolysis wavebands MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
! 05/03/2012 Changed to calculate backscattering coefficient, units
!             for extcoeff and backscatcoeff now Mm^-1
! 08/16/2012 Fixed backscattering coefficient, added submicron extinction & SSA
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added flag to decide between mixing rules
!            Flag = 0, core-in-shell (default)
!            Flag = 1, external mixture of BC and rest
!            Flag = 2, volume average dielectric constant mixing rule
!            Flag = 3, Maxwell-Garnett mixing rule
!
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE

        !Input Variables
        INTEGER :: Flag
	!Output Variables
	REAL*8, DIMENSION(NFastTUVWvl) :: ExtCoeff, SingScat, Assym, BackScatCoeff
	REAL*8, DIMENSION(NFastTUVWvl) :: SubExtCoeff, SubSingScat
        !Below are for external mixture of BC and "shell"
!	REAL*8, DIMENSION(NFastTUVWvl) :: ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep
!	REAL*8, DIMENSION(NFastTUVWvl) :: SubExtCoeffSep, SubSingScatSep

	!Internal variables
	INTEGER :: ib, numbin
        LOGICAL :: FASTTUV
	REAL*8, DIMENSION(NFastTUVWvl) :: ScatCoeff, SubScatCoeff
!	REAL*8, DIMENSION(NFastTUVWvl) :: ScatCoeffSep, SubScatCoeffSep
	REAL*8 :: AssymSum,xxx, hbs_frac, assym_sq, bc_vol_frac

        ! CMB: Gfortran doesn't like the additional definitions of the 
        !      public variables, nor does it like "Particle" being used
        !      both as a type and a variable name
	Type(Particle),POINTER :: p_Particle


        DO ib = 1, NFastTUVWvl
		ExtCoeff(ib) = 0.0
		ScatCoeff(ib) = 0.0
		BackScatCoeff(ib) = 0.0
		SubExtCoeff(ib) = 0.0
		SubScatCoeff(ib) = 0.0
                AssymSum = 0.0
!		ExtCoeffSep(ib) = 0.0
!		ScatCoeffSep(ib) = 0.0
!		BackScatCoeffSep(ib) = 0.0
!		SubExtCoeffSep(ib) = 0.0
!		SubScatCoeffSep(ib) = 0.0
!                AssymSumSep = 0.0
		p_Particle => Particles%First

		DO WHILE (ASSOCIATED(p_Particle))

!                     p_Particle%Temperature = TEMP
                     ! WRITE(*,*)    'cur%Temperature = ',TEMP   
                     !Recalculate optical parameters
                     !WRITE(*,*) "Call Optical, Particle #", I
                     CALL ShellRefIndAndRadFASTTUV(p_Particle)

                     IF (p_Particle%EffectiveRadius .gt. 0) THEN 
                       bc_vol_frac = (p_Particle%AbsCoreRad/&
                               p_Particle%EffectiveRadius)**3
                     ELSE
                       bc_vol_frac = 0.0
                     ENDIF
		     xxx = p_Particle%NumberofParticles*&
                        (PI*p_Particle%EffectiveRadius**2)

                     IF (FLAG .EQ. 0) THEN !Core in Shell  
			ExtCoeff(ib) = ExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
			AssymSum = AssymSum + &
                            xxx*p_Particle%ScaEff(ib)*p_Particle%AssymParam(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
                        ENDIF

                     ELSE IF (Flag .EQ. 1) THEN !External mixture
                        ExtCoeff(ib) = ExtCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ExtEffBC(ib)&
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
			AssymSum = AssymSum + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib)*&
                                p_Particle%AssymParamBC(ib) &
                                       +(1.0-bc_vol_frac)*xxx* &
                                       p_Particle%ScaEffShell(ib)*&
                                       p_Particle%AssymParamShell(ib) 
                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN	
                          SubExtCoeff(ib) = SubExtCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ExtEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
                        ENDIF
	              ELSE IF (FLAG .EQ. 2) THEN !Volume-averaged internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffVA(ib)*p_Particle%AssymParamVA(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
                        ENDIF
                      ELSE IF (FLAG .EQ. 3) THEN !Maxwell-Garnett internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffMG(ib)*&
                            p_Particle%AssymParamMG(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
                        ENDIF
                      ENDIF

                      p_Particle => p_Particle%Next
		END DO

		!Formula from p.327, eqn 9.103 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		SingScat(ib) = ScatCoeff(ib)/ExtCoeff(ib)
		SubSingScat(ib) = SubScatCoeff(ib)/SubExtCoeff(ib)                
!		SingScatSep(ib) = ScatCoeffSep(ib)/ExtCoeffSep(ib)
!		SubSingScatSep(ib) = SubScatCoeffSep(ib)/SubExtCoeffSep(ib)   
		!IF(ib .EQ. 419 .OR. ib .EQ. 420) WRITE(*,*) "SingScat: ", SingScat(ib), ScatCoeff(ib), ExtCoeff(ib) 
		
		!Formula from p.329, eqn 9.117 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		Assym(ib) = AssymSum/ScatCoeff(ib)
!		AssymSep(ib) = AssymSumSep/ScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "Assym: ", Assym(1)

		!Convert ExtCoeff and BackScatCoeff from 1/cm to 1/Mm
		ExtCoeff(ib) = 1.0e8*ExtCoeff(ib)
		SubExtCoeff(ib) = 1.0e8*SubExtCoeff(ib)
!		ExtCoeffSep(ib) = 1.0e8*ExtCoeffSep(ib)
!		SubExtCoeffSep(ib) = 1.0e8*SubExtCoeffSep(ib)


		!hbs_frac = hemispherical backscattering fraction
                !See Wiscombe and Grams (1976) J. Atmos. Sci. 33(12), pp. 2440-2451.
                !Integrate Eq. 23 using Henyey-Greenstein phase function from Eq. 16
                assym_sq = assym(ib)**2
                hbs_frac = 0.5*((assym_sq-1)/assym(ib))*(1.0/SQRT(assym_sq+2*assym(ib)+1)-1.0/SQRT(assym_sq+1))
		BackScatCoeff(ib) = hbs_frac*ScatCoeff(ib)
                BackScatCoeff(ib) = 1.0e8*BackScatCoeff(ib)

!                assym_sq = assymSep(ib)**2
!                hbs_frac = 0.5*((assym_sq-1)/assymSep(ib))*(1.0/SQRT(assym_sq+2*assymSep(ib)+1)-1.0/SQRT(assym_sq+1))
!		BackScatCoeffSep(ib) = hbs_frac*ScatCoeffSep(ib)
!                BackScatCoeffSep(ib) = 1.0e8*BackScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "ExtCoeff: ", ExtCoeff(1)

	END DO
	
        
	RETURN
END SUBROUTINE AerosolOptPropFASTTUV
