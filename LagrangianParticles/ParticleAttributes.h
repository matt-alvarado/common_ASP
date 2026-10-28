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

!! This subroutine calculates the volume-average refractive indices for
!! the reflective shell and then calls DMILAY for the efficienies and 
!! assymetry parameter
SUBROUTINE ShellRefIndAndRad(InParticle)
! 08/30/2010 Changed for the Photolyis Wavebands MJA
! 02/27/2012 Changed to wavelength-dependent refractive indices
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added Maxwell Garnett mixing rule
	USE ModelParameters
	USE Chemistry, ONLY : AqMolecularMass, &
			        SolidSaltDensity, AqPhaseChemicalNames, &
				OrgMolecularMass, OrgDensity, &
				AqOrgMolecularMass, AqOrgDensity, &
                                CationIonicRefraction, &
                                AnionIonicRefraction,  &
                                RefIndFlag

	USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, ERROR
	
	IMPLICIT NONE
	
	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!!External Subroutines 
        EXTERNAL  DMiLay

	!! Internal Variables
	INTEGER, PARAMETER  :: MaxAng = 100
 	INTEGER	:: I, JJ, NumAng, ib, w_ind_low
	REAL*8	:: SingleParticleVolume, VolumeConc, Factor,xxx
	REAL*8	:: Er, Ei, SumEr, SumEi, Wavelength, Sum
	REAL*8	:: SolnMolWeight, MolalRefraction, MolVol, & 
                   TotalVolConc, Rshell, ShellVolumeConc
        REAL*8  :: ErBC, EiBC, ErVA, EiVA, bcvolrat, ErMG, EiMG

	REAL*8, ALLOCATABLE	:: AqVolConc(:), &
				OrgVolConc(:), AqOrgVolConc(:), MolFrac(:)

        !These are refractive indices given as data statements
        !!Wavelength-dependent refractive indices from version 3.1a
        !!of Optical Properties of Aerosols and Clouds (OPAC) 
        !! M. Hess, P. Koepke, and I. Schult (1998): 
        !!Optical Properties of Aerosols and clouds: 
        !!The software package OPAC, Bull. Am. Met. Soc., 79, 831-844.
        !Wavelength-dependent complex refractive indices (WAVE, SOOT_REAL,
        !SOOT_IMAG, WATER_REAL, WATER_IMAG, WASO_REAL, WASO_IMAG, INSO_REAL,
        !INSO_IMAG, SEASALT_REAL, SEASALT_IMAG, SULF_REAL, SULF_IMAG,
        !DUST_REAL, DUST_IMAG) are module-level PARAMETER arrays declared in
        !LagrangianAerosols.f90 (OPAC v3.1a data; see that file for the full
        !citation) -- no longer rebuilt here on every call.
        REAL*8 :: a_w, b_w, n_real, n_imag, V_molal, R_molal

	COMPLEX :: SolarShellInd, SolarCoreInd, IRShellInd, IRCoreInd, VAInd
        COMPLEX :: Eshell, EBC, Factor2, EMG, MGInd
	REAL*8	:: Angles(MaxAng), CosAng(MaxAng), M1( MaxAng, 2), &
			M2( MaxAng, 2 ), S21( MaxAng, 2 ), &
		D21( MaxAng, 2 ), SolutionRefracIndex(2), SolutionVolConc



	!! BUILD N ARRAY for AqSpecies Volume Concentrations
	ALLOCATE (AqVolConc(HowManyAqChems+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "AqVolConc at size "//	&
	     TRIM(INT2STR(HowManyAqChems+HowManyAqCations+HowManyAqAnions)))

	!! BUILD N ARRAY for Org Species Volume Concentrations
	ALLOCATE (OrgVolConc(HowManyOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate OrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyOrgChems)))

	!! BUILD N ARRAY for Aq Org Species Volume Concentrations
	ALLOCATE (AqOrgVolConc(HowManyAqOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate AqOrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))

	!! BUILD N ARRAY for H2O and Ion Mol Fractions
	ALLOCATE (MolFrac(1+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate MolFrac at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))
	
	IF (InParticle%numberofparticles .GT. 1.0e-6) THEN
		SingleParticleVolume = (4.0*Pi/3.0)*(InParticle%effectiveradius)**3 !cm3
		
		VolumeConc = SingleParticleVolume*InParticle%numberofparticles !cm3/cm3
		
		!Convert aqueous concentrations to cm3 solution/cm3 air
		
		!Water
		AqVolConc(1)  = InParticle%AqChems(1) * &
                                InParticle%NumberOfParticles*AqMolecularMass(1)/ &
                                InParticle%SolutionDensity
		SolutionVolConc = AqVolConc(1)
		TotalVolConc = AqVolConc(1)
				
		!Solid Salts
		DO I = 2, HowManyAqChems
			AqVolConc(I)  = InParticle%AqChems(I) * &
                                        InParticle%NumberOfParticles*AqMolecularMass(I)/ &
                                        SolidSaltDensity(I)
			TotalVolConc = TotalVolConc + AqVolConc(I)
		END DO

		!Inorganic Ions
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(I)  = InParticle%AqChems(I) * &
                                        InParticle%NumberOfParticles*AqMolecularMass(I)/ &
                                        InParticle%SolutionDensity
			TotalVolConc = TotalVolConc + AqVolConc(I)
			SolutionVolConc = SolutionVolConc + AqVolConc(I)	
		END DO

		!Convert org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyOrgChems
			OrgVolConc(I)  = InParticle%OrgChems(I) * &
                                         InParticle%NumberOfParticles*OrgMolecularMass(I)/ &
                                         OrgDensity(I)
			TotalVolConc = TotalVolConc + OrgVolConc(I)	
		END DO
	
		!Convert hydrophilic org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(I)  = InParticle%AqOrgChems(I) * &
                                           InParticle%NumberOfParticles*AqOrgMolecularMass(I) &
                                           /AqOrgDensity(I)
			TotalVolConc = TotalVolConc + AqOrgVolConc(I)	
		END DO

		!Calculate shell volume conc (Total - BC)
		ShellVolumeConc = VolumeConc - OrgVolConc(1)
			
	ELSE    !No particles, use dummy numbers
		InParticle%ExtEff = 1.0
		InParticle%ScaEff = 1.0
		InParticle%BackScaEff = 1.0
		InParticle%AssymParam = 1.0
                RETURN          
	END IF

        !Calculate shell and core refractive index for each wavelength bin
        DO ib = 1, 451 !For photolysis between 250 and 700 nm, every 1 nm

           Wavelength = (250+(ib-1)*1) !in nm

           !Index of point on array WAVE closest to, but less than,
           !the wavelength for interpolation
           w_ind_low = INT(ib/50)+1

           !Calculate weightings of wavelength-dependent refractive index points
           a_w = (WAVE(w_ind_low+1)-Wavelength) &
                  /(WAVE(w_ind_low+1)-WAVE(w_ind_low))
           b_w = 1-a_w

           !Calculate Core Refractive Index
           n_real = a_w*SOOT_REAL(w_ind_low)+b_w*SOOT_REAL(w_ind_low+1)
           n_imag = -1.0*(a_w*SOOT_IMAG(w_ind_low)+b_w*SOOT_IMAG(w_ind_low+1))
	   SolarCoreInd = CMPLX(n_real, -1.0*n_imag)
           !IF (ib .eq. 419 .or. ib .eq. 420) WRITE(*,*) Wavelength, SolarShellInd, SolarCoreInd

           !We assume the core is BC, the shell everything else
           !BUT we need to handle the cases where the particle is all BC (no shell)
           ! or has no BC (all core). All core handled here.
           IF (ShellVolumeConc .LE. 0 .OR. InParticle%Dry) THEN
                !Shell and core refractive indices equal.
                SolarShellInd = SolarCoreInd
                !Calculate Core radius
	        InParticle%AbsCoreRad = OrgVolConc(1)/InParticle%NumberofParticles
	        InParticle%AbsCoreRad = (0.75*InParticle%AbsCoreRad/Pi)**(1.0/3.0) !In cm
                !Set Shell radius to above, and core a little less than shell
                InParticle%EffectiveRadius = InParticle%AbsCoreRad
                InParticle%AbsCoreRad = 0.99*InParticle%AbsCoreRad
           ELSE !Some or all shell, all shell handled at bottom
	     !Calculate shell refractive index
	     !Use the volume average dielectric constant mixing rule
	     !See Jacobson, "Fund. of Atm. Modeling", 2nd ed, p. 310
 	     SumEr = 0.
	     SumEi = 0.
 
             !Water refractive index and molal Refraction
             n_real = a_w*WATER_REAL(w_ind_low)+b_w*WATER_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WATER_IMAG(w_ind_low)+b_w*WATER_IMAG(w_ind_low+1))
             V_molal = 18.0674 !Tang, 1997: density = 0.9971 g/cm3, MW = 18.015 g/mol, V = MW/dens
	     R_molal = V_molal*(n_real**2 - 1.0)/(n_real**2 + 2) 

             !Calculate average refractive index for ionic solution
	     !Following Tang et. al., JGR V.102, D19, pp. 23269-23275, 1997. 
             !and Tang, JGR V.102, D2, pp. 1883-1893, 1997.)
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
                !WRITE(*,*) MolFrac(I), CationIonicRefraction(I-1)
	     END DO
	     DO I = 2+HowManyAqCations, 1+HowManyAqCations+HowManyAqAnions
	   	MolalRefraction = MolalRefraction + MolFrac(I)*AnionIonicRefraction(I-HowManyAqCations-1)
                !WRITE(*,*) MolFrac(I), AnionIonicRefraction(I-HowManyAqCations-1)
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

             !Solid Salts
             DO I = 2, HowManyAqChems
                IF(RefIndFlag(I) .EQ. 0) THEN
                   n_real = a_w*WATER_REAL(w_ind_low)+b_w*SULF_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*WATER_IMAG(w_ind_low)+b_w*SULF_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 1) THEN
                   n_real = a_w*SULF_REAL(w_ind_low)+b_w*SULF_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*SULF_IMAG(w_ind_low)+b_w*SULF_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 2) THEN
                   n_real = a_w*SEASALT_REAL(w_ind_low)+b_w*SEASALT_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*SEASALT_IMAG(w_ind_low)+b_w*SEASALT_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 3) THEN
                   n_real = a_w*DUST_REAL(w_ind_low)+b_w*DUST_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*DUST_IMAG(w_ind_low)+b_w*DUST_IMAG(w_ind_low+1))
                ENDIF
	        Er = n_real**2 - n_imag**2
	        Ei = 2*n_real*n_imag
                SumEr = SumEr + (AqVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqVolConc(I)/ShellVolumeConc)*Ei
	     END DO
         
             !Insoluble organics, skipping BC core
	     n_real = a_w*WASO_REAL(w_ind_low)+b_w*WASO_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WASO_IMAG(w_ind_low)+b_w*WASO_IMAG(w_ind_low+1))
	     Er = n_real**2 - n_imag**2
	     Ei = 2*n_real*n_imag
             DO I = 2, HowManyOrgChems
		SumEr = SumEr + (OrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (OrgVolConc(I)/ShellVolumeConc)*Ei
	     END DO

             !Hydrophilic Organics
	     n_real = a_w*WASO_REAL(w_ind_low)+b_w*WASO_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WASO_IMAG(w_ind_low)+b_w*WASO_IMAG(w_ind_low+1))
	     Er = n_real**2 - n_imag**2
	     Ei = 2*n_real*n_imag
             DO I = 1, HowManyAqOrgChems
		SumEr = SumEr + (AqOrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqOrgVolConc(I)/ShellVolumeConc)*Ei
	     END DO

            !!For core-in-shell mixing rule
	    Factor = SQRT(SumEr**2 + SumEi**2)
	    InParticle%ShellRealRefracInd(ib) = SQRT((Factor + SumEr)/2.0)
	    InParticle%ShellImagRefracInd(ib) = SQRT((Factor - SumEr)/2.0)
	    SolarShellInd = CMPLX(InParticle%ShellRealRefracInd(ib),-1.0*InParticle%ShellImagRefracInd(ib))

            !We assume the core is BC, the shell everything else
            !BUT we need to handle the cases where the particle is all BC (all core, no shell)
            ! or has no BC (all shell, no core). No core, all shell handeled here.

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
	  CALL DMiLay( InParticle%AbsCoreRad, &
			InParticle%EffectiveRadius, &
			Wavelength, &
			SolarCoreInd, & 
			SolarCoreInd, &
			CosAng, &
			NumAng, &
			InParticle%ExtEffBC(ib), &
			InParticle%ScaEffBC(ib), &
			InParticle%BackScaEffBC(ib), &
			InParticle%AssymParamBC(ib), & 
                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParamBC(ib) = &
                InParticle%AssymParamBC(ib)/InParticle%ScaEffBC(ib)

          !Now call for a hypothetical particle of
          !same effective size but all shell
!          WRITE(*,*) "Call 3"
	  CALL DMiLay( InParticle%AbsCoreRad, &
			InParticle%EffectiveRadius, &
			Wavelength, &
			SolarShellInd, & 
			SolarShellInd, &
			CosAng, &
			NumAng, &
			InParticle%ExtEffShell(ib), &
			InParticle%ScaEffShell(ib), &
			InParticle%BackScaEffShell(ib), &
			InParticle%AssymParamShell(ib), & 
                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParamShell(ib) = &
                InParticle%AssymParamShell(ib)/InParticle%ScaEffShell(ib)

!          WRITE(*,*) "Call 4"
          !Now call for volume average dielectric constant mixing rule
	  CALL DMiLay( InParticle%AbsCoreRad, &
			InParticle%EffectiveRadius, &
			Wavelength, &
			VAInd, & 
			VAInd, &
			CosAng, &
			NumAng, &
			InParticle%ExtEffVA(ib), &
			InParticle%ScaEffVA(ib), &
			InParticle%BackScaEffVA(ib), &
			InParticle%AssymParamVA(ib), & 
                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParamVA(ib) = &
                InParticle%AssymParamVA(ib)/InParticle%ScaEffVA(ib)

          !Now call for Maxwell-Garnett mixing rule
	  CALL DMiLay( InParticle%AbsCoreRad, &
			InParticle%EffectiveRadius, &
			Wavelength, &
			MGInd, & 
			MGInd, &
			CosAng, &
			NumAng, &
			InParticle%ExtEffMG(ib), &
			InParticle%ScaEffMG(ib), &
			InParticle%BackScaEffMG(ib), &
			InParticle%AssymParamMG(ib), & 
                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParamMG(ib) = &
                InParticle%AssymParamMG(ib)/InParticle%ScaEffMG(ib)
          
	END DO

        DEALLOCATE(AqVolConc,OrgVolConc,AqOrgVolConc, MolFrac, STAT = I)

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
