!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! LagrangianAerosols.f90
!! This is the main source file for particle structure, properties,
!! and aqueous-phase inorganic thermodynamics. It also includes the 
!! particle sorting routines needed for coagulation, condensation, 
!! and the moving-center distribution.
!! It mainly links all the particle header files together.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 09/07  2006   Matt Alvarado     Added array BoundaryParticles             !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 03/01  2012   Matt Alvarado     Added SetInitialDens to                   !!
!!                                  SUBROUTINE InitializeParticlesSectional  !!
!! 07/19  2013   Matt Alvarado     Moved IsNaN to InfrastructuralCode        !!  
!! 10/28  2026   Claude Sonnet 5.5 Commit 8cb0991: Phase 1.1: hoist OPAC     !!
!!                                 tables to module PARAMETERs               !!
!! 10/28  2026   Claude Sonnet 5.5 Commit a56950d: N6-N11: integrand caches, !!
!!                                 EqCoeff cache, tabulated Gibbs-Duhem      !!
!!                                 integral (N11 changes results ~2e-4)      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES				!!
!! 1. ModelParameters			!!
!! 2. InfrastructuralCode		!!
!! 3. GridPointFields			!!
!! 4. Time				!!
!! 5. Chemistry				!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES					            !!
!! 1. ParticleStructure.h					    !!
!! 2. InitializationAndSamplingRoutines.h			    !!
!! 3. ParticleAttributes.h					    !!
!! 4. Thermodynamics.h						    !!
!! 5. SortAerosol.h						    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE InitializeParticles()
!! 2. SUBROUTINE ClearParticles ()
!! 3. SUBROUTINE SCInitializeWaterParticles()
!! 4. SUBROUTINE ReInitializeParticles()
!! 5. FUNCTION Molality (ChemIndex, InParticle)
!! 6. FUNCTION CationMolality (ChemIndex, InParticle)
!! 7. FUNCTION AnionMolality (ChemIndex, InParticle)
!! 8. SUBROUTINE ChangeMolality (MolalityChange, ChemIndex, InParticle)
!! 9. SUBROUTINE ChangeCationMolality (MolalityChange, ChemIndex, InParticle)
!!10. SUBROUTINE ChangeAnionMolality (MolalityChange, ChemIndex, InParticle)
!!11. FUNCTION CopyLagrangianParticle (InParticle)
!!12. SUBROUTINE InitializeParticlesSectional()
!!13. SUBROUTINE DissociateAllElectrolytes()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Aerosols

IMPLICIT NONE

PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Now specify which subroutines may be accessed by the outside world	 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PUBLIC	:: AerosolModeNames,			&
           AerosolModes,                        &
           AerosolPH,				&
           CurvatureCorrection,			&
           OrgCurvatureCorrection,		&
           DumpParticleContentsAtError,         &
           FindElectrolyteEquilibrium,	        &
           FindAqElectrolyteEquilibriumForGridPoint, &
           HypotheticalElectrolyteConcentrations, &
           InitializeParticlesSectional,        &
           KusikMeissner,			&
           Molality,				&
           ParticleKnudsen,			&
           Particle,				&
           Particles, ParticleArray,	        &
           ParticleKnudsenForEnergy,	        &   
           ParticleMass,			&
           ReynoldsNumber,			&
           ReformElectrolyte,			&
           RecalculateRadius,			&
           SortAerosolAtGridPointForCondensation, &
           SortAerosolAtGridPointForCoagulation, &
           TerminalVelocity,			&
           UpdateWaterActivity,                 &
           RegridAerosol,                       & 
           EquilibriumWaterContentAmount,       &
           ReadBoundaryDistributionsOrganic,    &
           BoundaryParticles,                   &
           SortBoundaryPart,                    &
           AerosolOpticalDepth,                 &
           ShellRefIndAndRad,                   &
           AerosolOptProp,                      &
           PopulateParticlesSectionsRightAway,  &
           ReadDistributionsOrganic,            &
           InputFlagRatioOrMass

INCLUDE "ParticleStructure.h"

!! The flag from the input file AerosolModes.in
INTEGER :: InputFlagRatioOrMass

!! The windspeed for the O'Dowd 1997 Salt Parameterization 
!! (0 means don't use it)
!REAL*8  :: ODowdSaltWindSpeed
!INTEGER :: ODowdSaltNumbBins

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The Aerosol Modes array holds, in order (first index is mode identifier): !!
!!   Nt        (the number concentration of the mode per cc)		     !!
!!   Dbar_g_N  (the mean of the mode)					     !!
!!   Rho_g	   (the standard deviation)				     !!
!!   C(1...N)  (fractional mass concentration of all aq. chemicals,          !!
!!             using aqueous indexing)	                                     !!
!!   Corg(1..L) (fractional mass concentration of all org. chemicals,        !!
!!              using org. indexing)                                         !!
!!   IC(1..M)  (Number of each type of insoluble core aerosol per input)     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8, ALLOCATABLE :: AerosolModes(:,:)
INTEGER	            :: HowManyAerosolModes
LOGICAL             :: InsolCoresYesOrNo
CHARACTER (len=64), ALLOCATABLE  :: AerosolModeNames(:)

!!These are for the environmental (boundary) aerosol distribution
REAL*8, ALLOCATABLE :: EnvAerosolModes(:,:)
INTEGER		    :: HowManyEnvAerosolModes
CHARACTER (len=64), ALLOCATABLE  :: EnvAerosolModeNames(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The Particles array is the linked list                   !!
!! of particles observable to the outside world.	    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE (ParticleArray) :: Particles

!!BoundaryParticles is the environmental particle list
TYPE (ParticleArray) :: BoundaryParticles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Wavelength-dependent complex refractive indices used by ShellRefIndAndRad !!
!! (optical-properties workflow). These are simulation-invariant constants   !!
!! from OPAC v3.1a and are declared once here (PARAMETER, so compile-time    !!
!! constant and implicitly thread-safe) instead of being rebuilt via array   !!
!! constructor on every call to ShellRefIndAndRad (i.e. once per particle).  !!
!! Patch 1.1 of the OpticalPropertiesModel optimization plan.                !!
!! Source: M. Hess, P. Koepke, and I. Schult (1998), Optical Properties of   !!
!! Aerosols and Clouds: The software package OPAC, Bull. Am. Met. Soc.,      !!
!! 79, 831-844. All values at 0% RH.                                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL*8, PARAMETER, DIMENSION(61) :: WAVE = 1.0e3*(/ &
        2.500e-01,  3.000e-01,  3.500e-01,  4.000e-01,  4.500e-01, & 
        5.000e-01,  5.500e-01,  6.000e-01,  6.500e-01,  7.000e-01, &
        7.500e-01,  8.000e-01,  9.000e-01,  1.000e+00,  1.250e+00, &
        1.500e+00,  1.750e+00,  2.000e+00,  2.500e+00,  3.000e+00, &
        3.200e+00,  3.390e+00,  3.500e+00,  3.750e+00,  4.000e+00, &
        4.500e+00,  5.000e+00,  5.500e+00,  6.000e+00,  6.200e+00, &
        6.500e+00,  7.200e+00,  7.900e+00,  8.200e+00,  8.500e+00, &
        8.700e+00,  9.000e+00,  9.200e+00,  9.500e+00,  9.800e+00, &
        1.000e+01,  1.060e+01,  1.100e+01,  1.150e+01,  1.250e+01, &
        1.300e+01,  1.400e+01,  1.480e+01,  1.500e+01,  1.640e+01, &
        1.720e+01,  1.800e+01,  1.850e+01,  2.000e+01,  2.130e+01, &
        2.250e+01,  2.500e+01,  2.790e+01,  3.000e+01,  3.500e+01, &
        4.000e+01/)

        REAL*8, PARAMETER, DIMENSION(61) :: SOOT_REAL = (/ & !From file soot00 for black carbon
        1.620e+00,  1.740e+00,  1.750e+00,  1.750e+00,  1.750e+00, & 
        1.750e+00,  1.750e+00,  1.750e+00,  1.750e+00,  1.750e+00, &
        1.750e+00,  1.750e+00,  1.750e+00,  1.760e+00,  1.760e+00, &
        1.770e+00,  1.790e+00,  1.800e+00,  1.820e+00,  1.840e+00, &
        1.860e+00,  1.870e+00,  1.880e+00,  1.900e+00,  1.920e+00, &
        1.940e+00,  1.970e+00,  1.990e+00,  2.020e+00,  2.030e+00, &
        2.040e+00,  2.060e+00,  2.120e+00,  2.130e+00,  2.150e+00, &
        2.160e+00,  2.170e+00,  2.180e+00,  2.190e+00,  2.200e+00, &
        2.210e+00,  2.220e+00,  2.230e+00,  2.240e+00,  2.270e+00, &
        2.280e+00,  2.310e+00,  2.330e+00,  2.330e+00,  2.360e+00, &
        2.380e+00,  2.400e+00,  2.410e+00,  2.450e+00,  2.460e+00, &
        2.480e+00,  2.510e+00,  2.540e+00,  2.570e+00,  2.630e+00, &
        2.690e+00/) 
        REAL*8, PARAMETER, DIMENSION(61) :: SOOT_IMAG = (/ &
       -4.500e-01, -4.700e-01, -4.650e-01, -4.600e-01, -4.550e-01, &
       -4.500e-01, -4.400e-01, -4.350e-01, -4.350e-01, -4.300e-01, &
       -4.300e-01, -4.300e-01, -4.350e-01, -4.400e-01, -4.500e-01, &
       -4.600e-01, -4.800e-01, -4.900e-01, -5.100e-01, -5.400e-01, &
       -5.400e-01, -5.495e-01, -5.600e-01, -5.700e-01, -5.800e-01, &
       -5.900e-01, -6.000e-01, -6.100e-01, -6.200e-01, -6.250e-01, &
       -6.300e-01, -6.500e-01, -6.700e-01, -6.800e-01, -6.900e-01, &
       -6.900e-01, -7.000e-01, -7.000e-01, -7.100e-01, -7.150e-01, &
       -7.200e-01, -7.300e-01, -7.300e-01, -7.400e-01, -7.500e-01, &
       -7.600e-01, -7.750e-01, -7.900e-01, -7.900e-01, -8.100e-01, &
       -8.200e-01, -8.250e-01, -8.300e-01, -8.500e-01, -8.600e-01, &
       -8.700e-01, -8.900e-01, -9.100e-01, -9.300e-01, -9.700e-01, &
       -1.000e+00/) 

        REAL*8, PARAMETER, DIMENSION(61) :: WATER_REAL = (/ & !Refractive index from cuma00 file
        1.362e+00,  1.349e+00,  1.343e+00,  1.339e+00,  1.337e+00, &
        1.335e+00,  1.333e+00,  1.332e+00,  1.331e+00,  1.331e+00, &
        1.330e+00,  1.329e+00,  1.328e+00,  1.327e+00,  1.323e+00, &
        1.321e+00,  1.313e+00,  1.306e+00,  1.261e+00,  1.371e+00, &
        1.478e+00,  1.423e+00,  1.400e+00,  1.369e+00,  1.351e+00, &
        1.332e+00,  1.325e+00,  1.298e+00,  1.265e+00,  1.363e+00, &
        1.339e+00,  1.312e+00,  1.294e+00,  1.286e+00,  1.278e+00, &
        1.272e+00,  1.262e+00,  1.255e+00,  1.243e+00,  1.229e+00, &
        1.218e+00,  1.179e+00,  1.153e+00,  1.126e+00,  1.123e+00, &
        1.146e+00,  1.210e+00,  1.258e+00,  1.270e+00,  1.346e+00, &
        1.386e+00,  1.423e+00,  1.443e+00,  1.480e+00,  1.491e+00, &
        1.506e+00,  1.531e+00,  1.549e+00,  1.551e+00,  1.532e+00, &
        1.519e+00/) 
        REAL*8, PARAMETER, DIMENSION(61) :: WATER_IMAG = (/ &
       -3.350e-08, -1.600e-08, -6.500e-09, -1.860e-09, -1.020e-09, & 
       -1.000e-09, -1.960e-09, -1.090e-08, -1.640e-08, -3.350e-08, &
       -1.560e-07, -1.250e-07, -4.860e-07, -2.890e-06, -8.700e-06, &
       -2.000e-04, -1.000e-04, -1.100e-03, -1.740e-03, -2.720e-01, &
       -9.240e-02, -2.314e-02, -9.400e-03, -3.500e-03, -4.600e-03, &
       -1.340e-02, -1.240e-02, -1.160e-02, -1.070e-01, -8.800e-02, &
       -3.920e-02, -3.210e-02, -3.390e-02, -3.510e-02, -3.670e-02, &
       -3.790e-02, -3.990e-02, -4.150e-02, -4.440e-02, -4.790e-02, &
       -5.080e-02, -6.740e-02, -9.680e-02, -1.420e-01, -2.590e-01, &
       -3.050e-01, -3.700e-01, -3.960e-01, -4.020e-01, -4.270e-01, &
       -4.290e-01, -4.260e-01, -4.210e-01, -3.930e-01, -3.790e-01, &
       -3.700e-01, -3.560e-01, -3.390e-01, -3.280e-01, -3.360e-01, &
       -3.850e-01/) 

        !From file waso00, water-soluble sulfates, nitrates, and organic
        REAL*8, PARAMETER, DIMENSION(61) :: WASO_REAL = (/ & 
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.520e+00,  1.520e+00,  1.520e+00,  1.510e+00, &
        1.510e+00,  1.470e+00,  1.420e+00,  1.420e+00,  1.420e+00, &
        1.430e+00,  1.430e+00,  1.450e+00,  1.452e+00,  1.455e+00, &
        1.460e+00,  1.450e+00,  1.440e+00,  1.410e+00,  1.430e+00, &
        1.460e+00,  1.400e+00,  1.200e+00,  1.010e+00,  1.300e+00, &
        2.400e+00,  2.560e+00,  2.200e+00,  1.950e+00,  1.870e+00, &
        1.820e+00,  1.760e+00,  1.720e+00,  1.670e+00,  1.620e+00, &
        1.620e+00,  1.560e+00,  1.440e+00,  1.420e+00,  1.750e+00, &
        2.080e+00,  1.980e+00,  1.850e+00,  2.120e+00,  2.060e+00, &
        2.000e+00,  1.880e+00,  1.840e+00,  1.820e+00,  1.920e+00, &
        1.860e+00/) 
        REAL*8, PARAMETER, DIMENSION(61) :: WASO_IMAG = (/ &
       -3.000e-02, -8.000e-03, -5.000e-03, -5.000e-03, -5.000e-03, &
       -5.000e-03, -6.000e-03, -6.000e-03, -7.000e-03, -7.000e-03, &
       -8.500e-03, -1.000e-02, -1.300e-02, -1.550e-02, -1.900e-02, &
       -2.250e-02, -1.750e-02, -8.000e-03, -1.200e-02, -2.200e-02, &
       -8.000e-03, -7.050e-03, -5.000e-03, -4.000e-03, -5.000e-03, &
       -1.300e-02, -1.200e-02, -1.800e-02, -2.300e-02, -2.700e-02, &
       -3.300e-02, -7.000e-02, -6.500e-02, -1.000e-01, -2.150e-01, &
       -2.900e-01, -3.700e-01, -4.200e-01, -1.600e-01, -9.500e-02, &
       -9.000e-02, -7.000e-02, -5.000e-02, -4.700e-02, -5.300e-02, &
       -5.500e-02, -7.300e-02, -1.000e-01, -2.000e-01, -1.600e-01, &
       -2.420e-01, -1.800e-01, -1.700e-01, -2.200e-01, -2.300e-01, &
       -2.400e-01, -2.800e-01, -2.900e-01, -3.000e-01, -4.000e-01, &
       -5.000e-01/) 

       !From file ssam00, for sea salt
       REAL*8, PARAMETER, DIMENSION(61) :: SEASALT_REAL = (/ &
       1.510e+00,  1.510e+00,  1.510e+00,  1.500e+00,  1.500e+00, &
       1.500e+00,  1.500e+00,  1.490e+00,  1.490e+00,  1.490e+00, &
       1.490e+00,  1.480e+00,  1.480e+00,  1.470e+00,  1.470e+00, &
       1.460e+00,  1.450e+00,  1.450e+00,  1.430e+00,  1.610e+00, &
       1.490e+00,  1.480e+00,  1.480e+00,  1.470e+00,  1.480e+00, &
       1.490e+00,  1.470e+00,  1.420e+00,  1.410e+00,  1.600e+00, &
       1.460e+00,  1.420e+00,  1.400e+00,  1.420e+00,  1.480e+00, &
       1.600e+00,  1.650e+00,  1.610e+00,  1.580e+00,  1.560e+00, &
       1.540e+00,  1.500e+00,  1.480e+00,  1.480e+00,  1.420e+00, &
       1.410e+00,  1.410e+00,  1.430e+00,  1.450e+00,  1.560e+00, &
       1.740e+00,  1.780e+00,  1.770e+00,  1.760e+00,  1.760e+00, &
       1.760e+00,  1.760e+00,  1.770e+00,  1.770e+00,  1.760e+00, &
       1.740e+00/) 
       REAL*8, PARAMETER, DIMENSION(61) :: SEASALT_IMAG = (/ &
      -5.000e-06, -2.000e-06, -3.240e-07, -3.000e-08, -2.430e-08, &
      -1.550e-08, -1.000e-08, -1.600e-08, -4.240e-08, -2.000e-07, &
      -1.080e-06, -1.950e-06, -4.240e-05, -1.410e-04, -3.580e-04, &
      -5.700e-04, -7.620e-04, -1.000e-03, -4.000e-03, -1.000e-02, &
      -3.000e-03, -2.050e-03, -1.600e-03, -1.400e-03, -1.400e-03, &
      -1.400e-03, -2.500e-03, -3.600e-03, -1.100e-02, -2.200e-02, &
      -5.000e-03, -7.000e-03, -1.300e-02, -2.000e-02, -2.600e-02, &
      -3.000e-02, -2.800e-02, -2.620e-02, -1.800e-02, -1.600e-02, &
      -1.500e-02, -1.400e-02, -1.400e-02, -1.400e-02, -1.600e-02, &
      -1.800e-02, -2.300e-02, -3.000e-02, -3.500e-02, -9.000e-02, &
      -1.200e-01, -1.300e-01, -1.350e-01, -1.520e-01, -1.650e-01, &
      -1.800e-01, -2.050e-01, -2.750e-01, -3.000e-01, -5.000e-01, &
      -1.000e+00/) 

       !From file suso00, for acidic (stratospheric) H2SO4 droplets
       REAL*8, PARAMETER, DIMENSION(61) :: SULF_REAL = (/ &
        1.484e+00,  1.469e+00,  1.452e+00,  1.440e+00,  1.432e+00, &
        1.431e+00,  1.430e+00,  1.429e+00,  1.429e+00,  1.428e+00, &
        1.427e+00,  1.426e+00,  1.425e+00,  1.422e+00,  1.413e+00, &
        1.403e+00,  1.394e+00,  1.384e+00,  1.344e+00,  1.293e+00, &
        1.311e+00,  1.350e+00,  1.376e+00,  1.396e+00,  1.385e+00, &
        1.385e+00,  1.360e+00,  1.337e+00,  1.425e+00,  1.424e+00, &
        1.370e+00,  1.210e+00,  1.140e+00,  1.200e+00,  1.370e+00, &
        1.530e+00,  1.650e+00,  1.600e+00,  1.670e+00,  1.910e+00, &
        1.890e+00,  1.720e+00,  1.670e+00,  1.890e+00,  1.740e+00, &
        1.690e+00,  1.640e+00,  1.610e+00,  1.590e+00,  1.520e+00, &
        1.724e+00,  1.950e+00,  1.927e+00,  1.823e+00,  1.780e+00, &
        1.870e+00,  1.930e+00,  1.920e+00,  1.920e+00,  1.900e+00, &
        1.890e+00/) 
       REAL*8, PARAMETER, DIMENSION(61) :: SULF_IMAG = (/ &
       -1.000e-08, -1.000e-08, -1.000e-08, -1.000e-08, -1.000e-08, &
       -1.000e-08, -1.000e-08, -1.470e-08, -1.670e-08, -2.050e-08, &
       -7.170e-08, -8.630e-08, -2.550e-07, -1.530e-06, -6.940e-06, &
       -1.200e-04, -4.160e-04, -1.260e-03, -3.760e-03, -9.550e-02, &
       -1.350e-01, -1.578e-01, -1.580e-01, -1.310e-01, -1.260e-01, &
       -1.200e-01, -1.210e-01, -1.830e-01, -1.950e-01, -1.650e-01, &
       -1.280e-01, -1.760e-01, -4.880e-01, -6.450e-01, -7.550e-01, &
       -7.720e-01, -6.330e-01, -5.860e-01, -7.500e-01, -6.800e-01, &
       -4.550e-01, -3.400e-01, -4.850e-01, -3.740e-01, -1.980e-01, &
       -1.950e-01, -1.950e-01, -2.050e-01, -2.110e-01, -4.140e-01, &
       -5.900e-01, -4.100e-02, -3.025e-02, -2.352e-02, -2.925e-01, &
       -3.150e-02, -2.000e-01, -1.800e-01, -1.800e-01, -1.900e-01, &
       -2.200e-01/) 

       !From file miam00, for mineral dust (Quartz and clay particles)
       REAL*8, PARAMETER, DIMENSION(61) :: DUST_REAL = (/ &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.520e+00,  1.520e+00, &
        1.510e+00,  1.510e+00,  1.510e+00,  1.500e+00,  1.500e+00, &
        1.500e+00,  1.480e+00,  1.460e+00,  1.440e+00,  1.430e+00, &
        1.420e+00,  1.460e+00,  1.220e+00,  1.120e+00,  1.060e+00, &
        1.190e+00,  1.850e+00,  2.220e+00,  2.940e+00,  2.910e+00, &
        2.570e+00,  1.910e+00,  1.830e+00,  1.810e+00,  1.740e+00, &
        2.000e+00,  1.630e+00,  1.540e+00,  1.510e+00,  1.470e+00, &
        1.490e+00,  1.770e+00,  2.050e+00,  2.200e+00,  2.390e+00, &
        2.690e+00,  2.990e+00,  2.570e+00,  2.420e+00,  2.420e+00, &
        2.340e+00/) 
       REAL*8, PARAMETER, DIMENSION(61) :: DUST_IMAG = (/ &
       -3.000e-02, -2.500e-02, -1.700e-02, -1.300e-02, -8.500e-03, &
       -7.800e-03, -5.500e-03, -4.500e-03, -4.500e-03, -4.000e-03, &
       -4.000e-03, -4.000e-03, -4.000e-03, -4.000e-03, -5.000e-03, &
       -5.700e-03, -6.400e-03, -7.600e-03, -1.400e-02, -3.900e-02, &
       -2.400e-02, -1.925e-02, -1.800e-02, -1.200e-02, -6.700e-03, &
       -8.700e-03, -1.800e-02, -3.600e-02, -5.500e-02, -6.300e-02, &
       -5.200e-02, -1.300e-01, -8.900e-02, -1.200e-01, -2.100e-01, &
       -2.900e-01, -4.400e-01, -5.400e-01, -6.500e-01, -6.500e-01, &
       -5.000e-01, -2.500e-01, -2.000e-01, -3.500e-01, -5.000e-01, &
       -3.500e-01, -2.200e-01, -2.400e-01, -2.600e-01, -3.200e-01, &
       -3.700e-01, -4.700e-01, -5.700e-01, -8.200e-01, -9.400e-01, &
       -1.000e+00, -8.000e-01, -7.800e-01, -6.700e-01, -6.200e-01, &
       -7.000e-01/) 

!

	!Patch N7/N11 caches (see Thermodynamics.h).  NOT thread-safe: fill eagerly or make per-thread before OpenMP.
	REAL*8, ALLOCATABLE, SAVE :: EqCoeffCacheT(:), EqCoeffCacheV(:)
	!Patch C1: per-reaction name tests done once (string compares were ~6% of CondTest)
	LOGICAL, ALLOCATABLE, SAVE :: EqRxnAnionIsCO3(:), EqRxnIsLevitocite(:), EqRxnIsAmmSulfate(:)
	INTEGER, PARAMETER :: GDTabN = 256
	REAL*8, ALLOCATABLE, SAVE :: GDTabF(:,:), GDTabG(:,:), GDTabKey(:,:), GDTabSHi(:)
	LOGICAL, ALLOCATABLE, SAVE :: GDTabOK(:), GDTabBad(:)

CONTAINS

INCLUDE "InitializationAndSamplingRoutines.h"
INCLUDE "ParticleAttributes.h"
INCLUDE "Thermodynamics.h"
INCLUDE "SortAerosol.h"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Erase all particles so we can start again !!
SUBROUTINE ClearParticles ()
	
        ! CB: what happens if I comment this out?
	!USE Chemistry,       ONLY : HowManyAqChems, &
        !                            HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : HowManyBins
	USE InfrastructuralCode, ONLY : TRANSCRIPT, INT2STR

	IMPLICIT NONE

	INTEGER :: L
	TYPE(Particle), POINTER :: Leader, Current

	Current => Particles%First
	Leader  => Current%Next
	Particles%First => Null()

10	CALL TRANSCRIPT("Clearing particle "//TRIM(INT2STR(Current%ParticleID)))

	DEALLOCATE(Current%GammaMixed, stat=L)
	DEALLOCATE(Current%AqChems, stat=L)
	DEALLOCATE(Current%OrgChems, stat=L)
	DEALLOCATE(Current, stat=L)


	IF (ASSOCIATED(Leader)) THEN
		Current => Leader
		Leader  => Leader%Next
		GOTO 10
	END IF

	RETURN
END SUBROUTINE ClearParticles


!!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC!!
!! Make particles that fit a pre-specified set of sizes, all dry. !!
!!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC!!
SUBROUTINE SCInitializeWaterParticles()

	USE InfrastructuralCode, ONLY : ERROR
	USE ModelParameters,     ONLY : micron

	IMPLICIT NONE

	!! Local Variables
	INTEGER :: Allocation_Error,I,J,K,HowManyRadii
	REAL*8  :: Radius(1000)
	TYPE(Particle), POINTER :: CurrentParticle



	Particles%First => Null()
	
	HowManyRadii = 2000

        DO I = 1, HowManyRadii
	    Radius(I) = 2.*micron/2.
        END DO


	CurrentParticle        => SCTailorMakeWaterParticle(Radius(1))
	Particles%First => CurrentParticle

	DO I = 2,HowManyRadii

		CurrentParticle%Next => SCTailorMakeWaterParticle(Radius(I))
		CurrentParticle      => CurrentParticle%Next

	END DO

	RETURN
END SUBROUTINE SCInitializeWaterParticles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Molality is the unit of choice for calculating	  !!
!! anything that goes on in the aqueous phase.  But	  !!
!! the amount of each chemical in each particle is stored !!
!! in moles, to make conservation more exact, and so this !!
!! function converts moles to Molality.			  !!
!!							  !!
!! Molality is MOLES OF SOLUTE / KG of SOLVENT		  !!
!!							  !!
!! BUT THIS ISN'T SCALED BY MOLES, SO BE CAREFUL!         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! CMB (AER, Inc): None of these molality functions initialize the value of the 
!                 returned function, nor do they handle floating-point equality.
!                 Revisions made...
REAL*8 FUNCTION Molality (ChemIndex, InParticle)

	USE ModelParameters, ONLY : moles, grams,WaterMolecMass

	IMPLICIT NONE

	INTEGER		   :: ChemIndex
	TYPE(Particle), INTENT(IN) :: InParticle	!Patch N10

    ! CMB: mods for floating-point equality and function initialization
    Molality = 0.0
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	Molality = InParticle%AqChems(ChemIndex) / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

END FUNCTION Molality

!! -- Ion molalities are calculated in a pass-along format
REAL*8 FUNCTION CationMolality (ChemIndex, InParticle)

	USE ModelParameters, ONLY : HowManyAqChems

	IMPLICIT NONE

	INTEGER		   :: ChemIndex
	TYPE(Particle), INTENT(IN) :: InParticle	!Patch N10

    ! CMB: mods for floating-point equality and function initialization
    CationMolality = 0.0
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	CationMolality = Molality (ChemIndex + HowManyAqChems, InParticle)

END FUNCTION CationMolality

!! -- Ion molalities are calculated in a pass-along format
REAL*8 FUNCTION AnionMolality (ChemIndex, InParticle)

	USE ModelParameters, ONLY : HowManyAqChems, HowManyAqCations

	IMPLICIT NONE

	INTEGER		   :: ChemIndex
	TYPE(Particle), INTENT(IN) :: InParticle	!Patch N10

    ! CMB: mods for floating-point equality and function initialization
    AnionMolality = 0.0
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	AnionMolality = Molality (ChemIndex + HowManyAqChems + HowManyAqCations, InParticle)

END FUNCTION AnionMolality

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This function adds amount MolalityChange to the molality	   !!
!! for a certain chemical in a particular particle.  Chemicals !!
!! are stored in terms of moles, so conversion is necessary,   !!
!! hence this subroutine.									   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ChangeMolality (MolalityChange, ChemIndex, InParticle)

	USE ModelParameters, ONLY : moles, WaterMolecMass

	IMPLICIT NONE

	REAL*8		   :: MolalityChange
	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	!! Change the molality by the requested amount.
	InParticle%AqChems(ChemIndex) = InParticle%AqChems(ChemIndex) + MolalityChange / 1000. *	&
									(InParticle%AqChems(1) / moles * WaterMolecMass)

END SUBROUTINE  ChangeMolality

!! -- Changes in Ion molalities are calculated in a pass-along format
SUBROUTINE ChangeCationMolality (MolalityChange, ChemIndex, InParticle)

	USE ModelParameters, ONLY : HowManyAqChems

	IMPLICIT NONE

	REAL*8		   :: MolalityChange
	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

!print *, "Before Cation Molality:", MOLALITY(ChemIndex+ HowManyAqChems,InParticle)," + ",MolalityChange

	CALL ChangeMolality (MolalityChange, ChemIndex + HowManyAqChems, InParticle)

!print *, "After Cation Molality:", MOLALITY(ChemIndex+ HowManyAqChems,InParticle)


END SUBROUTINE ChangeCationMolality

!! -- Changes in Ion molalities are calculated in a pass-along format
SUBROUTINE ChangeAnionMolality (MolalityChange, ChemIndex, InParticle)

	USE ModelParameters, ONLY : HowManyAqChems, HowManyAqCations

	IMPLICIT NONE

	REAL*8		   :: MolalityChange
	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	CALL ChangeMolality (MolalityChange, ChemIndex + HowManyAqChems + HowManyAqCations, InParticle)

END SUBROUTINE ChangeAnionMolality 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Copy a Lagrangian Particle Exactly !!
FUNCTION CopyLagrangianParticle (InParticle)

	USE Chemistry,			 ONLY : HowManyAqChems,		&
						HowManyAqCations,	&
						HowManyAqAnions,	&
						AqMolecularMass,	&
						HowManyAqEqReactions

	USE InfrastructuralCode, ONLY : ERROR
	IMPLICIT NONE

	TYPE(Particle), POINTER :: CopyLagrangianParticle
	TYPE(Particle) :: InParticle

	INTEGER :: I, allocation_error

	ALLOCATE (CopyLagrangianParticle, stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle Failed in CopyLagrangianParticle()")
	ALLOCATE (CopyLagrangianParticle%AqChems(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%AqChems Failed in CopyLagrangianParticle()")
	ALLOCATE (CopyLagrangianParticle%GammaMixed(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%GammaMixed Failed in CopyLagrangianParticle()")

	CopyLagrangianParticle%ParticleID = InParticle%ParticleID
	CopyLagrangianParticle%ParticleDistribution = InParticle%ParticleDistribution

	CopyLagrangianParticle%Sectional = InParticle%Sectional
	CopyLagrangianParticle%NumberOfParticles = InParticle%NumberOfParticles
	CopyLagrangianParticle%Edges(1) = InParticle%Edges(1)
	CopyLagrangianParticle%Edges(2) = InParticle%Edges(2)

	CopyLagrangianParticle%velocity(1) = InParticle%velocity(1)
	CopyLagrangianParticle%velocity(2) = InParticle%velocity(2)
	CopyLagrangianParticle%velocity(3) = InParticle%velocity(3)

	CopyLagrangianParticle%Temperature = InParticle%Temperature
	CopyLagrangianParticle%EmbryoRadius = InParticle%EmbryoRadius
	CopyLagrangianParticle%InsolubleRadius = InParticle%InsolubleRadius
	CopyLagrangianParticle%EffectiveRadius = InParticle%EffectiveRadius

	CopyLagrangianParticle%IonicStr = InParticle%IonicStr
	CopyLagrangianParticle%WaterActivity = InParticle%WaterActivity
	CopyLagrangianParticle%SolutionDensity = InParticle%SolutionDensity
	CopyLagrangianParticle%ParticleDensity = InParticle%ParticleDensity
	CopyLagrangianParticle%SurfaceTension = InParticle%SurfaceTension

	DO I = 1, HowManyAqEqReactions
	CopyLagrangianParticle%GammaMixed(I) = InParticle%GammaMixed(I)
	END DO

	DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		CopyLagrangianParticle%AqChems(I) = InParticle%AqChems(I)
	END DO
	CopyLagrangianParticle%Next => Null()

	RETURN
END FUNCTION CopyLagrangianParticle

SUBROUTINE InitializeParticlesSectional()

	USE InfrastructuralCode, ONLY : ERROR
        USE ModelParameters,     ONLY : SetInitialDens

	IMPLICIT NONE

	!! Local Variables
	INTEGER :: Allocation_Error

	Particles%First => Null()
        CALL SetInitialDens(1.5000) !MJA 03-01-2012
        !Moved here so that the initial density can be set and reset at will

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Open the Particle Modes Input Deck and Read It In !!
	!! And then use it to populate the AerosolModes Array!!
	CALL ReadDistributionsOrganic() !!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Fill the Array with a Selection of Particles, !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!This puts all particles into sections
	Call PopulateParticlesSectionsRightAway(.FALSE.)
	RETURN
END SUBROUTINE InitializeParticlesSectional

SUBROUTINE DissociateAllElectrolytes()
	
	USE GridPointFields, ONLY : GetTemp,		        &
				    GetGridCellChemBurden,	&
				    ReplaceGridCellChemBurden,	&
				    GetRelativeHumidity

	USE Chemistry,	 ONLY : HowManyEvolveGasChems, &
				HowManyAqEqReactions, &
				AqEquilibriaList, &
				HowManyAqChems, &
				HowManyAqCations, &
				HowManyAqAnions

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterations,	&
				    AqThermoEquilibriumError, &
				    ThermoBulkMode, &
				    WaterContentPrecision

	USE Time, ONLY : CurrentTime, BeginningTime

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR

	IMPLICIT NONE

	!! Internal Variables
	REAL*8 :: II
	INTEGER :: I, J, K, L, M, Q, HowManyLinks, ReturnType, OrgRxnIndex, AqRxnIndex, StoreReturn
	TYPE(Particle),POINTER :: Current	
	
!!Matt Alvarado added this to dissociate all solid electrolytes at 
!!initialization of the particles, creating a metastable aerosol initially 
!!that will then dry to equilibrium. 

	Current => Particles%First
	IF (ASSOCIATED(Current)) THEN
20		DO I = 1, HowManyAqEqReactions
			IF (AqEquilibriaList(I,6) .EQ. 1.0 ) THEN 
                        !!Is a solid Eq. Rxn. 
		
			!! The change amount, two cases: 
                        !! equations with and without water
				IF (AqEquilibriaList(I,8) .EQ. 0) THEN
				  II = Current%AqChems(INT(AqEquilibriaList(I,1)))
				ELSE
				  II = MIN(Current%AqChems(INT(AqEquilibriaList(I,1))),	& ! water limited, or
					Current%AqChems(1)/AqEquilibriaList(I,8))	  ! electrolyte limited?
				END IF		
				
				!! Add to ions:
				!! CATION:
				Current%AqChems(HowManyAqChems+INT(AqEquilibriaList(I,2)))	=	&
				Current%AqChems(HowManyAqChems+INT(AqEquilibriaList(I,2))) +	&				
				II * AqEquilibriaList(I,4)

				!! ANION:
				Current%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(I,3))) =	&
				Current%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(I,3))) +	&
				II * AqEquilibriaList(I,5)

				IF(AqEquilibriaList(I,22) .NE. 0) THEN
					!!Levitocite case, extra anion
					Current%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(I,22))) =	&
					Current%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(I,22))) +	&
					II * AqEquilibriaList(I,23)
				END IF
				

				! And Adjust WATER CONTENT appropriately
				Current%AqChems(1) = Current%AqChems(1) - II * AqEquilibriaList(I,8)

				!! -- Then the Electrolyte
				!IF (AqEquilibriaList(I,1) .GT. 0.) THEN					
				Current%AqChems(INT(AqEquilibriaList(I,1)))	= &
				Current%AqChems(INT(AqEquilibriaList(I,1))) - II

			END IF
		END DO !!Solid Dissociation loop
		
		!Set the particle to WET initially, 
                !so dissolution reactions will happen
		Current%Dry = .FALSE. 
	
		IF (ASSOCIATED(Current%Next)) THEN
			Current => Current%Next
			GOTO 20
		END IF
	END IF

END SUBROUTINE DissociateAllElectrolytes

END MODULE Aerosols
