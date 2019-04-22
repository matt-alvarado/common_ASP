!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ChemicalPropertyRoutines.h
!! This calculates the molecular diffusion coefficient of all gases.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							    !!
!!									    !!
!! Month  Year   Name              Description				    !!
!! 07     2006   Matt Alvarado     Began Update History	                    !!
!! 10/15  2010   Matt Alvarado     Removed Optional Arguments		    !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP       !!
!!                                 a one-box model or subroutine.           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		     !!
!! 1. FUNCTION MolecularDiffusionCoefficient (ChemIndex)                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIFFUSIVITY OF TRACE CHEMICAL SPECIES                                     !!
!! as a function of temperature and pressure.                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FUNCTION Mole (ChemIndex)

!        Mole%r = 0.375
!        Mole%i = 0.0

!	RETURN
!END FUNCTION Mole


FUNCTION MolecularDiffusionCoefficient (ChemIndex)

	USE ModelParameters,     ONLY : cm,AirMolecMass,Angstrom,Rstar, &
                                        Avogadro,Pi,ChemScale
	USE InfrastructuralCode, ONLY : Warn, Error, INT2STR
	USE GridPointFields,	 ONLY : GetAirDensity, GetTemp

	USE XPLEXIFY
        USE MYTYPE
	IMPLICIT NONE

	!! External Variables
	INTEGER		:: ChemIndex
        TYPE (XPLEX)    :: MolecularDiffusionCoefficient
	!CHARACTER*(*), OPTIONAL :: ChemName

	!! Internal Variables
	INTEGER :: i 
        TYPE (XPLEX) :: temp

        i = ChemIndex


        temp = Rstar*GetTemp()*AirMolecMass*ChemScale/	&
	Avogadro/GasMolecularMass(i)*(GasMolecularMass(i)+		&
	Avogadro*AirMolecMass)/2/Pi
        
        temp = sqrt(temp)
        temp =  0.375/((4.5*Angstrom)*(4.5*Angstrom))/GetAirDensity()*temp
        MolecularDiffusionCoefficient = temp

	RETURN
END FUNCTION MolecularDiffusionCoefficient
