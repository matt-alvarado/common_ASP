!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Gibbs-DuhemEvaluationFunction.f90
!! Provides the needed integrad function for binary water activity
!! calculations.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES						        !!
!! 1. NONE									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. FUNCTION EvalIdLogGamma (IonicStr, Q, Stoich, Charge)         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A portion of the integrated form of the Gibbs-Duhem Equation is			!!
!! in integral from 0 to I of I d log_10 Gamma, which is integrated using	!!
!! DQXGS.  DQXGS requires a function that evaluates the integrand at a		!!
!! given value.  This is it.  This uses the same form of the equations as   !!
!! Resch, 1995.																!!
!! ----------------------------------------------------------------------   !!
!!    dlogG     q*B*I*(1+0.1 I)^(q-1)      0.5107(0.5*Sqrt[I]+0.069(C-1)I^4)!!
!!  I ----- = ------------------------  -  -----------------------------	!!
!!     dI     10ln10[1+B(1+0.1I)^q -B]            (1+C*Sqrt[I])^2			!!
!! ----------------------------------------------------------------------   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This can't be housed within the Module Aerosols and be made EXTERNAL as  !!
!! well.  So it's is floating on its own.									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION EvalIdLogGamma (IonicStr, Q, Stoich, Charge)

	!Patch N6: same equations as before, but (a) I**3. and I**4. (real exponents force pow())
	!are written as products, (b) EXP(-0.023*I**3), SQRT(I) and (1+0.1*I) are computed once
	!instead of once per Q term, (c) (1+0.1*I)**Q is computed once per Q term instead of twice,
	!and (d) the unused D = (1+0.1*I)**(Q-1) is no longer evaluated in the composite branch.
	IMPLICIT NONE

	REAL*8  :: IonicStr, Q(3), Stoich(4), Charge(4)		! Total Ionic Strength and the Kusik-Meissner Q Value for the I,J pair
	REAL*8  :: B,C,D,SQRTI, I3, I4, EXPT, ONEPL, PW, ARG, TERM
	INTEGER :: I

	SQRTI = SQRT(IonicStr)
	I3    = IonicStr*IonicStr*IonicStr
	I4    = I3*IonicStr
	EXPT  = EXP(-0.023*I3)
	ONEPL = 1. + 0.1*IonicStr

	IF ((Q(2) .EQ. 0.) .AND. (Q(3) .EQ. 0.) .AND. (Stoich(1) .EQ. 0.)) THEN
		B = 0.75-0.065*Q(1)
		C = 1. + 0.055*Q(1)*EXPT
		D = ONEPL**(Q(1)-1.)
		EvalIdLogGamma = Q(1)*B*IonicStr*D / (10.*LOG(10.)*(1.0+B*D*ONEPL-B)) &
						 - 0.5107*(0.5*SQRTI+0.069*I4*(C-1.))					   &
						 / (1.+C*SQRTI) / (1.+C*SQRTI)
	ELSE
		I = 1
		B = 0.75-0.065*Q(I)
		C = 1. + 0.055*Q(I)*EXPT
		PW  = ONEPL**Q(I)
		ARG = 1.+B*PW-B
		IF (ARG .LE. 0.0) THEN
			WRITE(*,*) "IonicStr", IonicStr
			WRITE(*,*) "B: ", B, "Q(1): ", Q(I), "Q(2): ", Q(2), "Q(3): ", Q(3)
			write(*,*) ARG
		END IF
		EvalIdLogGamma = Stoich(I) * Charge(I) *			        & !!-- Leading Factors
						  (LOG10(ARG) -	0.5107*SqrtI/(1.0+SqrtI*C))

		I = 2
		B = 0.75-0.065*Q(I)
		C = 1. + 0.055*Q(I)*EXPT
		ARG = 1.+B*ONEPL**Q(I)-B
		EvalIdLogGamma = EvalIdLogGamma + Stoich(I) * Charge(I) *		&
						 (LOG10(ARG) - 0.5107*SqrtI/(1.0+SqrtI*C))

		I = 3
		B = 0.75-0.065*Q(I)
		C = 1. + 0.055*Q(I)*EXPT
		ARG = 1.+B*ONEPL**Q(I)-B
		EvalIdLogGamma = EvalIdLogGamma - Stoich(I) * Charge(I) *		&
						 (LOG10(ARG) - 0.5107*SqrtI/(1.0+SqrtI*C))

		EvalIdLogGamma = EvalIdLogGamma / Charge(4) / Stoich(4)
	END IF

	RETURN
END FUNCTION EvalIdLogGamma
