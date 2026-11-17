!! ASP (c), 2004-2026, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! FEEReplay.f90
!!
!! Replay test for the electrolyte equilibrium solver.  Run it in a directory that holds
!! an ASP deck (InorganicInputDecks/, OutputFiles/) with the same chemistry as the run that
!! wrote the dump; the deck's solver settings select the solver that is replayed.
!!
!!   FEEReplay <dump file> [-out <file>] [-vdiag] [-repeat <n>]
!!
!!   -out <file>   write each replayed result (amounts and ReturnType) to <file>
!!   -vdiag        per-reaction diagnosis of the equilibrium violations in the results
!!   -repeat <n>   replay the dump n times (solver timing)
!!
!! Writes FEEReplay.txt (see ReplayFEE in FEEHarness.f90).  Dumps come from asp_feedump
!! (Tests/Makefile, Tests/README.md).

PROGRAM FEEReplay

	USE ModelParameters,      ONLY : SetDomainSize
	USE InfrastructuralCode,  ONLY : SetFileHandleCounter
	USE Chemistry,            ONLY : SetChemistryParams
	USE Aerosols,             ONLY : InitializeParticlesSectional
	USE Condensation,         ONLY : SetAllDissolution
	USE FEEHarness,           ONLY : ReplayFEE

	IMPLICIT NONE

	CHARACTER(LEN=512) :: DumpFile, ResultFile, Arg
	LOGICAL :: VDiag
	INTEGER :: NRep, I, NArgs

	DumpFile = "" ; ResultFile = "" ; VDiag = .FALSE. ; NRep = 1
	NArgs = COMMAND_ARGUMENT_COUNT ()
	I = 1
	DO WHILE (I .LE. NArgs)
		CALL GET_COMMAND_ARGUMENT (I, Arg)
		SELECT CASE (TRIM(Arg))
		CASE ("-out")
			I = I + 1 ; CALL GET_COMMAND_ARGUMENT (I, ResultFile)
		CASE ("-vdiag")
			VDiag = .TRUE.
		CASE ("-repeat")
			I = I + 1 ; CALL GET_COMMAND_ARGUMENT (I, Arg) ; READ (Arg, *) NRep
		CASE DEFAULT
			DumpFile = Arg
		END SELECT
		I = I + 1
	END DO
	IF (LEN_TRIM(DumpFile) .EQ. 0) STOP "usage: FEEReplay <dump file> [-out <file>] [-vdiag] [-repeat <n>]"

	!Model set-up as in the Lagrangian parcel run (chemistry, dissolution, particles)
	CALL SetDomainSize ()
	CALL SetFileHandleCounter (12)
	CALL ReadMainInputDeck ()
	CALL SetChemistryParams
	CALL SetAllDissolution
	CALL InitializeParticlesSectional

	CALL ReplayFEE (DumpFile, ResultFile, VDiag, NRep)
	WRITE (*,*) "FEEReplay done: see FEEReplay.txt"

END PROGRAM FEEReplay
