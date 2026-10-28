      LOGICAL FUNCTION TSTBAD( FUNCNM, RELERR )

c        Called by MIEV0 self-test (TESTMI) when a computed value differs
c        from the stored reference.  Reports the discrepancy and returns
c        .FALSE. so that MIEV0 treats the self-test as failed.
c        (Patch 1.3: the routine was referenced by MIEV0.f but not present.)

      IMPLICIT NONE
      CHARACTER*(*) FUNCNM
      REAL          RELERR

      WRITE( *, '(/,3A,1P,E11.2,A)' ) ' MIEV0 self-test: ', FUNCNM,
     &       ' differs from reference; relative error =', RELERR
      TSTBAD = .FALSE.
      RETURN
      END
