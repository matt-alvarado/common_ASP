/* Part of the FEE dump build (Tests/Makefile): attach the electrolyte-solver dump hook
   (FEEHarnessAttach in FEEHarness.f90) before the ASP main program starts, so that
   asp_feedump is the unchanged model plus the dump.  The hook only sets a pointer here;
   it reads its settings (ASP_FEE_DUMP...) at the first solver call. */
extern void fee_harness_attach(void);

__attribute__((constructor)) static void fee_dump_attach(void)
{
	fee_harness_attach();
}
