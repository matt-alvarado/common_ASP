# Electrolyte-solver replay test

Test programs for `FindElectrolyteEquilibrium` (`LagrangianParticles/Thermodynamics.h`). They are not part of
the model build; the model only carries the hook `FEEStateHook` (`LagrangianAerosols.f90`), which it never sets.

| File | Purpose |
|---|---|
| `FEEHarness.f90` | dump writer and reader, replay, equilibrium checks and the per-reaction violation diagnosis |
| `FEEReplay.f90` | replay program |
| `FEEDumpAttach.c` | attaches the dump writer to the model at start-up (for `asp_feedump`) |
| `Makefile` | builds `asp_feedump` and `FEEReplay` from the objects of a model build |

## Build

Build the model with the ASP makefile first (ASP repository, `build/ASPMakefile.linux.gfortran`), then

```
cd Tests
make ASP_BUILD=/path/to/ASP/build      # the directory with ObjectFiles/ and ModuleFiles/
```

## Dump

`asp_feedump` is the ASP model with the dump attached. Run it like `asp` in a run directory; with no variable set
it behaves exactly as the model.

```
ASP_FEE_DUMP=fee.bin ASP_FEE_DUMP_EVERY=100 /path/to/asp_feedump
```

- `ASP_FEE_DUMP=<file>`: write sampled solver calls to `<file>` (each record: the full particle state the solver
  reads, grid temperature, pressure and water vapour, then the state it returns).
- `ASP_FEE_DUMP_EVERY=<n>`: every n-th call (default 100). Calls of the first equilibration are not sampled.
- `ASP_FEE_DUMP_ALL=1`: also sample calls that update the activity coefficients (the ThermoTest scans need it).

## Replay

Run in a directory with a deck of the same chemistry (for example the same deck in Lagrangian parcel mode).
The deck's solver settings (`AqOuterSolver`, `AqInnerSolver`, ...) choose the solver that is replayed.

```
/path/to/FEEReplay fee.bin [-out results.bin] [-vdiag] [-repeat <n>]
```

`FEEReplay.txt` lists the number of records, how many results are bitwise identical to the dumped ones and have
the same `ReturnType`, the solver CPU time, the largest and mean difference relative to a particle's dissolved
moles, solver statistics per call, unstable solid sets and equilibrium violations. `-vdiag` adds, per reaction,
the violations in the replayed results and whether one single-reaction solve fixes them.

Replaying a dump with the solver that wrote it must give `bitwise same` = `records`.
