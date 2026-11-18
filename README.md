# SETUP

## Checking Out the Repository

```
% pwd
/rd47/scratch/ASP
% git clone git@lex-gitlab.aer.com:AQAC/common_ASP.git GitCommonASP
```

This repository will be used with other model repositories and not so much on its own. Currently, common_ASP is being utilized by the ASP and STILT-ASP repositories (the master and STILT-ASP branches, respectively), and we will refer the reader to their respective READMEs for a summary on how to set up each of those models with this common set of code. Note that with ASP and STILT-ASP, common_ASP is being cloned as a submodule and is thus cloned in a "detached head" state. If common_ASP is cloned by itself, the user will start on the master branch.

# Optional settings in ASPInputDeck.in

After the last required line of `ASPInputDeck.in` (the lights-off time), a deck may give solver and evaluation
settings as `Name = value` lines, in any order, until the end of the file. Names are not case sensitive; `!` starts a
comment as elsewhere in the deck. Settings not given keep their defaults, so decks without these lines run as before.

| Setting | Default | Meaning |
|---|---|---|
| `GDTableMode` | 1 | 1 = tabulated Gibbs-Duhem integral, 0 = full DGAUS8 integral on every call |
| `GDTableCells` | 256 | table cells per reaction (16 to 65536) |
| `AqInnerSolver` | 1 | single-reaction electrolyte solver: 1 = Newton, 0 = Jacobson mass-flux bisection |
| `AqOuterSolver` | 1 | electrolyte equilibrium: 1 = coupled Newton solve, 0 = worst-first loop |
| `AqPairRule` | -1 | solid-circuit precipitation rule: 1 = on, 0 = off, -1 = on with the coupled solve (`AqOuterSolver = 1`) only |
| `AqExitFix` | 0 | worst-first loop: 0 = early exit on the last reaction's deviation, 1 = no early exit, 2 = also check every reaction at the end |
| `AqDevCache` | 1 | worst-first loop deviation cache: 1 = on, 0 = off (same results) |
| `GasEvalLinked` | 0 | gas-phase ODE and Jacobian evaluation: 0 = flat arrays, 1 = linked lists (same results) |
| `ThermoScanUp` | 0 | ThermoTest RH scan: 0 = down 0.90 to 0.20, 1 = up 0.20 to 0.90 (use a deck with initial RH 0.20) |

Example:

```
360   !! When to turn off lights in minutes
AqOuterSolver = 0    ! worst-first loop
AqPairRule    = 1
```

Older decks may give the first four settings as positional lines instead (one value per line, in the order
`GDTableMode`, `GDTableCells`, `AqInnerSolver`, `AqOuterSolver`); positional lines must come before any `Name = value`
line. The model writes every setting it uses to `Transcript.TXT` and standard output at the start of the run.

These settings replace the environment variables `ASP_FEE_OUTER` (`AqOuterSolver`), `ASP_AQ_INNER` (`AqInnerSolver`),
`ASP_FEE_PAIRRULE`, `ASP_FEE_EXITFIX`, `ASP_FEE_DEVCACHE`, `ASP_GAS_LINKED` (`GasEvalLinked`) and `ASP_THERMO_SCAN`
(`ThermoScanUp`), which the model no longer reads. The electrolyte-solver dump and replay test is in `Tests/`
(see `Tests/README.md`).

# Documentation

| File | Contents |
|---|---|
| `ASP_LagrangianParcel_InformationFlow.md` | Lagrangian parcel model (gas + aerosol): call flow, data flow, key routines, cost |
| `ASP_CondTest_InformationFlow.md` | condensation test mode; condensation and the electrolyte solver in detail |
| `ASP_CoagTest_InformationFlow.md` | coagulation test mode; `StepSectionalCoagulationJacobson` and `RegridAerosol` in detail |
| `ASP_OpticalPropertiesModel_InformationFlow.md` | optical properties model; refractive indices and Mie calls |
| `ASP_Timbavati_Baseline_Comparison.md` | speed and results of the current code against the imported code on the Timbavati deck |
| `CHANGELOG.md` | dated list of changes since the 2026 import, with their effect on results and timings |
| `ASP_Patch_Summary.docx`, `ASP_CoagTest_Optimization_Summary.docx`, `ASP_CondTest_Optimization_Summary.docx` | detailed verification reports for the October 2026 changes (patches 00–19 in `CHANGELOG.md`) |
| `Tests/README.md` | electrolyte-solver dump and replay test |

The gas-phase rate-law flags (number of terms, formula, and the routine that evaluates each) are tabulated in
`Chemistry/ChemistryParametersInitRoutines.h`, above the gas-phase section of `MakeODESet`.

# Open items

- **`SurfaceTension` cost.** `RecalculateRadius` calls `SurfaceTension` (via `KusikMeissner` and `CalcLogGammaPure`)
  for every section; it is about half of a CoagTest run and about 9% of a parcel run. A cheaper evaluation changes how
  the thermodynamic state is used and must be tested on every run mode that calls coagulation. Not started.
- **Gas solver set-up.** `StepGasChemistry` calls `DLSODES` with `mf = 21` and `istate = 1` on every call while also
  loading `GasIa`/`GasJa` into `iwork`; whether LSODES uses the supplied structure has not been checked. `rwork` (at
  least 5,000,000 doubles) is allocated on every call. `DPRJS` zeroes a full column per Jacobian column (about 4% of a
  parcel run); clearing only the touched rows would need `GasJacobianEvaluator` to report them.
- **Speed-ups that would change results** (not taken): `x**y` with small integer `y` as multiplication, continuing
  `DLSODES` with `istate = 2`, and reusing the sparsity structure found by differencing.
- **Thread safety.** The `EqCoeff` cache, the Gibbs–Duhem tables, `PowMemo`, the `EquilibriumConstantsRatio` memo,
  the `DPREP` ordering cache and MIEV0's first-call state are module-level and filled lazily; none is safe under
  OpenMP as written.
- **Rate-flag checking.** `MakeODESet` does not reject an unknown 1- or 2-body rate flag (the reaction keeps all-zero
  rate terms); an unknown 3-body flag stops the run.
- **RemConMode1 coagulation step.** The deck uses a 30 s coagulation step because sections fell out of size order with
  60 s; `StepSectionalCoagulationJacobson` now regrids in that case, but the deck has not been re-run with 60 s.
- **Standard output in parcel runs.** The driver loop prints several progress lines per step (about 207,000 lines for a
  300-minute run); they are not under `StepPrints`.
- **Valgrind** has not been re-run since the October coagulation changes (an out-of-bounds read in the old
  `VolumeFractions` loop was removed by rewriting it).
