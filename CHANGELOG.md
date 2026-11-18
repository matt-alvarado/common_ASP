# Changelog

Changes to `common_ASP` since the 2026 import of the model (`34b7ec7`, M. Alvarado, 25 September 2026), newest
first. Each entry gives the change, its effect on model output, and the timings that justified it.

Conventions:

- **Byte-identical** means every numerical output file of the test runs compares equal with `cmp`; `Transcript.TXT`
  and standard output may differ in timing lines only.
- Test decks: Williams, Timbavati, Otavi and Alaska (`ASP/InorganicInputDecks`, Lagrangian parcel, 5 h, 40 bins
  unless stated), and in October also RemConMode1. Timings are from a shared two-core container (noise 5–10%) unless
  stated; instruction counts are from callgrind.
- "Patch NN" and "N1"–"N13" are the labels used in the commit messages and in the three Word reports, which hold the
  detailed verification of patches 00–19: `ASP_Patch_Summary.docx` (optical-properties mode, patches 00–10 including
  N1–N13 and "Phase 1.x"), `ASP_CoagTest_Optimization_Summary.docx` (11–15) and
  `ASP_CondTest_Optimization_Summary.docx` (16–19).
- The run-mode design documents (`ASP_*_InformationFlow.md`) describe the current code; this file holds the history.
  `ASP_Timbavati_Baseline_Comparison.md` compares the current code with the imported code on one deck.

## Changes that alter results

Every other change listed below is byte-identical on the test decks.

| Commit (label) | Change | Effect |
|---|---|---|
| `a56950d` (N11) | Tabulated Gibbs–Duhem integral | about 2e-4 relative (the imported `DGAUS8` integral's own tolerance error); `GDTableMode = 0` restores the full integral |
| `7737090` (Phase 1.3) | MIEV0 for the homogeneous-sphere Mie calls | at most about 1e-5 relative in the optics |
| `1f14880` (patch 10) | `OLDASP` initialised to `.FALSE.` | initial gas concentrations of decks whose `GasPhaseChems.in` starts with 0 (CO 1.0e4 ppb from the deck instead of 75.8 ppb); parcel particle number at 10 min 151 → 224 on Williams |
| `5b4647e` (patch 16) | `IonReset` initialised in `RegridAerosol` | removes NaN bins; Timbavati and RemConMode1 CoagTest and CondTest change, Williams does not |
| `508ba8d` | Gas rate flag 5 divides by the fall-off factor | OH + HNO₃ was 248 times too fast; Timbavati HNO₃ +202%, NO₂ −79% (see the Timbavati comparison) |
| `f3f1107` | 3-body rate flag 2 reads its own parameters | decks using 3-body flag 2 (not the baseline decks) |
| `2c3b2c4` (patch 28) | Newton single-reaction electrolyte solve | particle ions within 9e-4 and water within 0.6% relative, except Timbavati at 160 bins (chloride and water −0.9%) |
| `7ffe10e` (patch 31) | Coupled electrolyte solve is the default | particle water and chloride a few percent, NH₃ about 1% (Timbavati) |
| `87a863e`, `54bf450`, `5ae6004` (patches 32, 34, 36) | Coupled-solve fixes | ThermoTest results; parcel runs mostly within 1e-4, largest Timbavati particle water −25% at 20 min (patch 32) |
| `ac30032` (patch 35) | Dry-particle solve | ThermoTest scans only; parcel runs byte-identical |
| `e3faf76` (patch 37) | Regrid before coagulation when sections are out of order | runs in which sections fall out of order, mostly at 160 bins (Williams 160 bins within 0.56%) |
| `111426c` (patch 39) | Nitrate enhancement ratios from `NitrateGroups.in` | `EnhanceRatios.txt` dTotANs/dCO when the file is present |

## 2026-11-17 — input-deck settings, test harness, Fast-TUV optics table

- **Solver and evaluation switches as input-deck settings** (`11b100d`). The environment variables `ASP_FEE_OUTER`,
  `ASP_AQ_INNER`, `ASP_FEE_PAIRRULE`, `ASP_FEE_EXITFIX`, `ASP_FEE_DEVCACHE`, `ASP_GAS_LINKED` and `ASP_THERMO_SCAN` are
  no longer read; the settings `AqOuterSolver`, `AqInnerSolver`, `AqPairRule`, `AqExitFix`, `AqDevCache`,
  `GasEvalLinked` and `ThermoScanUp` replace them with the same defaults. `ReadOptionalDeckSettings` reads optional
  `Name = value` lines after the last required line of `ASPInputDeck.in`; the four older positional lines are still
  accepted; every setting is echoed to the transcript (see `README.md`). Numerical output byte-identical (10 test cases,
  three SAM-ASP-TUV decks); each former environment setting gives the same output as the equivalent deck setting
  (16 cases).
- **Electrolyte-solver dump/replay harness moved to `Tests/`** (`98117ae`). The model keeps one hook, `FEEStateHook`,
  around each `FindElectrolyteEquilibrium` solve. Byte-identical; see `Tests/README.md`.
- **Fast-TUV optics** (`3527f57`). The Mie lookup table of `ShellRefIndAndRadFASTTUV` is the `PARAMETER` array
  `MieLUT`; volume concentrations and shell/core refractive indices come from `OpticalVolumeConcentrations` and
  `ShellAndCoreIndices`, shared with `ShellRefIndAndRad`. Byte-identical (10 cases); a test program gives bit-identical
  results on 2,760 particle states.
- Later readability work (named constants, one rate-term parser, current-tense comments, these documents): no change
  in output; see `git log`.

## 2026-11-08 — gas-phase flat arrays, nitrate groups, quieter SAM runs

- **Flat arrays for the gas-phase derivative and Jacobian** (`bf6da0b`, patch 38). `BuildGasFlat` copies
  `GasPhaseODESet` and `GasPhaseJacobian` once into arrays in list order; `GasODEEvaluator` and `GasJacobianEvaluator`
  loop over them with the same operations in the same order, and the Jacobian column loop skips empty rows (5,929 of
  616² entries are non-empty in the baseline mechanism, 2,627 of 208² in CRACMM3). `GasEvalLinked = 1` selects the
  linked lists. Byte-identical (four decks and Williams_CRACMM3; also with the linked lists). Parcel wall time within
  noise; SAM-ASP-TUV asptest 66 → 54 s, bit-identical.
- **PN and AN enhancement ratios from species lists** (`111426c`, patch 39). `GasEnhanceOutput` sums dTotPNs/dCO and
  dTotANs/dCO over the species listed in the optional deck file `NitrateGroups.in` (`ReadNitrateGroups`); without the
  file it uses the old positional sums (33 species from PAN, 77 from CH3NO3) and writes -9999 instead of reading past
  the species list. The positional AN sum read past the end of the arrays in mechanisms with fewer species after
  CH3NO3 (CRACMM3), so that column varied from run to run. Without the file: byte-identical. With it: only the AN
  column changes on the baseline decks (Williams last row 1.54946E-03 → 1.54955E-03); CRACMM3 becomes deterministic.
- **`StepPrints`** (`3df702e`). Per-step progress messages of `StepASPOnce`, `SAM_wrapper`, `InputConcentrations` and
  `StepCondensationAll` print only when `StepPrints` is set (default on; SAM-ASP-TUV turns it off). `unifac.err` is
  opened once. Results unchanged.

## 2026-11-07 — SAM-ASP-TUV interface

- **SAM-ASP-TUV interface** (`1307434`). `StepASP4SAM.f90` takes TUV photolysis rates (`PhotolysisRates`) and returns
  7-bin Fast-TUV aerosol optics (`ShellRefIndAndRadFASTTUV`, `AerosolOptPropFASTTUV`); `M3_TO_CM3` is 1e6 (it was
  1e5); `InitializeASP` takes optional deck and output directories (`SetMainDeckDirs`). Other drivers unchanged.
- **Photolysis table sized from `PhotoRates.in`** (`d038ff2`): `HowManyPhotoRxns` columns instead of a fixed 51.
  Unchanged for 51-column decks.
- **Exports for SAM-ASP-TUV**: Fast-TUV wavelengths `NFastTUVWvl`, `FastTUVWavelengths` (`947edee`); initial plume CO
  `GasPhaseCOInit` (`a13a7c6`). No change for ASP.
- **`DPREP` ordering cache with 8 slots** (`8179e16`, extends patch 21). ASP alternates between ODE systems with
  different sparsity structures; asptest hit rate 81% → 99.99%. Byte-identical; asptest 233 → 197 s.
- **Optional cell-state reset** (`af6b9ee`). `SetASPCellReset(.TRUE.)` makes each `ASPInterface` call restore the saved
  particle state (`RestoreASPCellState`) and reset the Gibbs–Duhem tables (`ResetGibbsDuhemTables`), so a call's result
  depends only on its inputs (SAM with MPI). Off by default.

## 2026-11-04 — coupled-solve robustness, coagulation ordering

- **Patch 34** (`54bf450`, not byte-identical): fewer coupled-solve fallbacks. Water is held fixed in the coupled
  solve, amounts below 1e-200 count as zero, the water equilibrium is left to the validation pass on a singular
  Jacobian, and solids are dropped cumulatively (lowest DRH first) if single drops fail. ThermoTest replay fallbacks:
  down scan 85 → 10, up scan 855 → 46. 40-bin runs against patch 32: Otavi, Alaska, Timbavati totals within 2e-11,
  Williams within 5e-5.
- **Patch 35** (`ac30032`): `DrySolve` puts the electrolyte of a dry particle (water below 1e-3 of its electrolyte,
  moles) into solids by minimising the Gibbs energy of the solid assemblage (a linear programme, dense simplex), when
  the Newton solve fails. ThermoTest fallbacks 10 → 0 and 46 → 0; 40-bin parcel runs byte-identical (not reached).
  Tried and rejected: running `DrySolve` on every dry particle first (Timbavati nitrate up to 32%, chloride 48%, water
  81% different).
- **Patch 36** (`5ae6004`, not byte-identical): reactions with ratio up to 1e300 stay in the Newton active set,
  `PreAssociate` sets a strongly sub-equilibrium dissolved electrolyte directly, up to 100 Newton iterations, and a
  subset search over left-out solids (at most 2¹⁰ solves). ThermoTest and 160-bin parcel replays: no fallbacks and no
  wet violations. 40-bin runs against patch 35: Williams within 6e-5, Otavi 8e-7, Alaska 1.9e-3 (particle water),
  Timbavati 2e-6. ThermoTest aerosol fractions change by up to 0.08.
- **Patch 37** (`e3faf76`, not byte-identical): `StepSectionalCoagulationJacobson` runs `RegridAerosol` first (up to 10
  times) when the occupied sections are not in ascending order of effective radius. The semi-implicit scheme brackets
  each combined volume between neighbouring sections and needs that order; condensation can break it, and the volume
  check then stopped the run (RemConMode1 at 60 s with a 60 s coagulation step; Alaska at 160 bins). Chosen over
  sorting inside the step (M. Alvarado). Alaska 160 bins now completes with 60 s and 30 s steps. Against patch 36:
  Otavi (40, 160 bins), Alaska and Timbavati at 40 bins byte-identical; Williams 40 within 1.1e-6; Williams 160 within
  0.56%; Timbavati 160 chloride 3.0%, particle water 3.8%, other totals within 3.5e-4. CPU time unchanged.
- **Speed against the imported code after patch 36** (Lagrangian parcel, 5 h, CPU time, two alternating runs; the
  imported code plus the `BinEdges` and `IonReset` fixes):

  | Case, bins | Imported | Patch 36 | Faster |
  |---|---|---|---|
  | Timbavati, 160 | 519.1, 518.9 s | 84.8, 85.6 s | 6.1× |
  | Williams, 160 | 196.2, 200.9 s | 44.1, 45.3 s | 4.4× |
  | Otavi, 160 | 194.8, 189.9 s | 40.8, 40.9 s | 4.7× |
  | Alaska, 80 | 304.9 s | 47.4 s | 6.4× |
  | Alaska, 40 | 153.2 s | 25.4 s | 6.0× |

  The imported code stops Alaska at 160 bins after 4 min (the coagulation volume check that patch 37 addresses).

## 2026-11-03 — coupled electrolyte solve

- **Patch 29** (`f3ac62f`): deviation cache in the worst-first loop of `FindElectrolyteEquilibriumCore`: a reaction's
  deviation is re-evaluated only when one of the amounts it depends on (`ECRDepIdx`) or the water changed. Default on
  (`AqDevCache`). Byte-identical (four decks at 40 bins, Timbavati at 160). Replay solver CPU: Timbavati 2.47 → 2.34 s,
  Alaska 1.67 → 1.04 s.
- **Patch 30** (`04736cf`, off by default, byte-identical): `CoupledElectrolyteNewton` solves all movable reactions of
  a particle together for ln(ratio) = 0 (activity coefficients, water activity and temperature fixed; analytic
  Jacobian; mass and charge exact; solids by complementarity), with a validation pass by
  `EquilibrateAqueousDissociationReaction` and a fallback to the worst-first loop. Solid-circuit rule
  (`InitSolidCircuits`, `ResolveUnstableSets`): for linearly dependent sets of 2–4 solids below their DRH, the
  thermodynamically unstable side may not be completed by precipitation. Tried and rejected on the harness:
  Gauss–Seidel sweeps (cycle in Timbavati), fewer opening sweeps (mean difference 1.7e-3), early demotion of reactions.
- **Patch 31** (`7ffe10e`, not byte-identical, accepted by M. Alvarado): `AqOuterSolver` defaults to 1 (coupled solve).
  Replay at 160 bins, records violating equilibrium by more than twice the tolerance, worst-first → coupled: Williams
  1,166 → 0, Otavi 1,195 → 0, Alaska 192 → 0, Timbavati 3,762 → 0; the worst-first violations come from the loop's
  early exit on the last reaction's deviation. 5 h CPU, worst-first → coupled: 160 bins Williams 45.7 → 41.8 s, Otavi
  48.8 → 41.9 s, Alaska 103.8 → 97.7 s, Timbavati 105.7 → 86.0 s; 40 bins 13.2 → 12.7, 13.2 → 11.3, 29.2 → 25.7,
  30.1 → 23.2 s. Output: OA, NH₄, NO₃, SO₄ totals within 2e-5 (Otavi), 2e-3 (Williams), 7e-3 (Timbavati), 1.4%
  (Alaska); particle number −0.8% in Alaska; NH₃ within 0.24% except Timbavati (−1.3% / −1.1% at 160 / 40 bins);
  Timbavati chloride and water at 290 min +7.5% / +6.7% (160 bins), −3.3% / −3.6% (40 bins), alternating in sign
  between outputs.
- **Patch 32** (`87a863e`, not byte-identical), found with the ISORROPIA II ThermoTest comparison: (1) `AerosolOutput`
  never set the letovicite index, so letovicite was missing from `TotNH4` and `TotSO4` in the ThermoTest output
  (output only); (2) a solid circuit applies only if every solid on its stable side is present or supersaturated
  (`StableSideAvailable`), which fixes the letovicite composition (up to 32% of the ammonium had left as NH₃ at 313 K);
  (3) the neutralisation step of `EquilibrateAqueousDissociationReaction` uses at most 0.99 of the water (it could
  stop a nearly dry particle with negative molality; the imported code stopped this way at 6 of 9 temperatures);
  (4) the Newton single-reaction solve restores its last valid state if its bracket collapses. 40-bin parcel runs
  against patch 31: Otavi and Alaska byte-identical, Williams within 7e-5, Timbavati totals within 1.4e-3, chloride
  3.7%, particle water −25% at 20 min and within 2.5% after 60 min.
- **Patch 33** (`13397db`): ThermoTest can scan RH upward (now `ThermoScanUp = 1`). Down scans byte-identical.
- **Speed against the imported code after patch 31**, Timbavati, 160 bins, CPU, two runs: 590.2 / 580.2 s → 95.7 /
  94.1 s (6.2× faster). Along the way at 40 bins: 160.1 s imported; 72.7 s after N1–N11; 46.9 s after the 28–29
  October work (`258ee61`); 35.2 s after patch 28; 25.8 s after patch 31.

## 2026-11-02 — electrolyte solver speed

- **Patch 26** (`5ffb5c9`): `EquilibriumConstantsRatio` returns its stored result when every value it reads is
  bit-for-bit unchanged (35% hit rate). Byte-identical; Timbavati 160 bins 158.9 → 153.1 s.
- **Patch 27** (`d8fa49b`): dump/replay harness for `FindElectrolyteEquilibrium` (now in `Tests/`); the solver body is
  `FindElectrolyteEquilibriumCore`. Byte-identical; replay reproduces every record bit for bit.
- **Patch 28** (`2c3b2c4`, not byte-identical): `EquilibrateAqueousDissociationReaction` finds the reaction extent by
  safeguarded Newton on ln(ratio) inside Jacobson's bracket instead of bisection; `AqInnerSolver = 0` restores the
  bisection (byte-identical to patch 27). Replay at 160 bins: solver CPU 3.2× less (Timbavati) and 3.5× (Alaska), same
  return type on every record, largest change 3.7e-3 of a particle's dissolved moles. Timbavati 160 bins
  125.2 → 96.0 s (three runs each); 40 bins Alaska 45.5 → 27.2 s, Timbavati 48.3 → 31.8 s, Williams 15.7 → 12.1 s.
  Totals at 290 min: ions within 9e-4 and particle water within 0.6% relative, except Timbavati at 160 bins (chloride
  and water −0.9%); for scale, 80 instead of 160 bins changes Timbavati chloride by 14%.
- Before patch 26, wall time grew about linearly with bins (Williams 9.1 / 17.8 / 60.7 s, Timbavati 14.5 / 51.1 /
  158.9 s at 10 / 40 / 160 bins), and condensation was 71% of the Timbavati 160-bin profile.

## 2026-10-30 — particle-phase reactions and heterogeneous rate flags

- **`OrgPhaseReactions.in`** (`53d416e`): first-order particle-phase organic reactions (`StepOrgPhaseReactions`, called
  from `StepCondensationAll`); rates in 1/s or as a multiple of a gas photolysis rate. Without the file nothing changes.
- **1-body rate flags 18–21** (`e5f7593`, `5880438`, `e3e101a`, `0da9214`, `c506969`): acidity-dependent IEPOX/IMAE
  uptake (`ResetAcidUptakeRate`, CMAQ AERO6/7 form, particle acidity through `AcidStateHook`), NTR2 hydrolysis
  (`ResetNTR2Rate`), N₂O₅ uptake with the Davis et al. (2008) gamma (`ResetN2O5Rate`), NO₂ uptake on fine-particle
  area. Used by the CB6r5, CRACMM and SAPRC07 decks; no change for decks without them.

## 2026-10-29 — rate-law fixes and new mechanisms

- **Gas rate flag 5 fix** (`508ba8d`, not byte-identical): k = k₁ + k₃[M]/(1 + k₃[M]/k₂); the code multiplied by the
  fall-off factor. See `ASP_Timbavati_Baseline_Comparison.md`.
- **3-body flag 2 fix** (`f3f1107`): parameters were read one position early (the flag was used as the pre-factor).
- **New rate flags** for the CB6r5, CRACMM, SAPRC07, MCM and GEOS-Chem decks (`40c1760` to `ba61060`): 1-body 13, 14,
  17; 2-body 12–21; 3-body 4, 5; photolysis flag 0 with an optional scale term (`e4888c4`). 1-body flags 15/16 were
  added and then moved to 2-body 20/21 (`ba61060`). The table of all flags is in
  `Chemistry/ChemistryParametersInitRoutines.h` (`MakeODESet`). Existing flags unchanged.
- **Patch 25** (`9e33a67`): `EquilibrateAqueousDissociationReaction` reads its `AqEquilibriaList` entries once before
  the iteration loop. Byte-identical; 30-minute Williams run 53.9 → 50.9 G instructions.
- **Speed against the imported code at patch 25** (Williams, 60 min parcel run, three alternating runs, wall): 150.9 s
  → 39.8 s (3.8× faster; the imported code uses the full Gibbs–Duhem integral and uninitialised `OLDASP`).

## 2026-10-28 — first optimisation round and fixes

Optical-properties mode (patches 00–10, details in `ASP_Patch_Summary.docx`):

- **Build fix** (`09afe88`): guard the second `ALLOCATE` of `BinEdges` in `SetSectionalParameters` (gfortran aborts).
- **Phase 1.1, 1.2** (`8cb0991`, `3717b42`): OPAC refractive-index tables as module constants; category lookups
  hoisted. Byte-identical, no measurable gain (about 0.2% of the run).
- **N1–N5** (`bc80201`): redundant work removed in `StripToken`, `MakeJacobian`, `EquilibriumConstantsRatio`,
  `MeanActivity` and `EquilibrateAqueousDissociationReaction`. Byte-identical (45-run RH sweep); −26% to −44%.
- **N6–N11** (`a56950d`): Gibbs–Duhem integrand simplified (N6), `EqCoeff` cache (N7), `pow(x,1)` skipped (N8, not
  exponent 2), cheaper rate parsing (N9), `INTENT(IN)` (N10): all byte-identical. **N11** tabulates the Gibbs–Duhem
  integral per reaction and temperature (8-point Gauss–Legendre cells in √I, cubic Hermite interpolation, 256 cells):
  not byte-identical. The table matches a tightly converged integral to under 1e-6 relative; the imported `DGAUS8`
  integral (tolerance 1e-4) is off by up to 2.4% on small values, so N11 moves results by up to 2e-4 relative in the
  optics and 1–2% in a few size-distribution bins near a water-content threshold. Optical-mode times, imported →
  N1–N11: Timbavati 13.0 → 2.9 s, Williams 2.4 → 0.8 s, RemConMode1 3.4 → 0.6 s.
- **Phase 1.3** (`7737090`): the four homogeneous-sphere Mie calls of `ShellRefIndAndRad` use MIEV0
  (`MieHomogeneousSphere`); the core-shell call stays on `DMiLay`. At most 1e-5 relative against N1–N11; Williams −27%
  instructions.
- **N12, N13** (`9d89945`, `1a144fe`): rate flag parsed once per reaction; `SetChemSetupLevel` skips the ODE set and
  Jacobian in modes that do not integrate chemistry (level 2: optical, CoagTest; level 1: ThermoTest, CondTest).
  Numerical output byte-identical; `GasPhaseODESet.txt` and `GasPhaseJacobian.txt` are no longer written in those
  modes. −16% instructions.
- **Fixes 08–10** (`8064a92`, `21e319b`, `1f14880`): `MakeMonodisperse` sets `Edges`; new particles have zeroed
  optical fields; `OLDASP` initialised (changes initial gas concentrations, see the table above). After these, valgrind
  reports no errors on Williams and RemConMode1.

CoagTest (patches 11–15, `ASP_CoagTest_Optimization_Summary.docx`), all byte-identical:

- **11** (`24670d0`): single-particle kernel quantities once per section per step (`CoagPrecomputeParticle`,
  `CoagKernPairSectional`). **12** (`7c15ef6`): `MakeSection(1,J)` instead of an uninitialised index. **13**
  (`c20e470`): species-independent terms of the cubic loops collected once per step. **14** (`cf81a74`):
  `RegridAerosol(RadiiCurrent)` skips a repeated `RecalculateRadius` pass. **15** (`89d60a6`): `VolumeFractions` built
  with one bracket scan per (I, J). Wall time: Timbavati 4.40 → 1.12 s, Williams 3.43 → 0.63 s, RemConMode1
  1.31 → 0.38 s.

CondTest (patches 16–19, `ASP_CondTest_Optimization_Summary.docx`):

- **16** (`5b4647e`): `IonReset` initialised in the grow branch of `RegridAerosol`. The uninitialised flag made the
  receiving bin's aqueous species revert, giving water-free bins with NaN ionic strength; the imported code aborted
  Timbavati and RemConMode1 in SLATEC `DGAUS8`. Williams unchanged.
- **17** (`8766c40`): the Gibbs–Duhem table is used only for a finite non-negative ionic strength (otherwise
  `DGAUS8`, as before). No change on decks that ran.
- **18, 19** (`46b4665`, `86d3a78`): per-reaction species-name flags (`InitEqRxnNameFlags`) and `PowMemo`.
  Byte-identical; −11% to −13% instructions.
- **N11 switch** (`189da4e`, patch 20): `GDTableMode` and `GDTableCells` deck settings. Tables of 256, 1024 and 4096
  cells give identical Williams output.

Lagrangian parcel (patches 21–24), byte-identical on Williams, Timbavati and RemConMode1:

- **21** (`22dbf94`): `DPREP` reuses the minimum-degree ordering (`ODRV`, 28% of the 15-minute profile) when the
  sparsity structure is unchanged. **22** (`2f83092`): one `RecalculateRadius` instead of two in the parcel particle
  loop. **23** (`97e1d64`): `GasJacobianEvaluator` no longer zeroes `pdj` (DLSODES does). **24** (`f08134b`):
  `EquilibriumConstantsRatio` evaluated once at entry of `EquilibrateAqueousDissociationReaction`. Williams 300 min:
  41.57 → 37.86 s.
- Before patches 21–24 (patches 00–19), Williams took 41.9 s against 139.9 s for the imported code; Timbavati and
  RemConMode1 did not run on the imported code (`IonReset`).
- **Sensitivity of parcel runs** (found while checking N11): perturbing the input RH by 1e-4 relative changes the
  imported code's Williams results as much as N11 does (BC −9.8%, NH₄ −12% in the 60-minute sums), because a bin sits
  near a threshold in the equilibrium solver. Differences between builds in this mode therefore cannot be judged
  against the imported code after the first minutes.
