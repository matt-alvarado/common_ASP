# ASP Model — Information Flow for `DevelopmentDriver.f90` with the Lagrangian Parcel Model (gas + aerosol)

Source basis: `common_ASP` branch `asp_v_2.2.1`. Sections 1–6 were written against commit `fdf2da9` (patches 00–19) and their line numbers and profiles refer to that tree; section 7 is updated for patches 20–25 (branch tip `65d2ac9`, tag `v2.2.1`).

## 1. Scope and purpose

This document traces control and data flow when:

- `DevelopmentDriver.f90` (`PROGRAM ASP`) is the main program,
- the deck line "Lagrangian Parcel Model (1 = yes)" is `1` (all other run-mode lines `0`), and
- "Include aerosols?" is `1`, so the **gas + aerosol** branch (`ELSE IF (IncludeAerosols)`, line 567) runs. The gas-only branch (lines 400–566) is described briefly in section 3.9.

The parcel is a single well-mixed air parcel that is diluted with a background atmosphere. Every step it does gas-phase chemistry (stiff sparse ODE solve), mixes gases and particles with the background, and every 60 s it runs aerosol condensation/dissolution and coagulation on a sectional (binned) size distribution. Unlike the CoagTest and CondTest modes, nothing is held static: all processes interact, so cost is spread across the gas solver and the aerosol thermodynamics. Section 5 gives a measured profile; section 7 lists optimization hypotheses (none implemented).

## 2. High-level flow

```mermaid
flowchart TD
    A[PROGRAM ASP] --> B[SetDomainSize, ReadMainInputDeck]
    B --> C{LagrangianParcelModel and IncludeAerosols}
    C --> D[SetChemistryParams + ReadPhotolysisRates\nfull gas, aqueous and ODE-set tables]
    D --> D1[RecalculateReactionRates,\nRecalculatePhotoRates at UTCStartTime]
    D1 --> E[SetAllDissolution]
    E --> F[InitializeParticlesSectional\nsectional population from AerosolModes.in]
    F --> F1[ReadBoundaryDistributionsOrganic +\nPopulateParticlesSectionsRightAway TRUE\nbackground aerosol from EnvAerosolModes.in]
    F1 --> G[EquilibrateInternallyAtGridPoint water, RegridAerosol]
    G --> H[Hetero rates from size distribution,\nopen output files, write t = 0 output]
    H --> L{{Main loop: 600 + (RunLength-10)*60/MixStep iterations}}
    L --> L1[RecalculateReactionRates, RecalculatePhotoRates,\nHeteroNumberConc and radii, RecalculateHeteroRates]
    L1 --> L2[StepGasChemistry\nDLSODES, ChemTimestep]
    L2 --> L3{t multiple of CondTimeStep, t not 0}
    L3 -- yes --> L4[RegridAerosol, FindAqElectrolyteEquilibriumForGridPoint,\nStepCondensationAll]
    L3 -- no --> L5
    L4 --> L5{t multiple of CoagTimeStep, t not 0}
    L5 -- yes --> L6[StepSectionalCoagulationJacobson]
    L5 -- no --> L7
    L6 --> L7[Dilution of gases: T_mix = Slope*t + b]
    L7 --> L8[UpdateChemicalConcentrations]
    L8 --> L9[Sort particles, per bin: RecalculateRadius,\nTerminalVelocity, dilute + deposit number and mass,\nRecalculateRadius]
    L9 --> L10{t multiple of 600 s}
    L10 -- yes --> L11[ShellRefIndAndRad per bin, AerosolOptProp,\nwrite AeroOptCIS min file]
    L10 -- no --> L12
    L11 --> L12{t multiple of 60 s}
    L12 -- yes --> L13[GasOutput, GasEnhanceOutput,\nAerosolOutput for every bin]
    L12 -- no --> L14
    L13 --> L14[SizeDistOutput at t = 600 s multiples, 540 s, 2820 s]
    L14 --> L15[t = t + MixTimestep, ut = ut + MixTimestep/3600]
    L15 --> L
```

## 3. Step-by-step description

### 3.1 Program entry and branch selection

`ReadMainInputDeck` reads `ASPInputDeck.in`. The Lagrangian parcel branch is selected at line 400, and the aerosol sub-branch at line 567. The mixing parameter `b_mix` (deck section 3.b, 106.87 s in the three decks used) sets the dilution time.

### 3.2 Chemistry setup (full level)

`SetChemistryParams` builds all gas, aqueous and ODE-set tables (this branch does not use the reduced `SetChemSetupLevel(1)` of CondTest). `ReadPhotolysisRates` reads `PhotoRates.in`. `RecalculateReactionRates(IFGAS,IFAQ)` computes temperature-dependent rate constants and `RecalculatePhotoRates(ut, t, .TRUE.)` the initial photolysis rates at `UTCStartTime`.

### 3.3 Dissolution data and the particle population

- `SetAllDissolution`: inorganic and organic (hydrophobic and hydrophilic) dissolution data.
- `InitializeParticlesSectional` (line 592) builds the sectional population from `AerosolModes.in` (log-normal modes, section cutoffs). `RegridAerosol` follows.
- Lines 655–662: `ReadBoundaryDistributionsOrganic` and `PopulateParticlesSectionsRightAway(.TRUE.)` build the **background (boundary) aerosol** from `EnvAerosolModes.in`, then `SortBoundaryPart`. This is the aerosol the parcel mixes with.
- `EquilibrateInternallyAtGridPoint(EquilibrateWater=.TRUE., WaterOpenSystem=.TRUE.)` then `RegridAerosol` bring every bin to water equilibrium. The bins are then counted (`NumBins`, `NumAqChems`, `NumOrgChems`, `NumAqOrgChems`) and the work arrays (`NumConc`, `InitNumConc`, `EnvNumConc`, `HeteroNumberConc`, `HeteroRadii`) are allocated.
- Lines 745–767: the heterogeneous-reaction number concentration and radii are taken from the sections (`EffectiveRadius/100`, cm to m) and `RecalculateHeteroRates(..., .TRUE.)` is called.

### 3.4 Output at t = 0

`GasOut1/2/3.txt`, `EnhanceRatios.txt`, one `AerosolN.txt` per bin, and `SizeDist0min.txt`.

### 3.5 Main time loop (lines 841–1090)

The loop runs `INT(600 + (RunLength-10)*60/MixStep)` iterations. For the first 600 s the steps are 1 s (`ChemTimestep = MixTimestep = 1`); afterwards they are reset to the deck values (`ChemStep`, `MixStep`, 60 s in the test decks). Each iteration:

1. **Rates.** `GetTemp`, `RecalculateReactionRates`, `RecalculatePhotoRates(ut, t, .FALSE.)`. Heterogeneous number concentrations and radii are refreshed from the particle list, then `RecalculateHeteroRates`.
2. **Gas chemistry.** `StepGasChemistry(ChemTimestep, INT(MixTimestep/ChemTimestep))` (`Chemistry/EvolveGasChemistry.h`) integrates the gas ODEs with the sparse stiff solver `DLSODES` (`istate = 1` on every call, `mf = 21`, `iwork(6) = 2000`). The routine allocates `Y`, `rwork` (at least 5,000,000 doubles) and `iwork` on every call.
3. **Condensation (every `CondTimeStep`, not at t = 0).** `RegridAerosol`, `FindAqElectrolyteEquilibriumForGridPoint(UpdateThermo=.TRUE.)`, then `StepCondensationAll(CondTimeStep, 1)` (the same machinery as CondTest section 3.7 of `ASP_CondTest_InformationFlow.md`: LSODES for H2O/H2SO4 plus the inorganic equilibrium solver for each bin). `StepCondensationAll` calls `RegridAerosol` itself.
4. **Coagulation (every `CoagTimeStep`, not at t = 0).** `StepSectionalCoagulationJacobson(CoagTimeStep)` (same as CoagTest, `ASP_CoagTest_InformationFlow.md`). It stops the run with "Total Volume concentration is lost during step" if the total volume changes by more than 1e-4 (section 6).
5. **Gas dilution and deposition.** `T_mix = Slope*t + b` with `Slope = 2.0`, `b = b_mix`. For every gas species: `Chem(q) = (GridGasChem(q) + Background(q)*dt/T_mix) / (1 + dt/T_mix + dt*DepVel(q)/InvHgt)`. `UpdateChemicalConcentrations` stores the result.
6. **Particle dilution and deposition.** `SortAerosolAtGridPointForCoagulation`, then for each bin, with the matching background bin `Env`: `RecalculateRadius`, `TerminalVelocity`, the new number `N = (N0 + Nenv*dt/T_mix) / (1 + dt*(1/T_mix + TermV/InvHgt))`, the same formula for the mass of every aqueous, organic and aqueous-organic species (divided by the new number to get moles per particle), `Temperature = TEMP`, and a second `RecalculateRadius` (lines 943 and 1005).
7. **Optical properties (every 600 s).** `ShellRefIndAndRad` for each bin, `AerosolOptProp` (MIEV0/DMiLay Mie code), and `AeroOptCISNmin.txt` (451 wavelengths).
8. **Output (every 60 s).** `GasOutput`, `GasEnhanceOutput`, `AerosolOutput` for each bin; `SizeDistOutput` at multiples of 600 s, and at 540 s and 2820 s.
9. **Time update.** `t = t + MixTimestep`, `ut = ut + MixTimestep/3600`, `ChemTimestep` and `MixTimestep` reset after 600 s. The loop also writes many `WRITE(*,*)` lines per iteration ("In Loop", "Photo Okay", ..., "Before AerosolOutput, Bin #"); the Williams 300-minute run writes about 207,000 lines (9 MB) to stdout.

### 3.6 Key routines and where they live

| Stage | Routine | File |
|---|---|---|
| Gas ODE solve | `StepGasChemistry`, `GasODEEvaluator`, `GasJacobianEvaluator`, `DLSODES` | `Chemistry/EvolveGasChemistry.h`, `ComputationalSubroutines/opkdmain.f`, `opkda1.f` |
| Condensation | `StepCondensationAll`, `StepCondensation`, `EquilibrateGridPoint` | `Condensation/CondensationIntegrator.h`, `CondensationRelatedFunctions.h` |
| Inorganic equilibrium | `FindElectrolyteEquilibrium`, `EquilibrateAqueousDissociationReaction`, `EquilibriumConstantsRatio`, `MeanActivity`, `KusikMeissner` | `LagrangianParticles/Thermodynamics.h` |
| Particle properties | `RecalculateRadius`, `ParticleDensity`, `SurfaceTension` | `LagrangianParticles/ParticleAttributes.h` |
| Coagulation | `StepSectionalCoagulationJacobson` | `Coagulation/Coagulation.f90` |
| Regrid | `RegridAerosol`, `SortAerosolAtGridPointForCoagulation` | `LagrangianParticles/SortAerosol.h` |
| Optics | `ShellRefIndAndRad`, `AerosolOptProp`, `MIEV0`, `DMiLay` | `ParticleAttributes.h`, `MIEV0.f`, `DMiLay.f` |

### 3.9 Gas-only branch (not profiled here)

With "Include aerosols?" = 0 (lines 400–566), the same loop runs only the rates, `StepGasChemistry`, gas dilution and gas output, with a single heterogeneous bin from the deck (`hetero_num`, `hetero_rad`).

## 4. Inputs consumed and outputs produced

**Inputs** (in `InorganicInputDecks`):

| File | Used for |
|---|---|
| `ASPInputDeck.in` | run mode flags, T, P, RH, run length, the four time steps, inversion height, `b_mix`, photolysis settings |
| `AerosolModes.in` | modes (Nt, geometric mean diameter, sigma, composition), number of sections and cutoffs |
| `EnvAerosolModes.in` | background aerosol the parcel mixes with |
| `GasPhaseChems.in`, `GasChemicalMechanism.in`, `PhotoRates.in` | gas species (with initial and background values), reactions, photolysis |
| `AqPhaseChems.in`, `AqPhaseIons.in`, `AqEquilibriumReactions.in`, `AqChemicalMechanism.in`, `AerosolThermoKusikMeissner.in` | aqueous species, ions, equilibria and activity-coefficient data |
| `Dissolution.in`, `OrganicDissolution.in`, `HydrophilicOrganicDissolution.in`, `OrgPhaseChems.in`, `HydrophilicOrgChems.in` | gas-aerosol partitioning and organics |

**Outputs** (`OutputFiles`, 109 files for the 300-minute Williams run): `GasOut1/2/3.txt` and `EnhanceRatios.txt` (every 60 s), `Aerosol1..40.txt` (one per bin, every 60 s), 32 `SizeDistNmin.txt`, 30 `AeroOptCISNmin.txt`, `GasPhaseJacobian.txt`, `GasPhaseODESet.txt` and `Transcript.TXT`.

## 5. Baseline speed test and measured profile

### 5.1 Test setup

Three decks (Timbavati, Williams, RemConMode1; the optical-properties templates with Lagrangian Parcel = 1 and Optical Properties = 0). Common settings: T = 290 K, P = 1000 mb, RH = 0.7, run length 300 min, 60 s chemistry, mixing, condensation and coagulation steps, 40 bins. Builds: gfortran 16, `-O3 -g`, two-core container. Wall time, three runs per entry unless noted (noise about 5%).

### 5.2 Wall time

| Deck | Original (34b7ec7) | Patches 00–N5 (bc80201, 1 run) | Patches 00–09 (21e319b, 1 run) | After patch 10 (1f14880, 1 run) | Current (patches 00–19) |
|---|---|---|---|---|---|
| Williams | 139.9 s (139.4, 143.7, 136.7) | 102.3 s | 47.7 s | 48.1 s | 41.9 s (41.5, 42.4, 42.0) |
| Timbavati | aborts after about 1 s (SLATEC `DGAUS8`) | — | — | — | 46.8 s (46.8, 46.4, 47.1) |
| RemConMode1 | aborts after about 0.7 s (SLATEC `DGAUS8`) | — | — | — | 40.6 s with a 30 s coagulation step (stops after about 1.3 s with the 60 s step: coagulation volume check, section 6) |

So the baseline for Williams is about 140 s; the current branch is 3.3 times faster. Timbavati and RemConMode1 have no usable baseline because the original aborts on both (the `IonReset` bug, patch 16). The intermediate builds show where the Williams gain came from: N1–N5 save about 27%, patches N6–N11 and the Mie, ODE-set and run-mode changes cut a further 54%, and patches 11–19 about 13% more.

### 5.3 Output comparison with the baseline

- Patches 00 to N5 (bc80201): all 108 numerical output files on Williams are **byte-identical** to the original build.
- From N6–N11 (`a56950d`) on, output is no longer identical. At 10 min, total particle number, BC, water and K match the original to 1e-6 (SizeDist files), but by 60 min the sums differ by up to 10% in some columns (BC −10%, NH4 −12%) and by 290 min by up to 8%.
- **Cause of the drift (tested).** Building `a56950d` with only the N11 tabulated Gibbs–Duhem integral switched off (always DGAUS8) gives output byte-identical to the original on all 108 files, so N6–N10 are exact and N11 alone starts the drift. N11 changes the activity integral by about 2e-4, which shows up as 1e-6 to 4e-4 differences in individual bins at t = 0.
- **Why it grows (tested).** The same drift appears in the *original* build when the input RH is changed by only 1e-4 (relative; 0.7 to 0.70007, initial and background): same bin, same time, and the same 60-minute sums (BC −9.8%, NH4 −12%) as the N11 build. A change of 1e-6 instead gives differences of only about 4e-4 even at 290 min. In bin 30 the original shows a one-step composition event at t = 1500 s (solid (NH4)2SO4, K2SO4 and KNO3 jump, for example (NH4)2SO4 1.5e-22 against 1e-28 mol per particle) while particle number and radius are unchanged; after that K2SO4 stays about 70% above the perturbed runs. So the reference run sits near a threshold in the equilibrium solver, and any perturbation of about 1e-4 or more (N11, or a tiny RH change) tips it the other way. I did not identify which decision in the solver flips.
- Consequence: the drift is a sensitivity of the model, not evidence that N11 is wrong, but it means differences between builds in this mode cannot be judged by comparing files against the original after the first minutes. Gas and aerosol totals stay within a few percent of each other for 60 minutes, then individual bins differ at the order-one level where an event differs.
- **Switch (commit `189da4e`).** The last two optional lines of `ASPInputDeck.in` set `GDTableMode` (1 = N11 table, 0 = full DGAUS8 integral on every call, which reproduces the pre-N11 quadrature) and `GDTableCells` (default 256). Tables of 256, 1024 and 4096 cells give identical Williams output in the current build, so the table is already converged and the 2e-4 gap to the original comes from the original quadrature's own tolerance, not from table error. In the current build (after the OLDASP fix) table versus full integral agrees to 3e-5 in the enhancement ratios through 5 h, but individual bins still differ by up to 0.5% in water (bin 26, 60 min) and 1% in solid (NH4)2SO4 (bin 25, 5 h).
- Patch 10 (the approved OLDASP fix, `1f14880`) then changes the results a lot: total particle number in the SizeDist file at 10 min goes from 151 to 224 (+48%), presumably because `OLDASP` was uninitialised in the original and took a different path (I did not check the stack value). The patch is deliberate, but this shows its effect on parcel runs is far larger than in the optical and condensation modes.
- Consequence: a byte-identical check against the original build is only possible up to bc80201 in this mode. For later optimizations the reference must be the current branch.

### 5.4 Profile (Williams, first 15 minutes, callgrind, `-O3 -g`)

The run is 605 iterations: 600 with 1 s steps, 14 condensation and coagulation steps. Total 53.6 G instructions (inclusive shares):

| Part of the run | Instructions | Share |
|---|---|---|
| `StepGasChemistry` (605 calls, `DLSODES`) | 24.6 G | 45.9% |
| &nbsp;&nbsp;`DPREP`/`DIPREP` (815 calls, sparse set-up) | 15.9 G | 29.7% |
| &nbsp;&nbsp;&nbsp;&nbsp;`ODRV`/`MD` (ordering, minimum degree) | 14.9 G | 27.8% |
| &nbsp;&nbsp;`DPRJS` (4,329 calls, Jacobian) | 7.1 G | 13.2% |
| &nbsp;&nbsp;&nbsp;&nbsp;`GasJacobianEvaluator` (372,680 calls) | 3.6 G | 6.8% |
| &nbsp;&nbsp;`DSTODE` (integration steps) | 8.6 G | 16.0% |
| `StepCondensationAll` (14 calls) | 15.5 G | 28.9% |
| &nbsp;&nbsp;`FindElectrolyteEquilibrium` | 15.3 G | 28.6% |
| &nbsp;&nbsp;`EquilibrateAqueousDissociationReaction` | 12.8 G | 23.8% |
| &nbsp;&nbsp;`EquilibriumConstantsRatio` | 9.6 G | 17.9% |
| `RecalculateRadius` (24,200 calls from the dilution loop, two sites) | 4.9 G | 9.2% |
| &nbsp;&nbsp;`SurfaceTension` (all callers) | 5.9 G | 11.1% |
| `ShellRefIndAndRad` + Mie (`DMiLay`, `MIEV0`) | 1.9 G | 3.5% |
| `EquilibrateInternallyAtGridPoint` (set-up) | 2.9 G | 5.4% |
| `StepSectionalCoagulationJacobson` (14 calls) | 0.20 G | 0.4% |
| `RegridAerosol` | 0.6 G | 1.1% |
| `memset` (zeroing) | 4.0 G | 7.5% |

The `memset` total is almost all in two places: `DPRJS` (1.85 G) and `GasJacobianEvaluator` (1.84 G). The profile covers only the first 15 minutes; the shares over the full 300 minutes will differ (for example, the 1 s chemistry steps are all in the first 10 minutes) and I did not profile them.

## 6. Compatibility and issues found while preparing this document

1. **Original build aborts on Timbavati and RemConMode1** in the first condensation step (SLATEC `DGAUS8`, NaN ionic strength). This is the `IonReset` bug fixed in patch 16 (`ASP_CondTest_Optimization_Summary.docx`), which also hits this run mode through `RegridAerosol`.
2. **RemConMode1 stops at the first coagulation step (t = 60 s)** with "Total Volume concentration is lost during step" when the coagulation step is 60 s. The check is two-sided and here volume is slightly *gained*: the ratio is 1.000121 against a tolerance of 1e-4. It happens in builds without my coagulation patches 11–15, so those are not the cause. With the state at t = 60 s fixed, the mismatch is proportional to the coagulation step (60 s: 1.21e-4; 30 s: 6.0e-5; 15 s: 3.0e-5; 6 s: 1.2e-5), and every other call in the run conserves volume to about 1e-9. The first condensation step moves nearly all particles from bin 4 into bin 23 and a few near-empty bins (2, 3, 27, 28; about 1e-4 particles/cm³) end up with radii 2–4 times those of their neighbours, so the particle volumes handed to coagulation are not in increasing order. The partition coefficients in `StepSectionalCoagulationJacobson` assume increasing `SingleParticleVolume` (they bracket the combined volume between bins K and K+1), and the coagulation sort orders bins by their edges, not by their actual volume. I believe this is what trips the check but have not proved it. **Workaround adopted: the RemConMode1 deck uses a 30 s coagulation step** (the run then completes, 40.6 s, no errors). The proper fix, ordering by actual volume, is described in `ASP_Pending_Patch_Notes.md` and is not implemented.
3. **The Williams output no longer matches the original build** after N11 and the OLDASP fix (section 5.3).
4. **Stdout volume.** About 207,000 lines (9 MB) are written for the 300-minute Williams run, mostly "Okay" and "Before AerosolOutput, Bin #" progress messages. I did not measure their cost.
5. **Allocation.** `StepGasChemistry` allocates `rwork` of at least 5,000,000 doubles (40 MB) on each of its 605+ calls. The profile shows no large cost from this, but it is a candidate for a one-time allocation.
6. **LSODES options.** `mf = 21` with `istate = 1` on every call, while the code also loads `GasIa`/`GasJa` (sparsity structure) into `iwork`. I did not check which LSODES structure option `mf = 21` selects or whether the supplied structure is used; this matters for hypothesis 1.

## 7. Optimization status and next candidates

### 7.1 Implemented after this profile (byte-identical output)

| Patch | Commit | What | Evidence |
|---|---|---|---|
| 21 | `22dbf94` | `DPREP` caches the minimum-degree ordering (`ODRV`) and reuses it when the sparsity structure (`IA`, `JA`) is unchanged. Hypothesis 1 of the earlier list. | `ODRV` was 28% of the 15-minute profile. `ODRV` is deterministic in `N` and the structure, so the cached permutation is the one it would compute. |
| 22 | `2f83092` | Removed the first `RecalculateRadius` in the parcel particle loop. Hypothesis 4. | Instrumented build: 35,600 calls on each of Williams, Timbavati and RemConMode1, 0 bitwise differences in the 11 scalar fields it sets. |
| 23 | `97e1d64` | `GasJacobianEvaluator` no longer re-zeroes `pdj`, which `DPRJS` has already zeroed (LSODES contract). Partly hypothesis 2. | Only `DLSODES` calls this routine. |
| 24 | `f08134b` | `EquilibrateAqueousDissociationReaction` evaluates `EquilibriumConstantsRatio` once at entry and reuses it for the second check and the first convergence test (it was evaluated up to three times on an unchanged particle). Part of hypothesis 3. | `EquilibriumConstantsRatio` reads the particle only; reuse is skipped if the neutralisation step changed it. |
| 25 | `9e33a67` | `EquilibrateAqueousDissociationReaction` reads the `AqEquilibriaList` entries (columns 1–5, 8, 9, 22–25), the `INT` index arithmetic, the name flags and the half-iteration limit once before the iteration loop, and the loop uses the locals. `isnan(EqRatio)` is inlined as `EqRatio .NE. EqRatio`. Candidate 1 of the earlier list. | `AqEquilibriaList` is written only at initialisation; the locals are non-SAVE (the routine is `RECURSIVE`); `IsNaN(a)` is `a .ne. a`. |

All five give identical output on all 109 files of Williams, Timbavati and RemConMode1 (the latter with the 30 s coagulation step, section 6 item 2). Measured on Williams, quiet machine, three runs each: **41.57 s (42.3, 41.8, 40.7) before patch 21 to 37.86 s (37.9, 37.9, 37.8) after patch 24, 8.9% faster**. A 30-minute run falls from 59.4 G to 53.9 G instructions between patch 21 and patch 24 (patches 22–24, 9%), and to **50.9 G after patch 25 (a further 5.7%)**; Williams wall time is within noise of that (about 39 s on a busy machine).

### 7.2 Profile after patch 24 (30 minutes, 53.9 G instructions; a 90-minute profile before patches 22–24 gives the same ranking)

Over a longer run the mix shifts from gas chemistry to condensation: in the 90-minute profile `StepCondensationAll` is 69% and `DLSODES` 10% (it was 46% in the 15-minute profile, which is dominated by the first ten minutes of 1 s steps). After patch 24: `StepCondensationAll` 58%, `EquilibriumConstantsRatio` 31%, `MeanActivity` 18%, `KusikMeissner` 15%, libm `pow` 10%, `SurfaceTension` 9.4% (all callers), `DLSODES` 17%, `memset` 4.2%. The largest single self cost is `EquilibrateAqueousDissociationReaction` (15%).

### 7.3 Next candidates

1. **Done as patch 25** (see 7.1): 53.9 G to 50.9 G instructions, byte-identical.
2. **`SurfaceTension` (9.4%).** This is the option A parked earlier "until we test other run modes that call coagulation". Parcel mode is such a run mode. Waiting for your decision; any simplification changes how often the thermodynamic state is used and must be tested on this mode.
3. **Zeroing in `DPRJS` (memset 4.2%).** After patch 23 one full zeroing of `FTEM` per Jacobian column remains. Clearing only the entries the Jacobian wrote is not safe byte-for-byte: the sparsity pattern is found by differencing at the start of each call, and `GasJacobianEvaluator` can write entries outside it, so those would stay dirty. Needs `GasJacobianEvaluator` to report which rows it touched.
4. **Not byte-identical, listed for completeness:** replacing `x**y` for small integer exponents by multiplication (libm `pow` is 10% but the last bit can differ), continuing `DLSODES` with `istate = 2` (changes the step-size history), and re-using the sparsity structure found by differencing (it depends on the state).
5. **Per-call allocation in `StepGasChemistry`** (about 40 MB of `rwork` per call, no visible cost in the profile), **stdout writes** (not measured), **optics** (3.5%, 30 calls) and **coagulation** (0.4%): not worth attention in this mode.

### 7.4 Overall speed against the baseline (Williams, 60 min parcel run)

Alternating solo runs on the same machine, wall time:

| Build | Run 1 | Run 2 | Run 3 | Mean |
|---|---|---|---|---|
| Baseline (`master`, plus the one-line `BinEdges` gfortran build fix) | 149.0 s | 150.7 s | 153.0 s | **150.9 s** |
| Branch tip (patch 25) | 40.6 s | 40.9 s | 37.8 s | **39.8 s** |

That is 73.6% less wall time (3.8 times faster). The comparison is not between identical physics: the baseline uses the original `DGAUS8` Gibbs–Duhem integral (about 2e-4 away from the table, section 7.1 and the N11 bullet in section 6) and reads uninitialised memory for the initial gas concentrations (patch 10). Most of the gain is the tabulated integral plus the speed patches 16–25. The baseline aborts after about 1 s on RemConMode1 and Timbavati (SLATEC `DGAUS8`), so these decks have no baseline wall time. Instruction counts for the full baseline run were attempted and stopped as too slow under callgrind; the only exact instruction counts are the 30-minute figures in section 7.1 and 7.2 (59.4 G at patch 21, 53.9 G at patch 24, 50.9 G at patch 25).

### 7.5 Size-bin scaling (patches 26–28)

Baseline deck, Lagrangian parcel, 5 h, one core, `common_asp` c506969. Wall time grows about linearly with the number of bins (N^0.8–0.9 from 40 to 160 bins): Williams 9.1 / 17.8 / 60.7 s, Timbavati 14.5 / 51.1 / 158.9 s at 10 / 40 / 160 bins. In Timbavati at 160 bins condensation is 71% of the profile, almost all of it the electrolyte solver: `EquilibrateGridPoint` always runs its outer loop to the 11-iteration cap, and each `FindElectrolyteEquilibrium` call takes about 45 worst-reaction-first iterations, each re-evaluating every major reaction, with a bisection (mass-flux iteration) for each single-reaction solve.

| Patch | What | Evidence |
|---|---|---|
| 26 | `EquilibriumConstantsRatio` keeps, per reaction, the last result and the bit patterns of every value it reads (electrolyte, cation, anion and letovicite-anion amounts, water amount, number concentration, `Dry`, temperature, water activity, the two activity coefficients) and returns the stored result when all are unchanged. | Byte-identical (229 files, Timbavati 160 bins). Hit rate 35%; Timbavati 160 bins 158.9 s to 153.1 s (one run each). |
| 27 | Test harness for `FindElectrolyteEquilibrium`. With `ASP_FEE_DUMP=<file>` (and optionally `ASP_FEE_DUMP_EVERY=<n>`, default 100) every n-th condensation-path call writes the full particle state it reads and the state it returns, plus grid temperature, pressure and water vapour. With `ASP_FEE_REPLAY=<file>` a run initialises, re-solves every dumped state with the solver currently selected, times only the solver calls (CPU), writes `FEEReplay.txt` (records, bitwise matches, `ReturnType` matches, solver CPU, largest and mean difference per particle relative to its total dissolved moles) and stops. The solver itself is now `FindElectrolyteEquilibriumCore`; `FindElectrolyteEquilibrium` is a wrapper. | With no variable set the run is byte-identical. Replaying with the unchanged solver reproduces all 5,262 Timbavati and 5,254 Alaska records (160 bins) bit for bit, so the dump holds every input. |
| 28 | `EquilibrateAqueousDissociationReaction` finds the reaction extent by safeguarded Newton on ln(ratio) instead of bisection: the derivative is the analytic d ln Q / dX with activity coefficients and water activity fixed (stoichiometry² / amount for each ion and the dissolved electrolyte, plus the water-amount term), steps stay inside Jacobson's bracket [−Qnumer, Qdenom] with a bisection step whenever Newton would leave it, every trial state is starting amounts + stoichiometry × X (mass and charge conserved exactly), and an undersaturated solid that would need more than all of itself is dissolved completely. Convergence test and early returns are unchanged. **Not byte-identical; default on.** `AqInnerSolver` (optional `ASPInputDeck.in` line after `GDTableCells`: 1 = Newton, default; 0 = original bisection), overridden by `ASP_AQ_INNER=newton` or `jacobson`. | `AqInnerSolver = 0` (or `ASP_AQ_INNER=jacobson`) is byte-identical to patch 27 (4 cases × 109 files at 40 bins). Harness, 160 bins: solver CPU 3.2× less on Timbavati and 3.5× on Alaska, same `ReturnType` on every record, largest change 3.7e-3 of a particle's dissolved moles (mean 2.6e-5). See the tables below. |

Timbavati, 160 bins, three alternating runs each:

| Build | Run 1 | Run 2 | Run 3 | Mean |
|---|---|---|---|---|
| c506969 (before patch 26) | 125.3 s | 125.2 s | 125.2 s | **125.2 s** |
| Patch 28, Newton | 96.4 s | 95.0 s | 96.5 s | **96.0 s** |

That is 23% less wall time. At 40 bins (one run each, original solver → Newton): Williams 15.7 → 12.1 s, Otavi 12.5 → 12.5 s, Alaska 45.5 → 27.2 s, Timbavati 48.3 → 31.8 s.

Output change with Newton, totals over bins at 290 min (relative): OA and sulfate below 1e-4 in all cases; ammonium and nitrate below 7e-4; chloride below 9e-4 except Timbavati at 160 bins (−0.9%); particle water below 0.6% except Timbavati at 160 bins (−0.9%); end-of-run NH₃, HNO₃ and HCl within 0.2%; O₃ and NO₂ unchanged at output precision. For scale, Timbavati with 80 instead of 160 bins (original solver) changes chloride by 14% and water by 13%.

| Patch | What | Evidence |
|---|---|---|
| 29 | Deviation cache in the worst-first loop of `FindElectrolyteEquilibriumCore`. Each pass used to re-evaluate `|1 − ratio|` for every major reaction (about 1,550 ratio calls per solve in Timbavati at 160 bins, half of the solver time). A reaction's deviation depends only on the electrolyte, cation, anion and letovicite-anion amounts listed in `ECRDepIdx` and on water (activity coefficients, water activity and temperature are fixed in this loop), so the loop now keeps each deviation and re-evaluates it only when one of those amounts changed bit for bit since the previous pass (any change in water invalidates all). On by default; `ASP_FEE_DEVCACHE=0` turns it off. Harness additions: counters of outer passes, single-reaction solves and ratio evaluations per call, a histogram of passes per call, `ASP_FEE_REPLAY_REPEAT=<n>` and `ASP_FEE_REPLAY_OUT=<file>` (replayed results, for bitwise comparison). | Byte-identical to patch 28 (Williams, Otavi, Alaska, Timbavati at 40 bins, 109 files each; Timbavati at 160 bins, 229 files) and on every replayed record. Replay, 10 repeats: solver CPU 2.47 to 2.34 s (Timbavati), 1.67 to 1.04 s (Alaska). Full 160-bin runs, two each: Alaska 107.2 / 105.5 s to 100.4 / 96.2 s; Timbavati 99.2 / 109.8 s to 98.7 / 99.5 s (within noise). |

Outer-loop alternatives tried on the harness and not adopted (Timbavati and Alaska dumps, 160 bins): Gauss–Seidel sweeps (every out-of-tolerance reaction once per sweep) cycle among the potassium and ammonium salts in Timbavati and run to the sweep cap (mean difference 7e-2 of the dissolved moles against 2.6e-5); one opening sweep instead of three saves 19% on Timbavati but raises the mean difference to 1.7e-3; demoting a reaction as soon as equilibrating it changes nothing gives no speed-up and a mean difference of 2.2e-3. In Timbavati every solve takes 11–60 worst-first passes (none reaches the 100-pass cap) and about 171 single-reaction solves, 93 of them in the three opening sweeps.

| Patch | What | Evidence |
|---|---|---|
| 30 | Coupled Newton solve of the electrolyte equilibrium (`CoupledElectrolyteNewton`), replacing the opening sweeps and worst-first loop of `FindElectrolyteEquilibriumCore` after its first pass (not on the first equilibration). All movable reactions are solved together for ln(ratio) = 0 with activity coefficients, water activity and temperature fixed: analytic Jacobian, row-scaled dense elimination, amounts updated as n + ν·Δξ (mass and charge exact), steps kept positive; solids follow complementarity (active while present or supersaturated; a step that would make one negative stops at zero; an absent solid the step would dissolve is left out). Special cases (infinitely dissociating electrolytes, solids above DRH, neutralisation) stay with `EquilibrateAqueousDissociationReaction`, which runs a validation pass after each Newton cycle; on failure the particle is restored and the worst-first loop runs. **Solid-circuit rule:** `InitSolidCircuits` finds the minimal linearly dependent sets of 2–4 solids (36 in the baseline deck); pure solids have unit activity, so their metathesis has K = exp(Σ cᵢ ln Kᵢ(T)) and one side is unstable (e.g. NH₄NO₃ + K₂SO₄ against KNO₃ + (NH₄)₂SO₄, K = 0.39 at 288–298 K). An unstable side may not be completed by precipitation (coupled solve and single-reaction solver), an unstable side found present is dissolved lowest-DRH first (`ResolveUnstableSets`), and on a singular Jacobian active solids are dropped lowest-DRH first, accepted if the dropped solid ends undersaturated. A circuit applies only when all its solids are below their DRH. Selected by `AqOuterSolver` (optional `ASPInputDeck.in` line after `AqInnerSolver`; **0, off, in this patch**), overridden by `ASP_FEE_OUTER=newton` or `worst`; the rule goes with it (`ASP_FEE_PAIRRULE` overrides). Harness: `ASP_FEE_DUMP_ALL` (also calls that update activity coefficients, judged with the coefficients the solve used), `EqViolation` and `ViolationDiag` (`ASP_FEE_VDIAG`), counts of unstable solid sets, and `ASP_FEE_EXITFIX` (worst-first loop without its early exit on the last reaction's deviation). Developed on `global-newton-test`. | Byte-identical to patch 29 with the default (Williams, Otavi, Alaska, Timbavati at 40 bins, 109 files each; Timbavati at 160 bins, 229 files). |
| 31 | `AqOuterSolver` defaults to 1: the coupled solve with the solid-circuit rule is the electrolyte solver; the worst-first loop remains as `AqOuterSolver = 0` (or `ASP_FEE_OUTER=worst`). **Not byte-identical.** | Default output equals patch 30 with `ASP_FEE_OUTER=newton` bit for bit (four cases at 40 bins). See the tables below. Accepted by M. Alvarado. |

Equilibrium quality, replay at 160 bins (records violating the equilibrium conditions by more than twice the tolerance; water equilibrium and the special cases excluded; solver CPU):

| Case | Records | Worst-first (patch 29) | Coupled (patch 31) | Solver CPU, worst-first → coupled |
|---|---|---|---|---|
| Williams | 1,166 | 1,166 | 0 | 0.089 → 0.084 s |
| Otavi | 1,218 | 1,195 | 0 | 0.137 → 0.090 s |
| Alaska | 5,254 | 192 | 0 | 0.122 → 0.096 s |
| Timbavati | 5,262 | 3,762 | 0 | 0.248 → 0.156 s |

No coupled call falls back to the worst-first loop, and no result holds a complete unstable solid set (the worst-first results hold one in 661, 1,127, 2 and 3,529 records). Every worst-first violation is moved, and almost always fixed, by one single-reaction solve; they come from the loop's early exit on the deviation of the last reaction evaluated (`ASP_FEE_EXITFIX=1` removes it; that fixes Alaska at no cost but makes 67% of Timbavati calls run to the 100-pass cap).

Full 5 h runs, CPU time, one each (worst-first → coupled): 160 bins Williams 45.7 → 41.8 s, Otavi 48.8 → 41.9 s, Alaska 103.8 → 97.7 s, Timbavati 105.7 → 86.0 s; 40 bins 13.2 → 12.7, 13.2 → 11.3, 29.2 → 25.7, 30.1 → 23.2 s. Output change: OA, NH₄, NO₃ and SO₄ totals within 2e-5 (Otavi), 2e-3 (Williams), 7e-3 (Timbavati) and 1.4% (Alaska) at every 10-minute output; particle number −0.8% in Alaska, otherwise ≤ 0.06%; O₃ and NO₂ unchanged at output precision; HNO₃ and HCl within 0.12%; NH₃ within 0.24% except Timbavati (−1.3% / −1.1% at 160 / 40 bins). Timbavati particle chloride and water at 290 min change by +7.5% / +6.7% (160 bins) and −3.3% / −3.6% (40 bins); over the run the differences alternate in sign between 10-minute outputs (chloride, by then about 0.001 µg m⁻³, within ±4% from 80 to 280 min), and particle water is +43% at 20 min only, during early condensation.
| 32 | Fixes found with the ISORROPIA II ThermoTest comparison (seven bulk cases, 273.15–313.15 K, RH 0.90 → 0.20). (1) **Output:** `AerosolOutput` (ThermoTest output) declared but never set the letovicite index, so letovicite was missing from `TotNH4` and `TotSO4` (up to 99% of the sulfate at the letovicite composition, 313 K). Output only; present since before patch 0. (2) **Coupled-solve regression:** at NH₄:SO₄ = 1.5 below RH 0.40 the circuit (NH₄)₂SO₄ + NH₄HSO₄ → letovicite applied, but with ASP's per-electrolyte activity coefficients letovicite was undersaturated, so (NH₄)₂SO₄ was refused, letovicite did not form, about 1 µg m⁻³ of water stayed and up to 32% of the ammonium left as NH₃ at 313 K. A circuit now applies only if every solid on its stable side is present or supersaturated (`StableSideAvailable`); otherwise the unstable side may form, as with the worst-first loop. (3) **Crash:** the neutralisation step of `EquilibrateAqueousDissociationReaction` could consume more water than a nearly dry particle holds (negative molality, stop); it is now limited to 0.99 of the water. The pre-patch-0 code stopped this way at the letovicite composition at 6 of 9 temperatures. (4) The Newton single-reaction solve restores its last valid state if its bracket collapses without one. | ThermoTest: the letovicite case now dries out below RH 0.36 and keeps 0.90–1.00 of the ammonium, as the worst-first loop (0.906–1.00) and ISORROPIA (0.899–0.991); the other six cases are unchanged from patch 31; all 63 case–temperature runs complete and conserve every species to 0.2%. Replay at 160 bins: no unstable solid sets; violations Williams 1 (1 fallback), Otavi, Alaska, Timbavati 0. Full 40-bin runs against patch 31: Otavi and Alaska byte-identical; Williams within 7e-5; Timbavati OA, NH₄, NO₃, SO₄ within 1.4e-3, chloride within 3.7%, particle water −25% at 20 min, within ±2.5% after 60 min. |
| 33 | ThermoTest run mode can scan RH upward: `ASP_THERMO_SCAN=up` runs 0.20 → 0.90 (deliquescence branch; use a deck with initial RH 0.20); the default stays 0.90 → 0.20. | Down-scan ThermoTest output byte-identical to patch 32 (all 63 case–temperature runs, four output files each). The up scans complete and conserve every species to 0.2%. |
| 34 | Fewer coupled-solve fallbacks (singular Jacobians), found with the ThermoTest upward scans. (1) The water amount is held fixed in the coupled solve (water stoichiometry removed from the extents): a nearly dry particle cannot supply the water that NH₃ + H₂O ⇌ NH₄⁺ + OH⁻ or the water equilibrium consume, and the steps drove water to ~1e-28 and the molalities to ~1e14; the water amount is set by the water equilibrium. (2) Amounts below 1e-200 count as zero when choosing the active reactions (denormal NH₃ amounts gave infinite Jacobian entries). (3) On a singular Jacobian the water equilibrium is first left to the validation pass (with H⁺ or OH⁻ vanishingly small it is numerically dependent on the reactions sharing that ion). (4) If no single DRH-ordered solid drop works, solids are dropped cumulatively, lowest DRH first, until the solve converges with every dropped solid undersaturated. | Replay of every ThermoTest solve (7 cases × 9 T): fallbacks down scan 85 → 10 (letovicite composition), up scan 855 → 46 (42 at the letovicite composition, where the Newton and validation passes do not settle in 4 cycles, and 4 in Urban); records violating by > 2× tolerance down 129 → 13, up 1,014 → 1,178 (1,004 in K₂SO₄ + NH₄NO₃ at 278–283 K: (NH₄)₂SO₄ left supersaturated by the single-reaction solver's "not enough SO₄²⁻" rule, not moved by it). Lagrangian replays at 160 bins: Williams 1 fallback (a concentrated particle with six active solids), Otavi, Alaska, Timbavati 0. Full 40-bin runs vs patch 32: Otavi, Alaska and Timbavati totals unchanged (≤ 2e-11), Williams within 5e-5. ThermoTest: down scan within 0.039 in aerosol fraction (letovicite, 313 K); up scan changes at the deliquescence steps (single RH points up to 0.96, e.g. Urban 283 K now deliquesces at RH 0.80, as ISORROPIA); mean deliquescence RH changes by at most 0.04. |
| 35 | Dry particles in the coupled solve. Every remaining ThermoTest fallback was a dry particle (water ~1e-28 against ~1e-13 of electrolyte; molalities 1e14 and more), where neither the Newton solve nor the single-reaction solver settles. When the Newton solve fails on a particle whose water is below 1e-3 of its electrolyte (moles), `DrySolve` puts the electrolyte into solids by minimising the Gibbs energy of the solid assemblage: with pure solids at unit activity this is the linear programme min Σ nₛ ln Kₛ(T) over the reaction extents (mass and charge conserved; ions, dissolved electrolytes and solids above their DRH cost 1e4 per mole), solved by a dense simplex with Bland's rule; a 1e-9 trace of the incoming solution is kept for the water routines. Same pure-solid thermodynamics as the solid-circuit rule. Harness: dry results counted separately and violations counted for wet results only. Tried and not adopted: running `DrySolve` on every dry particle before the Newton solve (Timbavati 40 bins: nitrate up to 32%, chloride 48%, particle water 81% different), and treating the stable side of a circuit as available in dry particles (more fallbacks and wet violations at the letovicite composition). | Replay of every ThermoTest solve: fallbacks down 10 → 0, up 46 → 0; wet results violating by > 2× tolerance down 0, up 29 (worst-first loop: 1,361 and 1,472). All 126 scans complete. Lagrangian: 40-bin runs byte-identical to patch 34 (the dry solve is not reached); the Williams 160-bin replay keeps its one fallback (a concentrated wet particle with six active solids). ThermoTest outputs change only where the fallbacks were (4 downward and 10 upward scans: letovicite composition and Urban 313 K). |
| 36 | Last coupled-solve fallback and the wet violations. (1) The 29 wet violations on the ThermoTest upward scans were all HSO₄⁻ ⇌ H⁺ + SO₄²⁻ with HSO₄⁻ 10¹⁰–10¹⁶⁷ times below equilibrium: reactions with ratio ≥ 1e11 were left out of the Newton active set, and the validation pass moved HSO₄⁻ by less than its 1e-28 change threshold. The cut-off is now 1e300; before the Newton solve, a dissolved electrolyte that is below equilibrium is set directly to ratio × amount when that takes at most 1e-3 of its ions (`PreAssociate`); up to 100 Newton iterations (was 40). (2) The Williams fallback (a concentrated wet particle, water ≈ 1.5 × the dissolved moles, with KNO₃, KHSO₄, NH₄Cl, NH₄NO₃, (NH₄)₂SO₄ and letovicite active): when the single and cumulative DRH-ordered drops fail, every subset of the active solids is tried as the set left out (smallest first; at most 2¹⁰ solves) and the first that converges with every left-out solid undersaturated is accepted. | Replay of every ThermoTest solve (7 cases × 9 T, both scans): no fallbacks and no wet results violating by > 2× tolerance (worst-first loop: 1,353 / 1,490); all 126 scans complete. Lagrangian 160-bin replays: no fallbacks, no wet violations in Williams, Otavi, Alaska and Timbavati. Full 40-bin runs vs patch 35: Williams within 6e-5, Otavi within 8e-7, Alaska within 1.9e-3 (particle water; ions within 4.3e-4), Timbavati within 2e-6. ThermoTest: 28 scans per direction change by more than 1e-4 in aerosol fraction, largest 0.05 (letovicite NH₄, 298 K, down) and 0.08 (Urban NH₄, 313 K, up). Letovicite NH₃ gas–aerosol split against ISORROPIA stable: mean difference 0.001–0.025 by temperature, largest 0.11 (313 K). |

### 7.6 Speed against the pre-Claude code after patch 36 (Lagrangian parcel, 5 h, CPU time)

Baseline: `34b7ec7` (last commit before any Claude change) plus the `BinEdges` build guard (`09afe88`) and the `IonReset` fix (`5b4647e`), built with the makefile of that period; current: `asp_v2.2.2` patch 36 (`5ae6004`). Runs alternate on the same machine.

| Case, bins | Baseline | Patch 36 | Faster by |
|---|---|---|---|
| Timbavati, 160 | 519.1, 518.9 s | 84.8, 85.6 s | 6.1× |
| Williams, 160 | 196.2, 200.9 s | 44.1, 45.3 s | 4.4× |
| Otavi, 160 | 194.8, 189.9 s | 40.8, 40.9 s | 4.7× |
| Alaska, 160 | stops at 4 min (coagulation volume check) | 91.6 s | – |
| Alaska, 80 | 304.9 s | 47.4 s | 6.4× |
| Alaska, 40 | 153.2 s | 25.4 s | 6.0× |

The baseline Alaska run at 160 bins stops at t = 4 min with "Summed Total Volume concentration is lost during step" in `StepSectionalCoagulationJacobson` (the bin-ordering issue in `ASP_Pending_Patch_Notes.md`, section 1); with a 30 s coagulation step both builds stop at 1.5 min, so that workaround does not apply here. Williams and Otavi gain less because their electrolyte solves are fewer and simpler.

| Patch | What | Evidence |
|---|---|---|
| 37 | Coagulation volume check (pending-patch notes, section 1). The semi-implicit scheme in `StepSectionalCoagulationJacobson` needs the occupied sections in ascending order of their actual single-particle volume: the volume partition brackets vᵢ + vⱼ between neighbouring sections, and gains to section k use the already-updated smaller sections. Condensation can move a section's particles past its neighbour, and the two-sided volume check then stops the run. Now, after `SortAerosolAtGridPointForCoagulation`, if any occupied section is smaller than the one before it, `RegridAerosol` runs first (repeated up to 10 times, since it moves particles at most one section per call). Chosen over sorting by actual volume inside the step (M. Alvarado). | Alaska at 160 bins now completes with the 60 s and the 30 s coagulation step (patch 36 stopped at 1.5 min with 30 s; the pre-Claude code stops at 4 min with 60 s). Against patch 36: Otavi (40 and 160 bins), Alaska and Timbavati at 40 bins byte-identical (sections already in order); Alaska 160 totals within 1.6e-6; Williams 40 within 1.1e-6; Williams 160 within 0.56% at 290 min; Timbavati 160 OA, NH₄, NO₃, SO₄ within 3.5e-4, chloride 3.0% (190 min) and particle water 3.8% (20 min), largest change in one section's number 1.9% of the total; end-of-run NH₃ −4e-4. CPU time within run-to-run noise (Timbavati 160 88.2 s, Alaska 160 95.3 s). |
