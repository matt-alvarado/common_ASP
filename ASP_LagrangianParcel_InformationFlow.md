# ASP Model — Information Flow for `DevelopmentDriver.f90` with the Lagrangian Parcel Model (gas + aerosol)

Source basis: `common_ASP` branch `asp_v_2.2.1` (commit `fdf2da9`, patches 00–19 applied). Line numbers refer to `DevelopmentDriver.f90` in that tree.

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
| RemConMode1 | aborts after about 0.7 s (SLATEC `DGAUS8`) | — | — | — | stops after about 1.3 s: coagulation volume-loss error (section 6) |

So the baseline for Williams is about 140 s; the current branch is 3.3 times faster. Timbavati and RemConMode1 have no usable baseline because the original aborts on both (the `IonReset` bug, patch 16). The intermediate builds show where the Williams gain came from: N1–N5 save about 27%, patches N6–N11 and the Mie, ODE-set and run-mode changes cut a further 54%, and patches 11–19 about 13% more.

### 5.3 Output comparison with the baseline

- Patches 00 to N5 (bc80201): all 108 numerical output files on Williams are **byte-identical** to the original build.
- From N6–N11 (`a56950d`) on, output is no longer identical. At 10 min, total particle number, BC, water and K match the original to 1e-6 (SizeDist files), but by 60 min the sums differ by up to 10% in some columns (BC −10%, NH4 −12%) and by 290 min by up to 8%.
- **Cause of the drift (tested).** Building `a56950d` with only the N11 tabulated Gibbs–Duhem integral switched off (always DGAUS8) gives output byte-identical to the original on all 108 files, so N6–N10 are exact and N11 alone starts the drift. N11 changes the activity integral by about 2e-4, which shows up as 1e-6 to 4e-4 differences in individual bins at t = 0.
- **Why it grows (tested).** The same drift appears in the *original* build when the input RH is changed by only 1e-4 (relative; 0.7 to 0.70007, initial and background): same bin, same time, and the same 60-minute sums (BC −9.8%, NH4 −12%) as the N11 build. A change of 1e-6 instead gives differences of only about 4e-4 even at 290 min. In bin 30 the original shows a one-step composition event at t = 1500 s (solid (NH4)2SO4, K2SO4 and KNO3 jump, for example (NH4)2SO4 1.5e-22 against 1e-28 mol per particle) while particle number and radius are unchanged; after that K2SO4 stays about 70% above the perturbed runs. So the reference run sits near a threshold in the equilibrium solver, and any perturbation of about 1e-4 or more (N11, or a tiny RH change) tips it the other way. I did not identify which decision in the solver flips.
- Consequence: the drift is a sensitivity of the model, not evidence that N11 is wrong, but it means differences between builds in this mode cannot be judged by comparing files against the original after the first minutes. Gas and aerosol totals stay within a few percent of each other for 60 minutes, then individual bins differ at the order-one level where an event differs.
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
2. **RemConMode1 stops at the first coagulation step (t = 60 s)** with "Total Volume concentration is lost during step". It happens in a build with patches 00–10 plus the fixes 16 and 17 but without my coagulation patches 11–15, so the cause is not the coagulation optimizations. I have not found the cause. RemConMode1 runs fine in CoagTest and CondTest.
3. **The Williams output no longer matches the original build** after N11 and the OLDASP fix (section 5.3).
4. **Stdout volume.** About 207,000 lines (9 MB) are written for the 300-minute Williams run, mostly "Okay" and "Before AerosolOutput, Bin #" progress messages. I did not measure their cost.
5. **Allocation.** `StepGasChemistry` allocates `rwork` of at least 5,000,000 doubles (40 MB) on each of its 605+ calls. The profile shows no large cost from this, but it is a candidate for a one-time allocation.
6. **LSODES options.** `mf = 21` with `istate = 1` on every call, while the code also loads `GasIa`/`GasJa` (sparsity structure) into `iwork`. I did not check which LSODES structure option `mf = 21` selects or whether the supplied structure is used; this matters for hypothesis 1.

## 7. Optimization hypotheses (from the profile above; none implemented or tested yet)

1. **Sparse set-up repeated on every `DLSODES` call (30% of the run).** `DPREP` (815 calls, about 19.5 M instructions each) computes the sparse structure and minimum-degree ordering again at the start of each call because `istate = 1`. The sparsity pattern of the mechanism never changes. Options: reuse the ordering across calls (needs a change inside the OPK routines; continuing with `istate = 2` instead would change the step-size history and so is probably **not** byte-identical), or choose the LSODES structure option that takes the supplied `GasIa`/`GasJa` (check the LSODES documentation first, see section 6 item 6). The ordering result is deterministic, so a cached copy should give identical numbers; this is the largest single saving and also the riskiest.
2. **Zeroing in the Jacobian path (7% of the run).** `GasJacobianEvaluator` and `DPRJS` clear full arrays 372,680 and 4,329 times. Clearing only the non-zero positions (`GasIa`/`GasJa` give them) would produce identical results.
3. **Condensation thermodynamics (29%).** The hot spots are the same as in CondTest after patches 18 and 19: `EquilibriumConstantsRatio`, `MeanActivity`, `KusikMeissner`, `pow`. Further CondTest work carries over directly; avoiding re-solving unchanged reactions in `FindElectrolyteEquilibrium` is the biggest item.
4. **Duplicate `RecalculateRadius` in the dilution loop (9% total, possibly half avoidable).** The first call in the loop (line 943) acts on a state that the previous iteration's second call (line 1005), or `RegridAerosol` after condensation and coagulation, has already recomputed. Like CoagTest patch 14, it should be skippable on the 590 iterations that have no condensation or coagulation step. This must be checked with a bitwise probe; temperature is set to `TEMP` only before the second call.
5. **`SurfaceTension` inside `RecalculateRadius` (11%).** This is the option A parked earlier ("until we test other run modes that call coagulation"). Parcel mode is such a run mode. Caution: results in CoagTest ran with unequilibrated particles, but here the particles are equilibrated, so any simplification must be tested on this mode.
6. **Per-call allocation in `StepGasChemistry`.** Allocate `Y`, `rwork`, `iwork` once. Small but free.
7. **Stdout writes.** Commenting out the per-iteration progress messages would help only if I/O is significant. Not measured; the thermodynamics warnings should be left alone, as previously agreed.
8. **Optics every 10 minutes (3.5% in 15 minutes).** Only 30 calls in a full run. Not worth attention unless the optical output frequency changes.
9. **Coagulation (0.4%).** Not worth attention in this mode.
