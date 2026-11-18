# ASP information flow: Lagrangian parcel model (gas + aerosol)

How control and data move through `common_ASP` when `DevelopmentDriver.f90` (`PROGRAM ASP`) runs the Lagrangian
parcel model with aerosols. The history of the changes behind this design, with timings and their effect on results, is
in `CHANGELOG.md`.

## 1. Scope

The run mode is selected in `ASPInputDeck.in` by:

- "Lagrangian Parcel Model (1 = yes)" = 1, all other run-mode lines 0, and
- "Include aerosols?" = 1, so the driver takes the gas + aerosol branch (`ELSE IF (IncludeAerosols)` inside the
  `LagrangianParcelModel` block). The gas-only branch is summarised in section 3.9.

The parcel is one well-mixed air parcel diluted with a background atmosphere. Every step it integrates the gas-phase
chemistry (stiff sparse ODE solve), mixes gases and particles with the background, and at the condensation and
coagulation steps runs aerosol condensation/dissolution and coagulation on a sectional (binned) size distribution.
Nothing is held static: all processes interact, so the cost is spread over the gas solver and the aerosol
thermodynamics (section 6).

## 2. High-level flow

```mermaid
flowchart TD
    A[PROGRAM ASP] --> B[SetDomainSize, ReadMainInputDeck]
    B --> C{LagrangianParcelModel and IncludeAerosols}
    C --> D[SetChemistryParams + ReadPhotolysisRates\nfull gas, aqueous and ODE-set tables]
    D --> D1[RecalculateReactionRates,\nRecalculatePhotoRates at UTCStartTime]
    D1 --> E[SetAllDissolution]
    E --> F[InitializeParticlesSectional, RegridAerosol\nsectional population from AerosolModes.in]
    F --> F1[ReadBoundaryDistributionsOrganic +\nPopulateParticlesSectionsRightAway TRUE + SortBoundaryPart\nbackground aerosol from EnvAerosolModes.in]
    F1 --> G[EquilibrateInternallyAtGridPoint water, RegridAerosol]
    G --> H[Hetero rates from the size distribution,\nopen output files, write t = 0 output]
    H --> L{{Main loop: 600 + (RunLength-10)*60/MixStep iterations}}
    L --> L1[RecalculateReactionRates, RecalculatePhotoRates,\nhetero number and radii, RecalculateHeteroRates]
    L1 --> L2[StepGasChemistry\nDLSODES, ChemTimestep]
    L2 --> L3{t multiple of CondTimeStep, t not 0}
    L3 -- yes --> L4[RegridAerosol, FindAqElectrolyteEquilibriumForGridPoint,\nStepCondensationAll]
    L3 -- no --> L5
    L4 --> L5{t multiple of CoagTimeStep, t not 0}
    L5 -- yes --> L6[StepSectionalCoagulationJacobson]
    L5 -- no --> L7
    L6 --> L7[Dilution and deposition of gases: T_mix = Slope*t + b]
    L7 --> L8[UpdateChemicalConcentrations]
    L8 --> L9[SortAerosolAtGridPointForCoagulation; per bin:\nTerminalVelocity, dilute + deposit number and mass,\nRecalculateRadius]
    L9 --> L10{t multiple of 600 s}
    L10 -- yes --> L11[ShellRefIndAndRad per bin, AerosolOptProp,\nwrite AeroOptCIS<min> file]
    L10 -- no --> L12
    L11 --> L12{t multiple of 60 s}
    L12 -- yes --> L13[GasOutput, GasEnhanceOutput,\nAerosolOutput for every bin]
    L12 -- no --> L14
    L13 --> L14[SizeDistOutput at multiples of 600 s, and at 540 s and 2820 s]
    L14 --> L15[t = t + MixTimestep, ut = ut + MixTimestep/3600]
    L15 --> L
```

## 3. Step by step

### 3.1 Program entry and branch selection

`ReadMainInputDeck` reads `ASPInputDeck.in`, including the optional solver settings at its end
(`ReadOptionalDeckSettings`, see `README.md`). The driver's `IF / ELSE IF` chain enters the `LagrangianParcelModel`
block and, with aerosols on, its aerosol branch. The mixing parameter `b_mix` (deck section 3.b) sets the dilution time
scale.

### 3.2 Chemistry set-up (full level)

`SetChemistryParams` builds all gas, aqueous and ODE-set tables and the gas Jacobian (`SetChemSetupLevel` stays at its
default 0 in this mode). `MakeODESet` parses the reactions; the rate-law flags are listed in the table above its
gas-phase section in `Chemistry/ChemistryParametersInitRoutines.h`. `ReadPhotolysisRates` reads `PhotoRates.in`.
`RecalculateReactionRates(IFGAS, IFAQ)` evaluates the temperature-dependent rate constants
(`ResetGasPhaseReactionRates`) and `RecalculatePhotoRates(ut, t, .TRUE.)` the photolysis rates at `UTCStartTime`.

### 3.3 Dissolution data and the particle population

- `SetAllDissolution`: inorganic and organic (hydrophobic and hydrophilic) dissolution data.
- `InitializeParticlesSectional` builds the sectional population from `AerosolModes.in` (log-normal modes, section
  cutoffs); `RegridAerosol` follows.
- `ReadBoundaryDistributionsOrganic`, `PopulateParticlesSectionsRightAway(.TRUE.)` and `SortBoundaryPart` build the
  **background (boundary) aerosol** from `EnvAerosolModes.in`, the aerosol the parcel mixes with.
- `EquilibrateInternallyAtGridPoint(EquilibrateWater=.TRUE., WaterOpenSystem=.TRUE.)` and `RegridAerosol` bring every
  bin to water equilibrium. The bins are counted (`NumBins`, `NumAqChems`, `NumOrgChems`, `NumAqOrgChems`) and the work
  arrays (`NumConc`, `InitNumConc`, `EnvNumConc`, `HeteroNumberConc`, `HeteroRadii`) allocated.
- The heterogeneous-reaction number concentrations and radii are taken from the sections (`EffectiveRadius/100`, cm
  to m) and `RecalculateHeteroRates(..., .TRUE.)` sets the heterogeneous rates (`ResetHeteroRates`).

### 3.4 Output at t = 0

`GasOut1/2/3.txt`, `EnhanceRatios.txt`, one `AerosolN.txt` per bin, and `SizeDist0min.txt`.

### 3.5 Main time loop

The loop runs `INT(600 + (RunLength-10)*60/MixStep)` iterations. For the first 600 s the steps are 1 s
(`ChemTimestep = MixTimestep = 1`); after that they are the deck values (`ChemStep`, `MixStep`). Each iteration:

1. **Rates.** `GetTemp`, `RecalculateReactionRates`, `RecalculatePhotoRates(ut, t, .FALSE.)`; heterogeneous number
   concentrations and radii from the particle list, then `RecalculateHeteroRates`. The heterogeneous rate flags (1-body
   3, 18–21 and 2-body 6) are set in `ResetHeteroRates`; flag 18 reads per-particle acidity through `AcidStateHook`
   (filled by the Condensation module's `FillAcidState`), flags 19 and 20 the bulk fine-particle state through
   `BulkStateHook` (`FillBulkState`).
2. **Gas chemistry.** `StepGasChemistry(ChemTimestep, INT(MixTimestep/ChemTimestep))` (`Chemistry/EvolveGasChemistry.h`)
   integrates the gas ODEs with the sparse stiff solver `DLSODES` (`istate = 1` on every call, `mf = 21`,
   `iwork(6) = 2000`); `Y`, `rwork` (at least 5,000,000 doubles) and `iwork` are allocated on every call. The right-hand
   side and Jacobian come from `GasODEEvaluator` and `GasJacobianEvaluator`, which by default loop over flat copies of
   the ODE set and Jacobian (`BuildGasFlat`; `GasEvalLinked = 1` walks the linked lists instead, with the same result).
   `DPREP` keeps the minimum-degree orderings (`ODRV`) of the last 8 sparsity structures and reuses them.
3. **Condensation** (every `CondTimeStep`, not at t = 0). `RegridAerosol`,
   `FindAqElectrolyteEquilibriumForGridPoint(UpdateThermo=.TRUE.)`, then `StepCondensationAll(CondTimeStep, 1)`: LSODES
   for H₂O and H₂SO₄ plus the inorganic equilibrium of every bin, the organic condensation steps, the particle-phase
   organic reactions of `OrgPhaseReactions.in` if present, and `RegridAerosol` (details in
   `ASP_CondTest_InformationFlow.md`, section 3.7).
4. **Coagulation** (every `CoagTimeStep`, not at t = 0). `StepSectionalCoagulationJacobson(CoagTimeStep)` (details in
   `ASP_CoagTest_InformationFlow.md`). It first regrids if the occupied sections are out of size order, and stops the
   run with "Total Volume concentration is lost during step" if the total volume changes by more than 1e-4.
5. **Gas dilution and deposition.** `T_mix = Slope*t + b` with `Slope = 2.0`, `b = b_mix`. For every gas species,
   `Chem(q) = (Chem(q) + Background(q)*dt/T_mix) / (1 + dt/T_mix + dt*DepVel(q)/InvHgt)`;
   `UpdateChemicalConcentrations` stores the result.
6. **Particle dilution and deposition.** `SortAerosolAtGridPointForCoagulation`, then for each bin, with the matching
   background bin: `TerminalVelocity`; the new number `N = (N0 + Nenv*dt/T_mix) / (1 + dt*(1/T_mix + TermV/InvHgt))`;
   the same formula for the mass of every aqueous, organic and aqueous-organic species (divided by the new number to
   give moles per particle); `Temperature = TEMP`; `RecalculateRadius`.
7. **Optical properties** (every 600 s). `ShellRefIndAndRad` for each bin, `AerosolOptProp` (Mie code: `DMiLay` for the
   core-shell case, MIEV0 for the homogeneous spheres), and `AeroOptCIS<N>min.txt` (451 wavelengths).
8. **Output** (every 60 s). `GasOutput`, `GasEnhanceOutput`, `AerosolOutput` for each bin; `SizeDistOutput` at
   multiples of 600 s and at 540 s and 2820 s.
9. **Time update.** `t = t + MixTimestep`, `ut = ut + MixTimestep/3600`; `ChemTimestep` and `MixTimestep` are reset
   to the deck values after 600 s. The loop also writes several progress lines per iteration to standard output
   ("In Loop", "Photo Okay", ..., "Before AerosolOutput, Bin #"); a 300-minute Williams run writes about 207,000 lines
   (9 MB).

### 3.6 Electrolyte equilibrium

`FindElectrolyteEquilibrium` (a wrapper around `FindElectrolyteEquilibriumCore`, `LagrangianParticles/Thermodynamics.h`)
solves the aqueous ionic equilibrium of one particle. After a first pass of single-reaction solves
(`EquilibrateAqueousDissociationReaction`, by safeguarded Newton on ln(ratio) with `AqInnerSolver = 1`, or bisection
with 0), the default outer solver (`AqOuterSolver = 1`) is the coupled Newton solve `CoupledElectrolyteNewton`, which
solves all movable reactions together and follows the solid-circuit rule (`InitSolidCircuits`, `ResolveUnstableSets`,
`StableSideAvailable`); dry particles go to `DrySolve`. With `AqOuterSolver = 0` the worst-first loop is used, with an
optional deviation cache (`AqDevCache`). `EquilibriumConstantsRatio` (with `MeanActivity`, `EqCoeff`, `PowMemo` and a
by-value memo) gives each reaction's distance from equilibrium; the water activity uses the tabulated Gibbs–Duhem
integral (`GibbsDuhemIntegralTabulated`; `GDTableMode = 0` integrates with `DGAUS8` instead).

### 3.7 Key routines

| Stage | Routines | Files |
|---|---|---|
| Rate constants | `RecalculateReactionRates`, `ResetGasPhaseReactionRates`, `RecalculatePhotoRates`, `RecalculateHeteroRates`, `ResetHeteroRates` | `Chemistry/Chemistry.f90`, `EulerianEnvironment/GasPhaseChemistry.h` |
| Gas ODE solve | `StepGasChemistry`, `GasODEEvaluator`, `GasJacobianEvaluator`, `BuildGasFlat`, `DLSODES` | `Chemistry/EvolveGasChemistry.h`, `ComputationalSubroutines/opkdmain.f`, `opkda1.f` |
| Condensation | `StepCondensationAll`, `StepCondensation`, `EquilibrateGridPoint` | `Condensation/CondensationIntegrator.h`, `CondensationRelatedFunctions.h` |
| Inorganic equilibrium | `FindElectrolyteEquilibrium`, `CoupledElectrolyteNewton`, `EquilibrateAqueousDissociationReaction`, `EquilibriumConstantsRatio`, `MeanActivity`, `KusikMeissner` | `LagrangianParticles/Thermodynamics.h` |
| Particle properties | `RecalculateRadius`, `ParticleDensity`, `SurfaceTension` | `LagrangianParticles/ParticleAttributes.h` |
| Coagulation | `StepSectionalCoagulationJacobson` | `Coagulation/Coagulation.f90` |
| Regrid and sort | `RegridAerosol`, `SortAerosolAtGridPointForCoagulation` | `LagrangianParticles/SortAerosol.h` |
| Optics | `ShellRefIndAndRad`, `AerosolOptProp`, `MieHomogeneousSphere`, `DMiLay` | `LagrangianParticles/ParticleAttributes.h`, `MIEV0.f`, `DMiLay.f` |
| Output | `GasOutput`, `GasEnhanceOutput`, `AerosolOutput`, `SizeDistOutput` | `OutputRoutines/OutputRoutines.f90` |

### 3.8 Settings that select code paths

All are optional `Name = value` lines at the end of `ASPInputDeck.in` (`README.md`): `AqOuterSolver`, `AqInnerSolver`,
`AqPairRule`, `AqExitFix`, `AqDevCache` (electrolyte solver), `GDTableMode`, `GDTableCells` (Gibbs–Duhem integral),
`GasEvalLinked` (gas evaluator). `AqDevCache` and `GasEvalLinked` do not change results.

### 3.9 Gas-only branch

With "Include aerosols?" = 0 the same loop runs only the rates, `StepGasChemistry`, gas dilution and gas output, with a
single heterogeneous bin from the deck (`hetero_num`, `hetero_rad`).

## 4. Inputs and outputs

**Inputs** (in `InorganicInputDecks`):

| File | Used for |
|---|---|
| `ASPInputDeck.in` | run-mode flags, T, P, RH, run length, the four time steps, inversion height, `b_mix`, photolysis settings; optional solver settings at the end |
| `AerosolModes.in` | modes (number, geometric mean diameter, sigma, composition), number of sections and cutoffs |
| `EnvAerosolModes.in` | background aerosol the parcel mixes with |
| `GasPhaseChems.in`, `GasChemicalMechanism.in`, `PhotoRates.in` | gas species (initial and background values), reactions, photolysis |
| `AqPhaseChems.in`, `AqPhaseIons.in`, `AqEquilibriumReactions.in`, `AqChemicalMechanism.in`, `AerosolThermoKusikMeissner.in` | aqueous species, ions, equilibria and activity-coefficient data |
| `Dissolution.in`, `OrganicDissolution.in`, `HydrophilicOrganicDissolution.in`, `OrgPhaseChems.in`, `HydrophilicOrgChems.in` | gas–aerosol partitioning and organics |
| `OrgPhaseReactions.in` (optional) | first-order particle-phase organic reactions |
| `NitrateGroups.in` (optional) | species summed into the PN and AN enhancement ratios of `EnhanceRatios.txt` |

**Outputs** (`OutputFiles`; 109 files for a 300-minute 40-bin run): `GasOut1/2/3.txt` and `EnhanceRatios.txt` (every
60 s), `Aerosol1..N.txt` (one per bin, every 60 s), 32 `SizeDist<N>min.txt`, 30 `AeroOptCIS<N>min.txt`,
`GasPhaseJacobian.txt`, `GasPhaseODESet.txt` and `Transcript.TXT`.

## 5. Behaviour to know when comparing runs

- **Threshold sensitivity.** Parcel runs can sit near a threshold of the equilibrium solver: in Williams a relative
  change of 1e-4 in the input RH changes 60-minute size-distribution sums by about 10% (BC −9.8%, NH₄ −12%) through a
  one-step composition event in one bin. Differences between builds in this mode must therefore be judged on totals
  and over time, or byte for byte, not by small tolerances on individual bins.
- **Number of bins.** Results depend on the bin count (Timbavati chloride changes by 14% between 80 and 160 bins).
- **Coagulation step.** The RemConMode1 deck uses a 30 s coagulation step: with 60 s its sections fell out of size
  order after the first condensation step and the volume check stopped the run at t = 60 s (volume ratio 1.000121).
  `StepSectionalCoagulationJacobson` now regrids first in that situation; the deck has not been switched back.

## 6. Where the time goes

Measured with callgrind on the October 2026 code (details and later changes in `CHANGELOG.md`):

- In the first 15 minutes (1 s steps for 10 minutes) the gas solver dominated (`StepGasChemistry` 46%, of which the
  sparse ordering `ODRV` 28%, now cached) and condensation took 29%.
- Over longer runs condensation dominates: about 60–70% in `StepCondensationAll`, mostly the electrolyte solver
  (`EquilibriumConstantsRatio`, `MeanActivity`, `KusikMeissner`, libm `pow`), with `DLSODES` 10–17%, `SurfaceTension`
  about 9% and `memset` (zeroing in `DPRJS`) about 4%.
- Wall time grows about linearly with the number of bins (N^0.8–0.9 from 40 to 160 bins).
- Small: optics (about 3.5%), coagulation (0.4%), the per-call `rwork` allocation and the standard-output writes.

Open candidates are listed in `README.md` (open items).
