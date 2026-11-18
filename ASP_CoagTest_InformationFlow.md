# ASP information flow: coagulation test (`CoagTest`)

How control and data move through `common_ASP` when `DevelopmentDriver.f90` (`PROGRAM ASP`) runs with the deck flag
"Coagulation Test" set (`CoagTest = .TRUE.` in `ModelParameters`). The history of the changes behind this design, with
timings and their effect on results, is in `CHANGELOG.md`.

## 1. Scope

In this mode ASP runs **coagulation only**: no chemistry, condensation or transport. It builds a static sectional
(binned) aerosol population from `AerosolModes.in`, lets it coagulate for a fixed 10 hours in 1-minute steps, and writes
the number and radius of each bin every minute and a size-distribution snapshot every hour. Its purpose (per the comment
in the driver) is comparison with an analytic coagulation solution.

## 2. High-level flow

```mermaid
flowchart TD
    A[PROGRAM ASP] --> B[SetDomainSize, ReadMainInputDeck]
    B --> C{CoagTest ?}
    C -- yes --> D[SetChemSetupLevel 2 + SetChemistryParams\nspecies lists, equilibrium data only]
    D --> E[InitializeParticlesSectional]
    E --> E1[ReadDistributionsOrganic\nAerosolModes.in -> bins + modes]
    E1 --> E2[PopulateParticlesSectionsRightAway\none section per bin, lognormal modes summed]
    E2 --> F[SortAerosolAtGridPointForCoagulation + SpillBeans]
    F --> G[Count bins, open CoagTestOut.txt, write header]
    G --> H[RegridAerosol + SpillBeans]
    H --> I[SizeDistOutput t = 0]
    I --> J[Initial totals: number, volume, mass, K mass]
    J --> K{{Loop r = 1..600\nMixTimestep = CoagTimeStep = 60 s}}
    K --> K1[Write time + number and radius of each bin\nto CoagTestOut.txt]
    K1 --> K2[Recompute totals - diagnostic only]
    K2 --> K3[StepSectionalCoagulationJacobson 60 s]
    K3 --> K4[t = t + 60]
    K4 --> K5{t multiple of 3600 s ?}
    K5 -- yes --> K6[SizeDistOutput SizeDist NNNmin.txt]
    K5 -- no --> K
    K6 --> K
    K --> L[Final totals, DEALLOCATE, end of branch]
```

## 3. Step by step

### 3.1 Program entry and branch selection

`PROGRAM ASP` calls `SetDomainSize`, `SetFileHandleCounter(12)` and `ReadMainInputDeck()`. The deck sets the run-mode
flags, the environment (T, P, RH) and the coagulation settings of section 4. The driver's `IF / ELSE IF` chain enters
the `CoagTest` block and writes "Coagulation Test Selected." to the transcript.

### 3.2 Chemistry set-up (level 2)

`SetChemSetupLevel(2)`, then `SetChemistryParams`, which reads the gas, aqueous, organic and hydrophilic-organic species
definitions and the aqueous equilibrium reactions and sets the dissociation reactions; at level 2 it skips the gas- and
aqueous-phase ODE sets and Jacobians, which this mode never integrates. Results used later: species names, molecular
masses, densities (`AqMolecularMass`, `SolidSaltDensity`, `OrgDensity`, `AqOrgDensity`, ...), the equilibrium list
(`AqEquilibriaList`) and the counts (`HowManyAqChems`, `HowManyOrgChems`, ...). `SetAllDissolution` is not called.

### 3.3 Building the initial aerosol population

`InitializeParticlesSectional` (`LagrangianParticles/LagrangianAerosols.f90`):

1. Empties the particle list and sets the initial density (`SetInitialDens`).
2. `ReadDistributionsOrganic` parses `AerosolModes.in`: a Lagrangian/sectional flag, an input-type flag (relative molar
   abundance or mass burden), the sectional/Lagrangian size cutoff, the number of bins, the lowest-bin cutoff, then one
   line per aerosol mode (number concentration, size parameters, composition). It calls `SetSectionalParameters`
   (fills `BinEdges`, `HowManyBins`) and sets `ThermoBulkMode = .FALSE.`
3. `PopulateParticlesSectionsRightAway(.FALSE.)` creates one section per bin (`MakeSection(1, J)`: edges from
   `BinEdges`, number 0, composition and optical arrays allocated and zeroed), then adds each mode's lognormal number
   and mass share to each bin (bins with fewer than 1e-6 particles/cm³ are skipped) and fills `AqChems`, `OrgChems` and
   `AqOrgChems` (mol per particle). With the deck's `Monodisperse` flag it makes a single particle (`MakeMonodisperse`).

Result: a singly linked list `Particles%First -> ...` of `HowManyBins` sections holding `NumberOfParticles` (cm⁻³),
edges, composition and the derived radius and density.

### 3.4 Sorting, diagnostics, first regrid

- `SortAerosolAtGridPointForCoagulation` (`SortAerosol.h`) checks the list is ordered by population id, then by bin
  edges, and relinks it only when it is not (smallest to largest).
- `SpillBeans` prints a full-state dump to standard output (diagnostic only).
- The driver counts the bins (`NumBins`), gets `TEMP = GetTemp()`, opens `CoagTestOut.txt` and writes its header
  (`Time`, then a `Number` and `Eff. Radius` column pair per bin).
- `RegridAerosol` (section 3.7), `SpillBeans`, then `SizeDistOutput(0, "SizeDist0min.txt")`.
- `MixTimestep = 60` and `CoagTimeStep = 60` are **set in the driver**; the deck's coagulation step (`CoagStep`) is not
  used in this mode.
- Initial totals (number, volume, mass, potassium mass) are summed over sections. Potassium uses `FindChem("K+",1)` and
  `FindChem("K2SO4",1)`, so the test assumes those species exist. The totals are diagnostics only (their `WRITE`s are
  commented out).

### 3.5 Main time loop (600 iterations)

`DO r = 1, INT(10*60*60/MixTimestep)`, 600 passes of 60 s. Each pass:

1. **Time series.** `IF (MOD(t,60.0) <= 0.05)` (every pass with a 60 s step): writes `t/60`, calls
   `SortAerosolAtGridPointForCoagulation`, then writes `NumberOfParticles` and `EffectiveRadius` (µm) of every bin to
   `CoagTestOut.txt`.
2. **Totals.** Number, volume, K₂SO₄ and K masses over occupied sections (diagnostic only).
3. **`StepSectionalCoagulationJacobson(CoagTimeStep)`** (section 3.6).
4. `t = t + MixTimestep`.
5. **Hourly snapshot.** `IF (MOD(t,3600.0) <= 0.05)`: `SizeDistOutput(t, "SizeDist<minutes>min.txt")`; with t = 0,
   eleven files in total.

After the loop the driver recomputes the final totals, deallocates its arrays and leaves the `ELSE IF` chain.

### 3.6 `StepSectionalCoagulationJacobson` (`Coagulation/Coagulation.f90`)

Semi-implicit scheme of Jacobson (*Fundamentals of Atmospheric Modeling*, 2nd ed., sect. 15.2) on the sorted section
list. N = number of sections (`HowManyLinks`), S = number of aqueous species (`HowManyAqChems + HowManyAqCations +
HowManyAqAnions`), O and A the hydrophobic and hydrophilic organic counts.

| # | Stage | What happens | Cost |
|---|---|---|---|
| 0 | Guards and order | Returns if `DoCoagulationFlag` is false or `ThermoBulkMode` is true; sorts the list. If an occupied section's effective radius is smaller than the previous occupied section's, calls `RegridAerosol` and re-sorts (up to 10 times, since a regrid moves particles at most one section). Returns if all sections are empty. | O(N) |
| 1 | Allocate | `CoagArray(N,N)`, `VolumeFractions(N,N,N)`, `AqVolConc(N,S)`, `OrgVolConc(N,O)`, `AqOrgVolConc(N,A)` and length-N vectors, every call. | alloc per step |
| 2 | Monodisperse shortcut | One scalar kernel, explicit number update, concentration rescale, `RecalculateRadius`, return. | O(1) |
| 3 | **Kernel array** | Single-particle quantities (diffusion coefficient, mass, mean-free-path term, Reynolds number, terminal velocity) once per section (`CoagPrecomputeParticle`, environment read once). For every pair (I, J): zero if either is empty; `CoagKernel` if `ConstantKernel`; `CoagKernPairSectional` for two sectional particles; otherwise `CoagKernBrownianDiff + CoagKernGravCollection`. | N² pair evaluations, N precomputations |
| 4 | Volume bookkeeping | Average solution density over occupied sections. Per section: single-particle volume, number, volume concentration, and each species in cm³ per cm³ of air (water and ions with the average solution density, solid salts with `SolidSaltDensity`, organics with their densities). | O(N·(S+O+A)) |
| 5 | **Volume fractions** | `VolumeFractions(I,J,K)`: share of a merged I+J particle's volume assigned to section K, split linearly between the two sections that bracket it; one bracket scan over K per (I, J). | O(N³) |
| 6 | Term lists | The species-independent terms of the update (weights `VolumeFractions × CoagArray`, the K ranges, the gain/loss structure) are collected once into flat lists (`LowerK`, `KStart`, `EntI`, `EntJ`, `EntW`, `EntLast`) in the original I, J, K order. | O(N³) once |
| 7 | **Total volume and number** | Loss sum over J (`LowerSum`) and gain from smaller sections (`UpperSum`), implicit update `V_K = (V_K + dt·Upper)/(1 + dt·Lower)` in increasing K (gains use the already-updated smaller sections); number = V / particle volume. | O(N³) |
| 8 | **Per-species update** | The same update for every aqueous, hydrophobic and hydrophilic organic species, accumulated over the term lists. | O((S+O+A)·terms) |
| 9 | Conservation checks | Total volume and summed species volume must agree before and after to 1e-4, otherwise `ERROR` ("Total Volume concentration is lost during step"). Negative concentrations also raise `ERROR`. | O(N) |
| 10 | Write back | Sections with number > 0: `NumberOfParticles`, species back to mol/particle, `RecalculateRadius`. Empty sections are zeroed. | O(N) × `RecalculateRadius` |
| 11 | **`RegridAerosol(RadiiCurrent = .TRUE.)`** | Section 3.7; the radii are already current, so its `RecalculateRadius` pass is skipped for occupied sections. | O(N) (+ equilibrium on moves) |
| 12 | Free | Deallocates the work arrays. | |

**Kernels** (`Coagulation/CoagulationKernels.h`). `CoagKernBrownianDiff`: Brownian diffusion kernel with slip
correction, mean-free-path and transition-regime corrections and a Reynolds-number (ventilation) correction.
`CoagKernGravCollection`: gravitational collection, using `CollisionEfficiency` and the `TerminalVelocity` difference.
`CoagKernPairSectional` combines the precomputed quantities of two sections with the same expressions.

**`RecalculateRadius`** (`LagrangianParticles/ParticleAttributes.h`) returns early for empty sections, otherwise calls
`ParticleDensity` (volume and density from composition) and `SurfaceTension` (which calls `CalculateIonicStrength`,
`KusikMeissner` and `MeanActivity`).

### 3.7 `RegridAerosol` (`LagrangianParticles/SortAerosol.h`)

1. Re-sorts the list.
2. Calls `RecalculateRadius` on every section (skipped for occupied sections when `RadiiCurrent` is set).
3. For each occupied section: if its radius now exceeds its upper edge (and number > 2e-4), it moves essentially all of
   its particles and mass to the next larger section (leaving 1e-4 behind), recomputes `ParticleDensity` for both and
   runs `FindElectrolyteEquilibrium` on both. If its radius fell below its lower edge, the same happens towards the next
   smaller section (with diagnostic `WRITE`s).
4. Resets each section's temperature from the grid (`GetTempfromGrid`).
5. Checks the first particle for "number without mass".

### 3.8 Output routines

- `SizeDistOutput(t, name)` (`OutputRoutines/OutputRoutines.f90`): a header row and one row per bin with time, gas
  species, edges, effective radius, number and totals of the main aerosol components.
- `SpillBeans` writes RH, density and a full state dump to standard output.
- `CoagTestOut.txt` is written directly by the driver (section 3.5).

## 4. Inputs and outputs

**Inputs**

| File / source | Used for |
|---|---|
| `ASPInputDeck.in` | run-mode flag (`CoagTest`), T, P, RH, `DoCoagulationFlag`, `ConstantKernel` and `CoagKernel`, `Monodisperse`; `CoagStep` is read but not used in this mode |
| `AerosolModes.in` | bin definition (size cutoff, number of bins, lowest bin), lognormal modes and composition |
| Species and equilibrium decks (`AqPhaseChems.in`, `AqPhaseIons.in`, `OrgPhaseChems.in`, `HydrophilicOrgChems.in`, `AqEquilibriumReactions.in`, `AerosolThermoKusikMeissner.in`, gas-phase files) | `SetChemistryParams`: species lists and thermodynamic data |
| Set in the driver | 60 s step, 600 steps (10 h), hourly snapshots, potassium species names |

**Outputs** (13 files)

| File | Content |
|---|---|
| `CoagTestOut.txt` | header row; one row per minute: time (min), then number and effective radius (µm) of each bin |
| `SizeDist0min.txt`, `SizeDist60min.txt`, ..., `SizeDist600min.txt` | eleven size-distribution snapshots |
| `Transcript.TXT`, standard output | banner, `SpillBeans` dumps, `RegridAerosol` diagnostics ("Shrink", "Below 0!") when triggered |

## 5. Where the time goes

With the kernel and cubic loops reduced (patches 11–15 in `CHANGELOG.md`), most of the remaining cost is
`RecalculateRadius` in the write-back (about 60% of a RemConMode1 run when measured in October 2026), almost all of it
`SurfaceTension` → `KusikMeissner` → `CalcLogGammaPure`, plus the coagulation routine's own loops. A 10-hour run takes
well under a second on the test decks (Williams 0.63 s, Timbavati 0.83 s, RemConMode1 0.62 s, one run each). Not
pursued: kernel symmetry (the kernel is evaluated for all ordered pairs), per-step allocation, and the output.
