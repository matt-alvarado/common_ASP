# Timbavati: current code against the imported code

This report compares the model as optimised and corrected in October–November 2026 with the code imported into
`common_ASP` on 25 September 2026 (`34b7ec7`, M. Alvarado), on the Timbavati baseline deck. It answers two questions:
how much faster the current code is, and which changes move the results.

## Summary

- **One bug fix accounts for almost all of the change in results**: gas rate flag 5 multiplied where it should divide,
  which made OH + HNO₃ about 250 times too fast (`508ba8d`, section 4). Fixing it raises end-of-run HNO₃ by 202%,
  lowers NO₂ by 79% and raises O₃ by 14%; the extra HNO₃ forms ammonium salts that hold more water.
- **The performance work changes Timbavati results by less than 0.1%**, except the coupled electrolyte solve
  (`7ffe10e`), which changes particle water and chloride by about 3–4% and gas NH₃ by about 1%.
- **The current code runs the case about 6 times faster** (CPU time, 160 bins).

## 1. Method

**Case.** Timbavati baseline deck (`ASP/InorganicInputDecks/Timbavati`), Lagrangian parcel mode, 5 h, gfortran 13.3,
`-O3`, one core, at 40 and 160 size bins (line 61 of `AerosolModes.in`).

**Baseline.** The imported code cannot run this case: it stops after about 3 s in the first condensation step (SLATEC
`DGAUS8`, "ANS is probably insufficient", from a NaN ionic strength). The baseline is therefore `34b7ec7` plus the two
smallest changes needed to build and run it:

- `09afe88`: guard on the second `ALLOCATE` of `BinEdges` (needed by gfortran);
- `5b4647e`: initialise `IonReset` in the grow branch of `RegridAerosol` (the uninitialised flag causes the abort).

**Builds.** For each commit compared: `git archive <commit>` into `src/`; if the commit predates `09afe88` or `5b4647e`,
apply those with `git show <fix> | patch -p1`; build with the `ASP` `build/ASPMakefile.linux.gfortran` of the same
period (`ASP` `22ef539` for commits before `7737090`, with `-fallow-argument-mismatch` added for gfortran ≥ 10; the
optimisation flags are the same as in the current makefile).

**Comparison.** Particle totals: each column of `OutputFiles/SizeDist290min.txt` summed over the bins (µg m⁻³). Gas:
last line of `GasOut1–3.txt` (ppb).

## 2. Speed

Timbavati, 160 bins, CPU time, alternating runs on the same machine:

| Measured | Baseline | Current code | Faster |
|---|---|---|---|
| 3 November 2026, current = `7ffe10e` | 590.2, 580.2 s | 95.7, 94.1 s | 6.2× |
| 4 November 2026, current = `5ae6004` | 519.1, 518.9 s | 84.8, 85.6 s | 6.1× |

The machine was faster on 4 November (the baseline took 585 s on 3 November and 519 s on 4 November), hence the two
rows. On 3 November the mean wall times were 597.7 s and 97.7 s.

At 40 bins (one run each, CPU) the time fell in steps: baseline 160.1 s; 72.7 s after the first set-up and
thermodynamics optimisations (`a56950d`); 46.9 s at the end of the 28–29 October work (`258ee61`); 35.2 s with the
Newton single-reaction solve (`2c3b2c4`); 25.8 s with the coupled electrolyte solve (`7ffe10e`).

The other decks gain less because their electrolyte solves are fewer and simpler: at 160 bins Williams is 4.4 times and
Otavi 4.7 times faster; Alaska is 6.0–6.4 times faster at 40 and 80 bins (the baseline stops Alaska at 160 bins after
4 minutes in the coagulation volume check). Full table: `CHANGELOG.md`, 2026-11-04.

## 3. Output differences and the commits that cause them

Totals at 290 min (particles, µg m⁻³) and gas mixing ratios at the end of the run (ppb), 40 bins. Each row is the
change from the previous row; every row is a commit built and run with the same deck. "0" means less than 0.05%. The
baseline-to-current changes at 40 and 160 bins agree closely (for example NO₂ −79% at both, HNO₃ +202% at both,
particle water +149% and +153%).

| Step | OA | Water | NH₄ | Cl | NO₃ | SO₄ | O₃ | NO₂ | NH₃ | HNO₃ |
|---|---|---|---|---|---|---|---|---|---|---|
| Baseline value | 11.6 | 0.005 | 1.36 | 9.0e-5 | 2.93 | 3.87 | 75.3 | 5.25 | 1.62 | 1.99 |
| `09afe88`–`258ee61`: optimisations of 28–29 October and early fixes | 0 | +0.1% | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **`508ba8d`: gas rate flag 5 fix** | **+1.4%** | **+159%** | **+54%** | **+108%** | **+4.2%** | **+48%** | **+14%** | **−79%** | **−67%** | **+202%** |
| `f3f1107`–`c506969`: 3-body flag 2 fix, new rate flags, `OrgPhaseReactions.in` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `2c3b2c4`: Newton single-reaction solve | 0 | 0 | 0 | +0.1% | 0 | 0 | 0 | 0 | 0 | 0 |
| `7ffe10e`: coupled electrolyte solve | 0 | −3.6% | +0.2% | −3.3% | +0.4% | 0 | 0 | 0 | −1.1% | −0.1% |
| **Total, baseline → `7ffe10e`** | +1.4% | +149% | +54% | +101% | +4.7% | +48% | +14% | −79% | −68% | +202% |

Commits built: `34b7ec7`, `a56950d`, `1f14880`, `7c15ef6`, `8766c40`, `258ee61`, `508ba8d`, `f3f1107`, `e4888c4`,
`c506969`, `2c3b2c4`, `7ffe10e`.

Later changes to the electrolyte solver and coagulation (`87a863e` to `e3faf76`) move the Timbavati 40-bin totals by at
most 1.4e-3, except chloride (within 3.7%) and particle water early in the run (−25% at 20 min, within 2.5% after
60 min), all from `87a863e`; at 160 bins `e3faf76` changes chloride by up to 3.0% and particle water by up to 3.8%
(`CHANGELOG.md`).

## 4. The rate-flag 5 bug

Gas-phase rate flag 5 is the JPL/IUPAC form used for OH + HNO₃:

k = k₁ + k₃[M] / (1 + k₃[M]/k₂)

The imported code computed k = k₁ + k₃[M] · (1 + k₃[M]/k₂)
(`GasPhaseChemicalRates(L) = k1 + k3*Mgas*(1 + k3*Mgas/k2)`), multiplying by the fall-off factor instead of dividing
by it. In the Timbavati baseline deck the only flag-5 reaction is

HNO₃ + OH → NO₃ + H₂O; k₁ = 2.4e-14 exp(460/T), k₂ = 2.7e-17 exp(2199/T), k₃ = 6.5e-34 exp(1335/T)

At the case conditions (288 K, 900 mb, [M] = 2.26e19 molecule cm⁻³): k₁ = 1.19e-13, k₂ = 5.59e-14,
k₃[M] = 1.52e-12, so

- corrected: k = 1.7e-13 cm³ molecule⁻¹ s⁻¹ (the JPL value);
- imported code: k = 4.3e-11 cm³ molecule⁻¹ s⁻¹, **248 times too fast**.

With the imported code HNO₃ was destroyed by OH and recycled through NO₃ to NO₂ far too quickly. With the fix,
end-of-run HNO₃ is three times higher (+202%) and NO₂ 79% lower; O₃ is 14% higher. The extra HNO₃, taken up with NH₃
(gas NH₃ −67%), forms more ammonium salts (ammonium +54%, sulfate +48%, nitrate +4%), and the larger hygroscopic mass
holds more water (+159%) and chloride (+108%).

The fix applies to every deck that uses gas rate flag 5: results from the imported code for any such deck overestimate
HNO₃ loss by OH.

## 5. The other changes

- **Tabulated Gibbs–Duhem integral** (`a56950d`): particle water +0.1%, everything else less than 0.05%.
- **`OLDASP` initialisation (`1f14880`), uninitialised `ParticleDistribution` (`7c15ef6`), NaN guard on the
  Gibbs–Duhem table (`8766c40`):** no effect on this case.
- **3-body flag 2 fix (`f3f1107`) and the rate flags added for other mechanisms:** not used by the Timbavati baseline
  deck; no effect.
- **Newton single-reaction solve** (`2c3b2c4`): chloride +0.1%, everything else less than 0.05%.
- **Coupled electrolyte solve** (`7ffe10e`): particle water −3.6%, chloride −3.3%, NH₃ −1.1% at 40 bins. The
  worst-first loop it replaces stopped before reaching equilibrium in most Timbavati solves (3,762 of 5,262 replayed
  solves at 160 bins violated the equilibrium conditions by more than twice the tolerance; none with the coupled
  solve). Particle chloride is nearly depleted after the first hour, and its difference between the two solvers
  alternates in sign between output times.
