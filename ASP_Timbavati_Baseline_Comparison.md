# Timbavati: current `asp_v2.2.2` against the code before any Claude changes

**Main finding.** Almost all of the difference in model results between the current branch and the pre-Claude code comes from one bug fix, commit `508ba8d`: gas rate flag 5 multiplied where it should divide, which made OH + HNO₃ about 250 times too fast. The performance work (patches N1–N13, 0–31) changes the Timbavati results by less than 0.1%, except the coupled electrolyte solve (patch 31), which changes particle water and chloride by about 3–4% and gas NH₃ by about 1%. The current branch runs the case about 6 times faster.

Measured 3 November 2026 on the Timbavati baseline deck (`ASP/InorganicInputDecks/Timbavati`), Lagrangian parcel mode, 5 h, gfortran 13.3, `-O3`, one core.

## 1. Baseline used

The unmodified pre-Claude code, `common_asp` `34b7ec7` (M. Alvarado, 25 September 2026), cannot run this case: it stops after about 3 s in the first condensation step (SLATEC `DGAUS8`, "ANS is probably insufficient", from a NaN ionic strength). The baseline is therefore `34b7ec7` plus the two smallest changes needed to build and run it:

- `09afe88`: one-line guard on the re-`ALLOCATE` of `BinEdges` (needed by gfortran);
- `5b4647e` (patch 16): two lines that initialise `IonReset` in the grow branch of `RegridAerosol` (the uninitialised flag causes the abort).

Builds before `7737090` use the makefile of that period (`ASP` `22ef539`) with `-fallow-argument-mismatch` added for gfortran ≥ 10; the optimisation flags are the same as in the current makefile.

## 2. Speed

Timbavati, 160 size bins, alternating runs, CPU time:

| Build | Run 1 | Run 2 | Mean CPU | Mean wall |
|---|---|---|---|---|
| Baseline (`34b7ec7` + `09afe88` + `5b4647e`) | 590.2 s | 580.2 s | **585.2 s** | 597.7 s |
| Current `asp_v2.2.2` (`7ffe10e`) | 95.7 s | 94.1 s | **94.9 s** | 97.7 s |

The current branch uses 84% less CPU time (6.2 times faster).

Along the history (40 bins, one run each, CPU): baseline 160.1 s; after N1–N11 (`a56950d`) 72.7 s; end of the 28–29 October optimisation work (`258ee61`) 46.9 s; patch 28, Newton single-reaction solve (`2c3b2c4`) 35.2 s; patch 31, coupled electrolyte solve (`7ffe10e`) 25.8 s.

## 3. Output differences and the commits that cause them

Totals at 290 min (particles, µg m⁻³) and gas mixing ratios at the end of the run (ppb). The baseline-to-current changes at 40 and 160 bins agree closely (for example NO₂ −79% at both, HNO₃ +202% at both, particle water +149% and +153%); the table uses 40 bins. Each row is the change from the previous point; every point is a commit on `asp_v2.2.2`, built and run with the same deck. "0" means less than 0.05%.

| Step | OA | Water | NH₄ | Cl | NO₃ | SO₄ | O₃ | NO₂ | NH₃ | HNO₃ |
|---|---|---|---|---|---|---|---|---|---|---|
| Baseline value | 11.6 | 0.005 | 1.36 | 9.0e-5 | 2.93 | 3.87 | 75.3 | 5.25 | 1.62 | 1.99 |
| `09afe88`–`258ee61`: optimisation patches N1–N13, 0–25 and early fixes | 0 | +0.1% | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **`508ba8d`: fix gas rate flag 5** | **+1.4%** | **+159%** | **+54%** | **+108%** | **+4.2%** | **+48%** | **+14%** | **−79%** | **−67%** | **+202%** |
| `f3f1107`–`c506969`: 3-body flag 2 fix, new rate flags, `OrgPhaseReactions.in` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `2c3b2c4`: patch 28, Newton single-reaction solve | 0 | 0 | 0 | +0.1% | 0 | 0 | 0 | 0 | 0 | 0 |
| `7ffe10e`: patch 31, coupled electrolyte solve | 0 | −3.6% | +0.2% | −3.3% | +0.4% | 0 | 0 | 0 | −1.1% | −0.1% |
| **Total, baseline → current** | +1.4% | +149% | +54% | +101% | +4.7% | +48% | +14% | −79% | −68% | +202% |

Points built: `34b7ec7`, `a56950d`, `1f14880`, `7c15ef6`, `8766c40`, `258ee61`, `508ba8d`, `f3f1107`, `e4888c4`, `c506969`, `2c3b2c4`, `7ffe10e`.

## 4. The rate-flag bug

Gas-phase rate flag 5 is the JPL/IUPAC form used for OH + HNO₃:

k = k₁ + k₃[M] / (1 + k₃[M]/k₂)

The code before `508ba8d` computed

k = k₁ + k₃[M] · (1 + k₃[M]/k₂)

(`GasPhaseChemicalRates(L) = k1 + k3*Mgas*(1 + k3*Mgas/k2)`), multiplying by the fall-off factor instead of dividing by it. In the Timbavati baseline deck the only flag-5 reaction is

HNO₃ + OH → NO₃ + H₂O; k₁ = 2.4e-14 exp(460/T), k₂ = 2.7e-17 exp(2199/T), k₃ = 6.5e-34 exp(1335/T)

At the case conditions (288 K, 900 mb, [M] = 2.26e19 molecule cm⁻³): k₁ = 1.19e-13, k₂ = 5.59e-14, k₃[M] = 1.52e-12, so

- corrected: k = 1.7e-13 cm³ molecule⁻¹ s⁻¹ (the JPL value);
- old code: k = 4.3e-11 cm³ molecule⁻¹ s⁻¹, **248 times too fast**.

Consequences in Timbavati: with the old code HNO₃ was destroyed by OH and recycled through NO₃ to NO₂ far too quickly. With the fix, HNO₃ at the end of the run is three times higher (+202%) and NO₂ is 79% lower; O₃ is 14% higher. The extra HNO₃, taken up with NH₃ (gas NH₃ −67%), forms more ammonium salts (ammonium +54%, sulfate +48%, nitrate +4%), and the larger hygroscopic mass holds more water (+159%) and chloride (+108%).

The fix applies to every deck that uses gas rate flag 5. Results from the pre-fix code for any such deck overestimate HNO₃ loss by OH.

## 5. Other changes, for completeness

- **Gibbs–Duhem table (N11, `a56950d`):** particle water +0.1%, everything else less than 0.05%.
- **`OLDASP` initialisation (`1f14880`), uninitialised `ParticleDistribution` (`7c15ef6`), patch 17 (`8766c40`):** no effect on this case.
- **3-body flag 2 fix (`f3f1107`) and the rate flags added for the new mechanisms:** not used by the Timbavati baseline deck; no effect.
- **Patch 28, Newton single-reaction solve:** chloride +0.1%, everything else less than 0.05%.
- **Patch 31, coupled electrolyte solve:** particle water −3.6%, chloride −3.3%, NH₃ −1.1% at 40 bins. The worst-first loop it replaces stopped before reaching equilibrium in most Timbavati solves (flow doc `ASP_LagrangianParcel_InformationFlow.md`, section 7.5). In Timbavati, particle chloride is nearly depleted after the first hour and its difference between the two solvers alternates in sign between output times.

## 6. Reproducing

For each commit: `git archive <commit>` into `src/`; if the commit predates `09afe88` or `5b4647e`, apply those commits with `git show <fix> | patch -p1`; build with the `ASP` `build/ASPMakefile.linux.gfortran` of the same period; run the Timbavati deck with line 61 of `AerosolModes.in` set to the number of bins (40 or 160). The comparison sums each column of `OutputFiles/SizeDist290min.txt` over the bins and reads the last line of `GasOut1–3.txt`.
