# ASP pending patch notes (not implemented)

Items found while optimizing ASP that need a decision or a larger change. Nothing here is in the code yet. Branch `asp_v_2.2.1` of common_asp; master untouched.

## 1. Coagulation: order sectional bins by actual particle volume

**Symptom.** RemConMode1 in LagrangianParcel mode stops at t = 60 s (first coagulation call) with "Total Volume concentration is lost during step" when the coagulation step is 60 s. The measured volume ratio is 1.000121 (a small gain; tolerance 1e-4). Details and measurements: `ASP_LagrangianParcel_InformationFlow.md`, section 6, item 2.

**Current workaround.** The RemConMode1 deck uses `Coagulation Time Step = 30.0` s. The run then completes (300 min, 40.6 s). The deviation is first order in the step (1.2e-4 at 60 s, 6e-5 at 30 s, 3e-5 at 15 s, 1.2e-5 at 6 s).

**Likely mechanism (not proved).**
- `SortAerosolAtGridPointForCoagulation` (`LagrangianParticles/SortAerosol.h`) orders links by `ParticleDistribution` and `Edges`, that is by the bin's nominal size range.
- `StepSectionalCoagulationJacobson` (`Coagulation/Coagulation.f90`) takes each link's `SingleParticleVolume` from its *actual* `EffectiveRadius`, which condensation changes.
- The partition coefficients `VolumeFractions(I,J,K)` bracket the combined volume `BigV` between `SingleParticleVolume(K)` and `SingleParticleVolume(K+1)`. This only partitions volume correctly if `SingleParticleVolume` increases with K. If a bin is out of order, a combined volume can fall in no bracket (volume lost) or in two (volume gained).
- At t = 60 s in RemConMode1 several bins are out of order. At entry the volumes of bins 24, 25, 26, 27 are 1.7e-13, 3.5e-14, 4.8e-14 and 4.8e-13 cm³.

**How to sort by actual volume.**
1. In `StepSectionalCoagulationJacobson`, after `SingleParticleVolume(1:HowManyLinks)` is filled, test `ANY(SingleParticleVolume(2:HowManyLinks) .LT. SingleParticleVolume(1:HowManyLinks-1))`. If false (the usual case) do nothing, so all runs that pass today stay byte-identical.
2. If true, build a permutation `Order(1:HowManyLinks)` that sorts `SingleParticleVolume` ascending. `indexx.f` is already in the build and gives a rank index; ties must keep their original order (use the link index as a tie-breaker) so the result is deterministic. Treat empty links (`NumberOfParticles = 0`) as before.
3. Permute every per-link array into sorted order before the partition step: `SingleParticleVolume`, `NumConc`, `InitNumConc`, `VolumeConc`, `AqVolConc`, `OrgVolConc`, `AqOrgVolConc`, `SumVolConc`, and the coagulation kernel `CoagArray(I,J)` (rows and columns). Carry `Order` so results can be scattered back.
4. Run the existing partition and update code unchanged on the sorted arrays (`VolumeFractions`, `LowerK`, the `KStart`/`Ent*` term lists, the volume updates).
5. Scatter the updated arrays back to the original link order with `Order`, then continue with the existing code that writes volumes and number back to the particle list.
6. Keep the volume check at 1e-4; it should now pass because brackets exist for every combined volume.

**What changes.**
- Results change only in calls where bins are out of order, since in-order calls take the identity permutation. Williams may have such calls (it currently passes the check), so its output could change slightly; check by comparing before and after.
- A second option is to regrid particles into bins by actual radius before coagulation (as `RegridAerosol` does after condensation). That is a larger change because it moves particles between bins.

**Test plan.** (a) RemConMode1 with a 60 s coagulation step must complete, with the volume check passing. (b) Williams and Timbavati: compare all output files before and after; report which files differ and by how much. (c) CoagTest decks: confirm byte-identical output (their bins stay in order).
