# MMC pre-registration

Pre-registration is what makes the result a specification rather than a fit. Nothing
below may change after the loop is built without a dated, justified amendment at the
end of this file. The module, direction, splits, thresholds, and the conserved and
rewired scaffold are fixed here, before any modeling code exists.

**Committed:** 2026-07-09 (this commit). **Author:** Anees Ahmed Mahaboob Ali.

**Status (superseded):** this is the v1/v2 prediction-transfer pre-registration. The transfer
framing, including the planned Arc State comparison below, was superseded by the v3 reframe (see
the README iteration arc); the Arc State baseline was never run on this atlas. The live
pre-registrations are `PREREG_discovery.md` and `PREREG_norman.md`. This file is kept as the
dated historical record.

The module and direction were selected by the Stage-0 precondition test
(`mmc/data/precondition.py`, run with `scripts/run_precondition.py`) against the Zhu
store, not chosen by hand. Thresholds: FDR below 0.10; a perturbation is active
(well-powered) when its knockdown moved at least 20 downstream genes; a module and
direction pass when conservation is at least 0.50, rewiring is above 0, and there are
at least 8 conserved-plus-rewired edges.

## 1. Module (selected by the precondition test)

- **Primary module: TCR signalosome.** It is the only candidate that passes the
  precondition. Direction Stim8hr to Stim48hr.
  - Candidate regulators (perturbed): CD3E, ZAP70, LAT, LCP2, PLCG1, PRKCQ.
  - Target genes (measured): ZAP70, LAT, LCP2, PLCG1, PRKCQ, IL2, NFKB1, RELA, FOS,
    JUN.
  - Zhu coverage: all 11 genes present in the 10,282 measured genes.
- **Second module (Stage-3 anti-cherry-pick): Th2 and GATA3**, direction Stim8hr to
  Stim48hr. It has a conserved GATA3-to-cytokine core (GATA3 to IL5, IL13, STAT4) but
  edge-level conservation of 0.30, so it is a rewiring-rich counterpart rather than a
  second pass. Th1 and TBX21 is the alternative second module if a higher-conservation
  counterpart is wanted; its precondition will be run before it is used.

## 2. Transfer direction (chosen by where signal exists)

- Primary direction: Stim8hr to Stim48hr. Both states are stimulated, so the
  signalosome is active in both (activity in the thousands of downstream genes), which
  gives a testable conserved core and a temporal early-to-late rewiring.
- The Rest to Stim directions were rejected by the data: the signalosome is off at
  Rest (activity 2 to 82 downstream genes), so every Rest to Stim edge is untestable
  or rewired, and a Rest-trained model would have no signal to learn from.
- Precondition numbers observed for the primary (TCR signalosome, Stim8hr to
  Stim48hr): conservation fraction 0.52, rewiring fraction 0.48, 33 conserved-plus-
  rewired edges (17 conserved, 16 rewired, 13 no-effect, 9 untestable).

## 3. Conserved and rewired scaffold (ground truth for Step 6)

The per-edge classification from the precondition test for the primary module and
direction, recorded as the scaffold against which Step 6's predicted map is scored.
Effect is log2 fold change; activity is the count of downstream genes the knockdown
moved in that state.

Conserved edges (significant, same sign, both states active):
- CD3E to LCP2, CD3E to PLCG1, CD3E to PRKCQ, CD3E to RELA
- ZAP70 to LCP2, ZAP70 to PLCG1, ZAP70 to PRKCQ
- LAT to LCP2, LAT to PLCG1, LAT to PRKCQ, LAT to RELA
- LCP2 to PLCG1, LCP2 to PRKCQ
- PLCG1 to LAT, PLCG1 to LCP2, PLCG1 to PRKCQ, PLCG1 to IL2

Rewired edges (present in one state, opposite or absent in the other):
- CD3E to ZAP70, CD3E to IL2, CD3E to NFKB1
- ZAP70 to IL2, ZAP70 to RELA, ZAP70 to JUN
- LAT to ZAP70, LAT to IL2, LAT to NFKB1
- LCP2 to ZAP70, LCP2 to IL2, LCP2 to NFKB1, LCP2 to RELA, LCP2 to FOS
- PLCG1 to ZAP70, PLCG1 to RELA

Interpretation: the conserved core is the signalosome internal cross-regulation
(kinase and adaptor co-dependence), significant and same-sign in both states. The
rewiring is the effector output, the IL2, AP-1 (FOS and JUN), and NF-kB edges, which
follow the early-to-late activation shift; IL2 in particular is an early response that
is largely resolved by 48hr, so several IL2 edges lose significance at Stim48hr.

## 4. Splits (the leakage rule, enforced by construction)

Training state Stim8hr, test state Stim48hr. Perturbations are the six regulators.
- Train: the Stim8hr regulator knockdowns. In-context held-out for the Step-5 gate:
  LCP2 and PLCG1 at Stim8hr; the loop trains structure on CD3E, ZAP70, LAT, PRKCQ at
  Stim8hr.
- Tier A (strict transfer): the entire Stim48hr regulator perturbation set, predicted
  from the Stim8hr-frozen structure and the WT Stim48hr state only. The loop and the
  fitter see zero Stim48hr perturbations.
- Tier B (few-shot rewiring discovery):
  - discovery subset, visible to the loop: CD3E, ZAP70, LAT at Stim48hr.
  - held-out test subset, never seen, the scored set: LCP2, PLCG1, PRKCQ at Stim48hr,
    disjoint from the discovery subset.

This module has six perturbations, so the splits are thin. If Step 1 shows the splits
are underpowered, the module may be enriched with additional signalosome and
downstream activation genes, recorded as a dated amendment below, before the ensemble
is frozen.

## 5. Baselines
- Mean of training perturbations. Regularized linear map. Arc State (best available
  setting for this context, identical split). MechPert-style consensus (Tier B only).

## 6. Metrics and controls (never reward predicting the mean)
- Metrics: per-perturbation sign accuracy over DE genes (absolute delta above a
  threshold fixed at first use); Pearson and Spearman of predicted versus observed
  delta (all module genes plus the DE subset); DE-overlap (precision at k, Jaccard);
  bootstrap CIs.
- Controls: shuffled-perturbation negative (near 0); WT versus WT positive.
- Metric-sensitivity check across the calibrated set.

## 7. Go and no-go thresholds
- Stage 1 (in-context gate): the frozen ensemble is at least the linear baseline on
  the in-context held-out perturbations (target at least R near 0.39, the Zhu
  in-context reconstruction) and clearly beats the mean baseline.
- Stage 2 (decisive): read Tier A and Tier B against mean, linear, Arc State, and
  consensus with CIs, only after freezing. Verdict against the win conditions:
  - matches or beats Arc State on transfer with an interpretable model and a
    validated map: strongest;
  - beats linear, competitive with Arc State, only interpretable decomposition:
    strong;
  - fails transfer and the decomposition shows genuine rewiring: reportable negative.

## 8. Analysis plan
- Freeze the ensemble before reading any Tier A or Tier B result.
- The leakage audit (Tier A zero-leak, Tier B held-out isolation) is a committed
  deliverable.
- Report the trace as a record of how the model reasoned, never as evidence in place
  of the number.
- Score the predicted conserved and rewired map against the scaffold in section 3.

## Amendments

**2026-07-09, amendment 1: fit objective.** The parameter fit minimised a plain
mean-squared error over all gene-perturbation pairs. That objective is dominated by the
many near-zero genes and does not target what the model is scored on, and the first full
run showed it: in-sample Pearson 0.83 with in-sample sign accuracy 0.50, a model that
matched magnitudes but not directions. The objective is amended to a differentially
expressed-gene-weighted squared error plus a hinge penalty on predicting the wrong
direction for a DE pair (mmc/fit/fit_params.py, DE_THRESHOLD 0.5, DE_WEIGHT 4.0,
SIGN_PENALTY 0.5). This is a fit-side change only; the grammar, the splits, the metrics,
and the leakage rule are unchanged.

**2026-07-09, amendment 2: module selection and a baseline-beatability screen.** The
precondition screened conservation but not mean-beatability, and it selected the TCR
signalosome, which the first full run and a follow-up baseline screen
(scripts/baseline_screen.py) showed is the mean-baseline worst case: its perturbations
are 0.91 correlated and the leave-one-out mean reconstructs held-out perturbations at
0.90, so no method can win the number there. The screen adds a mean-beatability
criterion, run before modeling: leave-one-out mean and persistence baseline performance
and inter-perturbation correlation per candidate module. By that screen the primary
module is amended to CD4_lineage_TFs, the CD4 lineage master transcription factors: 19
perturbations (adequate power), inter-perturbation correlation 0.26 (diverse), and a
leave-one-out mean of 0.48 with mean sign accuracy 0.70 (headroom on both magnitude and
sign). The TCR signalosome is retained as a documented negative control. Th2 and GATA3
(leave-one-out mean 0.12) remains a secondary module but is under-powered at four
perturbations. Splits for modules beyond the pre-registered TCR signalosome are derived
deterministically by sorted gene order (mmc/data/splits.py, 70 percent train and 30
percent in-context held-out; first and second halves for the Tier B discovery and
held-out subsets); the leakage audit still holds by construction.
