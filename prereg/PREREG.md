# MMC pre-registration

Commit this, filled in, before any modeling code exists. Pre-registration is what
makes the result a specification rather than a fit. Nothing below may change after
the loop is built without a dated, justified amendment in this file.

**Committed:** FILL date. **Commit:** FILL hash. **Author:** Anees Ahmed Mahaboob Ali.

## 1. Module (from the Stage-0 precondition test)
- Primary module: Th2 and GATA3.
  - Targets (program genes): FILL (for example IL4, IL5, IL13).
  - Candidate regulators: GATA3, STAT6, TBX21, FILL.
- Second module (Stage-3 anti-cherry-pick): TCR signalosome (Stim8hr to Stim48hr)
  or Th1 and TBX21. FILL.
- Zhu coverage confirmed: every module gene in the 10,282 measured genes. FILL yes
  or list gaps.

## 2. Transfer direction (chosen by where signal exists)
- Primary direction: FILL (for example Rest to Stim8hr for Th2 and GATA3, since
  GATA3 carries Rest signal).
- Rationale and the precondition numbers observed:
  - conservation fraction: FILL. rewiring fraction: FILL. testable edges: FILL.

## 3. Conserved and rewired scaffold (ground truth for Step 6)
The per-edge classification from the precondition test (conserved, rewired,
untestable), recorded here as the scaffold against which Step 6's predicted map is
scored.
- FILL table: edge, train-condition effect (sign, FDR, r), test-condition effect,
  class.

## 4. Splits (the leakage rule, enforced by construction)
- Train: training-state perturbations; train and in-context-held-out split = FILL.
- Tier A (strict transfer): the entire test-state perturbation set, predicted from
  the WT test state only. The reasoning step and the fitter see zero test-state
  perturbations.
- Tier B (few-shot rewiring discovery):
  - discovery subset (visible to the reasoning step): FILL genes and size.
  - held-out test subset (never seen, the scored set): FILL genes and size,
    disjoint from the discovery subset.

## 5. Baselines
- Mean of training perturbations. Regularized linear map. Arc State (best available
  setting for this context, identical split). MechPert-style consensus (Tier B
  only).

## 6. Metrics and controls (never reward predicting the mean)
- Metrics: per-perturbation sign accuracy over DE genes (absolute delta above
  FILL); Pearson and Spearman of predicted versus observed delta (all module genes
  plus the DE subset); DE-overlap (precision at k, Jaccard); bootstrap CIs.
- Controls: shuffled-perturbation negative (near 0); WT versus WT positive.
- Metric-sensitivity check across the calibrated set: FILL.

## 7. Go and no-go thresholds
- Stage 1 (in-context gate): the frozen ensemble is at least the linear baseline on
  in-context held-out perturbations (target at least R near 0.39, the Zhu
  in-context reconstruction) and clearly beats the mean baseline.
- Stage 2 (decisive): read Tier A and Tier B against mean, linear, Arc State, and
  consensus with CIs, only after freezing. Verdict against the win conditions:
  - matches or beats Arc State on transfer with an interpretable model and a
    validated map: strongest;
  - beats linear, competitive with Arc State, only interpretable decomposition:
    strong;
  - fails transfer and the decomposition shows genuine rewiring: honest negative.

## 8. Analysis plan
- Freeze the ensemble before reading any Tier A or Tier B result.
- The leakage audit (Tier A zero-leak, Tier B held-out isolation) is a committed
  deliverable.
- Report the trace as a record of how the model reasoned, never as evidence in
  place of the number.
