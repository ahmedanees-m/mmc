# Model call ledger

Running record of model API calls and cost against the v4 budget. Committed as it accrues, per
`prereg/PREREG_v4.md` §12.

**Cap:** 1500 calls across all steps. **Reserve:** USD 20 held back for re-runs after review.

Calls routed to the NGC free tier carry no marginal cost and are tracked separately so the
Anthropic spend and the panel volume do not have to be disentangled later.

## Anthropic

| Date | Step | Model | Calls | Cached input | Cost (USD) | Note |
|---|---|---|---|---|---|---|
| 2026-08-17 | 1 (S1) | claude-opus-4-8 | ~4 | no | ~0.20 | Proposal loop on Cytokine_production, 3 iterations, 50 edges |
| 2026-08-17 | 1 (S1) | claude-opus-4-8 | ~3 | no | ~0.15 | Proposal loop on Th2_GATA3, 2 iterations, 14 edges |

**Running total:** approximately 7 calls, USD 0.35.

Counts are the loop's propose and repair calls inferred from the recorded iteration
history rather than read from a response header, so they are stated as approximate.
`claude-opus-4-8` is used for S1 rather than the models locked in PREREG_v4 section 7,
because S1 has to represent the proposer the existing corpus was built with; the
section 7 identifiers apply to the Step 6 panel. TCR_signalosome and CD4_lineage_TFs
reuse structures frozen earlier and cost nothing here.

## NGC hosted

| Date | Step | Model | Calls | Retry rate | Note |
|---|---|---|---|---|---|
| | | | | | |

**Running total:** 0 calls.

## Budget by step

Planned allocation from the v4 plan. Actuals fill in above as each step runs.

| Step | Planned calls | Planned USD |
|---|---|---|
| 1 structure-source comparator | 100 | 5 to 8 |
| 3 power curve | 250 | 12 to 20 |
| 4 Norman standard split | 150 | 8 to 12 |
| 5 anonymisation ablation | 360 | 18 to 25 |
| 6 proposer panel (Anthropic share) | 300 | 15 to 25 |
| 7 scale-up | 300 | 15 to 25 |
| 10 model-class robustness | 100 | 5 to 8 |
