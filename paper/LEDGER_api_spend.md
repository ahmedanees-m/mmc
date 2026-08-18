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
| 2026-08-17 | 5 anonymisation | claude-opus-4-8 | 48 | 0 | reconstructed: 24 runs over 2 modules, 4 arms, 3 seeds, 2 proposal iterations each |
| 2026-08-18 | 6 proposer panel | claude-opus-4-8 | 4 | not recorded | reconstructed from run records |
| 2026-08-18 | 6 proposer panel | claude-sonnet-5 | 4 | not recorded | reconstructed from run records |
| 2026-08-18 | 6 proposer panel | openai/gpt-oss-120b | 4 | not recorded | via NVIDIA |
| 2026-08-18 | 6 proposer panel | z-ai/glm-5.2 | 6 | not recorded | 2 in the first pass, 4 in the rerun |
| 2026-08-18 | 6 proposer panel | nvidia/nemotron-3-ultra-550b-a55b | 2 | not recorded | recovered after the timeout was raised to 900s |
| 2026-08-18 | 6 proposer panel | meta/llama-3.1-70b-instruct | 0 | all failed | no schema-valid response, recorded as a failure not a spend |
| 2026-08-18 | 1, 3, 4, 7 | none | 0 | n/a | these steps make no model calls |

**Running total: 68 calls**, against a hard cap of 1500.

**This table was reconstructed on 2026-08-19 and is not a contemporaneous record.** Section 12 requires the ledger to be committed as spend accrues, and it was not: the file sat empty while Steps 5 and 6 ran. The counts above are rebuilt from the run records in `results/step5.json`, `results/step6.json` and `results/step6b.json`, which record one entry per proposal iteration. Two limitations follow and are stated rather than smoothed over. Retry counts are not recoverable for Step 6, because the first pass predated the fix that routed Anthropic calls through the counting provider and recorded `calls = 0` for the models the panel was built on. And a call that failed schema validation and was retried inside the harness may be counted once rather than twice, so 68 is a lower bound on calls issued, though an accurate count of proposals obtained.

Going forward the counter in the provider is authoritative and each step appends its own row on completion.
