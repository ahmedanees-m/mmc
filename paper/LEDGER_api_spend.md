# Model call ledger

Running record of model API calls and cost against the v4 budget. Committed as it accrues, per
`prereg/PREREG_v4.md` §12.

**Cap:** 1500 calls across all steps. **Reserve:** USD 20 held back for re-runs after review.

Calls routed to the NGC free tier carry no marginal cost and are tracked separately so the
Anthropic spend and the panel volume do not have to be disentangled later.

## Anthropic

| Date | Step | Model | Calls | Cached input | Cost (USD) | Note |
|---|---|---|---|---|---|---|
| | | | | | | |

**Running total:** 0 calls, USD 0.00.

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
