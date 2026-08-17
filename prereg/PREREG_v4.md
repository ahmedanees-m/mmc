# MMC: v4 Pre-registration (identifiability ceiling and structure-source comparison)

*Supersedes the analysis plans in `PREREG.md`, `PREREG_discovery.md` and `PREREG_norman.md`, which stand as the record of the work already run. Nothing below reinterprets a completed run; the outcomes recorded in those files are fixed. This file constrains what is run next.*

**Committed:** 2026-08-17, before any v4 analysis code exists in the tree · **Author:** Anees Ahmed Mahaboob Ali · **Repository state:** cea4e22 · **Primary gate:** `mmc/eval/holdout.py`

---

## 0. What this document is for

The existing result is a characterised negative: across 76 distinct novel hypotheses on the Zhu CD4+ atlas, none appears in a model that beats a regularized linear map held-out, and the same holds for the combinatorial compose test on Norman. That negative is ambiguous in a way the earlier pre-registrations did not anticipate. It cannot currently be separated into "the proposer is weak" and "no structure in this grammar predicts on this data", because only one proposer has been tried and no upper bound on achievable performance has been measured.

Everything here follows from closing that ambiguity. Section 2 is the experiment that decides which paper gets written, and its branch rule is fixed below before the run. The remaining sections are conditioned on it.

## 1. Claim set

Claims the paper is permitted to make, with the rule that would retire each one. Nothing outside this list may be claimed without a dated amendment.

| # | Claim | Primary metric | Comparator | Retired if |
|---|---|---|---|---|
| C1 | No structure source achieves held-out predictive advantage over a linear map on the powered modules | held-out DE-overlap, paired Δ (§1.2) | S6 linear | any source, including S5, has Δ CI lower bound > 0 |
| C2 | The achievable ceiling is bounded well below the linear baseline, and the shortfall is a property of the data rather than of the proposer | Δ(S5, S6) and the ratio Δ(S1)/Δ(S5) | S5 oracle | Δ(S5, S6) lower bound > 0 while Δ(S1, S6) does not clear |
| C3 | Measurable properties of the response matrix, computed before any modelling, predict where the ceiling sits | out-of-sample calibration of the fitted relationship (§8) | intercept-only model | fewer than 4 of 5 prospective modules fall in the 80% prediction interval |
| C4 | The interventional data budget required to clear the baseline exceeds what current atlases provide for a single module | phase-diagram boundary in perturbation count and effect size | real Zhu module located on the same axes | calibration acceptance (§4.1) fails, in which case C4 is reported qualitatively or withdrawn |
| C5 | LLM-proposed structure is conditioned on the perturbation data rather than recalled from the gene names | edge Jaccard between arms against the within-arm seed ceiling (§6) | A1 vs A1' seed replicate | J(A1,A3) is not below J(A1,A1') with a CI excluding zero |
| C6 | The pattern in C1 and C5 is not specific to one model family | per-model Δ and pooled Wilson interval on the validated rate | across the panel (§7) | any panel model produces a validated hypothesis |
| C7 | Structure sources dissociate on annotation agreement and predictive advantage | CollecTRI precision, recall, Jaccard against Δ | across S1, S2, S4, S5 | the two axes rank the sources concordantly |
| C8 | The ceiling is not an artefact of the bounded sum-of-products grammar | Δ under four model classes (§11) | the original grammar | an alternative class clears the baseline where the original does not |

C2 and C5 are the two claims that turn the existing null into something positive. If Step 1 returns Branch B (§2.4) then C1 is withdrawn, C2 inverts, and the paper is a search-quality paper instead.

### 1.1 Primary and secondary metrics

Held-out DE-overlap is primary throughout. ACC_DEG is secondary and reported alongside in every table; it is never used to establish a claim on its own. Where a claim concerns a rate rather than a continuous score, the rate is primary and the interval is Wilson.

### 1.2 The advantage statistic

For a structure source S and a comparator B, evaluated under the same leave-one-perturbation-out folds:

> Δ(S, B) = mean over held-out folds of [ DE-overlap(S, fold) − DE-overlap(B, fold) ]

with a paired bootstrap over folds, 10,000 resamples, 95% percentile interval. **Advantage is established only when the lower bound of that interval exceeds zero.** Comparing two independently-computed marginal CIs and checking whether they overlap is not the test; the pairing is what gives the fold-level design its power, and the earlier runs left that on the table.

Where a result is reported as "separated CIs" for continuity with the earlier record, the paired Δ is reported next to it.

### 1.3 Multiple comparisons

Benjamini-Hochberg at q = 0.05, applied within these families and no others. Families are fixed here so they cannot be redrawn after the fact.

- **Family A.** The module-level test of Δ(S, S6) > 0, over all module-conditions entering the Step 7 scale-up, for a single structure source. One family per source; sources are not pooled.
- **Family B.** The Norman genetic-interaction subtype tests, secondary subtypes only. The primary hypothesis (§5.2) sits outside the family and is tested uncorrected.
- **Family C.** The per-model panel tests in §7.
- **Family D.** The diagnostic-to-ceiling regressions in §3, over the candidate diagnostics.

The two modules in Step 1 are not a family. They are two pre-specified primary modules and both are reported unadjusted, with the adjustment noted as inapplicable rather than omitted.

## 2. Step 1: the structure-source comparator

The point of this step is to put every available way of obtaining an edge set through one identical downstream path, so the only thing that varies is where the edges came from.

### 2.1 Modules

Two, fixed now: the cytokine-production module (54 genes, 456 DE entries, the powered module from `PREREG_discovery.md` §3) and Th2/GATA3. Th2 is underpowered by the earlier precondition and is included deliberately, as the contrast between a powered and an underpowered module is part of what §3 has to explain.

### 2.2 Sources

S1 Claude proposal, existing loop and current prompt, frozen structures reused where they exist. S2 textbook structure, hand-specified from the module's canonical biology and committed before the run. S3 a random ensemble of 1000 structures sampled uniformly over the grammar at edge count matched to S1, giving the null distribution. S4 interventional causal discovery through `scp_infer`: GIES, DCDI-G, DCDI-DSF, Bicycle, AVICI, GRNBoost2 and Mean Difference, with weighted outputs thresholded to binary edge sets and signs taken from the mean-difference direction. S5 the oracle, a greedy forward-backward edge search followed by simulated annealing, scored directly on held-out DE-overlap. S6 the regularized linear map. S7 the mean of training perturbations.

S5 is leaky by construction and is labelled a ceiling estimator in every table, figure and caption. A nested variant is run alongside it: structure selected on an inner split, scored on an outer split never touched during selection. Both numbers are reported wherever the ceiling is quoted, the leaky one as the upper bound and the nested one as the achievable estimate.

### 2.3 Reported per source

Δ against S6 and against S7 with paired CIs; ACC_DEG; percentile position within the S3 null; edge count; pairwise Jaccard against every other source. Seed stability of S5 from 10 seeds, reported as the spread of the ceiling estimate, since a ceiling that moves with the seed is not a ceiling.

### 2.4 Branch rule

Read on the cytokine module, which is the powered one. Th2 is reported but does not decide the branch.

**Branch A** if the lower bound of the Δ(S5, S6) interval is at or below zero. The claim becomes C1 and C2 as written: no structure in this grammar, from any source including one with access to the answer, predicts held-out perturbation responses better than a linear map on this data. The proposer-versus-data ambiguity dissolves.

**Branch B** if Δ(S5, S6) has a lower bound above zero while Δ(S1, S6) does not. The failure localises to proposal and search. C1 is withdrawn. The paper becomes a measurement of how far short LLM mechanistic reasoning falls from an achievable target, and §6 and §7 move from support to centrepiece. This is also the first real-data case in which the gate should grant certification, which the record has never yet shown.

**Branch C** if Δ(S1, S6) has a lower bound above zero. Given 0 of 76 this is unlikely, and it is pre-registered precisely so that it cannot be explained away if it happens. It would mean the earlier module-level results were underpowered rather than decisive, and the paper becomes a positive methods paper with the earlier negatives reinterpreted as a power problem. Honour it if it occurs.

Steps 3 through 11 do not begin until the branch is read and recorded in §14.

### 2.5 Known failure modes

DCDI does not scale past roughly 50 genes and ran to a 48-hour wall in the published benchmark. Modules stay under 40 genes for S4. If a method does not terminate inside 24 hours on the VM it is recorded as non-terminating at that module size, which is a result consistent with the published benchmark and is reported rather than dropped.

## 3. Step 2: identifiability diagnostics

Computed from the perturbation-response matrix alone, on training folds only, before any model is fit: effective rank (participation ratio) of the response matrix; fraction of held-out response variance carried by the leading training PC; the ratio of perturbation-specific to shared variance; structure-equivalence width, meaning the spread of held-out predictions across all structures within epsilon of the best training loss, taken from the S5 search trace; per-edge sign stability across near-optimal fits and across seeds; and the knockdown efficiency and effect-size distribution from the store.

Epsilon for the equivalence width is fixed at 1.05 times the best training loss. It is set now because it is exactly the kind of parameter that is tempting to tune after seeing which value separates the modules.

The regression of ceiling on diagnostics is underpowered at two modules and is not fitted until Step 7 supplies the module count. The diagnostics module and its tests are built now so that the numbers exist for every module from the start rather than being backfilled.

## 4. Step 3: power curve

Semi-synthetic, calibrated to the cytokine module's measured effect sizes, residual noise (bootstrapped, not assumed Gaussian), knockdown efficiency and shared-response component. Either the extended SERGIO shipped with `scp_infer` or the existing structural generator, whichever reproduces the real module's §3 diagnostics more closely.

### 4.1 Calibration acceptance criterion

The simulated module must match the real one on effective rank within 15% relative and on leading-PC variance fraction within 0.05 absolute. Both, not either.

If that fails, C4 drops to a qualitative boundary and the mismatch is itself reported: it would mean real Perturb-seq carries structure that current simulators do not reproduce, which is worth a paragraph and is not a reason to quietly rescope the figure.

### 4.2 Sweep

One axis at a time from the calibrated baseline point, never a full grid. Perturbation count 5, 10, 20, 40, 80, 160. Effect-size multiplier 0.25, 0.5, 1, 2, 4. Measurement noise 0.25, 0.5, 1, 2. Module size 10, 20, 40, 80. Knockdown efficiency 30%, 50%, 70%, 90%. Sampling depth per perturbation swept separately.

At each point: S5, S6, S7, the best-performing S4 method from Step 1, and S1 at subsampled points only, to stay inside the call budget. The deliverable is a phase diagram over perturbation count and effect size with the real Zhu and Norman modules located on it, and a stated answer of the form "clearing the linear baseline on a module of this character needs approximately N perturbations at effect size E, against the N and E the atlas actually provides."

## 5. Step 4: Norman under the field-standard split

The existing Norman null fit on singles only and predicted doubles. That is the strictly harder extrapolation, and the mechanistic argument in `PREREG_norman.md` §9 predicts failure there by construction. The field's standard additive split trains on all singles plus a subset of doubles. A reviewer who knows the benchmark will notice that the harder version was tested and a null reported. Both get run.

### 5.1 Splits

Standard additive: all 105 singles plus a random subset of the 131 doubles, predicting each held-out double, four folds. Compose-only: the existing protocol, reported alongside and labelled as the stricter test. Comparators on both: the structural model with gates, fitted-additive, mean-of-singles, control-mean, GEARS (§9), and the S5 oracle applied to the doubles task so it is known whether a null here is a search failure or a ceiling.

### 5.2 Primary hypothesis

**On the epistasis and suppression subsets, under the standard additive split, the structural model beats fitted-additive on held-out DE-overlap by the §1.2 rule.** These two subtypes are named in advance because they are the cases where one gene masks another, which is a hierarchy question the grammar can represent, rather than a novel-interaction question it cannot. Synergy, redundancy and neomorphism are secondary and fall in Family B.

H0 is that it does not, in which case the null is now reported under the field-standard protocol as well as the stricter one, stratified by the subtypes where mechanism had the best chance, with a ceiling estimate attached. That is a defensible null rather than a protocol artefact, and it is the expected outcome.

## 6. Step 5: anonymisation ablation

Four arms, one harness, same modules, same seeds.

| Arm | Symbols | Data |
|---|---|---|
| A1 | real | real |
| A2 | anonymised (G01, G02, ...) | real, identical |
| A3 | real | perturbation labels permuted |
| A4 | anonymised | permuted |

The alias mapping is randomised per seed so no consistent alias can leak across runs. Gene names are stripped from the residual reader output as well as from the module spec; leaving them in the repair step would defeat the whole arm. A redaction test that asserts no module gene symbol appears anywhere in the serialised prompt for A2 and A4 goes into the test suite and must pass before the run.

Per arm: in-sample fit, held-out DE-overlap, CollecTRI edge overlap, between-arm edge Jaccard, and distance to the Step 1 ceiling.

The comparison that carries C5 is between-arm Jaccard measured against the within-arm ceiling. J(A1, A1'), where A1' is A1 at a different seed, sets how much agreement to expect from an unperturbed replicate. Then J(A1, A2) below that ceiling means the names are doing work, and J(A1, A3) at that ceiling means the data is not. Both differences are bootstrapped over modules with 95% intervals.

Four outcomes are anticipated and all four are reportable. A2 close to A1 means proposals are data-conditioned and the ceiling result is not a recall artefact. A2 well below A1 means the proposals are substantially literature recall. A3 close to A1 is the strongest of the four and the least comfortable: it would mean the proposal barely depends on the data, and that the data-grounded rationale attached to STK11 was narration over a prior. High CollecTRI overlap in A1 with low overlap in A2 sharpens the existing finding, since it would locate the grounding in the literature rather than in the measurements.

## 7. Step 6: proposer panel

Two Anthropic tiers plus four NGC-hosted models from distinct families, one prompt for all of them.

Anthropic identifiers are locked now: `claude-sonnet-5` as the volume model, `claude-opus-5` on a matched subset. The four NGC identifiers are **not** locked in this commit, because the free-tier catalog turns over and locking a model that has been retired is worse than locking nothing. They are fixed by a dated amendment to this file, naming the exact identifiers, before the first panel call. The slots are: one large open frontier model, one Llama or Nemotron instruct model, one reasoning-specialised model, one further distinct family.

Per model: proposal validity rate, retry rate on structured output, novel-edge rate, edge-ablation necessity pass rate, held-out validated rate (primary), fraction of the Step 1 ceiling attained, in-sample fit against the S3 random control, and the A1-versus-A2 anonymisation sensitivity for at least three models. Pooled Wilson interval on the validated rate across the panel, plus per-model intervals.

No per-model prompt tuning. Tuning per model would make the panel a prompt-engineering result rather than a statement about the models, and the temptation to do it will be strongest on whichever model performs worst.

A monotone ordering across the panel without any model clearing the baseline is an informative outcome, not a wash: it would show the proposal step does real work while remaining capped, which supports rather than weakens the identifiability reading of the ceiling.

## 8. Step 7: scale-up and the fitted regime map

Around 20 to 25 module-conditions across Zhu's three states, plus at least five Replogle 2022 modules (K562 or RPE1) for a second cell type, to test whether the map is Zhu-specific.

Module selection is mechanical and fixed now: candidates from pathway annotation and from data-driven co-response clustering, filtered on perturbation count, coverage in Zhu's measured genes, and DE power by the existing precondition. **Every module extracted is reported, including those that fail power screening**, so that the set cannot be curated toward the conclusion. The screen's thresholds are those already in `scripts/baseline_screen.py` and are not adjusted mid-run.

Per module: the §3 diagnostics, S6, S7, the S5 ceiling and its nested variant, S1, and the best S4 structure. Then ceiling-over-linear is regressed on the diagnostics, with cross-validated predictive accuracy reported and the selected diagnostics named.

### 8.1 Prospective test

Five modules are held out of the fit entirely. Their ceilings are predicted from diagnostics alone, the predictions are committed, and only then are they run. **Success: at least 4 of 5 observed ceilings fall inside the 80% prediction interval.** With n = 5 this is a calibration check and is reported as such; it is not powered for a correlation coefficient and none is claimed as primary.

## 9. Step 8: GEARS

Own pinned environment, in its own container, sharing nothing. Run on the Norman standard additive split and on the Zhu module LOO splits where the setup permits, identical metrics and CIs. It appears as one row in the Step 1 table and one row in the Step 4 table, and is not given its own section.

Arc State is scoped out, for GPU budget and because GEARS is the comparator the Norman benchmarks use. Stating that is acceptable; leaving it unaddressed is not. If the GEARS environment cannot be built, the attempt and the failure mode are reported rather than the comparator silently disappearing.

## 10. Step 9: annotation agreement

Every edge set from S1, S2, S4 and S5 scored against CollecTRI by precision, recall and Jaccard, following the procedure Sprengel and Velten used so the numbers are comparable to theirs, and against the Moonen 2026 variant-to-CRE-to-gene links for the immune modules as a second, disease-relevant axis.

CollecTRI coverage of the specific module genes is checked and reported first. If coverage is thin, that is a stated limitation and the scores are reported with the coverage attached, not withheld.

The deliverable is a two-axis figure, annotation agreement against Δ, one point per structure source. A dissociation is expected: LLM edges high on annotation and low on prediction. If S5 scores low on annotation while scoring highest on prediction, that is worth stating plainly, since it would mean what predicts is not what is annotated.

## 11. Step 10: model-class robustness

Three alternatives to the bounded sum-of-products grammar, on the two Step 1 modules, S5 and S1 arms only: a signed-linear structural model with the same fixed-point solve and no gates; a Boolean model in the SCNS and CellNOpt lineage; and a per-node MLP inside the same fixed-point solve as an unbounded-expressivity control.

The MLP arm is deliberately uninterpretable and exists to separate two explanations that otherwise stay tangled. If it also fails to beat linear held-out, the ceiling is not about interpretability at all, which is the strongest form of C1. If it clears the baseline while the interpretable classes do not, that is a different finding, about the price of interpretability, and it gets reported as that rather than folded into the ceiling story.

Three classes, two modules, then stop. This is a robustness check and it is capped here because it is the step most likely to turn into its own project.

## 12. Stopping rules

Module-hunting for a positive is closed and does not reopen on any branch. The five-criterion discovery gate stands, but a hypothesis that clears it is a bonus reported as such, not the spine of the paper.

If Step 1 returns Branch A on the cytokine module, no additional immune modules are run to look for a different answer; the scale-up in §8 proceeds because it measures the relationship, not because it is another attempt.

If the §8.1 prospective test fails, the regime map is reported as descriptive and no further modules are added to rescue the fit.

Preprint goes out after Step 7 regardless of the state of Steps 8 through 10, which are added as a v2.

Hard cap of 1500 model calls across all steps, tracked in `paper/LEDGER_api_spend.md` and committed as it accrues. If a step would exceed its budget, module count or seed count comes down. The pre-registered protocol does not.

## 13. Out of scope

No novel validated regulatory biology. No prediction win in any regime unless §2.4 Branch C or §5.2 forces one under the stated rule, with correction applied. No claim beyond CD4+ T cells, K562 and RPE1. No claim about dynamical model classes on time-course data, which is untested here. No claim that mechanistic modelling is useless; the claim is bounded to this data regime, with the boundary quantified.

## 14. Amendments and outcomes

Amendments are dated, appended below, and state the reason. Silent changes are not acceptable and the value of this file depends entirely on that. The NGC identifiers in §7 are the one item known in advance to need an amendment.

Outcome sections are appended per step as each completes, in the manner of `PREREG_norman.md` §11, including the branch read in §2.4.

### Amendments

**A1, 2026-08-17. Module set and oracle search parameters, fixed before any comparator run.**

Nothing in the Step 1 comparator has been run. What has been done is descriptive: the modules were assembled from the store and their shape and section 3 diagnostics were computed, with no model fitted and no source compared. Three things that section 2.1 left loose or wrong are fixed here rather than after the fact.

*The cytokine module.* Section 2.1 described it as "54 genes, 456 DE entries". That is the candidate set from `PREREG_discovery.md` section 3, not the module the runs actually used. Commit 0aace0b trimmed it to a loop-tractable slice, and every frozen S1 structure and every number in the existing record refers to the trimmed version. Step 1 uses the trimmed module, which is 28 genes and 28 perturbations in Stim8hr, carrying 109 DE entries at FDR 0.10. Two reasons: S1 has to be the structures already proposed, or the comparison is not against the LLM that was actually tested; and section 2.5 caps S4 at 40 genes for DCDI, which the 54-gene set breaks and the 28-gene set does not. The 54-gene figure in section 2.1 should be read as superseded by this paragraph.

*Th2/GATA3 is not a viable second module, and this is stated before the run rather than discovered after it.* `Th2_GATA3` as pre-registered is 7 genes and 7 perturbations, and in Stim8hr it carries 5 DE entries in total, a mean of 0.71 per perturbation with a minimum of zero. Held-out DE-overlap is undefined on a fold with no observed DE gene, so most of its folds drop out and the paired statistic in section 1.2 would rest on three or four folds. It is retained and reported because section 2.1 committed to it and because the contrast between a powered and an unpowered module is part of what section 3 has to explain, but it cannot support a claim and it does not enter any Benjamini-Hochberg family. The branch in section 2.4 was already specified to be read on the cytokine module alone, so no branch condition changes.

*Two further modules are added,* because one informative module cannot say whether a ceiling is a property of the data regime or of one gene set. Both come from the already pre-registered `MODULES` registry and neither was chosen after seeing any comparator result. `CD4_lineage_TFs`, 33 genes and 32 perturbations with 70 DE entries, effective rank 5.75, is the diverse transcription-factor regime the baseline screen selected as having headroom. `TCR_signalosome`, 11 genes and 11 perturbations with 33 DE entries, is the redundant-cascade regime, and its leading-PC fraction of 0.65 is the highest of the five candidates, so it is where the shared-response explanation for a ceiling should bite hardest. Both sit under the 40-gene S4 cap. The primary module for the branch remains the cytokine module; these two are reported alongside it and both enter Family A. The condition is Stim8hr throughout, matching the existing record.

Measured shapes at Stim8hr, recorded here so they cannot be restated later:

| Module | Genes | Perturbations | DE entries | DE per perturbation | Effective rank | Leading-PC fraction |
|---|---|---|---|---|---|---|
| Cytokine_production | 28 | 28 | 109 | 3.89 | 3.64 | 0.276 |
| CD4_lineage_TFs | 33 | 32 | 70 | 2.19 | 5.75 | 0.202 |
| TCR_signalosome | 11 | 11 | 33 | 3.00 | 1.60 | 0.655 |
| Th2_GATA3 | 7 | 7 | 5 | 0.71 | 1.49 | 0.288 |

*Oracle search parameters.* Section 2.2 specified greedy forward-backward plus simulated annealing scored on held-out DE-overlap, without saying at what fold count or optimizer budget the search itself runs. Scoring every candidate under full leave-one-out at the reporting budget is roughly six minutes per candidate on the cytokine module, which puts the pre-registered ten seeds out of reach, so the search and the reporting are separated and both are fixed now.

The search objective is 5-fold cross-validated held-out DE-overlap, folds drawn by the seed, fits at one start and 80 iterations. Candidate batches are screened at 2 folds and only the top 20 go through the 5-fold objective, so screening changes which candidates are considered and never what a candidate is worth. Greedy runs over a 120-pair pool ranked by marginal association, capped at 30 edges; annealing then runs 800 proposals in batches of 24 with temperature falling from 0.02 to 0.001, drawing moves from the full edge space rather than the pool, and the annealing improvement over greedy is reported so that a pool that was silently capping the ceiling would show up. Ten seeds.

**The structure the search selects is then re-scored under the identical full leave-one-perturbation-out protocol at four starts and 250 iterations that every other source is held to, and that re-scored number is the reported ceiling.** The search budget is a nuisance parameter; the reported comparison remains like-for-like.

The nested variant of section 2.2 holds out 30 percent of perturbations as an outer split, runs the identical search on the inner 70 percent, and scores the selected structure on the outer split, which the selection never touched.

**A2, 2026-08-17. Correction to the primary metric: ties in DE-overlap.**

Found while validating the Step 1 harness on TCR_signalosome, before any comparator result was read. It affects the primary metric, so it is recorded here rather than in a footnote.

Held-out DE-overlap ranks genes by predicted absolute effect and takes the top k. The ranking used `argsort`, which breaks ties by array index, and array index is gene order. That is not harmless in this setting. A structural model only moves genes downstream of the perturbed one, so for any perturbation whose gene is not in the structure the model predicts exactly zero for every gene and its entire ranking is a single tie. Index-order tie-breaking then awarded it the first k genes of the module and scored it on whether those happened to be differentially expressed.

The size of the effect: on TCR_signalosome the zero baseline, which predicts no change at all, scored a held-out DE-overlap of 0.32. That is not a chance level, it is an artifact of the module's genes being sorted alphabetically with CD3E, FOS, IL2 and JUN at the front. Every sparse structure inherited the same bonus, and a sparser structure inherited more of it, which biases the comparator in the direction of whichever source proposes fewest edges.

Ties are now broken uniformly at random, with the score averaged over 64 draws to estimate the expected overlap, seeded from the prediction so the metric stays deterministic and reproducible. A prediction with no tie straddling the top-k boundary is arithmetically unchanged, so the linear, mean and other dense comparators keep exactly the numbers they had. The change is confined to predictions that are partly or wholly tied, which in practice means sparse structural models and the zero baseline.

Consequence for the existing record: the structural-model DE-overlap figures reported before this date, including the 0.18 against a linear baseline's 0.45 on the cytokine module, were computed under index-order tie-breaking and are not comparable to figures computed after it. The linear baseline's 0.45 is unaffected. Those runs are recomputed under the corrected metric and the recomputed values are what the paper reports; the earlier values are retained in the repository history rather than edited away. The direction of the earlier finding is not expected to change, since the artifact inflated the structural model rather than the baseline, but the recomputation is what settles that and it is not assumed here.

**A3, 2026-08-17. The S4 arm is reduced on Zhu, and the reason is the data.**

Section 2.2 specified S4 as GIES, DCDI-G, DCDI-DSF, Bicycle, AVICI, GRNBoost2 and Mean Difference, run through `scp_infer`. Five of the seven cannot be run on the Zhu modules, and this was established before any S4 result existed.

`scp_infer` takes an AnnData of cells by genes carrying a per-cell perturbation label, and GIES, DCDI in both variants, Bicycle, AVICI and SDCD all estimate their scores from a population of observations inside each intervention environment. The Zhu store does not hold that. Its manifest records it as derived from `GWCD4i.DE_stats.h5ad`, with log fold change, adjusted p value and z score for each of 33,983 perturbations against 10,282 genes: one summary row per intervention. GIES comes closest to running and still cannot, because its Gaussian score needs a non-singular covariance and one row per environment gives n equal to p. No per-cell Zhu data exists on the compute host or in the project archive, and the raw atlas would not fit the available storage.

Two methods survive the reduction and are run on Zhu: Mean Difference, which is the effect estimate the store already holds, and GRNBoost2, as per-target gradient-boosted regression over the response matrix with perturbations as samples. GRNBoost2 is implemented directly against the arboreto specification that `scp_infer` wraps, because at these matrix sizes the wrapper contributes a large dependency tree and nothing to the result; it is named as such wherever it is reported, and not described as having been run through `scp_infer`.

The other five are not dropped. Norman GSE133344 is present as per-cell counts, so the full suite runs there in Step 4 and is reported in that table. The Step 1 S4 row is explicitly a two-method row on Zhu, labelled with its scope.

This is worth a paragraph in the paper rather than a footnote. A comparator suite released specifically for interventional single-cell causal discovery cannot be applied to a published genome-scale atlas whose distributed form is differential-expression summaries. That is a statement about the practical reach of the method class, and it is exactly the kind of constraint the identifiability argument predicts will bind.

**A4, 2026-08-17. The oracle search maximises leave-one-out, not 5-fold.**

Amendment A1 set the S5 search objective to 5-fold cross-validated held-out DE-overlap, as a compute compromise, while the reported ceiling comes from leave-one-perturbation-out. Validating the harness on TCR_signalosome showed that compromise does not hold. The selected structure scored 0.4710 on the 5-fold objective it was chosen under and 0.2816 under the leave-one-out protocol it is reported at, which put it below the linear baseline and below a zero prediction.

For most sources a gap like that would be a result. For S5 it is a defect, because S5's only job is to be an upper bound: if the search maximises one quantity and the paper reports another, a low reported ceiling can mean the search failed rather than that the data has none to give, and Branch A would be reachable through a weak search. That is precisely the confound section 0.2 was written to remove.

The search objective is therefore leave-one-perturbation-out, the same fold scheme that is reported. The cheaper optimizer budget of one start and 80 iterations is kept, since the fit budget is a nuisance parameter and the selected structure is still re-scored at four starts and 250 iterations for the reported number. This became affordable because a fix to the backend, which was recompiling the loss once per fold, cut a full leave-one-out pass on the cytokine module from 336 to 83 seconds with fitted losses identical to 1e-10. Annealing is reduced from 800 proposals to 480 to hold the ten-seed budget roughly constant.

A1's other parameters stand: the 120-pair greedy pool, screening at 2 folds keeping the top 20, the 30-edge cap, and ten seeds.

**A5, 2026-08-17. The random null is swept over edge counts rather than fixed to one.**

Section 2.2 specified S3 as 1000 structures with the edge count matched to S1. Once the sources were actually built it became clear that one edge count cannot serve them. On the cytokine module the oracle selects 5 to 7 edges, the textbook structure holds 16, the mean-difference and GRNBoost2 arms 30, and the proposal arm 50. A random structure's expected held-out score moves with its size, so a single null placed all of them against a distribution built for one of them, and the percentile column would have compared a 50-edge structure against a 7-edge null.

The null is therefore swept: 200 structures at each of 5, 10, 16, 30 and 50 edges, holding the pre-registered total of 1000, and each source is placed against the band nearest its own edge count. The band nearest the oracle remains the headline null so the figure has one reference distribution. The sweep also answers a question the fixed version could not, which is whether the null itself rises with edge count; if random structures score better simply for being larger, that is worth knowing before any source is credited for a high score.

Ordering note: the cytokine module's null had already been reached under the single-band code when this was written, so that module carries a single band at the oracle's edge count and the swept null is run for it separately afterwards. Every module reports the swept version.

**A6, 2026-08-17. The four NGC panel identifiers are locked.**

Section 7 deferred these to a dated amendment because the free-tier catalogue turns over. Availability was established by calling the completions endpoint for each candidate rather than by reading the catalogue, and the two disagree sharply: the key lists 102 models, and of the 41 plausible panel candidates probed, **18 served a completion and 23 did not**, most returning HTTP 404 from the completions endpoint despite being listed. Notably every DeepSeek model failed, so the "DeepSeek-R1 class" slot in section 7 has to be substituted, which section 7 permits provided the substitution is recorded.

The four slots, all verified to return a completion on 2026-08-17:

| Slot (section 7) | Locked identifier | Family |
|---|---|---|
| Large open frontier model | `openai/gpt-oss-120b` | OpenAI open-weights |
| Llama or Nemotron large instruct | `meta/llama-3.1-70b-instruct` | Meta |
| Reasoning-specialised | `nvidia/nemotron-3-ultra-550b-a55b` | NVIDIA |
| Further distinct family | `z-ai/glm-5.2` | Z-AI |

Four distinct families, as section 7 requires. The reasoning slot substitutes an NVIDIA Nemotron reasoning model for the intended DeepSeek-R1 because no DeepSeek model is served on this key; `nvidia/nemotron-3-nano-omni-30b-a3b-reasoning` and `thinkingmachines/inkling` were the alternatives and are recorded here as the fallbacks if the locked one is withdrawn.

One implementation consequence, which the panel adapter has to handle or it will silently score several models as producing nothing: **the models differ in where they put their answer.** `openai/gpt-oss-120b`, `nvidia/nemotron-3-ultra-550b-a55b` and `thinkingmachines/inkling` return a populated `reasoning_content` field alongside `content`, and an early probe that read only `content` recorded `gpt-oss-120b` as returning an empty string. `meta/llama-3.1-70b-instruct` and `z-ai/glm-5.2` use `content` alone. The adapter reads both fields, and the proposal-validity rate in section 7 is reported only after that is verified per model, since a parsing failure would otherwise be indistinguishable from a model that cannot follow the schema.

**A7, 2026-08-17. Norman genetic-interaction subtypes, and two baseline changes to the standard split.**

Fixed before the standard split runs. Section 5.2 names epistasis and suppression as the primary hypothesis, and those are subtype labels rather than a continuous score, so they have to be obtainable before the hypothesis is testable.

*Subtypes.* Norman's published coefficient table is not shipped with the GEO matrix, which the earlier work already hit and worked around by recomputing non-additivity from the data. The subtype definitions themselves are recorded in `PREREG_norman.md` section 2: synergy is two large coefficients, suppression two small, epistasis an asymmetric pair in which one single accounts for the double while the other contributes almost nothing, and neomorphism a large deviation from any additive fit. They are therefore derived here from the additive coefficients fitted per pair, using quantiles of this dataset's own coefficient distribution rather than absolute thresholds, so that no cut can be tuned toward an outcome:

- **neomorphic**: additive-fit deviation in the top decile
- **epistasis**: coefficient asymmetry min(|c1|,|c2|)/max(|c1|,|c2|) in the bottom tertile, with the coefficient total at or above the lower tertile
- **suppression**: coefficient total in the bottom tertile
- **synergy**: coefficient total in the top tertile
- **additive**: everything else

The labels are operational reconstructions, not Norman's published assignments, and are described that way wherever they appear.

*Two baseline changes, both forced by the Step 1 result.* On the cytokine module 78 percent of everything the linear map achieved came from predicting the mean of the training perturbations. Any protocol that lets a model see some doubles has to be compared against simply averaging them, so **mean-of-training-doubles** joins the comparator set. Without it a model could look strong for reproducing a shared response it never had to reason about, which is precisely the failure the Step 1 decomposition exposed.

Second, the additive baseline is **fitted on the training doubles and applied to the held-out ones**. The compose test fitted the additive coefficients per pair on the pair being predicted, which sees the answer and is an oracle rather than a baseline. That per-pair version is retained under the name `fitted_additive_oracle` so the two cannot be confused, and the honest `additive_trained` is the reference the primary hypothesis is tested against.

**A8, 2026-08-17. Step 3 calibration fails its acceptance criterion, and the power curve is scoped down.**

Section 4.1 made the power curve conditional on the simulated module matching the real one on effective rank within 15 percent relative and leading-PC variance fraction within 0.05 absolute, both required, and stated that a failure is reported as a finding rather than worked around. It failed. This records the result and the resulting scope change before any sweep is run.

Real cytokine module: effective rank 3.639 of 28 perturbations, leading-PC fraction 0.2758.

| Generator | Effective rank | Leading-PC | Verdict |
|---|---|---|---|
| structural, 12 edges | 12.309 (+238%) | 0.0389 (err 0.237) | reject |
| structural, 24 edges | 12.025 (+231%) | 0.0519 (err 0.224) | reject |
| structural, 48 edges | 10.274 (+182%) | 0.0837 (err 0.192) | reject |
| structural, 84 edges (in-degree cap saturated) | 10.968 (+201%) | 0.0779 (err 0.198) | reject |
| structural + low-rank co-response, weight 2 | 4.511 (+24%) | 0.1955 (err 0.080) | reject |
| structural + low-rank co-response, weight 4 | 2.648 (+27%) | 0.3151 (err 0.039) | reject |
| low-rank plus noise, rank 2 | 3.483 (+4.3%) | 0.2779 (err 0.0022) | accept |
| low-rank plus noise, rank 6 | 3.578 (+1.7%) | 0.2749 (err 0.0008) | accept |

**No generator with a causal ground truth comes within tolerance, at any density or any mixing weight.** Density was swept precisely because a sparse graph routes each knockdown into a small, largely disjoint descendant set, which makes responses near-orthogonal and the rank high; saturating the grammar's in-degree cap at 84 edges moved the effective rank only from 12.3 to 11.0, nowhere near 3.6. Mixing in a low-rank co-response taken from the real data can hit either criterion separately but never both: at weight 2 the rank is close and the leading-PC is off by 0.080, at weight 4 the leading-PC is close and the rank is off by 27 percent.

The only generators that pass are low-rank-plus-noise surrogates with no causal structure at all, and they pass almost exactly, at 1.7 to 4.3 percent rank error. **Those cannot serve the power curve.** The sweep asks how much data is needed before structural modelling begins to win, and that question is meaningless against a generator with no structure to recover. Passing the criterion is necessary, not sufficient; the generator also has to contain the thing being searched for.

**Scope change.** Per section 4.1, the power curve is reported as a qualitative boundary rather than a numerical data budget, and claim C4 in section 1 is downgraded accordingly: the paper will not state "clearing the baseline requires approximately N perturbations at effect size E".

**What replaces it is stronger than the qualitative fallback suggests.** The mismatch is itself the finding, and it explains the Step 1 result rather than merely accompanying it. The Zhu response matrix is well described as a few shared response programs plus noise and badly described as arising from a sparse causal graph. That is why a linear map, which exploits low-rank structure directly, beats every sparse structural source including one selected on the held-out answer, and it is why the identifiability diagnostics of section 3 predict the ceiling. The claim moves from "here is the data budget you would need" to "the data does not have the shape this model class assumes, and here is the measurement showing it".

One correction to this script's own reporting is recorded for completeness. An earlier version mixed the co-response in as a constant vector added to every row and produced identical diagnostics at every weight. Both diagnostics are computed on the column-centred response matrix, so a constant offset is removed exactly by the centring and cannot change either. The component is now mixed per perturbation as a scaled low-rank term. This also sharpens the reading of the Step 1 decomposition: the low effective rank is not one offset shared by all perturbations but different perturbations producing scaled versions of a few common programs.

### Outcomes

**Step 1, primary module (Cytokine_production / Stim8hr), 2026-08-17. Branch A.**

Read on the powered module as section 2.4 requires. 28 genes, 28 perturbations, 109 DE entries, 19 folds carrying a DE gene. All figures are held-out DE-overlap under the corrected metric of amendment A2, with the paired statistic of section 1.2.

| Source | Edges | DE-overlap [95% CI] | Δ vs linear [95% CI] | Advantage | Selected from |
|---|---|---|---|---|---|
| linear (S6) | n/a | 0.4451 [0.3280, 0.5716] | reference | n/a | not a structure |
| mean (S7) | n/a | 0.3712 [0.2626, 0.4899] | -0.0739 [-0.2555, +0.1055] | no | not a structure |
| oracle (S5) | 7 | 0.2920 [0.1646, 0.4464] | -0.1531 [-0.2837, +0.0133] | no | held-out score |
| mean difference (S4) | 30 | 0.2468 [0.1530, 0.3529] | -0.1983 [-0.3532, -0.0646] | no | full response matrix |
| zero | n/a | 0.1252 [0.0958, 0.1549] | -0.3199 [-0.4549, -0.1967] | no | not a structure |
| textbook (S2) | 16 | 0.1211 [0.0734, 0.1768] | -0.3241 [-0.4610, -0.2020] | no | prior literature |
| GRNBoost2 (S4) | 30 | 0.1113 [0.0767, 0.1449] | -0.3338 [-0.4717, -0.2080] | no | full response matrix |

Random-structure null (S3), 1000 draws at the oracle's edge count: mean 0.1082.

**The branch is A.** The oracle does not clear the linear baseline; its paired interval against linear runs from -0.2837 to +0.0133 and therefore fails the section 1.2 rule, and its point estimate sits 0.15 below. Since S5 is an upper bound selected on the held-out answer itself, no structure in this grammar clears the bar on this module, and the limit is a property of the data rather than of any proposer. Section 2.4's Branch A language stands as written and C1 is not retired.

Supporting numbers, all pre-registered:

- **Ceiling seed stability (section 2.3).** Ten independent searches gave 0.1690 to 0.2920, sd 0.038. A ceiling that holds this steady is not a search artefact.
- **Nested honest variant (section 2.2).** Structure selected on an inner split, scored on an outer split selection never touched: 0.1333 [0.0797, 0.1828] over 8 folds. The gap to the leaky 0.2920 is the size of the selection advantage, and the honest number sits barely above the random null.
- **The oracle is finding something, and it is not enough.** At 0.2920 it is roughly 2.7 times the null mean of 0.1082, so the search is not merely fitting noise; it is simply far below a linear map.
- **Textbook structure performs at the floor.** 0.1211 against a null of 0.1082. Canonical immunology, compiled into this grammar, predicts held-out perturbation responses no better than a random structure of the same size.
- **Structure-equivalence width (section 3).** 107 structures lie within 5 percent of the best training loss; their held-out scores span 0.188 to 0.304, a width of 0.115. Fitting the training data equally well carries almost no information about predicting held-out responses, which is the mechanism behind the ceiling.
- **Sign stability (section 3).** Across that near-optimal class, 84 percent of the 50 distinct edges have a determined sign. `VARS2 -> IL5` appears activating in one search and repressing in another.
- **What the oracle selected.** `ATP5F1A -> IL22`, `PGS1 -> IL2`, `VARS2 -> IL5`, `INO80B -> IL22`, alongside `CD28 -> IL2`. Mostly metabolic and housekeeping genes rather than regulators, which is what an upper-bound search returns when there is little real signal to find.

One measurement worth recording against amendment A2: the linear baseline recomputed under the corrected metric is 0.4451, against 0.45 reported previously. The correction left the baseline where it was, as predicted, and moved only the sparse structural arms. The zero baseline fell to 0.1252 and now sits essentially on the random null, which is where a prediction of no change belongs.

Remaining Step 1 work before the branch is final for the paper: TCR_signalosome, Th2_GATA3 and CD4_lineage_TFs are running, the S1 proposal arm is evaluated in the second pass, and the swept null of amendment A5 replaces the single-band null above for the percentile column.

**Step 9, annotation agreement, 2026-08-17. The dissociation section 10 anticipated, and it is complete.**

CollecTRI via the OmniPath export: 62,404 interactions, 57,914 with an unambiguous sign. Coverage is reported before any score, as section 10 requires. Each source's precision is placed against a permutation chance level drawn from random edge sets of the same size over the same gene universe, because precision against a sparse annotation is not interpretable on its own.

Cytokine module, 6 annotated edges lying within it, 71 percent of its genes appearing in CollecTRI at all:

| Source | Edges | Hits | Precision | Chance | Ratio | p | Held-out Δ vs linear |
|---|---|---|---|---|---|---|---|
| textbook | 16 | 4 | 0.2500 | 0.0078 | **32.0** | **0.0005** | -0.3241 |
| mean difference | 30 | 4 | 0.1333 | 0.0080 | 16.6 | 0.0005 | -0.1983 |
| oracle | 7 | 0 | 0.0000 | 0.0084 | **0.00** | 1.0000 | -0.1531 |
| oracle nested | 5 | 0 | 0.0000 | 0.0074 | 0.00 | 1.0000 | n/a |
| GRNBoost2 | 30 | 0 | 0.0000 | 0.0080 | 0.00 | 1.0000 | -0.3338 |

**On this module the two axes are not merely uncorrelated, they are inverted.** The oracle carries the largest held-out advantage of any structure source, at 0.2920 against a null of 0.1082, and recovers **zero** annotated edges, at a permutation p of 1.0. The textbook structure predicts at the floor, 0.1211 against that same null, and is **32 times enriched** for annotated edges at p = 0.0005. What predicts is not what is annotated, and what is annotated does not predict.

Section 10 named this as the striking outcome to watch for and it has occurred, with a permutation test behind it rather than an eyeball comparison. It is the cleanest available statement of why validating AI-proposed mechanism on annotation agreement is miscalibrated: on this module, agreeing with the literature is anti-correlated with predicting held-out responses.

It also completes the account built up in the Step 1 and Step 3 outcomes. The data is low rank rather than sparse causal (A8), so the structures that predict best are those that happen to capture the dominant shared programs, and those are not regulatory edges. The oracle's selections bear this out directly: `ATP5F1A -> IL22`, `PGS1 -> IL2`, `VARS2 -> IL5`, metabolic and housekeeping genes rather than regulators.

TCR_signalosome, with 22 annotated edges and full gene coverage, shows no significant enrichment for any source: ratios run 0.23 to 2.59 with p from 0.082 to 0.99. Consistent with its low fold count, that module separates nothing, and it is reported as uninformative rather than as a weak version of the cytokine result.

**Step 4, standard additive split, 2026-08-17. H0, the pre-registered null.**

105 singles and 131 scorable doubles, four folds, DE-overlap at k = 50 on the held-out doubles. Subtypes derived per amendment A7: synergy 36, suppression 33, additive 33, epistasis 15, neomorphic 14. The trained-additive baseline's fitted coefficients were stable across folds (0.60 to 0.65 on each single).

| Method | DE-overlap@50 [95% CI] | Δ vs additive_trained | Advantage |
|---|---|---|---|
| fitted additive oracle | 0.5127 [0.4832, 0.5421] | +0.0153 [+0.0084, +0.0229] | yes, but it sees the answer |
| mean of singles | 0.4979 [0.4687, 0.5269] | +0.0005 [-0.0018, +0.0027] | no |
| sum of singles | 0.4979 [0.4687, 0.5269] | +0.0005 [-0.0018, +0.0027] | no |
| additive trained | 0.4974 [0.4682, 0.5264] | reference | n/a |
| structural | 0.4751 [0.4476, 0.5034] | -0.0223 [-0.0382, -0.0087] | no |
| mean of training doubles | 0.1879 [0.1745, 0.2009] | -0.3095 [-0.3368, -0.2831] | no |
| zero | 0.0047 [0.0031, 0.0067] | -0.4927 | no |

**Primary hypothesis (section 5.2), tested uncorrected as pre-registered:** the structural model does not beat the trained-additive baseline on either named subset. Epistasis, n = 15: 0.5560, Δ +0.0013 [-0.0147, +0.0147]. Suppression, n = 33: 0.4739, Δ -0.0424 [-0.0909, -0.0036], which is a significant deficit rather than an advantage. H0 stands, and the null is now reported under the field-standard protocol as well as the stricter compose-only one, stratified by the subtypes where mechanism had the best chance.

Two observations worth carrying into the write-up.

*The shared-response shortcut does not transfer.* On the Zhu cytokine module the mean of the training perturbations reached 78 percent of the linear baseline's achievable range. Here the mean of the training doubles reaches 0.1879 against an additive baseline's 0.4974, so it is a weak predictor. Norman's doubles are heterogeneous combinations of distinct gene pairs and have no dominant common response, whereas the Zhu module's perturbations do. The contrast supports reading the Zhu ceiling as a property of that data regime rather than of held-out perturbation prediction in general.

*An internal consistency check passed.* `sum_of_singles` and `mean_of_singles` score identically, which they must: DE-overlap at k is rank-based and a positive rescaling cannot change the ranking. That the two agree to four decimals confirms the metric is behaving as specified.
