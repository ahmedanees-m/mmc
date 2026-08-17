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

### Outcomes

*(pending)*
