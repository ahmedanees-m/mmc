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

Ordering note: the cytokine module's null had already been reached under the single-band code when this was written, so that module carries a single band at the oracle's edge count and the swept null was run for it separately afterwards.

**Exception recorded 2026-08-18: CD4_lineage_TFs carries the single-band null, not the swept one.** Its swept null was attempted twice and killed both times by host reboots, which occurred at 11:42 and 15:29 on 2026-08-18, roughly four hours apart. The job needs longer than that window on a 33-gene module. Rather than reduce the draws per band and quietly report a thinner null as though it were the pre-registered one, the second pass for that module was scoped to adding the proposal arm and reuses the 1000-draw single-band null from its first pass. The measured result on the cytokine module is that the null mean is essentially flat in edge count, 0.1076 to 0.1152 across 5 to 50 edges, so the single band is a close approximation; the widening upper tail is the part that is lost, and CD4's percentile column is reported with that stated.

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

**A9, 2026-08-17. Step 10's model classes are respecified, because as written they would not test what they were meant to test.**

Fixed before Step 10 runs. Section 11 named three alternative classes, a signed-linear structural model with no gates, a Boolean model, and a per-node MLP, with the MLP as the unbounded-expressivity control that would separate an interpretability cost from a data limit.

The Step 1 decomposition shows that specification does not do that. A structural prediction is the difference between the clamped and unclamped fixed points, so for any perturbation whose gene has no outgoing edges the prediction is identically zero for every gene, and for a gene that has edges only its descendants move. Seventy-eight percent of what the linear baseline achieves on the cytokine module comes from predicting the mean training response, and no structure can emit a response for a perturbation it does not contain. **All three named classes share the same fixed-point solve, the same `do(x = 0)` clamp and the same in-degree cap, so all three inherit that forfeit.** Running them as written would return the same ceiling for all four classes, which section 11 anticipated as the cleanest possible result, but the reason would be architectural rather than anything to do with interpretability or expressivity. The conclusion would look like "expressivity does not help" while actually demonstrating "none of these classes can represent a shared response", which is a different and much weaker statement.

Two classes are therefore added, both chosen to break the forfeit rather than to vary the logic:

1. **Structural plus offset.** The same structure and the same solve, with a fitted per-gene offset added to every perturbation's predicted delta. This hands the model exactly the component it currently cannot express, the stereotyped bulk response, and nothing else. It is the direct test of the section 2.2 mechanism: if structural-plus-offset closes the gap to linear, the ceiling is the forfeit; if it does not, the shortfall is elsewhere and the mechanism is wrong.
2. **Unbounded in-degree signed-linear.** A signed linear structural model with the three-regulator cap lifted, so a target may take any number of regulators. This tests whether the cap alone is binding, separately from the offset.

The three original classes are retained, so the comparison is over five classes plus the linear and mean baselines. The pre-registered cap on the step stands: these two modules only, and no further classes.

**The prediction this makes, recorded now so it can be wrong.** Structural-plus-offset is expected to land close to the mean baseline at 0.3712 and short of linear at 0.4451, because the offset supplies the shared response while the ridge map additionally exploits per-perturbation structure in the low-rank subspace. If instead it reaches linear, the ceiling is entirely the forfeit and the interpretable grammar is not itself the limitation. If it stays near the current oracle at 0.2920, the section 2.2 mechanism is wrong and the explanation has to be rebuilt.

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

**Step 1, all four modules, 2026-08-18. The comparison's answer is set by the module's regime, and the ceiling relative to chance is not.**

With CD4_lineage_TFs complete, the four modules can be read together. Held-out DE-overlap, with each module's own random-structure null.

| Module | Folds | Effective rank | Leading-PC | Specific/shared | Null | zero | mean | linear | oracle |
|---|---|---|---|---|---|---|---|---|---|
| TCR_signalosome | 7 | 1.60 | 0.654 | 0.989 | 0.2639 | 0.3150 | 0.4541 | 0.3929 | 0.4641 |
| Cytokine_production | 19 | 3.64 | 0.276 | 1.336 | 0.1082 | 0.1252 | 0.3712 | 0.4451 | 0.2920 |
| CD4_lineage_TFs | 11 | 5.75 | 0.202 | 1.924 | 0.1033 | 0.1220 | 0.0882 | 0.0706 | 0.3078 |
| Th2_GATA3 | 2 | 1.49 | 0.288 | 1.959 | 0.1626 | 0.2294 | 0.2500 | 0.2500 | 0.7500 |

Th2 carries two scoreable folds and is excluded from every reading below, per amendment A1.

**The baselines swing with the regime; the ceiling does not.** Expressed as a multiple of each module's own null:

| Module | linear / null | oracle / null |
|---|---|---|
| TCR_signalosome | 1.49 | 1.76 |
| Cytokine_production | 4.11 | 2.70 |
| CD4_lineage_TFs | 0.68 | 2.98 |

The linear baseline ranges over a factor of six, from four times chance on the cytokine module to below chance on CD4. The oracle ceiling ranges over a factor of 1.7, from 1.8 to 3.0 times chance. **The achievable structural ceiling is roughly invariant relative to chance while the linear baseline's performance is almost entirely a property of the module.**

The ordering that governs it is the perturbation-specific ratio, and its mirror the leading-PC fraction. Where responses are shared (TCR at 0.99 specific/shared, leading-PC 0.65; cytokine at 1.34 and 0.28) the mean and linear baselines are strong, because most of what there is to predict is the common response. Where responses are perturbation-specific (CD4 at 1.92, leading-PC 0.20) both collapse. This is the same mechanism the Step 1 decomposition and the Step 3 calibration identified, now visible as a gradient across modules rather than as a single module's property.

**Consequence for the branch, stated plainly.** On CD4 the oracle's paired advantage over linear is +0.2372 [+0.1114, +0.4154], which clears the section 1.2 rule. Read naively that is Branch B on this module. It should not be read that way: **the linear baseline on CD4 scores 0.0706 against a random-structure null of 0.1033 and a zero baseline of 0.1220, so it is performing below chance and below predicting nothing at all.** An advantage over a comparator that is itself below chance is not evidence that a structure predicts well. The oracle's 0.3078 does clear CD4's null at the 100th percentile, so it is finding real signal, but at 3.0 times chance it is in the same band as the cytokine module's 2.7 times, where the same comparison returned Branch A.

The branch was pre-registered to be read on the cytokine module alone, and it stands as Branch A. CD4 does not overturn it; it shows that the comparator, not the structure, is what varies.

**Verification of the CD4 linear anomaly, 2026-08-18. It is a property of the module, and the baseline is sound.**

The linear arm's low CD4 score was checked directly rather than assumed, since an advantage over a broken comparator would not be a result. Two candidate explanations were separated: the baseline is degenerate on this module, or the module is perturbation-specific enough that the typical response the ridge map learns is uninformative about any particular held-out perturbation.

The baseline is not degenerate. Across CD4's eleven scoreable folds the perturbed gene's response signature has median norm 1.03 and is non-zero in every fold, and the resulting prediction has median norm 1.68 and is non-zero in every fold. The map is well posed and produces substantial predictions; it simply ranks the wrong genes.

Counting hits against what a random ranking of the same size would give settles the reading:

| Module | Observed hits | Expected under random ranking | Ratio |
|---|---|---|---|
| Cytokine_production | 63 | 29.25 | 2.15 |
| CD4_lineage_TFs | 13 | 19.94 | 0.65 |

On the cytokine module the signature ridge is strongly informative, recovering more than twice the DE genes a random ranking would. On CD4 it recovers fewer than chance would. The deficit on CD4 is roughly 1.5 standard deviations on a Poisson scale and is **not** significantly below chance, so the correct statement is the weaker and safer one: **on CD4 the linear baseline carries no usable signal, rather than being actively anti-correlated.** The earlier wording of "below chance" overstates what 13 against 19.94 supports and is withdrawn.

This confirms the reading in the table above. CD4's oracle advantage over linear reflects a comparator with nothing to contribute on that module, not a structure that predicts well, and the module's conclusions continue to be drawn against its own null. It also strengthens the regime account rather than complicating it: the same baseline is highly informative where responses are shared and uninformative where they are perturbation-specific, which is the gradient the four-module table shows.

**A second finding from the same check, which changes amendment A9.** A one-hot ridge over which gene was perturbed was run on the identical folds as a differently-constructed linear comparator. Under leave-one-perturbation-out it is degenerate by construction: the held-out perturbation's indicator column never appears in the training rows, so its fitted coefficient is zero and the prediction is exactly zero for every gene. Its score duplicates the zero baseline exactly, 0.1220 on CD4 and 0.1252 on the cytokine module.

The `dense_linear` class added in amendment A9 for Step 10 is built this way and is therefore **not a usable model class under this evaluation protocol**. It is withdrawn from Step 10. The `structural + offset` class, which is the one that actually tests the section 2.2 mechanism, is unaffected and stands. If an unbounded in-degree comparison is still wanted it has to be built as a structural model with the cap lifted rather than as a one-hot regression, and that is left for the write-up to note as untested rather than added late.

**First co-response module, 2026-08-18. The two sources sit in different regimes on every diagnostic, and the co-response side is the better powered one.**

`coresponse_ACTR2` is the first module of the second source to complete, so the two can be compared directly for the first time. One module is not a basis for a claim about its source and none is made here; what follows is recorded because it bears on how the regression must be read, not as a result.

| | folds with DE | effective rank | leading PC | pert-specific ratio | random null |
|---|---|---|---|---|---|
| coresponse_ACTR2 | 19 | 1.69 | 0.562 | 3.00 | 0.1630 |
| the 10 regulon modules | 8 to 12 | 3.13 to 9.00 | 0.053 to 0.407 | 2.33 to 35.81 | 0.0461 to 0.1090 |

Two things follow that were anticipated and are now visible. The random null is markedly higher for the smaller module, 0.1630 against a regulon range topping out at 0.1090, which is what a rank agreement over 20 genes rather than 40 does to the chance level. Ceilings and nulls are therefore not on a common scale across the two sources and no cross-source comparison of either will be made. On effective rank and leading-PC fraction this module falls outside the entire regulon range, which is the non-overlap the confound entry above predicted; on the perturbation-specific ratio it does not, since regulon_YY1 at 2.33 sits below it.

The unexpected part is the fold count. This module carries 19 scoreable folds from 20 perturbations, where the 40-perturbation regulon modules carry 8 to 12. Nearly every perturbation in a co-response module produces detectable DE, which follows from how those modules were built, by clustering genes on their co-response across perturbations. The smaller source is thus the better determined one per module, roughly doubling the folds behind each estimate.

That matters for the prospective test. The five held-out modules are co-response, so the section 8.1 calibration check will be run in the better powered of the two regimes and against 13 co-response training modules of the same shape. The holdout composition, recorded as a limitation when it was drawn and then as an advantage on the confound axis, is on the power axis also the favourable side.

It also sets expectations for the ceiling. `coresponse_ACTR2` reaches 0.3388 against its null of 0.1630, a ratio of 2.1, which is within the range the regulon modules show, so the higher null does not by itself make the co-response modules look better or worse. Its nested honest ceiling of 0.0991 falls below its own null, which is the eleventh module for the pooled comparison and the reason that comparison is being kept.

**Module source is perfectly confounded with module size, recorded 2026-08-18 before the regression is fitted.**

Checking the module definitions behind the scale-up shows the two sources do not merely occupy different parts of the diagnostic range, they have no overlap in shape at all. All 10 regulon modules that passed the screen are 40 genes with 39 or 40 perturbations. All 18 co-response modules are 20 genes with 20 perturbations. There is no 20-gene regulon module and no 40-gene co-response module.

Source therefore cannot be separated from module size, perturbation count, or the fold count that follows from the perturbation count. A coefficient on the source indicator is a coefficient on all four at once. This matters for the metric as well as the model: DE-overlap is a rank agreement over the module's genes, so its chance level depends on how many genes there are, which is why the random nulls will not be on a common scale across the two sources either. Effective rank is bounded above by the smaller matrix dimension, so the regulon modules can reach values the co-response modules cannot attain at all, and only the normalised variant is comparable between them.

Three consequences are fixed now rather than after the numbers are in.

The within-source slopes remain interpretable and are the quantity to read. With 10 regulon modules and 13 co-response modules in the training set, both exceed the four-module minimum the fitting code enforces, so the regression is reported within each source as well as with the pre-specified source covariate. Where the covariate model and the within-source fits agree, the pooled slope can be quoted; where they disagree, the within-source fits are the ones that mean something.

The source coefficient is not interpretable as an effect of source and will not be described as one. It absorbs the size and perturbation-count difference in unknown proportion, and no design in this data set can decompose it.

The section 8.1 prospective test is unaffected on this axis, which is worth stating because the holdout being all co-response looked like a weakness when it was recorded. The five held-out modules are 20-gene co-response modules and 13 further co-response modules sit in the training set, so their predictions are interpolation within one regime rather than extrapolation across the confounded boundary. The test is narrower than intended and it is also cleaner than it would have been with a mixed holdout.

This is a property of how the modules were generated, not a defect introduced during the run: the regulon extractor takes 40-gene transcription-factor neighbourhoods and the co-response clusterer emits 20-gene clusters. It was not noticed when the generators were written and it should have been, because it constrains what the regime map can claim. Reporting it is the remedy available now.

**Step 7 analysis script corrected before any module data reached it, 2026-08-18.**

Reading `scripts/step2_regime_fit.py` against section 8 while the queue ran turned up a mismatch between what the pre-registration specifies and what the code did. Section 8 states that both module sources are carried into the fit and that source is reported as a covariate rather than pooled away. The fit was an ordinary univariate regression of ceiling advantage on one diagnostic at a time, with no source term; source appeared only as separate per-group summary counts, which is not a covariate.

The consequence is not cosmetic. Exercised on synthetic modules built with the structure section 8 describes, two sources sitting in different parts of the diagnostic range with a group offset between them, the univariate form returned a slope of -0.0024 where the generating slope was +0.0100, reversing the sign. Adding the source indicator recovered +0.0115 with a source term of +0.210 against a true offset of +0.200. Had the regression been run as written, the regime map could have reported the relationship backwards. The fit now carries the source indicator, which is what section 8 asked for.

That correction exposes a second problem it cannot solve. A covariate only separates from the diagnostic if the two groups overlap on it, and section 8's own extraction notes that regulon and co-response modules occupy different parts of the diagnostic range. On the synthetic case the perturbation-specific ratio correlated with source at -0.97 with zero range overlap, and the cross-validated fit degraded sharply because the indicator and the diagnostic are effectively one column entered twice. The script now measures this for every diagnostic, reporting the correlation with source and the overlap of the two ranges, and prints a statement that the slope and the source term are not separable whenever the correlation exceeds 0.8 or the overlap falls below 0.2. Where that fires, neither coefficient will be read on its own.

**Decision, 2026-08-18, recorded before the fit is run.** The pre-specified dependent variable is the ceiling's margin over the linear baseline and it stays primary. On modules where the linear arm sits at or below its own random null that margin carries the baseline's noise rather than a property of the module, and this has now been seen on CD4_lineage_TFs and on four of the nine completed Step 7 modules. The ceiling's margin over the module's own null is therefore fitted alongside it and reported, with the count of modules where linear is uninformative printed next to both. The margin over null is a secondary reported quantity and does not replace the pre-specified target. This is a different case from the nested ceiling withdrawn earlier: the null is estimated from 300 structures and is well determined, whereas the nested ceiling rested on 2 to 4 folds.

### Amendment A10, 2026-08-18. The section 8.1 holdout is designated mid-run, before the modules are run.

Section 8.1 requires five modules to be held out of the fit, their ceilings predicted from diagnostics alone, the predictions committed, and only then run. The scale-up queue was built to run all 28 modules in one pass, which would have left nothing to predict prospectively and reduced the primary test of claim C3 to a retrospective holdout. This was found with 9 modules complete and 19 not yet run, so the test is still available and is being set up now rather than abandoned.

The five are drawn uniformly without replacement from the 18 candidates that had not started, with the seed fixed to the date, and the pool, seed, rule and result recorded in `prereg/step7_holdout.json` and committed before the modules are blocked from running. No comparator result exists for any candidate in the pool, so the draw cannot have been informed by an outcome. The drawn modules are `coresponse_ELP2`, `coresponse_NELFCD`, `coresponse_SRSF10`, `coresponse_TADA1` and `coresponse_TUFM`.

One limitation follows from the timing and is recorded rather than worked around. The unrun pool at the moment of the draw was entirely co-response modules, because `regulon_AHR` is the last regulon candidate and it was already running. The prospective test therefore checks calibration only in the co-response regime, and says nothing about whether the map transfers to regulon modules. Had the holdout been designated before the queue started, it could have spanned both. That option was lost by building the queue without it, and reporting the test on the narrower basis is the honest remainder rather than a repair.

The order of operations from here is fixed: the remaining 13 non-holdout modules complete, the regression is fitted on those plus the 9 already done, the five predictions and their 80 percent intervals are committed to the repository, and only then are the five run. Success remains at least 4 of 5 observed ceilings inside the interval, and the section 12 rule stands that a failure makes the regime map descriptive with no further modules added to rescue it.

**Correction to the Step 7 interim entry, 2026-08-18, on two more modules.**

The entry below was written on seven modules and two of its statements do not survive the ninth. It said the honest ceiling is two to five times smaller on every module. On regulon_DNMT1 the ratio is 0.72, and on regulon_YY1 the honest figure of 0.3075 is larger than the leaky figure of 0.1749. Across the nine modules the ratio ranges from 0.21 to 1.76, so there is no uniform drop and the direction is not even fixed.

It also said the attainable ceiling lies between the two estimates. That framing is wrong and is withdrawn. The two are not bounds on a common quantity because they are not scored on the same data: the leaky ceiling is scored on every scoreable fold of the module, while the nested one is scored only on the held-out outer perturbations. A module whose outer perturbations happen to be the better determined ones can score higher held out than in sample, which is what regulon_YY1 shows. The two numbers do not bracket anything.

The reason is visible in the split sizes. Each module reserves 12 outer perturbations, but only 2 to 4 of those carry DE entries, so the honest ceiling rests on 2 to 4 folds. The intervals show the cost directly: regulon_DNMT1 gives 0.2244 [0.0000, 0.4488] on two folds and regulon_ETS1 gives 0.1503 [0.0000, 0.4508] on three, both including zero. Per-module honest ceilings at this precision cannot be ranked against each other, and the module-to-module variation in the ratio is consistent with split noise rather than with a property of the modules.

**Decision, 2026-08-18, superseding the one recorded below.** The Step 7 regression stays on the leaky ceiling with the pre-specified linear form, which is what the pre-registration meant by the ceiling and what the fold count supports. The plan to fit a second regression against the nested ceiling is withdrawn, because a 2 to 4 fold dependent variable would fit split noise. The nested estimator is instead reported pooled across modules, where roughly 28 modules at 2 to 4 folds each give enough folds for an aggregate statement about how much of the leaky ceiling survives held-out selection, with the per-module values shown alongside their intervals and explicitly not ranked. If the pooled comparison shows a systematic shortfall, that is a statement about the estimator applied across the module set, not a per-module measurement.

What stands from the entry below is narrower and still worth having. The leaky ceiling is measured with selection and scoring on the same folds and so cannot be read as evidence that structure generalises across perturbations. The pooled nested comparison is the test of that, and it is not yet done.

**Step 7 interim, 7 of 28 modules, 2026-08-18. The honest ceiling behaves very differently from the leaky one.**

Seven generated modules have completed the scale-up, each 40 genes and 40 perturbations with 8 to 11 scoreable folds. Two quantities are recorded per module: the leaky ceiling, where the oracle selects structure and is scored on the same folds, and the nested honest ceiling, where it selects on an inner set of perturbations and is scored on a held-out outer set.

| Module | Leaky ceiling | Honest ceiling | Null mean | Honest over null | Linear |
|---|---|---|---|---|---|
| regulon_HIF1A | 0.4558 | 0.1154 | 0.0628 | 1.8 | 0.0857 |
| regulon_NFE2L2 | 0.3315 | 0.0681 | 0.0710 | 0.96 | 0.0925 |
| regulon_EGR1 | 0.2948 | 0.1273 | 0.0487 | 2.6 | 0.0411 |
| regulon_STAT3 | 0.2432 | 0.0728 | 0.0461 | 1.6 | 0.1726 |
| regulon_ETS1 | 0.2431 | 0.1503 | 0.0602 | 2.5 | 0.0319 |
| regulon_CREB1 | 0.2147 | 0.1112 | 0.1090 | 1.0 | 0.0981 |
| regulon_USF1 | 0.1524 | 0.0623 | 0.0563 | 1.1 | 0.0239 |

The honest ceiling is two to five times lower than the leaky one on every module. Three of the seven land at or below their own random null once selection and scoring are separated, and the best of them reaches 2.6 times the null where the leaky figure reached 7.3. The ranking is not preserved either: HIF1A has the highest leaky ceiling and a middling honest one, while ETS1 is mid-table on the leaky measure and highest on the honest one.

Both estimators are biased and in opposite directions, so neither is the quantity of interest on its own. The leaky one is biased up because the structure is chosen using the folds it is then scored on. The nested one is biased down because it selects on fewer perturbations than the full set, so some of the drop reflects a smaller selection sample rather than leakage alone. The attainable ceiling lies between them, and the gap is wide enough that a claim resting on either alone would be unsafe.

The consequence for the regime map is that it must be fitted against both. **Decision, 2026-08-18: the Step 7 regression will be reported for the leaky and the nested ceiling separately, using the pre-specified linear form unchanged in both cases.** This adds a second dependent variable rather than altering the model, so no change of form is being made, and the pre-specified form stands as written. If the two regressions disagree in sign or in which diagnostic carries the relation, that disagreement is the result and will be reported as such rather than resolved by preferring one estimator.

Two further points are recorded now so that they are not read back into the data later. First, the leaky ceiling is a maximum over three search seeds and Step 1 used ten, so the two are not directly comparable and no Step 7 ceiling should be set beside a Step 1 ceiling. The spread across seeds is zero on four of the seven modules and large on the other three, reaching 0.099 to 0.456 on HIF1A, so the maximum is a stable statistic on some modules and an optimistic one on others. Second, the linear comparator is at or below the random null on four of these seven modules, which is the same failure already documented for CD4_lineage_TFs. On generated modules the linear arm is frequently uninformative, so the advantage rule against it carries little weight here and the null is the meaningful reference throughout Step 7.

**Defect found and corrected, 2026-08-18.** The scale-up skipped a module whose result file existed, and a result file appears at the first flush rather than at completion. The host reboot recorded below killed regulon_YY1 part-way through and the restart skipped it as though it were finished, leaving a parseable file with an empty comparison table. The guard now tests the content of the file, requiring a populated table and the random null it is judged against, rather than testing that the file exists. An audit of all Step 7 outputs against the corrected test found regulon_YY1 to be the only module affected, and it has been re-queued and is running. No completed module was disturbed.

**Step 1 complete, CD4_lineage_TFs second pass, 2026-08-18. All four modules, both passes.**

33 genes, 32 perturbations, 70 DE entries, 11 scoreable folds. Single-band null of 1000 draws at 13 edges, mean 0.1033, per the exception recorded in amendment A5.

| Source | DE-overlap [95% CI] | Delta vs linear [95% CI] | Clears the rule |
|---|---|---|---|
| mean difference | 0.4339 [0.2547, 0.6270] | +0.3633 [+0.1841, +0.5725] | yes |
| oracle | 0.3078 [0.1935, 0.4683] | +0.2372 [+0.1114, +0.4154] | yes |
| **proposal (claude)** | 0.1249 [0.0500, 0.2120] | +0.0542 [-0.0219, +0.1421] | no |
| zero | 0.1220 [0.0749, 0.1723] | +0.0514 [+0.0038, +0.0985] | yes |
| GRNBoost2 | 0.1056 [0.0545, 0.1611] | +0.0349 [-0.0128, +0.0805] | no |
| mean | 0.0882 [0.0312, 0.1550] | +0.0176 [-0.0160, +0.0533] | no |
| textbook | 0.0783 [0.0359, 0.1232] | +0.0076 [-0.0339, +0.0509] | no |
| linear | 0.0706 [0.0210, 0.1264] | reference | n/a |

**The zero baseline clears the rule against linear on this module.** A prediction of no change at all beats the ridge map with a separated interval. That is the plainest possible statement of the finding already verified for CD4: the linear comparator carries no usable signal here, recovering 13 DE genes where a random ranking of the same size would recover 19.9. Every "clears the rule" entry in this table has to be read against that, and the module's conclusions are drawn against its null of 0.1033 rather than against linear.

Read that way, the proposal arm at 0.1249 sits essentially on the zero baseline at 0.1220 and barely above the null at 0.1033. It finds nothing on this module. The oracle at 0.3078 is 3.0 times the null and is the only structure source with a real margin over chance, which matches its behaviour on the cytokine module at 2.7 times.

**Step 1 is now complete: four modules, seven structure sources, both passes.** Across all four, no proposer-produced structure clears the linear baseline anywhere the comparator itself is sound. The pre-registered branch, read on the cytokine module as specified, is Branch A, and nothing in the remaining three modules disturbs it: TCR is underpowered at seven scoreable folds, Th2 at two, and CD4's comparator is uninformative. The substantive content of the other three modules is not the branch but the regime gradient they establish, recorded above.

**Step 1 second pass, cytokine module, 2026-08-18. The complete table, with the proposal arm and the swept null.**

The second pass re-evaluated the four structures the first pass had already selected, added the S1 proposal arm, and replaced the single-band null with the swept version of amendment A5. Structures were reused rather than re-searched, so the oracle row is identical to the first pass by construction.

| Source | Edges | DE-overlap [95% CI] | Delta vs linear [95% CI] | Advantage | Null band | Percentile |
|---|---|---|---|---|---|---|
| linear (S6) | n/a | 0.4451 [0.3280, 0.5716] | reference | n/a | n/a | n/a |
| mean (S7) | n/a | 0.3712 [0.2626, 0.4899] | -0.0739 [-0.2555, +0.1055] | no | n/a | n/a |
| oracle (S5) | 7 | 0.2920 [0.1646, 0.4464] | -0.1531 [-0.2837, +0.0133] | no | 5 | 100.0 |
| mean difference (S4) | 30 | 0.2468 [0.1530, 0.3529] | -0.1983 [-0.3532, -0.0646] | no | 30 | 100.0 |
| **proposal, claude-opus-4-8 (S1)** | 50 | **0.1939 [0.1255, 0.2659]** | **-0.2513 [-0.3993, -0.1185]** | **no** | 50 | 100.0 |
| zero | n/a | 0.1252 [0.0958, 0.1549] | -0.3199 [-0.4549, -0.1967] | no | n/a | n/a |
| textbook (S2) | 16 | 0.1211 [0.0734, 0.1768] | -0.3241 [-0.4610, -0.2020] | no | 16 | 90.0 |
| GRNBoost2 (S4) | 30 | 0.1113 [0.0767, 0.1449] | -0.3338 [-0.4717, -0.2080] | no | 30 | 58.5 |

**Branch A is confirmed with every source present.** No structure source clears the linear baseline, and the ordering is stable: the leaky oracle above the algorithmic arms, those above the proposal arm, and the proposal arm above textbook and GRNBoost2.

The proposal arm sits at the 100th percentile of its own null band, so it is finding real signal rather than guessing, and it recovers about 43 percent of the oracle's above-chance advantage ((0.1939 - 0.1152) against (0.2920 - 0.1076)). It is nonetheless below the zero baseline's distance from linear and less than half of the mean baseline. The picture from the earlier record survives intact under the corrected metric and against a properly matched null: the proposals are better than chance and better than textbook structure, and nowhere near a linear map.

**The swept null answers the question amendment A5 raised, and the answer is that the concern was small.** Across 200 draws at each of five sizes the null mean is 0.1076, 0.1088, 0.1092, 0.1091 and 0.1152 for 5, 10, 16, 30 and 50 edges. A random structure's expected score is essentially flat in its edge count on this module, so the single-band null used in the first pass was not materially misleading. That is now measured rather than assumed, which is the point of having run it, and the percentile column is correct for each source's own size. The upper tail does widen with size, from a 95th percentile of 0.1174 at 5 edges to 0.1554 at 50, so a large structure has more room to score well by chance even though its average does not improve.

**Step 3 calibration on CD4_lineage_TFs, 2026-08-18. A calibrated generator with causal ground truth exists on this module, which partly reverses the scope-down in amendment A8.**

Amendment A8 scoped the power curve down to a qualitative boundary because no generator carrying a causal ground truth met the section 4.1 criterion on the cytokine module. That conclusion was drawn from one module. Repeating the gate on CD4_lineage_TFs gives a different answer.

Real CD4 module: effective rank 5.754, leading-PC fraction 0.2023.

| Generator | Effective rank | Leading-PC | Verdict |
|---|---|---|---|
| structural, 12 to 99 edges | 11.9 to 14.4 (+107% to +149%) | 0.023 to 0.067 | reject |
| structural + low-rank, weight 1.0 | 11.553 (+101%) | 0.1026 | reject |
| **structural + low-rank, weight 2.0** | **5.410 (+6.0%)** | **0.2049 (err 0.0027)** | **accept** |
| structural + low-rank, weight 4.0 | 2.679 (+53%) | 0.3346 | reject |
| low-rank plus noise, ranks 2 to 6 | 4.9 to 5.3 (+8% to +14%) | 0.205 to 0.227 | accept |

**A structural generator mixed with a low-rank co-response at weight 2.0 meets both criteria simultaneously on CD4**, at 6.0 percent rank error and 0.0027 leading-PC error. It carries a causal ground truth, so it can serve the sweep, which the low-rank surrogates cannot.

The reason the two modules differ is mechanical and worth stating. A sparse structural generator produces an effective rank of roughly 12 to 14 on both. CD4's real rank is 5.75 and the cytokine module's is 3.64, so the structural generator has less distance to travel on CD4, and a mixing weight exists that lands both diagnostics inside tolerance at once. On the cytokine module the weight that brought rank close (2.0) left the leading-PC fraction off by 0.080, and the weight that brought the leading-PC close (4.0) put rank off by 27 percent; no weight satisfied both.

**Consequence.** Amendment A8's scope-down stands for the cytokine module and is narrowed rather than withdrawn: the power curve is reported qualitatively there, and **can proceed quantitatively on CD4_lineage_TFs** using the structural-plus-low-rank generator at weight 2.0. Claim C4 is correspondingly restored in scope, limited to the higher-rank module and stated as such. The finding that a purely sparse causal generator cannot reproduce either module's identifiability signature is unaffected and remains the substantive result.

This also sharpens the regime account. The modules where a causal generator can be calibrated at all are the higher-rank ones, which are the same modules where the shared-response baselines are weak. Whether a simulator can be made to look like the data, and whether a linear map can predict it, are governed by the same property.

**Step 9 extended to every module with a proposal arm, 2026-08-18. The dissociation is sharpest where the annotation is densest.**

Supersedes the two-module version above, which is retained as the record of the first pass. CD4_lineage_TFs is absent because its second pass had not completed when this ran, and it is added when that lands.

Th2_GATA3, 17 annotated edges within the module, full gene coverage:

| Source | Edges | Hits | Precision | Chance | Ratio | p | Held-out DE-overlap |
|---|---|---|---|---|---|---|---|
| textbook | 9 | 8 | 0.8889 | 0.4056 | 2.19 | **0.0015** | 0.1000 |
| **proposal (claude)** | 14 | 11 | 0.7857 | 0.4031 | **1.95** | **0.0015** | **0.0000** |
| mean difference | 5 | 4 | 0.8000 | 0.3996 | 2.00 | 0.0730 | 0.4167 |
| oracle | 5 | 3 | 0.6000 | 0.3996 | 1.50 | 0.3133 | 0.7500 |
| GRNBoost2 | 21 | 9 | 0.4286 | 0.4029 | 1.06 | 0.4878 | 0.0000 |

**The proposal arm on Th2 recovers 11 of its 14 edges from the annotated regulon, at nearly twice chance and p = 0.0015, and scores exactly zero on held-out prediction.** That is grounded-but-non-predictive in its purest available form: a structure that agrees with the recorded literature to a degree that would pass any annotation-based validation, and that predicts nothing at all. It is the single clearest illustration the project has produced of why annotation agreement cannot substitute for held-out predictive advantage.

The cytokine module adds the mirror image. There the proposal arm recovers 1 of 50 edges, ratio 2.40 at p = 0.35, so it is not significantly enriched, while the textbook structure is enriched 32-fold at p = 0.0005 and the oracle recovers none at all. On TCR_signalosome nothing reaches significance for any source, ratios running 0.23 to 1.67, consistent with that module separating nothing.

Read together the three modules say that annotation agreement and predictive advantage vary independently: on Th2 the proposal is annotation-rich and prediction-empty, on the cytokine module the oracle is the best predictor and annotation-empty, and the textbook structure is annotation-rich and prediction-poor on both. Two of the four patterns section 10 anticipated appear in the same dataset, on different modules.

Th2's caveat stands throughout: two scoreable folds, so its held-out numbers carry no weight. That does not weaken the annotation reading, which does not depend on fold count, but the pairing of a significant enrichment with a zero held-out score should be presented as illustrative rather than as a measured effect size.

**Step 6, proposer panel, 2026-08-18. Zero validated across six model families, and schema validity separates them sharply.**

Cytokine module, two seeds per model, one prompt for every model with no per-model tuning, identical schema validation and retry path. Combines the original run with the rerun at the corrected timeout; where a model was rerun, the rerun supersedes.

| Model | Runs served | Schema failures | Validated | Held-out | Fraction of ceiling |
|---|---|---|---|---|---|
| claude-opus-4-8 | 2 | 0 | 0 | 0.1562 | 0.53 |
| claude-sonnet-5 | 2 | 0 | 0 | 0.1744 | 0.60 |
| openai/gpt-oss-120b | 2 | 0 | 0 | 0.1498 | 0.51 |
| meta/llama-3.1-70b-instruct | 0 | 2 | n/a | n/a | n/a |
| nvidia/nemotron-3-ultra-550b-a55b | 1 | 1 | 0 | 0.1946 | 0.67 |
| z-ai/glm-5.2 | 2 | 1 | 0 | 0.2211 | 0.76 |

**Pooled validated rate 0 of 9, Wilson 95% [0, 29.9%].**

Two findings, and they point in different directions.

*The reliability claim generalises across families.* No model in the panel produced a structure that beat the linear baseline, on a module where the leaky oracle does not either. Combined with Step 1, this closes the objection that the negative reflects one model's weakness: five families that served, spanning Anthropic, OpenAI open-weights, NVIDIA and Z-AI, all land between 0.15 and 0.22 against a linear baseline of 0.4451.

*The panel does not tighten the interval, and should not be presented as if it does.* At n = 9 the pooled Wilson interval is [0, 29.9%], which is far wider than the [0, 4.8%] the existing 76-hypothesis corpus already supports. The panel's contribution is breadth across families, not precision. Claiming a tighter bound from it would be wrong, and section 7's expectation that the pooled interval would tighten substantially is not met at this sample size.

*Schema validity separates the families where held-out score does not.* Four proposals across three models failed the ModelSpec schema after three attempts each. Both Anthropic tiers and `openai/gpt-oss-120b` produced a valid structure on every attempt; `meta/llama-3.1-70b-instruct` failed both of its; `nvidia/nemotron-3-ultra-550b-a55b` failed one of two; `z-ai/glm-5.2` failed one of three. The recurring error is a rule term naming a regulator with no corresponding edge, which is the grammar's coherence constraint rather than JSON validity. This is a concrete, reportable difference between families on a task the held-out metric cannot distinguish them on, and it is the proposal-validity column section 7 asked for.

One counter-intuitive detail worth stating rather than smoothing: among models that did serve, the non-Anthropic ones attained a *higher* fraction of the oracle ceiling (`glm-5.2` 0.76, `nemotron` 0.67) than the Anthropic tiers (0.51 to 0.60), while also being the ones that failed schema validation. Reliability of output format and quality of the structure produced are not the same axis, and on this evidence they are not even positively related.

**A methodological note that belongs in the write-up.** Two of these models were initially recorded as having failed to serve a completion. They had not: the client timeout was 180 seconds, set against an availability probe that asked each model to reply with the word "ok", while the real call is a 20,000-token structured generation. Raising the timeout to 900 seconds recovered `nemotron` entirely and revealed that `llama-3.1-70b`'s failure is schema validity rather than availability. An availability probe must exercise the real workload; verifying that an endpoint answers a trivial prompt establishes almost nothing.

**Step 5, anonymisation ablation, 2026-08-18. Complete on both modules. The structure is name-driven, and the data-conditioning step barely moves it.**

Four arms, three seeds, two modules, 24 runs, no failures. Every A2 and A4 prompt passed the redaction audit before its arm ran. Held-out scores cluster across all four arms in the 0.11 to 0.22 band, so the score does not separate them; edge agreement does, and section 6 named it as the statistic that carries claim C5.

Agreement is read against the within-arm seed ceiling, which is how much two runs of the same arm at different seeds agree.

| Comparison | Cytokine (ceiling 0.360) | TCR (ceiling 0.881) |
|---|---|---|
| A1 vs A2, names removed, data intact | 0.023 (ratio 0.06) | 0.112 (ratio 0.13) |
| **A1 vs A3, names intact, data destroyed** | **0.402 (ratio 1.12)** | **0.790 (ratio 0.90)** |
| A1 vs A4, both | 0.015 (ratio 0.04) | 0.064 (ratio 0.07) |

**Permuting the perturbation labels leaves the proposed structure at or above the level at which the unperturbed arm agrees with itself across seeds.** Removing the gene names collapses agreement to between 6 and 13 percent of that ceiling. Both modules give the same answer, and the cytokine ratio above 1.0 is itself informative: with the real data destroyed the proposals become *more* self-consistent than with it intact, which is what happens when the data contributes variation rather than signal.

Against the four patterns section 6 anticipated, this is the third and strongest: the proposal barely depends on the data.

**The architectural qualification, which must travel with the result.** The loop's initial proposal function receives the module's gene list and a one-line biological context and does not see the response data at all; data enters only through the repair step as a residual summary. With two loop iterations a structure is therefore one name-driven proposal followed by one or two data-informed repairs. The finding is consequently **not** that the model ignores data placed in front of it. It is that **the data-informed repair step barely moves the structure away from its name-driven starting point.**

That is the more useful statement of the two. It localises the failure to a specific stage of the loop rather than attributing it to the model in general, it is consistent with the Step 1 finding that the achievable ceiling is low regardless of proposer, and it is actionable: a loop of this design would have to be strengthened at the repair step, or restructured so the proposal sees data at all.

It also reframes what A2 measures. Removing the names does not merely withhold a hint, it removes the only input the initial proposal has, which is why aliased runs fail to agree even with each other. A2 is therefore a test of whether this loop can work from data alone, and on this evidence it cannot. Claim C5 is reported with that scope.

**Step 7, module extraction, 2026-08-17. Annotation-defined modules are systematically underpowered; data-defined ones are not.**

Extraction and screening only; the comparator runs on the passing set separately. Section 8 requires every extracted candidate be reported including screen failures, and all 59 are in `results/step7_modules.json`.

| Source | Extracted | Passed | Pass rate | Median scoreable folds | Median DE entries |
|---|---|---|---|---|---|
| CollecTRI regulon | 40 | 10 | 25% | 6 | 24 |
| Data-driven co-response | 19 | 18 | 95% | 16 | 52 |

Failure reasons across the 31 rejects: 28 for too few scoreable folds, 26 for too few DE entries, 1 for too few perturbations. The 28 passing modules carry 8 to 19 scoreable folds, median 16, which clears the target of roughly 20 to 25 module-conditions set in section 8.

**The contrast between the two sources is a result, not a nuisance.** A module defined by a transcription factor's annotated regulon usually does not carry enough differential expression to evaluate at all: `regulon_GATA3` assembles 40 genes and 40 perturbations and yields 6 scoreable folds. A module defined by clustering genes on their measured co-response almost always does, at a median of 16.

This is an independent line of evidence for the same conclusion the Step 1 decomposition, the Step 3 calibration and the Step 9 dissociation each reached from a different direction: the atlas's structure is shared response programs, not annotated regulatory neighbourhoods. Here it shows up before any model is fitted, in whether a module built on annotation even has signal to test.

It also has a practical consequence for the regime map. The passing set is dominated by co-response modules, so the fitted relationship between diagnostics and ceiling will be estimated mostly in that regime, and the ten passing regulon modules are the only counterweight. Both sources are carried into the fit and the source is reported as a covariate rather than pooled away.

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
