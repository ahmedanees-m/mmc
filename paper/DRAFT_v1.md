# Structure discovery, not identifiability, limits mechanistic prediction of perturbation responses

Anees Ahmed Mahaboob Ali

*Draft v3, 2026-08-22. Pre-registered as `prereg/PREREG_v4.md`, with twenty-one dated
amendments. Every number below is drawn from the recorded outcomes; sections that no measurement
supports are marked as open rather than filled in.*

*The title changed at v3. Extending the structure search from five seeds to twenty showed that
the ceiling it estimates sits above a regularized linear map on almost every module, which the
five-seed run had understated. The negative therefore is not that structure cannot help. It is
that nothing honest finds the structure that would. Section 3.2 carries the measurement and
section 9 records why the earlier reading was wrong.*

---

## Abstract

Mechanistic models of gene regulation are usually justified on the grounds that a sparse causal
structure, once recovered, should predict the effect of interventions the model was not fitted
on. This study tests that proposition directly on steady-state single-knockdown Perturb-seq
data, by holding the model class, the optimizer budget and the evaluation folds fixed and
varying only where the structure came from. Seven sources are compared, including structures
proposed by a large language model, curated textbook structures, two algorithmic inference
methods, random structures at matched edge count, and an oracle permitted to select its
structure using the held-out responses it is then scored on.

Across 27 gene modules, most comparisons are uninformative: on 14 the linear arm does not beat
its own random-structure null, and a structural advantage over a baseline performing at chance
means nothing. Twelve modules support the comparison.

On those twelve, with the search extended to twenty seeds and verified to have converged, the
oracle ceiling exceeds the linear map on **11 of 12** by point estimate and clears it under a
paired interval rule on 5. The information a structural model would need is therefore present in
the data. No honest procedure recovers it: the nested variant of the same search, scored only on
folds its selection never touched, exceeds the linear map on 5 of 12 and exceeds its own
random-structure null on 7 of 12, a median of 0.21 below the leaky arm. Language-model
proposals, textbook structure and two causal-discovery methods all fall further behind, and none
beats the linear map on any sound module.

One reason the honest search fails is a property of the data rather than of the model class. The
perturbation-by-gene response matrix on the primary module has an effective rank of 3.64 out of
28 perturbations and a leading-principal-component fraction of 0.276. Across seven modules, no purely structural
generator reproduces that signature at any edge density or in any of four topology families
including hub, scale-free and modular, while a structureless low-rank surrogate reproduces it
on every module, to within 1.7 percent on effective rank and 0.0008 on the leading-component
fraction on the primary one. On two of the seven a structural generator does reproduce it, but
only once a shared component four to eight times the structural contribution is added, which is
the component the sparse grammar cannot express.

Three objections are answered by measurement rather than by argument: that the metric is one
sparse models forfeit by construction, that the model class rather than the structure is the
binding constraint, and that a proposer shown only gene names is not being tested on anything
interesting. The first two do not survive; the third is partly correct, and the arm added to
address it shows that the proposer does condition on response data when the data is put in
front of it, without becoming competitive with the linear baseline.

Two results bound what any of this generalises to. The identifiability signature is stable
across all three states of this atlas but is absent in a second one: matched shape for shape,
a combinatorial CRISPRa atlas has two to three times the normalised effective rank. And on that
same second atlas, where the low-rank account does not apply, an untuned graph-prior method
still loses to a fitted-additive baseline by 0.27 in held-out DE-overlap. The difficulty of
beating simple baselines is therefore not a consequence of low rank; the two findings are
separate, and the low-rank account explains this atlas rather than the field.

---

## 1. Introduction

The argument for mechanistic modelling in transcriptomics rests on generalisation. A fitted
regression describes the perturbations it was trained on; a model carrying the right causal
structure should, in principle, predict the response to a perturbation it has never seen,
because the structure encodes how the system actually propagates change. That argument is
usually made rather than measured, and where it is measured, structure discovery and model
fitting are varied together, so a negative result cannot be attributed to either.

This study separates them. Every structure source supplies an edge set and nothing else. The
edge set is compiled into the same structural backend, fitted with the same optimizer budget,
and scored on the same leave-one-perturbation-out folds. The only quantity that varies between
rows of any comparison table is where the edges came from. A source that proposes better
structure should show up as a better held-out score, and nothing else can explain a difference.

The design includes a paired upper bound, and the pairing is what makes the result
interpretable. The oracle source runs a greedy forward-backward search followed by simulated
annealing over the full edge space. In its leaky form it selects its structure using the
held-out responses it is subsequently scored on, which measures what structure exists in the
grammar rather than what can be found. In its nested form the identical search selects on an
inner split and is scored on an outer split that selection never touched, which measures what
can be found. Same grammar, same budget, same folds; the only difference is access to the answer.

The gap between those two arms is the paper's main quantity. If the leaky arm fails, no honest
procedure in this grammar could succeed and the limitation is the model class. If the leaky arm
succeeds and the nested arm does not, the structure exists and the limitation is discovery. The
second is what happens.

The proposer under test is a large language model asked to propose regulatory edges within a
gene module. That is the practically interesting source, since it is the one being adopted, but
the design does not depend on it: the same comparison bounds every other source in the table.

---

## 2. Data and setup

**Data.** Steady-state single-knockdown Perturb-seq responses from an immune-cell atlas,
processed to per-perturbation response vectors against a pooled control. Perturbations are
grouped into 27 gene modules of three kinds: curated immunological modules, co-response modules
derived from response correlation, and regulon modules derived from annotated targets of a
single transcription factor.

**Evaluation.** Leave-one-perturbation-out. For each fold, the model is fitted without one
perturbation and scored on its held-out response. The primary metric is DE-overlap, the
precision at 50 between the predicted and observed differentially expressed genes, with ties
broken at random (amendment A2). Folds that carry no differentially expressed gene are not
scoreable and are excluded.

**Advantage statistic.** An earlier form of this work established advantage by checking whether
two independently bootstrapped intervals overlapped. That discards the fold-level pairing the
design provides, and perturbations differ enormously in how predictable they are, so the shared
between-fold variance dominates the comparison. Section 1.2 of the pre-registration fixes a
paired form instead: the mean over folds of the within-fold difference, with a paired bootstrap
over folds at 10,000 resamples, and advantage recorded only when the interval's lower bound
clears zero. Across families of modules, Benjamini-Hochberg correction is applied at q = 0.05.

**Structure sources.** S1 language-model proposal; S2 textbook; S3 random null at matched edge
count; S4 algorithmic, being GRNBoost2 and a mean-difference method; S5 oracle, in a leaky and
a nested honest variant; S6 regularized linear map; S7 mean of training perturbations. A zero
predictor is included as a floor.

**A necessary control.** Because a structural prediction is identically zero for any
perturbation of a gene with no outgoing edges, a comparison against a source that cannot
predict a fold is not a comparison at all. Every table below is read against the
random-structure null at matched edge count, and any module where the linear arm fails to beat
that null is treated as unable to support a comparison rather than as evidence.

---

## 3. No structure source beats a linear map where the baseline is sound

### 3.1 The primary module

The cytokine production module carries 28 genes, 28 perturbations, 109 differentially expressed
entries at FDR 0.10, and 19 folds carrying a differentially expressed gene.

| Source | Edges | DE-overlap [95% CI] | Delta vs linear [95% CI] | Advantage |
|---|---|---|---|---|
| linear (S6) | n/a | 0.4451 [0.3280, 0.5716] | reference | n/a |
| mean (S7) | n/a | 0.3712 [0.2626, 0.4899] | -0.0739 [-0.2555, +0.1055] | no |
| oracle (S5) | 7 | 0.2920 [0.1646, 0.4464] | -0.1531 [-0.2837, +0.0133] | no |
| mean difference (S4) | 30 | 0.2468 [0.1530, 0.3529] | -0.1983 [-0.3532, -0.0646] | no |
| zero | n/a | 0.1252 [0.0958, 0.1549] | -0.3199 [-0.4549, -0.1967] | no |
| textbook (S2) | 16 | 0.1211 [0.0734, 0.1768] | -0.3241 [-0.4610, -0.2020] | no |
| GRNBoost2 (S4) | 30 | 0.1113 [0.0767, 0.1449] | -0.3338 [-0.4717, -0.2080] | no |

The random-structure null over 1000 draws at the oracle's edge count is 0.1082.

Four properties make this a ceiling rather than a single unlucky run. Ten independent oracle
searches returned 0.1690 to 0.2920 with a standard deviation of 0.038, so the ceiling does not
move with the seed. The nested honest variant, selected on an inner split and scored on an
outer split that selection never touched, reaches 0.1333 [0.0797, 0.1828], barely above the
random null; the gap to the leaky 0.2920 measures the selection advantage directly. At 0.2920
the leaky oracle is 2.7 times the null mean, so the search is finding real signal and is still
far below a linear map. Textbook immunology compiled into this grammar reaches 0.1211 against a
null of 0.1082, which is to say it predicts held-out responses no better than a random
structure of the same size.

### 3.2 At scale: a reachable ceiling that nothing honest reaches

Across 27 modules, 14 have a linear arm that does not beat its own random-structure null. On
those the linear baseline is performing at chance, and a structural advantage over it carries no
information, so they cannot support the comparison. A thirteenth module drops out once its null
is estimated from twenty searches rather than five, its linear arm falling from just above the
null's 95th percentile to just below. **Twelve modules remain.**

**The search converges, which has to be shown rather than assumed.** The oracle ceiling is the
maximum over independent searches, and a maximum over k searches estimates a supremum only once
adding searches stops raising it. For each module and each k, the mean over all seed subsets of
size k of the maximum within the subset gives a saturation curve. At twenty seeds every module
has flattened: the twentieth seed raises the estimate by a median of 0.0008, and by at most
0.0045, buying a median 1.1 percent of the whole climb from a single seed. One module illustrates
why the check matters. Its curve reads 0.2373, 0.3121, 0.3619 at one to three seeds and 0.4177 at
five, an almost straight line that suggested no convergence at all; it then reaches 0.4533 at ten
and 0.4558 at both fifteen and twenty. A five-seed window had sat in the linear stretch of a
curve that flattens later.

**The ceiling sits above the linear map.**

| | Exceeds the linear map | Exceeds its own random null |
|---|---|---|
| oracle, selection sees the scored folds | **11 of 12** | 12 of 12 |
| oracle, nested so selection never sees them | **5 of 12** | 7 of 12 |

Under the paired interval rule the leaky arm clears the linear map on 5 of 12 and the nested arm
on none. The median gap between the two arms is 0.2077.

This is the paper's central measurement and it separates two questions the usual framing runs
together. Whether a structure exists in this grammar that predicts held-out responses better than
a linear map: yes, on almost every module. Whether any procedure that is not allowed to see the
answer finds it: no. The leaky arm and the nested arm are the same search under the same budget
on the same folds, differing only in whether selection may use the responses it is scored on, so
the gap between them is attributable to that and nothing else.

Every other source falls further behind the nested arm. On no sound module does a language-model
proposal, a textbook structure, GRNBoost2 or a mean-difference method beat the linear map.

**A note on summary statistics.** Reporting the ceiling as a median over searches rather than a
maximum gives 0 of 12 instead of 5. The maximum is the right statistic here because the claim
concerns what structure exists rather than what a typical search returns, and it is defensible
only because convergence was checked. An earlier version of this work reported the median-based
figure from a five-seed run, which understated the ceiling on exactly the modules where the
search had not yet converged; section 9 records the correction.

---

## 4. Why: the response matrix is low rank rather than sparse causal

If a sparse causal structure generated this data, a model in that class should be identifiable
from it. The identifiability diagnostics say it is not, and the reason is measurable before any
model is fitted.

**The diagnostics.** On the column-centred perturbation-by-gene response matrix, the effective
rank is computed as the participation ratio of the covariance spectrum, alongside the fraction
of variance in the leading principal component. The cytokine module gives an effective rank of
3.639 out of 28 perturbations and a leading-component fraction of 0.2758.

**The calibration test.** Structures were drawn from candidate generators, simulated through
the same fitting and simulation path as the real data, and their diagnostics compared against
the real module. Acceptance was fixed in advance at 15 percent relative on effective rank and
0.05 absolute on the leading component. All families were drawn at 48 edges.

| Generator | Effective rank | Relative error | Leading PC | Absolute error | Verdict |
|---|---|---|---|---|---|
| real module | 3.639 | | 0.2758 | | |
| uniform | 10.274 | 182.4% | 0.0837 | 0.1921 | reject |
| hub | 12.046 | 231.1% | 0.0445 | 0.2313 | reject |
| scale-free | 12.433 | 241.7% | 0.0283 | 0.2474 | reject |
| modular | 12.764 | 250.8% | 0.0284 | 0.2473 | reject |
| low-rank surrogate, rank 6 | 3.578 | 1.7% | 0.2749 | 0.0008 | accept |

No causal generator reproduces the signature at any density from 12 to 84 edges, or in any of
the four topology families. A structureless low-rank surrogate reproduces it closely.

**The topology families were added to challenge this claim.** The expectation was that
concentrating out-degree would produce a low-rank response matrix, since a few master regulators
driving most targets is the standard picture. On the primary module all three concentrated
families land further from the real module than the uniform generator, raising the effective
rank from 10.3 to between 12.0 and 12.8.

That ordering is module-specific and no general claim is made from it. Replicated across six
further modules, the direction is mixed on five and reversed on one, where all three
concentrated families sit closer to the real module than uniform. What holds everywhere is the
weaker and sufficient statement: every topology family rejects on every module tested, so the
signature is not reachable by concentrating out-degree, and which topology fails worst is not a
stable property.

**The sweep was replicated on six modules chosen before it ran.** The core result holds on four
of the six: no causal generator meets the acceptance criterion and the low-rank surrogate does.
On the remaining two, TCR_signalosome and coresponse_PIM1, a structural generator accepts, in
both cases only the variant carrying a shared component at four to eight times the structural
weight. A purely structural generator accepts nowhere, on any of the seven modules, at any
density or topology.

The claim is therefore narrowed to what the measurement supports. Sparse causal structure alone
does not reproduce the identifiability signature of this data. Sparse structure plus a dominant
shared component sometimes does, and that shared component is precisely what the grammar cannot
express and what section 5.2 supplies as a fitted offset.

**The diagnostics are not a thresholding artefact.** They are computed on the response matrix
rather than on a differential-expression call, so the differential-expression threshold cannot
reach them, and that matrix is 96.4 percent non-zero. An effective rank of 3.64 on an almost
fully populated 28 by 28 matrix is not a statement about sparsity.

**A prediction that did not hold.** If the response matrix is effectively low rank, then
truncating the ridge map near that rank should cost little held-out performance. That was
registered as a directional prediction with a rule fixed in advance, and on the 13 sound
modules it is refuted: the median ratio of reduced-rank to full-rank held-out DE-overlap is
0.865 against a pre-specified floor of 0.90.

The failure is informative rather than fatal. The ridge map does not clear the reduced-rank map
on a single one of the 13 modules, so no truncation loss is statistically established anywhere;
the refutation rests on a median of point estimates, which was the weaker of the two
pre-specified criteria. The modules that fail are the 40-perturbation regulon modules, which do
not recover the ridge's score until rank 24 to 32, while the modules the low-rank account was
built on behave as predicted, with the primary module at a ratio of 1.038 by rank 4.

The reason the bridge does not hold is worth stating, because it constrains how the diagnostic
should be read. The participation ratio measures how concentrated the variance is, while
DE-overlap ranks the top 50 genes, and a component carrying little variance can still change
which genes enter that ranking. The effective rank is therefore not the number of components
needed to predict, and this paper does not claim it is.

**Consistency with the discovery result.** This is what makes section 3.2's gap expected rather
than surprising: if many structures fit the training data equally well while predicting held-out
responses very differently, then a search scored on training loss has no way to choose among
them, and only a search allowed to see the held-out answer can. Within the near-optimal class of
structures, 107 lie within 5 percent of the best training loss, and their held-out scores span
0.188 to 0.304.
Fitting the training data equally well carries almost no information about predicting held-out
responses. Across that class, 84 percent of the 50 distinct edges have a determined sign, and
one edge appears activating in one search and repressing in another. The structure is not
identified by the data, which is what the diagnostics predicted.

### 4.3 The signature holds across cell state and does not hold on a second atlas

If the account above is right, the diagnostics should be a stable property of the data rather
than of one processing choice, and they should be measurable anywhere. Both were checked, and
the second check limits the claim.

**Across cell state, within this atlas, the signature is stable.** The same 27 modules measured
in all three states give a median normalised effective rank of 0.196 at rest, 0.182 at eight
hours of stimulation and 0.184 at forty-eight. Cell state does not move it.

**Across atlases it is not.** Comparing whole matrices would say nothing, since the modules here
are 11 to 40 genes while a combinatorial CRISPRa atlas is 105 perturbations by 20,421 genes, and
effective rank is not comparable across shapes. Submatrices were therefore drawn from that atlas
at the same shapes, with the readout set equal to the perturbed set exactly as the modules here
are constructed, 200 draws at each size.

| Shape | This atlas, normalised effective rank | Leading PC | Second atlas, normalised effective rank | Leading PC |
|---|---|---|---|---|
| 11 by 11 | 0.145 | 0.654 | 0.443 | 0.066 |
| 20 by 20 | 0.190 | 0.208 | 0.407 | 0.055 |
| 28 by 28 | 0.130 | 0.276 | 0.391 | 0.045 |
| 40 by 40 | 0.145 | 0.159 | 0.368 | 0.038 |

At every matched shape the second atlas has two to three times the normalised effective rank and
roughly a quarter of the leading-component fraction. Its response matrix is not low rank in the
sense that matters here.

This bounds the claim rather than breaking it. The account in this section explains why sparse
structure does not help *on this data*, it is measured rather than assumed, and it survives the
one within-dataset generalisation available. It is not a general law about perturbation
response matrices, and the diagnostic's value is that it can be computed on any new dataset
before a model is fitted, not that its value is known in advance. A plausible mechanism for the
difference, offered as interpretation and not as measurement, is that these are knockdowns,
whose responses collapse onto a shared stress-like axis, while the second atlas activates
transcription factors, each of which drives a distinct programme.

---

## 5. Three objections, answered by measurement

### 5.1 That the metric is one sparse models forfeit by construction

A structural prediction is identically zero for any perturbation of a gene with no outgoing
edges, so a sparse structure can lose folds it never had a chance at. The forfeit is real and
was quantified: searched structures cover a median of 0.49 of scoreable folds.

Two measurements answer the objection. First, conditioning the comparison on covered folds
only: the language model's own structures cover 95 and 100 percent of folds on the two curated
modules that carry them, and still lose to the linear map by 0.21 and 0.13 on the folds where
they do predict. The oracle's covered-fold gain is not evidence, because the oracle selected
the folds it covers.

Second, removing the shared component from the target outright. For each fold, the top-k right
singular vectors of the training responses are computed and both the observation and the
prediction are projected onto the orthogonal complement, so that what remains is the
perturbation-specific part of the response, the part the model class is supposed to be good at.
Seven modules were scored at the pre-specified k of 1.

| Module | Folds | DE mass removed | Residual null | Linear | Best structural | Delta | Status |
|---|---|---|---|---|---|---|---|
| Cytokine_production | 19 | 29.1% | 0.3787 | 0.5331 | oracle 0.4500 | -0.0831 | informative |
| TCR_signalosome | 7 | 48.5% | 0.4595 | 0.6168 | GRNBoost2 0.6276 | +0.0108 | informative |
| CD4_lineage_TFs | 11 | 15.4% | 0.2260 | 0.2701 | mean difference 0.3862 | +0.1160 | informative |
| coresponse_PIM1 | 18 | 20.4% | 0.3418 | 0.3468 | oracle 0.3617 | +0.0149 | informative, marginal |
| coresponse_HCCS | 17 | 3.5% | 0.1398 | 0.2157 | oracle 0.2184 | +0.0027 | does not test it |
| coresponse_MOV10 | 16 | 21.6% | 0.2113 | 0.1532 | oracle 0.3040 | +0.1508 | degenerate |
| Th2_GATA3 | 2 | 29.3% | 0.2798 | 0.5000 | oracle 0.4167 | -0.0833 | too few folds |

Three modules are excluded, each for a stated reason. HCCS lost 3.5 percent of its
differential-expression mass to the projection, so its residual is nearly the full response and
the comparison repeats the primary metric under another name. MOV10 has a linear arm at 0.1532
against a residual null of 0.2113, so its apparent structural gain of +0.1508 is measured
against a comparator performing below chance, which is the same failure the cross-tabulation of
section 3.2 was built to catch, reappearing in a new metric. Th2 carries two folds.

On the four modules where the measurement is sound, no structural source clears the linear map.
The cytokine module is negative. The other three are positive with every interval spanning
zero, and PIM1 is marginal twice over, since its linear arm beats its residual null by only
0.0050.

One methodological point came out of this sweep rather than from the hypothesis. A single
leading component removes between 3.5 and 48.5 percent of differential-expression mass
depending on the module, so the residual is not one quantity across modules and no score on it
should be pooled without that fraction reported alongside.

### 5.2 That the model class rather than the structure is the constraint

If the structural grammar simply cannot express a response shared across perturbations, the
comparison is unfair to it. Handing the structural model that shared response as a fitted
offset tests this. On the primary module it moves the oracle from 0.2920 to 0.4820, against the
linear map's 0.4451, so the grammar can reach a ridge map once it is given both a workable
structure and the component it could not express.

That is a real finding about the model class, and it does not rescue the proposer. The same
intervention is worth 0.0085 to the language model's own structure, and it costs between 0.017
and 0.087 on modules where structure already beats the mean baseline. The effect is conditional
on the module, and it helps where the shared response dominates. The binding constraint is
structure discovery.

### 5.3 That a proposer shown only gene names is not being tested on anything interesting

This objection is partly correct, and it is the one that required new work. The proposal step
was built blind to the data: it saw gene names and a context sentence. A negative result about
such a proposer is partly a property of the harness rather than of the model.

An arm was added in which the proposer is shown the response data. Across ten modules at five
seeds, 48 of 49 structure-agreement pairs between the data-blind and data-aware arms fall below
their own module's within-arm replicate median, while permuting which response belongs to which
perturbation leaves the proposal inside that replicate range. The proposer conditions on the
data, and on the correspondence between perturbation and response rather than on the presence
of numbers.

The extension was stratified, with the strata analysed separately and never pooled, and it
carried a prediction recorded before the run: that the advantage would be larger where the
data-blind arm has less to work with. Held-out, the generated stratum, where the blind arm sees
only a generic context sentence, averages +35.9 percent with five of six modules gaining, while
the curated stratum, which carries real biological context, averages +11.5 percent with one of
three. The prediction is confirmed and the gap widened rather than narrowed as the last module
landed.

This does not make the proposals competitive. The best data-aware score on any curated module
remains well below that module's linear baseline.

---

## 6. The combinatorial regime

Single-perturbation data cannot exercise non-additivity, which is the one regime with genuine
headroom over an additive baseline in principle. A separate pre-registered test on a
combinatorial CRISPRa atlas asks whether a structural model informed only by the single
perturbations predicts held-out doubles better than an additive baseline, on the pairs where
non-additivity is real.

It does not. On the non-additive tertile the structural model reaches 0.349 [0.303, 0.396]
against a fitted-additive oracle at 0.376 and a mean-of-singles at 0.370; on the additive
control it reaches 0.586 against 0.642 and 0.613. A positive control on the top decile of
non-additivity, the 13 strongest interactions, has the model beating both baselines on 0 of 13
pairs.

The headroom is real: even the fitted-additive oracle, which is allowed to see the double,
reaches only 0.38 on the non-additive pairs against 0.64 on the additive controls. But a model
informed only by the singles cannot reach it, because the pair-specific interaction is not
identifiable from single-perturbation marginals. A saturating logic gate imposes one particular
non-additive extrapolation, which is not the true interaction.

### 6.1 A graph-prior method does not beat an additive baseline either

The structural grammar used throughout this paper has no way to represent an external
gene-gene graph, and the obvious objection to section 3 is that a method carrying one would
succeed where it fails. GEARS is that method for this atlas. It was trained for 20 epochs on a
GPU, reaching a validation mean squared error of 0.0030.

Every arm is computed inside GEARS' own processed atlas and its own simulation split, because
that release is not the pseudobulk the rest of section 6 uses; the comparison is internally
consistent rather than pooled with the table above.

| Set | n | Non-additivity | GEARS | Fitted-additive | Mean-of-singles | GEARS minus fitted-additive |
|---|---|---|---|---|---|---|
| additive control | 23 | 0.090 | 0.5130 | 0.7939 | 0.7696 | -0.2809 [-0.3409, -0.2243] |
| non-additive | 23 | 0.292 | 0.3478 | 0.6304 | 0.6148 | -0.2826 [-0.3217, -0.2400] |
| all test doubles | 70 | 0.181 | 0.4457 | 0.7117 | 0.6949 | -0.2660 [-0.2966, -0.2346] |

The direction agrees with the published finding that additive baselines beat this method on
this atlas, which is a check on the implementation rather than a new result.

**It also loses on the metric it optimises.** Scoring only on DE-overlap would compare the
method against a target it never trains toward, so it was rescored on mean squared error, and
retrained for 60 epochs, across all 70 test doubles.

| Metric | GEARS | Fitted-additive | Mean-of-singles |
|---|---|---|---|
| DE-overlap, higher is better | 0.4426 | 0.7117 | 0.6949 |
| MSE, all genes, lower is better | 0.00504 | 0.00123 | 0.00275 |
| MSE, top-20 DE genes, lower is better | 0.28169 | 0.04593 | 0.21405 |

Against the fitted-additive arm the paired intervals are -0.2691 [-0.2989, -0.2377] on
DE-overlap, -0.00381 [-0.00458, -0.00314] on overall MSE and -0.23576 [-0.28268, -0.19298] on
top-20 MSE. The fitted-additive arm sees the double it is scored on; the mean-of-singles arm does
not, has no fitted parameters at all, and still wins on every metric. Validation error bottoms
near epoch 20 and drifts upward by epoch 42 while training error stays flat, so more training was
not what was missing.

The remaining caveat is hyperparameter tuning, which was not attempted, so the statement is
scoped to an untuned model at default settings. What is closed is the objection that the gap was
an artefact of scoring on the wrong metric.

**One consequence matters for the argument as a whole.** This is the atlas that does *not* show
the low-rank signature, by section 4.3, and an additive baseline still wins on it. So the
difficulty of beating simple baselines is not downstream of low rank. The two findings are
separate: section 4 explains section 3 on this project's atlas, and it does not explain this
result.

### 6.2 A widely reused genetic-interaction stratification is not reproducible from the published record

This began as a blocked step and is reported because the blockage is the finding.

The pre-registered analysis referred to the published genetic-interaction subtypes of the
source atlas, the categories named synergy, suppression, redundancy, neomorphism and epistasis.
Those categories are reused downstream, including by later methods papers that benchmark on
this dataset. They are not recoverable from the published record.

**What was checked.** The categories are not in the GEO deposit, whose cell-identity table
carries only barcode, guide identity, read and UMI counts, coverage and cell number. They are
not in either loader of the standard perturbation-data package, nor anywhere in the
distribution of the best-known method trained on this atlas, which ships no data files at all.
They are not in the supplementary material: a byte-level search of every XML part of all nine
supplementary workbooks for the five category names returns zero matches, and the supplementary
text confirms those nine are the complete set. The supplementary PDF mentions the concepts only
in prose and never as an assignment per gene pair.

**What is published is the feature set, not the labels.** One supplementary table gives 125
double perturbations with sixteen continuous columns: the fitness interaction score, the
Theil-Sen coefficients for each single, their magnitude and asymmetry, differential-expression
counts, and a family of distance correlations. Every column is numeric and none is categorical.

**The categories also cannot be regenerated from a published rule.** The supplementary methods
describe the grouping as an adaptation of a dimensionality-reduction method using UMAP over
those features, in which the number of clusters was fixed by qualitative judgement, in the
authors' description a tradeoff between interpretability and granularity, and the projection
shown was selected from 10,000 random-seed iterates. A grouping produced by a seed-dependent
embedding with a hand-chosen cluster count has no deterministic procedure behind it, so it
cannot be reconstructed from the published data even in principle. The one fully specified
numeric criterion in the methods, an asymmetry threshold combined with a minimum interaction
score, orients arrows in a figure rather than assigning subtypes, and applied to the published
table it selects 1 of 125 doubles.

**Why this matters beyond one analysis.** A categorical stratification that is cited and reused
is, in practice, whatever each downstream group re-derives. Counts quoted for these subtypes in
the literature do not match any published table, and re-derivations following the same prose
description need not agree with one another or with the original figure. The general point is
that a stratification is part of a paper's evidence, and publishing the features it was derived
from is not the same as publishing the stratification.

**Consequence here, stated plainly.** The analysis in section 6 cannot be run against the
published assignments, so its pair selection uses a non-additivity measure recomputed from the
pseudobulk. That is a deviation from the pre-registration and it is recorded as one. No
alternative label source was substituted, because deriving categories from the published
coefficients would reproduce exactly the problem described above. The published quantitative
table would allow the selection to be re-run on published quantities rather than a self-derived
proxy, and that is open work rather than a result claimed here.

---

## 7. What was not established

**Whether module properties predict where the ceiling sits was not tested.** Every
cross-validated R squared for the ceiling-versus-diagnostics regression is negative or
indistinguishable from zero, and an intercept-only model predicts held-out ceilings better than
the fitted one. A prospective check on five modules held out before the regression was fitted
passes, but it passes because the prediction intervals are wide enough to be honest about how
little the model knows, not because the point predictions are good. A power calculation shows
that a slope of the size actually estimated would have been missed more often than not by this
design. The correct statement is untested, not refuted, and the regime map is reported as a
limitation. Testing it properly needs a different design: more seeds on the ceiling, the margin
over the module null rather than over the linear map as the target, and roughly double the
modules.

**Three claims from the original plan were withdrawn during execution.** That what predicts is
anti-correlated with what is annotated: this does not survive scaling to 30 modules, where the
within-module rank correlation averages -0.050 and only 2 modules carry enough annotated edges
to test it. That the identifiability diagnostics predict the ceiling: see above. That the
proposal is name-driven in a way that says something general about language models: the
proposal step was blind to the data, so the original finding was partly a property of the
harness, and section 5.3 addresses the harder question directly.

---

## 8. Limitations

**The comparison itself rests on one atlas and one cell type.** Sections 3, 4 and 5 are all
measured on a single steady-state atlas. The identifiability diagnostics have now been computed
on a second one, and they do not hold there (section 4.3), so the low-rank account is
established for this data rather than for perturbation data generally. What was not repeated on
a second atlas is the structure-source comparison itself, which is the more expensive half:
running the full seven-source table elsewhere is the natural next study and would establish
whether claim 1 travels even where claim 2's explanation does not apply.

**The graph-prior comparator is run but not tuned.** GEARS was trained for 20 epochs on default
hyperparameters and scored on this paper's ranking metric, which it does not optimise, having
been trained on mean squared error. It loses to a fitted-additive baseline by 0.27 (section
6.1), in the same direction as the published result for that atlas. That is enough to answer
the objection that a method with external pairwise priors would obviously succeed, and it is
not enough to characterise the method's ceiling. A tuned comparison, and one scored on the
metric the method optimises, is not attempted here.

**Steady state.** Nothing here bears on time-resolved data, where a structural model has
dynamics to exploit and the identifiability picture may differ.

**The oracle bounds this grammar, not all structural models.** The upper bound is an upper
bound within the bounded grammar and its optimizer budget. Annealing gains measured between
0.000 and 0.017 indicate the candidate pool was not silently capping the ceiling, but a
different model class is a different question.

**A retention gap, closed for the modules that matter and open elsewhere.** Per-seed per-fold
scores were not stored during the original comparator runs. They are now retained, and the 13
modules the headline count depends on have been re-run with them, which is what makes the seed
median count in section 3.2 exact rather than reconstructed. The other 14 modules keep the gap,
as do the per-pair names in the combinatorial analysis.

**The ceiling is a supremum over one grammar and one search budget.** Convergence is
established within that budget: twenty seeds, a 120-pair candidate pool ranked by marginal
association, and simulated annealing drawing moves from the whole edge space. Annealing gains of
0.000 to 0.017 indicate the pool was not capping the result, and the saturation curves indicate
the seed count was not either. Neither rules out that a different grammar, or a search of a
different kind, would find more. The claim is bounded to what was searched.

**The honest arm is one honest procedure, not all of them.** The nested oracle is the strongest
structure-discovery method tested here and it is still a greedy search with annealing under a
particular split. That it fails to reach the ceiling is evidence about this class of procedure
on this data, not a proof that no procedure could. What makes it informative is that it is the
same search as the leaky arm under the same budget, so the comparison isolates access to the
answer rather than confounding method with capability.

---

## 9. Pre-registration and deviations

The study was pre-registered before the comparator runs and amended twenty times, each
amendment dated and recorded before the work it governs. The amendments that matter to the
reading of this draft are: the paired advantage statistic replacing an interval-overlap check;
the random tie-break in the differential-expression overlap metric; the substituted
non-additivity measure in section 6 and its consequences; the residual decomposition; the
topology generators; the stratified data-aware arm with its advance prediction; the replication
of the generator sweep on six modules named in advance; the external diagnostics; and the
re-run of the thirteen sound modules with per-seed per-fold scores retained.

**One amendment reversed a conclusion this work had already drawn, and the reversal is the
reason the paper's title changed.** A saturation check on the five-seed data showed the oracle
maximum still climbing on 10 of 13 modules, with one module near-linear in seed count. Read on
its own that says the search never converged and no ceiling was estimated, which would have
forced every claim to rescope to the structures actually searched. Rather than rescope on an
inference from a short run, the search was extended to twenty seeds with the reading rule fixed
in advance: a final-seed gain below 0.005 makes the maximum a supremum estimate. All 13 modules
cleared it. The module that had looked worst reaches the same value at fifteen and twenty seeds.
The five-seed window had been sitting in the linear stretch of a curve that flattens later, so
the apparent non-convergence was an artefact of the range examined. The corrected reading raises
the count and moves the paper's central claim from identifiability to discoverability, which is
a stronger result than the one it replaced.

**One amendment carried a directional prediction and it was refuted.** If the response matrix
is effectively low rank, a reduced-rank map truncated near that rank should recover the full
ridge map's held-out score. The rule was fixed before the run: the ridge clears the reduced-rank
map on at most 4 of the 13 sound modules, and the median ratio is at least 0.90. The first held
completely, since the ridge clears on none of the 13, and the second failed at 0.865, so by the
stated rule the prediction is refuted. Two faults in how that rule was written are recorded with
it: the deciding criterion was a median of point estimates carrying no uncertainty, which was
the weaker of the two, and the rule as written left its own stated middle category unreachable.
The result is reported as refuted rather than reinterpreted, and its lesson is in section 4: the
participation ratio measures how concentrated the variance is, while the metric ranks the top
fifty genes, so the effective rank is not the number of components needed to predict.

Corrections made during execution are recorded in full in the project record, including several
where a figure was stated and then revised on further data. The most consequential were the
discovery that eight of ten corrected advantages were vacuous, the withdrawal of the annotation
inversion, the replacement of three-seed held-out figures that proved to be noise, and the
withdrawal of a mechanism for why concentrated network topologies fail, which held on the
primary module and did not replicate.

---

## 10. Reproduction

All computation runs in Docker. Structure sources are interchangeable at a single interface,
and every comparison table is produced by the same scoring path. Scripts, pre-registration,
amendments and recorded outcomes are in the repository; the extracted published tables used in
section 6 are under `results/norman_supp/`.

---

## Open items in this draft

- The structure-source comparison has not been repeated on a second atlas. The diagnostics have
  been, and they do not carry over, which makes the repeat more interesting rather than less.
- The two modules where a structural generator accepts deserve their own short treatment rather
  than a sentence, since they are the closest thing in this study to a positive result for
  sparse structure.
- The five modules where the converged ceiling clears the linear map are reported as a count and
  not examined. What structure the search found there, and whether any of it is annotated, is the
  obvious next question and the one most likely to say something biological.
- The nested honest arm exceeds its own random null on only 7 of 12 modules. Whether that is the
  split costing it power or the search genuinely finding little is not separated here.
- Section 8's GEARS gap is the most likely reviewer objection and has no measurement behind it.
- Section 6's pair selection can be re-run on published quantities and has not been.
- Figures are not yet drawn. The candidates are the calibration table as a two-panel scatter,
  the coverage-conditioned comparison, and the residual sweep with differential-expression mass
  removed on one axis.
