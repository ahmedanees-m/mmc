# MMC: Track B Pre-registration (Norman epistasis capability)

*Commit before building the Norman pipeline. Same discipline as `PREREG_discovery.md`: the held-out rule is fixed in advance, and the null is a reportable, unifying outcome, not a failure.*

**Committed:** 2026-07-10 (this commit, before the Norman pipeline is built) · **Author:**
Anees Ahmed Mahaboob Ali · **Backend:** structural steady-state (JAX), commit f668781 (with
the activation operator, §3) · **Gate:** held-out prediction of doubles (mmc/eval/holdout.py).

---

## 0. Role, stated first
This is the **supporting capability pillar**, not the centerpiece and not a rescue for the paper. The paper is already decided and floored (method + limit-map + honest immune negative + engineer behavior). Track B either adds *"mechanism wins where non-additivity is real"* or adds *"the fit-not-predict boundary holds combinatorially too."* **Both strengthen the same honest paper; neither is load-bearing.** Non-immune (K562), the immune/disease story is deliberately **not** claimed for this piece.

## 1. Hypothesis (pre-registered)
On **strongly non-additive** pairs, the logic-gate structural model predicts **held-out double perturbations** better than the fitted-additive and mean-of-singles baselines; on **additive** pairs, it shows **no advantage**. The *contrast* (advantage concentrated on the non-additive pairs, absent on additive ones) is the claim: it shows the sum-of-products logic gates capture *exactly* the non-additivity that additive models structurally cannot, rather than a generic edge.
**H0 (reportable, unifying):** the model does not beat additive baselines on held-out doubles even for non-additive pairs, i.e. it can *represent* epistasis but not *predict* it out-of-sample. This would unify with Track A into a single clean boundary (mechanism fits but does not predict, single **and** combinatorial), a stronger finding than a lone immune negative.

## 2. Dataset
Norman et al. 2019, *Science* (GEO **GSE133344**): K562, **CRISPRa (gain-of-function)**, ~200k cells, single- and dual-gene perturbations of growth/differentiation genes. Non-additivity is quantified in the paper as **deviation from a fitted additive model** (double ≈ c₁·single₁ + c₂·single₂), measured by **distance correlation d**; GI subtypes follow the coefficients (synergy = two large c; suppression = two small; epistasis = asymmetric; **neomorphic = high d**, a near-silent single amplifying another). Their model parameters + d per pair are in the paper's table of GI coefficients.

## 3. Engineering change (verify before the run)
Norman is **activation**, not knockdown. The structural backend's perturbation operator must support **overexpression**: `do(x_g = high)` (clamp the node high / drive its production up), **not** the `do(x_g = 0)` knockdown clamp used for Zhu. Implement and **verify on a toy** (an activation of a known activator raises its targets) before fitting. Everything else in the backend (structural fixed-point, gradient fit, the gate) is unchanged.

## 4. Interacting-pair selection (mechanical, pre-registered)
- **Non-additive set (where mechanism should win):** pairs in the **top tertile of distance
  correlation d** (high deviation from additive), including the synergy and neomorphic
  subtypes.
- **Additive control set (where no advantage is expected):** pairs in the **bottom tertile
  of d** (well-explained by the additive model).
- Both sets fixed from the Norman GI table **before** any modeling. Record n per set.

## 5. Model + baselines
- **Model:** the logic-gate structural model (sum-of-products), structure over the pair's genes + their readout, fit on the singles.
- **Baselines:** (a) **fitted-additive**: the Norman linear model c₁·s₁ + c₂·s₂ (the harder, meaningful bar); (b) **mean-of-singles**: 0.5·s₁ + 0.5·s₂ (the simple additive floor). (Optional foundation comparator: GEARS on the identical split, not required for the capability claim.)

## 6. Held-out design (the compose test)
Fit structure + parameters on the **single** perturbations only; predict the **held-out doubles** (both singles seen, the double unseen). **Leakage rule:** no double is seen at fit time. This tests whether mechanism *composes* singles into the correct non-additive double, the honest analogue of Track A's held-out gate.

## 7. Metrics
DE-overlap (**primary**) + ACC_DEG on the held-out doubles, bootstrap CIs (same harness pattern as Track A). Report per set (non-additive vs additive).

## 8. Success rule (pre-registered)
**Win =** model **> fitted-additive AND > mean-of-singles** on **held-out DE-overlap (separated CIs)** on the **non-additive** set, **AND** the advantage is **specific** (present on non-additive, absent/smaller on additive, the §1 contrast). Beating additive on non-additive pairs while *also* beating it on additive pairs would suggest a generic (non-epistasis) effect and weakens the mechanistic claim; report it as such.

## 9. Honest caution (the Track-A lesson, applied here)
**Expressivity ≠ generalization.** The logic gates *representing* non-additivity does not guarantee *predicting* held-out doubles; that exact gap sank Track A. Run expecting the possibility that it does not clear even here. A pre-registered **null unifies the paper** (fit-not-predict, single and combinatorial) and is reported as a finding, not a failure. Do not relax the rule after seeing results.

## 10. Scope
Capability validation only; K562, CRISPRa, non-immune. Every claim scoped accordingly. This pillar backs the *method*; the disease relevance lives entirely in the immune limit-map (Track A / C).
