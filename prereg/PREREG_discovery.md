# MMC Part-4 Discovery Pre-registration (cytokine-production module)

Committed before running the loop on the cytokine module. Per v3 Part 4, a discovery
chosen after seeing the loop's output is not a discovery. This file fixes what counts,
and the decision rule, in advance.

**Committed:** 2026-07-10. **Author:** Anees Ahmed Mahaboob Ali. **Backend:** structural
steady-state (JAX), commit f668781. **Gate:** mmc/eval/holdout.py. **Module builder:**
scripts/cytokine_module.py (commit 22b2a7e).

## 0. Stopping rule (stated first)

This is the first properly-powered test in the favorable TF-to-cytokine regime, and it is
the immune positive's last module. After it:
- Clean win (model greater than linear on held-out DE-overlap, with separated CIs, on a
  powered module) then harvest and validate discovery candidates.
- No win despite adequate power then it is decisive that mechanism fits but does not
  predict on single-knockdown steady-state immune data; the positive centerpiece moves to
  the Norman epistasis capability, and the immune work becomes the method plus the
  limit-map.
- Either way, stop module-hunting. No fourth immune module.

## 1. Hypothesis

**H (discovery):** the loop implicates at least one dark-candidate regulator
(data-supported, non-canonical, not named by Zhu for any program) in cytokine control,
adding mechanistic structure beyond Zhu's correlational DE (a signed edge, an intermediate,
a context-specificity), that is held-out-required, seed-stable, and orthogonally validated.
**H0:** every held-out-required edge is canonical or Zhu-named, or no dark edge is
held-out-required, or none validates. H0 is a reportable outcome.

## 2. Module (scripts/cytokine_module.py)

- **Readout panel (DE columns):** 21 cytokines and chemokines (IL2/3/4/5/9/10/13/21,
  IL17A/F, IL22, IFNG, TNF, LTA, CSF2, IL16, CCL3/4/5, CXCL8, TGFB1). Breadth is power.
- **Known backbone (18):** the strongest canonical TFs and textbook TCR signalosome and
  cytokine receptors by cytokine effect. Scaffold and power; excluded as discovery sources.
- **Dark candidates (15):** regulators with FDR-significant cytokine effects that are not
  canonical, not textbook signaling, and not Zhu-named. Selected within a breadth window (a
  floor removes failed-knockdown artefacts, a top-decile ceiling removes the global-stress
  confound) and ranked by cytokine-regulatory power. The only admissible discovery sources.
- **Edge regime:** regulator-to-cytokine edges (the high-signal regime); a small set of
  known regulator-to-regulator scaffold edges permitted. No dense TF-to-TF module.
- **Not-reported-by-Zhu, operational:** admissible only if not in CANONICAL_REGULATORS,
  KNOWN_SIGNALING, or ZHU_NAMED_REGULATORS. ZHU_NAMED is complete: the SAGA and Mediator
  complexes and the arrayed and bulk-RNAseq validated genes (ATP2A2, CYB5R4, ELOB, MEN1,
  KDM1A, SGF29, MED24) from the abstract and the Zhu analysis repository, plus the top one
  percent (coef_rank at least 0.99, 279 genes) of Zhu's polarization and aging
  regulator-coefficient tables, since a discovery is void if Zhu named it for any program.

## 3. Power precondition (locked)

- Total DE entries: **456** (at least 40 required). **PASS.**
- Dark candidates included: **15**. Module gene count: **54** (21 readouts + 18 backbone +
  15 dark). Th2 failed here at 5 entries; this module is powered.

## 4. Splits, baselines, metrics (the held-out gate)

- **Split:** leave-one-perturbation-out over the module's regulators (mmc/eval/holdout.py).
  The held-out perturbation is never seen at fit time.
- **Baselines:** mean-of-training, **linear** (ridge reconstruction, the real bar), zero.
- **Metrics:** ACC_DEG and **DE-overlap** (primary; sign accuracy saturates), bootstrap CIs.
- **Success = model greater than linear on DE-overlap, held-out, with separated CIs.**
  Beating the mean while tying linear (the Th2 outcome) is not a win; it is the
  Ahlmann-Eltze signature.

## 5. Discovery criteria (all five; cheapest first)

A held-out-required dark edge is a discovery only if it clears all of:
1. **Data-supported** (FDR-significant cytokine effect, adequate effect size).
2. **Held-out-required** (section 6): removing it lowers held-out DE-overlap beyond
   threshold, seed-stably. In-sample fit improvement does not count.
3. **Non-textbook** (not canonical; passes the deleted-knowledge test: the edge is in the
   data pattern, surfaced from the numbers independent of Claude's priors).
4. **Not reported by Zhu** (section 2) and MMC adds mechanistic structure beyond Zhu's DE.
5. **Orthogonally validated.** For a dark edge, validation is Moonen CRE-gene links, TF
   ChIP-seq or motif evidence, or an independent arrayed or FACS screen. Zhu's own arrayed
   validation table validates the genes Zhu named, which are excluded, so a dark candidate
   will be absent from it; that absence is a negative check confirming it is not Zhu's
   finding, not evidence that it is unvalidated. Do not read absence from the Zhu table as
   unvalidated.

A candidate failing any gate is logged as proposed-but-unvalidated, not a discovery.

## 6. Edge necessity (mmc/eval/holdout.py edge_ablation_holdout)

An edge is held-out-required iff, removed from the frozen structure and LOO re-run, the
module's held-out DE-overlap drops by at least **0.10**, and the drop is seed-stable
(seed_stability std below **0.05**). Necessity is measured on held-out prediction, never on
in-sample fit. Th2 showed why this needs power: 5 DE entries gave a uniform, uninformative
ablation.

## 7. Harvest and validation

Freeze the loop's structure on training only; keep the full trace and the engineer-behavior
log (textbook / novel-but-known / novel-and-unreported per proposal). Confirm the held-out
gate; if no win despite power, invoke section 0 and pivot to Norman. If a win, enumerate the
held-out-required dark edges and run each through the five criteria; prefer clean-mechanism
candidates (a kinase such as STK11 or SIK3, a TF cofactor such as CBFB) over chromatin
candidates, which are conceptually close to Zhu's SAGA and Mediator theme and are held to a
higher bar. A candidate clearing all five is the primary result; none clearing all five is a
null result plus the behavior log.

## 8. Anti-theater guardrails

- Pre-registered before the run (this file, committed with a timestamp).
- ZHU_NAMED complete from the paper and repository (cytokine and program).
- Necessity and success on held-out, never in-sample (the corrected Th2 lesson).
- Orthogonal validation required; the trace is a record of reasoning, not evidence in place
  of validation.
- Every claim scoped to this atlas, this module, CD4+ T cells.
- Norman epistasis runs in parallel as the discovery-independent insurance positive.
