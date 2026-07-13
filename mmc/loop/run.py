"""The propose, read, repair, ensemble discovery loop on a module in one state.

The loop proposes an executable structure, fits it, reads the structural residuals
from the diagnostic gate, and repairs the structure, iterating to a residual plateau
or an iteration budget while tracking the ensemble of equally-fitting structures.
Only the training-state data is used; held-out and test-state data never enter the
loop, which the trace records.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..data import module_extract
from ..shared.trace import Step, Trace
from . import propose as propose_mod
from . import repair as repair_mod
from . import residuals
from .ensemble import Ensemble
from .propose import model_hash

DEFAULT_CONTEXT = {
    "CD4_lineage_TFs": (
        "These genes are the lineage-defining transcription factors of CD4+ T-helper "
        "differentiation and their cytokine output in stimulated primary human CD4+ T "
        "cells. GATA3 (Th2), TBX21 (Th1), RORC (Th17), FOXP3 (Treg), BCL6 (Tfh), and the "
        "STATs upstream of them cross-regulate: GATA3 and TBX21 are mutually antagonistic, "
        "STAT6 drives GATA3, STAT4 and STAT1 drive TBX21, STAT3 drives RORC. Their output "
        "is the cytokines IL4, IL5, IL13 (Th2), IFNG (Th1), IL17A, IL17F (Th17), IL10, IL21, "
        "and IL2. Unlike a single signaling cascade, these regulators do different things "
        "when knocked down, some activating and some repressing their targets. Knockdown "
        "effects are log2 fold changes at 8 hours of stimulation."
    ),
    "TCR_signalosome": (
        "These genes are the T-cell-receptor signalosome and its immediate transcriptional "
        "output in stimulated primary human CD4+ T cells. CD3E, ZAP70, LAT, LCP2, PLCG1, and "
        "PRKCQ form the proximal signaling complex; their activity converges on IL2, NF-kB "
        "(NFKB1, RELA), and AP-1 (FOS, JUN). Knockdown effects are log2 fold changes at 8 "
        "hours of stimulation."
    ),
}


@dataclass
class LoopResult:
    module: str
    condition: str
    ensemble: Ensemble
    trace: Trace
    history: list[dict] = field(default_factory=list)


def discover(module: str, condition: str, context: str | None = None, *,
             perts: list[str] | None = None, start_spec=None, backend: str = "structural",
             max_iters: int = 4, n_starts: int = 8, max_iter: int = 300,
             ensemble_k: int = 5, ensemble_margin: float = 0.2,
             trace: Trace | None = None, replay_log: list | None = None) -> LoopResult:
    # The structural steady-state backend (v3) is the default; the ODE backend is kept
    # for the earlier experiments and the tests.
    if backend == "structural":
        from ..fit import fit_structural as fitmod
    else:
        from ..fit import diagnose as fitmod
    trace = trace or Trace(module=f"{module}:{condition}")
    context = context or DEFAULT_CONTEXT.get(module, f"Module {module} in CD4+ T cells.")
    observed = module_extract.observed_deltas(module, condition)   # training state only
    if perts is not None:
        # restrict to the split's visible perturbations, so a held-out knockdown never
        # enters the loop even though it is measured in this state.
        keep = set(perts)
        observed = {p: d for p, d in observed.items() if p in keep}
    genes = module_extract.model_genes(module)
    ens = Ensemble(k=ensemble_k, margin=ensemble_margin)
    history: list[dict] = []

    if start_spec is not None:
        # Tier B adaptation begins from the frozen structure rather than a fresh proposal,
        # so the repairs are the rewiring the discovery subset forces on a fixed prior.
        spec = start_spec
        trace.log(Step(kind="propose", model_hash=model_hash(spec),
                       rationale="adapt from the frozen structure", edit=None,
                       train_score=None))
        step_rationale, step_edit = "adapt from the frozen structure", None
    else:
        spec, step_rationale = propose_mod.propose(genes, context, trace=trace)
        step_edit = None
    for it in range(max_iters):
        best, labels, stats = fitmod.diagnose(spec, observed, n_starts=n_starts, max_iter=max_iter)
        ens.add(spec, best["loss"], best["params"])
        items = residuals.structural_items(spec, best, labels, observed,
                                           residuals_fn=fitmod.residuals)
        trace.log(Step(kind="fit", model_hash=model_hash(spec),
                       rationale=f"{len(items)} structural residual(s)",
                       edit=None, train_score=best["loss"]))
        history.append({"hash": model_hash(spec), "loss": round(best["loss"], 4),
                        "n_structural": len(items)})
        if replay_log is not None:
            # A per-iteration record for the deterministic loop replay: the structure at this
            # step, the residuals the reasoning step was shown, its rationale, the edit that
            # produced this structure, and the training-fit diagnostics.
            replay_log.append({
                "n": it + 1,
                "genes": list(genes),
                "edges": [{"src": e.regulator, "dst": e.target, "sign": int(e.sign)}
                          for e in spec.edges],
                "residuals": items,
                "rationale": step_rationale,
                "edit": step_edit,
                "fit": round(float(best["loss"]), 4),
                "n_structural": len(items),
            })
        if not items or it == max_iters - 1:
            break
        summary = residuals.summary_text(items, stats)
        spec, step_rationale = repair_mod.repair(spec, context, summary, trace=trace)
        step_edit = trace.steps[-1].edit

    top = ens.best()
    trace.log(Step(kind="freeze",
                   model_hash=model_hash(top.spec) if top else "",
                   rationale=(f"froze an ensemble of {len(ens.members)} structure(s) "
                              f"within {ens.margin} loss of the best; training data only "
                              f"({condition}), no held-out or test-state data seen"),
                   edit=None, train_score=top.loss if top else None))
    return LoopResult(module, condition, ens, trace, history)
