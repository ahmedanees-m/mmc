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
from ..fit.diagnose import diagnose
from ..shared.trace import Step, Trace
from . import propose as propose_mod
from . import repair as repair_mod
from . import residuals
from .ensemble import Ensemble
from .propose import model_hash

DEFAULT_CONTEXT = {
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
             perts: list[str] | None = None, start_spec=None,
             max_iters: int = 4, n_starts: int = 8, max_iter: int = 60,
             ensemble_k: int = 5, ensemble_margin: float = 0.2,
             trace: Trace | None = None) -> LoopResult:
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
    else:
        spec, _ = propose_mod.propose(genes, context, trace=trace)
    for it in range(max_iters):
        best, labels, stats = diagnose(spec, observed, n_starts=n_starts, max_iter=max_iter)
        ens.add(spec, best["loss"], best["params"])
        items = residuals.structural_items(spec, best, labels, observed)
        trace.log(Step(kind="fit", model_hash=model_hash(spec),
                       rationale=f"{len(items)} structural residual(s)",
                       edit=None, train_score=best["loss"]))
        history.append({"hash": model_hash(spec), "loss": round(best["loss"], 4),
                        "n_structural": len(items)})
        if not items or it == max_iters - 1:
            break
        summary = residuals.summary_text(items, stats)
        spec, _ = repair_mod.repair(spec, context, summary, trace=trace)

    top = ens.best()
    trace.log(Step(kind="freeze",
                   model_hash=model_hash(top.spec) if top else "",
                   rationale=(f"froze an ensemble of {len(ens.members)} structure(s) "
                              f"within {ens.margin} loss of the best; training data only "
                              f"({condition}), no held-out or test-state data seen"),
                   edit=None, train_score=top.loss if top else None))
    return LoopResult(module, condition, ens, trace, history)
