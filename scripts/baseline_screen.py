"""Baseline-only mean-beatability screen (no modeling).

Before building anything, measure how strong the mean and persistence baselines already
are on each candidate module. A module where the mean predicts held-out perturbations
almost perfectly has no headroom for any method to add value (the Ahlmann-Eltze regime);
a module whose perturbations are diverse leaves the mean weak, and that is where a
mechanistic model can win. For each candidate and the transfer direction, report:

  perts        the number of regulator knockdowns measured in each state (power)
  loo-mean     leave-one-out mean baseline in the train state (the in-context ceiling)
  persistence  train-state effect used as the test-state prediction (transfer ceiling)
  xfer-mean    train mean used as the test-state prediction
  redundancy   mean pairwise correlation of the perturbation responses (1 = identical)

The number that decides the ceiling is loo-mean pooled Pearson: if it is near 1, no model
can win the in-context number on that module.
"""
from __future__ import annotations

import numpy as np

from mmc.baselines import linear as linear_bl
from mmc.baselines import mean as mean_bl
from mmc.eval import metrics
from mmc.shared import store

TRAIN, TEST = "Stim8hr", "Stim48hr"
DE = 0.5

# Candidate modules, defined here so the screen stays independent of the pre-registered
# set. The screen filters each to the regulators that actually have knockdown data.
CANDIDATES: dict[str, dict[str, list[str]]] = {
    "TCR_signalosome": {
        "regulators": ["CD3E", "ZAP70", "LAT", "LCP2", "PLCG1", "PRKCQ"],
        "targets": ["ZAP70", "LAT", "LCP2", "PLCG1", "PRKCQ",
                    "IL2", "NFKB1", "RELA", "FOS", "JUN"],
    },
    "Th2_GATA3": {
        "regulators": ["GATA3", "STAT6", "TBX21", "STAT4"],
        "targets": ["GATA3", "STAT6", "TBX21", "STAT4", "IL4", "IL5", "IL13", "IFNG"],
    },
    "CD4_lineage_TFs": {
        "regulators": ["GATA3", "TBX21", "STAT6", "STAT4", "STAT1", "STAT3", "RORC",
                       "FOXP3", "IRF4", "BATF", "BCL6", "MAF", "PRDM1", "RUNX3", "AHR",
                       "IKZF2", "FOXO1", "TCF7", "LEF1", "ID2"],
        "targets": ["GATA3", "TBX21", "STAT6", "STAT4", "STAT1", "STAT3", "RORC",
                    "FOXP3", "IRF4", "BATF", "IL4", "IL5", "IL13", "IFNG", "IL17A",
                    "IL17F", "IL2", "IL10", "IL21", "IL9", "TNF", "CSF2", "CTLA4"],
    },
    "CD4_TF_network": {
        "regulators": ["GATA3", "TBX21", "STAT6", "STAT4", "STAT1", "STAT3", "RORC",
                       "FOXP3", "IRF4", "BATF", "BCL6", "PRDM1", "RUNX3", "MAF"],
        "targets": ["GATA3", "TBX21", "STAT6", "STAT4", "STAT1", "STAT3", "RORC",
                    "FOXP3", "IRF4", "BATF", "BCL6", "PRDM1", "RUNX3", "MAF"],
    },
}


def observed(regulators: list[str], targets: list[str], condition: str):
    genes = list(dict.fromkeys(regulators + targets))
    regs = [g for g in regulators if g in genes]
    df = store.module_effects(regs, genes, condition)
    out: dict[str, dict[str, float]] = {}
    for _, r in df.iterrows():
        p, g = r["perturbation"], r["target_gene"]
        if p == g:
            continue
        out.setdefault(p, {})[g] = float(r["effect_size"])
    return out, genes


def loo_mean(obs: dict, genes: list[str]) -> dict:
    perts = list(obs)
    preds, obss = {}, {}
    for p in perts:
        others = {q: obs[q] for q in perts if q != p}
        if not others:
            continue
        preds[p] = mean_bl.fit(others, genes)
        obss[p] = obs[p]
    return metrics.score_set(preds, obss, genes, DE) if preds else {"pooled": None}


def redundancy(obs: dict, genes: list[str]) -> float | None:
    perts = list(obs)
    if len(perts) < 2:
        return None
    M = np.array([[obs[p].get(g, 0.0) for g in genes] for p in perts], dtype=float)
    corrs = []
    for i in range(len(perts)):
        for j in range(i + 1, len(perts)):
            a, b = M[i], M[j]
            if np.ptp(a) == 0 or np.ptp(b) == 0:
                continue
            corrs.append(np.corrcoef(a, b)[0, 1])
    return float(np.mean(corrs)) if corrs else None


def _p(score: dict, key: str = "pearson"):
    pooled = score.get("pooled") if score else None
    return pooled.get(key) if pooled else None


def _fmt(x) -> str:
    return "n/a" if x is None or (isinstance(x, float) and x != x) else f"{x:.3f}"


def screen(name: str, spec: dict) -> dict:
    tr_obs, genes = observed(spec["regulators"], spec["targets"], TRAIN)
    te_obs, _ = observed(spec["regulators"], spec["targets"], TEST)
    common = [p for p in te_obs if p in tr_obs]

    loo = loo_mean(tr_obs, genes)
    # transfer: predict each test-state perturbation from the train state.
    pers = linear_bl.persistence(tr_obs, common)
    xmean = mean_bl.predict(mean_bl.fit(tr_obs, genes), common)
    te_common = {p: te_obs[p] for p in common}
    pers_s = metrics.score_set(pers, te_common, genes, DE) if common else {"pooled": None}
    xmean_s = metrics.score_set(xmean, te_common, genes, DE) if common else {"pooled": None}

    return {
        "module": name, "n_train": len(tr_obs), "n_test": len(te_obs), "n_common": len(common),
        "loo_mean_r": _p(loo), "loo_mean_sign": _p(loo, "sign_accuracy"),
        "persistence_r": _p(pers_s), "persistence_sign": _p(pers_s, "sign_accuracy"),
        "xfer_mean_r": _p(xmean_s), "redundancy": redundancy(tr_obs, genes),
    }


def main() -> None:
    rows = [screen(name, spec) for name, spec in CANDIDATES.items()]
    print(f"{'module':18s} {'perts':>10s} {'loo-mean r':>11s} {'loo sign':>9s} "
          f"{'persist r':>10s} {'xfer-mean':>10s} {'redundancy':>11s}")
    for r in rows:
        perts = f"{r['n_train']}/{r['n_test']}"
        print(f"{r['module']:18s} {perts:>10s} {_fmt(r['loo_mean_r']):>11s} "
              f"{_fmt(r['loo_mean_sign']):>9s} {_fmt(r['persistence_r']):>10s} "
              f"{_fmt(r['xfer_mean_r']):>10s} {_fmt(r['redundancy']):>11s}")
    print("\nheadroom = 1 - loo-mean r: a module with a low loo-mean r is where a model can win.")


if __name__ == "__main__":
    main()
