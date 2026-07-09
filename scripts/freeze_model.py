"""Freeze the model on the training state and save it for evaluation.

The loop runs on the training state (Stim8hr) restricted to the split's training
perturbations, so the in-context held-out and every test-state perturbation stay unseen.
The frozen structure, its fitted parameters, and the ensemble are written to JSON so the
Step 5 and Step 6 evaluation loads a fixed model rather than re-running the loop.
"""
from __future__ import annotations

import json
import sys

from mmc.data.splits import build
from mmc.loop.propose import model_hash
from mmc.loop.run import discover


def _params_json(p: dict | None) -> dict | None:
    if p is None:
        return None
    return {"basal": [float(v) for v in p["basal"]],
            "decay": [float(v) for v in p["decay"]],
            "terms": p["terms"]}


def main() -> None:
    module = sys.argv[1] if len(sys.argv) > 1 else "TCR_signalosome"
    max_iters = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    n_starts = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    max_iter = int(sys.argv[4]) if len(sys.argv) > 4 else 150
    out_path = sys.argv[5] if len(sys.argv) > 5 else "/app/frozen_model.json"

    s = build(module)
    result = discover(module, s.train_state, perts=list(s.train_perts),
                      max_iters=max_iters, n_starts=n_starts, max_iter=max_iter)
    top = result.ensemble.best()

    print(f"frozen structure {model_hash(top.spec)}  training loss {top.loss:.4f}")
    print(f"trained on {s.train_state} perturbations {list(s.train_perts)}")
    print("training-fit history:")
    for h in result.history:
        print("  ", h)

    out = {
        "module": module, "state": s.train_state, "train_perts": list(s.train_perts),
        "genes": top.spec.genes,
        "best": {"hash": model_hash(top.spec), "loss": top.loss,
                 "spec": json.loads(top.spec.to_json()),
                 "params": _params_json(top.params)},
        "ensemble": [{"hash": model_hash(m.spec), "loss": m.loss,
                      "spec": json.loads(m.spec.to_json()),
                      "params": _params_json(m.params)} for m in result.ensemble.members],
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"saved {out_path} ({len(result.ensemble.members)} ensemble members)")


if __name__ == "__main__":
    main()
