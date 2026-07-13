"""Capture a replay artifact of the discovery loop for the demo.

Runs the propose, read, repair loop on a module in one condition and writes a per-iteration
replay JSON: the structure at each step, the residuals the reasoning step was shown, its
verbatim rationale, the structural edit, the training loss, and the structural-residual count.
The demo (demo/app.py) replays this file deterministically. Requires MMC_ZHU_STORE and
ANTHROPIC_API_KEY.

    python scripts/capture_loop_trace.py TCR_core Stim8hr results/trace/loop_replay.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from mmc.data import module_extract
from mmc.loop.run import discover


def main() -> None:
    module = sys.argv[1] if len(sys.argv) > 1 else "TCR_core"
    condition = sys.argv[2] if len(sys.argv) > 2 else "Stim8hr"
    out = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("results/trace/loop_replay.json")

    # Preflight: confirm the store actually has this module in this condition before spending
    # any model calls, so a missing slice fails fast with a clear message.
    observed = module_extract.observed_deltas(module, condition)
    if not observed:
        raise SystemExit(f"No observed knockdown deltas for {module} in {condition}; "
                         "check MMC_ZHU_STORE points at the full store.")
    print(f"preflight: {len(observed)} perturbation(s) for {module} in {condition}: "
          f"{sorted(observed)}")

    replay: list = []
    # The loop's rigorous defaults: enough restarts and optimizer iterations for the
    # fit-versus-structure gate to separate a trapped optimizer from a genuine structural error
    # before any residual is handed to the reasoning step, iterating to a residual plateau.
    result = discover(module, condition, replay_log=replay)

    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "module": f"{module}:{condition}",
        "iterations": replay,
        "ensemble_size": len(result.ensemble.members),
    }
    out.write_text(json.dumps(payload, indent=2))

    print(f"wrote {out} with {len(replay)} iteration(s)")
    for r in replay:
        print(f"  iteration {r['n']}: loss {r['fit']}, structural residuals {r['n_structural']}, "
              f"{len(r['edges'])} edges")
        print(f"    rationale: {r['rationale']}")
        if r["edit"]:
            print(f"    edit: {r['edit']}")


if __name__ == "__main__":
    main()
