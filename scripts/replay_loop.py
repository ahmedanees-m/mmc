"""Replay the captured discovery-loop trace in the terminal.

A deterministic, offline replay of demo/loop_replay.json. For each iteration of a captured
run it prints the structure, the residuals the reasoning step was shown, its verbatim
rationale, the structural edit, and the training loss. This replays a real recorded run; it
does not call the model and adds no result of its own. The pacing between lines is for
readability only.

    python scripts/replay_loop.py
    python scripts/replay_loop.py --pace 1.0
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "demo", "loop_replay.json")

PATTERN = {"missing_effect": "effect missing", "wrong_sign": "sign wrong",
           "spurious_effect": "spurious effect", "magnitude": "magnitude"}


def edit_summary(edit):
    if not edit:
        return None
    removed = {e.rsplit("(", 1)[0]: e for e in edit.get("removed", [])}
    added = {e.rsplit("(", 1)[0]: e for e in edit.get("added", [])}
    parts = []
    for key, val in added.items():
        if key in removed and "(-1)" in val and "(+1)" in removed[key]:
            src, dst = key.split("->")
            parts.append(f"flip {src} -> {dst} to repressive")
        else:
            parts.append(f"add {val}")
    for key, val in removed.items():
        if key not in added:
            parts.append(f"remove {val}")
    return "; ".join(parts) if parts else "structure revised"


def main():
    ap = argparse.ArgumentParser(description="Terminal replay of the captured discovery loop.")
    ap.add_argument("--trace", default=DEFAULT, help="path to loop_replay.json")
    ap.add_argument("--pace", type=float, default=0.8,
                    help="seconds between lines (readability pacing; use 0 for instant)")
    args = ap.parse_args()

    # The captured rationale carries the model's arrow and Greek notation; emit UTF-8 so it
    # prints on terminals whose default encoding cannot (Windows cp1252).
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

    with open(args.trace, encoding="utf-8") as fh:
        data = json.load(fh)
    iters = data["iterations"]

    def line(text=""):
        print(text, flush=True)
        if args.pace:
            time.sleep(args.pace)

    module = data.get("module", "module")
    bar = "=" * 70
    print(bar, flush=True)
    line(f"MMC discovery loop   deterministic replay of a captured run   {module}")
    line("The reasoning step sets structure; the optimizer sets magnitudes; the")
    line("fit-versus-structure gate surfaces the structural residuals.")
    print(bar, flush=True)
    line()

    prev = None
    for it in iters:
        n = it["n"]
        verb = "proposes the initial structure" if it["edit"] is None else "revises the structure"
        line(f"[iteration {n}] the reasoning step {verb} ...")
        line(f"              {len(it['genes'])} genes, {len(it['edges'])} edges")
        line(f"              rationale: {it['rationale']}")
        summary = edit_summary(it["edit"])
        if summary:
            line(f"              structural edit: {summary}")
        line(f"[iteration {n}] fitting (multi-start) and reading residuals ...")
        line(f"              training loss = {it['fit']:.2f}   "
             f"structural residuals = {it['n_structural']}")
        for r in it["residuals"][:5]:
            tag = PATTERN.get(r["pattern"], r["pattern"])
            line(f"                {r['perturbation']} KD -> {r['target']}: "
                 f"model {r['predicted']:+.2f}, observed {r['observed']:+.2f}   [{tag}]")
        if prev is not None:
            line(f"              training loss change: {it['fit'] - prev:+.2f}")
        prev = it["fit"]
        line()

    print(bar, flush=True)
    line("structure frozen. This structure is what the held-out gate evaluates next.")
    line(f"training loss {iters[0]['fit']:.2f} -> {iters[-1]['fit']:.2f}   "
         f"structural residuals {iters[0]['n_structural']} -> {iters[-1]['n_structural']}")
    print(bar, flush=True)


if __name__ == "__main__":
    main()
