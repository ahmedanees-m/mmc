"""Export each module's response matrices so the analyses can run without the store.

For every module in the definitions file and each condition, writes one compressed archive
holding the effect-size matrix, the FDR matrix, and the gene and perturbation labels. A
directory of these is what MMC_MODULE_MATRICES points at, and it is enough to re-run every
analysis on this atlas.

The generated modules are not in the static registry, so each is registered from its gene
list before its matrices are read. `model_genes` returns the ordered union of regulators and
targets, so registering a module's gene list as both reproduces that list exactly, which is
what the comparator was given.

A module with no perturbed gene in a condition has no matrix to write; it is reported as
skipped rather than failing the run, because the screen outcome for such a module is itself
recorded in the definitions file.

    python scripts/export_module_matrices.py --definitions data/module_definitions.json \
        --conditions Rest,Stim8hr,Stim48hr --out data/module_matrices

Reads the store named by MMC_ZHU_STORE.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from mmc.data import module_extract
from mmc.data.module_data import store_matrices


def load_modules(path: Path) -> list[dict]:
    """The module entries, from whichever key the definitions file uses."""
    d = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(d, list):
        entries = d
    else:
        for key in ("modules", "candidates"):
            if isinstance(d.get(key), list):
                entries = d[key]
                break
        else:
            raise SystemExit(f"{path} has no module list; keys are {sorted(d)}")
    missing = [e for e in entries if not e.get("name") or not e.get("genes")]
    if missing:
        raise SystemExit(f"{len(missing)} entries lack a name or a gene list")
    return entries


def modules_from_results(dirs: str) -> list[dict]:
    """Module entries recovered from comparator result files.

    The curated modules are not in the definitions file, which holds the mechanically
    assembled candidates. Each result records the gene list it was run on under every
    structure source in `specs`, so the module is recoverable from its own output. Where a
    module appears in more than one result the gene lists must agree, or it is not one
    module.
    """
    found: dict[str, dict] = {}
    for d in (Path(p) for p in dirs.split(",") if p):
        for f in sorted(d.glob("*.json")):
            r = json.loads(f.read_text(encoding="utf-8"))
            name = r.get("module")
            specs = r.get("specs") or {}
            lists = {tuple(v["genes"]) for v in specs.values()
                     if isinstance(v, dict) and v.get("genes")}
            if not name or len(lists) != 1:
                continue
            genes = list(next(iter(lists)))
            if name in found and found[name]["genes"] != genes:
                raise SystemExit(f"{name} has two different gene lists across results")
            found[name] = {"name": name, "genes": genes,
                           "condition": r.get("condition"), "source": "curated"}
    return list(found.values())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--definitions", default="data/module_definitions.json",
                    help="one or more definitions files, comma separated")
    ap.add_argument("--conditions", default="Rest,Stim8hr,Stim48hr")
    ap.add_argument("--out", default="data/module_matrices")
    ap.add_argument("--only-passing", action="store_true",
                    help="skip the candidates that failed the power screen")
    ap.add_argument("--force", action="store_true",
                    help="rewrite matrices that already exist in the output directory")
    ap.add_argument("--from-results", default="",
                    help="comma-separated result directories to recover further modules "
                         "from, for the curated modules the definitions file omits")
    args = ap.parse_args()

    entries, seen = [], set()
    for part in args.definitions.split(","):
        for e in load_modules(Path(part)):
            if e["name"] not in seen:
                seen.add(e["name"])
                entries.append(e)
    if args.only_passing:
        entries = [e for e in entries
                   if e.get("passed_screen") or e.get("source") == "curated"]
    if args.from_results:
        known = seen
        extra = [e for e in modules_from_results(args.from_results)
                 if e["name"] not in known]
        print(f"{len(extra)} further modules recovered from result files: "
              f"{', '.join(sorted(e['name'] for e in extra))}")
        entries = extra + entries
    conditions = [c for c in args.conditions.split(",") if c]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    print(f"{len(entries)} modules by {len(conditions)} conditions", flush=True)

    written, existing, skipped = 0, 0, []
    for i, entry in enumerate(entries, 1):
        module, genes = entry["name"], entry["genes"]
        # the ordered union of these two is the gene list the comparator used
        module_extract.register_module(module, genes, genes)
        for condition in conditions:
            dest = out / f"{module}__{condition}.npz"
            if dest.exists() and not args.force:
                existing += 1
                continue
            try:
                g, perts, obs, q = store_matrices(module, condition)
            except Exception as exc:                       # noqa: BLE001
                skipped.append(f"{module}/{condition}: {exc}")
                continue
            np.savez_compressed(
                dest,
                genes=np.array(g, dtype=str),
                perts=np.array(perts, dtype=str),
                effect_size=obs,
                fdr=q,
            )
            written += 1
        print(f"  [{i}/{len(entries)}] {module}: {written} written so far", flush=True)

    print(f"\nwrote {written} matrix files to {out}")
    if skipped:
        print(f"{len(skipped)} module-condition pairs had nothing to write:")
        for s in skipped:
            print(f"  {s}")
    if written == 0 and existing == 0:
        print("nothing was written", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
