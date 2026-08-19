"""Re-derive the ceiling on the seed median instead of the seed maximum.

Review 2 item 4: the ceiling reported throughout is a maximum over search seeds, which is
an upward-biased order statistic whose bias grows with the seed count. Step 1 used ten
seeds and Step 7 used three, so the two are not on one scale. The median is the stable
summary and should be primary, with the maximum reported as an explicit upper bound.

What this can and cannot do is bounded by what the runs stored. Per-seed fold scores were
not recorded, only a per-seed summary with its own bootstrap interval, so the paired
permutation test of Family A cannot be recomputed exactly on the median seed without
re-running. What is available is each seed's own mean and interval, which supports a
comparison against linear and against the module null, and that is what is reported here.
"""
import glob
import json
import os
import statistics as st

RES = "G:/My Drive/VERDICT/MMC/run_archive/results"
AUDIT = "G:/My Drive/VERDICT/mmc-repo/results/step7_branch_audit.json"

audit = {r["module"]: r for r in json.load(open(AUDIT))["rows"]}

rows = []
for pattern in (RES + "/step7/step1_*.json", RES + "/step1_*_full.json",
                RES + "/step1_CD4_lineage_TFs.json"):
    for f in sorted(glob.glob(pattern)):
        d = json.load(open(f))
        if not (d.get("table") and d.get("random_null")):
            continue
        by = [v for v in ((d.get("oracle_ceiling_spread") or {}).get("by_seed") or [])
              if isinstance(v, (int, float))]
        if not by:
            continue
        tab = {r["source"]: r for r in d["table"]}
        name = d["module"]
        a = audit.get(name, {})
        rows.append({
            "module": name,
            "n_seeds": len(by),
            "max": max(by),
            "median": st.median(by),
            "mean": st.mean(by),
            "linear": tab["linear"]["de_overlap_mean"],
            "null_mean": d["random_null"]["mean"],
            "null_p95": a.get("null_p95"),
            "linear_beats_null": a.get("linear_beats_null"),
            "bh_reject_on_max": a.get("oracle_clears_linear"),
        })

rows.sort(key=lambda r: -r["max"])
print("Ceiling on the seed maximum against the seed median, %d modules\n" % len(rows))
print("%-24s %5s %8s %8s %8s %8s %9s" % ("module", "seeds", "max", "median", "drop",
                                          "linear", "null p95"))
drops = []
for r in rows:
    drop = r["median"] - r["max"]
    drops.append(drop)
    print("%-24s %5d %8.4f %8.4f %+8.4f %8.4f %9s"
          % (r["module"], r["n_seeds"], r["max"], r["median"], drop, r["linear"],
             ("%.4f" % r["null_p95"]) if r["null_p95"] is not None else "-"))

print("\nmedian ceiling is lower than the max on %d of %d modules; mean drop %.4f, worst %.4f"
      % (sum(1 for d in drops if d < 0), len(drops), st.mean(drops), min(drops)))

# how the cross-tab moves when the ceiling is the median seed
print("\nCross-tab of section 2.22, recomputed on the median-seed ceiling.")
print("A module counts as the oracle clearing linear only if the median-seed ceiling")
print("exceeds linear at all; this is weaker than the paired test and is an upper bound")
print("on how many would survive it.\n")
cells = {}
for r in rows:
    if r["linear_beats_null"] is None:
        continue
    clears_max = bool(r["bh_reject_on_max"])
    clears_med = r["median"] > r["linear"]
    cells.setdefault((clears_med, r["linear_beats_null"]), []).append(r["module"])
tt = cells.get((True, True), [])
tf = cells.get((True, False), [])
ft = cells.get((False, True), [])
ff = cells.get((False, False), [])
print("%-38s %22s %16s" % ("", "linear carries signal", "linear is noise"))
print("%-38s %22d %16d" % ("median ceiling exceeds linear", len(tt), len(tf)))
print("%-38s %22d %16d" % ("median ceiling does not", len(ft), len(ff)))
print("\nmodules where the contest is real and the median ceiling still exceeds linear (%d):"
      % len(tt))
for m in tt:
    r = next(x for x in rows if x["module"] == m)
    print("  %-24s median %.4f  linear %.4f  margin %+.4f"
          % (m, r["median"], r["linear"], r["median"] - r["linear"]))

n_real = sum(1 for r in rows if r["linear_beats_null"])
print("\nEffective N: %d modules ran, %d have a linear arm that beats its own null,"
      % (len(rows), n_real))
print("so %d are modules where the comparison carries information." % n_real)
