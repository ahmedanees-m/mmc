"""Run the Stage-0 precondition test against the Zhu store and print the report.

Prints module-gene coverage, a conservation and rewiring summary per module and
transfer direction, and the per-edge classification (the scaffold for PREREG.md) for
each module and direction that passes. Requires MMC_ZHU_STORE.
"""
from __future__ import annotations

from mmc.data import precondition as pre


def _fmt(t) -> str:
    if t is None:
        return "not measured"
    e, f, r, nd = t
    e = float("nan") if e is None else float(e)
    f = float("nan") if f is None else float(f)
    act = "-" if nd is None else str(int(nd))
    return f"e={e:+.2f} fdr={f:.1e} act={act}"


def main() -> None:
    print("=== module-gene coverage ===")
    for name in pre.MODULES:
        cov = pre.coverage(name)
        print(f"{name}: {len(cov['present'])} present, missing {cov['missing'] or 'none'}")

    print(f"\n=== precondition summary "
          f"(pass = conservation >= 0.5, rewiring > 0, real edges >= {pre.N_MIN}) ===")
    results = pre.run_all()
    print(f"{'module':16s} {'direction':18s} {'cons':>6s} {'rew':>6s} "
          f"{'real':>5s} {'C':>3s} {'R':>3s} {'noeff':>6s} {'untest':>7s}  pass")
    for r in results:
        c = r["counts"]
        cons = f"{r['conservation']:.2f}" if r["conservation"] is not None else "-"
        rew = f"{r['rewiring']:.2f}" if r["rewiring"] is not None else "-"
        direction = f"{r['train']}->{r['test']}"
        print(f"{r['module']:16s} {direction:18s} {cons:>6s} {rew:>6s} "
              f"{r['real_edges']:5d} {c['conserved']:3d} {c['rewired']:3d} "
              f"{c['no_effect']:6d} {c['untestable']:7d}  {'YES' if r['passed'] else '.'}")

    print("\n=== per-edge classification for candidate module and direction pairs "
          "(real edges >= 4) ===")
    for r in sorted(results, key=lambda x: -x["real_edges"]):
        if r["real_edges"] < 4:
            continue
        cons = f"{r['conservation']:.2f}" if r["conservation"] is not None else "-"
        rew = f"{r['rewiring']:.2f}" if r["rewiring"] is not None else "-"
        print(f"\n## {r['module']}  {r['train']} -> {r['test']}  "
              f"(conservation {cons}, rewiring {rew}, "
              f"{'PASS' if r['passed'] else 'no pass'})")
        for e in r["edges"]:
            if e.label in ("conserved", "rewired"):
                print(f"  {e.regulator:6s} -> {e.target:6s}  {e.label:9s}  "
                      f"train[{_fmt(e.train)}]  test[{_fmt(e.test)}]")


if __name__ == "__main__":
    main()
