"""Reproduce the regulator-activity precondition fragment from the public Zhu
summary table (no big object needed). Downloads the 4.8MB suppl table and reports
cross-state conservation of regulator activity + candidate-module breadth.

Usage: python scripts/precondition_summary.py
"""
import io
import urllib.request

import pandas as pd
from scipy.stats import spearmanr

URL = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
       "marson2025_data/suppl_tables/DE_stats.suppl_table.csv")
CONDS = ["Rest", "Stim8hr", "Stim48hr"]


def main():
    raw = urllib.request.urlopen(URL).read()
    df = pd.read_csv(io.BytesIO(raw))
    piv = (df.pivot_table(index="target_contrast_gene_name", columns="culture_condition",
                          values="n_downstream", aggfunc="first").reindex(columns=CONDS))
    sub = piv.dropna()
    print(f"genes in all 3 conditions: {len(sub)}")
    for a, b in [("Rest", "Stim8hr"), ("Rest", "Stim48hr"), ("Stim8hr", "Stim48hr")]:
        print(f"  breadth Spearman {a} vs {b}: {spearmanr(sub[a], sub[b]).correlation:.3f}")
    for mod, genes in {
        "Th2/GATA3": ["GATA3", "STAT6", "TBX21"],
        "TCR": ["CD3E", "ZAP70", "LAT", "PLCG1"],
    }.items():
        print(f"\n{mod} breadth by condition:")
        print(piv.reindex(genes))


if __name__ == "__main__":
    main()
