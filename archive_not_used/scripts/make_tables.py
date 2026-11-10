#!/usr/bin/env python3
"""
make_tables.py

Writes the tables of the paper in a form that can be pasted straight into Word,
so that no number and no rank is ever retyped by hand. The rank discrepancies in
Table 1 of the first version came from manual transcription.

Output, in --out:

    table_criterion_correlation.csv   Spearman correlation of the nine criteria
    table_criterion_correlation.txt   the same, tab separated, for Word
    table_index.csv                   composite index, tiers, rank ranges
    table_criteria.csv                observed values and ranks per criterion

To turn a .txt file into a Word table: paste it, select it, then
Insert > Table > Convert Text to Table, separator = tab.

Usage
-----
    python3 reproduce_all.py --data data --out results
    python3 make_tables.py --results results/results.json --out tables
"""

import argparse, csv, json, os, shutil

CRIT_LABEL = {"dmean": "mean", "dvar": "variance", "dskew": "skewness",
              "dexc": "excess kurtosis", "dH": "entropy", "chi": "chi2 GoF",
              "dks": "KS distance", "chiser": "chi2 serial", "r1": "lag-1 autocorr."}


def correlation_table(res, out):
    cc = res["diagnostics"]["criterion_correlations"]
    crits, M = cc["criteria"], cc["spearman"]
    lab = [CRIT_LABEL.get(c, c) for c in crits]

    with open(os.path.join(out, "table_criterion_correlation.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([""] + lab)
        for i, name in enumerate(lab):
            w.writerow([name] + [f"{M[i][j]:.3f}" for j in range(len(lab))])

    # lower triangle only: the upper half is redundant and wastes column width
    lines = ["\t".join([""] + lab)]
    for i, name in enumerate(lab):
        row = [name] + [f"{M[i][j]:.3f}" if j <= i else "" for j in range(len(lab))]
        lines.append("\t".join(row))
    txt = "\n".join(lines)
    with open(os.path.join(out, "table_criterion_correlation.txt"), "w") as fh:
        fh.write(txt + "\n")

    print(txt)
    hi = [(lab[i], lab[j], M[i][j])
          for i in range(len(lab)) for j in range(i)
          if abs(M[i][j]) > 0.9]
    if hi:
        print("\ncriterion pairs with |rho| > 0.9 under H0:")
        for a, b, r in hi:
            print(f"  {a} and {b}: rho = {r:.4f}")
        print("these do not carry independent information and must not be counted"
              "\nas separate evidence in the text.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/results.json")
    ap.add_argument("--out", default="tables")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    with open(args.results) as fh:
        res = json.load(fh)
    correlation_table(res, args.out)

    src = os.path.dirname(args.results)
    for name in ("table_index.csv", "table_criteria.csv", "table_scores.csv"):
        p = os.path.join(src, name)
        if os.path.exists(p):
            shutil.copy(p, os.path.join(args.out, name))

    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
