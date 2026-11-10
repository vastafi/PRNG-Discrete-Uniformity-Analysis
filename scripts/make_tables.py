#!/usr/bin/env python3
"""
make_tables.py

Copies the CSV tables written by reproduce_all.py into tables/, so that no
number or rank in the paper is retyped by hand.

    table_index.csv      composite index, tiers, rank ranges
    table_criteria.csv   observed values and ranks per criterion
    table_scores.csv     score of every source on every criterion

Usage
-----
    python3 scripts/reproduce_all.py --data data --out results/main
    python3 scripts/make_tables.py --results results/main/results.json --out tables
"""

import argparse, os, shutil


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/main/results.json")
    ap.add_argument("--out", default="tables")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    src = os.path.dirname(args.results)
    for name in ("table_index.csv", "table_criteria.csv", "table_scores.csv"):
        p = os.path.join(src, name)
        if os.path.exists(p):
            shutil.copy(p, os.path.join(args.out, name))

    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
