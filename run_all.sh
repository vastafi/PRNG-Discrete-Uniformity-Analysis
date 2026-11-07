#!/usr/bin/env bash
# Reproduces every number, table and figure of the paper from the digit files.
# Run from the repository root:   bash run_all.sh
# Total time: several hours on one core. The commands marked (parallel) are
# independent and can be started at the same time in separate terminals.
set -euo pipefail

# 1. Main analysis: nine criteria, null calibration, composite index, per-run indices
python3 reproduce_all.py --data data --out results
python3 make_figures.py  --results results/results.json --out figures
python3 make_tables.py   --results results/results.json --out tables

# 2. Multi-seed experiment and extraction-path control
python3 multiseed.py --seeds 200 --n 2614690 --null 120000 --seq-null 4000 --out results_multiseed
python3 extraction_audit.py --n 2614690 --null 120000 --seq-null 4000 --out results_extraction

# 3. Additional analyses
OUT=results_extra
mkdir -p "$OUT"
# G statistic, ANOVA without pi, sample-size sensitivity (seconds)
python3 quick.py --results results/results.json --out "$OUT/quick.json"
# positive control, 1,000 replicates (parallel; ~25 min each)
python3 pc1000.py --reps 500 --start 0   --out "$OUT/pc_a.json"
python3 pc1000.py --reps 500 --start 500 --out "$OUT/pc_b.json"
# joint null on complete sequences (parallel)
python3 joint_null.py --N 2614690  --B 10000         --out "$OUT/joint_run.npz"   # ~8 min
python3 joint_null.py --N 15688140 --B 2000          --out "$OUT/joint_p0.npz"    # ~13 min
python3 joint_null.py --N 15688140 --B 1500 --tag 1  --out "$OUT/joint_p1.npz"    # ~10 min
python3 joint_null.py --N 15688140 --B 1500 --tag 2  --out "$OUT/joint_p2.npz"    # ~10 min
# analyses (the order of the pooled files matters)
python3 joint_analysis.py --joint-run "$OUT/joint_run.npz" \
    --joint-pooled "$OUT/joint_p0.npz" "$OUT/joint_p1.npz" "$OUT/joint_p2.npz" \
    --out "$OUT/joint_analysis.json"
python3 final.py --results results/results.json --multiseed results_multiseed/multiseed.json \
    --pc "$OUT/pc_a.json" "$OUT/pc_b.json" --joint-run "$OUT/joint_run.npz" \
    --joint-pooled "$OUT/joint_p0.npz" "$OUT/joint_p1.npz" "$OUT/joint_p2.npz" \
    --out "$OUT/final.json"

echo "done"
