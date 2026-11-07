#!/usr/bin/env bash
# Reproduces every number added in version 3 of the paper.
# Run from the repository root, after reproduce_all.py and multiseed.py:
#   bash run.sh
# Total time: about 1.5 h on one core; the commands marked (parallel) can run at once.
set -euo pipefail
OUT=results
mkdir -p "$OUT"

# 1. G statistic, ANOVA without pi, sample-size sensitivity (seconds)
python3 quick.py --results results/results.json --out "$OUT/quick.json"

# 2. Positive control, 1,000 replicates (parallel; ~25 min each)
python3 pc1000.py --reps 500 --start 0   --out "$OUT/pc_a.json"
python3 pc1000.py --reps 500 --start 500 --out "$OUT/pc_b.json"

# 3. Joint null on complete sequences (parallel)
python3 joint_null.py --N 2614690  --B 10000         --out "$OUT/joint_run.npz"   # ~8 min
python3 joint_null.py --N 15688140 --B 2000          --out "$OUT/joint_p0.npz"    # ~13 min
python3 joint_null.py --N 15688140 --B 1500 --tag 1  --out "$OUT/joint_p1.npz"    # ~10 min
python3 joint_null.py --N 15688140 --B 1500 --tag 2  --out "$OUT/joint_p2.npz"    # ~10 min

# 4. Analyses (the order of the pooled files matters)
python3 joint_analysis.py --joint-run "$OUT/joint_run.npz" \
    --joint-pooled "$OUT/joint_p0.npz" "$OUT/joint_p1.npz" "$OUT/joint_p2.npz" \
    --out "$OUT/joint_analysis.json"
python3 final.py --results results/results.json --multiseed results_multiseed/multiseed.json \
    --pc "$OUT/pc_a.json" "$OUT/pc_b.json" --joint-run "$OUT/joint_run.npz" \
    --joint-pooled "$OUT/joint_p0.npz" "$OUT/joint_p1.npz" "$OUT/joint_p2.npz" \
    --out "$OUT/final.json"
echo "done: results in $OUT/"
