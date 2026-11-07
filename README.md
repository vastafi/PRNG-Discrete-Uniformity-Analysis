# Reproduction package

The digit files are in the companion repository
https://github.com/vastafi/PRNG-Random-Digits-0-9-Dataset (ten sources, six runs of
2,614,690 digits each).

## Files

| file | what it does |
|---|---|
| `reproduce_all.py` | main pipeline: nine criteria, Monte Carlo calibration, composite index, diagnostics, CSV tables |
| `make_figures.py` | Figures 1–4 from `results.json`, vector PDF + 600 dpi PNG |
| `make_tables.py` | tables in Word-ready form, including the criterion-correlation table |
| `positive_control.py` | false-alarm rate and power of the index on deliberately corrupted streams |
| `multiseed.py` | the same analysis over many independent seeds, with the variance decomposition |
| `extraction_audit.py` | control experiment: one algorithm, five digit-extraction paths |
| `requirements.txt` | dependencies | — |


## Run order

```bash
pip install -r requirements.txt

python3 reproduce_all.py --data data --out results
python3 make_figures.py  --results results/results.json --out figures
python3 make_tables.py   --results results/results.json --out tables
python3 positive_control.py --data data --reps 50 --out results_control
python3 multiseed.py --seeds 200 --n 2614690 --null 120000 --seq-null 4000 --out results_multiseed
python3 extraction_audit.py --n 2614690 --null 120000 --seq-null 4000 --out results_extraction
```
## Additional analyses

These scripts sit in the repository root next to `reproduce_all.py` and read only its
outputs and those of `multiseed.py`.

| file | what it does |
|---|---|
| `common.py` | nine criteria, likelihood-ratio statistic G, null calibration, two-sided scores, additive and geometric index (shared by the scripts below) |
| `quick.py` | G for every source, Pearson–G agreement, ANOVA without pi, sensitivity of the required number of runs |
| `joint_null.py` | all nine criteria computed together on complete uniform sequences (joint null) |
| `joint_analysis.py` | joint null vs. separately calibrated blocks; tests of block independence |
| `pc1000.py` | positive control over 1,000 replicates, additive and geometric index |
| `final.py` | criterion correlation table, G-based index variant, Wilson intervals, Monte Carlo check of the power formula |
| `power_curves.py` | detection rate vs. defect size for six rules (additive, geometric, minimum, Fisher, X² alone, serial alone) |
| `reversal_null.py` | null distribution of the disagreement between one-sided and two-sided orderings, and where the observed one falls |
| `run_all.sh` | runs the whole pipeline in order |

Outputs are written to `results_extra/`.

### Scoring of the moment criteria

The four moment criteria (mean, variance, skewness, excess kurtosis) are scored two-sided
on the **signed** deviation. Their absolute deviation has its null mode at zero, so a
two-sided score on it would treat a deviation close to zero as atypical. Pearson's
chi-square, the entropy deviation and the Kolmogorov–Smirnov distance are scored
two-sided on the statistic itself, where a value far below the null centre is genuine
under-dispersion; the lag-1 autocorrelation uses the upper tail of |Z|.

### Run everything

```bash
bash run_all.sh
```

All seeds are fixed (20260908 and seeds derived from it), so repeated runs give identical
numbers, except the secrets row of the multi-seed experiment, which cannot be seeded.

### Note on `results_control/positive_control.json`

That file comes from an earlier 50-replicate run with a pair-based null of 200 sequences
of 3×10⁶ digits and does not match the paper. The positive control reported in the paper
is `results_extra/pc_a.json` + `pc_b.json` (1,000 replicates, pair-based null of 4,000
sequences of 10⁶ digits).
