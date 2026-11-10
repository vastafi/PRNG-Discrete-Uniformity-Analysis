# Reproduction package

Code and results for the paper *Concordance with the discrete uniform distribution of ten
commonly used random number sources and the limits of a composite index* (Astafi, Leahu,
Ciorbă).

The digit files are in the companion repository
https://github.com/vastafi/PRNG-Random-Digits-0-9-Dataset (ten sources, six runs of
2,614,690 digits each). Clone it into `data/`, which is not tracked here:

```bash
git clone https://github.com/vastafi/PRNG-Random-Digits-0-9-Dataset data
pip install -r requirements.txt
bash run_all.sh
```

## Layout

```
scripts/            all Python code (run from the repository root)
results/main/       results.json: criteria, null calibration, composite index, per-run indices
results/multiseed/  multi-seed experiment (200 seeds × 6 generators) and its ANOVA
results/extraction/ extraction-path control
results/extra/      additional analyses: positive control, joint null, power curves,
                    Proposition 2, reversal null
figures/            the two figures of the paper (PDF + PNG)
tables/             CSV tables of criterion values, scores and index
run_all.sh          runs the whole pipeline in order
archive_not_used/   older scripts, figures and tables not used in the paper (see its README)
```

## Scripts

| script | what it does |
|---|---|
| `reproduce_all.py` | main pipeline: nine criteria, Monte Carlo calibration, composite index, diagnostics, CSV tables |
| `common.py` | nine criteria, likelihood-ratio statistic G, null calibration, two-sided scores, additive and geometric index (shared by the scripts below) |
| `make_figures.py` | Figure 1 (per-run index) from `results.json` |
| `make_tables.py` | copies the CSV tables into `tables/` |
| `multiseed.py` | the same analysis over many independent seeds, with the variance decomposition |
| `extraction_audit.py` | control experiment: one algorithm, five digit-extraction paths |
| `quick.py` | G for every source, Pearson–G agreement, ANOVA without π, sensitivity of the required number of runs |
| `joint_null.py` | all nine criteria computed together on complete uniform sequences (joint null) |
| `joint_analysis.py` | joint null vs. separately calibrated blocks; tests of block independence |
| `pc1000.py` | positive control over 1,000 replicates, additive and geometric index |
| `final.py` | criterion-correlation table, G-based index variant, Wilson intervals, Monte Carlo check of the power formula |
| `power_curves.py` | detection rate vs. defect size for six rules (additive, geometric, minimum, Fisher, X² alone, serial alone) |
| `reversal_null.py` | null distribution of the disagreement between one-sided and two-sided orderings |
| `ceiling_bound.py` | Proposition 2: the power ceiling of the additive index, P0(I_A < 0.287) |
| `make_figure_power.py` | Figure 2 (power curves) from `pc_tilt.json` and `pc_block.json` |

## Where each result of the paper comes from

| paper | file |
|---|---|
| Figure 1 (per-run index) | `figures/figure1_per_run.pdf` (`make_figures.py`) |
| Figure 2 (power curves) | `figures/figure2_power.pdf` (`make_figure_power.py`) |
| Tables of criterion values, scores and index | `results/main/results.json`, `tables/` |
| Criterion-correlation table | `results/extra/final.json` (`corr_table`) |
| Positive-control table | `results/extra/pc_a.json`, `pc_b.json` |
| Proposition 2 (power ceiling, 0.241) | `results/extra/ceiling_bound.json` |
| One-sided vs two-sided null (0.30, 0.16) | `results/extra/reversal_null.json` |
| Multi-seed ANOVA | `results/multiseed/multiseed.json` |
| Extraction-path control | `results/extraction/extraction_audit.json` |

All seeds are fixed (20260908 and seeds derived from it), so repeated runs give identical
numbers, except the `secrets` row of the multi-seed experiment, which cannot be seeded.
The values in the paper are those of the run published here.

Software versions, hardware and the reproducibility notes are in
[`README_computational_environment.md`](README_computational_environment.md).

## Scoring of the moment criteria

The four moment criteria (mean, variance, skewness, excess kurtosis) are scored two-sided
on the **signed** deviation. Their absolute deviation has its null mode at zero, so a
two-sided score on it would treat a deviation close to zero as atypical. Pearson's
chi-square, the entropy deviation and the Kolmogorov–Smirnov distance are scored
two-sided on the statistic itself, where a value far below the null centre is genuine
under-dispersion; the lag-1 autocorrelation uses the upper tail of |Z|.
