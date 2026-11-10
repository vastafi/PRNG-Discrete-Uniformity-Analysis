# Not used in the paper

Files kept for the record. They are not needed to reproduce the paper and are not
run by `run_all.sh`. Nothing here is read by the scripts in `scripts/`.

| file | why it is here |
|---|---|
| `scripts/positive_control.py` | first version of the positive control (50 replicates, pair-based null of 200 sequences of 3×10⁶ digits). The paper uses `scripts/pc1000.py` (1,000 replicates). |
| `scripts/make_figures.py` | earlier version that also drew Figures 1, 3 and 4 below; the current one draws only the per-run figure used in the paper. |
| `scripts/make_tables.py` | earlier version that also wrote the criterion-correlation table below. |
| `figures/figure1.*` | composite index by source with the null band (not in the paper). |
| `figures/figure3.*` | score of every source on every criterion (heat map, removed from the paper). |
| `figures/figure4.*` | Spearman correlation of the nine criteria under H0, from an earlier run. The correlation table of the paper is `corr_table` in `results/extra/final.json`. |
| `tables/table_criterion_correlation.*` | the same correlation matrix as `figure4`, from the earlier run. |
| `results_main_csv/` | copies of the CSV tables that `reproduce_all.py` writes next to `results.json`; identical to `tables/`. |
