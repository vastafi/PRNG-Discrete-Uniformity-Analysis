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
| `requirements.txt` | dependencies |
| `common.py` | nine criteria, likelihood-ratio statistic G, null calibration, two-sided scores, additive and geometric index | all |
| `quick.py` | G for every source, Pearson–G agreement, ANOVA without pi, sensitivity of the required number of runs | Sections 3.2, 5; Tables 2, 5 |
| `joint_null.py` | all nine criteria computed together on complete uniform sequences (joint null) | Section 3.1 |
| `joint_analysis.py` | joint null vs block-product null; tests of block independence | Section 3.1 |
| `pc1000.py` | positive control over 1,000 replicates, additive and geometric index | Section 6.3, Table 6 |
| `final.py` | correlation table, G-based index variant, Wilson intervals, Monte Carlo check of the power formula | Tables 1, 6; Section 5 |
| `run.sh` | runs all of the above in order | — |


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


### Run

```bash
python3 reproduce_all.py --data data --out results          # version-2 pipeline, if not yet run
python3 multiseed.py --seeds 200 --n 2614690 --null 120000 --seq-null 4000 --out results_multiseed
bash run.sh                                            # about 1.5 h on one core
```

Run each command on its own line. `--data` expects one subfolder per source and
one file per run; folder and file names do not matter.
