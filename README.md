# Reproduction and revision package

## Files

| file | what it does |
|---|---|
| `reproduce_all.py` | main pipeline: nine criteria, Monte Carlo calibration, composite index, tiers, diagnostics, CSV tables |
| `make_figures.py` | Figures 1–4 from `results.json`, vector PDF + 600 dpi PNG |
| `extraction_audit.py` | control experiment: one algorithm, five digit-extraction paths |
| `multiseed.py` | the same analysis over many independent seeds, with the variance decomposition |
| `make_tables.py` | tables in Word-ready form, including the criterion-correlation table |
| `requirements.txt` | pinned dependencies |

## Run order

```bash
pip install -r requirements.txt

python3 reproduce_all.py --data data --out results
python3 make_figures.py  --results results/results.json --out figures
python3 extraction_audit.py --n 15688140 --out results_extraction
python3 make_tables.py  --results results/results.json --out tables
python3 multiseed.py --seeds 200 --n 2614690 --out results_multiseed
```

`--data` expects one subfolder per source and one file per run, exactly as in
version 1; folder and file names do not matter.

Runtime on a current laptop, at the sample size of the paper: about 2 minutes
per source for the descriptive part, plus roughly 4 minutes for the
sequence-based null at `--seq-null 200`. `multiseed.py` with 200 seeds and
2.6 million digits per seed takes a few hours; reduce `--n` first, not
`--seeds`.
