# Reproduction package

Analysis code for *A Multi-Criterion Statistical Assessment of Discrete Uniformity in
Random and Pseudorandom Number Sources* (Astafi, Leahu, Ciorba).

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

Run each command on its own line. `--data` expects one subfolder per source and
one file per run; folder and file names do not matter.
