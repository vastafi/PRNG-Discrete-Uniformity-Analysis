# Reproduction and revision package

Analysis code for *A Multi-Criterion Statistical Assessment of Discrete Uniformity in
Random and Pseudorandom Number Sources* (Astafi, Leahu, Ciorba).

Version 2 of the pipeline. It reproduces every number of version 1 and adds the
corrections required by the statistical objections raised in review.

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

## What changed, and why

**Two-sided calibration.** Version 1 scored each criterion with
`u = P(T_null > T_obs)`, which returns `u → 1` when the observed statistic is
unusually *small*. For chi-square, KS or the entropy deficit an unusually small
value means under-dispersion, not quality, so the published index rewards
agreement that is too good to be true. Version 2 computes
`t = 2·min(P(T_null > T_obs), P(T_null < T_obs))`, which is ~U(0,1) under H0 and
flags both tails. The version-1 score is still computed and reported as
`crqi_v1_oneside`, so the published table can be regenerated.

**Monte Carlo null for the sequence-based criteria.** The serial chi-square on
overlapping pairs is slightly under-dispersed relative to chi2(81) (simulated
variance ≈ 152 against 162, false positive rate 4.25% instead of 5%), so the
asymptotic p-value of version 1 is mildly conservative. Version 1 also built the
null distribution of the index with `U["r1"] = rng.random(B)` and
`U["chiser"] = rng.random(B)`, that is, it assumed both scores are exactly
uniform under H0. Version 2 simulates them.

**Tiers without a magic constant.** Version 1 decided the "top tier" inside
`make_figures.py` with `TOP_CUT = 0.60`, a number that appears nowhere in the
manuscript. Version 2 groups sources by overlap of the 95% confidence interval
of their per-run index, so the tier structure states what the data can resolve.
The rule must be described in the paper.

**The entropy criterion is not independent.** `reproduce_all.py` now verifies
that `H − Ĥ = χ²/(2N ln2)` on the null replicates. At the sample size of the
paper the Spearman correlation between the entropy deviation and the
goodness-of-fit chi-square is exactly 1 and the identity holds to a relative
error below 10⁻³. The entropy criterion therefore carries no information beyond
the chi-square, and the index effectively gives weight 0.30 to a single
statistic. An alternative weighting without entropy is computed alongside
(`index_twoside_no_entropy`), and `make_tables.py` writes the dependence
structure as a Word-ready table.

**Entropy has the wrong target.** The plug-in estimator is biased by
`−(k−1)/(2N ln2) ≈ −4.1×10⁻⁷` at N = 15,688,140, so a perfect generator gives
`|Ĥ − H| ≈ 4.1×10⁻⁷`, not 0. Ranking that column towards zero is close to
meaningless; the calibrated score handles it correctly, the raw table did not.

**Diagnostics the referees will ask for**, all written to `results.json`:
the variance decomposition of the index into within-source and between-source
components with an F test and the ICC; the joint null probability of the
ThreadLocalRandom configuration (excellent low-order moments together with a
poor goodness-of-fit), which is the one anomaly in the paper that is actually
significant; and the Spearman correlation matrix of the nine criteria under H0.

**CSV tables.** `table_criteria.csv`, `table_index.csv` and `table_scores.csv`
are written by the script. The rank discrepancies in Table 1 of the manuscript
came from transcribing ranks by hand; paste these files instead.

## The two experiments that still have to be run

`extraction_audit.py` holds the algorithm and the seed fixed and varies only the
code that turns generator output into decimal digits. On a short test run the
composite index already spans more than 0.3 across the five paths, which is
larger than the between-source range reported in the paper. Until this is run at
the full sample size and the extraction path of each of the ten sources is
documented, differences between sources cannot be attributed to the generators.

`multiseed.py` repeats the whole analysis over many independent seeds. This is
what converts "this stream" into "this generator". If the F test of the
generator effect is not significant, that is the result to report.

## Data

Deposit both the digit streams and this code on Zenodo and cite the DOI. A
GitHub URL is not accepted as persistent archiving by Elsevier, Springer or
Taylor & Francis, and the uncompressed streams exceed GitHub's per-file limit.
