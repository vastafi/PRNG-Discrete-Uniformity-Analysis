## Computational environment

The digit files were generated with:

| Source | Software | Machine |
|---|---|---|
| Microsoft Excel | Excel 16.91 (Office LTSC 2021 for Mac) | MacBook Pro (Apple M2) |
| Wolfram Mathematica | Mathematica 12 | Intel Core i9-12900K desktop |
| R | R 4.2.3 | MacBook Pro (Apple M2) |
| Python (NumPy, secrets) | Python 3.9.6 | MacBook Pro (Apple M2) |
| C++ `rand()` | Apple Clang 16.0.0 | MacBook Pro (Apple M2) |
| Java SecureRandom, ThreadLocalRandom | OpenJDK 21.0.3 (Corretto) | MacBook Pro (Apple M2) |

All statistics, simulations, tables and figures are computed in Python with NumPy and SciPy (versions in `requirements.txt`).

## Notes on two sources

- **C++ `rand()`** (from `<cstdlib>`, platform-dependent, not `std::mt19937`): each run is a separate execution seeded with `srand(time(0))` and writing `rand() % 10`. The modulo bias is negligible at this sample size: it raises the expected goodness-of-fit statistic at 15,688,140 digits by about 0.2 even for `RAND_MAX = 32767`. The six files were checked to be distinct and non-overlapping; the seeds were not recorded.
- **Random.org**: integers 0–9 generated from atmospheric noise, collected in batches of at most 10,000 (at least 262 requests per run).

## Reproducibility

Three levels of reproducibility should be distinguished.

1. **Analysis.** Every number reported in the paper is computed from the published digit files by the scripts in this repository, with fixed random seeds (20260908 for the main analysis and seeds derived from it elsewhere). This covers the null calibrations, the joint-null check, the positive control, the power curves and the multi-seed experiment. Repeated runs give identical results, with one exception: the Python `secrets` row of the multi-seed experiment, because `secrets` cannot be seeded.
2. **Generation of the digit files.** This is not reproducible for every source. The C++ seeds were derived from the system clock and not recorded; the Java generators, Excel and Mathematica were seeded by their runtime; the `secrets` module cannot be seeded; Random.org and the digits of π are external data. For these sources the published files in [PRNG-Random-Digits-0-9-Dataset](https://github.com/vastafi/PRNG-Random-Digits-0-9-Dataset) are the primary data. Only the R and NumPy runs can be regenerated from their seeds.
3. **Tables and figures.** Every table and figure can be regenerated from the two repositories alone (`run_all.sh`). The authors regenerated them on macOS, and an independent regeneration on Linux reproduced every number to floating-point precision.
