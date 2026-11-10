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

## Reproducibility

Three levels of reproducibility should be distinguished.

1. **Analysis.** Every number reported in the paper is computed from the published digit files by the scripts in this repository, with fixed random seeds (20260908 for the main analysis and seeds derived from it elsewhere). This covers the null calibrations, the joint-null check, the positive control, the power curves and the multi-seed experiment. Repeated runs give identical results, with one exception: the Python `secrets` row of the multi-seed experiment, because `secrets` cannot be seeded.
2. **Generation of the digit files.** This is not reproducible for every source. The C++ seeds were derived from the system clock and not recorded; the Java generators, Excel and Mathematica were seeded by their runtime; the `secrets` module cannot be seeded; Random.org and the digits of π are external data. For these sources the published files in [PRNG-Random-Digits-0-9-Dataset](https://github.com/vastafi/PRNG-Random-Digits-0-9-Dataset) are the primary data. Only the R and NumPy runs can be regenerated from their seeds.
3. **Tables and figures.** Every table and figure can be regenerated from the two repositories alone (`run_all.sh`). The authors regenerated them on macOS, and an independent regeneration on Linux reproduced every number to floating-point precision.
