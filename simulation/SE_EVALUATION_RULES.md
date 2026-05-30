# SE Evaluation Rules

This note defines when it is reasonable to reduce the simulation settings and how to judge whether the reported standard errors are stable enough.

## Short Answer

Reducing `rep` and `k.runs` is statistically acceptable for pilot work and for choosing a computational budget. It is not automatically acceptable for final claims.

The smaller setting is acceptable only if the Monte Carlo error induced by the smaller setting is below a pre-specified tolerance and the substantive conclusions are unchanged.

## What Each Parameter Controls

`rep` controls the Monte Carlo precision of the simulation summaries:

- mean estimate
- bias and relative bias
- empirical SE, meaning `sd(Est)` across simulated datasets
- coverage

`k.runs` controls the Monte Carlo precision of bootstrap SE within each simulated dataset.

These are different sources of uncertainty. A small `k.runs` can make `SD.BT`, `SD.BT1`, and `SD.BT2` unstable even if the main model fit is stable.

## Key Formulas

Approximate relative Monte Carlo error of an empirical SE estimated from `R` simulation repetitions:

```text
rel_MCSE(empirical_SE) ~= 1 / sqrt(2 * (R - 1))
```

Approximate relative Monte Carlo error of a bootstrap SE estimated from `B` bootstrap samples:

```text
rel_MCSE(bootstrap_SE) ~= 1 / sqrt(2 * (B - 1))
```

Examples:

| Count | Relative MCSE |
|---:|---:|
| 2 | 70.7% |
| 5 | 35.4% |
| 10 | 23.6% |
| 20 | 16.2% |
| 25 | 14.4% |
| 50 | 10.1% |
| 100 | 7.1% |

This table applies to `B` for bootstrap SE and approximately to `R` for empirical SE, using the formulas above.

## Recommended Decision Thresholds

For exploratory work:

- `rel_MCSE(empirical_SE) <= 25%` is usable for screening.
- `rel_MCSE(bootstrap_SE) <= 25%` is usable for rough comparison.
- `se_ratio = avg_reported_SE / empirical_SE` should usually be between `0.67` and `1.50`.

For a result we are willing to report:

- `rel_MCSE(empirical_SE) <= 15%`, preferably `<= 10%`.
- `rel_MCSE(bootstrap_SE) <= 15%`, preferably `<= 10%`.
- `se_ratio` should usually be between `0.80` and `1.25`.
- The conclusion should be stable across at least two independent seed batches.
- Bootstrap valid-retention after bias filtering should be at least 80% for the SE variant being reported.

For coverage claims:

- Do not trust coverage from very small `rep`.
- With `rep=20`, a nominal 95% coverage estimate has a worst-case Monte Carlo SE around 11 percentage points.
- With `rep=100`, the worst-case Monte Carlo SE is still around 5 percentage points.
- Coverage needs larger `rep` than rough SE ranking.

## Practical Recommendation For This Project

The existing pilots with `k.runs=1` or `k.runs=2` are useful for timing and checking code paths. They are not enough to settle bootstrap SE.

If we want a reduced but defensible setting, use a staged plan:

1. Run `rep=20`, `k.runs=10` as a screening batch.
2. If `SD.BT` conclusions matter, increase to at least `k.runs=25`; `k.runs=10` has about 23.6% relative Monte Carlo error for bootstrap SE.
3. If empirical SE or coverage conclusions matter, increase `rep` to at least 50, preferably 100.
4. Compare two independent batches with different seeds. If key SE ratios and conclusions differ by more than 10-15%, the reduced setting is too small.

## How To Run The Evaluator

After combining parallel outputs:

```powershell
Rscript .\simulation\combine_parallel_results.R "E:\path\to\parallel_output"
```

Run the SE evaluator:

```powershell
Rscript .\simulation\evaluate_se_stability.R "E:\path\to\parallel_output" --out "E:\path\to\parallel_output\se_evaluation"
```

You can pass multiple combined results:

```powershell
Rscript .\simulation\evaluate_se_stability.R `
  "E:\path\to\Setting_I_batch" `
  "E:\path\to\Setting_II_batch" `
  --out "E:\path\to\se_evaluation_all"
```

The script writes:

- `estimate_metrics.csv`
- `se_metrics.csv`
- `run_diagnostics.csv`
- `SE_EVALUATION_SUMMARY.md`

## Interpretation Of Main Columns

`empirical_se`:
The Monte Carlo empirical standard deviation of estimates across simulated datasets.

`avg_reported_se`:
The root-mean-square reported SE from the fitted method, such as `SD`, `SD.BT`, `SD.BT1`, or `SD.BT2`.

`se_ratio`:
`avg_reported_se / empirical_se`. Values near 1 are better. Large deviations mean the reported SE is not matching empirical variation.

`rel_mcse_empirical_se`:
Monte Carlo uncertainty of `empirical_se` due to finite `rep`.

`bootstrap_se_rel_mcse_nominal`:
Monte Carlo uncertainty of bootstrap SE due to finite `k.runs`.

`coverage`:
Empirical 95% interval coverage. Use only as secondary evidence unless `rep` is sufficiently large.
