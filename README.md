# Fitting the Background Intensity in ETAS Models Using Regression-Based Poisson Point Process Models

This repository contains the R code used to produce the simulation study, the
empirical application, and the additional analyses carried out in response
to reviewer comments, for the manuscript:

> Tarantino, M., D'Angelo, N., Chiodi, M., Adelfio, G. "Fitting the
> Background Intensity in ETAS Models Using Regression-Based Poisson Point
> Process Models." Submitted to *Stochastic Environmental Research and Risk
> Assessment*.

## Overview

The proposed method (ETAS-P) replaces the kernel-based background component
of the classical ETAS model with a regression-based Poisson point process
background, fitted through an alternating estimation procedure (ETAS
parameter update → probabilistic declustering → Poisson background
regression). Two competing estimators are used as benchmarks throughout:
the classical kernel-based ETAS estimator, and its FLP-enhanced version
(ETAS-FLP).

The repository is organised into three parts, reflecting the structure of
the manuscript:

1. **Core estimation code** (`R/01_core_estimation/`) — defines the three
   estimators, the ETAS catalogue simulator, and the main simulation plan.
2. **Application code** (`R/02_application/`) — fits the three models to
   the Italian earthquake catalogue.
3. **Reviewer-response analyses** (`R/03_reviewer_response_analyses/`) —
   seven additional scripts, each written to address a specific reviewer
   comment (see the table below).

## Requirements

- R (>= 4.2)
- Packages: `mgcv`, `etasFLP`, `parallel`, `ggplot2`

```r
install.packages(c("mgcv", "ggplot2"))
install.packages("etasFLP")
```

## Data

The empirical application uses the Italian earthquake catalogue
`catalog.withcov`, distributed with the `etasFLP` R package:

```r
library(etasFLP)
data(catalog.withcov)
```

No external data files are required or provided in this repository; the
catalogue is loaded automatically once the package is attached.

## File-by-file guide

### `R/01_core_estimation/`

| File | Produces | Seeds / tolerances / quadrature |
|---|---|---|
| `etas_parametric_final_plan_no_smooth.R` | Defines `etasclass()` (classical ETAS and ETAS-FLP, via the `flp` argument) and `etasclass.par()` (ETAS-P); runs the main simulation plan and the four original secondary sensitivity analyses (catalogue size, background misspecification, local background irregularity, extreme triggering). Produces Tables S1–S8 of the Supplementary Information and Figures 1–6 of the main text. | Per-replicate simulation seeds are set explicitly in the plan orchestration section (see `SEED_BASE` and the seed construction inside the replicate loop). Optimisation tolerances: `epsmax`, `iterlim`. Berman–Turner quadrature construction for the background regression: `mult.bg`, `ncube.bg`; number of angular sectors for the triggering kernel: `ntheta`. |
| `add_etas_flp_to_existing_results_NO_OVERWRITE_v3.R` | Adds the ETAS-FLP fit (`flp=TRUE`) to catalogues already simulated by the main plan, without overwriting existing ETAS/ETAS-P results. | Reads the `.rds` outputs of the main plan; does not simulate new catalogues, so no additional seeds are introduced. |

### `R/02_application/`

| File | Produces | Notes |
|---|---|---|
| `studio_simulazione_etas_parametrico.R` | Fits ETAS, ETAS-FLP, and ETAS-P to the Italian catalogue (`catalog.withcov`). Produces Table 2 (ETAS parameters) and Table 3 (background component summary) of the main text, and Figures 7–8. | Uses `magn.threshold = 2.5` and `magn.threshold.back = 3.9`; the latter is verified in Supplementary Information Sec. S3 to have no effect on the reported fits, given how `w`/`hdef` are supplied in these calls. |

### `R/03_reviewer_response_analyses/`

| File | Reviewer comment addressed | Produces |
|---|---|---|
| `etas_nonstationary_background_robustness_check.R` | Reviewer 1, points 1–2 (scope of the proposed model relative to genuinely nonstationary ETAS formulations; robustness to an unmodelled linear trend in the background rate) | Supplementary Information, Sec. S2.5, Tables S9–S11 |
| `etas_coverage_bootstrap.R` | Reviewer 1, point 3 (coverage of conditional vs. parametric-bootstrap confidence intervals for the background regression coefficients, validated on simulated data with known truth) | Supplementary Information, Sec. S2.6, Table S12 |
| `catalog_completeness_analysis.R` | Reviewer 1, point 5 / Reviewer 3 (completeness magnitude $M_c$ vs. the threshold $M_0$ used for estimation) | Supplementary Information, Sec. S3.1, Table S13 |
| `m0_sensitivity_application.R` | Reviewer 1, point 5 / Reviewer 3 (sensitivity of the empirical application to the completeness threshold $M_0$) | Supplementary Information, Sec. S3.1, Tables S14–S15 |
| `etas_application_bootstrap.R` | Reviewer 1, point 3 (parametric bootstrap of the complete alternating procedure, applied directly to the empirical fit reported in the main text) | Supplementary Information, Sec. S3.2, Table S16 |
| `explore_nonstationary_results.R` | — (exploration/reporting script) | Console tables and diagnostic plots summarising the output of `etas_nonstationary_background_robustness_check.R` |
| `explore_coverage_bootstrap_results.R` | — (exploration/reporting script) | Console tables and diagnostic plots summarising the output of `etas_coverage_bootstrap.R` |

## Reproducibility notes

- All stochastic simulations use explicit, fixed seeds (see the table above
  and in-code comments in each script); re-running a script with the same
  seed reproduces the corresponding table or figure, subject to R and
  package version differences.
- Optimisation bounds, starting values, and convergence tolerances
  (`epsmax`, `iterlim`, `params.ind`) are set explicitly at each call to
  `etasclass()`/`etasclass.par()`, following the conventions of the
  `etasFLP` package (Chiodi and Adelfio, 2017).
- Scripts in `R/03_reviewer_response_analyses/` checkpoint their progress
  incrementally (partial results are saved to disk after each replicate or
  threshold) and are resumable: re-running a script after an interruption
  skips work already completed rather than starting over.
- Several of these analyses are computationally intensive (parametric
  bootstraps with 30–150 replicates, each requiring a full refit of the
  alternating ETAS-P procedure) and were run on a multi-core machine.
  Parallelisation is configurable via environment variables documented at
  the top of each script (typically `*_N_CORES` and `*_MODE`); setting
  `*_N_CORES = 1` reproduces the fully sequential behaviour.

## Citation

If you use this code, please cite the manuscript above (full reference to
be updated upon acceptance) and, where applicable, this repository via its
Zenodo DOI: **[to be added once the Zenodo release is created]**.

## License

This repository is released under the Apache License 2.0. See `LICENSE`
for details.
