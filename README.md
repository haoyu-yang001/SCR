# SCR: Sequential Covariate-Adjusted Randomization

<!-- badges: start -->
[![R-CMD-check](https://github.com/haoyu-yang001/SCR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haoyu-yang001/SCR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

An R implementation of the sequential covariate-adjusted randomization (SCR) procedure from

> Yang, H., Qin, Y., Li, Y., and Hu, F. (2024). Sequential covariate-adjusted randomization via hierarchically minimizing Mahalanobis distance and marginal imbalance. *Biometrics*, 80(2), ujae047. [doi:10.1093/biomtc/ujae047](https://doi.org/10.1093/biomtc/ujae047)

SCR randomizes participants **individually and sequentially**, as they arrive. It separates the
marginal imbalance (the difference in group sizes) from the covariate imbalance (a modified
Mahalanobis distance that uses covariate sums instead of means) and minimizes them in an explicit
hierarchical order. Covariates are balanced only while the group sizes differ by less than `d`.

## Installation

SCR requires R (>= 4.0.0). Its only dependency outside base R is **MASS**, which ships with R.

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("haoyu-yang001/SCR")
```

To also build the vignette (needs **knitr** and **rmarkdown**):

```r
remotes::install_github("haoyu-yang001/SCR", build_vignettes = TRUE)
vignette("SCR")
```

Or install from a local clone:

```sh
git clone https://github.com/haoyu-yang001/SCR.git
R CMD build SCR
R CMD INSTALL SCR_*.tar.gz
```

## Usage

```r
library(SCR)

set.seed(1)
n <- 500; p <- 10
x <- matrix(rnorm(n * p), n, p)

# Allocate a whole sequence of participants (Algorithm 1)
fit <- SCR(x, d = 5, q1 = 0.75, q2 = 0.85)
fit$assignment            # 1 = treatment, 0 = control
fit$Mahalanobis_Distance

# Randomize a newly arrived participant given those already enrolled
scr_next(x[101, ], covariate = x[1:100, ], assignment = fit$assignment[1:100])

# Continuous + categorical covariates: W(n) = w * M(n) + (1 - w) * ||D_n||^2
z <- data.frame(sex = sample(c("F", "M"), n, TRUE), site = sample(1:4, n, TRUE))
SCR(x, categorical = z, w = 0.7)

# Imbalance measures, randomization test, multi-arm extension
scr_imbalance(x, fit$assignment)
y <- 0.3 * fit$assignment + rowSums(x) + rnorm(n)
scr_test(y, x, fit$assignment, B = 500)
smart(x, K = 3, d = 5, q = 0.75)
```

See `vignette("SCR")` for a full walk-through (install with `build_vignettes = TRUE`).

| Function          | Purpose                                                              |
|-------------------|----------------------------------------------------------------------|
| `SCR()`           | Allocate a sequence of units (two arms)                              |
| `scr_next()`      | Allocate one newly arrived unit given the units already enrolled     |
| `scr_imbalance()` | Mahalanobis, modified Mahalanobis, within-stratum and marginal imbalance |
| `scr_test()`      | Randomization test and treatment effect estimates                    |
| `smart()`         | Multi-arm version                                                    |

## Parameters

* `d`: marginal imbalance constraint. Covariates are balanced while `|n1 - n0| < d`. The paper recommends a moderate value such as 5 or 10. `d = n` balances covariates only.
* `q1`: biased-coin probability of the covariate-adaptive step (default 0.75).
* `q2`: biased-coin probability of the marginal step (default 0.85).
* `w`: weight of the Mahalanobis part when categorical covariates are included (default 0.7).

## Citation

If you use SCR in your work, please cite the paper:

> Yang, H., Qin, Y., Li, Y., and Hu, F. (2024). Sequential covariate-adjusted randomization via hierarchically minimizing Mahalanobis distance and marginal imbalance. *Biometrics*, 80(2), ujae047. https://doi.org/10.1093/biomtc/ujae047

BibTeX:

```bibtex
@article{yang2024scr,
  title   = {Sequential Covariate-Adjusted Randomization via Hierarchically Minimizing {M}ahalanobis Distance and Marginal Imbalance},
  author  = {Yang, Haoyu and Qin, Yichen and Li, Yang and Hu, Feifang},
  journal = {Biometrics},
  year    = {2024},
  volume  = {80},
  number  = {2},
  pages   = {ujae047},
  doi     = {10.1093/biomtc/ujae047}
}
```

The same reference is available in R:

```r
citation("SCR")
```

## License

GPL (>= 2)
