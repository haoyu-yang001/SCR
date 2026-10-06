# SCR 2.0.0

A rewrite of SCR 1.0.0, the package published in the supplementary material of Yang, Qin, Li and Hu (2024, *Biometrics*).

## New features

* `scr_next()` assigns a single newly arrived unit. Calling it repeatedly gives the same assignments as `SCR()` under the same seed.
* `SCR()` supports categorical covariates through `categorical` and `w`, using the combined imbalance measure W(n) of Section 2.2 of the paper.
* `scr_imbalance()` computes the Mahalanobis distance, modified Mahalanobis distance, within-stratum imbalance, W(n) and marginal imbalance.
* `scr_test()` provides a randomization test of the treatment effect and both treatment effect estimators. Its p-value is (1 + #{|t_b| >= |t_obs|}) / (B + 1), so it is never 0.
* `SCR()` now returns an object of class `"SCR"` that also holds the allocation probability and the step of Algorithm 1 used for each unit.

## Changes from 1.0.0

* When the two potential imbalances are equal, the unit is assigned with probability 0.5, as in Algorithm 1. Version 1.0.0 favoured treatment 1.
* `d = 0` no longer leaves units unassigned. All units after the first two then use the marginal step, and the probability is 0.5 when the groups are equal.
* The first two units are assigned to the two arms in random order. Version 1.0.0 always used the same order.
* The marginal constraint is still `|n1 - n0| < d`, as in 1.0.0.
* Input and output assignments are both coded 0/1. Version 1.0.0 took 1/2 as input and returned 0/1.
* `Mahalanobis_Distance` is n (x1bar - x0bar)' cov(x)^-1 (x1bar - x0bar), as defined in the paper. Version 1.0.0 used n/2.
* The algorithm updates its statistics incrementally and is much faster: about 0.05 s instead of 2 s for n = 1000 and p = 10.
* The dependency on tidyverse, dplyr and readr is removed. The package now imports only MASS.
* `smart()` (multi-arm) is kept. It now uses the same tie handling, and the first K units go to the K arms in random order.
