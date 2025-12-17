# Multi-Step AIC

<!-- badges: start -->
<!-- badges: end -->

## Overview

`Multi-Step AIC` implements a **multi-path forward model selection** procedure.
Traditional stepwise model selection commits to a single best model at each iteration,
which often leads to instability when many models have very similar AIC values.
This package instead maintains a *frontier* of near-optimal models at each step,
allowing multiple promising paths to be evaluated and compared. Unlike classical 
forward stepwise procedures (e.g. stepAIC), which follow a single greedy path, 
this approach explicitly preserves uncertainty in early selection decisions.

This enables:

- exploration of multiple plausible model sequences,
- identification of unstable variable-selection decisions,
- and post-selection screening through stability-based filtering.

A typical workflow consists of:

1. Build a multi-path model-selection forest using `build_paths()`  
2. Assess selection stability through resampling using `stability()`  
3. Filter and prioritize models by AIC closeness and empirical stability using `plausible_models()`  


## Installation

```r
install.packages("remotes")

remotes::install_github(
  "CooperNiebuhr/Multi-Step-AIC",
  dependencies = TRUE
)
```

---

## Example usage: Gaussian Linear Regression

```r
set.seed(1)
# Simulate a linear regression dataset
n <- 120; p <- 8
X <- matrix(rnorm(n * p), n, p)
beta <- c(2, -1.5, 0, 0, 1, rep(0, p - 5))
y <- as.numeric(X %*% beta + rnorm(n, sd = 1))
colnames(X) <- paste0("x", 1:p)
```

### 1. Build multi-path selection forest
```r
forest <- build_paths(X, y, family = "gaussian")
```
Advanced controls (branching tolerance, frontier size, etc.) are available; see ?build_paths.

### 2. Compute model stability via bootstrap resampling

```r
stab <- stability(
  X             = X,
  y             = y,
  family        = "gaussian",
  B             = 50,
  resample_type = "bootstrap")
```
Advanced branching controls are available; see ?stability.

### 3. Identify plausible models
```r
plaus <- plausible_models(
  forest = forest,
  stab   = stab,
  Delta  = 2,
  tau    = 0.6
)

head(plaus)
#> model_key   size   AIC   avg_stability
#> x1+x5        2     94.1       0.78
#> x1+x3+x5     3     92.6       0.74
```

Models near the AIC optimum that also exhibit high stability are natural candidates for downstream interpretation or validation.


For extended examples and diagnostics, see package vignettes:

- `vignette("diabetes_progression")`
- `vignette("branching_behavior")`


### Relation to Classical Stepwise Selection

Classical forward stepwise procedures select a single variable at each step based on marginal AIC improvement, which can lead to unstable model choices when predictors are correlated or effects are weak.  
Multi-Step AIC generalizes this idea by retaining multiple near-optimal candidates at each step, allowing analysts to assess model uncertainty rather than committing to a single path.

