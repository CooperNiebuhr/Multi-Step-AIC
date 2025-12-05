# final_project_group2

<!-- badges: start -->
<!-- badges: end -->

## Overview

`finalprojectgroup2` implements a **multi-path forward model selection** procedure.
Traditional stepwise model selection commits to a single best model at each iteration,
which often leads to instability when many models have very similar AIC values.
This package instead maintains a *frontier* of near-optimal models at each step,
allowing multiple promising paths to be evaluated and compared.

This enables:

- exploration of multiple plausible model sequences,
- identification of unstable variable-selection decisions,
- and post-selection validation through stability-based filtering.

A typical workflow consists of:

1. Build a multi-path model-selection forest using `build_paths()`  
2. Assess selection stability through resampling using `stability()`  
3. Filter and prioritize models by AIC closeness and empirical stability using `plausible_models()`  

For extended examples and diagnostics, see package vignettes:

- `vignette("diabetes_progression")`
- `vignette("branching_behavior")`

## Installation

```r
remotes::install_github("R-4-Data-Science/finalprojectgroup2")
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
forest <- build_paths(
  X      = X,
  y      = y,
  family = "gaussian",
  K      = min(p, 10),
  eps    = 1e-6,
  delta  = 1,
  L      = 50
)
```

### 2. Compute model stability via bootstrap resampling

```r
stab <- stability(
  X             = X,
  y             = y,
  family        = "gaussian",
  B             = 50,
  resample_type = "bootstrap",
  K             = min(p, 10),
  eps           = 1e-6,
  delta         = 1,
  L             = 50
)
```

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

---

## Example usage: Binomial Logistic Regression

```r
set.seed(2)

## Simulate a logistic regression dataset
n <- 200
p <- 8

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)

beta <- c(1.5, -1.0, 0.8, 0, 0, 0, 0, 0)
eta  <- X %*% beta
prob <- 1 / (1 + exp(-eta))

y <- rbinom(n, size = 1, prob = prob)
```

### 1. Multi-path forward selection

```r
forest <- build_paths(
  X      = X,
  y      = y,
  family = "binomial",
  K      = min(p, 10),
  eps    = 1e-6,
  delta  = 1,
  L      = 50
)
```

### 2. Stability assessment

```r
stab <- stability(
  X             = X,
  y             = y,
  family        = "binomial",
  B             = 30,
  resample_type = "bootstrap",
  K             = min(p, 10),
  eps           = 1e-6,
  delta         = 1,
  L             = 50
)
```

### 3. Extract stable and AIC-compatible models

```r
plaus <- plausible_models(
  forest = forest,
  stab   = stab,
  Delta  = 2,
  tau    = 0.6
)

head(plaus)
```

This identifies models that are simultaneously competitive in AIC and empirically stable across resampled datasets.
