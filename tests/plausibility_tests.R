set.seed(1)
n <- 120; p <- 8
X <- matrix(rnorm(n*p), n, p)
beta <- c(2, -1.5, 0, 0, 1, rep(0, p-5))
y <- X %*% beta + rnorm(n, sd = 1)

forest <- build_paths(X = X, y = as.numeric(y), family = "gaussian", K = min(ncol(X), 10), eps = 1e-6, delta = 1, L = 50)
stab <- stability(X = X, y = as.numeric(y), B = 50, resample_type = "bootstrap", ...)


# compute AIC cutoff manually
models <- forest$aic_by_model
best_aic <- min(models$aic)
Delta <- 2
cutoff <- best_aic + Delta
cutoff

# get plausible models
pm <- plausible_models(forest, stab, Delta = Delta, tau = 0)

# check if all models retain AIC constraint, should be TRUE
all(pm$aic <= cutoff)

# assert correct exclusion
any(models$aic > cutoff & models$model_id %in% pm$model_id)

tau <- 0.6
pm <- plausible_models(forest, stab, Delta = Delta, tau = tau)
# check if models pass threshold, should be TRUE
all(pm$pi_bar >= tau)

# check that no excluded model should have passed, should be FALSE
candidates <- forest$aic_by_model[
  forest$aic_by_model$aic <= best_aic + Delta, ]

any(candidates$pi_bar > tau & !(candidates$model_id %in% pm$model_id))

# check both Delta and tau together
pm <- plausible_models(forest, stab, Delta = Delta, tau = tau)

all(pm$aic <= best_aic + Delta)        # AIC condition
all(pm$pi_bar >= tau)                  # Stability condition



