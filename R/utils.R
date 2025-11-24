#' Internal utility: compute AIC for a given variable set
#'
#' @keywords internal
fit_aic <- function(vars_idx, df, var_names, family) {
  # vars_idx: integer indices into var_names (i.e., columns of X)
  # df: data.frame with response column y and all predictors
  # var_names: character vector of predictor names (colnames(X))
  # family: "gaussian" or "binomial"
  
  if (length(vars_idx) == 0L) {
    form <- y ~ 1
  } else {
    vars_char <- var_names[vars_idx]
    form <- as.formula(
      paste("y ~", paste(vars_char, collapse = " + "))
    )
  }
  
  if (family == "gaussian") {
    fit <- stats::lm(formula = form, data = df)
  } else if (family == "binomial") {
    fit <- stats::glm(formula = form,
                      family = stats::binomial(),
                      data = df)
  } else {
    stop("Unsupported family: ", family)
  }
  
  aic_val <- stats::AIC(fit)
  list(aic = aic_val, vars = vars_idx)
}

#' Internal utility: build a key string for a model
#'
#' @keywords internal
model_key <- function(vars_idx) {
  if (length(vars_idx) == 0L) return("")
  paste(sort(vars_idx), collapse = ",")
}
