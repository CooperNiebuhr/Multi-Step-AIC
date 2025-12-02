#' Step 3: Selecting Plausible Models
#'
#' Combine the full-data multi-path search (Step 1) with stability scores
#' from resampling (Step 2) to obtain a set of final plausible models.
#'
#' Following Algorithm 3 in the project description, we:
#' \enumerate{
#'   \item Find the lowest AIC across all models from the full-data search.
#'   \item Keep all models whose AIC is within \code{Delta} of this minimum.
#'   \item For each retained model \(S\), compute its average stability
#'         \eqn{\pi(S) = (1 / |S|) \sum_{j \in S} \pi_j}.
#'   \item Keep only those models with \eqn{\pi(S) \ge \tau}.
#' }
#'
#' @param forest An object of class \code{"path_forest"} produced by
#'   \code{build_paths()} on the full data.
#' @param stab An object of class \code{"path_stability"} produced by
#'   \code{stability()} on the same predictors.
#' @param Delta Numeric AIC tolerance \eqn{\Delta}. Models with
#'   \code{AIC <= AIC_min + Delta} are considered close to the best model.
#'   Typical value: \code{2}.
#' @param tau Numeric stability threshold \eqn{\tau} in \eqn{[0, 1]}.
#'   Models must have average stability at least \code{tau}. Typical value:
#'   \code{0.6}.
#'
#' @return A data frame of plausible models with one row per model and
#'   at least the columns:
#'   \itemize{
#'     \item \code{model_id}, \code{step}, \code{size}, \code{aic},
#'           \code{key}, \code{vars} (as in \code{forest$aic_by_model}).
#'     \item \code{pi_bar}: average stability \eqn{\pi(S)} for the model.
#'   }
#'   Rows are sorted by \code{aic} (ascending). The returned object has
#'   class \code{"plausible_models"}.
#'
#' @export
plausible_models <- function(forest,
                             stab,
                             Delta = 2,
                             tau   = 0.6) {
  if (!inherits(forest, "path_forest")) {
    stop("`forest` must be an object of class 'path_forest'.")
  }
  if (!inherits(stab, "path_stability")) {
    stop("`stab` must be an object of class 'path_stability'.")
  }
  
  if (!is.numeric(Delta) || length(Delta) != 1L || Delta < 0) {
    stop("`Delta` must be a single non-negative number.")
  }
  if (!is.numeric(tau) || length(tau) != 1L || tau < 0 || tau > 1) {
    stop("`tau` must be a single number in [0, 1].")
  }
  
  models <- forest$aic_by_model
  if (is.null(models) || nrow(models) == 0L) {
    warning("No models found in `forest$aic_by_model`.")
    out <- models
    class(out) <- c("plausible_models", class(out))
    attr(out, "Delta") <- Delta
    attr(out, "tau")   <- tau
    return(out)
  }
  if (is.null(models$vars)) {
    stop("`forest$aic_by_model` must contain a `vars` list-column.")
  }
  
  pi_vec <- stab$pi
  if (is.null(pi_vec)) {
    stop("`stab` must contain a `pi` component with stability scores.")
  }
  
  # 1–3: AIC window around the best model
  best_aic <- min(models$aic, na.rm = TRUE)
  in_window <- models$aic <= best_aic + Delta
  candidates <- models[in_window, , drop = FALSE]
  
  if (nrow(candidates) == 0L) {
    warning("No models within the AIC window defined by `Delta`.")
    out <- candidates
    class(out) <- c("plausible_models", class(out))
    attr(out, "Delta") <- Delta
    attr(out, "tau")   <- tau
    return(out)
  }
  
  # 4: For each model S, compute average stability pi(S)
  n_models <- nrow(candidates)
  pi_bar   <- numeric(n_models)
  
  for (i in seq_len(n_models)) {
    vars_i <- candidates$vars[[i]]
    
    if (length(vars_i) == 0L) {
      # Empty model: average stability is undefined; mark as NA so it's dropped
      pi_bar[i] <- NA_real_
    } else {
      if (any(vars_i < 1L | vars_i > length(pi_vec))) {
        stop("Model ", i, " contains invalid variable indices.")
      }
      pis <- pi_vec[vars_i]
      pi_bar[i] <- mean(pis)
    }
  }
  
  candidates$pi_bar <- pi_bar
  
  # 5: Keep only models with pi(S) >= tau
  keep <- !is.na(candidates$pi_bar) & candidates$pi_bar >= tau
  out <- candidates[keep, , drop = FALSE]
  
  if (nrow(out) > 0L) {
    out <- out[order(out$aic), , drop = FALSE]
  }
  
  attr(out, "Delta") <- Delta
  attr(out, "tau")   <- tau
  class(out) <- c("plausible_models", class(out))
  
  out
}