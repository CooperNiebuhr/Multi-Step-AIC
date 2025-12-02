#' Resampling-based variable stability
#'
#' Runs the multi-path search on many resamples and computes, for each
#' predictor j, the proportion of models it appears in on each resample,
#' then averages those proportions across resamples.
#'
#' @param X Predictor matrix or data frame (n x p).
#' @param y Response vector of length n.
#' @param family Model family, "gaussian" (lm) or "binomial" (glm).
#' @param B Number of resamples.
#' @param resample_type Type of resampling: "bootstrap" or "subsample".
#' @param m Subsample size if \code{resample_type = "subsample"}.
#'   If \code{NULL}, defaults to \code{floor(sqrt(n))}.
#' @param K,eps,delta,L Multi-path parameters passed to \code{build_paths()}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{"path_stability"} with components:
#' \itemize{
#'   \item \code{pi}: numeric vector of length \code{p}, stability scores
#'         between 0 and 1 for each predictor.
#'   \item \code{z}: B x p matrix of per-resample proportions
#'         \eqn{z_j^{(b)}}.
#'   \item \code{meta}: list with resampling and model parameters.
#' }
#'
#' @export
stability <- function(X, y,
                      family = c("gaussian", "binomial"),
                      B = 50,
                      resample_type = c("bootstrap", "subsample"),
                      m = NULL,
                      K = NULL,
                      eps = 1e-6,
                      delta = 0,
                      L = 50,
                      seed = NULL) {
  family <- match.arg(family)
  resample_type <- match.arg(resample_type)
  
  X <- as.data.frame(X)
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(y) != n) {
    stop("Length of y must match nrow(X).")
  }
  
  if (is.null(m) && resample_type == "subsample") {
    m <- floor(sqrt(n))
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Matrix to hold z_j^(b) values
  z_mat <- matrix(0, nrow = B, ncol = p)
  
  for (b in seq_len(B)) {
    # Draw indices for this resample
    idx <- switch(
      resample_type,
      "bootstrap" = sample.int(n, size = n, replace = TRUE),
      "subsample" = sample.int(n, size = m, replace = FALSE)
    )
    
    X_b <- X[idx, , drop = FALSE]
    y_b <- y[idx]
    
    # Run multi-path search on the resample
    forest_b <- build_paths(
      X = X_b,
      y = y_b,
      family = family,
      K = K,
      eps = eps,
      delta = delta,
      L = L
    )
    
    # Extract all models (across all levels) on this resample
    models_b <- forest_b$aic_by_model$vars
    n_models_b <- length(models_b)
    
    # Defensive: in principle there is always at least the empty model
    if (n_models_b == 0L) {
      next
    }
    
    counts <- integer(p)
    for (k in seq_len(n_models_b)) {
      vars_k <- models_b[[k]]
      if (length(vars_k) > 0L) {
        counts[vars_k] <- counts[vars_k] + 1L
      }
    }
    
    z_mat[b, ] <- counts / n_models_b
  }
  
  # Average over resamples to get π_j
  pi_vec <- colMeans(z_mat)
  
  meta <- list(
    family        = family,
    B             = B,
    resample_type = resample_type,
    m             = m,
    K             = K,
    eps           = eps,
    delta         = delta,
    L             = L,
    n             = n,
    p             = p,
    call          = match.call()
  )
  
  out <- list(
    pi   = pi_vec,
    z    = z_mat,
    meta = meta
  )
  class(out) <- "path_stability"
  out
}
