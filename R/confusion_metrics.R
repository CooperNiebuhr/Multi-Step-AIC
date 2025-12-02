#' Confusion matrix and diagnostic metrics for logistic models
#'
#' Convenience function to compute a confusion matrix and standard
#' classification metrics at a fixed probability cutoff (default 0.5).
#'
#' @param y_true Numeric or logical vector of true binary outcomes
#'   (0/1, FALSE/TRUE, or equivalent).
#' @param p_hat Numeric vector of predicted probabilities for the
#'   positive class (Pr(Y = 1 | X)).
#' @param cutoff Probability threshold for classifying observations
#'   as positive. Default is 0.5 as in the project outline.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{confusion}: 2x2 matrix with rows = true class (0,1),
#'         columns = predicted class (0,1), in the order
#'         \code{c("0","1")}.
#'   \item \code{metrics}: named numeric vector with
#'         \code{prevalence}, \code{accuracy},
#'         \code{sensitivity}, \code{specificity},
#'         \code{FDR}, and \code{DOR}.
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- rnorm(n)
#' p <- 1 / (1 + exp(-(1 + 2 * x)))
#' y <- rbinom(n, 1, p)
#' p_hat <- p + rnorm(n, sd = 0.05)  # noisy probs
#' p_hat <- pmin(pmax(p_hat, 0), 1)
#' cm <- confusion_metrics(y, p_hat, cutoff = 0.5)
#' cm$confusion
#' cm$metrics
#'
#' @export
confusion_metrics <- function(y_true, p_hat, cutoff = 0.5) {
  if (length(y_true) != length(p_hat)) {
    stop("`y_true` and `p_hat` must have the same length.")
  }
  
  # Coerce outcomes to 0/1 numeric
  y <- as.numeric(y_true)
  if (any(is.na(y))) {
    stop("`y_true` must be coercible to numeric 0/1.")
  }
  
  # Predicted class at cutoff
  if (!is.numeric(p_hat)) {
    stop("`p_hat` must be numeric predicted probabilities.")
  }
  if (any(p_hat < 0 | p_hat > 1, na.rm = TRUE)) {
    warning("Some probabilities are outside [0,1]. Proceeding anyway.")
  }
  
  pred <- ifelse(p_hat >= cutoff, 1L, 0L)
  
  # Confusion counts
  TP <- sum(pred == 1L & y == 1L, na.rm = TRUE)
  TN <- sum(pred == 0L & y == 0L, na.rm = TRUE)
  FP <- sum(pred == 1L & y == 0L, na.rm = TRUE)
  FN <- sum(pred == 0L & y == 1L, na.rm = TRUE)
  
  n <- TP + TN + FP + FN
  if (n == 0L) {
    stop("No valid observations after removing NAs.")
  }
  
  confusion <- matrix(
    c(TN, FP,
      FN, TP),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      truth = c("0", "1"),
      pred  = c("0", "1")
    )
  )
  
  # Metrics
  prevalence  <- (TP + FN) / n
  
  accuracy    <- (TP + TN) / n
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  
  FDR <- if ((TP + FP) > 0) FP / (TP + FP) else NA_real_
  
  # Diagnostic odds ratio: (TP/FN) / (FP/TN) = (TP * TN) / (FP * FN)
  if (TP > 0 && TN > 0 && FP > 0 && FN > 0) {
    DOR <- (TP * TN) / (FP * FN)
  } else {
    DOR <- NA_real_
  }
  
  metrics <- c(
    prevalence  = prevalence,
    accuracy    = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    FDR         = FDR,
    DOR         = DOR
  )
  
  list(
    confusion = confusion,
    metrics   = metrics
  )
}
