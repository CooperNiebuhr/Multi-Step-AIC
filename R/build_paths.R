#' Multi-path forward selection using AIC
#'
#' @param X Predictor matrix or data frame (n x p).
#' @param y Response vector of length n.
#' @param family Model family, "gaussian" (lm) or "binomial" (glm).
#' @param K Maximum number of steps (maximum model size). If NULL, uses min(p, 10).
#' @param eps Minimum AIC improvement required to expand a parent model.
#' @param delta AIC tolerance for keeping near-ties among children of the same parent.
#' @param L Optional cap on the number of models kept per level (after dedup).
#'
#' @return An object of class \code{"path_forest"}:
#' \itemize{
#'   \item \code{frontiers}: list of data frames, one per level (step 0 = empty model).
#'   \item \code{aic_by_model}: data frame summarizing all models across all steps.
#'   \item \code{meta}: list of metadata (family, K, eps, delta, L, variable names, etc.).
#' }
#'
#' Each frontier data frame has columns:
#' \itemize{
#'   \item \code{model_id}: integer id (unique across all steps).
#'   \item \code{step}: step index (0 for empty model).
#'   \item \code{size}: number of predictors in the model.
#'   \item \code{aic}: AIC value.
#'   \item \code{key}: character key encoding the variable set (sorted indices).
#'   \item \code{vars}: list-column of integer indices of included predictors.
#' }
#'
#' @export
build_paths <- function(X, y,
                        family = c("gaussian", "binomial"),
                        K = NULL,
                        eps = 1e-6,
                        delta = 0,
                        L = 50) {
  family <- match.arg(family)
  
  # Coerce X to data frame and give default names if needed
  X <- as.data.frame(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }
  
  p <- ncol(X)
  n <- nrow(X)
  
  if (length(y) != n) {
    stop("Length of y must match nrow(X).")
  }
  
  if (is.null(K)) {
    K <- min(p, 10L)
  }
  
  # Data for modeling and predictor names
  df <- data.frame(y = y, X, check.names = FALSE)
  var_names <- colnames(X)
  
  # Storage
  frontiers <- list()
  all_rows <- list()  # for aic_by_model
  model_id <- 0L
  
  # Step 0: empty model
  empty_info <- fit_aic(integer(0L), df = df, var_names = var_names, family = family)
  model_id <- model_id + 1L
  
  frontier0 <- data.frame(
    model_id = model_id,
    step     = 0L,
    size     = 0L,
    aic      = empty_info$aic,
    key      = model_key(integer(0L)),
    stringsAsFactors = FALSE
  )
  frontier0$vars <- list(integer(0L))
  
  frontiers[[1L]] <- frontier0
  all_rows[[1L]] <- frontier0
  
  current_frontier <- frontier0
  
  # Main forward-selection loop
  for (k in seq_len(K)) {
    parent_frontier <- current_frontier
    
    # For each parent, we will generate candidate children
    n_parents <- nrow(parent_frontier)
    
    # Store all children temporarily as a list
    children_list <- list()
    child_counter <- 0L
    
    # Track best child AIC per parent
    parent_best_aic <- rep(Inf, n_parents)
    parent_children_indices <- vector("list", n_parents)
    
    for (i in seq_len(n_parents)) {
      parent_vars <- parent_frontier$vars[[i]]
      parent_aic  <- parent_frontier$aic[i]
      
      available <- setdiff(seq_len(p), parent_vars)
      if (length(available) == 0L) next
      
      for (j in available) {
        child_vars <- sort(c(parent_vars, j))
        info <- fit_aic(child_vars, df = df, var_names = var_names, family = family)
        
        child_counter <- child_counter + 1L
        children_list[[child_counter]] <- list(
          parent_row = i,
          vars       = info$vars,
          aic        = info$aic,
          size       = length(info$vars),
          key        = model_key(info$vars)
        )
        
        parent_children_indices[[i]] <-
          c(parent_children_indices[[i]], child_counter)
        
        if (info$aic < parent_best_aic[i]) {
          parent_best_aic[i] <- info$aic
        }
      }
    }
    
    # No children at all (e.g., no remaining variables)
    if (child_counter == 0L) {
      break
    }
    
    # Filter children per parent using eps and delta
    keep_child <- rep(FALSE, child_counter)
    
    for (i in seq_len(n_parents)) {
      idxs <- parent_children_indices[[i]]
      if (length(idxs) == 0L) next
      
      parent_aic  <- parent_frontier$aic[i]
      best_child  <- parent_best_aic[i]
      
      # Check if best child improves by at least eps
      if (is.finite(best_child) && best_child <= parent_aic - eps) {
        # Keep children within delta of the parent's best child
        for (idx in idxs) {
          if (children_list[[idx]]$aic <= best_child + delta) {
            keep_child[idx] <- TRUE
          }
        }
      }
      # else: this parent contributes no children
    }
    
    if (!any(keep_child)) {
      # No parent achieved sufficient improvement
      break
    }
    
    kept_children <- children_list[keep_child]
    if (length(kept_children) == 0L) {
      break
    }
    
    # Deduplicate by key, keeping the lowest AIC for each key
    keys <- vapply(kept_children, function(z) z$key, character(1))
    unique_keys <- unique(keys)
    
    dedup_list <- vector("list", length(unique_keys))
    names(dedup_list) <- unique_keys
    
    for (idx in seq_along(kept_children)) {
      ch <- kept_children[[idx]]
      kkey <- ch$key
      
      existing <- dedup_list[[kkey]]
      if (is.null(existing) || ch$aic < existing$aic) {
        dedup_list[[kkey]] <- ch
      }
    }
    
    dedup_children <- Filter(Negate(is.null), dedup_list)
    
    # If too many models, keep best L by AIC globally
    if (!is.null(L) && length(dedup_children) > L) {
      aics_vec <- vapply(dedup_children, function(z) z$aic, numeric(1))
      ord <- order(aics_vec)
      dedup_children <- dedup_children[ord[seq_len(L)]]
    }
    
    n_new <- length(dedup_children)
    if (n_new == 0L) {
      break
    }
    
    # Build frontier data frame for step k
    frontier_k <- data.frame(
      model_id = model_id + seq_len(n_new),
      step     = rep(k, n_new),
      size     = vapply(dedup_children, function(z) z$size, integer(1)),
      aic      = vapply(dedup_children, function(z) z$aic, numeric(1)),
      key      = vapply(dedup_children, function(z) z$key, character(1)),
      stringsAsFactors = FALSE
    )
    vars_list <- lapply(dedup_children, function(z) z$vars)
    frontier_k$vars <- vars_list
    
    model_id <- model_id + n_new
    
    frontiers[[k + 1L]] <- frontier_k
    all_rows[[length(all_rows) + 1L]] <- frontier_k
    
    current_frontier <- frontier_k
  }
  
  aic_by_model <- do.call(rbind, all_rows)
  
  meta <- list(
    family    = family,
    K         = K,
    eps       = eps,
    delta     = delta,
    L         = L,
    n         = n,
    p         = p,
    var_names = var_names
  )
  
  out <- list(
    frontiers    = frontiers,
    aic_by_model = aic_by_model,
    meta         = meta
  )
  class(out) <- "path_forest"
  out
}
                        #' Multi-path forward selection using AIC
#'
#' @param X Predictor matrix or data frame (n x p).
#' @param y Response vector of length n.
#' @param family Model family, "gaussian" (lm) or "binomial" (glm).
#' @param K Maximum number of steps (maximum model size). If NULL, uses min(p, 10).
#' @param eps Minimum AIC improvement required to expand a parent model.
#' @param delta AIC tolerance for keeping near-ties among children of the same parent.
#' @param L Optional cap on the number of models kept per level (after dedup).
#'
#' @return An object of class \code{"path_forest"}:
#' \itemize{
#'   \item \code{frontiers}: list of data frames, one per level (step 0 = empty model).
#'   \item \code{aic_by_model}: data frame summarizing all models across all steps.
#'   \item \code{meta}: list of metadata (family, K, eps, delta, L, variable names, etc.).
#' }
#' @export
build_paths <- function(X, y,
                        family = c("gaussian", "binomial"),
                        K = NULL,
                        eps = 1e-6,
                        delta = 0,
                        L = 50) {
  family <- match.arg(family)
  
  # Data preparation
  X <- as.data.frame(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }
  
  p <- ncol(X)
  n <- nrow(X)
  
  if (length(y) != n) stop("Length of y must match nrow(X).")
  if (is.null(K)) K <- min(p, 10L)
  
  # Combined data for fitting
  df <- data.frame(y = y, X, check.names = FALSE)
  var_names <- colnames(X)
  
  # Storage
  frontiers <- list()
  all_rows <- list() 
  model_id <- 0L
  
  # --- Step 0: Empty Model ---
  empty_info <- fit_aic(integer(0L), df = df, var_names = var_names, family = family)
  model_id <- model_id + 1L
  
  frontier0 <- data.frame(
    model_id = model_id,
    step     = 0L,
    size     = 0L,
    aic      = empty_info$aic,
    key      = model_key(integer(0L)),
    stringsAsFactors = FALSE
  )
  frontier0$vars <- list(integer(0L))
  
  frontiers[[1L]] <- frontier0
  all_rows[[1L]] <- frontier0
  
  current_frontier <- frontier0
  
  # --- Main Loop: Steps 1 to K ---
  for (k in seq_len(K)) {
    parent_frontier <- current_frontier
    n_parents <- nrow(parent_frontier)
    
    # Temporary storage for this step's children
    children_list <- list()
    child_counter <- 0L
    
    # Track best AIC per parent to apply 'delta' logic
    parent_best_aic <- rep(Inf, n_parents)
    parent_children_indices <- vector("list", n_parents)
    
    # 1. Expand every parent
    for (i in seq_len(n_parents)) {
      parent_vars <- parent_frontier$vars[[i]]
      parent_aic  <- parent_frontier$aic[i]
      
      available <- setdiff(seq_len(p), parent_vars)
      if (length(available) == 0L) next
      
      for (j in available) {
        child_vars <- sort(c(parent_vars, j))
        info <- fit_aic(child_vars, df = df, var_names = var_names, family = family)
        
        child_counter <- child_counter + 1L
        children_list[[child_counter]] <- list(
          parent_row = i,
          vars       = info$vars,
          aic        = info$aic,
          size       = length(info$vars),
          key        = model_key(info$vars)
        )
        
        parent_children_indices[[i]] <- c(parent_children_indices[[i]], child_counter)
        
        if (info$aic < parent_best_aic[i]) {
          parent_best_aic[i] <- info$aic
        }
      }
    }
    
    if (child_counter == 0L) break
    
    # 2. Filter children based on eps and delta
    keep_child <- rep(FALSE, child_counter)
    
    for (i in seq_len(n_parents)) {
      idxs <- parent_children_indices[[i]]
      if (length(idxs) == 0L) next
      
      parent_aic <- parent_frontier$aic[i]
      best_child <- parent_best_aic[i]
      
      # Rule: Best child must improve parent by at least 'eps'
      if (is.finite(best_child) && best_child <= parent_aic - eps) {
        # Rule: Keep any child within 'delta' of that best child
        for (idx in idxs) {
          if (children_list[[idx]]$aic <= best_child + delta) {
            keep_child[idx] <- TRUE
          }
        }
      }
    }
    
    if (!any(keep_child)) break
    
    kept_children <- children_list[keep_child]
    
    # 3. Deduplicate children by variable set key
    # If duplicates exist (reached via different parents), keep the one with better AIC (rare but possible)
    keys <- vapply(kept_children, function(z) z$key, character(1))
    unique_keys <- unique(keys)
    
    dedup_list <- vector("list", length(unique_keys))
    names(dedup_list) <- unique_keys
    
    for (ch in kept_children) {
      kkey <- ch$key
      existing <- dedup_list[[kkey]]
      if (is.null(existing) || ch$aic < existing$aic) {
        dedup_list[[kkey]] <- ch
      }
    }
    dedup_children <- Filter(Negate(is.null), dedup_list)
    
    # 4. Prune if > L models
    if (!is.null(L) && length(dedup_children) > L) {
      aics_vec <- vapply(dedup_children, function(z) z$aic, numeric(1))
      ord <- order(aics_vec)
      dedup_children <- dedup_children[ord[seq_len(L)]]
    }
    
    if (length(dedup_children) == 0L) break
    
    # 5. Finalize Step Data Frame
    n_new <- length(dedup_children)
    frontier_k <- data.frame(
      model_id = model_id + seq_len(n_new),
      step     = rep(k, n_new),
      size     = vapply(dedup_children, function(z) z$size, integer(1)),
      aic      = vapply(dedup_children, function(z) z$aic, numeric(1)),
      key      = vapply(dedup_children, function(z) z$key, character(1)),
      stringsAsFactors = FALSE
    )
    frontier_k$vars <- lapply(dedup_children, function(z) z$vars)
    
    model_id <- model_id + n_new
    frontiers[[k + 1L]] <- frontier_k
    all_rows[[length(all_rows) + 1L]] <- frontier_k
    current_frontier <- frontier_k
  }
  
  out <- list(
    frontiers = frontiers,
    aic_by_model = do.call(rbind, all_rows),
    meta = list(family = family, K = K, eps = eps, delta = delta, L = L)
  )
  class(out) <- "path_forest"
  out
}

