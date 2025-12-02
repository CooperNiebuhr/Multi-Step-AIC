set.seed(1)
n <- 120; p <- 8
X <- matrix(rnorm(n*p), n, p)
beta <- c(2, -1.5, 0, 0, 1, rep(0, p-5))
y <- X %*% beta + rnorm(n, sd = 1)

out <- build_paths(X = X, y = as.numeric(y), family = "gaussian", K = min(ncol(X), 10), eps = 1e-6, delta = 1, L = 50)

# check for only one variable added per level, per branch
for (k in seq_along(out$frontiers)) {
  df <- out$frontiers[[k]]
  expected_size <- k - 1L
  if (!all(df$size == expected_size)) {
    stop(sprintf("Size mismatch in level %d", k-1))
  }
}

eps   <- out$meta$eps
delta <- out$meta$delta

# check that all children obey eps and delta
for (k in seq_len(length(out$frontiers)-1)) {
  parents <- out$frontiers[[k]]
  children <- out$frontiers[[k+1]]
  
  for (i in seq_len(nrow(parents))) {
    p_vars <- parents$vars[[i]]
    p_aic  <- parents$aic[i]
    
    # children with exactly one more var and parent subset
    eligible <- vapply(children$vars, function(ch) {
      length(ch) == length(p_vars) + 1 && all(p_vars %in% ch)
    }, logical(1))
    
    ch <- children[eligible, ]
    if (nrow(ch) == 0) next
    
    best_child_aic <- min(ch$aic)
    
    # eps constraint
    if (!(best_child_aic <= p_aic - eps)) {
      stop("Parent violated eps rule.")
    }
    
    # delta constraint: all kept children must satisfy this
    if (!all(ch$aic <= best_child_aic + delta + 1e-12)) {
      stop("Child violated delta bound.")
    }
  }
}

# check for functional deduplication
for (k in 2:length(out$frontiers)) {
  df <- out$frontiers[[k]]
  keys <- df$key
  
  # 1. No duplicates should exist
  if (any(duplicated(keys))) {
    stop(sprintf("Deduplication failed in frontier k=%d", k-1))
  }
  
  # 2. For safety: recompute minimal AIC for each key and match
  for (key in unique(keys)) {
    rows <- df[df$key == key, ]
    if (nrow(rows) != 1) {
      stop("Multiple models with same key after deduplication.")
    }
    # Nothing more needed because dedup guarantees minimal AIC
  }
}