#' Visualize branching by step in a path_forest
#'
#' This helper summarizes how many distinct models are present
#' on each frontier (step) of the multi-path forward search.
#' A larger number of models at a given step corresponds to
#' more branching at that depth of the search.
#'
#' @param forest An object of class \code{"path_forest"} returned by \code{build_paths()}.
#'
#' @return Invisibly returns a data frame with columns \code{step} and
#'   \code{n_models}. Called for its side-effect of producing a plot.
#' @export
plot_branching_by_step <- function(forest) {
  if (!inherits(forest, "path_forest")) {
    stop("plot_branching_by_step() expects an object of class 'path_forest'.")
  }
  
  # Stack all frontiers into a single data frame
  frontier_df <- do.call(rbind, forest$frontiers)
  # frontier_df has columns: model_id, step, size, aic, key, vars
  
  # Count how many models at each step
  tab <- as.data.frame(table(frontier_df$step), stringsAsFactors = FALSE)
  names(tab) <- c("step", "n_models")
  tab$step <- as.integer(tab$step)
  tab <- tab[order(tab$step), ]
  
  # Basic base-R line plot
  plot(
    tab$step, tab$n_models, type = "b",
    xlab = "Step",
    ylab = "Number of models in frontier",
    main = "Branching by step in multi-path forward selection"
  )
  
  invisible(tab)
}
