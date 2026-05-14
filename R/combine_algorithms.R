#' Combine two phylogenetic tree pruning algorithms
#'
#' Runs two tree pruning algorithms sequentially.
#'
#' @param tree_file character; a phylogenetic tree in .nwk, .newick or .tree format.
#' @param algorithms character vector of length 2. Allowed values are "CPA", "IQR", and "PSFA".
#' @param total_percent_remaining numeric; final percentage of tips to retain after both algorithms.
#' @param prune_percentages numeric vector of length 2; percentage of original tips to prune by each algorithm.
#'   For example, c(4, 6) means prune 4 percent with the first algorithm and 6 percent with the second.
#' @param algorithm_args named list; optional algorithm-specific arguments.
#'   See all possible arguments in the help section of each argument.
#' @param output character; desired name of the final output file.
#'
#' @return A list containing the final tree, percentage retained compared to the original, and step-wise results.
#'   The step-wise results for the first show the percentage of remaining tips compared to the original.
#'   For the second step, it shows the percentage of remaining tips compared to the first step.
#'
#' @export



combine_algorithms <- function(
    tree_file,
    algorithms,
    total_percent_remaining = 90,
    prune_percentages,
    algorithm_args = list(),
    output = "combined_tree.nwk"
) {
  if (!is.character(tree_file) || length(tree_file) != 1) {
    stop("`tree_file` must be a single character string.")
  }

  if (!file.exists(tree_file)) {
    stop("`tree_file` does not exist: ", tree_file)
  }

  if (!is.character(algorithms) || length(algorithms) != 2) {
    stop("`algorithms` must be a character vector of length 2.")
  }

  algorithms <- toupper(algorithms)

  allowed_algorithms <- c("CPA", "IQR", "PSFA")

  if (!all(algorithms %in% allowed_algorithms)) {
    stop(
      "`algorithms` must contain only: ",
      paste(allowed_algorithms, collapse = ", ")
    )
  }

  if (!is.numeric(total_percent_remaining) || length(total_percent_remaining) != 1) {
    stop("`total_percent_remaining` must be a single numeric value.")
  }

  if (total_percent_remaining <= 0 || total_percent_remaining > 100) {
    stop("`total_percent_remaining` must be > 0 and <= 100.")
  }

  if (!is.numeric(prune_percentages) || length(prune_percentages) != 2) {
    stop("`prune_percentages` must be a numeric vector of length 2.")
  }

  if (any(prune_percentages < 0)) {
    stop("`prune_percentages` cannot contain negative values.")
  }

  expected_remaining <- 100 - sum(prune_percentages)

  if (!isTRUE(all.equal(expected_remaining, total_percent_remaining))) {
    stop(
      "`prune_percentages` must sum to 100 - `total_percent_remaining`.\n",
      "sum(prune_percentages) = ", sum(prune_percentages), "\n",
      "expected total_percent_remaining = ", expected_remaining
    )
  }

  if (!is.list(algorithm_args)) {
    stop("`algorithm_args` must be a named list.")
  }

  first_remaining <- 100 - prune_percentages[1]

  thresholds <- c(
    first_remaining,
    total_percent_remaining / first_remaining * 100
  )

  names(thresholds) <- paste0(algorithms, "_step", seq_along(algorithms))


  algorithm_functions <- list(
    CPA = prune_tree_CPA,
    IQR = prune_tree_IQR,
    PSFA = prune_tree_PSFA
  )

  intermediate_file <- tempfile(fileext = ".nwk")

  on.exit({
    if (file.exists(intermediate_file)) {
      unlink(intermediate_file)
    }
  }, add = TRUE)

  step_results <- vector("list", length = 2)
  names(step_results) <- paste0(algorithms, "_step", seq_along(algorithms))

  current_tree_file <- tree_file

  for (i in seq_along(algorithms)) {
    alg <- algorithms[i]
    fun <- algorithm_functions[[alg]]

    step_output <- if (i == 1) intermediate_file else output

    args <- list(
      tree_file = current_tree_file,
      threshold = thresholds[i],
      output = step_output
    )

    extra_args <- algorithm_args[[alg]]

    if (is.null(extra_args)) {
      extra_args <- list()
    }

    extra_args <- validate_algorithm_args(alg, extra_args)

    args <- utils::modifyList(args, extra_args)

    step_results[[i]] <- do.call(fun, args)

    if (is.null(step_results[[i]]$tree)) {
      stop("Algorithm ", alg, " did not return a `tree` element.")
    }

    ape::write.tree(step_results[[i]]$tree, file = step_output)

    current_tree_file <- step_output
  }

  original_tree <- ape::read.tree(tree_file)
  original_n_tips <- length(original_tree$tip.label)

  final_n_tips <- length(step_results[[2]]$tree$tip.label)

  percent_remaining_original <- final_n_tips / original_n_tips * 100

  list(
    tree = step_results[[2]]$tree,
    percent_remaining_original = percent_remaining_original,
    #percent_remaining_last_step = step_results[[2]]$percent_remaining,
    #algorithms = algorithms,
    #thresholds = thresholds,
    #prune_percentages_original = prune_percentages,
    step_results = step_results,
    output = output
  )
}
