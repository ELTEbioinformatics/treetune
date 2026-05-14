test_that("combine_algorithms validates basic inputs", {
  tree_file <- tempfile(fileext = ".nwk")
  tree <- ape::rtree(20)
  ape::write.tree(tree, file = tree_file)

  expect_error(
    combine_algorithms(
      tree_file = tree_file,
      algorithms = "CPA",
      prune_percentages = c(4, 6)
    ),
    "`algorithms` must be a character vector of length 2",
    fixed = TRUE
  )

  expect_error(
    combine_algorithms(
      tree_file = tree_file,
      algorithms = c("CPA", "BAD"),
      prune_percentages = c(4, 6)
    ),
    "`algorithms` must contain only",
    fixed = TRUE
  )

  expect_error(
    combine_algorithms(
      tree_file = tree_file,
      algorithms = c("CPA", "IQR"),
      total_percent_remaining = 90,
      prune_percentages = c(5, 6)
    ),
    "`prune_percentages` must sum to 100 - `total_percent_remaining`",
    fixed = TRUE
  )

  expect_error(
    combine_algorithms(
      tree_file = "missing_file.nwk",
      algorithms = c("CPA", "IQR"),
      prune_percentages = c(4, 6)
    ),
    "`tree_file` does not exist",
    fixed = TRUE
  )
})



test_that("combine_algorithms passes algorithm-specific arguments", {
  tree_file <- tempfile(fileext = ".nwk")
  output_file <- tempfile(fileext = ".nwk")

  tree <- ape::rtree(100)
  ape::write.tree(tree, file = tree_file)

  observed_cpa_args <- list()
  observed_psfa_args <- list()

  mock_cpa <- function(
    tree_file,
    threshold,
    output,
    root_to_node_ratio = 0.1,
    beta = 25,
    radius_ratio = 0
  ) {
    observed_cpa_args$root_to_node_ratio <<- root_to_node_ratio
    observed_cpa_args$beta <<- beta
    observed_cpa_args$radius_ratio <<- radius_ratio

    tree <- ape::read.tree(tree_file)
    n_tips <- length(tree$tip.label)

    n_keep <- ceiling(n_tips * threshold / 100)
    tips_to_drop <- tree$tip.label[(n_keep + 1):n_tips]

    if (length(tips_to_drop) > 0) {
      tree <- ape::drop.tip(tree, tips_to_drop)
    }

    ape::write.tree(tree, file = output)

    list(
      tree = tree,
      percent_remaining = length(tree$tip.label) / n_tips * 100
    )
  }

  mock_psfa <- function(
    tree_file,
    threshold,
    output,
    longest_to_average = 9
  ) {
    observed_psfa_args$longest_to_average <<- longest_to_average

    tree <- ape::read.tree(tree_file)
    n_tips <- length(tree$tip.label)

    n_keep <- ceiling(n_tips * threshold / 100)
    tips_to_drop <- tree$tip.label[(n_keep + 1):n_tips]

    if (length(tips_to_drop) > 0) {
      tree <- ape::drop.tip(tree, tips_to_drop)
    }

    ape::write.tree(tree, file = output)

    list(
      tree = tree,
      percent_remaining = length(tree$tip.label) / n_tips * 100
    )
  }

  testthat::local_mocked_bindings(
    prune_tree_CPA = mock_cpa,
    prune_tree_PSFA = mock_psfa,
    prune_tree_IQR = function(tree_file, threshold, output) {
      stop("IQR should not be called in this test.")
    }
  )

  combine_algorithms(
    tree_file = tree_file,
    algorithms = c("CPA", "PSFA"),
    total_percent_remaining = 90,
    prune_percentages = c(4, 6),
    algorithm_args = list(
      CPA = list(
        root_to_node_ratio = 0.2,
        beta = 30,
        radius_ratio = 5
      ),
      PSFA = list(
        longest_to_average = 8
      )
    ),
    output = output_file
  )

  expect_equal(observed_cpa_args$root_to_node_ratio, 0.2)
  expect_equal(observed_cpa_args$beta, 30)
  expect_equal(observed_cpa_args$radius_ratio, 5)

  expect_equal(observed_psfa_args$longest_to_average, 8)
})
