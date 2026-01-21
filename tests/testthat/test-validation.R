test_that("Validation errors are informative (PSFA)", {
  expect_error(prune_tree_PSFA(0), "must be a non-empty")
  expect_error(prune_tree_PSFA("nope.nwk"), "does not exist")
  tf <- tempfile(fileext = ".nwk"); file.create(tf)
  expect_error(prune_tree_PSFA(tf, threshold = 200), "numeric value in")
  expect_error(prune_tree_PSFA(tf, longest_to_average = -1), "numeric value in")
})

test_that("Validation errors are informative (IQR)", {
  expect_error(prune_tree_IQR(0), "must be a non-empty")
  expect_error(prune_tree_IQR("missing.nwk"), "does not exist")
  tf <- tempfile(fileext = ".nwk"); file.create(tf)
  expect_error(prune_tree_IQR(tf, threshold = -5), "numeric value in")
})

test_that("Validation errors are informative (CPA)", {
  tf <- tempfile(fileext = ".nwk"); file.create(tf)
  expect_error(prune_tree_CPA(0), "must be a non-empty")
  expect_error(prune_tree_CPA(tf, root_to_node_ratio = -0.1), "numeric value in")
  expect_error(prune_tree_CPA(tf, min_num_of_roots = 0), "numeric value in")
  expect_error(prune_tree_CPA(tf, M_n = "function"), "must be a function")
  expect_error(prune_tree_CPA(tf, beta = 200), "numeric value in")
  expect_error(prune_tree_CPA(tf, show_plot = "yes"), "must be either TRUE or FALSE")
  expect_error(prune_tree_CPA(tf, safe_tips = 1:3), "character vector")
})
