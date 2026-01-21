test_that("CPA prunes or retains and writes output", {
  f_in  <- .tmp_treefile(70)
  f_out <- tempfile(fileext = ".nwk")
  original_n <- .read_ntip(f_in)

  res <- prune_tree_CPA(
    tree_file = f_in,
    root_to_node_ratio = 0.1,
    min_num_of_roots = 10,
    threshold = 90,
    beta = 20,
    radius_ratio = 0,        # allow pruning
    safe_tips = character(0),
    show_plot = FALSE,
    show_pruned_tips = TRUE,
    output = f_out
  )

  .expect_phylo_with_percent(res, original_n)
  expect_true(file.exists(f_out))
  expect_type(res$pruned_tips, "character")
})


test_that("CPA prunes or retains and writes output when M_n is supplied", {
  f_in  <- .tmp_treefile(70)
  f_out <- tempfile(fileext = ".nwk")
  original_n <- .read_ntip(f_in)

  res <- prune_tree_CPA(
    tree_file = f_in,
    root_to_node_ratio = 0.1,
    min_num_of_roots = 10,
    M_n = function(x) {x/2},
    threshold = 90,
    beta = 20,
    radius_ratio = 0,        # allow pruning
    safe_tips = character(0),
    show_plot = FALSE,
    show_pruned_tips = TRUE,
    output = f_out
  )

  .expect_phylo_with_percent(res, original_n)
  expect_true(file.exists(f_out))
  expect_type(res$pruned_tips, "character")
})


test_that("CPA respects safe_tips (never prunes them)", {
  # Build a tree, mark a couple of tips as 'safe'
  tr <- ape::rtree(50)
  f_in <- tempfile(fileext = ".nwk")
  ape::write.tree(tr, f_in)

  safe <- tr$tip.label[1:3]

  res <- prune_tree_CPA(
    tree_file = f_in,
    show_plot = FALSE,
    show_pruned_tips = TRUE,
    safe_tips = safe,
    output = tempfile(fileext = ".nwk")
  )

  expect_false(any(safe %in% res$pruned_tips))
  expect_true(all(safe %in% res$tree$tip.label))
})

test_that("CPA with radius_ratio too strict results in no pruning", {
  f_in <- .tmp_treefile(60)
  res  <- prune_tree_CPA(
    tree_file = f_in,
    radius_ratio = 100,   # extremely strict
    show_plot = FALSE,
    show_pruned_tips = TRUE,
    output = tempfile(fileext = ".nwk")
  )
  # percent_remaining should be 100 when no pruning occurs
  expect_equal(res$percent_remaining, 100)
})



test_that("CPA produces a plot", {
  f_in <- .tmp_treefile(60)
  expect_silent(
    prune_tree_CPA(
    tree_file = f_in,
    show_plot = T,
    show_pruned_tips = TRUE,
    output = tempfile(fileext = ".nwk")
    )
  )
})



test_that("CPA does not show pruned tips if asked", {
  f_in <- .tmp_treefile(60)
  res  <- prune_tree_CPA(
      tree_file = f_in,
      show_pruned_tips = F,
      output = tempfile(fileext = ".nwk")
      )
  expect_type(res$pruned_tips, "NULL")
})
