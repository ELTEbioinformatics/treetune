test_that("IQR orchestrator returns valid tree and file", {
  f_in  <- .tmp_treefile(80)
  f_out <- tempfile(fileext = ".nwk")
  original_n <- .read_ntip(f_in)

  res <- prune_tree_IQR(
    tree_file = f_in,
    threshold = 80,                     # large budget (20% remain target)
    output = f_out
  )

  .expect_phylo_with_percent(res, original_n)
  expect_true(file.exists(f_out))
  expect_s3_class(ape::read.tree(f_out), "phylo")
})

test_that("IQR handles very small trees without crashing", {
  # With <4 tips the IQR fence returns Inf; should no-op safely
  f_in <- .tmp_treefile(3)
  res  <- prune_tree_IQR(f_in, threshold = 90, output = tempfile(fileext = ".nwk"))
  expect_s3_class(res$tree, "phylo")
  expect_equal(ape::Ntip(res$tree), .read_ntip(f_in))  # likely unchanged
})
