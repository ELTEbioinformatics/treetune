test_that("PSFA runs and returns a pruned tree (or same) with sane percent", {
  f_in  <- .tmp_treefile(60)
  f_out <- tempfile(fileext = ".nwk")
  original_n <- .read_ntip(f_in)

  res <- prune_tree_PSFA(
    tree_file = f_in,
    threshold = 10,
    longest_to_average = 5,
    output = f_out
  )

  .expect_phylo_with_percent(res, original_n)
  expect_true(file.exists(f_out))

  # output file is valid Newick
  tr_out <- ape::read.tree(f_out)
  expect_s3_class(tr_out, "phylo")
})

test_that("PSFA threshold bounds behave", {
  f_in  <- .tmp_treefile(40)
  original_n <- .read_ntip(f_in)

  # threshold = 0 → shouldn't prune due to budget (may still early-stop)
  res0 <- prune_tree_PSFA(f_in, threshold = 0, output = tempfile(fileext = ".nwk"))
  expect_true(ape::Ntip(res0$tree) <= original_n)

  # threshold = 90 gives room to prune
  res90 <- prune_tree_PSFA(f_in, threshold = 90, output = tempfile(fileext = ".nwk"))
  expect_true(ape::Ntip(res90$tree) <= original_n)
})
