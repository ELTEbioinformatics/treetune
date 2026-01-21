# Helper to write a random tree to a temp Newick file and return the path
.tmp_treefile <- function(n = 40) {
  tr <- ape::rtree(n)
  f <- tempfile(fileext = ".nwk")
  ape::write.tree(tr, f)
  f
}

# Read a Newick file and return Ntip
.read_ntip <- function(path) {
  ape::Ntip(ape::read.tree(path))
}

# Expect a valid phylo object with sensible tip counts
.expect_phylo_with_percent <- function(x, original_n) {
  testthat::expect_s3_class(x$tree, "phylo")
  testthat::expect_true(is.numeric(x$percent_remaining))
  testthat::expect_true(x$percent_remaining > 0)
  testthat::expect_true(x$percent_remaining <= 100)
  # Internal consistency: percent ~ Ntip
  testthat::expect_equal(
    x$percent_remaining,
    100 * ape::Ntip(x$tree) / original_n,
    tolerance = 1e-8
  )
}
