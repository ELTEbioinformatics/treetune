test_that("IQR fence returns Inf for <4 and a number otherwise", {
  expect_equal(calculate_iqr_threshold(c(1,2,3)), Inf)
  th <- calculate_iqr_threshold(1:10)
  expect_true(is.finite(th))
})

test_that("ensure_node_labels fills node labels", {
  tr <- ape::rtree(20)
  tr$node.label <- NULL
  tr2 <- ensure_node_labels(tr)
  expect_equal(length(tr2$node.label), tr2$Nnode)
  expect_true(all(nzchar(tr2$node.label)))
})
