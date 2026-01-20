#' Prune a Phylogenetic Tree Using the IQR Algorithm
#' This function implements a two-part tree pruning algorithm based on branch length
#' and root-to-tip distance outlier detection using the Inter Quartile Range (IQR) method.
#' @param tree_file character; a phylogenetic tree in .nwk, .newick or .tree format.
#' @param threshold numeric; the minimum percentage of the leaves that the algorithm must keep. Defaults to 90.
#' @param output character; the desired name of the output file. Defaults to "iqr_tree.nwk".
#' @return A list with two elements:
#' \describe{
#'   \item{tree}{The pruned tree in phylo format.}
#'   \item{percent_remaining}{A numeric value indicating the percentage of leaves retained after pruning.}
#' }
#' @details
#' Abbreviations:
#'
#' - IQR: the upper fence of the Inter Quartile Range of a vector: Q3 + 3 * (Q3 - Q1) = extreme outlier threshold.
#' - R2T: root-to-tip distance on a midpoint rooted tree.
#'
#' First part of the algorithm: Pruning the unrooted tree based on the upper fence of IQR branch lengths.
#'
#' Calculating the minimum number of tips for each branch (unrooted bifurcating tree, both directions are looked up, than the least tip number is chosen to represent the branch).
#'
#' Calculating the IQR for branch lengths.
#'
#' Excluding those extreme outlier branches containing less than 0.05 (or a given) proportion of the tips.
#'
#' Second part of the algorithm: Pruning the midpoint rooted tree based on root-to-tip distances. Excluding tips based on the R2T IQRs, while iteratively midpoint rooting and excluding the top greatest extreme outlier.
#'
#' Midpoint rooting the tree (after pruning with method described in the first part).
#'
#' Calculating the IQR for root-to-tip distances.
#'
#' Excluding the most extreme outlier tip based on the IQR for root-to-tip distances.
#'
#' Repeating the previous three steps until there are no more extreme outlier IQR tip is found.
#'
#' Unrooting the pruned tree.
#' @export


prune_tree_IQR <- function(tree_file, threshold = 90, output = "iqr_tree.nwk") {

  # Validation
  .check_string(tree_file, "tree_file")
  .check_file_exists(tree_file)
  .check_number(threshold, "threshold", 0, 100)
  .check_string(output, "output")


  tree <- ape::read.tree(tree_file)
  original_num_leaves <- ape::Ntip(tree)

  res1 <- prune_by_branch_length(tree, threshold)
  tree1 <- res1$tree
  percent_left <- res1$percent_left_budget

  tree2 <- prune_by_root_to_tip(tree1, percent_left, original_num_leaves)

  ape::write.tree(tree2, file = output)
  list(tree = tree2, percent_remaining = 100 * ape::Ntip(tree2) / original_num_leaves)
}
