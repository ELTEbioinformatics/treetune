#' Prune a Phylogenetic Tree Using the Primitive Straight-Forward Approach (PSFA)
#'
#' This functions Identifies the longest edge and calculates the sizes of the two components resulting from its removal.
#' If the size of the smaller component is within the tolerance range and the edge was excessively long, removes it.
#' @param tree_file character; a phylogenetic tree in .nwk, .newick or .tree format.
#' @param threshold numeric; the minimum percentage of the leaves that the algorithm must keep. Defaults to 90.
#' @param longest_to_average numeric; determines how many times longer should a branch be than the average, to consider it too long. Defaults to 9.
#' @param output character; the desired name of the output file. Defaults to "psfa_tree.nwk".
#' @return A list with two elements:
#' \describe{
#'   \item{tree}{The pruned tree in phylo format.}
#'   \item{percent_remaining}{A numeric value indicating the percentage of leaves retained after pruning.}
#' }
#' @details You have to define criteria for considering an edge as "excessively long".
#' You also have to specify the tolerance range.
#' @export


prune_tree_PSFA <- function(tree_file, threshold = 90, longest_to_average = 9, output = "psfa_tree.nwk") {

  # Validation
  .check_string(tree_file, "tree_file")
  .check_file_exists(tree_file)
  .check_number(threshold, "threshold", 0, 100)
  .check_number(longest_to_average, "longest_to_average", 0)
  .check_string(output, "output")

  tree <- ape::read.tree(tree_file)
  original_num_leaves <- length(tree$tip.label)
  percent_left <- 100 - threshold

  repeat {
    result <- one_primitive_step(tree, original_num_leaves, percent_left, longest_to_average)
    tree <- result$tree
    percent_left <- result$percent_left
    if (!result$PRUNE) break
  }

  ape::write.tree(tree, file = output)
  current_num_leaves <- length(tree$tip.label)
  list(tree = tree, percent_remaining = 100 * current_num_leaves / original_num_leaves)
}
