#' Prune a Phylogenetic Tree Using the Circular Pruning Algorithm (CPA)
#'
#' This function implements a more complex (and often better) algorithm than PSFA but it also has a longer runtime.
#' It aims to prune the tree to make it roughly circular
#' @param tree_file character; a phylogenetic tree in .nwk, .newick or .tree format.
#' @param root_to_node_ratio numeric; determines how many nodes does the algorithm try out as "roots".
#' By increasing this parameter, the precision and the runtime becomes greater too.
#' The default value means that the algorithm will try out a randomly selected 10% of the non-leaf nodes as "roots".
#' Apart from these random nodes the algorithm will always try the "midpoint root" of the tree, as it usually gives a good result.
#' The parameter must be between 0 and 1, although the algorithm can handle the other cases as well (as long as it is a number). Defaults to 0.1.
#' @param min_num_of_roots integer; the minimum number of nodes that the algorithm tries out as "roots".
#' This parameter only matters if the root_to_node_ratio chooses less nodes than min_num_of_roots.
#' It must be less than the number of non-leaf nodes, although the algorithm can handle the other cases as well (as long as it is an integer).
#' Defaults to 15.
#' @param M_n function; a function, which you should define before calling the algorithm if you intend to use this parameter.
#' During the CPA the leaves are being placed in brackets based on their distances from the "root".
#' The algorithm will delete the leaves in the last couple of brackets if certain conditions are met.
#' M_n determines the number of brackets: the total number of brackets will be M_n(num_leaves), where num_leaves is the number of leaves on the tree.
#' The default value means that the user did not give any function to the argument. In this case M_n(num_leaves) = num_leaves.
#' Assigning anything to M_n other than a function will result in an error message.
#' @param threshold numeric; the minimum percentage of the leaves that the algorithm must keep. Defaults to 90.
#' @param beta numeric; determines where to prune the tree. After placing the leaves in the brackets, the first couple of brackets will be safe, until the total number of leaves in them reaches alpha*num_leaves.
#' When the CPA finds this threshold, it is ready for pruning. The algorithm will start the pruning if a bracket contains beta percent less leaves than the previous one.
#' After this point the CPA will delete every leaf. Defaults to 25.
#' @param radius_ratio numeric;  controls whether pruning occurs based on the tree's radius. If the longest root-to-tip distance in the pruned tree is not at least radius_ratio percent smaller than in the original tree (measured from the same root), no pruning is performed.
#' Defaults to 0.
#' @param safe_tips list; a list where you can input the names of the tips that you do not want to be cut off.
#' After the algorithm chose the best root and started the pruning, it checks which leaves are in the list and leaves them on the tree.
#' The default value means that there are no safe leaves. It is important that the names in the list must be characters.
#' @param show_plot bool; the default value means that after running the code it will not show any plot. By setting this value to TRUE, the brackets will be shown to the user.
#' These brackets correspond to the "best root" (the center of the outscribed circle of the tree). On the x axis you can see the distances from the root, while on the y axis you can see how many leaves were placed in each bracket.
#' Blue indicates the preserved brackets, while red indicates the ones which have been deleted.
#' @param show_pruned_tips bool; responsible for indicating which leaves have been cut from the tree. By setting this value to TRUE, the function will return a list that contains the names of the pruned tips, in addition to the default return values.
#' This list will be the third return value.
#' @param output character; the desired name of the output file. Defaults to "cpa_tree.nwk".
#' @return A list with the following elements:
#' \describe{
#'   \item{tree}{The pruned tree in phylo format.}
#'   \item{percent_remaining}{A numeric value indicating the percentage of leaves retained after pruning.}
#'   \item{pruned_tips}{(Optional) A character vector containing the names of the pruned tips. This is returned only if `show_pruned_tips = TRUE`.}
#' }
#' @details You have to define criteria for considering an edge as "excessively long".
#' You also have to specify the tolerance range.
#' @export


prune_tree_CPA <- function(
    tree_file,
    root_to_node_ratio = 0.1,
    min_num_of_roots = 15,
    M_n = 0,
    threshold = 90,
    beta = 20,
    radius_ratio = 0,
    safe_tips = character(0),
    show_plot = FALSE,
    show_pruned_tips = FALSE,
    output = "cpa_tree.nwk"
) {

  # Validation
  .check_string(tree_file, "tree_file")
  .check_file_exists(tree_file)
  .check_number(root_to_node_ratio, "root_to_node_ratio", 0, 1)
  .check_number(min_num_of_roots, "min_num_of_roots", 1)
  .check_function(M_n, "M_n")
  .check_number(threshold, "threshold", 0, 100)
  .check_number(beta, "beta", 0, 100)
  .check_number(radius_ratio, "radius_ratio", 0, 100)
  .check_character_vector(safe_tips, "safe_tips")
  .check_boolean(show_plot, "show_plot")
  .check_boolean(show_pruned_tips, "show_pruned_tips")
  .check_string(output, "output")



  tree <- ape::read.tree(tree_file)
  original_num_leaves <- length(tree$tip.label)

  if (is.function(M_n)) {
    M_n_int <- max(1, floor(M_n(original_num_leaves)))
  } else {
    M_n_int <- ifelse(M_n == 0, original_num_leaves, max(1, floor(M_n)))
  }

  hang <- function(tree, root_node) {
    dist_nodes <- ape::dist.nodes(tree)
    tips <- which(tree$edge[, 2] <= length(tree$tip.label))
    dist_to_root <- dist_nodes[root_node, tips]
    as.numeric(dist_to_root)
  }

  internal_nodes <- unique(tree$edge[, 1])
  num_non_leaf_nodes <- length(internal_nodes)
  num_roots <- min(max(floor(root_to_node_ratio * num_non_leaf_nodes), min_num_of_roots), num_non_leaf_nodes)
  set.seed(123)
  roots <- sample(internal_nodes, num_roots)
  roots <- unique(c(roots, ape::Nnode(tree) + 1))

  best_p_v <- Inf
  best_data <- NULL
  best_root <- NULL
  best_stop <- NULL
  best_original_radius <- NULL

  for (root in roots) {
    data <- hang(tree, root)
    original_radius <- max(data)
    bin_edges <- seq(min(data), original_radius, length.out = M_n_int + 1)
    hist_counts <- graphics::hist(data, breaks = bin_edges, plot = FALSE)$counts

    cum_freq <- cumsum(hist_counts)
    total_freq <- sum(hist_counts)
    threshold_index <- which(cum_freq >= (threshold / 100) * total_freq)[1]
    stop <- M_n_int - 1

    if (threshold_index == M_n_int - 1) {
      stop <- threshold_index
    } else if (threshold_index == M_n_int - 2) {
      if (hist_counts[threshold_index + 1] < (beta / 100) * hist_counts[threshold_index]) {
        stop <- threshold_index
      }
    } else {
      for (i in threshold_index:(M_n_int - 2)) {
        if (hist_counts[i + 1] < (beta / 100) * hist_counts[i]) {
          stop <- i
          break
        }
      }
    }

    p_v <- bin_edges[stop + 2]
    if (p_v < best_p_v) {
      best_p_v <- p_v
      best_root <- root
      best_stop <- stop
      best_original_radius <- original_radius
    }
  }

  if (show_plot && !is.null(best_root)) {
    data <- hang(tree, best_root)
    bin_edges <- seq(min(data), max(data), length.out = M_n_int + 1)
    hist_obj <- graphics::hist(data, breaks = bin_edges, plot = FALSE)
    df <- data.frame(
      mids = hist_obj$mids,
      counts = hist_obj$counts,
      color = ifelse(seq_along(hist_obj$counts) > best_stop, "red", "gray")
    )
    plot(ggplot2::ggplot(df, ggplot2::aes(x = mids, y = counts, fill = color)) +
           ggplot2::geom_bar(stat = "identity", color = "black") +
           ggplot2::scale_fill_identity() +
           ggplot2::labs(title = "CPA", x = "Distance from root", y = "Frequency"))
  }

  if ((100 - 100 * best_p_v / best_original_radius) >= radius_ratio) {
    pruned <- c()
    retained <- c()
    tip_ids <- 1:ape::Ntip(tree)
    tip_distances <- hang(tree, best_root)
    for (i in seq_along(tip_ids)) {
      name <- tree$tip.label[i]
      dist <- tip_distances[i]
      if (dist <= best_p_v || name %in% safe_tips) {
        retained <- c(retained, name)
      } else {
        pruned <- c(pruned, name)
      }
    }
    tree <- ape::keep.tip(tree, retained)
    ape::write.tree(tree, file = output)

    if (show_pruned_tips) {
      return(list(tree = tree, percent_remaining = 100 * length(tree$tip.label) / original_num_leaves, pruned_tips = pruned))
    } else {
      return(list(tree = tree, percent_remaining = 100 * length(tree$tip.label) / original_num_leaves))
    }
  }

  if (show_pruned_tips) {
    return(list(tree = tree, percent_remaining = 100, pruned_tips = c()))
  }
  list(tree = tree, percent_remaining = 100)
}



