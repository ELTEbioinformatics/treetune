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
  ntip <- ape::Ntip(tree)
  root_id <- ntip + 1L

  if (is.function(M_n)) {
    M_n_int <- max(1L, as.integer(floor(M_n(original_num_leaves))))
  } else {
    M_n_int <- if (identical(M_n, 0)) original_num_leaves else max(1L, as.integer(floor(M_n)))
  }

  hang <- function(tree, root_node) {
    dist_nodes <- ape::dist.nodes(tree)
    as.numeric(dist_nodes[root_node, 1:ape::Ntip(tree)])
  }

  edge_len_between <- function(tree, u, v) {
    e <- tree$edge
    el <- tree$edge.length
    k <- which(e[, 1] == u & e[, 2] == v)
    if (length(k) == 1) return(el[k])
    k <- which(e[, 1] == v & e[, 2] == u)
    if (length(k) == 1) return(el[k])
    stop("Could not find edge between nodes on path (unexpected).")
  }

  midpoint_internal_node <- function(tree) {
    ntip <- ape::Ntip(tree)
    if (ntip < 2) return(NA_integer_)

    d <- ape::cophenetic.phylo(tree)
    d[lower.tri(d, diag = TRUE)] <- -Inf
    ij <- which(d == max(d, na.rm = TRUE), arr.ind = TRUE)[1, ]
    tip_i <- ij[1]; tip_j <- ij[2]
    diam <- d[tip_i, tip_j]
    half <- diam / 2

    path_nodes <- ape::nodepath(tree, from = tip_i, to = tip_j)[[1]]
    if (length(path_nodes) < 2) return(NA_integer_)

    cum <- 0
    for (k in seq_len(length(path_nodes) - 1)) {
      u <- path_nodes[k]
      v <- path_nodes[k + 1]
      seg <- edge_len_between(tree, u, v)

      if (cum + seg >= half) {
        t_from_u <- half - cum
        t_from_v <- seg - t_from_u

        cand <- c(u, v)
        dist_to_mid <- c(t_from_u, t_from_v)

        is_internal <- cand > ntip
        if (!any(is_internal)) return(NA_integer_)


        internal_idx <- which(is_internal)
        best <- internal_idx[which.min(dist_to_mid[internal_idx])]
        chosen <- cand[best]

        if (chosen <= ntip) return(NA_integer_)
        return(as.integer(chosen))
      }

      cum <- cum + seg
    }

    chosen <- path_nodes[length(path_nodes)]
    if (chosen <= ntip) return(NA_integer_)
    as.integer(chosen)
  }

  children <- split(tree$edge[, 2], tree$edge[, 1])

  bfs_nodes <- function(start) {
    q <- start
    out <- integer(0)
    while (length(q) > 0) {
      v <- q[1]
      q <- q[-1]
      out <- c(out, v)
      ch <- children[[as.character(v)]]
      if (!is.null(ch)) q <- c(q, ch)
    }
    out
  }

  order_all <- bfs_nodes(root_id)
  internal_nodes <- order_all[order_all > ntip]
  num_non_leaf_nodes <- length(internal_nodes)

  num_roots <- min(
    max(floor(root_to_node_ratio * num_non_leaf_nodes), min_num_of_roots),
    num_non_leaf_nodes
  )

  set.seed(123)
  roots <- sample(internal_nodes, num_roots)

  mid_root <- midpoint_internal_node(tree)
  if (!is.na(mid_root)) roots <- unique(c(roots, mid_root))

  best_p_v <- Inf
  best_root <- NULL
  best_stop0 <- NULL
  best_original_radius <- NULL

  for (root in roots) {
    data <- hang(tree, root)
    original_radius <- max(data)

    a <- min(data)
    b <- original_radius
    n <- M_n_int


    bin_edges <- a + (0:n) * (b - a) / n


    idx <- findInterval(data, bin_edges, rightmost.closed = TRUE)
    idx <- idx[idx >= 1 & idx <= n]
    hist_counts <- tabulate(idx, nbins = n)

    cum_freq <- cumsum(hist_counts)
    total_freq <- cum_freq[length(cum_freq)]


    threshold_index0 <- which(cum_freq >= (threshold / 100) * total_freq)[1] - 1L


    stop0 <- n - 1L

    if (threshold_index0 == n - 1L) {
      stop0 <- threshold_index0
    } else if (threshold_index0 == n - 2L) {
      if (hist_counts[threshold_index0 + 2L] < (beta / 100) * hist_counts[threshold_index0 + 1L]) {
        stop0 <- threshold_index0
      }
    } else {
      for (i0 in threshold_index0:(n - 2L)) {
        if (hist_counts[i0 + 2L] < (beta / 100) * hist_counts[i0 + 1L]) {
          stop0 <- i0
          break
        }
      }
    }


    p_v <- bin_edges[stop0 + 2L]


    if (p_v < best_p_v) {
      best_p_v <- p_v
      best_root <- root
      best_stop0 <- stop0
      best_original_radius <- original_radius
    }
  }

  if (show_plot && !is.null(best_root)) {
    data <- hang(tree, best_root)
    a <- min(data)
    b <- max(data)
    n <- M_n_int
    bin_edges <- a + (0:n) * (b - a) / n

    idx <- findInterval(data, bin_edges, rightmost.closed = TRUE)
    idx <- idx[idx >= 1 & idx <= n]
    counts <- tabulate(idx, nbins = n)

    mids <- (bin_edges[-1] + bin_edges[-(n + 1L)]) / 2
    df <- data.frame(
      mids = mids,
      counts = counts,
      color = ifelse(seq_len(n) > (best_stop0 + 1L), "red", "gray")
    )

    plot(
      ggplot2::ggplot(df, ggplot2::aes(x = mids, y = counts, fill = color)) +
        ggplot2::geom_bar(stat = "identity", color = "black") +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = "CPA", x = "Distance from root", y = "Frequency")
    )
  }

  if ((100 - 100 * best_p_v / best_original_radius) >= radius_ratio) {
    pruned <- character(0)
    retained <- character(0)

    tip_distances <- hang(tree, best_root)
    for (i in seq_len(ntip)) {
      nm <- tree$tip.label[i]
      dist <- tip_distances[i]
      if (dist <= best_p_v || nm %in% safe_tips) {
        retained <- c(retained, nm)
      } else {
        pruned <- c(pruned, nm)
      }
    }

    tree <- ape::keep.tip(tree, retained)
    ape::write.tree(tree, file = output)

    percent_remaining <- 100 * length(tree$tip.label) / original_num_leaves
    if (show_pruned_tips) {
      return(list(tree = tree, percent_remaining = percent_remaining, pruned_tips = pruned))
    }
    return(list(tree = tree, percent_remaining = percent_remaining))
  }

  if (show_pruned_tips) {
    return(list(tree = tree, percent_remaining = 100, pruned_tips = character(0)))
  }
  list(tree = tree, percent_remaining = 100)
}






