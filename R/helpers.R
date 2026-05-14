#' @importFrom utils globalVariables
NULL

#' @importFrom stats cor.test lm optimize runif sd
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("mids", "counts", "color"))
}



# --- Validators ---
#' @keywords internal
.check_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower || x > upper) {
    stop(sprintf("`%s` must be a numeric value in [%s, %s].", name, lower, upper), call. = FALSE)
  }
}

#' @keywords internal
.check_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || x == "") {
    stop(sprintf("`%s` must be a non-empty character string.", name), call. = FALSE)
  }
}

#' @keywords internal
.check_file_exists <- function(path) {
  if (!file.exists(path)) stop(sprintf("File '%s' does not exist.", path), call. = FALSE)
}

#' @keywords internal
.check_character_vector <- function(x, name = deparse(substitute(x))) {
  if (!(is.character(x) && is.vector(x))) {
    stop(sprintf("`%s` must be a character vector.", name), call. = FALSE)
  }
  if (any(is.na(x))) {
    stop(sprintf("`%s` contains missing values.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.check_boolean <- function(x, name = deparse(substitute(x))) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be either TRUE or FALSE (a single logical value).", name),
         call. = FALSE)
  }
  invisible(TRUE)
}


#' @keywords internal
.check_function <- function(x, name = deparse(substitute(x))) {
  if (is.function(x)) {
    return(invisible(TRUE))
  }
  if (is.numeric(x) && length(x) == 1L && !is.na(x) && x == 0) {
    return(invisible(TRUE))
  }
  stop(sprintf("`%s` must be a function.", name),
       call. = FALSE)
}

# --- PSFA ---
#' @keywords internal
one_primitive_step <- function(tree, original_num_leaves, percent_left, longest_to_average) {
  edges <- tree$edge.length
  longest_branch_length <- max(edges)
  longest_branch_index <- which.max(edges)

  average_branch_length <- mean(edges)
  tip_ids <- 1:ape::Ntip(tree)
  num_leaves <- length(tip_ids)

  edge_node <- tree$edge[longest_branch_index, 2]
  leaves_on_one_side <- phangorn::Descendants(tree, edge_node, type = "tips")[[1]]
  other_side <- setdiff(tip_ids, leaves_on_one_side)

  PRUNE <- FALSE
  if (longest_branch_length > longest_to_average * average_branch_length) {
    if ((percent_left * original_num_leaves) / 100 > length(leaves_on_one_side)) {
      tree <- ape::drop.tip(tree, leaves_on_one_side)
      percent_left <- percent_left - 100 * length(leaves_on_one_side) / original_num_leaves
      PRUNE <- TRUE
    } else if ((percent_left * original_num_leaves) / 100 > length(other_side)) {
      tree <- ape::drop.tip(tree, other_side)
      percent_left <- percent_left - 100 * length(other_side) / original_num_leaves
      PRUNE <- TRUE
    }
  }
  list(tree = tree, PRUNE = PRUNE, percent_left = percent_left)
}

#' @keywords internal
calculate_iqr_threshold <- function(values) {
  values <- as.numeric(values)
  if (length(values) < 4L) return(Inf)
  qs <- as.numeric(stats::quantile(values, probs = c(0.25, 0.75), type = 7, names = FALSE))
  q1 <- qs[1]; q3 <- qs[2]
  q3 + 3 * (q3 - q1)
}

#' @keywords internal
root_like_ete <- function(tr) {
  if (ape::is.rooted(tr)) return(tr)
  root_node <- ape::Ntip(tr) + 1L
  ape::root(tr, node = root_node, resolve.root = TRUE)
}

#' @keywords internal
ensure_node_labels <- function(tr) {
  if (is.null(tr$node.label) || length(tr$node.label) != tr$Nnode) {
    tr$node.label <- paste0("N", seq_len(tr$Nnode))
  }
  tr
}

#' @keywords internal
label_to_node_number <- function(tr) {
  labs <- c(tr$tip.label, tr$node.label)
  structure(seq_along(labs), names = labs)
}

# children labels of a node label (in current tree/order)
#' @keywords internal
children_labels <- function(tr, parent_label) {
  lab2num <- label_to_node_number(tr)
  pnum <- unname(lab2num[parent_label])
  if (is.na(pnum)) return(character(0))
  rows <- which(tr$edge[,1] == pnum)
  if (!length(rows)) return(character(0))
  cld <- tr$edge[rows, 2]
  labs <- c(tr$tip.label, tr$node.label)
  labs[cld]
}

# ---- the dynamic-preorder branch-length IQR pruner ----
#' @keywords internal
prune_by_branch_length <- function(tree, threshold) {
  stopifnot(inherits(tree, "phylo"))
  if (is.null(tree$edge.length)) stop("phylo has no edge.length; cannot use branch-length thresholding.")

  get_desc_tip_labels <- function(tr, node_id) {
    ntip <- ape::Ntip(tr)
    children_map <- split(tr$edge[, 2], tr$edge[, 1])

    out <- character(0)
    stack <- node_id
    while (length(stack) > 0) {
      x <- stack[[1]]
      stack <- stack[-1]
      if (x <= ntip) {
        out <- c(out, tr$tip.label[x])
      } else {
        kids <- children_map[[as.character(x)]]
        if (!is.null(kids) && length(kids) > 0) stack <- c(kids, stack)
      }
    }
    unique(out)
  }


  levelorder_node_ids <- function(tr) {
    tr <- ape::reorder.phylo(tr, order = "cladewise")
    d <- ape::node.depth(tr)
    parent_depth <- d[tr$edge[, 1]]
    child_depth  <- d[tr$edge[, 2]]
    tr$edge[order(parent_depth, child_depth, seq_len(nrow(tr$edge))), 2]
  }


  upper_fence <- calculate_iqr_threshold(tree$edge.length)
  total_tips  <- ape::Ntip(tree)
  percent_left <- 100 - threshold


  tree <- ape::reorder.phylo(tree, order = "cladewise")

  repeat {
    pruned_this_pass <- FALSE


    node_order <- levelorder_node_ids(tree)

    for (node_id in node_order) {
      ntip <- ape::Ntip(tree)
      max_node_id <- ntip + tree$Nnode
      if (node_id < 1 || node_id > max_node_id) next

      all_leaves <- tree$tip.label
      num_leaves <- length(all_leaves)

      leaves_on_one_side <- get_desc_tip_labels(tree, node_id)
      one_side <- length(leaves_on_one_side)
      other_side <- num_leaves - one_side
      min_tips <- min(one_side, other_side)

      edge_idx <- which(tree$edge[, 2] == node_id)
      if (length(edge_idx) == 0) next
      node_dist <- tree$edge.length[edge_idx[1]]

      if (node_dist > upper_fence && min_tips < (percent_left * total_tips) / 100) {


        if (one_side > other_side) {
          tips_to_remove <- setdiff(all_leaves, leaves_on_one_side)
          tree <- ape::drop.tip(tree, tips_to_remove)
          percent_left <- percent_left - 100 * other_side / total_tips
        } else {
          tips_to_remove <- leaves_on_one_side
          tree <- ape::drop.tip(tree, tips_to_remove)
          percent_left <- percent_left - 100 * one_side / total_tips
        }

        tree <- ape::reorder.phylo(tree, order = "cladewise")
        pruned_this_pass <- TRUE
        break
      }
    }

    if (!pruned_this_pass) break
  }

  list(tree = tree, percent_left_budget = percent_left)
}





# ---------- Root-to-tip phase (unchanged; ape-only) ----------
#' @keywords internal
root_to_tip_ape <- function(tr) ape::node.depth.edgelength(tr)[1:ape::Ntip(tr)]


#' @keywords internal
midpoint_root <- function(tr) {
  return(phangorn::midpoint(tr))
}


#' @keywords internal
prune_by_root_to_tip <- function(tree, percent_left, original_num_leaves) {
  stopifnot(inherits(tree, "phylo"))

  cutoff <- (percent_left * original_num_leaves) / 100

  tr0 <- midpoint_root(tree)
  dist0 <- root_to_tip_ape(tr0)
  names(dist0) <- tr0$tip.label

  if (length(dist0) >= 4L) {
    fence0 <- calculate_iqr_threshold(unname(dist0))
    outliers0 <- which(dist0 > fence0)
    if (length(outliers0) > cutoff) return(tree)
  } else {
    return(tree)
  }

  removed_counter <- 0L
  while (removed_counter < cutoff) {
    tr0 <- midpoint_root(tree)
    dist_rt <- root_to_tip_ape(tr0)
    if (length(dist_rt) < 4L) break
    names(dist_rt) <- tr0$tip.label

    fence <- calculate_iqr_threshold(unname(dist_rt))


    extreme_name <- names(dist_rt)[which.max(dist_rt)]
    if (dist_rt[extreme_name] > fence) {
      tree <- ape::drop.tip(tree, extreme_name)
      removed_counter <- removed_counter + 1L
      percent_left <- percent_left - 100 / original_num_leaves
      cutoff <- (percent_left * original_num_leaves) / 100
    } else {
      break
    }
  }

  ape::unroot(tree)
}



validate_algorithm_args <- function(algorithm, args) {
  allowed_args <- list(
    CPA = c(
      "root_to_node_ratio",
      "min_num_of_roots",
      "M_n",
      "beta",
      "radius_ratio",
      "safe_tips",
      "show_plot",
      "show_pruned_tips"
    ),
    PSFA = c(
      "longest_to_average"
    ),
    IQR = character(0)
  )

  unknown_args <- setdiff(names(args), allowed_args[[algorithm]])

  if (length(unknown_args) > 0) {
    stop(
      "Unknown argument(s) for ", algorithm, ": ",
      paste(unknown_args, collapse = ", ")
    )
  }

  args
}
