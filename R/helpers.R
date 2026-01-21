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
  # 0) set up original reference & fence
  tree0 <- ensure_node_labels(root_like_ete(tree))
  upper_fence0 <- calculate_iqr_threshold(tree0$edge.length)
  total_tips0  <- ape::Ntip(tree0)
  percent_left <- 100 - threshold

  # 1) work on a live copy (keep node labels so we can refind nodes)
  tr <- tree0
  tr <- ensure_node_labels(tr)

  # define root label once, based on current tr
  root_label <- c(tr$tip.label, tr$node.label)[ape::Ntip(tr) + 1L]

  # Queue for levelorder (breadth-first)
  queue <- list(root_label)

  dequeue <- function() {
    x <- queue[[1]]
    queue <<- queue[-1]
    x
  }

  enqueue <- function(labels) {
    if (length(labels)) queue <<- c(queue, as.list(labels))
  }

  while (length(queue)) {
    node_lbl <- dequeue()
    lab2num <- label_to_node_number(tr)
    node_num <- unname(lab2num[node_lbl])
    if (is.na(node_num)) next  # may have been pruned

    # skip root
    if (node_num != (ape::Ntip(tr) + 1L)) {
      hit <- which(tr$edge[,2] == node_num)
      if (length(hit)) {
        node_dist <- tr$edge.length[hit[1L]]
        if (!is.na(node_dist) && node_dist > upper_fence0) {
          tip_ids <- 1:ape::Ntip(tr)
          leaves_on_one_side <- phangorn::Descendants(tr, node_num, "tips")[[1]]
          leaves_on_other_side <- setdiff(tip_ids, leaves_on_one_side)
          one_side <- length(leaves_on_one_side)
          other_side <- length(leaves_on_other_side)
          min_tips <- min(one_side, other_side)

          if (min_tips < (percent_left * total_tips0) / 100) {
            to_remove <- if (one_side <= other_side) leaves_on_one_side else leaves_on_other_side
            if (length(to_remove) > 0) {
              tr <- ape::drop.tip(tr, to_remove, collapse.singles = TRUE)
              removed <- length(to_remove)
              percent_left <- percent_left - 100 * removed / total_tips0
              tr <- ensure_node_labels(tr)
            }
          }
        }
      }
    }

    # enqueue children after processing current node
    ch_labs <- children_labels(tr, node_lbl)
    enqueue(ch_labs)
  }

  list(tree = tr, percent_left_budget = percent_left)
}

# ---------- Root-to-tip phase (unchanged; ape-only) ----------
#' @keywords internal
root_to_tip_ape <- function(tr) ape::node.depth.edgelength(tr)[1:ape::Ntip(tr)]


#' @keywords internal
midpoint_root <- function(tr) {
  return(phytools::midpoint.root(tr))
}

#' @keywords internal
prune_by_root_to_tip <- function(tree, percent_left, original_num_leaves) {
  tr0 <- midpoint_root(tree)
  # Identify the two child edges of the new root
  root_edges <- which(tr0$edge[,1] == ape::Ntip(tr0) + 1L)
  root_mean <- mean(tr0$edge.length[root_edges])

  tr0$edge.length[root_edges] <- rep(root_mean[1], length(root_edges))


  dist0 <- root_to_tip_ape(tr0)
  names(dist0) <- tree$tip.label
  out0 <- if (length(dist0) >= 4L) which(dist0 > calculate_iqr_threshold(dist0)) else integer(0)
  #out0 <- names(dist0)[dist0 > calculate_iqr_threshold(dist0)]
  if (length(out0) > (percent_left * original_num_leaves) / 100) {
    return(tree)  # early exit, like Python
  }

  removed_counter <- 0L
  while (removed_counter < (percent_left * original_num_leaves) / 100) {
    tr0 <- midpoint_root(tree)
    root_edges <- which(tr0$edge[,1] == ape::Ntip(tr0) + 1L)
    root_mean <- mean(tr0$edge.length[root_edges])

    tr0$edge.length[root_edges] <- rep(root_mean[1], length(root_edges))
    dist_rt <- root_to_tip_ape(tr0)
    if (length(dist_rt) < 4L) break
    fence <- calculate_iqr_threshold(dist_rt)
    out <- which(dist_rt > fence)
    if (length(out) == 0) break
    extreme <- out[ which.max(dist_rt[out]) ]
    tree <- ape::drop.tip(tree, tree$tip.label[extreme])
    removed_counter <- removed_counter + 1L
    percent_left <- percent_left - 100 / original_num_leaves
  }
  ape::unroot(tree)
}
