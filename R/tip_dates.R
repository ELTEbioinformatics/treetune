#' Get valid and ordered tip dates for rooting
#'
#' This function takes a phylogenetic tree and a data frame containing tip
#' information, and returns a vector of tip dates ordered according to the
#' tree's tip order. Character strings representing ambiguous dates are
#' coerced to Date format using the specified strategy.
#'
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param df A data frame containing tip information.
#' @param label_var Name of the column in \code{df} containing tip labels.
#' @param date_var Name of the column in \code{df} containing dates.
#' @param strategy Strategy for handling ambiguous date strings: "lower" 
#' (earliest possible date), "middle" (midpoint), "upper" (latest possible 
#' date), or "random" (random date within range). Default is "middle". See 
#' \code{date_lower}, \code{date_middle}, \code{date_upper}, and 
#' \code{date_runif} for details.
#' @return A vector of Date objects ordered to match \code{tree$tip.label}.
#' @export
tip_dates <- function(tree, df, label_var, date_var, strategy = "middle") {
  assert(tree, "phylo")
  assert(df, "data.frame")
  strategy <- match.arg(
    strategy, choices = c("lower", "middle", "upper", "random")
  )
  if (!label_var %in% names(df)) {
    stop(paste0("'df' must contain a column called '", label_var, "'."))
  }
  assert(df[[label_var]], "character")
  # Check that all tree tips are present in df
  missing_tips <- setdiff(tree$tip.label, df[[label_var]])
  if (length(missing_tips) > 0) {
    stop("The following tips are missing from 'df': ", 
         paste(missing_tips, collapse = ", "))
  }
  # Filter df to only include rows matching tree tips
  df <- df[which(df[[label_var]] %in% tree$tip.label), ]
  # Reorder df rows to match tree$tip.label order
  df <- df[match(tree$tip.label, df[[label_var]]), ]
  # Check dates
  if (!date_var %in% names(df)) {
    stop(paste0("'df' must contain a column called '", date_var, "'."))
  }
  tip_dates <- df[[date_var]]
  # Fix dates, if necessary
  if (inherits(tip_dates, "Date")) {
    # already a Date vector; no conversion needed
  } else if (inherits(tip_dates, c("POSIXct", "POSIXt"))) {
    tip_dates <- as.Date(tip_dates)
  } else if (inherits(tip_dates, "character")) {
    if (strategy == "lower") {
      tip_dates <- date_lower(tip_dates)
    } else if (strategy == "middle") {
      tip_dates <- date_middle(tip_dates)
    } else if (strategy == "upper") {
      tip_dates <- date_upper(tip_dates)
    } else if (strategy == "random") {
      tip_dates <- date_runif(tip_dates)
    }
  } else {
    stop("'date_var' must be of class Date, POSIXct, POSIXt, or character.")
  }
  return(tip_dates)
}