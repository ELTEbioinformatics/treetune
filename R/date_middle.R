#' Convert date to numeric format for phylogenetic dating
#'
#' This is a utility function that converts a date to a numeric value. This
#' format is commonly used in phylogenetic dating. If date is incomplete the
#' function returns the middle of the interval.
#' @param dates character; a full or partial date in "YYYY-MM_DD" like format.
#' @param out_format character; output format, either \code{"decimal"} or
#' \code{"date"}.
#' @examples
#' date_middle("2021-18-11")
#' date_middle("1988-03")
#' date_middle("2000")
#' @export
date_middle <- function(dates, out_format = "date"){
  
  out_format <- match.arg(out_format, choices = c("decimal", "date"))

  dates <- as.character(dates)
  
  foo <- function(x){
    if(is.na(x)) return(NA)
    lower <- date_lower(x, out_format = "decimal")
    upper <- date_upper(x, out_format = "decimal")
    middle <- mean(lower)
    return(middle)
  }
  out <- sapply(dates, foo)
  
  if (out_format == "date") {
    out <- lubridate::date_decimal(out) |> unname() |> as.Date()
  }
  
  return(out)
}
