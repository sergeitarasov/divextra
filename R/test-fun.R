#' Report p-value
#'
#' Reports a p-value in APA style.
#'
#' @param p The p-value
#' @param digits The number of digits to round to (default = 3)
#'
#' @return A string with the format "p = .040" or "p < .001"
#' @examples
#'
#' report_p(0.02018) # returns "p = .020"
#' report_p(0.00028) # returns "p < .001"
#'
#' @export
report_p <- function(p, digits = 3) {
  if (p < .001) return("p < .001")

  p_round <- round(p, digits) %>%
    as.character() %>%
    # omit leading zero for APA-style
    stringr::str_replace("0.", ".") %>%
    # pad right with zeros
    stringr::str_pad(digits+1, "right", 0)

  p_string <- paste("p =", p_round)

  return(p_string)
}
