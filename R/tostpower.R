# Some useful keyboard shortcuts for package authoring:
#' TOST power and sample size calculations.
#'
#' Calculate the power to detect equivalence or non-inferiority for either
#' paired or unpaired experiments.
#'
#' @param n The sample size: for \code{paired == TRUE}, the number differences;
#' for \code{paired == FALSE} (the default), the number in each group.
#' @param power The power to detect equivalence (or non-inferiority).
#' @param true_dif The true (ie, unknown) difference in population mean (as a
#'  proportion of the reference mean).
#' @param threshold The maximum difference between the two means that can be
#'  considered equivalent (as a proportion of the reference mean).
#' @param sd The standard devation (as a proportion of the reference mean).
#' @param sig_level The significance (i.e. alpha) level of the test.
#' @param paired Is this for paired comparisons?
#' @param type Is this designed to show equivalence or simply non-inferiority?
#'
#' @return Power to detect equivalence
#'
#' @author Dennis L. Malandro, \email{dennismalandro@@gmail.com}
#'
#' @references
#' Kem F. Phillips (1990) Power of the Two One-Sided Tests Procedure in
#' Bioequivalence, \emph{J. Pharmacokinetics and Biopharmaceutics} Vol 18, No. 2
#'
#' Donald J. Schuirmann (1987) A Comparison of the Two One-Sided Tests Procedure and the
#' Power Approach for Assessing the Equivalence of Average Bioavailability,
#'  \emph{J. Pharmacokinetics and Biopharmaceutics} Vol 15, No. 6
#'
#' @seealso \code{\link{power.t.test}}
#'
#' @examples
#' tost_power(9)
#'
#' tost_ss(paired = TRUE, type = 'non_inferiority')
#'
#' @keywords hypothesis-test power sample-size
#'
#' @export
tost_power <- function(n, true_dif = 0, threshold = 0.2, sd = 0.2,
  sig_level = 0.05, paired = FALSE, type = c('equivalence', 'non_inferiority')) {

  type <- match.arg(type)

  m <- 2 - paired # 1 for paired, 2 for !paired

  # for paired, n and se correspond to *differences*
  # for !paired, n is number in *each* group, 2n is total ss
  std_error <- sd * sqrt(m / n) # weird way to do it
  nu <- (n - 1) * m

  qu <- qt(sig_level, nu, lower.tail = FALSE)

  if(length(threshold) == 1) threshold <- c(-abs(threshold), abs(threshold))
  ncp <- (true_dif - threshold) / std_error

  if (type == 'equivalence') {
#    pmvt(lower = c(qu, -Inf), upper = c(Inf, -qu), df = floor(nu), delta = ncp)
    p <- pmvt(lower = c(qu, -Inf), upper = c(Inf, -qu), df = nu, delta = ncp)
    attributes(p) <- NULL
  } else {
    p <- pt(qu, df = nu, ncp = ncp[1], lower.tail = FALSE)
  }
#   list(power = p, n = n, true_dif = true_dif,
#     threshold = paste0('(', threshold[1], ', ', threshold[2], ')'), sd = sd,
#     sig_level = sig_level, paired = paired, type = type)
  p

}

#' @export
#'
#' @rdname tost_power
#'
tost_ss <- function(power = 0.8, true_dif = 0, threshold = 0.2, sd = 0.2,
  sig_level = 0.05, paired = FALSE, type = c('equivalence', 'non_inferiority')) {

  type <- match.arg(type)
  if(length(threshold) == 1) threshold <- c(-abs(threshold), abs(threshold))

  u <- uniroot(
    function(n) tost_power(n, true_dif = true_dif, threshold = threshold,
      sd = sd, sig_level = sig_level, paired = paired, type = type) - power,
    c(2, 999), extendInt = 'no')
  ceiling(u$root)
}

