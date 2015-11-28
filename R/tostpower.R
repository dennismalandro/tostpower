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
#' @return Power to detect equivalence/non-inferiority
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
tost_power <- function(n, true_dif = 0, threshold = 0.2,
  sd = 0.2, sig_level = 0.05, paired = FALSE,
  type = c('equivalence', 'non_inferiority')) {

  type <- match.arg(type)

  m <- 2 - paired  # m = 1 for paired data, 2 for unpaired
  std_error <- sd * sqrt(m / n)
  nu <- (n - 1) * m
  # for non-paired data, n is number in *each* group

  qu <- qt(sig_level, nu, lower.tail = FALSE)

  if(length(threshold) == 1)
    threshold <- abs(threshold) * c(-1, 1)

  ncp <- (true_dif - threshold) / std_error

  if (type == 'equivalence') {
    p <- mvtnorm::pmvt(
      lower = c(qu, -Inf),
      upper = c(Inf, -qu),
      df = nu,
      delta = ncp  # ncp is a length-two vector
    )
    attributes(p) <- NULL
  } else {# univariate t-dist for non-inferiority
    p <- pt(qu, df = nu, ncp = ncp[1], lower.tail = FALSE)
  }
  p

}
