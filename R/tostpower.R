#' TOST power and sample size calculations.
#'
#' Calculate the power to detect equivalence or non-inferiority for either
#' paired or unpaired experiments.
#'
#' @param n The number in each group for \code{!paired} (the default);
#'  the number of pairs for \code{paired}
#' @param sigma The standard devation: total for \code{!paired}, within for \code{paired}.
#' @param true_diff The true (unknown) difference of population means.
#' @param eqv_interval The interval within which \code{true_diff} is _equivalant_
#'  considered equivalent (as a proportion of the reference mean).
#' @param alpha The significance level of the test.
#' @param paired Is this a paired comparison?
#' @param df Degress of freedom, (usually not set independent of \code{n}, but can be, see details).
#'
#' @return Power for both equivalence and non-inferiority alternatives
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
#' # Reproduce of n=9 part of power curve in Phillips, Fig. 1 on p. 139
#'
#' delta_mu <- c(0, 5, 10, 15)
#'
#' sapply(delta_mu, function(x) {
#'   tost_power(true_diff = x,
#'     n = 9, sigma = 10,
#'     eqv_interval = c(-20, 20), df = 7)[['eq']]
#'   }
#' )
#'
#' @keywords hypothesis-test power sample-size
#'
#' @export
tost_power <- function(
  n = 10, sigma = 1,
  true_diff = 0,
  eqv_interval = c(-1, 1),
  alpha = 0.05,
  paired = FALSE,
  df = NULL) {

  # set df ind of n if needed cuzza other params (eg x-over)
  if(is.null(df)) df <- ifelse(paired, n - 1, 2 * n - 2)


  # rejection region under H0: univariate central t twice
  rej_region <- qt(c(1 - alpha, alpha), df = df)


  # bivariate noncentrality parameter
  ncp <- (true_diff - eqv_interval) / sigma * sqrt(n / 2)


  if(true_diff < eqv_interval[1]) {
    stop("mean difference can't be below equivalence region")

  } else {

    # equivalence alternative is bivariate noncentral t
    eqv_power <- mvtnorm::pmvt(
      lower = c(rej_region[1], -Inf),
      upper = c(Inf, rej_region[2]),
      df = df,
      delta = ncp,
      algorithm = mvtnorm::GenzBretz(abseps = 0.00001))
    attributes(eqv_power) <- NULL

    # noninf alt is univariate noncentral
    noninf_power <- pt(rej_region[1], df, ncp = ncp[1], lower = F)}

  # power for both equivalence and noninferiority alternatives
  c(eq = eqv_power, noninf = noninf_power)}
