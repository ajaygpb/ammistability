#' AMMI Stability Index
#'
#' \code{ASI.AMMI} computes the AMMI Stability Index (ASI) (Jambhulkar et al.,
#' 2014; Jambhulkar et al., 2015; Jambhulkar et al., 2017) considering the first
#' two interaction principal components in the AMMI model. Using ASI, the Yield
#' stability Index (YSI) is also calculated.
#'
#' The AMMI Stability Index (\eqn{ASI}) is computed as follows:
#'
#' \deqn{ASI = \sqrt{\left [ PC_{1}^{2} \times \theta_{1}^{2} \right ]+\left [
#' PC_{2}^{2} \times \theta_{2}^{2} \right ]}}
#'
#' Where, \eqn{PC_{1}} and \eqn{PC_{2}} are the scores of 1st and 2nd
#' interaction principal components (IPCs) respectively; and \eqn{\theta_{1}} and
#' \eqn{\theta_{2}} are percentage sum of squares explained by the 1st and 2nd
#' principal component interaction effect respectively.
#'
#' The Yield Stability Index (\eqn{YSI}) is computed as follows:
#'
#' \deqn{YSI = R_{ASI} + R_{Y}}
#'
#' Where, \eqn{R_{ASI}} is the ASI rank of the genotype and \deqn{R_{Y}} is
#' the mean yield rank of the genotype.
#'
#' @inheritParams MASV.AMMI
#'
#' @return
#'
#' @importFrom agricolae tapply.stat
#' @importFrom agricolae AMMI
#' @export
#'
#' @references
#'
#' \insertRef{jambhulkar_ammi_2014}{AMMIStbP}
#'
#' \insertRef{jambhulkar_genotype_2015}{AMMIStbP}
#'
#' \insertRef{jambhulkar_stability_2017}{AMMIStbP}
#'
#' @examples
ASI.AMMI <- function(model, n, alpha = 0.05) {

  # Check model class
  if (!is(model, "AMMI")) {
    stop('"model" is not of class "AMMI"')
  }

  # Check alpha value
  if (!(0 < alpha && alpha < 1)) {
    stop('"alpha" should be between 0 and 1 (0 < alpha < 1)')
  }

  # Find number of significant IPCs according to F test
  if (missing(n) | is.null(n)) {
    n = sum(model$analysis$Pr.F <= alpha, na.rm = TRUE)
  }

  # Check for n
  if (n %% 1 != 0 && length(n) != 1) {
    stop('"n" is not an integer vector of unit length')
  }

  A <- model$biplot[, 1:4]
  A <- A[A[, 1] == "GEN", -c(1, 2)]

  th1 <- model$analysis["PC1",]$percent/100
  th2 <- model$analysis["PC2",]$percent/100

  ASI <- sqrt(((A[,"PC1"]^2) + (th1^2)) + ((A[,"PC2"]^2) + (th2^2)))

  rk <- rank(ASI)
  B <- model$means
  W <- tapply.stat(B[,3],B[,2],function(x) mean(x,rm.na = TRUE))
  Rx <- rank(-W[,2])
  YSI_ASI <- rk + Rx
  ranking <- data.frame(ASI, YSI_ASI, rASI = rk, rY = Rx, means = W[,2])

  return(ranking)

}
