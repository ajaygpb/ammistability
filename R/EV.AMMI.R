#' Averages of the Squared Eigenvector Values
#'
#' \code{EV.AMMI} computes the Sums of the Averages of the Squared Eigenvector
#' Values (EV) (Zobel, 1994) all significant interaction principal components in
#' the AMMI model. Using EV, the Yield stability Index (YSI) is also calculated.
#'
#' The Averages of the Squared Eigenvector Values (\eqn{EV}) is computed as
#' follows:
#'
#' \deqn{EV = \sum_{n=1}^{N'}\frac{\gamma_{in}^2}{N'}}
#'
#' Where, \eqn{N'} is the number of significant IPCAs (number of IPC that were
#' retained in the AMMI model via F tests); and \eqn{\gamma_{in}} is the
#' eigenvector value for \eqn{i}th genotype.
#'
#' The Yield Stability Index (\eqn{YSI}) is computed as follows:
#'
#' \deqn{YSI = R_{EV} + R_{Y}}
#'
#' Where, \eqn{R_{EV}} is the EV rank of the genotype and \eqn{R_{Y}} is the
#' mean yield rank of the genotype.
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
#' \insertRef{zobel_stress_1994}{AMMIStbP}
#'
#' @examples
EV.AMMI <- function(model, n, alpha = 0.05) {

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


# EV1<-model$biplot[,3]^2
# EV2<-model$biplot[,4]^2
# EV3<-model$biplot[,5]^2
# EV<-EV1+EV2+EV3
# rEV<-rank(EV)
# rEV<-data.frame(EV1,EV2,EV3,EV,rEV)
# rEV

  # GxE matrix
  ge <- array(model$genXenv, dim(model$genXenv), dimnames(model$genXenv))
  # SVD
  svdge <- svd(ge)
  gamma <- svdge$u[,1:n]

  EV <- rowSums(gamma^2/n)

  rk <- rank(EV)
  B <- model$means
  W <- tapply.stat(B[,3],B[,2],function(x) mean(x,rm.na = TRUE))
  Rx <- rank(-W[,2])
  YSI_EV <- rk + Rx
  ranking <- data.frame(MASV, YSI_EV, rEV = rk, rY = Rx, means = W[,2])

  return(ranking)

}
