#' Sum Across Environments of GEI Modelled by AMMI
#'
#' \code{AMGE.AMMI} computes the Sum Across Environments of GEI Modelled by AMMI
#' (AMGE) (Sneller et al., 1997) considering all significant interaction
#' principal components in the AMMI model. Using AMGE, the Yield stability Index
#' (YSI) is also calculated.
#'
#' The Sum Across Environments of GEI Modelled by AMMI (\eqn{AMGE}) is computed
#' as follows:
#'
#' \deqn{AMGE = \sum_{j=1}^{E} \sum_{n=1}^{N'} \lambda_{n} \gamma_{in}
#' \delta_{jn}}
#'
#'
#' Where, \eqn{N'} is the number of significant IPCAs (number of IPC that were
#' retained in the AMMI model via F tests); \eqn{\lambda_{n}} is the is the
#' singular value for IPC \eqn{n} and correspondingly \eqn{\lambda_{n}^{2}}  is
#' its eigen value; \eqn{\gamma_{in}} is the eigenvector value for \eqn{i}th
#' genotype; and \eqn{\delta_{jn}}  is the eigenvector value for \eqn{j}th
#' environment.
#'
#' The Yield Stability Index (\eqn{YSI}) is computed as follows:
#'
#' \deqn{YSI = R_{AMGE} + R_{Y}}
#'
#' Where, \eqn{R_{AMGE}} is the SIPC rank of the genotype and \eqn{R_{Y}} is the
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
#' \insertRef{sneller_repeatability_1997}{AMMIStbP}
#'
#' @examples
AMGE.AMMI <- function(model, n, alpha = 0.05) {

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

  # GxE matrix
  ge <- array(model$genXenv, dim(model$genXenv), dimnames(model$genXenv))
  # SVD
  svdge <- svd(ge)
  lambda.n <- svdge$d[1:n]
  gamma.n <- svdge$u[,1:n]
  delta.n <- svdge$v[,1:n]

  ge.n <- gamma.n %*% diag(lambda.n) %*% t(delta.n)

  AMGE <- rowSums(ge.n)

  rk <- rank(AMGE)
  B <- model$means
  W <- tapply.stat(B[,3],B[,2],function(x) mean(x,rm.na = TRUE))
  Rx <- rank(-W[,2])
  YSI_AMGE <- rk + Rx
  ranking <- data.frame(AMGE, YSI_AMGE, rAMGE = rk, rY = Rx, means = W[,2])

  return(ranking)

}
