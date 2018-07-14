#' AMMI Based Stability Parameter
#'
#' \code{ASTAB.AMMI} computes the AMMI Based Stability Parameter (ASTAB) (Rao
#' and Prabhakaran, 2005) considering all significant interaction principal
#' components (IPCs) in the AMMI model. Using ASTAB, the Yield stability Index
#' (YSI) is also calculated.
#'
#' The AMMI Based Stability Parameter value (\eqn{ASTAB}) is computed as
#' follows:
#'
#' \deqn{ASTAB = \sum_{n=1}^{N'}\lambda_{n}\gamma_{in}^{2}}
#'
#' Where, \eqn{N'} is the number of significant IPCAs (number of IPC that were
#' retained in the AMMI model via F tests); \eqn{\lambda_{n}} is the is the
#' singular value for IPC \eqn{n} and correspondingly \eqn{\lambda_{n}^{2}}  is
#' its eigen value; and \eqn{\gamma_{in}} is the eigenvector value for \eqn{i}th
#' genotype.
#'
#' The Yield Stability Index (\eqn{YSI}) is computed as follows:
#'
#' \deqn{YSI = R_{ASTAB} + R_{Y}}
#'
#' Where, \eqn{R_{ASTAB}} is the ASTAB rank of the
#' genotype and \eqn{R_{Y}} is the mean yield rank of the genotype.
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
#' \insertRef{rao_use_2005}{AMMIStbP}
#'
#' @examples
ASTAB.AMMI <- function(model, n, alpha = 0.05) {

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
  gamma.n <- svdge$u[,1:n]
  lambda.n <- svdge$d[1:n]

  ASTAB <- rowSums(lambda.n*((gamma.n)^2))

  rk <- rank(ASTAB)
  B <- model$means
  W <- tapply.stat(B[,3],B[,2],function(x) mean(x,rm.na = TRUE))
  Rx <- rank(-W[,2])
  YSI_ASTAB <- rk + Rx
  ranking <- data.frame(ASTAB, YSI_ASTAB, rASTAB = rk, rY = Rx, means = W[,2])

  return(ranking)

}
