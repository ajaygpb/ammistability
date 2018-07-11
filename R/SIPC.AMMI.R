#' Sums of the Absolute Value of the IPC Scores
#'
#' \code{SIPC.AMMI} computes the Sums of the Absolute Value of the IPC Scores
#' (ASI) (Sneller et al., 1997) considering all significant interaction
#' principal components (IPCs) in the AMMI model. Using SIPC, the Yield
#' stability Index (YSI) is also calculated.
#'
#' The Sums of the Absolute Value of the IPC Scores (\eqn{SIPC}) is computed as
#' follows:
#'
#' \deqn{SIPC = \sum_{n=1}^{N'} \left | \lambda_{n}^{0.5}\gamma_{in} \right |}
#'
#' OR
#'
#' \deqn{SIPC = \sum_{n=1}^{N'}\left | PC_{n} \right |}
#'
#'
#' Where, \eqn{N'} is the number of significant IPCs (number of IPC that were
#' retained in the AMMI model via F tests); \eqn{\lambda_{n}} is the is the
#' singular value for IPC \eqn{n} and correspondingly \eqn{\lambda_{n}^{2}}  is
#' its eigen value; \eqn{\gamma_{in}} is the eigenvector value for \eqn{i}th
#' genotype; and \eqn{PC_{1}}, \eqn{PC_{2}}, \eqn{\cdots}, \eqn{PC_{n}} are the
#' scores of 1th, 2th, ..., and \eqn{n}th IPC.
#'
#' The closer the SIPC scores are to zero, the more stable the genotypes are
#' across test environments.
#'
#' The Yield Stability Index (\eqn{YSI}) is computed as follows:
#'
#' \deqn{YSI = R_{SIPC} + R_{Y}}
#'
#' Where, \eqn{R_{SIPC}} is the SIPC rank of the genotype and \eqn{R_{Y}} is the
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
SIPC.AMMI <- function(model, n, alpha = 0.05) {

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

# cova <- cov(model$genXenv)
# values <- eigen(cova)
# SIPC1 <- (sqrt(values$values[1])*model$biplot[,3])
# SIPC2 <- (sqrt(values$values[2])*model$biplot[,4])
# SIPC3 <- (sqrt(values$values[3])*model$biplot[,5])
# SIPC <- SIPC1+SIPC2+SIPC3
# rs <- rank(SIPC)
# rSIPC <- data.frame(SIPC1,SIPC2,SIPC3,SIPC,rs)
# rSIPC

  # # GxE matrix
  # ge <- array(model$genXenv, dim(model$genXenv), dimnames(model$genXenv))
  # # SVD
  # svdge <- svd(ge)
  # lambda.n <- svdge$d[1:n]
  # gamma.n <- svdge$u[,1:n]
  # A <- sqrt(lambda.n)*gamma.n

  A <- model$biplot
  A <- A[A[,1] == "GEN",-c(1,2)]
  A <- A[,1:n] # Fetch only n IPCs

  SIPC <- unname(rowSums(apply(A, 2, FUN = abs)))

  rk <- rank(SIPC)
  B <- model$means
  W <- tapply.stat(B[,3],B[,2],function(x) mean(x,rm.na = TRUE))
  Rx <- rank(-W[,2])
  YSI_SIPC <- rk + Rx
  ranking <- data.frame(SIPC, YSI_SIPC, rSIPC = rk, rY = Rx, means = W[,2])

  return(ranking)

}
