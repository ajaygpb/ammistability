#' Stability Measure Based on Fitted AMMI Model
#'
#' \code{FA.AMMI} computes the Stability Measure Based on Fitted AMMI Model (FA)
#' (Raju, 2002) considering all significant interaction principal components
#' (IPCs) in the AMMI model. Using FA, the Simultaneous Selection Index for
#' Yield and Stability (SSI) is also calculated according to the argument
#' \code{ssi.method}.
#'
#' The Absolute value of the Relative Contribution of IPCs to the Interaction
#' (\eqn{D_{Za}}) is computed as follows:
#'
#' \deqn{FA = \sum_{n=1}^{N'}\lambda_{n}^{2}\gamma_{in}^{2}}
#'
#' Where, \eqn{N'} is the number of significant IPCs (number of IPC that were
#' retained in the AMMI model via F tests); \eqn{\lambda_{n}} is the is the
#' singular value for IPC \eqn{n} and correspondingly \eqn{\lambda_{n}^{2}}  is
#' its eigen value; and \eqn{\gamma_{in}} is the eigenvector value for \eqn{i}th
#' genotype.
#'
#' When \eqn{N'} is replaced by 1 (only first IPC axis is considered for
#' computation), then the parameter \eqn{FP} can be estimated.
#'
#' \deqn{FP = \lambda_{1}^{2}\gamma_{i1}^{2}}
#'
#' @inheritParams MASV.AMMI
#'
#' @return
#'
#' @importFrom agricolae AMMI
#' @export
#'
#' @references
#'
#' \insertRef{raju_study_2002}{AMMIStbP}
#'
#' \insertRef{zali_evaluation_2012}{AMMIStbP}
#'
#' @seealso \code{\link[AMMIStbP]{SSI}}
#'
#' @examples
FA.AMMI <- function(model, n, alpha = 0.05,
                    ssi.method = c("farshadfar", "rao"), a = 1) {

  # Check model class
  if (!is(model, "AMMI")) {
    stop('"model" is not of class "AMMI"')
  }

  # Check alpha value
  if (!(0 < alpha && alpha < 1)) {
    stop('"alpha" should be between 0 and 1 (0 < alpha < 1)')
  }

  # Find number of significant IPCs according to F test
  if (missing(n) || is.null(n)) {
    n = sum(model$analysis$Pr.F <= alpha, na.rm = TRUE)
  }

  # Check for n
  if (n %% 1 != 0 && length(n) != 1) {
    stop('"n" is not an integer vector of unit length')
  }

  ssi.method <- match.arg(ssi.method)



# cova<-cov(model$genXenv)
# values<-eigen(cova)
# D1<-sqrt((values$values[1]*model$biplot[,3])^2)
# D2<-sqrt((values$values[2]*model$biplot[,4])^2)
# D3<-sqrt((values$values[3]*model$biplot[,5])^2)
# D<-D1+D2+D3
# FA<-D1^2+D2^2+D3^2
# rf<-rank(FA)
# rFA<-data.frame(FA,rf)
# rFA

  # GxE matrix
  ge <- array(model$genXenv, dim(model$genXenv), dimnames(model$genXenv))
  # SVD
  svdge <- svd(ge)
  gamma.n <- svdge$u[,1:n]
  lambda.n <- svdge$d[1:n]

  FA <- rowSums((lambda.n^2) * (gamma.n^2))

  B <- model$means
  W <- aggregate(B$Yield, by = list(model$means$GEN), FUN = mean, na.rm = TRUE)
  SSI_FA <- SSI(y = W$x, sp = FA, gen = W$Group.1,
                method = ssi.method, a = a)
  ranking <- SSI_FA
  colnames(ranking) <- c("FA", "SSI", "rFA", "rY", "means")

  return(ranking)

}
