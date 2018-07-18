#' Zhang's D Parameter
#'
#' \code{DZ.AMMI} computes the Zhang's D Parameter values or AMMI statistic
#' coefficient or AMMI distance or AMMI stability index
#' (\eqn{\textrm{D}_{\textrm{z}}}) (Zhang, 1998) considering all significant
#' interaction principal components (IPCs) in the AMMI model. It is the distance
#' of IPC  point  from  origin  in space. Using \eqn{\textrm{D}_{\textrm{z}}},
#' the Simultaneous Selection Index for Yield and Stability (SSI) is also
#' calculated according to the argument \code{ssi.method}.
#'
#' The Zhang's D Parameter value (\eqn{D_{z}}) is computed as follows:
#'
#' \deqn{D_{z} = \sqrt{\sum_{n=1}^{N'}\gamma_{in}^{2}}}
#'
#' Where, \eqn{N'} is the number of significant IPCAs (number of IPC that were
#' retained in the AMMI model via F tests); and \eqn{\gamma_{in}} is the
#' eigenvector value for \eqn{i}th genotype.
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
#' \insertRef{zhang_analysis_1998}{AMMIStbP}
#'
#' @seealso \code{\link[AMMIStbP]{SSI}}
#'
#' @examples
DZ.AMMI <- function(model, n, alpha = 0.05,
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

  # GxE matrix
  ge <- array(model$genXenv, dim(model$genXenv), dimnames(model$genXenv))
  # SVD
  svdge <- svd(ge)
  gamma.n <- svdge$u[,1:n]

  DZ <- sqrt(rowSums((gamma.n)^2))

  B <- model$means
  W <- aggregate(B$Yield, by = list(model$means$GEN), FUN = mean, na.rm = TRUE)
  SSI_DZ <- SSI(y = W$x, sp = DZ, gen = W$Group.1,
                method = ssi.method, a = a)
  ranking <- SSI_DZ
  colnames(ranking) <- c("DZ", "SSI", "rDZ", "rY", "means")

  return(ranking)

}
