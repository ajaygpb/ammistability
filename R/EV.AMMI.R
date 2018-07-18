#' Averages of the Squared Eigenvector Values
#'
#' \code{EV.AMMI} computes the Sums of the Averages of the Squared Eigenvector
#' Values (EV) (Zobel, 1994) considering all significant interaction principal
#' components (IPCs) in the AMMI model. Using EV, the Simultaneous Selection
#' Index for Yield and Stability (SSI) is also calculated according to the
#' argument \code{ssi.method}.
#'
#' The Averages of the Squared Eigenvector Values (\eqn{EV}) is computed as
#' follows:
#'
#' \deqn{EV = \sum_{n=1}^{N'}\frac{\gamma_{in}^2}{N'}}
#'
#' Where, \eqn{N'} is the number of significant IPCs (number of IPC that were
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
#' \insertRef{zobel_stress_1994}{AMMIStbP}
#'
#' @seealso \code{\link[AMMIStbP]{SSI}}
#'
#' @examples
EV.AMMI <- function(model, n, alpha = 0.05,
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
  gamma.n <- svdge$u[,1:n]

  EV <- rowSums(gamma.n^2/n)

  B <- model$means
  W <- aggregate(B$Yield, by = list(model$means$GEN), FUN = mean, na.rm = TRUE)
  SSI_EV <- SSI(y = W$x, sp = EV, gen = W$Group.1,
                method = ssi.method, a = a)
  ranking <- SSI_EV
  colnames(ranking) <- c("EV", "SSI", "rEV", "rY", "means")

  return(ranking)

}
