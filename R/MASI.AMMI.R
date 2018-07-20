#' Modified AMMI Stability Index
#'
#' \code{MASI.AMMI} computes the Modified AMMI Stability Index (MASI)
#' considering all significant interaction principal components (IPCs) in the
#' AMMI model. Using MASI, the Simultaneous Selection Index for Yield and
#' Stability (SSI) is also calculated according to the argument
#' \code{ssi.method}.
#'
#' The Modified AMMI Stability Index (\eqn{MASI}) is computed as follows:
#'
#' \deqn{MASI = \sqrt{ \sum_{n=1}^{N'} PC_{n}^{2} \times \theta_{n}^{2}}}
#'
#' Where, \eqn{PC_{n}} are the scores of nth IPC; and \eqn{\theta_{n}} is the
#' percentage sum of squares explained by the nth principal component
#' interaction effect.
#'
#' @inheritParams MASV.AMMI
#'
#' @return A data frame with the following columns:  \item{MASI}{The MASI
#'   values.} \item{SSI}{The computed values of simultaneous selection index for
#'   yield and stability.} \item{rMASI}{The ranks of MASI values.} \item{rY}{The
#'   ranks of the mean yield of genotypes.} \item{means}{The mean yield of the
#'   genotypes.}
#'
#'   The names of the genotypes are indicated as the row names of the data
#'   frame.
#'
#' @importFrom agricolae AMMI
#' @export
#'
#' @seealso \code{\link[AMMIStbP]{MASI}}, \code{\link[AMMIStbP]{SSI}}
#'
#' @examples
#' library(agricolae)
#' data(plrv)
#'
#' # AMMI model
#' model <- with(plrv, AMMI(Locality, Genotype, Rep, Yield, console = FALSE))
#'
#' # ANOVA
#' model$ANOVA
#'
#' # IPC F test
#' model$analysis
#'
#' # Mean yield and IPC scores
#' model$biplot
#'
#' # G*E matrix (deviations from mean)
#' array(model$genXenv, dim(model$genXenv), dimnames(model$genXenv))
#'
#' # With default n (N') and default ssi.method (farshadfar)
#' MASI.AMMI(model)
#'
#' # With n = 4 and default ssi.method (farshadfar)
#' MASI.AMMI(model, n = 4)
#'
#' # With default n (N') and ssi.method = "rao"
#' MASI.AMMI(model, ssi.method = "rao")
#'
#' # Changing the ratio of weights for Rao's SSI
#' MASI.AMMI(model, ssi.method = "rao", a = 0.43)
#'
#' # ASI.AMMI same as MASI.AMMI with n = 2
#'
#' a <- ASI.AMMI(model)
#' b <- MASI.AMMI(model, n = 2)
#'
#' identical(a$ASI, b$MASI)
#'
MASI.AMMI <- function(model, n, alpha = 0.05,
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

  A <- model$biplot
  A <- A[A[,1] == "GEN",-c(1,2)]
  A <- A[,1:n] # Fetch only n IPCs

  thn <-  model$analysis[1:n,]$percent/100

  MASI <- sqrt(rowSums(as.matrix(A^2) %*% (diag(thn^2))))

  B <- model$means
  W <- aggregate(B$Yield, by = list(model$means$GEN), FUN = mean, na.rm = TRUE)
  SSI_MASI <- SSI(y = W$x, sp = MASI, gen = W$Group.1,
                 method = ssi.method, a = a)
  ranking <- SSI_MASI
  colnames(ranking) <- c("MASI", "SSI", "rMASI", "rY", "means")

  return(ranking)

}
