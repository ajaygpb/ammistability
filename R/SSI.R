#' Simultaneous Selection Indices for Yield and Stability
#'
#' \code{SSI} computes the Simultaneous Selection Index for Yield and Stability
#' (SSI) according to the methods specified in the argument \code{method}.
#'
#' The SSI proposed by Rao and Prabhakaran (2005) (\eqn{I_{i}}) is computed as
#' follows:
#'
#' \deqn{I_{i} = \frac{\overline{Y}_{i}}{\overline{Y}_{..}} + \alpha
#' \frac{\frac{1}{SP_{i}}}{\frac{1}{T}\sum_{i=1}^{T}\frac{1}{SP_{i}}}}
#'
#' Where \eqn{SP_{i}} is the stability measure of the \eqn{i}th genotype under
#' AMMI procedure;  \eqn{Y_{i}} is mean performance of \eqn{i}th genotype;
#' \eqn{Y_{..}} is the overall mean; \eqn{T} is the number of genotypes under
#' test and \eqn{\alpha} is the ratio of the weights given to the stability
#' components (\eqn{w_{2}}) and yield (\eqn{w_{1}}) with a restriction that
#' \eqn{w_{1} + w_{2} = 1}. The weights can be specified as required.
#'
#' \tabular{rrr}{ \strong{\eqn{\alpha}} \tab \strong{\eqn{w_{1}}} \tab
#' \strong{\eqn{w_{2}}}\cr 1.00 \tab 0.5 \tab 0.5\cr 0.67 \tab 0.6 \tab 0.4\cr
#' 0.43 \tab 0.7 \tab 0.3\cr 0.25 \tab 0.8 \tab 0.2\cr }
#'
#' The SSI proposed by Farshadfar (2008) is called the Genotype stability index
#' (GSI) or Yield stability index (YSI) (Farshadfar et al., 2011) and is computed
#' by summation of the ranks of the stability index/parameter and the ranks of
#' the mean yields.
#'
#' \deqn{GSI = YSI = R_{SP} + R_{Y}}
#'
#' Where, \eqn{R_{SP}} is the stability parameter/index rank of the genotype and
#' \eqn{R_{Y}} is the mean yield rank of the genotype.
#'
#' @param y A numeric vector of the mean yield/performance of genotypes.
#' @param sp A numeric vector of the stability paramter/index of the genotypes.
#' @param gen A character vector of the names of the genotypes.
#' @param method The method for the computation of simultaneous selection index.
#'   Either \code{"farshadfar"} or \code{"rao"} (See \strong{Details}).
#' @param a The ratio of the weights given to the stability components for
#'   computation of SSI when \code{method = "rao"} (See \strong{Details}).
#'
#' @return A data frame with the following columns:  \item{SP}{The stability
#'   parameter values.} \item{SSI}{The computed values of simultaneous selection
#'   index for yield and stability.} \item{rSP}{The ranks of the stability
#'   parameter.} \item{rY}{The ranks of the mean yield of genotypes.}
#'   \item{means}{The mean yield of the genotypes.}
#'
#'   The names of the genotypes are indicated as the row names of the data frame.
#'
#' @export
#'
#' @references
#'
#' \insertRef{rao_use_2005}{AMMIStbP}
#'
#' \insertRef{farshadfar_incorporation_2008}{AMMIStbP}
#'
#' \insertRef{farshadfar_ammi_2011}{AMMIStbP}
#'
#' @seealso \code{\link[AMMIStbP]{AMGE.AMMI}}, \code{\link[AMMIStbP]{ASI.AMMI}},
#'   \code{\link[AMMIStbP]{ASTAB.AMMI}}, \code{\link[AVAMGE]{AMGE.AMMI}},
#'   \code{\link[AMMIStbP]{DA.AMMI}}, \code{\link[AMMIStbP]{DZ.AMMI}},
#'   \code{\link[AMMIStbP]{EV.AMMI}}, \code{\link[AMMIStbP]{FA.AMMI}},
#'   \code{\link[AMMIStbP]{MASV.AMMI}}, \code{\link[AMMIStbP]{SIPC.AMMI}},
#'   \code{\link[AMMIStbP]{ZA.AMMI}}
#'
#'  @examples
#' library(agricolae)
#' data(plrv)
#' model <- with(plrv, AMMI(Locality, Genotype, Rep, Yield, console=FALSE))
#'
#' yield <- aggregate(model$means$Yield, by= list(model$means$GEN),
#'                FUN=mean, na.rm=TRUE)[,2]
#' stab <- DZ.AMMI(model)$DZ
#' genotypes <- rownames(DZ.AMMI(model))
#'
#' # With default ssi.method (farshadfar)
#' SSI(y = yield, sp = stab, gen = genotypes)
#'
#' # With  ssi.method = "rao"
#' SSI(y = yield, sp = stab, gen = genotypes, method = "rao")
#'
#' # Changing the ratio of weights for Rao's SSI
#' SSI(y = yield, sp = stab, gen = genotypes, method = "rao", a = 0.43)
#'
SSI <- function(y, sp, gen, method = c("farshadfar", "rao"), a = 1) {
  # Check if argument y is of type numeric
  if (!is.numeric(y)) {
    stop("'y' should be a numeric vector")
  }

  # Check if argument sp is of type numeric
  if (!is.numeric(sp)) {
    stop("'sp' should be a numeric vector")
  }

  # Check if y, sp and gen are of equal length
  if (length(y) != length(sp) | length(y) != length(gen)) {
    stop("'y', 'sp' and 'gen' lengths differ")
  }

  method <- match.arg(method)

  rk <- rank(sp)
  Rx <- rank(-y)

  if (method == "rao") {
    if (!is.numeric(a) || length(a) != 1) {
      stop("'a' should be a numeric vector of length 1")
    }
    SSI <- (y/mean(y)) + (a * ((1/sp)/mean(1/sp))) # Rao's I
  }

  if (method == "farshadfar") {
  SSI <- rk + Rx #YSI
  }

  out <- data.frame(SP = sp, SSI, rSP = rk, rY = Rx, means = y)
  rownames(out) <- gen

  return(out)

}
