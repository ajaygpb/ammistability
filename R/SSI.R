#' Simultaneous Selection Indices for Yield and Stability
#'
#' \code{SSI} computes the Simultaneous Selection Index for Yield and Stability
#' (SSI) according to the methods specified in the argument \code{method}.
#'
#' The SSI according to Rao and Prabhakaran (2005)
#' (\ifelse{html}{\out{<i>I<sub>i</sub></i>}}{\eqn{I_{i}}}) is computed as
#' follows:
#'
#' \ifelse{html}{\out{<p style="text-align: center;"><em>I<sub>i</sub> =
#' <big>[</big><sup><span style="text-decoration:overline">Y</span><sub>i</sub>
#' </sup>&frasl;<sup> </sup><sub><span
#' style="text-decoration:overline">Y</span></sub><sub>..</sub><big>]</big> +
#' &alpha; <big>[</big><sup>(1/SP<sub>i</sub>)</sup> <big>&frasl;</big>
#' <sub>((1/T) <big>&sum;</big><sup>T</sup><sub style="line-height: 1.8;
#' margin-left: -1ex;">i=1</sub>
#' (1/SPi))</sub><big>]</big></em></p>}}{\deqn{I_{i} =
#' \frac{\overline{Y}_{i}}{\overline{Y}_{..}} + \alpha
#' \frac{\frac{1}{SP_{i}}}{\frac{1}{T}\sum_{i=1}^{T}\frac{1}{SP_{i}}}}}
#'
#' Where \ifelse{html}{\out{<i>SP<sub>i</sub></i>}}{\eqn{SP_{i}}} is the
#' stability measure of the \ifelse{html}{\out{<i>i</sub></i>}}{\eqn{i}}th
#' genotype under AMMI procedure; \ifelse{html}{\out{<i><span
#' style="text-decoration:overline">Y</span><sub>i</sub></i>}}{\eqn{\overline{Y}_{i}}}
#' is mean performance of \ifelse{html}{\out{<i>i</sub></i>}}{\eqn{i}}th
#' genotype; \ifelse{html}{\out{<i><span
#' style="text-decoration:overline">Y</span><sub>..</sub></i>}}{\eqn{\overline{Y}_{..}}}
#' is the overall mean; \ifelse{html}{\out{<i>T</sub></i>}}{\eqn{T}} is the
#' number of genotypes under test and
#' \ifelse{html}{\out{<i>&alpha;</sub></i>}}{\eqn{\alpha}} is the ratio of the
#' weights given to the stability components
#' (\ifelse{html}{\out{<i>w<sub>2</sub></i>}}{\eqn{w_{2}}}) and yield
#' (\ifelse{html}{\out{<i>w<sub>1</sub></i>}}{\eqn{w_{1}}}) with a restriction
#' that \ifelse{html}{\out{<i>w<sub>1</sub> + w<sub>2</sub> = 1</i>}}{\eqn{w_{1}
#' + w_{2} = 1}}. The weights can be specified as required.
#'
#' \tabular{rrr}{ \strong{\ifelse{html}{\out{<i>&alpha;</i>}}{\eqn{\alpha}}}
#' \tab \strong{\ifelse{html}{\out{<i>w<sub>1</sub></i>}}{\eqn{w_{1}}}} \tab
#' \strong{\ifelse{html}{\out{<i>w<sub>2</sub></i>}}{\eqn{w_{2}}}}\cr 1.00 \tab
#' 0.5 \tab 0.5\cr 0.67 \tab 0.6 \tab 0.4\cr 0.43 \tab 0.7 \tab 0.3\cr 0.25 \tab
#' 0.8 \tab 0.2\cr }
#'
#' The SSI proposed by Farshadfar (2008) is called the Genotype stability index
#' (\ifelse{html}{\out{GSI}}{\eqn{GSI}}) or Yield stability index
#' (\ifelse{html}{\out{YSI}}{\eqn{YSI}}) (Farshadfar et al., 2011) and is
#' computed by summation of the ranks of the stability index/parameter and the
#' ranks of the mean yields.
#'
#' \ifelse{html}{\out{<p style="text-align: center;"><em>GSI = YSI =
#' R<sub>SP</sub> + R<sub>Y</sub></em></p>}}{\deqn{GSI = YSI = R_{SP} + R_{Y}}}
#'
#' Where, \ifelse{html}{\out{<i>R<sub>SP</sub></i>}}{\eqn{R_{SP}}} is the
#' stability parameter/index rank of the genotype and
#' \ifelse{html}{\out{<i>R<sub>Y</sub></i>}}{\eqn{R_{Y}}} is the mean yield rank
#' of the genotype.
#'
#' @param y A numeric vector of the mean yield/performance of genotypes.
#' @param sp A numeric vector of the stability parameter/index of the genotypes.
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
#'   The names of the genotypes are indicated as the row names of the data
#'   frame.
#'
#' @export
#'
#' @references
#'
#' \insertRef{rao_use_2005}{ammistability}
#'
#' \insertRef{farshadfar_incorporation_2008}{ammistability}
#'
#' \insertRef{farshadfar_ammi_2011}{ammistability}
#'
#' @seealso \code{\link[ammistability]{AMGE.AMMI}},
#'   \code{\link[ammistability]{ASI.AMMI}},
#'   \code{\link[ammistability]{ASTAB.AMMI}},
#'   \code{\link[ammistability]{AVAMGE.AMMI}},
#'   \code{\link[ammistability]{DA.AMMI}}, \code{\link[ammistability]{DZ.AMMI}},
#'   \code{\link[ammistability]{EV.AMMI}}, \code{\link[ammistability]{FA.AMMI}},
#'   \code{\link[ammistability]{MASV.AMMI}},
#'   \code{\link[ammistability]{SIPC.AMMI}},
#'   \code{\link[ammistability]{ZA.AMMI}}
#'
#' @examples
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
    SSI <- (y / mean(y)) + (a * ((1 / sp)/mean(1 / sp))) # Rao's I
  }

  if (method == "farshadfar") {
  SSI <- rk + Rx #YSI
  }

  out <- data.frame(SP = sp, SSI, rSP = rk, rY = Rx, means = y)
  rownames(out) <- gen

  return(out)

}
