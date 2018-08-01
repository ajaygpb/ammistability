#' Sum Across Environments of GEI Modelled by AMMI
#'
#' \code{AMGE.AMMI} computes the Sum Across Environments of Genotype-Environment
#' Interaction (GEI) Modelled by AMMI (AMGE) (Sneller et al., 1997) considering
#' all significant interaction principal components (IPCs) in the AMMI model.
#' Using AMGE, the Simultaneous Selection Index for Yield and Stability (SSI) is
#' also calculated according to the argument \code{ssi.method}.
#'
#' The Sum Across Environments of GEI Modelled by AMMI
#' (\ifelse{html}{\out{<i>AMGE</i>}}{\eqn{AMGE}}) is computed as follows:
#'
#' \ifelse{html}{\out{<p style="text-align: center;"><i>AMGE =
#' <big>&sum;</big><sup>E</sup><sub style="line-height: 1.8; margin-left:
#' -3ex;">j=1</sub><big>&sum;</big><sup>N'</sup><sub style="line-height: 1.8;
#' margin-left: -3ex;">n=1</sub>
#' &lambda;<sub>n</sub>&gamma;<sub>in</sub>&delta;<sub>jn</sub></i></p>}}{\deqn{AMGE
#' = \sum_{j=1}^{E} \sum_{n=1}^{N'} \lambda_{n} \gamma_{in} \delta_{jn}}}
#'
#' Where, \ifelse{html}{\out{<i>N'</i>}}{\eqn{N'}} is the number of significant
#' IPCs (number of IPC that were retained in the AMMI model via F tests);
#' \ifelse{html}{\out{<i>&lambda;<sub>n</sub></i>}}{\eqn{\lambda_{n}}} is the
#' singular value for \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}th IPC and correspondingly
#' \ifelse{html}{\out{<i>&lambda;<sup>2</sup><sub style="line-height: 1.8;
#' margin-left: -1ex;">n</sub></i>}}{\eqn{\lambda_{n}^{2}}} is its eigen value;
#' \ifelse{html}{\out{<i>&gamma;<sub>in<sub></i>}}{\eqn{\gamma_{in}}} is the
#' eigenvector value for \ifelse{html}{\out{<i>i</i>}}{\eqn{i}}th genotype; and
#' \ifelse{html}{\out{<i>&delta;<sub>jn<sub></i>}}{\eqn{\delta{jn}}} is the
#' eigenvector value for the \ifelse{html}{\out{<i>j</i>}}{\eqn{j}}th environment.
#'
#' @inheritParams MASV.AMMI
#'
#' @return A data frame with the following columns:  \item{AMGE}{The AMGE
#'   values.} \item{SSI}{The computed values of simultaneous selection index for
#'   yield and stability.} \item{rAMGE}{The ranks of AMGE values.} \item{rY}{The
#'   ranks of the mean yield of genotypes.} \item{means}{The mean yield of the
#'   genotypes.}
#'
#'   The names of the genotypes are indicated as the row names of the data
#'   frame.
#'
#' @importFrom methods is
#' @importFrom stats aggregate
#' @importFrom agricolae AMMI
#' @export
#'
#' @references
#'
#' \insertRef{sneller_repeatability_1997}{ammistability}
#'
#' @seealso \code{\link[agricolae]{AMMI}}, \code{\link[ammistability]{SSI}}
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
#' AMGE.AMMI(model)
#'
#' # With n = 4 and default ssi.method (farshadfar)
#' AMGE.AMMI(model, n = 4)
#'
#' # With default n (N') and ssi.method = "rao"
#' AMGE.AMMI(model, ssi.method = "rao")
#'
#' # Changing the ratio of weights for Rao's SSI
#' AMGE.AMMI(model, ssi.method = "rao", a = 0.43)
#'
AMGE.AMMI <- function(model, n, alpha = 0.05,
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

  # Check if n > N
  if (n > nrow(model$analysis)) {
    stop('"n" is greater than the number of IPCs in "model"')
  }

  ssi.method <- match.arg(ssi.method)

  # GxE matrix
  ge <- array(model$genXenv, dim(model$genXenv), dimnames(model$genXenv))
  # SVD
  svdge <- svd(ge)
  lambda.n <- svdge$d[1:n]
  gamma.n <- svdge$u[,1:n]
  delta.n <- svdge$v[,1:n]

  ge.n <- gamma.n %*% diag(lambda.n) %*% t(delta.n)

  AMGE <- rowSums(ge.n)

  B <- model$means
  W <- aggregate(B$Yield, by = list(model$means$GEN), FUN = mean, na.rm = TRUE)
  SSI_AMGE <- SSI(y = W$x, sp = AMGE, gen = W$Group.1,
                  method = ssi.method, a = a)
  ranking <- SSI_AMGE
  colnames(ranking) <- c("AMGE", "SSI", "rAMGE", "rY", "means")

  return(ranking)

}
