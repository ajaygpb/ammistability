#' Stability Measure Based on Fitted AMMI Model
#'
#' \code{FA.AMMI} computes the Stability Measure Based on Fitted AMMI Model (FA)
#' (Raju, 2002) considering all significant interaction principal components
#' (IPCs) in the AMMI model. Using FA, the Simultaneous Selection Index for
#' Yield and Stability (SSI) is also calculated according to the argument
#' \code{ssi.method}.
#'
#' The Stability Measure Based on Fitted AMMI Model
#' (\ifelse{html}{\out{<i>FA</i>}}{\eqn{FA}}) is computed as follows:
#'
#' \ifelse{html}{\out{<p style="text-align: center;"><i>FA =
#' <big>&sum;</big><sup>N'</sup><sub style="line-height: 1.8; margin-left:
#' -2ex;">n=1</sub> &lambda;<sup>2</sup><sub style="line-height: 1.8;
#' margin-left: -1ex;">n</sub>&gamma;<sup>2</sup><sub style="line-height: 1.8;
#' margin-left: -1ex;">in</sub></i></p>}}{\deqn{FA =
#' \sum_{n=1}^{N'}\lambda_{n}^{2}\gamma_{in}^{2}}}
#'
#' Where, \ifelse{html}{\out{<i>N'</i>}}{\eqn{N'}} is the number of significant
#' IPCs (number of IPC that were retained in the AMMI model via F tests);
#' \ifelse{html}{\out{<i>&lambda;<sub>n</sub></i>}}{\eqn{\lambda_{n}}} is the
#' singular value for \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}th IPC and
#' correspondingly \ifelse{html}{\out{<i>&lambda;<sup>2</sup><sub
#' style="line-height: 1.8; margin-left:
#' -1ex;">n</sub></i>}}{\eqn{\lambda_{n}^{2}}} is its eigen value; and
#' \ifelse{html}{\out{<i>&gamma;<sub>in<sub></i>}}{\eqn{\gamma_{in}}} is the
#' eigenvector value for \ifelse{html}{\out{<i>i</i>}}{\eqn{i}}th genotype.
#'
#' When \ifelse{html}{\out{<i>N'</i>}}{\eqn{N'}} is replaced by 1 (only first
#' IPC axis is considered for computation), then the parameter
#' \ifelse{html}{\out{<i>FP</i>}}{\eqn{FP}} can be estimated (Zali et al.,
#' 2012).
#'
#' \ifelse{html}{\out{<p style="text-align: center;"><i>FP
#' =&lambda;<sup>2</sup><sub style="line-height: 1.8; margin-left:
#' -1ex;">1</sub>&gamma;<sup>2</sup><sub style="line-height: 1.8; margin-left:
#' -1ex;">i1</sub></i></p>}}{\deqn{FP = \lambda_{1}^{2}\gamma_{i1}^{2}}}
#'
#' When \ifelse{html}{\out{<i>N'</i>}}{\eqn{N'}} is replaced by 2 (only first
#' two IPC axes are considered for computation), then the parameter
#' \ifelse{html}{\out{<i>B</i>}}{\eqn{B}} can be estimated  (Zali et al., 2012).
#'
#' \ifelse{html}{\out{<p style="text-align: center;"><i>B =
#' <big>&sum;</big><sup>2</sup><sub style="line-height: 1.8; margin-left:
#' -1ex;">n=1</sub> &lambda;<sup>2</sup><sub style="line-height: 1.8;
#' margin-left: -1ex;">n</sub>&gamma;<sup>2</sup><sub style="line-height: 1.8;
#' margin-left: -1ex;">in</sub></i></p>}}{\deqn{B =
#' \sum_{n=1}^{2}\lambda_{n}^{2}\gamma_{in}^{2}}}
#'
#' When \ifelse{html}{\out{<i>N'</i>}}{\eqn{N'}} is replaced by
#' \ifelse{html}{\out{<i>N</i>}}{\eqn{N}} (All the IPC axes are considered for
#' computation), then the parameter estimated is equivalent to Wricke's
#' ecovalence (\ifelse{html}{\out{<i>W<sub>(AMMI)</sub></i>}}{\eqn{W_{(AMMI)}}})
#' (Wricke, 1962; Zali et al., 2012).
#'
#' \ifelse{html}{\out{<p style="text-align: center;"><i>W<sub>(AMMI)</sub> =
#' <big>&sum;</big><sup>N</sup><sub style="line-height: 1.8; margin-left:
#' -1ex;">n=1</sub> &lambda;<sup>2</sup><sub style="line-height: 1.8;
#' margin-left: -1ex;">n</sub>&gamma;<sup>2</sup><sub style="line-height: 1.8;
#' margin-left: -1ex;">in</sub></i></p>}}{\deqn{W_{(AMMI)} =
#' \sum_{n=1}^{N}\lambda_{n}^{2}\gamma_{in}^{2}}}
#'
#' @inheritParams MASV.AMMI
#'
#' @return A data frame with the following columns:  \item{FA}{The FA values.}
#'   \item{SSI}{The computed values of simultaneous selection index for yield
#'   and stability.} \item{rFA}{The ranks of FA values.} \item{rY}{The ranks of
#'   the mean yield of genotypes.} \item{means}{The mean yield of the
#'   genotypes.}
#'
#'   The names of the genotypes are indicated as the row names of the data
#'   frame.
#'
#' @importFrom agricolae AMMI
#' @export
#'
#' @references
#'
#' \insertRef{wricke_method_1962}{AMMIStbP}
#'
#' \insertRef{raju_study_2002}{AMMIStbP}
#'
#' \insertRef{zali_evaluation_2012}{AMMIStbP}
#'
#' @seealso \code{\link[AMMIStbP]{SSI}}
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
#' FA.AMMI(model)
#'
#' # With n = 4 and default ssi.method (farshadfar)
#' FA.AMMI(model, n = 4)
#'
#' # With default n (N') and ssi.method = "rao"
#' FA.AMMI(model, ssi.method = "rao")
#'
#' # Changing the ratio of weights for Rao's SSI
#' FA.AMMI(model, ssi.method = "rao", a = 0.43)
#'
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

  FA <- rowSums((gamma.n^2) %*% diag(lambda.n^2))

  B <- model$means
  W <- aggregate(B$Yield, by = list(model$means$GEN), FUN = mean, na.rm = TRUE)
  SSI_FA <- SSI(y = W$x, sp = FA, gen = W$Group.1,
                method = ssi.method, a = a)
  ranking <- SSI_FA
  colnames(ranking) <- c("FA", "SSI", "rFA", "rY", "means")

  return(ranking)

}
