#' Modified AMMI Stability Value
#'
#' \code{MASV.AMMI} computes the Modified AMMI Stability Value (MASV) from a
#' modified formula of AMMI Stability Value (ASV) (Purchase et al. 1997). This
#' formula calculates AMMI stability value considering all significant principal
#' components. Using MASV, the Modified Yield stability Index (MYSI) is also
#' calculated.
#'
#' The Modified AMMI Stability Value (\deqn{MASV}) is computed as follows:
#'
#' \deqn{MASV = \sqrt{\sum\limits_{n=1}^{N-1} \frac{SSIPC_n}{SSIPC_{n+1}} \times
#' (PC_n)^2 + (PC_{n+1})^2}}
#'
#' Where, \eqn{PC_n}
#'
#' @param model
#'
#' @return
#'
#' @import agricolae tapply.stat
#' @import agricolae AMMI
#' @export
#'
#' @examples
MASV.AMMI <- function(model) {

  A <- model$biplot
  A <- A[A[,1] == "GEN",-c(1,2)]

  MASV <- rep(0, nrow(A))

  for (i in seq_along(A)) {
    pc <- model$analysis[i,4]/model$analysis[(i + 1),4]
    MASV <- MASV + apply(A, 1, function(x) pc * (x[i])^2 + (x[i + 1])^2)

    if ((i + 1) == max(seq_along(A))) break()
  }

  MASV <- sqrt(MASV)

  rk <- rank(MASV)
  B <- model$means
  W <- tapply.stat(B[,3],B[,2],function(x) mean(x,rm.na = TRUE))
  Rx <- rank(-W[,2])
  YSI <- rk + Rx
  ranking <- data.frame(MASV, MYSI, rMASV = rk, rMYSI = Rx, means = W[,2])
  invisible(ranking)
}
