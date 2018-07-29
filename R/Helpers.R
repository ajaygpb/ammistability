#' Title
#'
#' @param df
#' @param names
#' @param group
#' @param force.grouping
#' @param line.size
#' @param line.alpha
#' @param line.col
#' @param point.size
#' @param point.alpha
#' @param point.col
#' @param legend.position
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
rankslopegraph <- function(df, names, group, force.grouping = TRUE,
                           line.size = 1, line.alpha = 0.5, line.col = NULL,
                           point.size = 1, point.alpha = 0.5, point.col = NULL,
                           text.size = 2, legend.position = "bottom"){
  # check if names and group present in df

  if (missing(group) || is.null(group)) {
    ids <- names
  } else {
    ids <- c(names, group)
  }

  dfmelt <- reshape2::melt(df, id.vars = ids)


  # # consolidate duplicated ranks
  # dfmelt2 <- aggregate(. ~value+variable, data=dfmelt,
  #                     function(x) paste(unique(x), collapse = "\n"))

  dfmelt <- transform(dfmelt, lab1 = ave(get(names), variable, value,
                                         FUN = toString))
  dfmelt$lab1 <- gsub(", ", "\n", dfmelt$lab1)
  # dfmelt$lab1 <- dfmelt[, names]
  dfmelt$lab2 <- dfmelt$value

  dpcheck <- dfmelt[,c("variable","value")]
  dfmelt[duplicated(dpcheck), ]$lab1 <- NA
  dfmelt[duplicated(dpcheck), ]$lab2 <- NA

  if (missing(group) || is.null(group)) {
    if (!force.grouping) {
      gp <- NULL
    } else{
      gp <- names
    }
  } else {
    gp <- group
    dfmelt[, gp] <- as.factor(dfmelt[, gp])
  }

  slopeg <- ggplot(data = dfmelt, aes_string(x = "variable", y = "value",
                                             group = names))

  if (missing(line.col) || is.null(line.col)) {
    slopeg <- slopeg +
      geom_line(aes_string(color = gp), size = line.size, alpha = line.alpha)
  } else {
    slopeg <- slopeg +
      geom_line(aes_string(color = gp), size = line.size, alpha = line.alpha,
                colour = line.col)
  }

  if (missing(point.col) || is.null(point.col)) {
    slopeg <- slopeg +
      geom_point(aes_string(color = gp), size = point.size, alpha = point.alpha)
  } else {
    slopeg <- slopeg +
      geom_point(aes_string(color = gp), size = point.size, alpha = point.alpha,
                 colour = point.col)
  }

  slopeg <- slopeg +
    geom_text(aes_string(label = "lab1"), size = text.size,
              vjust = 0 + 0, na.rm = TRUE, nudge_y = point.size/5) +
    geom_text(aes_string(label = "lab2"), size = text.size,
              vjust = 1 - 0, na.rm = TRUE, nudge_y = -point.size/10) +
    scale_y_reverse(breaks = 1:max(dfmelt$value)) +
    ylab("Rank") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          legend.position = legend.position)

  return(slopeg)
}



#' Title
#'
#' @param df
#' @param increasing
#' @param decreasing
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
rankdf <- function(df, increasing, decreasing, ...) {
  df[, decreasing] <- lapply(df[, decreasing, drop=FALSE],
                             function(x) rank(-x, ...))
  df[, increasing] <- lapply(df[, increasing, drop=FALSE],
                             function(x) rank(x, ...))
  return(df)
}
