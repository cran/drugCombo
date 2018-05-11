
#' Inverse 4-parameter log-logistic (response-dose) function 
#' 
#' @param y Response. 
#' @param h Hill's coefficient (slope of the curve).
#' @param b Baseline effect (at zero dose).
#' @param m Maximal/asymptote effect (at infinite dose).
#' @param logEC50 Point of inflection (in logarithmic terms).
#' @return Dose level. 
invL4 <- function(y, h, b, m, logEC50) {
  
  exp(logEC50) * ((y-b)/(m-y))^(1.0/abs(h))
  
}

# Plot the data (NB: works only for checkerboard-design data)
plotData <- function(data, side = c("d1", "d2")) {
    
  side <- match.arg(side)
  otherSide <- setdiff(c("d1", "d2"), side)
  
  # check names
  if (!all(c("d1", "d2") %in% names(data))) {
    if (all(c("A", "B") %in% names(data))) {
      data$d1 <- data$A
      data$d2 <- data$B
    } else {
      stop("Please name doses either `d1` and `d2` or `A` and `B`.")
    }
  }
  if (!"effect" %in% names(data)) {
    if ("response" %in% names(data)) {
      data$effect <- data$response
    } else {
      stop("Please name response either `effect` or `response`.")
    }
  }
  
  p <- ggplot(data = data, 
          aes_string(x = paste0("log10T(", side, ")"), y = "effect")) + 
      geom_point() +
      stat_summary(fun.y = "mean", color = "gray", geom = "line") +
      facet_wrap(otherSide, nrow = 1) + 
      scale_x_continuous(breaks = unique(log10T(data[[side]])), 
          labels = unique(data[[side]]), minor_breaks = NULL) +  
      ggtitle(paste("Observed data by levels of", otherSide)) +
      xlab(side) + 
      theme_bw() + 
      theme(
          legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  
  p
}

getData <- function(fit, columns = NULL, ...) {
  stopifnot(inherits(fit, "HarbronFit")) # nls is enough?
  if (is.null(columns))
    columns <- setdiff(ls(fit$m$getEnv()), names(coef(fit)))
  
  colList <- mget(columns, envir = fit$m$getEnv(), inherits = FALSE, 
      ifnotfound = list(NULL))
  # remove columns not found
  colList <- Filter(Negate(is.null), colList)
  
  data.frame(colList, ...)
  
}