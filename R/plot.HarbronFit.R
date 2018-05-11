#' Plot method for \code{HarbronFit} objects
#' 
#' Diagnostic and goodness-of-fit plots (in 2D and 3D).
#' 
#' @param x A \code{HarbronFit} object returned by \code{\link{fitModel}}.
#' @param y A (optional) second \code{HarbronFit} object returned by 
#' \code{\link{fitModel}}. If provided, \code{which} argument is ignored and 
#' a 2d-plot comparing two model fits is produced. Note that the two models 
#' should have been fitted on the same data. Note that this argument can also 
#' be used as \code{which}. See examples.
#' @param which Whether to show default \code{\link[nlme]{plot.nls}} ("nls"), 
#' a 'slice' plot with fitted curves overlaid on top of the observed data ("2d")
#' or a 3d-plot with fitted surface overlaid on top of the observed data ("3d").
#' @param ... Further arguments passed to \code{\link[nlme]{plot.nls}}, 
#' \code{\link{fitPlot2d}} or \code{\link{fitPlot3d}} depending on \code{which}.
#' @return Output from \code{\link[nlme]{plot.nls}}, a \code{ggplot2} object if 
#' \code{which = "2d"}, or a 3d \code{rgl} plot if \code{which = "3d"}. 
#' @author Maxim Nazarov
#' @import nlme
#' @importFrom stats terms predict
#' @examples
#' \donttest{
#' data("checkerboardData", package = "drugCombo")
#' data1 <- checkerboardData[checkerboardData$exp == 1, ]
#' fitUniform <- fitModel(data1, model = "uniform")
#' fitLinear <- fitModel(data1, model = "linear1")
#' plot(fitUniform, fitLinear)
#' plot(fitLinear, "2d")  # here 2nd argument is interpreted as `which`
#' }
#' @export
plot.HarbronFit <- function(x, y = NULL, which = c("nls", "2d", "3d"), ...) {
  
  # allow interpretation of the second argument as 'which' too
  if (!is.null(y) && !inherits(y, "HarbronFit") && is.character(y)) {
    which <- y
    y <- NULL
  }
  
  if (inherits(y, "HarbronFit")) # for now we can only compare model fits in 2d 
    which <- "2d"
  
  which <- match.arg(which)
  
  if (which == "nls") {
    NextMethod(...)
  } else if (which == "2d") {
    fitPlot2d(x, fit2 = y, ...)
  } else { #3d
    fitPlot3d(x, ...)
  }
}


## Transform function for doses
log10T <- function(z, offset = NULL) 
  log10(z + 0.5 * if (!is.null(offset)) offset else min(z[z != 0]))
# NB: this is a proper inverse only when there were no 0's originally, 
#     which should be the case for tauSurface
invLog10T <- function(t, offset = NULL)  # NB: not for 0's
  10^t - if (!is.null(offset)) 1/2 * offset else 1/3 * min(10^t) 


#' Plot 2d surface (slices) of observations and model fit
#' 
#' @param fit A \code{HarbronFit} object returned by \code{\link{fitModel}}.
#' @param fit2 An optional \code{HarbronFit} object returned by 
#' \code{\link{fitModel}}. If provided, a 2d-plot comparing two model fits is 
#' produced. Note that the two models should have been fitted on the same data. 
#' Note that this argument can also be used as \code{side}. 
#' @param side Which side ("d1", "d2" or "total" for the sum of d1 and d2) to 
#' use as x-axis.
#' @param useFineGrid Whether to use fine grid for plotting fitted curves 
#' (default), or calculate predictions only at the observed data points.
#' @param modelNames Model names to use for the plot legend in the case of model 
#' comparison (i.e. when \code{fit2} is provided).  
#' @return ggplot2 object
#' @author Maxim Nazarov
#' @export
fitPlot2d <- function(fit, fit2 = NULL, side = c("d1", "d2", "total"), 
    useFineGrid = TRUE, modelNames = NULL) {

  # are we plotting two models?
  compare <- FALSE
  if (inherits(fit2, "HarbronFit"))
    compare <- TRUE
  
  # allow interpretation of the second argument as 'side' too, TODO: needed?
  if (!is.null(fit2) && !inherits(fit2, "HarbronFit") && is.character(fit2))
    side <- fit2
   
  side <- match.arg(side)
  
  if (side != "total")
    otherSide <- setdiff(c("d1", "d2"), side)
  
  # extract data
  defaultCols <- c("d1", "d2", "effect")
  
  extraCols <- extraCols2 <- c()
  if (fit$tauSpec == "symbolic") 
    extraCols <- setdiff(all.vars(fit$tauFormula), defaultCols) # TODO better?
  
  if (compare && fit2$tauSpec == "symbolic") { 
    extraCols2 <- setdiff(all.vars(fit2$tauFormula), defaultCols) # TODO better?
    extraCols <- union(extraCols, extraCols2)
  }
  
  plotData <- getData(fit, c(defaultCols, extraCols)) 
  if (compare) {
    data2 <- getData(fit2, c(defaultCols, extraCols))
    # check compatibility, can we support different data?
    # we check only common columns, might be dangerous
    commonCols <- intersect(names(plotData), names(data2))
    if (!isTRUE(all.equal(plotData[commonCols], data2[commonCols], 
            check.attributes = FALSE)))
      stop("Model fits not compatible, please compare models fitted to the same data.")
  }
  
  if (side == "total" && !("total" %in% names(plotData)))
    plotData$total <- plotData$d1 + plotData$d2
  
  # TODO: what is the proper way to fo this?
  factorsInFormula <- grepl("factor", deparse(attr(terms(fit$tauFormula), "variables")))
    
  if (useFineGrid && length(extraCols) == 0 && side != "total" && 
      fit$tauSpec != "grouping" && !factorsInFormula) {
    # FIXME: less strict?
    # get finer grid for smooth prediction lines
    # can't do it if extra columns (not d1/d2) are present
    n <- 100
    
    nonZeros <- plotData[[side]][plotData[[side]] != 0]
    offset <- min(nonZeros)
    
    fineGrid <- expand.grid(
        c(if (min(plotData[[side]]) == 0) 0 else NULL, 
            invLog10T(offset = offset, seq(from = min(log10T(nonZeros, offset)), 
                    to = max(log10T(nonZeros, offset)), length.out = n))),
        unique(plotData[[otherSide]])
    )
    names(fineGrid) <- c(side, otherSide)
    
    predictedData <- cbind(fineGrid, predict = predict(fit, newdata = fineGrid))
  } else 
    predictedData <- cbind(plotData, predict = predict(fit, newdata = plotData))
  
  if (compare) {
    factorsInFormula2 <- grepl("factor", deparse(attr(terms(fit2$tauFormula), "variables")))

    # FIXME: can we avoid repetition?
    if (useFineGrid && length(extraCols2) == 0 && side != "total" && 
        fit2$tauSpec != "grouping" && !factorsInFormula2) {
      # get finer grid for smooth prediction lines
      # can't do it if extra columns (not d1/d2) are present
      n <- 100
      
      nonZeros <- plotData[[side]][plotData[[side]] != 0]
      offset <- min(nonZeros)
      
      fineGrid <- expand.grid(
          c(if (min(plotData[[side]]) == 0) 0 else NULL, 
              invLog10T(offset = offset, seq(from = min(log10T(nonZeros, offset)), 
                      to = max(log10T(nonZeros, offset)), length.out = n))),
          unique(plotData[[otherSide]])
      )
      names(fineGrid) <- c(side, otherSide)
      
      predictedData2 <- cbind(fineGrid, predict = predict(fit2, newdata = fineGrid))
    } else 
      predictedData2 <- cbind(plotData, predict = predict(fit2, newdata = plotData))
    
  }
  
  
  p <- ggplot(data = plotData, 
          aes_string(x = paste0("log10T(", side, ")"), y = "effect")) + 
      geom_point(color = "lightgray") +
      stat_summary(fun.y = "mean", color = "black", geom = "point") +
      geom_line(data = predictedData, aes(x = log10T(get(side)), y = predict, 
              color = "one"), show.legend = isTRUE(compare))
  
  if (compare) {
    p <- p + geom_line(data = predictedData2, aes(x = log10T(get(side)), y = predict, 
            color = "two")) 
    # model names
    p <- p + scale_color_manual(name = "", 
        values = c("one" = "red", "two" = "darkgreen"), 
        label = if (!is.null(modelNames) && length(modelNames) == 2) modelNames else c("model 1", "model 2"))
  }
  
  p <- p + xlab(side) + 
      theme_bw() + 
      theme(
          legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )

  # if too many points do not show all labels, but use default axis
  maxLabels <- 15
  if (length(unique(log10T(plotData[[side]]))) <= maxLabels)
    p <- p + scale_x_continuous(breaks = unique(log10T(plotData[[side]])), 
            labels = unique(plotData[[side]]), minor_breaks = NULL)
      else
        p <- p + scale_x_continuous() + xlab(bquote(log[10](.(side))))
  
  if (length(extraCols) > 0) {
    p <- p + facet_wrap(extraCols) +
        ggtitle(paste("Observed and predicted data split by levels of", 
                paste(extraCols, collapse = ", ")))
  } else if (side != "total") { 
    maxFacetCols <- 10 # to avoid unreadable plots
    p <- p + facet_wrap(otherSide, ncol = maxFacetCols) + 
        ggtitle(paste("Observed and predicted data split by levels of", 
                otherSide)) 
  } else {
    p <- p + ggtitle("Observed and predicted data")
  }
        
  p
}

#' Plot 3d surface of observations and model fit
#' 
#' @param fit A \code{HarbronFit} object returned by 
#' \code{\link{fitModel}}.
#' @param logScale Whether to use log-scale for x and y axes.
#' @param useFineGrid Whether to use fine grid for plotting the fitted surface 
#' (default), or calculate predictions only at the observed data points.
#' @param showMesh Whether to show 'mesh' on the plot. Note: currently doesn't 
#' play nicely with the \code{useFineGrid} argument, so the latter is ignored.
#' @param widget Whether to return a "htmlwidget" object instead of plotting on 
#' a 3d device.
#' @return A 3d plot is shown or an object of class "htmlwidget" as returned by 
#' \code{\link[rgl]{rglwidget}}.
#' @import rgl
#' @author Maxim Nazarov
#' @export
fitPlot3d <- function(fit, logScale = TRUE, useFineGrid = TRUE, 
    showMesh = TRUE, widget = FALSE) {
  
  
  defaultCols <- c("d1", "d2", "effect")
  
  extraCols <- NULL
  if (fit$tauSpec == "symbolic") 
    extraCols <- setdiff(all.vars(fit$tauFormula), defaultCols) # TODO better?
  
  data <- getData(fit, c(defaultCols, extraCols)) 

  # fix offset to have the same for finegrid and dose grid 
  offset <- min(data[data$d1 != 0, "d1"], data[data$d2 != 0, "d2"])
    
  transformF <- if (logScale) function(z) log10T(z, offset = offset) else function(z) z
  invTransformF <- if (logScale) function(z) invLog10T(z, offset = offset) else function(z) z
  
  if (widget) {  # do not open 3d device
    open3d(useNULL = TRUE)
  } else {
    open3d()
  }
  
  # all points
  colorPoints = c("black", "sandybrown", "brown", "white")
  with(data, plot3d(effect ~ transformF(d1) + transformF(d2), 
          type = "s", radius = 0.05, col = colorPoints[1 + 1 * (d2 == 0) + 2 * (d1 == 0)],
          xlab = "d1", ylab = "d2", zlab = "effect"))
  
  
  uniqueDoses <- with(data, list("d1" = sort(unique(d1)), "d2" = sort(unique(d2))))
  
  # Axis labels
  xlab <- uniqueDoses$d1
  ylab <- uniqueDoses$d2
  xat <- transformF(xlab)
  yat <- transformF(ylab)
  
  # Add predicted surface
  if (showMesh || !useFineGrid) {
    # do not use fine grid, otherwise it looks strange
    doseGrid <- expand.grid(uniqueDoses)
    
    predictedGrid <- cbind(doseGrid, 
        predict = predict(fit, newdata = doseGrid))
    
    # sort manually
    predictedGrid <- predictedGrid[order(predictedGrid$d2, predictedGrid$d1), ]
    
    predictMatrix <- matrix(predictedGrid[["predict"]], 
        nrow = length(uniqueDoses$d1), ncol = length(uniqueDoses$d2))
    
  } else { # get finer grid for smooth prediction surface
    n <- 100
    fineGrid <- expand.grid(
        d1 = c(
            if (min(data$d1) == 0) 0 else NULL,
            invTransformF(
                seq(
                    from = min(transformF(data$d1[data$d1 != 0])), 
                    to = max(transformF(data$d1[data$d1 != 0])), 
                    length.out = n
                )
            )
        ),
        d2 = c(
            if (min(data$d2) == 0) 0 else NULL, 
            invTransformF(
                seq(
                    from = min(transformF(data$d2[data$d2 != 0])), 
                    to = max(transformF(data$d2[data$d2 != 0])), 
                    length.out = n
                )
            )
        )
    )
    
    predictedGrid <- cbind(fineGrid, 
        predict = predict(fit, newdata = fineGrid))
    
    # sort manually
    predictedGrid <- predictedGrid[order(predictedGrid$d2, predictedGrid$d1), ]
    
    uniqueDoses <- with(predictedGrid, 
        list("d1" = sort(unique(d1)), "d2" = sort(unique(d2))))

    predictMatrix <- matrix(predictedGrid[["predict"]], 
        nrow = length(uniqueDoses$d1), ncol = length(uniqueDoses$d2))
    
  }  

  # plot surface
  persp3d(transformF(uniqueDoses$d1), transformF(uniqueDoses$d2), predictMatrix,
      add = TRUE, col = "grey70", alpha = 0.5, lit = FALSE, lwd = 2)
#  bbox3d(xat = xat, yat = yat, xlab = xlab, ylab = ylab, 
#      front = "lines", back = "cull")
  
  # add axis labels
  rgl.bbox(xlen = 0, ylen = 0, zlen = 0, # no tickmarks
      color = "#000000", front = "lines", back = "cull")
  axis3d(edge = "x--", at = xat, labels = xlab)
  axis3d(edge = "y--", at = yat, labels = ylab)
  axis3d(edge = "z--")
  
  # add mesh
  if (showMesh)
    persp3d(transformF(uniqueDoses$d1), transformF(uniqueDoses$d2), predictMatrix,
        add = TRUE, col = "grey70", alpha = 0.5, lit = FALSE, lwd = 2)
  
  aspect3d(1)  # show the bounding box as a cube
  
  if (widget)
    rglwidget()
}