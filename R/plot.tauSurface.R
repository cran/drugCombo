#' Plot method for "tauSurface" objects
#'
#' 2D slice and 3D surface plots for the interaction index surface.
#' 
#' @param x A \code{tauSurface} object returned by \code{\link{getTauSurface}}.
#' @param y An optional second \code{tauSurface} object returned by 
#' \code{\link{getTauSurface}}. If provided, \code{which} argument is ignored 
#' and a 2d-plot comparing two tau surfaces is produced. Note that the two 
#' estimates should have been calculated on the same \code{HarbronFit} object. 
#' Note that this argument can also be used as \code{which}. See examples.
#' @param which Whether to show a 3d plot (surface plot) or a 2d plot 
#' (slice plot). Correspondingly \code{\link{tauPlot3d}} or 
#' \code{\link{tauPlot2d}} is called.
#' @param ... Further arguments passed to \code{\link{tauPlot3d}} or 
#' \code{\link{tauPlot2d}}.
#' @seealso \code{\link{tauPlot3d}}, \code{\link{tauPlot2d}} for the underlying 
#' functions and their arguments. \code{\link{contour.tauSurface}} for another 
#' visual representation of the estimated interaction indices.
#' @importFrom methods hasArg
#' @examples
#' \donttest{
#' data("checkerboardData", package = "drugCombo")
#' data1 <- checkerboardData[checkerboardData$exp == 1, ]
#' fitUniform <- fitModel(data1, model = "uniform")
#' tauUniform <- getTauSurface(fitUniform)
#' fitLinear <- fitModel(data1, model = "linear1")
#' tauLinear <- getTauSurface(fitLinear)
#' plot(tauUniform)
#' plot(tauLinear, which = "2d", side = "d2", facetBy = "d1")
#' plot(tauLinear, which = "3d")
#' plot(tauUniform, tauLinear, tauNames = c("uniform", "linear"))
#' plot(tauUniform, tauLinear, continuous2 = FALSE)
#' }
#' @export
plot.tauSurface <- function(x, y = NULL, which = c("2d", "3d"), ...) {
  
  # allow interpretation of the second argument as 'which' too
  if (!is.null(y) && !inherits(y, "tauSurface") && is.character(y)) {
    which <- y
    y <- NULL
  }
  
  compare <- FALSE
  if (inherits(y, "tauSurface")) {# for now we can only compare tau surfaces in 2d 
    which <- "2d"
    compare <- TRUE
    
    # check compatibility
    uniqueTau <- function(tauSurface) {
      res <- tauSurface$data
      if (length(notUsed <- setdiff(c("d1", "d2"), all.vars(tauSurface$tauFormula))))
        res <- res[!duplicated(res[, setdiff(names(res), notUsed)]), ]
      
      res
    }
    
    tauData1 <- uniqueTau(x)
    tauData2 <- uniqueTau(y)
    
    tauCols <- c("tau", "se.est", "lower", "upper")
    dataCols1 <- setdiff(colnames(tauData1), tauCols)
    dataCols2 <- setdiff(colnames(tauData2), tauCols)
    
    if (!setequal(dataCols1, dataCols2))
      stop("Tau estimates not compatible, please compare estimates calculated on the same model fit.")
    
  }
    
  which <- match.arg(which)

  funs <- NULL
  if (!is.null(x$surfaceFormula))
    funs <- c(funs, tau = x$surfaceFormula)
  if (!is.null(x$ciFormula))
    funs <- c(funs, 
        lower = x$ciFormula(side = 'lower'), 
        upper = x$ciFormula(side = 'upper'))
  
  if (which == "3d") {
    
    if (hasArg("continuous")) {
      tauPlot3d(x$data, funs = funs, ...)
    } else {
      tauPlot3d(x$data, continuous = isTRUE(!is.null(funs)), 
          funs = funs, ...)
    }
    
  } else {
    
    if (compare) {
      
      # check if can 'simplify' the plots
      plotData  <- unique(x)
      plotData2 <- unique(y)
      defaultCols <- c("se.est", "lower", "upper")
      sameSide <- setequal(
          setdiff(names(plotData), defaultCols), 
          setdiff(names(plotData2), defaultCols)
      )
      
      if (!sameSide) {
        # restore all data
        plotData <- x$data
        plotData2 <- y$data
      }
      
      funs2 <- NULL
      if (!is.null(y$surfaceFormula))
        funs2 <- c(funs2, tau = y$surfaceFormula)
      if (!is.null(y$ciFormula))
        funs2 <- c(funs2, 
            lower = y$ciFormula(side = 'lower'), 
            upper = y$ciFormula(side = 'upper'))
      
      # FIXME: this is ugly :(
      if (hasArg("continuous")) {
        if (hasArg("continuous2")) {
          tauPlot2d(plotData, plotData2, funs = funs, funs2 = funs2, ...)
        } else {
          tauPlot2d(plotData, plotData2, continuous2 = !is.null(funs2), 
              funs = funs, funs2 = funs2, ...)
        }
      } else {
        if (hasArg("continuous2")) {
          tauPlot2d(plotData, plotData2, continuous = !is.null(funs),
              funs = funs, funs2 = funs2, ...)
        } else {
          tauPlot2d(plotData, plotData2, continuous = !is.null(funs),
              continuous2 = !is.null(funs2), funs = funs, funs2 = funs2, ...)
        }
      }

    } else {
      
      plotData <- unique(x)

      if (hasArg("continuous")) {
        tauPlot2d(plotData, funs = funs, ...)
      } else {
        tauPlot2d(plotData, continuous = !is.null(funs), funs = funs, ...)
      }
    
    }
  }
}

#' Plot estimated interaction index surface slice along one of the doses
#' 
#' @details The function returns a 2d plot for the interaction index (tau) 
#' estimates as a function of one of the two doses in a checkerboard design, or 
#' rays in a ray design. Pointwise confidence intervals are displayed as error 
#' bars. In addition to plotting tau estimates from one \code{tauSurface} object, 
#' the function can be used to compare two \code{tauSurface} objects. 
#' This can be used, for example, to see the difference between Wald-type and 
#' boostrap-based confidence intervals for tau.
#' Although the function can be used 'manually', typically one calls the 
#' \code{\link{plot.tauSurface}} method, which then calls this function when
#' \code{which = "2d"}.
#' 
#' @param tauSurface A \code{tauSurface} object returned by 
#' \code{\link{getTauSurface}}.
#' @param tauSurface2 An optional second \code{tauSurface} object returned by 
#' \code{\link{getTauSurface}}. If provided, a 2d-plot comparing the two tau 
#' surfaces is produced. Note that the estimates should have been calculated on 
#' the same \code{HarbronFit} object.
#' @param side Data column to use as x-axis: "d1", "d2", "total" or another 
#' variable from the data in the \code{tauSurface} object.
#' @param groupBy Data column to use as grouping. Note that if comparison of two 
#' surfaces is performed, this will be ignored.
#' @param colorBy Data column to use for coloring.
#' @param facetBy Whether to facet plots by extra variables used in the tau 
#' formula ("auto") or manually provided data column(s) to facet by.
#' @param continuous Whether continuous type of plot is requested (for 
#' \code{tauSurface1}). This is automatically detected if used via 
#' \code{\link{plot.tauSurface}}, but can be overwritten. 
#' @param continuous2 Whether continuous type of plot is requested (for 
#' \code{tauSurface2}). This is automatically detected if used via 
#' \code{\link{plot.tauSurface}}, but can be overwritten. 
#' @param addLine Whether to connect tau estimates for subsequent doses.
#' @param funs A list with functions to compute tau surface and confidence 
#' bands (for \code{tauSurface1}). These are returned by the 
#' \code{\link{getTauSurface}} and are automatically used when the 
#' \code{\link{plot.tauSurface}} is called.
#' @param funs2 A list with functions to compute tau surface and confidence 
#' bands (for \code{tauSurface2}). These are returned by the 
#' \code{\link{getTauSurface}} and are automatically used when the 
#' \code{\link{plot.tauSurface}} is called.
#' @param title Plot title.
#' @param tauNames Tau surface names to use for the plot legend in the case of 
#' comparison of estimates (i.e. when \code{tauSurface2} is provided).  
#' @param digits Number of digits used in \code{\link[base]{format}} for the 
#' dose labels.
#' @param facetOpts Arguments passed to \code{\link[ggplot2]{facet_wrap}}
#' @return A ggplot2 object.
#' @author Maxim Nazarov
#' @seealso \code{\link{plot.tauSurface}}, \code{\link{tauPlot3d}}
#' @export
tauPlot2d <- function(tauSurface, tauSurface2 = NULL, side = "d1", 
    groupBy = NULL, colorBy = groupBy, facetBy = "auto",
    continuous = FALSE, continuous2 = FALSE, addLine = continuous || continuous2, 
    funs = NULL, funs2 = NULL, title = NULL, tauNames = NULL, digits = 4,
    facetOpts = NULL) {

  tauCompareColumn <- "tau.est"
  compare <- FALSE

  # allow interpretation of the second argument as 'side' too
  if (!is.null(tauSurface2) && !inherits(tauSurface2, "data.frame") & is.character(tauSurface2))
    side <- tauSurface2
  
  tauSurface1 <- tauSurface  # keep for separate plotting
  
  if (inherits(tauSurface2, "data.frame")) {
    
    compare <- TRUE

    # add any of the default cols (se.est, lower, upper) that are missing
    # this can happen if a tau ests has 'boot' method or if addCI = FALSE was 
    # requested
    defaultCols <- c("se.est", "lower", "upper")
    if (length(missingCols1 <- setdiff(defaultCols, names(tauSurface1))))
      tauSurface1[, missingCols1] <- NA_real_
    if (length(missingCols2 <- setdiff(defaultCols, names(tauSurface2))))
      tauSurface2[, missingCols2] <- NA_real_

    # add comparison names
    tauNames <- if (!is.null(tauNames) && length(tauNames) == 2) tauNames else 
          c("model 1", "model 2")
    tauSurface1[[tauCompareColumn]] <- tauNames[1]
    tauSurface2[[tauCompareColumn]] <- tauNames[2]
    
    # combine the two
    tauSurface <- rbind(tauSurface1, tauSurface2)  
    
    if (!is.null(groupBy))
      warning("When comparing two tau fits, `groupBy` is ignored", 
          call. = FALSE)
    
    groupBy <- tauCompareColumn # FIXME: more flexible?
  
  }

  otherSide <- NULL
  
  if (side %in% c("d1", "d2")) {
    
    otherSide <- setdiff(c("d1", "d2"), side)
    
    if (!side %in% names(tauSurface)) {
      sideName <- side # for axis label
      side <- NA
      continuous <- FALSE
    }
    
    if (!otherSide %in% names(tauSurface))
      otherSide <- NULL
    
    if (is.na(side) && !is.null(otherSide))
      warning("Plotted data values are the same for all doses of ", sideName, ".",
          " Consider plotting with `side = \"", otherSide, "\"`.", 
          call. = FALSE)
    
  } else if (side == "total" && all(c("d1", "d2") %in% names(tauSurface))) {
    tauSurface$total <- tauSurface$d1 + tauSurface$d2
  } else if (!side %in% names(tauSurface)) {
    stop("No column called '", side, "' found for the 'side' argument.")
  }
  
  groupVar <- otherSide
  if (!is.null(groupBy)) {
    if (!groupBy %in% names(tauSurface))
      stop("No column called '", groupBy, "' found for the 'groupBy' argument.")
    
    if (compare) {
      groupVar <- c(groupVar, groupBy) 
    } else {
      groupVar <- groupBy
    }
          
  }
  colorVar <- otherSide
  if (!is.null(colorBy)) {
    if (!colorBy %in% names(tauSurface))
      stop("No column called '", colorBy, "' found for the 'colorBy' argument.")
    
    colorVar <- colorBy
  }
  
  # make sure both 1 or 2 group variables are supported
  groupVarStr <- if (!is.null(groupVar))
    paste0("interaction(`", paste0(groupVar, collapse = "`, `"), "`)") else NULL
  
  p <- ggplot(data = tauSurface, aes_string(x = side, y = "tau",
          group = groupVarStr, col = colorVar, fill = colorVar))
  
  # adding smooth CIs
  getFineGrid <- function(tauSurface, side, groupVar, funs, n = 51) {
    fineGrid <- expand.grid(
        c(list(10^seq(
            from = log10(min(tauSurface[[side]], na.rm = TRUE)), 
            to = log10(max(tauSurface[[side]], na.rm = TRUE)), 
            length.out = n)),
        if (!is.null(groupVar)) unique(tauSurface[groupVar]) else NA)
    )
    names(fineGrid) <- c(side, groupVar)
    
    fineGrid$tau <- if ("tau" %in% names(funs)) funs$tau(fineGrid$d1, fineGrid$d2) else NA_real_
    fineGrid$lower <- if ("lower" %in% names(funs)) funs$lower(fineGrid$d1, fineGrid$d2) else NA_real_
    fineGrid$upper <- if ("upper" %in% names(funs)) funs$upper(fineGrid$d1, fineGrid$d2) else NA_real_
    
    fineGrid
  }
  
  fineGrid <- if (continuous && !is.null(funs) && !is.na(side)) {
        getFineGrid(tauSurface1, side, groupVar, funs)
      } else tauSurface1
  
  fineGrid2 <- if (continuous2 && !is.null(funs2) && !is.na(side)) {
        getFineGrid(tauSurface2, side, groupVar, funs2)
      } else tauSurface2
  
#  if (compare) { # subset for potentially different display (continuous/discrete)
#    fineGrid  <- fineGrid[fineGrid[[tauCompareColumn]] == tauNames[1], ]
#    fineGrid2 <- fineGrid2[fineGrid2[[tauCompareColumn]] == tauNames[2], ]
#  }
  
  if (is.numeric(tauSurface[[side]])) {
    minDiff <- min(abs(diff(log10(tauSurface[[side]]))))
    dodgeWidth <- min(0.2, minDiff/1.5)
  } else { 
    dodgeWidth <- 0.2
  }
  errorbarWidth <- 0.75*dodgeWidth # this must be > than errorbarWidth
  
  # add confidence intervals/bands
  if (all(c("lower", "upper") %in% names(tauSurface))) {

    p <- p + aes_string(ymin = "lower", ymax = "upper")
    
    # accommodate different combinations of discrete/continuous intervals
    if (!continuous) { 
      
      if (compare) { 
        if (!continuous2) { 
          p <- p + geom_errorbar(width = errorbarWidth, position = position_dodge(width = dodgeWidth), na.rm = TRUE)
        } else {
          p <- p + geom_ribbon(data = fineGrid2, alpha = 0.4, color = NA)
          p <- p + geom_errorbar(data = fineGrid, width = errorbarWidth, na.rm = TRUE)
        }
      } else {
        p <- p + geom_errorbar(width = errorbarWidth, na.rm = TRUE)
      }
      
    } else {

      p <- p + geom_ribbon(data = fineGrid, alpha = 0.4, color = NA)
      
      if (compare) { 
        if (!continuous2) { 
          p <- p + geom_errorbar(data = fineGrid2, width = 0.15, na.rm = TRUE)
        } else {
          p <- p + geom_ribbon(data = fineGrid2, alpha = 0.4, color = NA)
        }
      }
      
    }
    
  }
  
  # connect the points
  if (addLine && !is.na(side)) {
    if (!continuous) {
      warning("Adding line to discretely specified tau-function, interpret with caution.", 
          call. = FALSE, immediate. = TRUE)
      p <- p + geom_line(data = if (compare) tauSurface1 else tauSurface, linetype = 'dotted')
    } else
      p <- p + geom_line(data = fineGrid)
    
    if (compare) {
      if (!continuous2) {
        warning("Adding line to discretely specified tau-function, interpret with caution.", 
            call. = FALSE, immediate. = TRUE)
        p <- p + geom_line(data = if (compare) tauSurface2 else tauSurface, linetype = 'dotted')
      } else 
        p <- p + geom_line(data = fineGrid2)
    }
  }
  
  # add additivity (tau = 1) line
  p <- p + geom_point(position = if (compare && !continuous && !continuous2) 
                position_dodge(width = dodgeWidth) else "identity") + 
      geom_hline(yintercept = 1)
  
  # add custom labels, if not too many
  maxLabels <- 12
  
  if (!is.na(side)) {
    if (is.numeric(tauSurface[[side]])) { # we allow also character/factor variables for 'side'
      sideLabels <- unique(tauSurface[[side]])
      fewSideLabels <- (length(sideLabels) <= maxLabels)
      p <- p + scale_x_log10(
          breaks = if (fewSideLabels) sideLabels else waiver(),
          labels = if (fewSideLabels) format(sideLabels, digits = digits) else waiver(),
          minor_breaks = if (fewSideLabels) NULL else waiver())
    }
  } else
    p <- p + scale_x_discrete(labels = "", name = if (is.null(groupVar)) "" else sideName)
  
  if (!is.null(colorVar)) {
    colorLabels <- unique(tauSurface[[colorVar]])
    fewColorLabels <- (length(colorLabels) <= maxLabels)
    
    if (is.numeric(tauSurface[[colorVar]]))
      p <- p + scale_colour_continuous(trans = 'log',
              breaks = if (fewColorLabels) colorLabels else waiver()) +
          scale_fill_continuous(trans = 'log',
              breaks = if (fewColorLabels) colorLabels else waiver())
    else 
      p <- p + scale_colour_discrete(breaks = colorLabels) +
          scale_fill_discrete(breaks = colorLabels)
  }
  
  if (!is.null(facetBy)) {
    if (length(facetBy) == 1 && facetBy == "auto") {
      defaultCols <- c("d1", "d2", "tau", "se.est", "lower", "upper", "total", tauCompareColumn)
      facetBy <- setdiff(names(tauSurface), union(defaultCols, side))
    } else if (length(missingCols <- setdiff(facetBy,  names(tauSurface))) > 0) {
      stop("No column", if (length(missingCols) > 1) "s" , 
          " called '", paste(missingCols, collapse = "', '"), 
          "' found for the 'facetBy' argument.")
    }
    if (length(facetBy) > 0)
      p <- p + do.call(facet_wrap, list(facets = facetBy, facetOpts))
  }
    
  if (!is.null(title))
    p <- p + ggtitle(title)
  
  p + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}


#' Plot 3d surface of interaction index estimates
#' 
#' @details The function returns a 3d plot for the interaction index (tau) 
#' estimates as a function of the doses of the two drugs. Pointwise confidence 
#' intervals are displayed as error bars.
#' Although the function can be used 'manually', typically one calls the 
#' \code{\link{plot.tauSurface}} method, which then calls this function when 
#' \code{which = "3d"}.
#' 
#' @param tauSurface A \code{tauSurface} object returned by 
#' \code{\link{getTauSurface}}.
#' @param logScale Whether to use log-scale for x and y axes.
#' @param continuous Whether continuous type of plot is requested. This is 
#' automatically detected if used via \code{\link{plot.tauSurface}}, but can be 
#' overwritten. 
#' @param funs A list with functions to compute tau surface and confidence 
#' bands. These are returned by the \code{\link{getTauSurface}} and are 
#' automatically used when the \code{\link{plot.tauSurface}} is called.
#' @param addPlane Whether to add estimated tau plane.
#' @param colorPoints Whether to color points by synergy/antagonism. Blue color 
#' is used for points deemed synergistic (if the confidence interval lies 
#' entirely below 1), red for points deemed antagostic' (if the confidence 
#' interval lies entirely above 1). Other points are colored white.
#' @param widget Whether to return a "htmlwidget" object instead of plotting on 
#' 3d device.
#' @return A 3d plot is shown or an object of the class "htmlwidget" as returned by 
#' \code{\link[rgl]{rglwidget}}.
#' @import rgl
#' @author Maxim Nazarov
#' @seealso \code{\link{plot.tauSurface}}, \code{\link{tauPlot2d}}
#' @export
tauPlot3d <- function(tauSurface, logScale = TRUE, 
    continuous = FALSE, funs = NULL, 
    addPlane = continuous, colorPoints = TRUE, widget = FALSE) {
  
  ## Transform function for doses
  log10T <- function(z) log10(z + 0.5 * min(z[z != 0]))
  transformF <- if (logScale) log10T else function(z) z
  # NB: this is a proper inverse only when there were no 0's originally, 
  #     which should be the case for tauSurface
  invLog10T <- function(t) 10^t - 1/3 * min(10^t) 
  invTransformF <- if (logScale) invLog10T else function(z) z
  
  synColors <- if (colorPoints && all(c("lower", "upper") %in% names(tauSurface)))  
        1 + 1*(tauSurface$lower > 1) + 2*(tauSurface$upper < 1) else 1
  colorPallette <- c("grey70", "red", "blue")
  synColors <- colorPallette[synColors]
  uniqueDoses <- with(tauSurface, list("d1" = sort(unique(d1)), "d2" = sort(unique(d2))))
  doseGrid <- expand.grid(uniqueDoses)
  plotMatrices <- lapply(merge(doseGrid, tauSurface, sort = FALSE), matrix, 
      nrow = length(uniqueDoses$d1), ncol = length(uniqueDoses$d2))
  
  if (widget) {  # do not open 3d device
    open3d(useNULL = TRUE)
  } else {
    open3d()
  }
  
  # points
  plot3d(tauSurface$tau~transformF(tauSurface$d1)+transformF(tauSurface$d2), 
      type = "s", size = 1, col = synColors,
      xlab = "d1", ylab = "d2", zlab = "tau")
  
  # Add axes
  xlab <- uniqueDoses$d1
  xat <- transformF(xlab)
  ylab <- uniqueDoses$d2
  yat <- transformF(ylab)
  bbox3d(xat = xat, yat = yat, 
      xlab = xlab, ylab = ylab, 
      front = "line", back = "line")
  
  # draw error-bars
  if (all(c("lower", "upper") %in% names(tauSurface)))
    segments3d(scale = FALSE,
        rep(transformF(tauSurface$d1), each = 2), 
        rep(transformF(tauSurface$d2), each = 2), 
        as.vector(t(cbind(tauSurface$lower, tauSurface$upper))),
        col = rep(synColors, each = 2),
        lwd = 2)
    
  # add additivity surface
  planes3d(a = 0, b = 0, c = 1, d = -1, col = "gray", alpha = 0.5)
  
  # add plane connecting points
  if (addPlane) {
    
    if (!isTRUE(continuous)) {
      warning("Adding a plane to discretely specified tau-function, interpret with caution.", 
          call. = FALSE, immediate. = TRUE)
      
      persp3d(transformF(uniqueDoses$d1), transformF(uniqueDoses$d2), plotMatrices$tau,
          add = TRUE, col = "gray70", alpha = 0.5)
      
    } else {
      
      if (!is.null(funs) && "tau" %in% names(funs)) {
        
        # continuous functions
        persp3d(function(d1, d2) funs$tau(invLog10T(d1), invLog10T(d2)),
            xlim = range(xat), ylim = range(yat), 
            add = TRUE, col = "gray70", alpha = 0.5, lit = FALSE)
        
      } else { # use default doses
        
        persp3d(transformF(uniqueDoses$d1), transformF(uniqueDoses$d2), plotMatrices$tau,
            add = TRUE, col = "gray70", alpha = 0.5)
        
      }
      
    }
  }
  
  aspect3d(1)  # show the bounding box as a cube 
  
  if (widget)
    rglwidget()
}


#' Contour-plot method for "tauSurface" objects
#'
#' @details The function returns a heatmap-like plot displaying the interaction 
#' index (tau) estimates and pointwise confidence intervals numerically. 
#' Blue and red colours are used to indicate areas of synergy (if confidence 
#' interval lies entirely below 1) or antagonism (if confidence interval lies 
#' entirely above 1). Note that color intensity is determined by the absolute 
#' values of the interaction indices.
#' 
#' @param x A \code{tauSurface} object returned by \code{\link{getTauSurface}}.
#' @param digits Number of digits used in \code{\link[base]{format}} for the 
#' dose labels.
#' @param ... Further parameters, currently not used.
#' @return a ggplot2 object
#' @import ggplot2
#' @author Maxim Nazarov
#' @export
contour.tauSurface <- function(x, digits = 4, ...) {
  
  plotData <- x$data
  
  ## make manual color scale
  # use RColorBrewer's brewer.pal(8, "Reds") and brewer.pal(8, "Blues")
  redPalette <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", 
      "#EF3B2C", "#CB181D", "#99000D")
  bluePalette <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", 
      "#4292C6", "#2171B5", "#084594")
  allColors <- c(rev(bluePalette), "#FFFFFF", redPalette)
  nCols <- length(allColors) # should be odd
  allColors <- setNames(allColors, seq_len(nCols))

  # determine color based on distance from the closest CI limit to 1,
  # also limit range to 0 - 2 (since antagonism can go higher)
  # handle 0 (additivity) individually
  plotData$dist <- pmin(1, pmax(0, plotData$lower-1)) - pmin(1, pmax(0, 1-plotData$upper))
  plotData$cut <- as.character(cut(plotData$dist, breaks = seq(-1, 1, length.out = nCols), 
      include.lowest = TRUE, labels = seq_len(nCols)[-(nCols+1)/2]))
  plotData[plotData$dist == 0, "cut"] <- as.character((nCols+1)/2)
  
  # prepare fill legend
  synCalls <- c("none", "antagonism", "synergy") 
  plotData$synLabel <- factor(
      synCalls[1 + 1*(plotData$lower > 1) + 2*(plotData$upper < 1)],
      levels = synCalls
  )
  legendColors <- setNames(allColors[c((nCols+1)/2, nCols-1, 2)], nm = synCalls)
  # subset to only the colors that are present in the data
  legendColors <- legendColors[names(legendColors) %in% as.character(unique(plotData$synLabel))]

  # text to show
  plotData$label <- if (all(c("lower", "upper") %in% names(plotData))) {
        sprintf("%.2f\n(%.2f, %.2f)", plotData$tau, plotData$lower, plotData$upper)
      } else {
        sprintf("%.2f", plotData$tau)
      }
  
  # show doses on equidistant grid
  plotData$d1 <- factor(plotData$d1)
  plotData$d2 <- factor(plotData$d2)
    
  
  p <- ggplot(data = plotData, aes_string(x = "d1", y = "d2")) +
      geom_tile(aes_string(fill = "cut"), color = "grey") + 
      geom_text(aes_string(label = "label"), show.legend = FALSE, size = 3) +
      # invisible points, used only for labels
      geom_point(aes_string(color = "synLabel"), alpha = 0) +
      # round dose labels to digits
      scale_x_discrete(labels = format(as.numeric(levels(plotData$d1)), digits = digits)) + 
      scale_y_discrete(labels = format(as.numeric(levels(plotData$d2)), digits = digits)) +
      scale_fill_manual(breaks = seq_len(nCols), values = allColors, 
          guide = "none") + 
      scale_color_manual( # for a nicer legend
          values = setNames(1:3, nm = synCalls),
          limits = force,
          guide = guide_legend(title = "call:", 
              override.aes = list(alpha = 1, shape = 22, size = 8, color = "grey", 
                  fill = legendColors))
      ) + 
      theme_minimal() + 
      theme(
          panel.grid.major = element_blank(), 
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1)
      )
  
  p
  
}
