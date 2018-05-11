#' Compute interaction index surface and confidence intervals 
#' 
#' Computes estimates and confidence intervals for the interaction surface for 
#' all dose combinations.
#' @param data Dose combinations to compute interaction index for. If NULL 
#' (default), taken from the \code{fit} object. 
#' @param addCI Whether confidence intervals need to be computed.
#' @param method Which method to use to calculate confidence intervals: 
#' "default" for Wald-type or "boot" for nonparametric bootstrap.
#' @param level The confidence level required for the confidence intervals 
#' (default is 0.95).
#' @inheritParams getBootTaus
#' @param ... Further parameters that are passed to \code{\link{getBootTaus}}.
#' @return An object of class "tauSurface" which is essentially a list with 
#' the following components: data frame with interaction index (tau) estimates,
#' standard errors and pointwise confidence intervals, formulas for computing 
#' tau at any given dose (only for models with continuous functions used to 
#' define tau), and details on the tau specification from the \code{fit}. In 
#' addition, if the "boot" method was used, all the bootstrap estimates are 
#' returned and can be accessed with \code{\link{bootstrapCoefs}}.
#' @importFrom stats model.matrix coef quantile setNames vcov confint.default
#' @author Maxim Nazarov
#' @seealso \code{\link{plot.tauSurface}}, \code{\link{contour.tauSurface}} for
#' visual representation of the tau surface.
#' @export
getTauSurface <- function(fit, data = NULL, addCI = TRUE, 
    method = c("default", "boot"), level = 0.95, niter = 100, 
    resample = c("all", "mono", "stratified"), seed = NULL, ...) {
  
  # checks
  if (!inherits(fit, "HarbronFit")) {
    stop("Unrecognized `fit` object provided. Function accepts only results from `fitModel`.")
  }
  if(!"tauSpec" %in% names(fit))
    stop("Can't get the surface since tau specification is not defined. Please run `fitModel` first.")
  
  formulaCols <- c("d1", "d2")
  if (fit$tauSpec == "symbolic") 
    formulaCols <- union(formulaCols, all.vars(fit$tauFormula))
  
  # define data
  if (is.null(data))
    data <- data.frame(mget(formulaCols, envir = fit$m$getEnv(), inherits = FALSE),
        stringsAsFactors = FALSE)

  data <- droplevels(unique(data[data$d1*data$d2 != 0, formulaCols])) #c("d1", "d2")])
  rownames(data) <- NULL
  
  cis <- NULL
  surfaceFormula <- NULL
  ciFormula <- NULL
  
  method <- match.arg(method)
  
  tauB <- NULL
  errorsB <- NULL
  if (method != "boot")
    niter <- NULL
  
  if (fit$tauSpec == "symbolic") {
    
    tres <- tauDeltaMethod(fit = fit, data = data, level = level)
    tau <- tres[["tau"]] 
        
    # add formula for continuous plotting
    # currently do this if we have only d1 and/or d2 in the formula and they are not factors
    # FIXME: is there a better way to do this?
    factorsInFormula <- grepl("factor", deparse(attr(terms(fit$tauFormula), "variables")))
    continuous <- isTRUE(setequal(formulaCols, c("d1", "d2")) && !factorsInFormula)
    
    if (continuous)
      surfaceFormula <- function(d1, d2) {

        tauIdx <- grepl("tau", names(coef(fit)))
        tauEstimates <- coef(fit)[tauIdx]
        # add fixed tau values & reorder
        tauEstimates <- c(tauEstimates, fit$fixedTau)
        tauOrder <- order(as.numeric(gsub("[a-z]", "", names(tauEstimates))))
        tauEstimates <- tauEstimates[tauOrder]
        # tauLog
        if (fit$tauLog) {
          tauEstimates <- exp(tauEstimates)
        }
        
        if (length(tauEstimates))  
          as.vector(
              model.matrix(
                  fit$tauFormula,
                  as.data.frame(Filter(Negate(is.null), list(d1 = d1, d2 = d2))) 
                  # ^ this is to remove NULLs to avoid error when converting to data frame
              )  %*% tauEstimates
          )
        else 
          1
      }
    
    if (addCI) {
      if (method == "boot") {
        
        bootRes <- getBootTaus(fit, niter = niter, resample = resample, seed = seed, ...)
        tauB <- bootRes$tauB
        errorsB <- bootRes$errors
    
        # tauLog
        if (fit$tauLog) {
          tauB <- exp(tauB)
        }
        
        tauBoots <- model.matrix(fit$tauFormula, data) %*% t(tauB)
        
        cis <- as.data.frame(t(apply(tauBoots, 1, function(row) {
                      setNames(c(quantile(row, probs = c((1 - level)/2, (1 + level)/2))), 
                          c("lower", "upper"))
                    })))
        
        # add formula for continuous plotting
        # currently can do this if we have only d1/d2 in the formula
        if (continuous)
          ciFormula <- function(side = c("lower", "upper"), level = 0.95) {
            side <- match.arg(side)
            function(d1, d2) {
              tauBoots <- apply(tauB, 1, function(row) {
                    model.matrix(fit$tauFormula, Filter(Negate(is.null), list(d1 = d1, d2 = d2))) %*% row
                  })
              
              prob <- switch(side, "lower" = (1 - level)/2, "upper" = (1 + level)/2)
              apply(tauBoots, 1, quantile, probs = prob, names = FALSE)
            }
          }
        
      
      } else {
        
        cis <- tres[, c("se.est", "lower", "upper")]
        
        # add formula for continuous plotting
        if (continuous)
          
          ciFormula <- function(side = c("lower", "upper"), level = 0.95) {
            side <- match.arg(side)
            function(d1, d2) 
              tauDeltaMethod(fit, data = Filter(Negate(is.null), list(d1 = d1, d2 = d2)), level = level)[[side]]
          } 
        
      }
    }
  
  } else if (fit$tauSpec == "literal") {
    
    tres <- tauDeltaMethod(fit = fit, data = data, level = level)
    tau <- tres[["tau"]]
    
    surfaceFormula <- function(d1, d2) { 
      
      tauNames <- grep("tau", names(coef(fit)), value = TRUE)
      tauEstimates <- mget(tauNames, envir = fit$m$getEnv(), inherits = FALSE)
      if (fit$tauLog) 
        tauEstimates <- lapply(tauEstimates, exp)
      
      eval(expr = do.call(substitute, list(fit$tauFormula[[2]], tauEstimates)),
          envir = list(d1 = d1, d2 = d2)
      )
    }
    
    if (addCI) {
      
      if (method == "boot") {
        
        bootRes <- getBootTaus(fit, niter = niter, resample = resample, seed = seed, ...)
        tauB <- bootRes$tauB
        errorsB <- bootRes$errors
        
        # tauLog
        if (fit$tauLog) {
          tauB <- exp(tauB)
        }
        
        tauBoots <- apply(tauB, 1, function(row) {
              eval(fit$tauFormula[[2]], data, enclos = list2env(as.list(row)))
            })
        
        cis <- as.data.frame(t(apply(tauBoots, 1, function(row) {
                      setNames(c(quantile(row, probs = c((1 - level)/2, (1 + level)/2))), 
                          c("lower", "upper"))
                    })))
        
        # add formula for continuous plotting
        ciFormula <- function(side = c("lower", "upper"), level = 0.95) {
          side <- match.arg(side)
          function(d1, d2) {
            tauBoots <- apply(tauB, 1, function(row) {
                  eval(fit$tauFormula[[2]], list(d1 = d1, d2 = d2), enclos = list2env(as.list(row)))
                })
            
            prob <- switch(side, "lower" = (1 - level)/2, "upper" = (1 + level)/2)
            apply(tauBoots, 1, quantile, probs = prob, names = FALSE)
          }
        }
      
      } else {
        cis <- tres[, c("se.est", "lower", "upper")]
        
        # add formula for continuous plotting
        ciFormula <- function(side = c("lower", "upper"), level = 0.95) {
          side <- match.arg(side)
          function(d1, d2) 
            tauDeltaMethod(fit, data = list(d1 = d1, d2 = d2), level = level)[[side]]
        }
        
      }
    }
  } else {
    stop("Unknown tau specification type.")
  }
  
  if (fit$stage == 2 && addCI && method != "boot") {
    warning("Confidence intervals computed don't take into account error propagation from the 1st stage (mono). Please use method = 'boot' for this.")
  }
  
  res <- cbind(data, tau = c(tau))
  if (!is.null(cis))
    res <- cbind(res, cis)
  
  out <- list(data = res[], tauSpec = fit$tauSpec, tauFormula = fit$tauFormula, 
      surfaceFormula = surfaceFormula, ciFormula = ciFormula, method = method)
  
  if (method == "boot")
    out <- c(out, 
        list(niter = niter, bootstrapCoefs = tauB, bootstrapErrors = errorsB))
  
  class(out) <- append("tauSurface", class(out))
  out
  
}


#' Show estimated model parameters from all bootstrap iterations
#' 
#' @param tauSurface A \code{tauSurface} object returned by 
#' \code{\link{getTauSurface}}.
#' @return matrix of parameter estimates from the bootstrap iterations as 
#' returned by \code{\link{getBootTaus}}
#' 
#' @author Maxim Nazarov
#' @export
bootstrapCoefs <- function(tauSurface) {
  
  if (!inherits(tauSurface, "tauSurface") || tauSurface$method != "boot" || is.null(tauSurface$bootstrapCoefs)) {
    warning("No bootstrap coefficients available")
    return(NULL)
  }
 
  tauSurface$bootstrapCoefs
  
}

#' Unique method for "tauSurface" objects
#' 
#' "unique" method to extract unique interaction index estimates from a 
#' "tauSurface" object in a tabular format.
#' 
#' @param x Output of \code{\link{getTauSurface}}.
#' @param ... Further arguments, currently not used.
#' @export
unique.tauSurface <- function(x, ...) {
  if (x$tauSpec != "grouping") {
    if (length(notUsed <- setdiff(c("d1", "d2"), all.vars(x$tauFormula))))
      toPrint <- unique(x$data[, setdiff(names(x$data), notUsed)])
    else
      toPrint <- x$data
  } else {
    toPrint <- x$data[!duplicated(eval(x$tauFormula, x$data)), ]
  }
  rownames(toPrint) <- NULL
  
  toPrint
}

#' Print method for "tauSurface" objects
#' 
#' @param x Output of \code{\link{getTauSurface}}.
#' @param ... Further arguments, currently not used.
#' @export
print.tauSurface <- function(x,  ...) {
  
  toPrint <- x$data
  
  rownames(toPrint) <- NULL
  
  cat("A \"tauSurface\" object")
  if ("lower" %in% names(x$data))
    cat(", conf.intervals method:", deparse(x$method), "\n")
  print(toPrint)
  invisible(toPrint)
}


# Calculate interaction index estimates and Wald-type confidence intervals 
# 
# The delta method is used for literal tau formulas.
# 
# @param fit Object of class "HarbronFit".
# @param data Data to calculate confidence intervals for (list or data frame).
# @param level Confidence level.
#' @importFrom Deriv Deriv
#' @importFrom stats coef vcov qnorm
# @return data frame
# @keywords internal
# @author Maxim Nazarov
tauDeltaMethod <- function(fit, data, level = 0.95) { 
  
  # whether fitting was done on log-transformed taus
  tauLog <- isTRUE(fit$tauLog) 

  tauNames <- grep("tau", names(coef(fit)), value = TRUE)
  
  tauEstimates <- as.list(coef(fit)[tauNames])
  # add fixed tau values
  tauEstimates <- c(tauEstimates, fit$fixedTau)
  if (tauLog) 
    tauEstimates <- lapply(tauEstimates, exp)
  
  tauVar <- vcov(fit)[tauNames, tauNames, drop = FALSE]
  
  if (fit$tauSpec == "literal") {
    
    tau <- eval(fit$tauFormula[[2]], envir = data, enclos = list2env(tauEstimates))

    # if tauLog = TRUE, we substitute all 'tauN' into 'exp(tauN)' in the formula
    dFun <- if (tauLog) {
          do.call(
              substitute, 
              list(fit$tauFormula, 
                  setNames(nm = tauNames, lapply(tauNames, function(x) substitute(exp(x), list(x = as.name(x)))))
              )
          )
        } else fit$tauFormula
    
    # gradient
    gd <- as.matrix(eval(envir = data, enclos = fit$m$getEnv(),
        expr = Deriv::Deriv(f = dFun, x = tauNames, combine = "cbind")
    ))
    
  } else if (fit$tauSpec == "symbolic") {
    
    # handle fixed tau
    fixedTauNames <- names(fit$fixedTau)
    # add 0 variance for fixed tau values
    if (length(fixedTauNames)) {
      # adjust tauVar
      extraColMat <- matrix(0, nrow = nrow(tauVar), ncol = length(fixedTauNames), 
          dimnames = list(c(), fixedTauNames))
      tauVar <- cbind(tauVar, extraColMat)
      
      extraRowMat <- matrix(0, nrow = length(fixedTauNames), ncol = ncol(tauVar), 
          dimnames = list(fixedTauNames, c()))
      tauVar <- rbind(tauVar, extraRowMat)
      
      # order tauN's by N
      tauOrder <- order(as.numeric(gsub("[a-z]", "", colnames(tauVar))))
      
      tauEstimates <- tauEstimates[tauOrder]
      tauVar <- tauVar[tauOrder, tauOrder, drop = FALSE]
    }
    
    gd <- model.matrix(fit$tauFormula, as.data.frame(data))
    
    tau <- if (length(tauEstimates)) gd %*% unlist(tauEstimates) else 1
    
    # FIXME: understand this fully
    if (tauLog)
      tauVar <- tauVar * unlist(tauEstimates) %*% t(unlist(tauEstimates))
    
  } else 
    stop("Unsupported tau specification for CI calculation.")

  se.est <- sqrt(diag(gd %*% tauVar %*% t(gd)))
  
  lower <- tau - qnorm((1 + level)/2) * se.est
  upper <- tau + qnorm((1 + level)/2) * se.est
  
  data.frame(tau, se.est, lower, upper, row.names = NULL)
  
}
