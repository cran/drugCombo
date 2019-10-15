#' Run nonparametric bootstrap on the interaction index model
#'
#' Function to run nonparametric bootstrap on the interaction index model.
#' It is usually called from \code{\link{getTauSurface}}.
#'  
#' @param fit A \code{HarbronFit} object returned by 
#' \code{\link{fitModel}}.
#' @param niter Number of bootstrap samples to use.
#' @param resample Resampling method for bootstrap. Either "all" (default) for 
#' resampling from all data, "mono" for separately resampling monotherapy and 
#' combination data, or "stratified" for resampling at each dose combination 
#' separately. Note that the latter method is not meaningful if there are no 
#' replicates in the data. 
#' @param seed Random seed to use for bootstrap
#' @param verbose Whether to show progress information.
#' @param ... Further arguments passed to the \code{\link{fitModel}} calls.
#' @return A matrix of interaction index estimates based on the bootstrap samples.
#' @importFrom stats coef
#' @author Maxim Nazarov
#' @export
getBootTaus <- function(fit, niter = 100, 
    resample = c("all", "mono", "stratified"), seed = NULL, verbose = FALSE, ...) {
  
  if (!is.null(seed))
    set.seed(seed)
  
  taus <- errors <- list()
  
  # original data - do we need an arg?
  dataCols <- c("d1", "d2", "effect")
  if (fit$tauSpec == "symbolic") 
    dataCols <- union(dataCols, all.vars(fit$tauFormula))
  
  data0 <- data.frame(mget(dataCols, envir = fit$m$getEnv(), inherits = FALSE))
  
  resample <- match.arg(resample)
  
  for (b in seq_len(niter)) {
    # resample
    dataB <- switch(resample,
        "all" = data0[sample(seq_len(nrow(data0)), size = nrow(data0), replace = TRUE), ],
        "mono" = { # resample mono & combi data separately
          monoIdx <- which(data0$d1 * data0$d2 == 0)
          otherIdx <- which(data0$d1 * data0$d2 != 0)
          rbind(
              data0[sample(monoIdx,  size = length(monoIdx),  replace = TRUE), ],
              data0[sample(otherIdx, size = length(otherIdx), replace = TRUE), ]
          )
        },
        "stratified" = {
          do.call(rbind, 
              lapply(split(data0, interaction(data0$d1, data0$d2)), function(d) {
                    d[sample(seq_len(nrow(d)), size = nrow(d), replace = TRUE), ]
                  })
          )
        }
    )
    rownames(dataB) <- NULL
    
    monoB <- try(fitMarginals(dataB, fixed = getConstraintsAsFixed(fit$mono)), 
        silent = TRUE)
    if (inherits(monoB, "try-error")) next
    
    bootSep <- suppressMessages(try({
              
              fitModel(dataB, monoB, 
                  model = fit$tauModel, 
                  tauFormula = fit$originalTauFormula, # use this, as the fit$tauFormula is modified if any taus are fixed
                  tauLog = fit$tauLog, 
                  tauStart = fit$tauStart,
                  stage = fit$stage,
                  fixed = c(fit$fixedMono, fit$fixedTau),
                  inactiveIn = fit$inactiveIn,
                  ...) # TODO: also keep and pass ... arguments from fitModel
              
            }, silent = TRUE))
    
    if (inherits(bootSep, "try-error"))
      errors[[b]] <- attr(bootSep, "condition")$message
    
    if (!inherits(bootSep, "try-error") && isTRUE(bootSep$convInfo$isConv)) {
    
      coefs <- coef(bootSep)
      # add fixed
      coefs <- c(coefs, bootSep$fixedTau)
      monoCoefNames <- c("h1", "h2", "b", "m1", "m2", "e1", "e2")
      # order by increasing N in tauN's
      tauCoefs <- grep("tau\\d+", names(coefs), value = TRUE)
      tauCoefs <- tauCoefs[order(as.numeric(gsub("[a-z]", "", tauCoefs)))]
      coefOrder <- c(intersect(names(coefs), monoCoefNames), tauCoefs)
      
      # check that number of coefficients match original
      if (bootSep$nTaus != fit$nTaus) {
        errors[[b]] <- "Non-conformable number of tau-parameters during bootstrap. Consider another resampling method."   
      } else {      
        taus[[b]] <- coefs[coefOrder]  # TODO: save more info?
      }
      
    }
    
    if (verbose)
      cat("iteration: ", b, " / ", niter, ", converged: ", 
          b-sum(vapply(taus[seq_len(b)], is.null, logical(1))), " / ", b, 
          "\r", sep = "")
    
  }
  
  if (length(taus) == 0 && length(errors) > 0)
    stop("None of the bootstrap iterations succeeded. The errors are:\n - ", 
        paste(errors[!sapply(errors, is.null)], collapse = "\n - "))
  
  notConverged <- vapply(taus, is.null, logical(1))
  if (sum(notConverged) > 0)
    warning("The model fit didn't converge for ", sum(notConverged), "/", niter, 
        " iterations during bootstraping", call. = FALSE, immediate. = TRUE)
# non-null
  ntaus <- do.call(rbind, taus[!notConverged])
  ntaus <- ntaus[, grepl("tau", colnames(ntaus)), drop = FALSE]
  
  list(tauB = ntaus, errors = errors)

}
