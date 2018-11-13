# calculation of d/D in the additivity/synergy equation
doseRatio <- function(response, d, h, b, m, e) {
  
  stopifnot(length(response) == length(d))
  D <- rep(0, length(d))
  
  lower <- min(b, m)
  upper <- max(b, m)
  
  # In case of different asymptotes, response can be outside of [lower, upper] 
  # range, this needs to be handled separately (otherwise formula gives NaN)
  validIdx <- response > lower & response < upper
  
  D[validIdx] <- invL4(response[validIdx], h, b, m, e)
  D[response <= lower] <- if (b < m) 0 else Inf
  D[response >= upper] <- if (b < m) Inf else 0
  
  dPart <- ifelse(d == 0, 0, d/D)
  dPart
}

# Model specification. Used inside \code{\link{fitModel}}.
# @author Maxim Nazarov, Nele Goeyvaerts based on C. Harbron's original code
tauModelFormula <- function(d1, d2, h1 = NA, h2 = NA, e1 = NA, e2 = NA, 
    m1 = NA, m2 = NA, b = NA, b2 = b, ..., gp = NULL,
    formula = ~1, tauSpec = NULL, tauLog = FALSE, inactiveIn = 0, 
    niterations = 50, eps = 1e-14) {
  
  n <- length(d1)

  tau <- list(...)
  
  # back-transform (exp)
  if (!is.null(tau) && tauLog) 
    tau <- lapply(tau, exp)
  
  mono <- (d1 == 0) | (d2 == 0)
  
  if (inactiveIn == 0)  # separated, as we need all arguments below to be defined
    if ((b < m1 && b2 > m2) || (b > m1 && b2 < m2)) 
      stop("Agonist/antagonist combination is not currently supported.")
  
  ## define taus
  vars <- all.vars(formula)
  # define formula type
  if (all(vars %in% c("d1", "d2"))) { 
    # type 1: lm ('symbolic')
    mm <- model.matrix(formula, data.frame(d1 = d1, d2 = d2))
    mm <- mm[!mono, colSums(mm[!mono, , drop = FALSE]) != 0, drop = FALSE]
    stopifnot(ncol(mm) == length(tau))
    tauComp <- if (length(tau)) mm %*% unlist(tau) else 1 # ~= tau[gp]
  } else if (tauSpec == "symbolic") {
    # parse custom formula, depending on other columns in the data (e.g. ray)
    mm <- model.matrix(formula)
    # here we exclude columns (=levels) that are actually not used       
    mm <- mm[!mono, colSums(mm[!mono, , drop = FALSE]) != 0, drop = FALSE]
    
    stopifnot(ncol(mm) == length(tau))
    tauComp <- if (length(tau)) mm %*% unlist(tau) else 1 # ~= tau[gp]
  } else {
    # type 2: nls ('literal')
    # here we need to foresee possible tauLog
    tauComp <- eval(formula[[2]],  
        envir  = data.frame(d1 = d1[!mono], d2 = d2[!mono]), 
        enclos = list2env(tau))
  }
  
  # avoid taking ratios 
  tauVector <- rep(1, n)
  tauVector[!mono] <- tauComp
  
  # define region and initialize for binary search
  U <- rep(max(b, b2, m1, m2, na.rm = TRUE), n)
  L <- G <- rep(min(b, b2, m1, m2, na.rm = TRUE), n)               
  dPart1 <- dPart2 <- rep(0, n)
  
  # perform binary search
  for(i in 1:niterations) {
    
    response <- (L + U)/2.0
    
    if (inactiveIn != 1)  # first compound is active
      dPart1 <- doseRatio(response, d1, h1, b, m1, e1) 
    
    if (inactiveIn != 2)   # second compound is active
      dPart2 <- doseRatio(response, d2, h2, b2, m2, e2)
    
    G <- dPart1 + dPart2 - tauVector
    
    G[abs(G) < eps] <- 0

    # find approximate response(y) by "binary search" 
    L[G == 0] <- U[G == 0] <- response[G == 0]   # found response value
    
    if ((inactiveIn == 0 && (b > m1 || b2 > m2)) || 
        (inactiveIn == 1 && b2 > m2) || 
        (inactiveIn == 2 && b > m1)) { # slope is negative  
      # increase response on next iteration by increasing L [as response = (L+U)/2]
      L[G < 0] <- response[G < 0]   
      # decrease response on next iteration by decreasing U
      U[G > 0] <- response[G > 0]    
    } else {  # otherwise, switch L and U 
      U[G < 0] <- response[G < 0]
      L[G > 0] <- response[G > 0]
    }
  }
  
  return(response)
}



defaultModels <- function() {
  list(
      "additive"    = ~ 0, 
      "uniform"     = ~ 1,
      "linear1"     = ~ log10(d1),
      "separate1"   = ~ 0 + factor(d1),
      "linear2"     = ~ log10(d2),
      "separate2"   = ~ 0 + factor(d2),
      "separate12"  = ~ 0 + factor(d1):factor(d2), 
      "zhao"        = ~ exp(tau1 + tau2*log(d1) + tau3*log(d2) + tau4*log(d1)*log(d2) + tau5*log(d1)^2 + tau6*log(d2)^2)
  )
}


#' Fit drug interaction index model according to Harbron's framework
#' 
#' @description This is the main function to fit an interaction index model to 
#' drug combination data based on the Loewe additivity model. The interaction 
#' index can be specified in a flexible way as a function of doses and other 
#' variables.
#' 
#' @param data A (long) data frame to fit the model to. Required columns are 
#' "d1", "d2" and "effect". Other variables are allowed and can be used in 
#' \code{tauFormula}.
#' @param mono An optional "MarginalFit" object, obtained from 
#' \code{\link{fitMarginals}}. 
#' @param model A pre-defined model to use for the interaction index tau. One of 
#' "additive", "uniform", "linear1", "separate1", "linear2", "separate2", 
#' "separate12" or "zhao". See details. 
#' @param tauFormula A formula to define the interaction index tau, using either 
#' 'literal' (as in \code{nls}) or 'symbolic' (as in \code{lm}) specification.
#' @param tauLog Whether to fit the model using log-transformed tau parameters.
#' This is mostly useful for "separate"-type tau models for better convergence.
#' Note that if TRUE, tau cannot be negative, which may be not approriate
#' for some models, such as "linear1" and "linear2".
#' Note that this affects the coefficient names in the result 
#' ("logtau1", "logtau2", ... instead of "tau1", "tau2", ...), so if 
#' \code{fixed} argument is used, this should be taken into account. 
#' @param tauStart Vector of starting values for tau parameters, either of 
#' length 1 or of the same length as the total number of tau parameters.
#' @param stage Whether to run a 1-stage or 2-stage estimation.
#' @param fixed Constraints on monotherapy and/or tau parameters as a vector 
#' of the form `name = value`, if NULL (default), taken from \code{mono}. Note 
#' that the tau parameters should be named "tau1", "tau2", ... if 
#' \code{tauLog = FALSE} (default), or "logtau1", "logtau2", ... if 
#' \code{tauLog = TRUE}.
#' @param inactiveIn which compound is inactive (1 or 2), or 0 (default) when 
#' both compounds are active.
#' @param verbose Whether to show extra information useful for debugging.
#' @param ... Further arguments passed to the \code{\link[minpack.lm]{nlsLM}} 
#' call. For example, \code{trace = TRUE} is useful to see the trace of the 
#' fitting process and may help identify issues with convergence.
#' @return Fitted object of class "HarbronFit" which is an \code{nls}-like 
#' object with extra elements.
#' 
#' @details
#' There are different ways to specify a model for the interaction index tau:
#' \itemize{
#' \item Using one of the pre-defined models as specified in the \code{model} 
#' argument:
#' 	\itemize{
#' 		\item "additive", for additivity model,
#' 		\item "uniform", one overall value for tau,
#' 		\item "linear1", linear dependency on log10 dose of the first compound,
#' 		\item "linear2", linear dependency on log10 dose of the second compound,
#' 		\item "separate1", different tau for each dose of the first compound,
#' 		\item "separate2", different tau for each dose of the second compound, 
#' 		\item "separate12", different tau for each combination of doses of the two compounds,
#' 		\item "zhao", quadratic response surface model following Zhao et al 2012.
#' 	}
#' \item Using a literal or symbolic formula. Note that for the monotherapies, 
#' tau is assumed to be equal to 1. Therefore, continuous models may entail 
#' discontinuities in the interaction index when d1 and d2 approach 0.
#' }
#' 
#' @importFrom minpack.lm nlsLM
#' @importFrom BIGL fitMarginals
#' @importFrom stats setNames as.formula
#' @author Maxim Nazarov
#' @examples
#' \donttest{
#' data("checkerboardData", package = "drugCombo")
#' data1 <- checkerboardData[checkerboardData$exp == 1, ]
#' mono1 <- fitMarginals(data1, fixed = c(b = 1))
#' # all three ways below are equivalent
#' fitLin1 <- fitModel(data = data1, mono = mono1, model = "linear1")
#' fitLin1b <- fitModel(data1, mono1, tauFormula = ~ log10(d1))
#' fitLin1c <- fitModel(data1, mono1, tauFormula = ~ tau1+tau2*log10(d1))
#' }
#' @export
fitModel <- function(data, mono = NULL, 
    model = NULL, tauFormula = NULL, tauLog = FALSE, tauStart = 1 * (!tauLog), 
    stage = 1, fixed = NULL, inactiveIn = 0, verbose = FALSE, 
    ...) {
  
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
  
  ## monotherapies
  if (is.null(mono)) {  # default mono-fitting
    message("Running default monotherapy fitting. Note that if constraints are needed, run mono fitting manually first")
    mono <- fitMarginals(data)
  } else if (!inherits(mono, "MarginalFit")) {
    # TODO: or do we still support other ways of specifying mono?
    stop("Please provide a 'marginalFit' object as `mono` parameter, or use `mono = NULL` to apply default monotherapy fitting.")
  }

  if (inherits(mono, "MarginalFit")) {
    monoCoefs <- mono$coef
  } else {
    stop("Error while fitting monotherapy, please try fitting it manually and pass as `mono` argument.")
  }
  
  monoNamesConstrained <- NULL
  
  if (stage == 1) {
    
    monoStart <- monoCoefs
    monoPars <- setNames(nm = c("h1", "h2", "e1", "e2", "b", "m1", "m2")) 
    
    ## use constraints from mono
    # parse constraints into fixed
    if (is.null(fixed) && !is.null(mono$model$constraints)) {
      fixed <- getConstraintsAsFixed(mono)
    }
    
    if (!is.null(fixed)) {
      # TODO: support BIGL linear constraints?
      if (length(fixed) > 0) {

        monoNamesConstrained <- intersect(names(fixed), names(monoPars))
        monoPars[monoNamesConstrained] <- fixed[monoNamesConstrained]
        monoStart <- monoStart[!(names(monoStart) %in% monoNamesConstrained)]
      
      }
    }
  } else if (stage == 2) {
    
    monoStart <- c()
    monoPars <- monoCoefs
    
    if (!is.null(fixed) && any(names(monoCoefs) %in% names(fixed)))
      stop("For stage = 2, monotherapy constraints should be defined inside the `mono` object. ", 
          "Please run monotherapy fitting separately and pass as `mono` parameter.")
        
  } else {
    stop("`stage` must be either 1 or 2.")
  }
  
  # if `inactiveIn` is specified, exclude corresponding mono coefs
  if (length(inactiveIn) == 1 && (inactiveIn == 1 || inactiveIn == 2)) {
    monoPars <- monoPars[!grepl(inactiveIn, names(monoPars))]
    monoStart <- monoStart[!grepl(inactiveIn, names(monoStart))]
  } else { # silently ignore other cases
    inactiveIn <- 0
  }
  
  monoPars <- c(d1 = "d1", d2 = "d2", monoPars)
  
  monoParsString <- makeString(monoPars) #paste(names(monoPars), "=", monoPars, collapse = ", ")
  
  model <- match.arg(model, names(defaultModels()))
  
  tauFormula <- if (is.null(tauFormula)) defaultModels()[[model]] else tauFormula
  
  ## determine type of formula supplied
  vars <- all.vars(tauFormula)
  
  tauSpec <- "symbolic"
  
  if (all(vars %in% names(data))) {
    # type 1: lm
    if (verbose)
      message("Assuming symbolic formula")
    mm <- model.matrix(tauFormula, data)
    monoIdx <- (data$d1 * data$d2 == 0)
    mm <- mm[!monoIdx, colSums(mm[!monoIdx, , drop = FALSE]) != 0, drop = FALSE]
    nTaus <- ncol(mm)
    
    tauNames <- if (nTaus > 0) setNames(nm = paste0(if (tauLog) "log", "tau", seq_len(nTaus))) else character(0)
    
  } else {
    # type 2: nls
    tauSpec <- "literal"
    nTaus <- sum(grepl("tau\\d+", vars)) #NB: not robust...
    if (verbose)
      message("Assuming literal formula with ", nTaus, " taus")
    #model <- NULL
    
    tauNames <- setNames(nm = grep("tau\\d+", vars, value = TRUE))
    
  }
  
  # take possible tauX constraints into account
  tauNamesFixed <- intersect(names(fixed), names(tauNames))
  tauNames[tauNamesFixed] <- fixed[tauNamesFixed]
  # FIXME: ugly?
  tauNamesFree <- tauNames[!(names(tauNames) %in% tauNamesFixed)]
  nFreeTaus <- length(tauNamesFree)
  
  # Here we put tau constraints into the formula, and make sure it is still 
  # a formula, this is similar to calling 
  # pryr::substitute_q(tauFormula, as.list(fixed))
  # TODO: not sure it works correctly for symbolic formulas
  origFormula <- tauFormula # keep original formula for bootstrap
  tauFormula <- as.formula(
      do.call(substitute, list(tauFormula, as.list(fixed[tauNamesFixed]))), 
      env = attr(tauFormula, ".Environment")
  )
  
  modelFormula <- as.formula(paste("effect ~ tauModelFormula(",
          monoParsString, ",", 
          if (nTaus > 0) paste(makeString(tauNames), ",") else NULL, 
          "formula = ", deparse(tauFormula, width.cutoff = 500),
          ", tauSpec = ", shQuote(tauSpec),
          ", tauLog = ", tauLog,
          ", inactiveIn = ", inactiveIn,
          ")"))
  
  tauStart <- if (nFreeTaus > 0) {
    if (length(tauStart) == 1)
      setNames(rep_len(tauStart, nFreeTaus), names(tauNamesFree))
    else if (length(tauStart) == nFreeTaus)
      setNames(tauStart, names(tauNamesFree))
    else
      stop("`tauStart` must be a vector of length 1 or length equal to the number of tau parameters (", nFreeTaus, ")")
  } else NULL
  
  tauStart <- tauStart[!(names(tauStart) %in% tauNamesFixed)]

  # inform user
  if (length(c(monoNamesConstrained, tauNamesFixed)) > 0)
    message("Using constraints: ", 
        makeString(fixed[c(monoNamesConstrained, tauNamesFixed)])
    )
  
  startingValues <- c(monoStart, tauStart)
  if (verbose) {
    cat("nls formula:\n  ")
    print(modelFormula)
    cat("starting values:\n")
    print(startingValues)
  }    
  fit <- nlsLM(formula = modelFormula, data = data, start = startingValues, ...)
  
  # TODO: better structure?
  fit$tauFormula <- tauFormula
  fit$originalTauFormula <- origFormula
  fit$tauSpec <- tauSpec
  fit$tauLog <- tauLog
  fit$tauStart <- if (length(unique(tauStart)) == 1) unique(tauStart) else tauStart
  fit$tauModel <- model
  fit$stage <- stage
  fit$mono <- mono
  # if 'fixed' has some names that are not used, we remove them
  fit$fixedMono <- fixed[monoNamesConstrained]  
  fit$fixedTau  <- fixed[tauNamesFixed]
  fit$inactiveIn <- inactiveIn
  fit$nTaus <- nTaus
  
  class(fit) <- append("HarbronFit", class(fit))
  fit
}

#' @importFrom BIGL fitMarginals
#' @export
BIGL::fitMarginals


getConstraintsAsFixed <- function(mono) {
  if (is.null(mono$model$constraints))
    return(NULL)
  
  parNames <- apply(mono$model$constraints$matrix, 1, 
      function(row) { 
        if (setequal(unique(row), 0:1)) 
          mono$model$order[as.logical(row)] 
        else NULL
      })
  validNames <- !vapply(parNames, is.null, logical(1))
  if (any(!validNames))
    warning("General linear constraints are ignored, only simple ones will be used.")
  
  fixedVector <- mono$model$constraints$vector[validNames]
  names(fixedVector) <- parNames[validNames]
  
  fixedVector
}

makeString <- function(x) { 
  if (is.null(names(x))) 
    names(x) <- x
  paste(names(x), "=", x, collapse = ", ")
}