###############################################################################
#                        Package: mnlogit                                     #
#                                                                             #
# Multinomial logit, maximum likelihood estimation by Newton-Raphson method   #
#                                                                             #
#               Scientific Computing Group, Sentrana Inc.                     #
###############################################################################

###############################################################################
#                 Main user function                                          #
# Args:  							              #
#   formula     - same format as mlogit. See help(formula)                    #
#   data        - input data (as a data.frame object) in "long" format        #
#   choiceVar   - the data column containing alternative names                #
#   maxiter     - maximum number of Newton-Raphson iterations to run          #
#   ftol        - function tolerance.                                         #
#                 Difference of two consecutive function evaluation           #
#                 Criteria of terminating Newton's Iterative process          #
#   gtol        - gradient norm tolerance.                                    #
#   weights     - an optional vector of positive frequency weights.           #
#   ncores      - number of processors allowed to use                         #
#   na.rm       - if FALSE then stop(), else remove rows of data with NA      #
#   print.level - increase from 0 to progressively print more runing info     #
#   linDepTol   - tolerance with which linear dep among cols is detected      #
#   start       - initial vector of coefficients                              #
#   alt.subset  - subset of alternatives to perform estimation on             #
#   ...         - currently unused                                            #
#                                                                             #
# Output:                                                                     #
#   mnlogit object                                                            #
#                                                                             #
# Note:                                                                       #
#   'ftol', 'gtol' & 'maxiter' specify Newton-Raphson termination criteria    #
###############################################################################

modified_mnlogit <- function (formula, data, choiceVar = NULL, maxiter = 50, ftol = 1e-8,
                              gtol = 1e-8, weights = NULL, ncores = 1, na.rm = TRUE, print.level = 0,
                              linDepTol = 1e-8, start = NULL, alt.subset = NULL, ...)
{
  #browser()
  # startTime <- proc.time()[3]
  # initcall <- match.call()
  # if (!is.data.frame(data))
  #   stop("data must be a data.frame in long format or a mlogit.data object")
  # if (ncores < 1) {
  #   ncores <- 1
  #   warning("Setting ncores equal to: 1")
  # }
  # if (!is.null(choiceVar) && is.factor(data[[choiceVar]])) {
  #   warning(paste("Column", choiceVar, "in data will NOT be treated as",
  #                 "factor, but as character string!"))
  # }
  # if (is.null(choiceVar) && !any(class(data) == "mlogit.data"))
  #   stop("Arg data MUST be a mlogit.data object when arg choiceVar = NULL")
  # if (is.null(choiceVar)) {
  #   choiceVar <- "_Alt_Indx_"
  #   data[[choiceVar]] <- attr(data, "index")$alt
  # }
  formula <- parseFormula(f = formula)
  response <- attr(formula, "response")
  interceptOn <- attr(formula, "Intercept")
  # csvChVar <- attr(formula, "csvChCoeff")
  # indspVar <- attr(formula, "indSpVar")
  covariates <- attr(formula, "indSpVar")
  # csvGenVar <- attr(formula, "csvGenCoeff")
  # covariates <- c(csvChVar, indspVar, csvGenVar)
  varNames <- attr(formula, "varNames")
  # if (is.null(covariates) && !interceptOn)
  #   stop("Error! Predictor variable(s) must be specified")
  # if (is.null(response))
  #   stop("Error! Alternative variable must be specified")
  # if (!is.null(alt.subset)) {
  #   if (sum(unique(data[[choiceVar]]) %in% alt.subset) <
  #       2)
  #     stop("Error! Atleast 2 alternatives in data must be in alt.subset")
  #   keepRows <- data[[choiceVar]] %in% alt.subset
  #   if (sum(keepRows) <= 0)
  #     stop("Error! No altrnative in 'alt.subset' is in data.")
  #   data <- data[keepRows, , drop = FALSE]
  # }
  choice.set <- c("A", "C", "G", "T")
  # choice.set <- unique(data[[choiceVar]])
  K <- length(choice.set)
  if (nrow(data)%%K)
    stop("Mismatch between number of rows in data and number of choices.")
  N <- nrow(data)/K
  if (!is.null(weights) && length(weights) != N)
    stop("Length of 'weights' arg must match number of observations in data.")
  #if (!is.null(weights) && !all(weights > 0))
    #stop("All entries in 'weights' must be strictly positive.")
  #if (!is.null(weights))
  #weights <- weights * N/sum(weights)
  # data <- data[, c(varNames, choiceVar)]
  # na.rows <- c()
  # for (col in 1:ncol(data)) na.rows <- union(na.rows, which(is.na(data[[col]])))
  # Ndropped <- 0
  # if (length(na.rows) > 0) {
  #   if (!na.rm)
  #     stop("NA present in input data.frame with na.rm = FALSE.")
  #   keepRows <- rep(TRUE, nrow(data))
  #   keepRows[na.rows] <- FALSE
  #   for (i in 1:N) {
  #     if (!all(keepRows[((i - 1) * K + 1):(i * K)]))
  #       keepRows[((i - 1) * K + 1):(i * K)] <- FALSE
  #   }
  #   data <- data[keepRows, , drop = FALSE]
  #   if (!is.null(weights)) {
  #     weights <- weights[keepRows[seq(1, N * K, K)]]
  #   }
  #   N <- nrow(data)/K
  #   Ndropped <- (length(keepRows) - sum(keepRows))/K
  # }
  # if (print.level && Ndropped > 0)
  #   cat(paste("Num of dropped observations (due to NA)  =",
  #             Ndropped, "\n"))
  ind <- sort_ind(data[[choiceVar]], nrow(data)) + 1 #sort_ind only sort ATCG
  data <- data[ind[, 1], ]
  # choice.set <- unique(data[[choiceVar]])
  respVec <- data[[attr(formula, "response")]]
  # if (is.factor(respVec))
  #   respVec <- droplevels(respVec)
  # respVec <- as.numeric(respVec)
  # min.respVec <- min(respVec)
  # spread <- max(respVec) - min.respVec
  # if (spread != 1) {
  #   stop(paste("Response variable", attr(formula, "response"),
  #              "must be a factor with exactly 2 levels."))
  # }
  # respVec <- respVec - min.respVec
  # freq.choices <- colSums(matrix(respVec, nrow = N, ncol = K))/N
  # loFreq <- min(freq.choices)
  # loChoice <- choice.set[which(loFreq == freq.choices)]
  # names(freq.choices) <- choice.set
  # if (loFreq < 1e-07) {
  #   cat("Frequencies of alternatives in input data:\n")
  #   print(prop.table(freq.choices), digits = 4)
  #   stop(paste("Frequency, in response, of choice:", loChoice,
  #              "< 1e-7."))
  # }

  X <- formDesignMat_r(data = data, N = N, formula = formula, varVec = attr(formula, "indSpVar"),
                       includeIntercept = attr(formula, "Intercept"))
  #X <- formDesignMat(dat = data, N = N)
  #colnames(X) <- c("(Intercept)", "read_pos", "ref_pos", "qua", "hap_nucC", "hap_nucG", "hap_nucT",
                   #"qua:hap_nucC", "qua:hap_nucG", "qua:hap_nucT")
  # X <- if (!is.null(X))
  # X <- X[1:N, , drop = FALSE]
  # Y <- formDesignMat(varVec = attr(formula, "csvChCoeff"),
  #                    includeIntercept = FALSE)
  # Z <- formDesignMat(varVec = attr(formula, "csvGenCoeff"),
  #                    includeIntercept = FALSE)
  # badColsList <- list(indSpVar = NULL, csvChCoeff = NULL, csvGenCoeff = NULL)
  # badColsList$indSpVar <- getNullSpaceCols(X, tol = linDepTol)
  badColsList <- getNullSpaceCols(X, tol = linDepTol)
  # for (i in 1:K) {
  #   init <- (i - 1) * N + 1
  #   fin <- i * N
  #   badColsList$csvChCoeff <- union(badColsList$csvChCoeff,
  #                                   getNullSpaceCols(Y[init:fin, , drop = FALSE], tol = linDepTol))
  # }
  # badColsList$csvGenCoeff <- getNullSpaceCols(Z, tol = linDepTol)
  # badVarsList <- list()
  # badVarsList$indSpVar <- colnames(X[, badColsList$indSpVar,
  #                                    drop = FALSE])
  badVarsList <- colnames(X[, badColsList, drop = FALSE])
  # badVarsList$csvChCoeff <- colnames(Y[, badColsList$csvChCoeff,
  #                                      drop = FALSE])
  # badVarsList$csvGenCoeff <- colnames(Z[, badColsList$csvGenCoeff,
  #                                       drop = FALSE])
  badCoeffNames <- makeCoeffNames(badVarsList, choice.set)
  #if (!is.null(X))
  if (!is.null(badColsList))
    X <- X[, setdiff(1:ncol(X), badColsList), drop = FALSE]
  # if (!is.null(Y))
  #   Y <- Y[, setdiff(1:ncol(Y), badColsList$csvChCoeff),
  #          drop = FALSE]
  # if (!is.null(Z))
  #   Z <- Z[, setdiff(1:ncol(Z), badColsList$csvGenCoeff),
  #          drop = FALSE]
  # varNamesList <- list()
  # varNamesList$indSpVar <- colnames(X)
  varNamesList <- colnames(X)
  # varNamesList$csvChCoeff <- colnames(Y)
  # varNamesList$csvGenCoeff <- colnames(Z)
  coeffNames <- makeCoeffNames(varNamesList, choice.set)
  if (!is.null(start))
    names(start) <- coeffNames
    #start[coeffNames] <- start
  baseChoiceName <- choice.set[1]
  # if (!is.null(Z)) {
  #   for (ch_k in 2:K) {
  #     Z[((ch_k - 1) * N + 1):(ch_k * N), ] <- Z[((ch_k -
  #                                                   1) * N + 1):(ch_k * N), , drop = FALSE] - Z[1:N,
  #                                                                                               , drop = FALSE]
  #   }
  # }
  # Z <- Z[(N + 1):(K * N), , drop = FALSE]
  respVec <- respVec[(N + 1):(K * N)]
  # t1 <- proc.time()[3]
  # gc()
  # prep.time <- t1 - startTime
  # if (print.level > 1) {
  #   cat(paste0("Base alternative is: ", baseChoiceName))
  #   cat(paste0("\nPreprocessing data for estimation took ",
  #              round(prep.time, 3), " sec.\n"))
  # }
  result <- newtonRaphson(respVec, X, NULL, NULL, K, maxiter, gtol,
                          ftol, ncores, print.level, coeffNames, weights = weights,
                          start = start)
  # result$est.stats$prepTimeSecs <- prep.time
  colnames(result$hessMat) <- coeffNames
  # rownames(result$hessMat) <- coeffNames
  names(result$grad) <- coeffNames
  od <- reordering(varNamesList, choice.set)
  coeffNames <- makeCoeffNames(varNamesList, choice.set)
  coefficients <- c(result$coeff, if (is.null(badCoeffNames)) NULL else rep(NA,
                                                                            length(badCoeffNames)))
  names(coefficients) <- c(coeffNames, badCoeffNames[reordering(badVarsList,
                                                                choice.set)])
  reordered_coeff <- c(result$coeff[od], if (is.null(badCoeffNames)) NULL else rep(NA,
                                                                                   length(badCoeffNames)))
  names(reordered_coeff) <- c(coeffNames[od], badCoeffNames[reordering(badVarsList,
                                                                       choice.set)])
  colnames(result$probability) <- choice.set
  if (maxiter > 0)
    colnames(result$residual) <- choice.set
  # result$model.size$intercept <- interceptOn
  attributes(formula) <- NULL
  logLik <- structure(-result$loglikelihood, class = "logLik")
  AIC <- 2 * ((K - 1)*ncol(X) + result$loglikelihood)
  # index <- data.frame(chid = rep(1:result$model.size$N, result$model.size$K),
  #                     alt = data[[choiceVar]])
  # attr(data, "index") <- index
  fit <- structure(list(coefficients = coefficients, logLik = logLik,
                        gradient = -result$grad, hessian = result$hessMat, est.stat = result$est.stats,
                        fitted.values = 1 - attr(result$residual, "outcome"),
                        probabilities = result$probability, residuals = result$residual, AIC = AIC), class = "mnlogit")
  # if (print.level)
  #   cat(paste0("\nTotal time spent in mnlogit = ", round(proc.time()[3] -
  #                                                          startTime, 3), " seconds.\n"))
  return(fit)
}

# Makes names of model coefficients
# Sets ordering of the Hessian and gradient rows.
# Ensures order is in confirmation with design matrices
makeCoeffNames <- function (varNames, choices)
{
  if (length(varNames) == 0) return(NULL)
  choices <- as.vector(choices)
  # coeffName <- c(outer(varNames$indSpVar, choices[-1], paste, sep=":"),
  #                outer(varNames$csvChCoeff, choices, paste, sep=":"),
  #                varNames$csvGenCoeff)
  coeffName <- c(outer(varNames, choices[-1], paste, sep=":"))
}

# Generate a re-ordering of coeff names (group choices together)
reordering <- function(varList, choices)
{
  if (length(varList) == 0) return(NULL)
  K <- length(as.vector(choices))
  p <- length(as.vector(varList))
  # p <- length(as.vector(varList$indSpVar))
  # f <- length(as.vector(varList$csvChCoeff))
  # d <- length(as.vector(varList$csvGenCoeff))
  orig <-  c(if (p > 0) rep(1:p, K-1) else NULL)
  # orig <-  c(if (p > 0) rep(1:p, K-1) else NULL,
  #            if (f > 0) rep((p+1):(p+f), K) else NULL,
  #            if (d > 0) (p+f+1):(p+f+d) else NULL)
  order(orig)
}

getNullSpaceCols <- function(mat, tol = 1e-07) {
  if (is.null(mat))
    return(NULL)
  if (ncol(mat) == 1)
    return(NULL)
  qrdecomp <- qr(mat, tol = tol)
  rank <- qrdecomp$rank
  if (rank == ncol(mat))
    return(NULL)
  nullSpCols <- qrdecomp$pivot[(rank + 1):ncol(mat)]
  return(nullSpCols)
}

formDesignMat_r <- function(data, N, formula, varVec = NULL, includeIntercept = TRUE) {
  # if (is.null(varVec) && !includeIntercept)
  #   return(NULL)
  fm <- paste(attr(formula, "response"), "~")
  if (!is.null(varVec))
    fm <- paste(fm, paste(varVec, collapse = "+"))
  if (!includeIntercept)
    fm <- paste(fm, "-1 ")
  else fm <- paste(fm, "+1 ")
  return(model.matrix(as.formula(fm), data[1:N, ]))
}
