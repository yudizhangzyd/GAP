
# @description fit mnlogit and update beta
# @param dat expanded data
# @param par parameters
# @param ncores number of ncores to run mnlogit
# @return betas and CE_llk

Mstep <- function(dat, read_length = NULL, par, weight_id, formula, num_cat = 4, ncores, weights = TRUE, given_weights = NULL) {
  #give a initial value for optimization
  if (weights) {
    if(is.null(given_weights)) { ## for old_tpphase
      weights <- data.table::rbindlist(foreach(i = 1:ncol(par$wic)) %dopar% {
        data.frame(wic = rep.int(par$wic[, i], read_length[i]))
      })
      if(!is.null(weight_id))
        weights <- weights[-c(weight_id), ]
      weights <- weights$wic
    } else {
      weights <- given_weights
      if(!is.null(weight_id))
        weights <- weights[-c(weight_id)]
    }

    start <- par$beta %>% as.vector
    fit <- modified_mnlogit(formula = formula,
                            data = dat, weights = weights, choiceVar = "nuc", ncores = ncores, start = start)
  }
  else {
    fit <- modified_mnlogit(formula = formula,
                            data = dat, choiceVar = "nuc", ncores = ncores)
  }

  #A <- fit$coefficients[which(str_detect(attr(fit$coefficients, "names"), ":A") == 1)]
  C <- fit$coefficients[which(str_detect(attr(fit$coefficients, "names"), ":C") == 1)]
  G <- fit$coefficients[which(str_detect(attr(fit$coefficients, "names"), ":G") == 1)]
  T <- fit$coefficients[which(str_detect(attr(fit$coefficients, "names"), ":T") == 1)]

  mat <- matrix(cbind(C, G, T), ncol = 3)
  beta <- mat
  logLik <- fit$logLik

  res <- list()
  res$logLik <- logLik
  res$beta <- beta

  return(res)
}

# @description prepare data and call mnlogit
# @return logLik of mnlogit and betas

m_beta <- function(res, weight_id, data, id, formula, reads_lengths = NULL, ncores, old_version, weight) {
  par <- list()

  if(length(res$excluded_id) != 0) {
    data_rm <- data %>% filter(mode == 1)
    weight_id <- c(weight_id, which(data_rm$id %in% res$excluded_id))
    data <- data %>% filter(!(id %in% res$excluded_id))
  }

  if(old_version)
    par$wic <- t(res$param$w_ic)
  par$beta <- res$param$beta #beta from last step as starting value
  data <- data[, !names(data) %in% c("id")] # (modified_moligit exclude first column: id)

  Mpar <- Mstep(dat = data, read_length = reads_lengths, weight_id = weight_id, formula = formula, par = par,
                ncores = ncores, weights = TRUE, given_weights = weight)
  par$beta <- Mpar$beta

  if(old_version) {
    par$wic <- res$param$w_ic ## For the use of update haplotype
    par$eta <- res$param$mixture_prop
    par$del_rate <- res$param$del_rate
    par$ins_rate <- res$param$ins_rate
    par$excluded_read <- res$param$excluded_read
  }

  results <- list()
  results$par <- par
  results$CEllk <- Mpar$logLik

  return(results)
}



