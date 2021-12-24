
ini_par <- function(dat, n_observation, formula, n_class, weight_id, num_cat, ncores, old_version) {
  par <- list()
  #set.seed(seed)
  #par$eta <- runif(n_class, 0, 1)
  #par$eta <- par$eta/sum(par$eta)
  par$eta <- rep(1/num_cat, n_class)
  par$wic <- matrix(1/num_cat, nrow = n_observation, ncol = n_class)
  if(old_version) {
    par$ins_rate <- 1e-5
    par$del_rate <- 1e-5
    par$excluded_read <- rep(0, n_observation)
  }

  Mpar <- Mstep(dat, read_length = NULL, par, weight_id = weight_id, formula, num_cat, ncores, weights = FALSE)
  par$beta <- Mpar$beta

  return(par)
}
