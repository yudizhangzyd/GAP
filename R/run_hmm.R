opts <- new.env()
assign("majority", 0.9, envir = opts)
assign("cut_off", 0, envir = opts)
assign("three_hap", 0.62, envir = opts)
assign("third_nuc", 0, envir = opts)
assign("two_hap", 0.45, envir = opts)
assign("left_range", 1.2, envir = opts)
assign("emprical", 1, envir = opts)
assign("sampling", 0, envir = opts)
assign("n_sample", 10, envir = opts)
assign("lower_ab", 0.95, envir = opts)
assign("single", 0, envir = opts)
assign("no_emi_update", 1, envir = opts)
assign("dbheter", 0, envir = opts)
assign("db_penal", 0.5, envir = opts)

setOpt <- function(...) {
  opts <- getOpt()
  args <- list(...)
  if(length(args)==1 && is.list(args[[1]])) {
    args <- args[[1]]
  }
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) { # class() OK here, since all opts are simple objects with single classes
      if(class(getOpt(opnm)) == class(args[[opnm]])) {
        assign(opnm, args[[opnm]], envir=opts)
      } else {
        warning(paste0(opnm, " not set, value provided has different class (", class(args[[opnm]]),
                       ") then current option value (", class(getOpt(opnm)), ")"))
      }
    } else {
      warning(opnm, " is not a valid option.")
    }
  }
}

getOpt <- function(option = NULL) {
  if(is.null(option)) option <- ls(opts)

  if(!all(option %in% ls(opts))) {
    warning("Tried to get an invalid option: ", option[!(option %in% ls(opts))])
    option <- option[option %in% ls(opts)]
  }

  ropts <- lapply(option, function(x) get(x, envir=opts))
  names(ropts) <- option
  if(length(ropts) == 1) ropts <- ropts[[1]]  # If just one option requested, return it alone
  return(ropts)
}

#' HMM for genotyping and phasing for alleotetraploid individual
#' @param alignment Universal reference genome (Generated from calling `call_aln`)
#' @param datafile Input data file path (Generated from calling `call_aln`)
#' @param res_file Output the results in a file. Default: NULL
#' @param seed Set seed. Default:0
#' @param tol Convergence tolerance. Default:1e-05
#' @param max_iter Iteration in HMM. Default:20
#' @param ncores Cores to use when calling `mnlogit`
#' @param verbose Print process. Default: FALSE
#' @param formula Formula pass to `mnlogit`
#' @param ... Set options
#'
#' @import tibble
#' @import stringr
#' @import dplyr
#' @import purrr
#' @import foreach
#' @import mnlogit
#' @import doParallel
#' @import Formula
#' @import RcppArmadillo
#' @import Rcpp
#' @importFrom magrittr %>%
#' @importFrom stats as.formula formula median model.matrix terms
#' @importFrom data.table rbindlist
#'
#' @export altraphase
#' @return A list contain the phased haplotypes, SNPs and SNP locations and other auxiliary values.
#'
#' @examples
#' datFile <- system.file("extdata", "out.txt", package = "GAP")
#' refFile <- system.file("extdata", "uni.fa", package = "GAP")
#' res <- altraphase(datafile = datFile, alignment = refFile, max_iter = 10)
altraphase <- function(datafile = NULL, alignment = NULL, res_file = NULL,
                       formula = mode~1|read_pos + qua + hap_nuc + qua:hap_nuc, max_iter = 20,
                       seed = 0, tol = 1e-05, ncores = 2, verbose = FALSE, ...)  {
  registerDoParallel(cores = ncores)
  call <- sys.call(1)
  # Read in default opts and then replace with any that were passed in to the function
  opts <- getOpt()
  args <- list(...)
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) {
      opts[[opnm]] <- args[[opnm]]
    } else {
      warning(opnm, " is not a valid option.")
    }
  }
  set.seed(seed)
  res <- list()
  start.time <- Sys.time()

  use_MC = !opts$sampling
  universal <- readLines(alignment)[2] %>% str_split("") %>% unlist()
  dat_info <- read_data(datafile, old_v = 0)
  HMM <- hmm_info(dat_info = dat_info, uni_alignment = universal, opt = opts)
  pre_hs <- get_hidden_state(hmm_info = HMM, dat_info = dat_info,
                               uni_alignment = universal, opt = opts)
  HMM <- append(HMM, pre_hs)
  rm(pre_hs)
  ########################## baum-welch (iterate until converge)

  ## initialization
  ##### use linkage info to limit some unlikely happened transition
  hap_length <- dat_info$ref_length_max - dat_info$ref_start
  original_nuc <- dat_info$nuc ## record nuc
  linkage_in <- linkage_info(dat_info = dat_info, undecided_pos = HMM$undecided_pos,
                             uni_alignment = universal, nuc_count = HMM$nuc_count, nuc_unique = HMM$nuc_unique)
  mutated_nuc <- linkage_in$obs
  dat_info$nuc <- linkage_in$obs

  ini_s <- Sys.time()
  if(opts$emprical == 1) {
    overlap_info <- get_overlap(del_ = linkage_in$deletion_snp, coverage_ = linkage_in$coverage,
                                hmm_ = HMM, pos_possi_ = HMM$pos_possibility, hap_min_pos = dat_info$ref_start)
    HMM$num_states <- overlap_info$num_states
    linkage_in <- linkage_in$link
    hap_full_info <- full_hap_new(HMM, linkage_in, overlap_info, hap_length, dat_info$ref_start,
                                  use_MC = use_MC, db_heter = opts$dbheter)
  } else {
    overlap_info <- get_overlap(del_ = linkage_in$deletion_snp, coverage_ = linkage_in$coverage,
                                hmm_ = HMM, pos_possi_ = NULL, hap_min_pos = dat_info$ref_start)
    linkage_in <- linkage_in$link
    hap_pool <- determine_hidden(overlap_info = overlap_info, link_in = linkage_in, opt = opts,
                                 hmm_info = HMM, uni_alignment = universal, hap_min_pos = dat_info$ref_start)
    HMM[names(hap_pool)] <- hap_pool
    rm(hap_pool)
    hap_full_info <- full_hap_samp2(hmm_info = HMM, overlap_info = overlap_info, hap_length = hap_length,
                                    hap_min_pos = dat_info$ref_start, lower_ab = opts$lower_ab, linkage = linkage_in)
  }
  ### if each hidden state only has one possibility, then get the output directly
  if(all(hap_full_info$new_num_states == 1)) {
    haplotypes <- connect_hap(hmm_info = HMM, dat_info = dat_info, hap_info = hap_full_info$full_hap,
                              overlap_info = overlap_info)
    end.time <- Sys.time()
    res <- list()
    snp_location <- HMM$undecided_pos + 1
    snps <- haplotypes[, snp_location]
    snp_location <- snp_location + dat_info$ref_start
    res$start_pos <- dat_info$ref_start
    res$haplotypes$hap_final <- haplotypes
    res$snp_location <- snp_location
    res$snps <- snps
    res$combination <- hap_full_info$combination
    res$loci <- overlap_info$location
    res$undecided_pos <- HMM$undecided_pos
    res$cov_record <- HMM$cov_record
    res$hap_block <- overlap_info$start_t
    res$time <- end.time - start.time
    res$iter <- 0

    if(!is.null(res_file)) {
      saveRDS(res, res_file)
      return()
    } else
      return(res)
  }
  ### some states and corresponding values should be removed! e.g.: when state 0 has site 0 and state 1
  ### has nothing but state 2 has 01, then transition might be wrong
  de <- find_deleted(hmm_info = HMM, overlap_info = overlap_info)
  left <- NULL
  if(length(de$delete_t) != 0) {
    removal <- remove_deletion(de = de, hmm = HMM, overlap = overlap_info, hap = hap_full_info)
    HMM[intersect(names(HMM), names(removal))] <- removal[intersect(names(HMM), names(removal))]
    hap_full_info[intersect(names(hap_full_info), names(removal))] <- removal[intersect(names(hap_full_info), names(removal))]
    overlap_info[intersect(names(overlap_info), names(removal))] <- removal[intersect(names(overlap_info), names(removal))]
    left <- removal$del
    rm(removal)
  }
  rm(de)
  HMM$num_states <- hap_full_info$new_num_states

  merge_info <- merge_no_connection(hmm_info = HMM, overlap_info = overlap_info, full_hap = hap_full_info$full_hap)
  if(length(merge_info$targets) != 0) {
    removal <- merge_res(merge = merge_info, hmm = HMM, overlap = overlap_info, hap = hap_full_info)
    HMM[intersect(names(HMM), names(removal))] <- removal[intersect(names(HMM), names(removal))]
    hap_full_info[intersect(names(hap_full_info), names(removal))] <- removal[intersect(names(hap_full_info), names(removal))]
    overlap_info[intersect(names(overlap_info), names(removal))] <- removal[intersect(names(overlap_info), names(removal))]
    rm(removal)
    rm(merge_info)
  } else
    HMM$time_pos <- as.list(HMM$time_pos)

  merge_info <- merge_states(hmm_info = HMM, overlap_info = overlap_info, full_hap = hap_full_info$full_hap)
  if(length(merge_info$targets) != 0) {
    removal <- merge_res(merge = merge_info, hmm = HMM, overlap = overlap_info, hap = hap_full_info)
    HMM[intersect(names(HMM), names(removal))] <- removal[intersect(names(HMM), names(removal))]
    hap_full_info[intersect(names(hap_full_info), names(removal))] <- removal[intersect(names(hap_full_info), names(removal))]
    overlap_info[intersect(names(overlap_info), names(removal))] <- removal[intersect(names(overlap_info), names(removal))]
    rm(removal)
    rm(merge_info)
  }
  ###indicate which transfer could happen
  hap_full <- hap_full_info$full_hap
  trans_indicator <- trans_permit(num_states = HMM$num_states, combination = hap_full_info$combination,
                                  overlap_info = overlap_info, t_max = HMM$t_max)
  ### if double heter, then mark the transitions
  trans_constraint = 0
  if(opts$dbheter)
    trans_constraint <- trans_const(overlap_info = overlap_info, combination = hap_full_info$combination,
                                    db_sites = HMM$db_loci, num_states = HMM$num_states, t_max = HMM$t_max)

  ### start initializing
  ## initialize hap
  dat_info$nuc <- original_nuc
  hap_info <- sample_hap2(hmm_info = HMM, hap_length = hap_length, hap_min_pos = dat_info$ref_start)
  hapinit <- hap_info$haplotype
  data <- format_data(dat_info = dat_info, haplotype = hapinit)
  if(hap_info$gap_in) {
    data_rm <- data %>% filter(mode == 1)
    data <- data %>% filter(hap_nuc != -1) # mnlogit only takes data without indels in read or in haplotypes
  }
  data$nuc <- to_char_r(data$nuc)
  data$hap_nuc <- to_char_r(data$hap_nuc)
  id <- data["id"]
  data <- data[, !names(data) %in% c("id")]
  par <- list()
  n_class = 4
  num_cat = 4
  par <- ini_par(dat = data, n_observation = dat_info$n_observation, formula = formula, old_version = 0,
                 n_class = n_class, num_cat = num_cat, ncores = ncores, weight_id = NULL)
  weight_id <- NULL

  cat("initialization: \n");
  data_new <- format_data2(hmm_info = HMM, d_info = dat_info, hap_info = hap_full)
  data <- data_new$df_new
  bw <- baum_welch_init(hmm_info = HMM, data_info = dat_info, hap_info = hap_full, par = par,
                        PD_LENGTH = nrow(par$beta), trans_indicator = trans_indicator,
                        hash_idx = data_new$idx, trans_constraint = trans_constraint,
                        db = opts$dbheter, penality = opts$db_penal)

  if(hap_info$gap_in) {
    data_rm <- data %>% filter(mode == 1)
    weight_id <- which(data_rm$hap_nuc == -1)
    if(length(weight_id) == 0)
      weight_id <- NULL
    data <- data %>% filter(hap_nuc != -1) # mnlogit only takes data without indels in read or in haplotypes
  }
  data$nuc <- to_char_r(data$nuc)
  data$hap_nuc <- to_char_r(data$hap_nuc)
  id <- data["id"]

  ini_e <- Sys.time()
  res$ini_time <- ini_e - ini_s
  start.hmm <- Sys.time()
  for (m in (1:max_iter)) {
    cat("iter: ", m, "\n")
    full_llk <- bw$par_hmm_bf$full_llk
    par_hmm_old <- bw$par_hmm
    phi_old <- bw$par_hmm$phi
    if(m > 1 & opts$no_emi_update) {
      update = 1
    } else {
      cat("no update\n")
      tmp <- m_beta(res = bw$par_aux, id = id, weight_id = weight_id, data = data, formula = formula,
                    ncores = ncores, old_version = 0, weight = bw$par_aux$weight)
      par$beta <- tmp$par$beta
      update = 0
    }
    ## estimation other parameters
    dat_info$nuc <- mutated_nuc
    bw <- baum_welch_iter(hmm_info = HMM, par_hmm = bw, data_info = dat_info, hap_info = hap_full,
                          beta = par$beta, PD_LENGTH = nrow(par$beta), hash_idx = data_new$idx,
                          no_emi_upt = update)
    # cat(bw$par_aux$eta, "\n")
    if (abs(bw$par_hmm_bf$full_llk - full_llk) < tol |
        ((all(abs(exp(bw$par_hmm$phi) - exp(phi_old)) < tol) == TRUE) &&
         all(compare_par(new = bw$par_hmm, old = par_hmm_old, name = "emit", tol) == TRUE) &&
         all(compare_par(new = bw$par_hmm, old = par_hmm_old, name = "trans", tol) == TRUE)))
      break;
  }
  end.hmm <- Sys.time()
  HMM_time <- end.hmm - start.hmm
  res$HMM_time <- HMM_time
  cat("HMM time: ", HMM_time, "\n")
  ### viterbi decoding
  cat("viterbi decoding\n");

  hap <- viterbi(hmm_info = HMM, dat_info = dat_info, hap_info = hap_full,
                 par_hmm = bw$par_hmm, left_ = left, overlap_info = overlap_info)
  haplotypes <- matrix(to_char_r(hap$hap_final), nrow = n_class)
  snp_location <- HMM$undecided_pos + 1
  snps <- haplotypes[, snp_location]
  snp_location <- snp_location + dat_info$ref_start

  end.time <- Sys.time()

  rl <- reads_llk(hmm_info = HMM, dat_info = dat_info, hap_info = hap_full, beta = par$beta,
                  eta = bw$par_aux$eta, PD_LENGTH = 10, chosed_state = hap$chosed_state)
  res$rl <- rl
  res$time <- end.time - start.time
  res$start_pos <- dat_info$ref_start
  res$haplotypes <- hap
  res$snp_location <- snp_location
  res$snps <- snps
  res$gamma <- bw$par_hmm_bf$gamma
  res$combination <- hap_full_info$combination
  res$loci <- overlap_info$location
  res$undecided_pos <- HMM$undecided_pos
  res$cov_record <- HMM$cov_record
  res$hap_block <- overlap_info$start_t
  res$iter <- m

  if(!is.null(res_file))
    saveRDS(res, res_file)
  else
    return(res)
}






