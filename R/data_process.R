
to_char_r <- function(x) {
  as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G", "4" = "N", "-1" = "N")[as.character(x)])
}

to_xy_r <- function(x) {
  as.numeric(c("A" = "0", "T" = "2", "C" = "1", "G" = "3")[t(x)])
}

compare_par <- function(new, old, name, tol) {
  mapply(FUN = function(A, B) {
    abs(exp(A) - exp(B)) < tol
  }, A = new[[name]], B = old[[name]]) %>% flatten_lgl()
}

remove_deletion <- function(de, hmm, overlap, hap) {
  del <- list()
  res <- list()
  deletion <- de$delete_t + 1
  del$full_hap <- hap$full_hap[deletion]
  del$time_pos <- hmm$time_pos[deletion]
  del$p_tmax <- hmm$p_tmax[deletion]
  res$del <- del

  len <- length(deletion)
  res$t_max <- hmm$t_max - len
  res$n_t <- hmm$n_t[-deletion]
  res$n_in_t <- hmm$n_in_t[-deletion]
  res$time_pos <- hmm$time_pos[-deletion]
  res$p_tmax <- hmm$p_tmax[-deletion]
  res$num_states <- hmm$num_states[-deletion]
  res$location <- overlap$location[-deletion]
  res$overlapped <- overlap$overlapped[-deletion]
  res$overlapped_id <- de$overlapped_id[-deletion]
  res$full_hap <- hap$full_hap[-deletion]
  res$new_num_states <- hap$new_num_states[-deletion]
  res$combination <- hap$combination[-deletion]
  res$start_t <- de$new_start

  return(res)
}

merge_res <- function(merge, hmm, overlap, hap) {
  res <- list()
  deletion <- merge$targets + 1

  len <- length(deletion)
  res$t_max <- hmm$t_max - len
  res$n_t <- merge$n_t[-deletion]
  res$n_in_t <- merge$n_in_t[-deletion]
  res$time_pos <- merge$new_time_pos[-deletion]
  res$p_tmax <- merge$p_tmax[-deletion]
  res$num_states <- hmm$num_states[-deletion]
  res$location <- overlap$location[-deletion]
  res$overlapped <- overlap$overlapped[-deletion]
  res$overlapped_id <- merge$overlapped_id[-deletion]
  res$full_hap <- merge$full_hap[-deletion]
  res$new_num_states <- hap$new_num_states[-deletion]
  res$combination <- hap$combination[-deletion]
  res$start_t <- merge$start_t

  return(res)
}

