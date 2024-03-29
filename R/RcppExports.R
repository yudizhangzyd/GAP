# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

get_overlap <- function(del_, coverage_, hmm_, pos_possi_, hap_min_pos) {
    .Call('_GAP_get_overlap', PACKAGE = 'GAP', del_, coverage_, hmm_, pos_possi_, hap_min_pos)
}

full_hap_new <- function(hmm_info, linkage_info, overlap_info, hap_length, hap_min_pos, use_MC = 0L, db_heter = 0L) {
    .Call('_GAP_full_hap_new', PACKAGE = 'GAP', hmm_info, linkage_info, overlap_info, hap_length, hap_min_pos, use_MC, db_heter)
}

reads_llk <- function(hmm_info, dat_info, hap_info, beta, eta, PD_LENGTH, chosed_state) {
    .Call('_GAP_reads_llk', PACKAGE = 'GAP', hmm_info, dat_info, hap_info, beta, eta, PD_LENGTH, chosed_state)
}

format_data2 <- function(hmm_info, d_info, hap_info) {
    .Call('_GAP_format_data2', PACKAGE = 'GAP', hmm_info, d_info, hap_info)
}

baum_welch_init <- function(hmm_info, data_info, hap_info, PD_LENGTH, par, trans_indicator, hash_idx, trans_constraint, db, penality) {
    .Call('_GAP_baum_welch_init', PACKAGE = 'GAP', hmm_info, data_info, hap_info, PD_LENGTH, par, trans_indicator, hash_idx, trans_constraint, db, penality)
}

baum_welch_iter <- function(hmm_info, par_hmm, data_info, hap_info, beta, PD_LENGTH, hash_idx, no_emi_upt) {
    .Call('_GAP_baum_welch_iter', PACKAGE = 'GAP', hmm_info, par_hmm, data_info, hap_info, beta, PD_LENGTH, hash_idx, no_emi_upt)
}

trans_permit <- function(num_states, overlap_info, combination, t_max) {
    .Call('_GAP_trans_permit', PACKAGE = 'GAP', num_states, overlap_info, combination, t_max)
}

trans_const <- function(overlap_info, combination, db_sites, num_states, t_max) {
    .Call('_GAP_trans_const', PACKAGE = 'GAP', overlap_info, combination, db_sites, num_states, t_max)
}

find_deleted <- function(hmm_info, overlap_info) {
    .Call('_GAP_find_deleted', PACKAGE = 'GAP', hmm_info, overlap_info)
}

merge_no_connection <- function(hmm_info, overlap_info, full_hap) {
    .Call('_GAP_merge_no_connection', PACKAGE = 'GAP', hmm_info, overlap_info, full_hap)
}

merge_states <- function(hmm_info, overlap_info, full_hap) {
    .Call('_GAP_merge_states', PACKAGE = 'GAP', hmm_info, overlap_info, full_hap)
}

read_data <- function(path, old_v) {
    .Call('_GAP_read_data', PACKAGE = 'GAP', path, old_v)
}

format_data <- function(dat_info, haplotype, time_pos = -1L) {
    .Call('_GAP_format_data', PACKAGE = 'GAP', dat_info, haplotype, time_pos)
}

to_char <- function(nuc) {
    .Call('_GAP_to_char', PACKAGE = 'GAP', nuc)
}

hmm_info <- function(dat_info, uni_alignment, opt, sbs = 1L) {
    .Call('_GAP_hmm_info', PACKAGE = 'GAP', dat_info, uni_alignment, opt, sbs)
}

get_hidden_state <- function(hmm_info, dat_info, uni_alignment, opt) {
    .Call('_GAP_get_hidden_state', PACKAGE = 'GAP', hmm_info, dat_info, uni_alignment, opt)
}

linkage_info <- function(dat_info, undecided_pos, uni_alignment, nuc_count, nuc_unique) {
    .Call('_GAP_linkage_info', PACKAGE = 'GAP', dat_info, undecided_pos, uni_alignment, nuc_count, nuc_unique)
}

sample_hap <- function(dat_info, start, idx, hap_deletion_len) {
    .Call('_GAP_sample_hap', PACKAGE = 'GAP', dat_info, start, idx, hap_deletion_len)
}

sample_hap2 <- function(hmm_info, hap_length, hap_min_pos) {
    .Call('_GAP_sample_hap2', PACKAGE = 'GAP', hmm_info, hap_length, hap_min_pos)
}

sort_ind <- function(input, n) {
    .Call('_GAP_sort_ind', PACKAGE = 'GAP', input, n)
}

formDesignMat_c <- function(dat, N, PD_LENGTH) {
    .Call('_GAP_formDesignMat_c', PACKAGE = 'GAP', dat, N, PD_LENGTH)
}

determine_hidden <- function(overlap_info, link_in, opt, hmm_info, uni_alignment, hap_min_pos) {
    .Call('_GAP_determine_hidden', PACKAGE = 'GAP', overlap_info, link_in, opt, hmm_info, uni_alignment, hap_min_pos)
}

full_hap_samp2 <- function(hmm_info, overlap_info, linkage, hap_length, hap_min_pos, lower_ab) {
    .Call('_GAP_full_hap_samp2', PACKAGE = 'GAP', hmm_info, overlap_info, linkage, hap_length, hap_min_pos, lower_ab)
}

viterbi <- function(hmm_info, dat_info, hap_info, overlap_info, par_hmm, left_) {
    .Call('_GAP_viterbi', PACKAGE = 'GAP', hmm_info, dat_info, hap_info, overlap_info, par_hmm, left_)
}

connect_hap <- function(hmm_info, dat_info, hap_info, overlap_info) {
    .Call('_GAP_connect_hap', PACKAGE = 'GAP', hmm_info, dat_info, hap_info, overlap_info)
}

