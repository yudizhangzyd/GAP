#ifndef BAUMWELCH_H
#define BAUMWELCH_H

#include <Rcpp.h>
using namespace Rcpp;

double site_likelihood(unsigned int i, unsigned int K, NumericMatrix beta, unsigned int PD_LENGTH,
                        List dat_info, IntegerMatrix haplotype, unsigned int time_pos);
List find_combination(IntegerVector undecided_pos, IntegerVector pos_possibility,
                      unsigned int p_tmax, unsigned int time_pos, int hap_min_pos);
IntegerMatrix fill_all_hap(List hidden_states, unsigned int hap_length, IntegerVector n_row);
List get_overlap(Nullable<List> del_ = R_NilValue, Nullable<List> coverage_ = R_NilValue,
                 Nullable<List> hmm_ = R_NilValue, Nullable<IntegerVector> pos_possi_ = R_NilValue, int hap_min_pos = 0);
#endif
