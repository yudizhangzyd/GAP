#ifndef DATA_FORMAT_H
#define DATA_FORMAT_H

#include <Rcpp.h>
using namespace Rcpp;
DataFrame format_data(List dat_info, IntegerMatrix haplotype, int time_pos = -1);
IntegerMatrix format_data_simple(List dat_info, IntegerMatrix haplotype, int time_pos = -1);
int top_n_map(List unique_map);
int uni_sum(List unique_map, unsigned int cut_off);
List linkage_info(List dat_info, IntegerVector undecided_pos);
#endif
