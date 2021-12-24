#ifndef SAMPLING_H
#define SAMPLING_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

void sampling_link(vector<int> &res, int in_id, vector<vector<int> > reads_pool,
                   vector<vector<int> > start_end, int last_pos, int flag);
void cat_reads_forward(vector<int> &partial, vector<int> res, vector<vector<int> > start_end,
                       int in_id, int s, IntegerMatrix uni_read);
void cat_reads_backward(vector<int> &partial, vector<int> res, vector<vector<int> > start_end,
                        int in_id, int e, IntegerMatrix uni_read);
IntegerMatrix sample_hap(std::vector<string> map, IntegerMatrix uni_read,
                         List all_id, int sample);
List limit_comb_samp(IntegerMatrix combination, List hidden_states, IntegerVector location,
                     List hap_block, List block_sites, int t, unsigned int num_states,
                     double lower_ab, IntegerMatrix linkage);
IntegerMatrix determine_hap(IntegerMatrix hap_pool, double lower_ab, IntegerMatrix linkage);
#endif
