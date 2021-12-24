#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "data_format.h"
#include "baumwelch.h"
#include "hmm_state.h"
#include "utils.h"
using namespace Rcpp;
using namespace std;

#define MLOGIT_CLASS 4
#define NUM_CLASS 4

// [[Rcpp::depends(RcppArmadillo)]]
// if double heter exists, then penalize the initialize state (low probability)
List ini_hmm (unsigned int t_max,  IntegerVector num_states,
              List trans_indicator) {
  NumericVector phi(num_states[0]);
  List trans(t_max - 1);
  // List emit(t_max);
  unsigned int w, m, t;

  // change this to log scale
  for(w = 0; w < num_states[0]; ++w)
    phi[w] = log(1/double(num_states[0]));

  // keep some trans to 0
  for(t = 0; t < t_max - 1; ++t) {
    NumericMatrix transition(num_states[t], num_states[t + 1]);
    // Rcout << t << " " << num_states[t] <<":" <<  num_states[t + 1] << "\n";
    // if(trans_indicator.isNotNull()) {
    IntegerMatrix trans_ind = trans_indicator[t];
    for (w = 0; w < num_states[t]; ++w) {
      int new_num = num_states[t + 1];
      for (m = 0; m < num_states[t + 1]; ++m)
        if (trans_ind(w, m)) // 1 means not the same, so cannot b transferred
          new_num--;
        for (m = 0; m < num_states[t + 1]; ++m)
          if (!trans_ind(w, m))
            transition(w, m) = log(1/double(new_num));
          else
            transition(w, m) = R_NegInf;
    }
    trans(t) = transition;
  }

  List par_hmm = List::create(
    Named("phi") = phi,
    Named("trans") = trans);
  return(par_hmm);
}

List init_penal (List par_hmm, List trans_constraint, IntegerVector num_states,
                 unsigned int t_max, double penality) {
  NumericVector phi = par_hmm["phi"];
  List trans = par_hmm["trans"];
  unsigned int w, m, t;

  IntegerVector constraint = trans_constraint[0];
  // change this to log scale
  int count = 0;
  double norm = 0;
  for (w = 0; w < phi.size(); ++w)
    if (constraint[w]) {
      phi[w] += log(penality);
      norm += phi[w];
      count++;
    }

  if(count) {
    norm = log(1 - exp(norm)) - log(phi.size() - count);
    for (w = 0; w < phi.size(); ++w)
      if (!constraint[w])
        phi[w] = norm;
  }

  for(t = 0; t < t_max - 1; ++t) {
    NumericMatrix transition = trans[t];
    IntegerMatrix constraints = trans_constraint[t + 1];
    for (w = 0; w < num_states[t]; ++w) {
      count = 0;
      norm = 0;
      for (m = 0; m < num_states[t + 1]; ++m)
        if (constraints(w, m)) {
          count++;
          transition(w, m) += log(penality);
          norm += transition(w, m);
        }
        if(!count)
          continue;
        norm = log(1 - exp(norm)) - log(num_states[t + 1] - count);
        for (m = 0; m < num_states[t + 1]; ++m)
          if (!constraints(w, m))
            transition(w, m) = norm;
    }
    trans(t) = transition;
  }
  List par_new = List::create(
    Named("phi") = phi,
    Named("trans") = trans);
  return(par_new);
}


//  to avoid underflow, use log sum exp
List forward_backward(List par_hmm, unsigned int t_max, IntegerVector num_states)
{
  NumericVector phi = par_hmm["phi"];
  List trans = par_hmm["trans"];
  List emit = par_hmm["emit"];

  List alpha(t_max);
  List beta_wt(t_max);
  unsigned int w, t, m;
  double max_penality;

  for (t = 0; t < t_max; ++t) {
    NumericVector alp(num_states[t]);
    NumericVector emission = emit(t);
    NumericVector alp_last;
    if (t == 0)
      for (w = 0; w < num_states[t]; ++w)
        alp(w) = phi[w] + emission(w);
    else {
      alp_last = alpha(t - 1);
      NumericMatrix transition = trans(t - 1);
      for (w = 0; w < num_states[t]; ++w) {
        max_penality = R_NegInf;
        for (m = 0; m < num_states[t - 1]; ++m) {
          double value = alp_last(m) + transition(m, w);
          if (value > max_penality)
            max_penality = value;
        }
        for (m = 0; m < num_states[t - 1]; ++m)
          alp(w) += exp(alp_last(m) + transition(m, w) - max_penality); // Log-Sum-Exp Trick
        alp(w) = log(alp(w)) + emission(w) + max_penality;
      }
    }
    alpha(t) = alp;
  }

  for (t = t_max; t --> 0;) {
    NumericVector beta(num_states[t]);
    NumericVector beta_after;
    if(t == t_max - 1) {
      for (w = 0; w < num_states[t]; ++w)
        beta(w) = 0;
    } else {
      NumericVector emission = emit(t + 1);
      NumericMatrix transition = trans(t);
      beta_after = beta_wt(t + 1);
      for (w = 0; w < num_states[t]; ++w) {
        max_penality = R_NegInf;
        for (m = 0; m < num_states[t + 1]; ++m) {
          double value = beta_after(m) + transition(w, m) + emission(m);
          if (value > max_penality)
            max_penality = value;
        }
        for (m = 0; m < num_states[t + 1]; ++m)
          beta(w) += exp(beta_after(m) + transition(w, m) + emission(m) - max_penality);
        beta(w) = log(beta(w)) + max_penality;
      }
    }
    beta_wt(t) = beta;
  }

  List gamma(t_max);
  double sum = 0;
  for (t = 0; t < t_max; ++t) {
    NumericVector beta = beta_wt(t);
    NumericVector alp = alpha(t);
    sum = 0;
    NumericVector gam(num_states[t]);
    max_penality = R_NegInf;
    for (w = 0; w < num_states[t]; ++w) {
      gam(w) = alp(w) + beta(w);
      if (gam(w) > max_penality)
        max_penality = gam(w);
    }
    for (w = 0; w < num_states[t]; ++w)
      sum += exp(gam(w) - max_penality); // log sum exp trick
    for (w = 0; w < num_states[t]; ++w)
      gam(w) = gam(w) - (log(sum) + max_penality);
    gamma(t) = gam;
  }

  List xi(t_max - 1);
  //arma::Cube<double> xi(total_state, total_state, t_max);
  for (t = 0; t < t_max - 1; ++t) {
    NumericVector emission = emit(t + 1);
    NumericMatrix transition = trans(t);
    NumericVector beta = beta_wt(t + 1);
    NumericVector alp = alpha(t);
    NumericMatrix x(num_states[t], num_states[t + 1]);
    sum = 0;
    max_penality = R_NegInf;
    for (w = 0; w < num_states[t]; ++w)
      for (m = 0; m < num_states[t + 1]; ++m) {
        x(w, m) = alp(w) + transition(w, m) + beta(m) + emission(m);
        if (x(w, m) > max_penality)
          max_penality = x(w, m);
      }
      for (w = 0; w < num_states[t]; ++w)
        for (m = 0; m < num_states[t + 1]; ++m)
          sum += exp(x(w, m) - max_penality);
    for (w = 0; w < num_states[t]; ++w)
      for (m = 0; m < num_states[t + 1]; ++m)
        x(w, m) = x(w, m) - (log(sum) + max_penality);
    xi(t) = x;
  }
  // compute full likelihood
  NumericVector alp = alpha[t_max - 1];
  double full_llk = 0;
  max_penality = max(alp);
  for (w = 0; w < num_states[t_max - 1]; ++w)
    full_llk += exp(alp[w] - max_penality);
  full_llk = log(full_llk) + max_penality;
  Rcout << "full log likelihood: " << full_llk << "\n";

  List par_hmm_out = List::create(
    Named("phi") = phi,
    Named("trans") = trans,
    Named("emit") = emit,
    Named("alpha") = alpha,
    Named("beta_wt") = beta_wt,
    Named("gamma") = gamma,
    Named("xi") = xi,
    Named("full_llk") = full_llk);
  return(par_hmm_out);
}

NumericVector update_phi(List par_hmm, unsigned int num_states) {
  List gamma = par_hmm["gamma"];
  NumericVector phi(num_states);
  NumericVector gam = gamma[0];

  for (unsigned int w = 0; w < num_states; ++w)
    phi[w] = gam(w);
  return phi;
}

List update_trans(List par_hmm, unsigned int t_max, IntegerVector num_states) {
  List gamma = par_hmm["gamma"];
  List xi = par_hmm["xi"];
  List trans(t_max - 1);
  unsigned int t, m, w;
  for (t = 0; t < t_max - 1; ++t) {
    NumericMatrix transition(num_states[t], num_states[t + 1]);
    NumericVector gam = gamma[t];
    NumericMatrix x = xi[t];
    for (w = 0; w < num_states[t]; ++w)
      for (m = 0; m < num_states[t + 1]; ++m)
        transition(w, m) = x(w, m) - gam(w);
    trans[t] = transition;
  }
  return(trans);
}

double site_likelihood (unsigned int i, unsigned int K, NumericMatrix beta, unsigned int PD_LENGTH,
                        List dat_info, IntegerMatrix haplotype, unsigned int time_pos) {
  unsigned int j, l;
  double qua_in, read_pos_in, ref_pos_in;
  double tail, sum;
  double xb = 0;

  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];
  IntegerVector length = dat_info["length"];
  IntegerVector index = dat_info["start_id"];
  IntegerMatrix ref_index = dat_info["ref_idx"];

  NumericVector hap_nuc(MLOGIT_CLASS);
  NumericVector hnuc_qua(MLOGIT_CLASS);
  NumericVector pred_beta(MLOGIT_CLASS - 1);
  double read_llk = 0.;
  arma::vec predictor(PD_LENGTH);
  //NumericVector site_llk(length[i]);

  /* non-indel positions for both reads and haplotypes */
  for(j = 0; j < length[i]; ++j) {
    if(obs[index[i] + j] == -1)
      continue;

    // Rprintf("j %d, position %d\n", j, index[i] + j);
    read_pos_in = read_pos[index[i] + j];
    qua_in = qua[index[i] + j];
    ref_pos_in = ref_pos[index[i] + j]; // this is different since the sub-hap is counted 0 from the time t

    // don't take indels of genomes, apply for any K
    if(haplotype(0, ref_pos_in - time_pos) == -1 || haplotype(NUM_CLASS - 1, ref_pos_in - time_pos) == -1)
      continue;

    for (l = 0; l < MLOGIT_CLASS; ++l) {
      hap_nuc[l] = 0;
      hnuc_qua[l] = 0;
    }
    // Rcout << "haplotype " << haplotype(K, ref_pos_in - time_pos) << "\t";

    if(haplotype(K, ref_pos_in - time_pos) == 1) {
      hap_nuc[0] = 1;
      hnuc_qua[0] = qua_in;
    } else if(haplotype(K, ref_pos_in - time_pos) == 3) {
      hap_nuc[1] = 1;
      hnuc_qua[1] = qua_in;
    } else if(haplotype(K, ref_pos_in - time_pos) == 2) {
      hap_nuc[3] = 1;
      hnuc_qua[3] = qua_in;
    }

    predictor = {1, read_pos_in, qua_in, hap_nuc[0], hap_nuc[1],
                 hap_nuc[3], hnuc_qua[0], hnuc_qua[1], hnuc_qua[3]};

    arma::mat beta_ar = as<arma::mat>(beta);
    arma::vec pb = beta_ar.t() * predictor;
    pred_beta = as<NumericVector>(wrap(pb));

    // Rcout << "predictor : " << predictor.t() << "\n";
    //Rcout << "beta_ar : " << beta_ar << "\n";
    //Rcout << "pred_beta : " << pred_beta << "\n";

    sum = 0.0;
    for (l = 0; l < MLOGIT_CLASS - 1; ++l)
      sum += exp(pred_beta[l]);
    tail = log(1/(1 + sum));

    if(obs[index[i] + j] == 1) {
      xb = pred_beta[0];
    } else if(obs[index[i] + j] == 3) {
      xb = pred_beta[1];
    } else if(obs[index[i] + j] == 2) {
      xb = pred_beta[2];
    } else if(obs[index[i] + j] == 0) {
      xb = 0;
    }
    read_llk += xb + tail; /* Notice here use log likelihood, not likelihood */
  }
  //Rcout << "read_llk: " << read_llk << "\t";
  return read_llk;
} /* site_likelihood */

// NOTE: Rcpp doesn't accept comparison between int and unsigned int!!!
// notice: if the covered location is deletion, and it is the only linked read then a new start
// [[Rcpp::export]]

List get_overlap(Nullable<List> del_, Nullable<List> coverage_,  Nullable<List> hmm_,
                 Nullable<IntegerVector> pos_possi_, int hap_min_pos)
{
  List hmm_info(hmm_);
  IntegerVector num_states = hmm_info["num_states"]; // this need to be updated if any new start found
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  List n_in_t = hmm_info["n_in_t"];
  unsigned int t_max = hmm_info["t_max"];

  IntegerVector start_t(t_max);
  // get the first t which has variation
  for (unsigned int t = 0; t < t_max; ++t)
    if(num_states[t] != 1) {
      start_t[0] = t;
      break;
    }
  List overlapped(t_max);
  List location(t_max);
  IntegerVector overlapped_idx(t_max);
  int begin, end, end1, min;
  unsigned int s_t = 1;

  unordered_map<int, vector<int>> mp;
  unordered_map<int, vector<int>> cover;
  unordered_map<int, int> update_state;
  if(del_.isNotNull()) {
    List deletion_snp(del_);
    List coverage(coverage_);
    for(int i = 0; i < deletion_snp.size(); ++i) {
      vector<int> tmp = deletion_snp[i];
      int site = tmp.back();
      tmp.pop_back();
      mp[site] = tmp;
    }
    for (int j = 0; j < undecided_pos.size(); ++j) {
      vector<int> tmp = coverage[j];
      cover[undecided_pos[j]] = tmp;
    }
  }
  if(pos_possi_.isNotNull()) {
    IntegerVector pos_possibility(pos_possi_);
    for (int j = 0; j < undecided_pos.size(); ++j)
      update_state[undecided_pos[j]] = pos_possibility[j];
  }

  for (unsigned int t = 0; t < t_max; ++t) {
    if (num_states[t] == 1) {
      overlapped[t] = -1;
      overlapped_idx[t] = -1;
      location[t] = -1;
      continue;
    }
    int ss = 0;

    if(t == start_t[0]) {
      overlapped[t] = -1;
      overlapped_idx[t] = -1;
      begin = time_pos[start_t[0]] - hap_min_pos;
      end = time_pos[start_t[0]] + p_tmax[start_t[0]] - hap_min_pos;
      int num = 0;
      IntegerVector location_t(undecided_pos.size());
      for (unsigned int m = 0; m < undecided_pos.size(); ++m)
        if (undecided_pos[m] >= begin && undecided_pos[m] < end)
          location_t(num++) = undecided_pos[m];
      // check if the last site is a deletion
      if(num > 1) {
        vector<int> be;
        IntegerVector reads = n_in_t[t];
        for(int i = 0; i < num - 1; ++i) {
          int first = location_t[i];
          int second = location_t[i + 1];
          if(mp.find(first) != mp.end() || mp.find(second) != mp.end()) {
            IntegerVector index = wrap(mp[first]); // read_id
            IntegerVector index2 = wrap(mp[second]); // read_id
            IntegerVector allreads1 = wrap(cover[first]);
            IntegerVector allreads2 = wrap(cover[second]);
            IntegerVector mutual = intersect(allreads1, allreads2);
            // if no mutual reads, then break the linkage
            // Rcout << first << ",mutual: " << mutual << "|" << index << "\n";
            if(mutual.size() == 0) {
              IntegerVector if_in = intersect(reads, index);
              IntegerVector if_in2 = intersect(reads, index2);
              if(if_in.size() != 0)
                be.push_back(i);
              else if(if_in2.size() != 0)
                be.push_back(i + 1);
            }
          }
        }
        if(!be.empty()) {
          Rcout << t << " linkage got broken\n";
          if(be.size() == 1) {
            if(be[0] == 0)
              ss = be[0] + 1;
            else
              num = be[0];
          } else {
            for(int i = be.size() - 1; i >= 1; --i)
              if(be[i] - be[i - 1] > 1) {
                num = be[i] - 1;
                ss = be[i - 1] + 1;
              }
          }
          if(pos_possi_.isNotNull()) {
            num_states[t] = 1;
            // update the num_states
            for(int i = ss; i < num; ++i)
              num_states[t] *= update_state[location_t[i]];
          }
        }
      }
      // Rcout << ss << " " << num - 1 << "\n";
      location[start_t[0]] = location_t[Range(ss, num - 1)];
      continue;
    }
    begin = time_pos[t] - hap_min_pos;
    end = time_pos[t] + p_tmax[t] - hap_min_pos;
    // store location
    int num = 0;
    IntegerVector location_t(undecided_pos.size());
    for (unsigned int m = 0; m < undecided_pos.size(); ++m)
      if (undecided_pos[m] >= begin && undecided_pos[m] < end)
        location_t(num++) = undecided_pos[m];

    if(num > 1) {
      vector<int> be;
      IntegerVector reads = n_in_t[t];
      // Rcout << reads << "\n";
      for(int i = 0; i < num - 1; ++i) {
        int first = location_t[i];
        int second = location_t[i + 1];
        if(mp.find(first) != mp.end() || mp.find(second) != mp.end()) {
          IntegerVector index = wrap(mp[first]); // read_id
          IntegerVector index2 = wrap(mp[second]); // read_id
          IntegerVector allreads1 = wrap(cover[first]);
          IntegerVector allreads2 = wrap(cover[second]);
          IntegerVector mutual = intersect(allreads1, allreads2);
          // if no mutual reads, then break the linkage
          // Rcout << first << " " << second << ",mutual: " << mutual << "|" << index2 << "|" << index << "\n";
          if(mutual.size() == 0) {
            IntegerVector if_in = intersect(reads, index);
            IntegerVector if_in2 = intersect(reads, index2);
            if(if_in.size() != 0)
              be.push_back(i);
            else if(if_in2.size() != 0)
              be.push_back(i + 1);
          }
        }
      }
      if(!be.empty()) {
        // Rcout << t << " linkage got broken\n";
        // print_c_intvec(be);
        // Rcout << ss << " " << num - 1 << "\n";
        if(be.size() == 1) {
          if(be[0] == 0){
            ss = 1;
          } else // break in the middle or at last
            num = be[0];
        } else {
          for(int i = be.size() - 1; i >= 1; --i)
            if(be[i] - be[i - 1] > 1) {
              num = be[i] - 1;
              ss = be[i - 1] + 1;
            }
        }
        if(pos_possi_.isNotNull()) {
          num_states[t] = 1;
          for(int i = ss; i < num; ++i)
            num_states[t] *= update_state[location_t[i]];
        }
          // Rcout << ss << " " << num - 1 << "\n";
      }
    }
      //
    location[t] = location_t[Range(ss, num - 1)];

    int len = -1;
    int id_t = 0;
    int index = 0;

    int flag = 0; //indicate if this has overlapped
    for (unsigned int t1 = 0; t1 < t; ++t1) {
      if(num_states[t1] == 1)
        continue;
      IntegerVector last_location = location[t1];
      // begin1 = time_pos[t1] - hap_min_pos;
      end1 = last_location[last_location.size() - 1];
      // Rcout << "end1 " << end1 << "location_t[0] " << location_t[0] << "\n";
      if(location_t[ss] <= end1) {
        flag = 1;
        min = end1;
        // minimum overlapped region
        if(location_t[num - 1] < end1)
          min = location_t[num - 1];
        // find the time t which has the longest coverage
        if(min - location_t[ss] > len) {
          len = min - location_t[ss];
          id_t = min;
          index = t1;
        }
      }
    }
    if(!flag) {
      start_t[s_t++] = t;
      overlapped[t] = -1;
      overlapped_idx[t] = -1;
    } else {
      int count = 0;
      IntegerVector position(undecided_pos.size());
      for (unsigned int m = 0; m < undecided_pos.size(); ++m)
        if (undecided_pos[m] >= begin && undecided_pos[m] <= id_t)
          position(count++) = undecided_pos[m];
        overlapped[t] = position[Range(0, count - 1)];
        overlapped_idx[t] = index;
    }
  }

  List overlap = List::create(
    Named("location") = location,
    Named("overlapped") = overlapped,
    Named("overlapped_id") = overlapped_idx,
    Named("start_t") = start_t[Range(0, s_t - 1)],
                              Named("num_states") = num_states);
  return(overlap);
}


// [[Rcpp::export]]
List full_hap_new (List hmm_info, IntegerMatrix linkage_info, List overlap_info, unsigned int hap_length,
                   int hap_min_pos, unsigned int use_MC = 0, unsigned int db_heter = 0) {
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector n_row = hmm_info["n_row"];
  IntegerVector pos_possibility = hmm_info["pos_possibility"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int t, m;
  List full_hap(t_max);
  List comb(t_max);
  IntegerMatrix hap = fill_all_hap(hidden_states, hap_length, n_row);
  IntegerVector new_num_states(t_max);

  List overlapped = overlap_info["overlapped"];
  IntegerVector overlapped_id = overlap_info["overlapped_id"];
  IntegerVector start_t = overlap_info["start_t"];
  List loci = overlap_info["location"];
  int n_start = start_t.size();
  //start t info(store this for t that needs to use this info)
  // Rcout << "start t " << start_t << "\n";

  List full_hap_t;

  for(t = 0; t < n_start; ++t) {
    // Rcout << start_t[t] << " ";
    int count = 0;
    IntegerVector location = loci[start_t[t]];
    List comb_info_t0 = find_combination(location, undecided_pos, pos_possibility);

    IntegerMatrix combination = comb_info_t0["combination"];
    unsigned int num = comb_info_t0["num"];
    // find index in linkage matrix
    unsigned int start_idx = 0;
    for(m = 0; m < undecided_pos.size(); ++m)
      if(location[0] == undecided_pos[m]) {
        start_idx = m;
        break;
      }

    List t0 = limit_comb_t0(combination, hidden_states, location, linkage_info, num,
                            start_idx, num_states[start_t[t]], use_MC);
    IntegerVector exclude = t0["exclude"];
    // Rcout << "\n"<< exclude << "\n";
    new_num_states[start_t[t]] = t0["num_states"];
    IntegerMatrix new_comb(new_num_states[start_t[t]], combination.ncol());

    for(m = 0; m < num_states[start_t[t]]; ++m)
      if(!exclude[m])
        new_comb(count++, _) = combination(m, _);
      // dereplicate
      IntegerMatrix final_comb = dereplicate_states(new_comb, hidden_states, location, combination.ncol(),
                                                    new_num_states[start_t[t]], db_heter);
      comb[start_t[t]] = final_comb;
      new_num_states[start_t[t]] = final_comb.nrow();
  }

  IntegerVector exclude_last;
  IntegerMatrix comb_in;
  // get the states
  for(t = 0; t < t_max; ++t) {
    if(num_states[t] != 1 && overlapped_id[t] != -1) {
      // int count = 0;
      // Rcout << t << "\t";
      // int identical = 0;
      int last_t = overlapped_id[t];
      IntegerVector overlapped_t = overlapped[t];
      IntegerVector loci_lastt = loci[last_t];
      IntegerVector loci_currt = loci[t];

      if(loci_lastt[0] <= loci_currt[0] && loci_lastt[loci_lastt.size() - 1] >= loci_currt[loci_currt.size() - 1]) {
        if(loci_lastt.size() > loci_currt.size()) { // if current is in its overlap
          // get the unique overlapped combination from the last t
          IntegerMatrix comb_in = comb[last_t];
          // Rcout << "contained in\n";
          IntegerMatrix new_comb = unique_overlap(overlapped_t, comb_in, loci_lastt, new_num_states[last_t]);
          new_num_states[t] = new_comb.nrow();
          comb[t] = new_comb;
        } else if (loci_lastt.size() == loci_currt.size()) {// if current is same as its overlap
          IntegerMatrix new_comb = comb[last_t];
          comb[t] = new_comb;
          new_num_states[t] = new_num_states[last_t];
          // }
        }
      }
      else {
        // new variable site in this t, need to get the new combination while making sure it can be transferred to the next t
        IntegerMatrix comb_in = comb[last_t];
        IntegerMatrix new_comb = new_combination(hmm_info, loci_currt, overlapped_t, comb_in,
                                                 loci_lastt, linkage_info, new_num_states[last_t], use_MC);
        new_num_states[t] = new_comb.nrow();
        comb[t] = new_comb;
      }
    }
    else if(num_states[t] == 1) {
      new_num_states[t] = 1;
      comb[t] = -1;
    }
  }

  for(t = 0; t < t_max; ++t) {
    int count = 0;

    if(num_states[t] == 1) {
      full_hap_t = List(1);
      full_hap_t[0] = hap(_, Range(time_pos[t] - hap_min_pos, time_pos[t] + p_tmax[t] - hap_min_pos - 1));
    } else {
      IntegerMatrix new_comb = comb[t];
      IntegerVector loci_currt = loci[t];
      full_hap_t = List(new_num_states[t]);
      for(m = 0; m < new_num_states[t]; ++m) {
        IntegerMatrix haplotype = make_hap(hidden_states, hap, loci_currt, p_tmax[t], new_comb(m, _), time_pos[t], loci_currt.size(), hap_min_pos);
        full_hap_t(count++) = haplotype;
      }
    }
    full_hap[t] = full_hap_t;
  }

  List out = List::create(
    Named("full_hap") = full_hap,
    Named("new_num_states") = new_num_states,
    Named("combination") = comb);

  return(out);
}
// [[Rcpp::export]]
List reads_llk(List hmm_info, List dat_info, List hap_info, NumericMatrix beta,
               NumericVector eta, int PD_LENGTH, IntegerVector chosed_state)
{
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  List n_in_t = hmm_info["n_in_t"];
  List time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int i, k, t;
  double sum_emit_prob;
  int n_obs = dat_info["n_observation"];
  NumericVector weight_llk(NUM_CLASS);
  // compute emission based on the initial value of eta and beta
  IntegerVector ass(n_obs);
  NumericVector rl(n_obs);
  int count = 0;
  for(t = 0; t < t_max; ++t) {
    // give h_t, each t has many possible combinations
    List full_hap_t = hap_info(t);
    sum_emit_prob = 0;
    IntegerVector idx = n_in_t[t];
    IntegerVector tp_read = time_pos[t];
    int tp = tp_read[0];
    // Rcout << "t: " << t << "\n";
    // NumericVector read_likelihood(n_t[t]);
    NumericMatrix read_class_llk(n_t[t], NUM_CLASS);
    IntegerMatrix haplotype = full_hap_t(chosed_state[t]);
    for (i = 0; i < n_t[t]; ++i) {
      unsigned int id = idx[i];
      int max_id = 0;
      double max = -INFINITY;
      // Rcout << "id" << id << ": " << tp << "\n";
      for (k = 0; k < NUM_CLASS; ++k) {
        read_class_llk(i, k) = site_likelihood(id, k, beta, PD_LENGTH, dat_info, haplotype, tp);
        if(read_class_llk(i, k) > max) {
          max = read_class_llk(i, k);
          max_id = k;
        } else if(read_class_llk(i, k) == max) {
          double r = ((double) rand() / (RAND_MAX));
          if(r < 0.5)
            max_id = k;
        }
        // weight_llk[k] = eta[k] * exp(read_class_llk(i, k));
        // read_likelihood[i] += weight_llk[k];
      }
      rl[count] =  max;
      ass[count++] = max_id;
    }
  }

  List out = List::create(
    Named("rl") = rl,
    Named("ass") = ass);

  return(out);
}

List compute_emit(List hmm_info, List dat_info, List hap_info, NumericMatrix beta, NumericVector eta, int PD_LENGTH)
{
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  List n_in_t = hmm_info["n_in_t"];
  List time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int i, k, m, t;
  double sum_emit_prob;
  List w_ic(t_max);
  List emit(t_max);
  NumericVector weight_llk(NUM_CLASS);
  // compute emission based on the initial value of eta and beta
  for(t = 0; t < t_max; ++t) {
    // give h_t, each t has many possible combinations
    List full_hap_t = hap_info(t);
    List w_icm(num_states[t]);
    NumericVector emission(num_states[t]);
    sum_emit_prob = 0;
    IntegerVector idx = n_in_t[t];
    IntegerVector tp_read = time_pos[t];
    int tp = tp_read[0];
    // Rcout << "t: " << t << "\n";
    for(m = 0; m < num_states[t]; ++m) {
      //Rcout << m;
      NumericVector read_likelihood(n_t[t]);
      NumericMatrix read_class_llk(n_t[t], NUM_CLASS);
      NumericMatrix w_tic(n_t[t], NUM_CLASS);
      IntegerMatrix haplotype = full_hap_t(m);
      for (i = 0; i < n_t[t]; ++i) {
        unsigned int id = idx[i];
        // Rcout << "id" << id << ": " << tp << "\n";
        for (k = 0; k < NUM_CLASS; ++k) {
          read_class_llk(i, k) = site_likelihood(id, k, beta, PD_LENGTH, dat_info, haplotype, tp);
          weight_llk[k] = eta[k] * exp(read_class_llk(i, k));
          read_likelihood[i] += weight_llk[k];
        }
        for (k = 0; k < NUM_CLASS; ++k)
          w_tic(i, k) = weight_llk[k]/read_likelihood[i];
        emission[m] += log(read_likelihood[i]);
      }
      emission[m] = exp(emission[m]);
      sum_emit_prob += emission[m];
      w_icm(m) = w_tic;
      // Rcout << w_tic;
    }
    w_ic(t) = w_icm;
    // scale the emission prob (log scale)
    for(m = 0; m < num_states[t]; ++m)
      emission[m] = log(emission[m]) - log(sum_emit_prob);
    if(num_states[t] == 1)
      emission[0] = 0;
    emit(t) = emission;
  }
  List par_hmm_out = List::create(
    Named("emit") = emit,
    Named("w_ic") = w_ic);
  return(par_hmm_out);
}

NumericVector update_eta(List w_ic, List gamma, IntegerVector num_states, IntegerVector n_t, unsigned int t_max, unsigned int n_observation) {
  double inner;
  unsigned int i, k, m, t;
  NumericVector eta_new(NUM_CLASS);

  for(k = 0; k < NUM_CLASS; ++k) {
    for(t = 0; t < t_max; ++t) {
      List w_icm = w_ic(t);
      NumericVector gam = gamma(t);
      inner = 0;
      for(m = 0; m < num_states[t]; ++m) {
        NumericMatrix w_tic = w_icm(m);
        double sum_wic = 0;
        for (i = 0; i < n_t[t]; ++i)
          sum_wic += w_tic(i, k);
        inner += exp(gam(m)) * sum_wic;
      }
      eta_new[k] += inner;
    }
  }
  for(k = 0; k < NUM_CLASS; ++k)
    eta_new[k] = eta_new[k]/n_observation;

  return(eta_new);
}

List sub_sample(List hmm_info, List dat_info) {
  unsigned int i, j, t;
  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector obs_index = dat_info["id"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];
  IntegerVector n_t = hmm_info["n_t"];
  IntegerVector index = dat_info["start_id"];
  IntegerVector length = dat_info["length"];
  IntegerVector hap_len = hmm_info["p_tmax"];
  unsigned int t_max = hmm_info["t_max"];
  List n_in_t = hmm_info["n_in_t"];

  List dat_out(t_max);
  unsigned int id, len, total;
  for(t = 0; t < t_max; ++t) {
    len = 0;
    List dat_info_t;
    IntegerVector idx = n_in_t[t];
    for(i = 0; i < n_t[t]; ++i) {
      id = idx[i];
      len += length[id];
    }
    IntegerVector qua_out(len);
    IntegerVector obs_out(len);
    IntegerVector obs_index_out(len);
    IntegerVector ref_pos_out(len);
    IntegerVector read_pos_out(len);
    total = 0;
    for(i = 0; i < n_t[t]; ++i) {
      id = idx[i];
      for(j = 0; j < length[id]; ++j) {
        qua_out(total) =  qua[index[id] + j];
        obs_out(total) = obs[index[id] + j];
        obs_index_out(total) = obs_index[index[id] + j];
        ref_pos_out(total) = ref_pos[index[id] + j];
        read_pos_out(total++) = read_pos[index[id] + j];
      }
    }
    dat_info_t["total"] = total;
    dat_info_t["ref_length_max"] = hap_len[t];
    dat_info_t["qua"] = qua_out;
    dat_info_t["nuc"] = obs_out;
    dat_info_t["id"] = obs_index_out;
    dat_info_t["ref_pos"] = ref_pos_out;
    dat_info_t["read_pos"] = read_pos_out;
    dat_out(t) = dat_info_t;
  }
  return(dat_out);
}

/*
 * format the data for mnlogit, note here ref_pos should be shifted according the strating t
 */
// [[Rcpp::export]]
List format_data2(List hmm_info, List d_info, List hap_info) {
  IntegerVector num_states = hmm_info["num_states"];
  List time_pos = hmm_info["time_pos"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int t, m, i, j, k;

  List subsample = sub_sample(hmm_info, d_info);

  unsigned int all = 0;
  for(t = 0; t < t_max; ++t) {
    List data_info = subsample(t);
    unsigned int tt = data_info["total"];
    all += num_states[t] * tt * NUM_CLASS;
  }
  // Rcout << all << "\n";
  IntegerMatrix full_dat(all, 6);
  unsigned int numerate = 0;
  for(t = 0; t < t_max; ++t) {
    List full_hap_t = hap_info(t);
    List data_info = subsample(t);
    unsigned int num = data_info["total"];
    unsigned int len = num * NUM_CLASS;
    // Rcout << "t :" << t << " ";
    // unsigned int new_len = num * NUM_CLASS;
    for(m = 0; m < num_states[t]; ++m) {
      IntegerMatrix haplotype = full_hap_t(m);
      IntegerVector tp = time_pos[t];
      IntegerMatrix df = format_data_simple(data_info, haplotype, tp[0]);
      // group the same record together
      for(j = 0; j < len; ++j) {
        full_dat(numerate, 0) = df(j, 0);
        full_dat(numerate, 1) = df(j, 1);
        full_dat(numerate, 2) = df(j, 2);
        full_dat(numerate, 3) = df(j, 3);
        full_dat(numerate, 4) = df(j, 4);
        full_dat(numerate++, 5) = df(j, 5);
      }
    }
  }
  // Rcout << "all_len" << numerate << " " << all << "\n";
  // hash the matrix
  List hashed_dat = hash_mat(full_dat);
  List all_id = hashed_dat["all_id"];
  IntegerVector idx =  hashed_dat["idx"];
  // make a new dataframe based on the hash table, i.e. repeat each record 4 times with changing the nuc
  IntegerMatrix subset_dat = ss(full_dat, idx);
  IntegerVector ref_pos = subset_dat(_, 1);
  IntegerVector read_pos = subset_dat(_, 2);
  IntegerVector qua = subset_dat(_, 4);
  IntegerVector obs = subset_dat(_, 3);
  IntegerVector hap_nuc = subset_dat(_, 5);
  IntegerVector obs_index = subset_dat(_, 0);

  int input_arr[] = {0, 1, 2, 3};
  unsigned int total = idx.size() * MLOGIT_CLASS;

  IntegerVector r_ref_pos(total);
  IntegerVector r_read_pos(total);
  IntegerVector r_qua(total);
  IntegerVector r_obs(total);
  IntegerVector r_hap_nuc(total);
  IntegerVector r_mode(total);
  IntegerVector r_id(total);

  for (k = 0; k < idx.size(); ++k)
    for (i = 0; i < MLOGIT_CLASS; ++i) {
      r_ref_pos[k * MLOGIT_CLASS + i] = ref_pos[k];
      r_read_pos[k * MLOGIT_CLASS + i] = read_pos[k];
      r_qua[k * MLOGIT_CLASS + i] = qua[k];
      r_id[k * MLOGIT_CLASS + i] = obs_index[k];
      r_hap_nuc[k * MLOGIT_CLASS + i] = hap_nuc[k];
      r_obs[k * MLOGIT_CLASS + i] = input_arr[i];
    }
  for (i = 0; i < idx.size(); ++i)
    for (k = 0; k < MLOGIT_CLASS; ++k)
      if (r_obs[k + MLOGIT_CLASS * i] == obs[i])
        r_mode[k + MLOGIT_CLASS * i] = 1;

  DataFrame df_new = DataFrame::create(
    Named("id") = r_id,
    Named("mode") = r_mode,
    Named("read_pos") = r_read_pos,
    Named("ref_pos") = r_ref_pos,
    Named("qua") = r_qua,
    Named("nuc") = r_obs,
    Named("hap_nuc") = r_hap_nuc);

  return List::create(_["idx"] = all_id, _["df_new"] = df_new, _["subset_dat"] = subset_dat);
  // return List::create(_["subset_dat"] = subset_dat, _["hashed_dat"] = hashed_dat, _["full_dat"] = full_dat);
  // return List::create( _["full_dat"] = full_dat);
}

NumericVector make_weight(List wic, List gamma, List hmm_info, List dat_info) {
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  IntegerVector length = dat_info["length"];
  List n_in_t = hmm_info["n_in_t"];
  unsigned int t_max = hmm_info["t_max"];

  unsigned int t, m, i, k, j;
  unsigned int id;
  unsigned int num, count = 0;
  IntegerVector n_set;
  for(t = 0; t < t_max; ++t) {
    num = 0;
    n_set = n_in_t[t];
    for (i = 0; i < n_t[t]; ++i) {
      id = n_set[i];
      num += length[id];
    }
    count += num * NUM_CLASS * num_states[t];
  }
  // Rcout << "weight length:" << count << "\n";
  NumericVector weight(count);
  count = 0;
  for(t = 0; t < t_max; ++t) {
    List w_icm = wic(t);
    NumericVector gam = gamma(t);
    n_set = n_in_t[t];
    for(m = 0; m < num_states[t]; ++m) {
      NumericMatrix w_tic = w_icm(m);
      for (i = 0; i < n_t[t]; ++i) { // todo: maybe better to make a read set for each t
        id = n_set[i];
        for(j = 0; j < length[id]; ++j)
          for (k = 0; k < NUM_CLASS; ++k)
            weight(count++) = w_tic(i, k) * exp(gam(m));
      }
    }
  }
  return(weight);
}

NumericVector comb_weight(NumericVector weight, List hash_idx) {
  int len = hash_idx.size();
  NumericVector new_weight(len);
  for (int i = 0; i < len; ++i) {
    // subset weight
    IntegerVector idx = hash_idx[i];
    NumericVector weig = weight[idx];
    for(int j = 0; j < weig.size(); ++j) {
      new_weight[i] += weig[j];
    }
  }
  return(new_weight);
}
// [[Rcpp::export]]
List baum_welch_init(List hmm_info, List data_info, List hap_info, int PD_LENGTH, List par,
                     List trans_indicator, List hash_idx, List trans_constraint,
                     int db, double penality)
{
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  NumericMatrix beta = par["beta"];
  NumericVector eta = par["eta"];
  unsigned int t_max = hmm_info["t_max"];

  List par_hmm;
  /* initialize hmm par */
  par_hmm = ini_hmm(t_max, num_states, trans_indicator);
  if (db)
    par_hmm = init_penal(par_hmm, trans_constraint, num_states, t_max, penality);
  NumericVector phi = par_hmm["phi"];
  List trans = par_hmm["trans"];
  List for_emit = compute_emit(hmm_info, data_info, hap_info, beta, eta, PD_LENGTH);
  par_hmm["emit"] = for_emit["emit"];
  List w_ic = for_emit["w_ic"];
  List par_hmm_bf = forward_backward(par_hmm, t_max, num_states);
  List gamma = par_hmm_bf["gamma"];
  /* prepare weight for beta */
  NumericVector wei = make_weight(w_ic, gamma, hmm_info, data_info);
  NumericVector weight = comb_weight(wei, hash_idx);
  //store the parmaeters for calling mnlogit
  List par_aux = List::create(
    Named("beta") = beta,
    Named("w_ic") = w_ic,
    Named("weight") = weight);

  List ls = List::create(
    Named("par_hmm_bf") = par_hmm_bf,
    Named("par_aux") = par_aux,
    Named("par_hmm") = par_hmm);

  return(ls);
}
// [[Rcpp::export]]
List baum_welch_iter(List hmm_info, List par_hmm, List data_info, List hap_info,
                     NumericMatrix beta, int PD_LENGTH, List hash_idx, int no_emi_upt)
{
  List par_aux = par_hmm["par_aux"];
  List par_hmm_bf = par_hmm["par_hmm_bf"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  List w_ic = par_aux["w_ic"];
  List gamma = par_hmm_bf["gamma"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int n_observation = data_info["n_observation"];
  List par_hmm_new;
  List wic_new;

  /* update eta */
  NumericVector eta_new = update_eta(w_ic, gamma, num_states, n_t, t_max, n_observation);

  /* update emit */
  if(no_emi_upt) {
    List par_hmm_bf = par_hmm["par_hmm"];
    par_hmm_new["emit"] = par_hmm_bf["emit"];
    wic_new = par_aux["w_ic"];
  } else {
    List for_emit_new = compute_emit(hmm_info, data_info, hap_info, beta, eta_new, PD_LENGTH);
    par_hmm_new["emit"] = for_emit_new["emit"];
    wic_new = for_emit_new["w_ic"];
  }

  /* update hmm par (except for emit) */
  NumericVector phi_new = update_phi(par_hmm_bf, num_states[0]);
  List trans_new = update_trans(par_hmm_bf, t_max, num_states);
  par_hmm_new["phi"] = phi_new;
  par_hmm_new["trans"] = trans_new;

  List par_hmm_bf_new = forward_backward(par_hmm_new, t_max, num_states);
  List gamma_new = par_hmm_bf_new["gamma"];

  /* prepare weight for beta */
  NumericVector weight;
  if(!no_emi_upt) {
    NumericVector wei = make_weight(wic_new, gamma_new, hmm_info, data_info);
    weight = comb_weight(wei, hash_idx);
  }

  List par_aux_out = List::create(
    Named("beta") = beta,
    Named("eta") = eta_new,
    Named("w_ic") = wic_new,
    Named("weight") = weight);

  List ls = List::create(
    Named("par_aux") = par_aux_out,
    Named("par_hmm") = par_hmm_new,
    Named("par_hmm_bf") = par_hmm_bf_new);

  return(ls);
}

// find which states at t can transfer to the next states at t+1
// NOTICE: EVEN IF TWO HEIGHBOURING STATES DO NOT HAVE OVERLAP! THERE ARE STILL POSSIBLY THAT THEY HAVE OVERLAP THROUGH COMMOM STATES!
// undecided_pos starts from 0; but time_pos starts from int hap_min_pos = dat_info["ref_start"];

// [[Rcpp::export]]
List trans_permit(IntegerVector num_states, List overlap_info, List combination, int t_max) {
  List trans_permits(t_max - 1);
  List loci = overlap_info["location"];
  IntegerVector overlapped_id = overlap_info["overlapped_id"];
  // List overlapped = overlap_info["overlapped"];
  IntegerVector start_t = overlap_info["start_t"];

  unsigned int t, j, m, w;
  IntegerVector start(t_max);
  if(start_t.size() > 1)
    for(t = 0; t < t_max; ++t)
      for(j = 0; j < start_t.size(); ++j)
        if(t == start_t[j]) {
          start[t] = 1;
          break;
        }

  for(t = 1; t < t_max; ++t) {
    // Rcout << t << "t " << num_states[t] << "\t" << num_states[t - 1] << "\t start\t" << start[t] << "\n";
    if(num_states[t] != 1 && num_states[t - 1] != 1 && start[t] != 1) {
      IntegerMatrix trans(num_states[t - 1], num_states[t]);
      IntegerMatrix comb_t1 = combination[t - 1];
      IntegerMatrix comb_t2 = combination[t];

      IntegerVector location_t1 = loci[t - 1];
      IntegerVector location_t2 = loci[t];
      int len = comb_t2.ncol() + comb_t1.ncol();
      IntegerVector index(len);

      // Rcout << location_t1 << "\n";
      // Rcout << location_t2 << "\n";
      // get the overlapped_id region, this might be different from the overlapped_id states we had
      int id = -1;
      for(j = 0; j < location_t1.size(); ++j)
        if(location_t1[j] == location_t2[0]) {
          id = j;
          break;
        }
      int end = id;
      for(j = id + 1; j < location_t1.size(); ++j)
        if(location_t1[j] == location_t2[j - id])
          end = j;

        if(id != -1) {
          // Rcout << t << ":overlapped_id " << id << "\t" << end << "\n";
          for(m = 0; m < num_states[t - 1]; ++m) {
            IntegerVector hap_t1 = comb_t1(m, _);
            for(w = 0; w < num_states[t]; ++w) {
              IntegerVector hap_t2 = comb_t2(w, _);
              // Rcout << hap_t1 << "|| " << hap_t2 << "\n";
              for(j = id; j < end + 1; ++j) {
                if (hap_t1[j] != hap_t2[j - id]) {
                  trans(m, w) = 1; // represents m cannot transfer to w
                  break;
                }
              }
            }
          }
          trans_permits(t - 1) = trans;
        } else {// find the transition between t2 and t0, then map to the common between t1 and t0
          IntegerVector location_t0 = loci[overlapped_id[t]];
          IntegerMatrix comb_t0 = combination[overlapped_id[t]]; // overlapped_id t
          int count = 0;
          int id_t1 = -1;
          // Rcout << t << " not overlapped with last\n";
          for(j = 0; j < location_t0.size(); ++j)
            if(location_t0[j] == location_t1[location_t1.size() - 1]) {
              index[count++] = j;
              id_t1 = j;
              break;
            }
          if(id_t1 == -1) // t1 and t2 have no connection
            continue;

          // for(j = id_t1 + 1; j < location_t0.size(); ++j)
          //   if(location_t0[j] == location_t1[j - id_t1])
          //     index[count++] = j;

          for(j = 0; j < location_t0.size(); ++j)
            if(location_t0[j] == location_t2[0]) {
              index[count++] = j;
              id = j;
              break;
            }
          for(j = id + 1; j < location_t0.size(); ++j)
            if(location_t0[j] == location_t2[j - id])
              index[count++] = j;

          index.erase(count, len);
          // Rcout << index << "\n";
          int new_len = comb_t2.ncol() + 1;
          for(m = 0; m < num_states[t - 1]; ++m) {
            IntegerVector hap_t1 = comb_t1(m, _);
            for(w = 0; w < num_states[t]; ++w) {
              int flag = 0;
              IntegerVector hap_t2 = comb_t2(w, _);
              // get this combination
              IntegerVector new_v(new_len);
              hap_t2.insert(0, hap_t1[hap_t1.size() - 1]);

              for(int i = 0; i < num_states[overlapped_id[t]]; ++i) {
                IntegerVector hap_t0 = comb_t0(i, _);
                IntegerVector tmp = hap_t0[index]; // get the sites has connection in t0
                // Rcout << tmp << "|| ";
                count = 0;
                for(j = 0; j < tmp.size(); ++j)
                  if(tmp[j] == hap_t2[j])
                    count++;
                  if(count == tmp.size()) {
                    flag = 1;
                    break;
                  }
              }
              if(!flag)
                trans(m, w) = 1;
            }
          }
          trans_permits(t - 1) = trans;
        }
    } else {
      IntegerMatrix temp(num_states[t - 1], num_states[t]);
      trans_permits(t - 1) = temp;
    }
  }
  return(trans_permits);
}

// [[Rcpp::export]]
List trans_const(List overlap_info, List combination, IntegerVector db_sites,
                 IntegerVector num_states, int t_max) {
  unsigned int t, j, m, w;
  List trans_const(t_max);
  List loci = overlap_info["location"];
  IntegerVector start_t = overlap_info["start_t"];

  IntegerVector start(t_max);
  if(start_t.size() > 1)
    for(t = 0; t < t_max; ++t)
      for(j = 0; j < start_t.size(); ++j)
        if(t == start_t[j]) {
          start[t] = 1;
          break;
        }

  unordered_set<int> db(db_sites.begin(), db_sites.end());
  if(num_states[0] == 1) {
    IntegerMatrix trans(num_states[0], num_states[1]);
    trans_const(0) = trans;
  } else {
    IntegerMatrix comb = combination[0];
    IntegerVector location = loci[0];
    IntegerVector tran(num_states[0]);
    // mark single -> db/ db -> single, should I also mark db to db?
    for(m = 0; m < comb.nrow(); ++m)
      for(j = 0; j < location.size(); ++j)
        if(db.find(location[j]) != db.end())
          if(comb(m, j) > 1) { // the index of db class
            tran[m] = 1;
            break;
          }
    trans_const[0] = tran;
  }

  for(t = 1; t < t_max; ++t) {
    // Rcout << "t " << t << "\n";

    IntegerMatrix trans(num_states[t - 1], num_states[t]);
    if(num_states[t] == 1) {
      trans_const(t) = trans;
      continue;
    }

    IntegerVector location_t2 = loci[t];
    IntegerMatrix comb_t2 = combination[t];
    if(start[t] == 1 && num_states[t - 1] == 1) {
      for(m = 0; m < comb_t2.nrow(); ++m)
        for(j = 0; j < location_t2.size(); ++j)
          if(db.find(location_t2[j]) != db.end())
            if(comb_t2(m, j) > 1) { // the index of db class
              trans(0, m) = 1;
              break;
            }
            trans_const(t) = trans;
            continue;
    }

    IntegerMatrix comb_t1 = combination[t - 1];
    IntegerVector location_t1 = loci[t - 1];
    //
    // print_intmat(comb_t1);
    // print_intmat(comb_t2);
    // Rcout << location_t1 << ", " << location_t2 << "\n";
    if((num_states[t] != 1 && num_states[t - 1] == 1) ||
       (start[t] == 1 && num_states[t - 1] != 1)) {
      for(m = 0; m < comb_t2.nrow(); ++m) {
        int flag = 0;
        for(j = 0; j < location_t2.size(); ++j) {
          if(comb_t2(m, j) > 1) {
            if(db.find(location_t2[j])!= db.end())
              for(w = 0; w < comb_t1.nrow(); ++w) {
                trans(w, m) = 1;
                flag = 1;
              }
          }
          if(flag)
            break;
        }
      }
      if(start[t] == 1) {
        for(m = 0; m < comb_t1.nrow(); ++m) {
          int flag = 0;
          for(j = 0; j < location_t1.size(); ++j) {
            if(db.find(location_t1[j])!= db.end()) {
              if(comb_t2(m, j) > 1)
                for(w = 0; w < comb_t2.nrow(); ++w) {
                  trans(m, w) = 1;
                  flag = 1;
                }
            }
            if(flag)
              break;
          }
        }
      }
    }
    trans_const(t) = trans;
  }
  return(trans_const);
}

/*
 * Find the hidden states need to be deleted, start need to be adjusted, new nums is the num of states after constraint
 */
// [[Rcpp::export]]
List find_deleted(List hmm_info, List overlap_info) {

  unsigned int t, i;
  IntegerVector overlapped_id = overlap_info["overlapped_id"];
  IntegerVector start_t = overlap_info["start_t"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector end_t(start_t.size());
  unsigned int t_max = hmm_info["t_max"];
  start_t.push_back(t_max - 1);
  // end t for each block

  int count = 0;
  for(i = 0; i < start_t.size() - 1; ++i) {
    int max = 0;
    for(t = start_t[i]; t < start_t[i + 1]; ++t) {
      if(overlapped_id[t] != -1 && num_states[t] != 1) {
        if(t > max)
          max = t;
      }
    }
    end_t[count++] = max;
  }
  start_t.erase(start_t.end() - 1);

  // Rcout << end_t << "\n";

  vector<int> delete_t;
  // delete states with hap block that has no variation

  // start_t need to be changed so as the overlapped id
  IntegerVector new_start(start_t.size());
  for(t = 0; t < start_t.size(); ++t)
    new_start[t] = start_t[t];

  for(t = 0; t < start_t.size(); ++t) {
    count = 0;
    for(i = start_t[t]; i <= end_t[t]; ++i)
      if(num_states[i] == 1) {
        count++;
        delete_t.push_back(i);
      }
      if(count) {
        for(i = t + 1; i < start_t.size(); ++i)
          new_start[i] -= count;
      }
  }

  IntegerVector new_overlapped_id(overlapped_id.size());
  if(!delete_t.empty()) {
    for(t = 0; t < overlapped_id.size(); ++t)
      new_overlapped_id[t] = overlapped_id[t];
    for(t = 0; t < delete_t.size(); ++t)
      for(i = delete_t[0] + 1; i < overlapped_id.size(); ++i)
        if(overlapped_id[i] != -1)
          if(delete_t[t] < i)
            new_overlapped_id[i]--;
  }
  // Rcout << "overlapped_id " << overlapped_id << "\n";
  // Rcout << "new_start " << new_start << "\n";
  List del = List::create(
    Named("delete_t") = delete_t,
    Named("overlapped_id") = new_overlapped_id,
    Named("new_start") = new_start);

  return(del);
}

// if do not use merge, then use this

List remove_non_variant(List hmm_info, List overlap_info) {
  unsigned int t, i;
  unsigned int t_max = hmm_info["t_max"];
  IntegerVector overlapped_id = overlap_info["overlapped_id"];
  List location = overlap_info["location"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector start_t = overlap_info["start_t"];
  // some new states is 1, and it lies after non-one num of states state
  vector<int> a;

  for(t = 1; t < t_max; ++t) {
    IntegerVector cur = location[t];
    if(cur[0] == -1)
      continue;
    if(num_states[t] == 1)
      a.push_back(t);
  }
  int length = 1;
  vector<int> targets;

  // Traverse the array from first position
  for(i = 1; i < a.size(); i++) {
    if (i == a.size() -1 || a[i] - a[i - 1] != 1) {
      if (length == 1) {
        // Rcout << a[i - length] << "\n";
        // Rcout << num_states[a[i - length] - 1] << "\n";
        // Rcout << num_states[a[i - length] + 1] << "\n";
        int l = a[i - length] - 1;
        int r = a[i - length] + 1;
        // states cannot be 1 and they need to have overlap
        if(num_states[l] != 1 && num_states[r] != 1) {
          IntegerVector left = location[l];
          IntegerVector right = location[r];
          IntegerVector inter = intersect(left, right);
          if(inter.size() != 0)
            targets.push_back(a[i - length]);
        }
      }
      length = 1;
    } else {
      length++;
    }
  }

  // if one location is removed then remomber to remove this in hmm info
  if(!targets.empty())
    for(t = 0; t < targets.size(); ++t) {
      for(i = targets[0] + 1; i < overlapped_id.size(); ++i)
        if(overlapped_id[i] != -1)
          if(targets[t] < i)
            overlapped_id[i]--;
          for(i = 0; i < start_t.size(); ++i)
            if(targets[t] < start_t[i])
              start_t[i]--;
    }

    List del = List::create(
      Named("delete_t") = targets,
      Named("overlapped_id") = overlapped_id,
      Named("start_t") = start_t);

  return(del);
}

/*
 * Find the hidden states satisfies: A: 1, 2 , B: 3, C:1, 2, and B has no connection with A and C
 * merge C to A,
 */
// [[Rcpp::export]]
List merge_no_connection(List hmm_info, List overlap_info, List full_hap) {

  unsigned int t, i;
  unsigned int t_max = hmm_info["t_max"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector overlapped_id = overlap_info["overlapped_id"];
  IntegerVector start_t = overlap_info["start_t"];
  List location = overlap_info["location"];
  List n_in_t = hmm_info["n_in_t"];
  IntegerVector n_t = hmm_info["n_t"];

  List new_time_pos(t_max);
  vector<int> targets;
  unordered_map<string, vector<int>> seen_map; // record the same location and their index
  for(t = 0; t < t_max; ++t) {
    new_time_pos[t] = time_pos[t];
    IntegerVector cur = location[t];
    if(cur[0] == -1)
      continue;
    string s = "";
    for(i = 0; i < cur.size(); ++i)
      s += to_string(cur[i]);
    seen_map[s].push_back(t);
  }
  map<int, vector<int>> ordered_id;
  // put the first loci index: t, and the other t have the same loci in it
  for(auto &m: seen_map) {
    vector<int> id = m.second;
    ordered_id[id[0]] = id;
  }

  for(auto l = ordered_id.begin(); l != ordered_id.end(); ++l) {
    vector<int> id = l->second;
    if(id[0] >= start_t[start_t.size() - 1])
      continue;
    // the index should be in the same hap block, i.e. start_t1 < id < start_t2
    set<int> tmp; // this tmp also include the B
    int which_t = 0;
    for(t = 0; t < start_t.size() - 1; ++t) // can use binary search??
      if(id[0] >= start_t[t] && start_t[t + 1] > id[0]) {
        which_t = start_t[t + 1];
        break;
      }
    // Rcout << "next: " << which_t << "\n";
    for(i = id.size() - 1; i >= 1; --i)
      if(id[i] >= which_t)
        tmp.insert(id[i]);

    if(!tmp.empty()) {
      // merge reads in time t
      // for(auto &m:tmp)
      //   Rcout << m << " ";
      // Rcout << "are merged to " << id[0] << "\n";
      int merge_to = id[0];
      int max_len = p_tmax[merge_to] + time_pos[merge_to];
      int max_id = 0;
      int start_pos = time_pos[merge_to];
      vector<int> rid = as<vector<int>>(n_in_t[merge_to]);
      IntegerVector all_pos;
      for(i = 0; i < rid.size(); ++i)
        all_pos.push_back(start_pos);
      for(auto &m:tmp) {
        int bound = p_tmax[m] + time_pos[m];
        if(bound > max_len) {
          max_id = m;
          max_len = bound;
        }
        vector<int> rid2 = as<vector<int>>(n_in_t[m]);
        int len = rid2.size();
        n_t[merge_to] += len;
        rid.insert(rid.end(), rid2.begin(), rid2.end());
        for(i = 0; i < len; ++i)
          all_pos.push_back(time_pos[m]);
      }
      if(max_id) { // merge haplotypes
        List hap_ori = full_hap[merge_to];
        List hap_new = full_hap[max_id];
        IntegerMatrix extra = hap_new[0];
        int tail = max_len - (p_tmax[merge_to] + time_pos[merge_to]);
        // Rcout << "extra: " << max_len << "-" << p_tmax[merge_to] << "-" << time_pos[merge_to] << "\n";
        IntegerMatrix ex = extra(_, Range(extra.ncol() - tail, extra.ncol() - 1));
        arma::Mat<int> m1 = as<arma::Mat<int>>(ex);
        for(i = 0; i < hap_ori.size(); ++i) {
          arma::Mat<int> m2 = as<arma::Mat<int>>(hap_ori[i]);
          arma::Mat<int> out = join_rows(m2, m1);
          hap_ori[i] = wrap(out);
        }
        full_hap[merge_to] = hap_ori;
      }

      new_time_pos[merge_to] = all_pos;
      n_in_t[merge_to] = rid;
      p_tmax[merge_to] = max_len - time_pos[merge_to];
      // Rcout << "merge_to: " << merge_to << ": max len: " << p_tmax[merge_to] << "\n";
      // print_c_intvec(rid);
      targets.insert(targets.end(), tmp.begin(), tmp.end());
    }
  }

  if(!targets.empty())
    for(t = 0; t < targets.size(); ++t) {
      for(i = targets[0] + 1; i < overlapped_id.size(); ++i)
        if(overlapped_id[i] != -1)
          if(targets[t] < i)
            overlapped_id[i]--;
          for(i = 0; i < start_t.size(); ++i)
            if(targets[t] < start_t[i])
              start_t[i]--;
    }

  List ls = List::create(
    Named("targets") = targets,
    Named("n_in_t") = n_in_t,
    Named("n_t") = n_t,
    Named("p_tmax") = p_tmax,
    Named("new_time_pos") = new_time_pos,
    Named("full_hap") = full_hap,
    Named("start_t") = start_t,
    Named("overlapped_id") = overlapped_id);

  return(ls);
}

// merge states appeared in another hap block which should be in the previous one
// [[Rcpp::export]]
List merge_states(List hmm_info, List overlap_info, List full_hap) {

  unsigned int t, i, j, k;
  unsigned int t_max = hmm_info["t_max"];
  List time_po = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector overlapped_id = overlap_info["overlapped_id"];
  IntegerVector start_t = overlap_info["start_t"];
  List loci = overlap_info["location"];
  List n_in_t = hmm_info["n_in_t"];
  IntegerVector n_t = hmm_info["n_t"];
  List new_time_pos = time_po;

  vector<int> targets;
  start_t.push_back(t_max);
  vector<std::set<int>> index_set;
  for(t = 0; t < start_t.size() - 1; ++t) {
    std::set<int> index_set1;
    int s = start_t[t];
    int e = start_t[t + 1];
    for(j = s; j < e; ++j) {
      IntegerVector location = loci[j];
      if(location[0] == -1)
        continue;
      for(k = 0; k < location.size(); ++k)
        index_set1.insert(location[k]);
    }
    index_set.push_back(index_set1);
  }
  for(t = 0; t < index_set.size() - 1; ++t) {
    std::set<int> index_set1 = index_set[t];
    std::set<int> index_set2 = index_set[t + 1];
    for(int id:index_set2) {
      if(index_set1.find(id) != index_set1.end()) {// merge this one
        Rcout << "merge " << id << " ";
        int s = start_t[t + 1];
        int e = start_t[t + 2];
        for(j = s; j < e; ++j) {
          IntegerVector location = loci[j];
          if(location[0] == -1)
            continue;
          vector<int> loc = as<vector<int>>(location);
          if(find(loc.begin(), loc.end(), id) != loc.end()) {
            int s1 = start_t[t];
            int e1 = start_t[t + 1];
            for(i = e1 - 1; i >= s1; --i) {
              location = loci[i];
              if(location[0] == -1)
                continue;
              vector<int> loc0 = as<vector<int>>(location);
              if(find(loc0.begin(), loc0.end(), id) != loc0.end()) {
                Rcout << "from " << j << " to " << i << "\n";
                int merge_to = i;
                IntegerVector time_pos = time_po[merge_to];
                int start_pos = time_pos[0];
                IntegerVector pos = new_time_pos[merge_to];
                IntegerVector all_pos;
                vector<int> rid2 = as<vector<int>>(n_in_t[j]);
                vector<int> rid = as<vector<int>>(n_in_t[merge_to]);
                if(pos.size() == 1) { // meaning this has not been merged
                  if(time_pos.size() == 1)
                    for(i = 0; i < rid.size(); ++i)
                      all_pos.push_back(start_pos);
                  else
                    all_pos = time_pos;
                } else
                  all_pos = pos;
                IntegerVector time_pos_from = time_po[j];
                if(time_pos_from.size() == 1)
                  for(i = 0; i < rid2.size(); ++i)
                    all_pos.push_back(time_pos_from[0]);
                else
                  for(i = 0; i < time_pos_from.size(); ++i)
                    all_pos.push_back(time_pos_from[i]);
                rid.insert(rid.end(), rid2.begin(), rid2.end());

                List hap_ori = full_hap[merge_to];
                List hap_new = full_hap[j];
                IntegerMatrix extra = hap_new[0];
                int max_len = p_tmax[j] + time_pos_from[0];
                int tail = max_len - (p_tmax[merge_to] + start_pos);
                IntegerMatrix ex = extra(_, Range(extra.ncol() - tail, extra.ncol() - 1));
                arma::Mat<int> m1 = as<arma::Mat<int>>(ex);
                for(i = 0; i < hap_ori.size(); ++i) {
                  arma::Mat<int> m2 = as<arma::Mat<int>>(hap_ori[i]);
                  arma::Mat<int> out = join_rows(m2, m1);
                  hap_ori[i] = wrap(out);
                }
                full_hap[merge_to] = hap_ori;

                new_time_pos[merge_to] = all_pos;
                n_in_t[merge_to] = rid;
                n_t[merge_to] = rid.size();
                p_tmax[merge_to] = max_len - start_pos;
                break;
              }
            }
            targets.push_back(j);
          }
        }
      }
    }
  }

  start_t.erase(start_t.size() - 1);
  if(!targets.empty())
    for(t = 0; t < targets.size(); ++t) {
      for(i = targets[0] + 1; i < overlapped_id.size(); ++i)
        if(overlapped_id[i] != -1)
          if(targets[t] < i)
            overlapped_id[i]--;
          for(i = 0; i < start_t.size(); ++i)
            if(targets[t] < start_t[i])
              start_t[i]--;
    }

  List ls = List::create(
    Named("targets") = targets,
    Named("n_in_t") = n_in_t,
    Named("n_t") = n_t,
    Named("p_tmax") = p_tmax,
    Named("start_t") = start_t,
    Named("new_time_pos") = new_time_pos,
    Named("overlapped_id") = overlapped_id,
    Named("full_hap") = full_hap);

  return(ls);
}
