#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "hmm_state.h"
#include "utils.h"

#define NUM_CLASS 4

using namespace Rcpp;
using namespace std;

void sampling_link(vector<int> &res, int in_id, vector<vector<int> > reads_pool,
                   vector<vector<int> > start_end, int last_pos, int flag);
void cat_reads_forward(vector<int> &partial, vector<int> res, vector<vector<int> > start_end,
                       int in_id, int s, IntegerMatrix uni_read);
void cat_reads_backward(vector<int> &partial, vector<int> res, vector<vector<int> > start_end,
                        int in_id, int e, IntegerMatrix uni_read);

void sampling_link(vector<int> &res, int in_id, vector<vector<int> > reads_pool,
                   vector<vector<int> > start_end, int last_pos, int flag) {
  res.push_back(in_id);
  int pos = start_end[in_id][flag];
  if(pos != last_pos) {// keeps sampling until reach 0
    // reads_pool need to be updated
    vector<int> sub_reads = reads_pool[in_id];
    IntegerVector seeds = sample_int(wrap(sub_reads), 1, true);
    int sub_choosed = seeds[0];
    sampling_link(res, sub_choosed, reads_pool, start_end, last_pos, flag);
  }
}
// sample to the left
void cat_reads_forward(vector<int> &partial, vector<int> res, vector<vector<int> > start_end,
                       int in_id, int s, vector<vector<int> > uni_read) {
  if(s != 0) {
    int subs = start_end[res[in_id]][0];
    vector<int> wr = uni_read[res[in_id]];
    //subset the read [begin of current read and the last read]
    vector<int> read = vector<int>(wr.begin() + subs, wr.begin() + s);
    for(int m = read.size() - 1; m >= 0; --m)
      partial.insert(partial.begin(), read[m]);
    cat_reads_forward(partial, res, start_end, in_id + 1, subs, uni_read);
  }
}
// sample to the right
void cat_reads_backward(vector<int> &partial, vector<int> res, vector<vector<int> > start_end,
                        int in_id, int e, vector<vector<int> > uni_read) {
  if(in_id != res.size()) {
    int sube = start_end[res[in_id]][1];
    vector<int> wr = uni_read[res[in_id]];
    // Rcout << in_id<< " " << res.size() << "end:" << e + 1 << " " << sube << "\n";
    //subset the read [end of last read and the current read]
    vector<int> read = vector<int>(wr.begin() + e + 1, wr.begin() + sube + 1);
    for(int m = 0; m < read.size(); ++m)
      partial.push_back(read[m]);
    cat_reads_backward(partial, res, start_end, in_id + 1, sube, uni_read);
  }
}


// we go to either or both ends of this read and randomly sample any read that
// covers the next (missing?) site
// We discard any part of the new read overlapping the already sampled part

// thoughts: some unlikely nuc should be treated as missing? Since in simulation
// it is likely we pick such sample
// and simulate haplotype block instead of two consuctive haplotype, then determine
// the Ht???
// and how to choose from 4 haplotypes??? link_in should be sites at least more than 4 reads connect
IntegerMatrix sample_hap(vector<string> map, IntegerMatrix uni_read,
                         List all_id, int sample) {
  //remove non-covered reads
  int i, j, k;
  int num = uni_read.ncol();
  int n_reads = uni_read.nrow();
  vector<int> n_sreads(n_reads);
  int count = 0;
  for(k = 0; k < n_reads; ++k) {
    IntegerVector same_reads = all_id[k];
    n_sreads[k] = same_reads.size();
    count += n_sreads[k];
  }
  IntegerMatrix hap_pool;
  vector<vector<int> > samples_hap;

  if(num == 2) {
    unordered_map<string, vector<int>> unique_count; // store the fully connected reads and counts
    for(k = 0; k < n_reads; ++k)
      if(uni_read(k, 1) != -1 && uni_read(k, 0) != -1) {
        string s = "";
        s = to_string(uni_read(k, 0)) + to_string(uni_read(k, 1));
        unique_count[s].push_back(k);
      }
      // find which read covers which loci
      unordered_map<int, vector<int>> loci_nuc;
      for(k = 0; k < n_reads; ++k) {
        int loci = 1;
        if(uni_read(k, 1) == -1)
          loci = 0;
        loci_nuc[loci].push_back(uni_read(k, loci));
      }

      int iter = 0;
      while(iter < sample) {
        for(k = 0; k < n_reads; ++k) { // change this to sample a k
          int loci;
          if(uni_read(k, 1) == -1)
            loci = 0;
          else if(uni_read(k, 0) == -1)
            loci = 1;
          else {
            for(j = 0; j < n_sreads[k]; ++j) {
              IntegerVector tp = uni_read(k, _);
              vector<int> tmp = as<vector<int>>(tp);
              samples_hap.push_back(tmp);
            }
            continue;
          }
          vector<int> pool;
          // Rcout << "read:" << loci << "|"<< uni_read(k, 0) << uni_read(k, 1) << "\n";
          for(auto &c:unique_count) {
            string s = c.first;
            int nuc = s[loci] - '0';
            if(nuc == uni_read(k, loci)) {
              vector<int> tmp = c.second;
              pool.insert(pool.end(), tmp.begin(), tmp.end());
            }
          }
          // if no connection, can be anything
          int flag = 0;
          if(pool.empty()) {
            pool = loci_nuc[!loci];
            flag = 1;
          }
          // print_c_intvec(pool);
          int n_same_reads = n_sreads[k];
          // Rcout << "choosed: ";
          IntegerVector ret = sample_int(wrap(pool), n_same_reads, true);
          // Rcout << ret << "\n";
          for(j = 0; j < n_sreads[k]; ++j) {
            if(flag) {
              if(loci)
                samples_hap.push_back({ret[j], uni_read(k, loci)});
              else
                samples_hap.push_back({uni_read(k, loci), ret[j]});
            }
            else{
              IntegerVector tp = uni_read(ret[j], _);
              samples_hap.push_back(as<vector<int>>(tp));
            }
          }
        }
        iter++;
      }
  } else {
    // store the start & end position of each unique read
    vector<vector<int> > start_end;
    // get the start and end position first
    // deal with the missing later, find reads cover the left/right of a read
    for(k = 0; k < n_reads; ++k) {
      vector<int> se = {0, num - 1};
      if(uni_read(k, 0) == -1) // left missing, find where start
        for(i = 1; i < num; ++i)
          if(uni_read(k, i) != -1) {
            se[0] = i;
            break;
          }
          if(uni_read(k, num - 1) == -1)
            for(i = num - 2; i >= 0; --i)
              if(uni_read(k, i) != -1) {
                se[1] = i;
                break;
              }
              start_end.push_back(se);
              // Rcout << start_end[k][0] << " " << start_end[k][1] <<"\n";
    }
    vector<int> excluded(n_reads, 0);

    // if there is deletion (or caused by unwanted nuc), then sample the deletions first
    for(k = 0; k < n_reads; ++k) {
      int redo = 0;
      for(i = 0; i < num; ++i)
        if(i > start_end[k][0] && uni_read(k, i) == -1 && i < start_end[k][1]) {
          redo = 1;
          for(j = 0; j < n_reads; ++j) {
            if(j == k || uni_read(j, i) == -1)
              continue;
            // check if th overlaps are the same
            int flag = 1;
            int index = i - 1;
            while(index >= 0) {
              if(uni_read(k, index) != uni_read(j, index)) {
                flag = 0;
                break;
              }
              index--;
            }
            if(flag) {
              // Rcout << "Replace missing:before" << uni_read(k, i) << "after" << uni_read(j, i) << "\n";
              uni_read(k, i) = uni_read(j, i);
              break;
            }
          }
        }
        // check if this read already exits
        if(redo) {
          string s = "";
          for(int j = 0; j < num; ++j)
            s += to_string(uni_read(k, j));
          int find = 0;
          for (auto &itr: map) {
            if(s == itr) {
              excluded[k] = 1;
              // Rcout << "find" << find << " " << n_sreads[find] << "\n";
              n_sreads[find]++;
              break;
            }
            find++;
          }
        }
    }

    vector<vector<int> > unique_reads;
    for(k = 0; k < n_reads; ++k)
      if(!excluded[k]) {
        // start_end2.push_back(start_end[k]);
        IntegerVector tp = uni_read(k, _);
        vector<int> tmp = as<vector<int>>(tp);
        // print_c_intvec(tmp);
        unique_reads.push_back(tmp);
      }
      for(k = n_reads - 1; k >= 0; --k) {
        if(excluded[k]) {
          start_end.erase(start_end.begin() + k);
          n_sreads.erase(n_sreads.begin() +k);
        }
      }
      n_reads = unique_reads.size();
      unordered_map<int, vector<int> > pos_covers;
      vector<vector<int> > l_overlap;
      l_overlap.resize(n_reads);
      vector<vector<int> > r_overlap;
      r_overlap.resize(n_reads);
      // find the reads covers at each position
      for(i = 0; i < num; ++i)
        for(k = 0; k < n_reads; ++k)
          if(unique_reads[k][i] != -1)
            pos_covers[i].push_back(k);

          // IntegerVector exclude(n_reads);
          // record the first missing pos for all reads, and record if a read only has one nuc
          for(k = 0; k < n_reads; ++k) {
            vector<int> se = start_end[k];
            // Rcout << "K:" << k << "| ";
            // print_c_intvec(unique_reads[k]);
            // Rcout << "\n";
            if(se[0] != 0) {
              vector<int> same_cover = pos_covers[se[0]];
              vector<int> left_cover = pos_covers[se[0] - 1];
              // get the mutual reads of these two sets
              vector<int> v_intersection;
              set_intersection(same_cover.begin(), same_cover.end(),
                               left_cover.begin(), left_cover.end(),
                               back_inserter(v_intersection));
              // Rcout << "l_intersection: \n";
              // print_c_intvec(v_intersection);
              // make sure the reads in this intersection has the same context as k,
              // and make sure this context has to be longer than 1
              // if none of them have the same context, then include all
              // give more weight to the reads contains more common nuc
              vector<int> record;
              vector<int> weight(v_intersection.size(), 1);
              int k_end = se[1];
              int flag = 0;
              while(k_end >= se[0]) {
                for(j = 0; j < v_intersection.size(); ++j) {
                  int cnt = 0;
                  vector<int> sub_se = start_end[v_intersection[j]];
                  int end = min(sub_se[1], k_end);
                  for(i = se[0]; i < end + 1; ++i) {
                    // Rcout << v_intersection[j] << ":"<< unique_reads[v_intersection[j]][i] << "\n";
                    if((unique_reads[k][i] != unique_reads[v_intersection[j]][i]) && (unique_reads[k][i] != -1)
                         && (unique_reads[v_intersection[j]][i] != -1)) {
                      // Rcout << "exclude\n";
                      record.push_back(j);
                      break;
                    }
                    cnt++;
                  }
                  weight[j] = cnt;
                }
                if(record.size() == v_intersection.size()) {
                  record.clear();
                  k_end--;
                  flag = 1;
                  // Rcout << "all exclude" << se[0] << "," << k_end <<"\n";
                } else
                  break;
              }
              // if all been excluded, then require a shorter linkage, not the entire overlap, but only a portion
              // if all excluded then give equal weight (no overlapped at all)
              if(record.empty() && flag) {
                // Rcout << "all exclude\n";
                for(i = 0; i < left_cover.size(); ++i)
                  l_overlap[k].push_back(left_cover[i]);
                int nuc;
                int next = se[0] - 1; // next has to have equal likely different nuc?
                unordered_map<int, int> nuc_cnt;
                for(j = 0; j < left_cover.size(); ++j) {
                  nuc = unique_reads[left_cover[j]][next];
                  nuc_cnt[nuc]++;
                }
                for(auto &m:nuc_cnt) {
                  if(m.first != nuc) {
                    int sum1 = m.second;
                    int sum2 = nuc_cnt[nuc];
                    // Rcout << "exclude all, " << m.first << ":" << nuc << "," << sum1 << ":" << sum2 << "\n";
                    int diff;
                    int marker;
                    if(sum1 > sum2){
                      diff = sum1 - sum2;
                      marker = nuc;
                    } else {
                      marker = m.first;
                      diff = sum2 - sum1;
                    }
                    // if(diff > 0)
                    //   Rcout << "balance the connection\n";
                    while (diff > 0) {
                      for(i = 0; i < left_cover.size(); ++i)
                        if(unique_reads[left_cover[i]][next] == marker && diff > 0) {
                          diff--;
                          l_overlap[k].push_back(left_cover[i]);
                        }
                    }
                  }
                }
              }
              if(record.size() != v_intersection.size() && !record.empty()) {
                for(i = record.size() - 1; i >= 0; --i) {
                  v_intersection[record[i]] = v_intersection.back();
                  v_intersection.pop_back();
                  weight[record[i]] = weight.back();
                  weight.pop_back();
                }
              }
              // Rcout << "weight ";
              // print_c_intvec(weight);
              for(i = 0; i < v_intersection.size(); ++i)
                for(j = 0; j < weight[i]; ++j)
                  l_overlap[k].push_back(v_intersection[i]);
              // Rcout << "left overlapped reads\n";
              // for(j = 0; j < l_overlap[k].size(); ++j)
              //   Rcout << l_overlap[k][j] << " ";
              // Rcout << "\n";
            }
            if(se[1] != num - 1) {
              vector<int> same_cover = pos_covers[se[1]];
              vector<int> left_cover = pos_covers[se[1] + 1];
              // get the mutual reads of these two sets
              vector<int> v_intersection;
              set_intersection(same_cover.begin(), same_cover.end(),
                               left_cover.begin(), left_cover.end(),
                               back_inserter(v_intersection));
              // Rcout << "r_intersection: \n";
              // print_c_intvec(left_cover);
              vector<int> record;
              vector<int> weight(v_intersection.size(), 1);
              int k_start = se[0];
              int flag = 0;
              while(k_start <= se[1]) {
                // Rcout << k_start << ":\n";
                for(j = 0; j < v_intersection.size(); ++j) {
                  // Rcout << v_intersection[j] << ":\n";
                  int cnt = 0;
                  vector<int> sub_se = start_end[v_intersection[j]];
                  int end = max(sub_se[0], k_start);
                  // Rcout << v_intersection[j] << ": " << end << "," << se[1] + 1 << "\n";
                  for(i = end; i < se[1] + 1; ++i) {
                    if((unique_reads[k][i] != unique_reads[v_intersection[j]][i]) && (unique_reads[k][i] != -1)
                         && (unique_reads[v_intersection[j]][i] != -1)) {
                      record.push_back(j);
                      break;
                    }
                    cnt++;
                  }
                  weight[j] = cnt;
                }
                if(record.size() == v_intersection.size()){
                  record.clear();
                  k_start++;
                  flag = 1;
                  // Rcout << "all exclude" << k_start << "," << se[1] <<"\n";
                } else
                  break;
              }
              if(record.empty() && flag) {
                // Rcout << "all exclude\n";
                for(i = 0; i < left_cover.size(); ++i)
                  r_overlap[k].push_back(left_cover[i]);
                int nuc;
                int next = se[1] + 1; // next has to have equal likely different nuc?
                unordered_map<int, int> nuc_cnt;
                for(j = 0; j < left_cover.size(); ++j) {
                  nuc = unique_reads[left_cover[j]][next];
                  nuc_cnt[nuc]++;
                }
                for(auto &m:nuc_cnt) {
                  if(m.first != nuc) {
                    int sum1 = m.second;
                    int sum2 = nuc_cnt[nuc];
                    // Rcout << "exclude all, "  << m.first << ":" << nuc << "," << sum1 << ":" << sum2 << "\n";
                    int diff;
                    int marker;
                    if(sum1 > sum2){
                      diff = sum1 - sum2;
                      marker = nuc;
                    } else {
                      marker = m.first;
                      diff = sum2 - sum1;
                    }
                    // if(diff > 0)
                    //   Rcout << "balance the connection\n";
                    while (diff > 0) {
                      for(i = 0; i < left_cover.size(); ++i)
                        if(unique_reads[left_cover[i]][next] == marker && diff > 0) {
                          diff--;
                          r_overlap[k].push_back(left_cover[i]);
                        }
                    }
                  }
                }
              }
              // int flag = 0;
              if(record.size() != v_intersection.size() && !record.empty()) { // remove the ones overlap only 1 and overlaps do not consistent
                for(i = record.size() - 1; i >= 0; --i) {
                  v_intersection[record[i]] = v_intersection.back();
                  v_intersection.pop_back();
                  weight[record[i]] = weight.back();
                  weight.pop_back();
                }
              }
              // Rcout << "weight ";
              // print_c_intvec(weight);
              for(i = 0; i < v_intersection.size(); ++i)
                for(j = 0; j < weight[i]; ++j)
                  r_overlap[k].push_back(v_intersection[i]);
              // Rcout << "right overlapped reads\n";
              // for(i = 0; i <  r_overlap[k].size(); ++i)
              //   Rcout <<  r_overlap[k][i] << " ";
              // Rcout << "\n";
            }
          }
          // for each read, start sampling
          int iter = 0;
          while(iter < sample) {
            // Rcout << "iter:" << iter << "\n";
            for(k = 0; k < n_reads; ++k) {
              vector<int> exists_full = unique_reads[k];
              // if(exclude[k])
              //   continue;
              // Rcout << "K:" << k << "| ";
              // print_c_intvec(exists_full);
              vector<int> l_reads = l_overlap[k];
              vector<int> r_reads = r_overlap[k];
              int n_same_reads = n_sreads[k];
              vector<vector<int> > sample_reads_r(n_same_reads);
              vector<vector<int> > sample_reads_l(n_same_reads);
              NumericVector indicator = Rcpp::rbinom(n_same_reads, 1, 0.5);
              if(!l_reads.empty()) {
                int s = start_end[k][0];
                // random sample n_same_reads from the pool (need to be adjusted to the reads share same haps)
                IntegerVector ret = sample_int(wrap(l_reads), n_same_reads, true);
                vector<int> choosed = as<vector<int>>(ret);
                // Rcout << ret << "\n";
                // then sample of them until the end
                for(i = 0; i < n_same_reads; ++i) {
                  if(indicator[i] == 0)
                    continue;
                  int in_id = choosed[i];
                  vector<int> res;
                  sampling_link(res, in_id, l_overlap, start_end, 0, 0);
                  // Rcout << "l choosed:\n";
                  // for(int m = 0; m < res.size(); ++m)
                  //   Rcout << res[m] << " ";
                  // Rcout << "\n";
                  vector<int> partial;
                  if(!res.empty()) { // cat the reads
                    // Rcout << "filling reads\n";
                    cat_reads_forward(partial, res, start_end, 0, s, unique_reads);
                    // for(int m = 0; m < partial.size(); ++m)
                    //   Rcout << partial[m] << " ";
                    // Rcout << "\n";
                  }
                  for(int m = 0; m < partial.size(); ++m)
                    sample_reads_l[i].push_back(partial[m]);
                }
              }
              if(!r_reads.empty()) {
                int e = start_end[k][1];
                IntegerVector ret = sample_int(wrap(r_reads), n_same_reads, true);
                vector<int> choosed = as<std::vector<int>>(ret);
                // Rcout << ret << "\n";
                // then sample of them until the end
                for(i = 0; i < n_same_reads; ++i) {
                  if(indicator[i] == 0)
                    continue;
                  int in_id = choosed[i];
                  vector<int> res;
                  sampling_link(res, in_id, r_overlap, start_end, num - 1, 1);
                  // Rcout << "r choosed:\n";
                  // for(int m = 0; m < res.size(); ++m)
                  //   Rcout << res[m] << " ";
                  // Rcout << "\n";
                  vector<int> partial;
                  if(!res.empty()) { // cat the reads
                    // Rcout << "filling reads\n";
                    cat_reads_backward(partial, res, start_end, 0, e, unique_reads);
                    // for(int m = 0; m < partial.size(); ++m)
                    //   Rcout << partial[m] << " ";
                    // Rcout << "\n";
                  }
                  for(int m = 0; m < partial.size(); ++m)
                    sample_reads_r[i].push_back(partial[m]);
                }
              }
              // formulate the fully connected sample
              for(i = 0; i < n_same_reads; ++i) {
                if(indicator[i] == 0)
                  continue;
                vector<int> exists = vector<int>(exists_full.begin() + start_end[k][0], exists_full.begin() + start_end[k][1] + 1);
                if(!l_reads.empty()) {
                  vector<int> partial = sample_reads_l[i];
                  exists.insert(exists.begin(), partial.begin(), partial.end());
                }
                if(!r_reads.empty()) {
                  vector<int> partial = sample_reads_r[i];
                  exists.insert(exists.end(), partial.begin(), partial.end());
                }
                // Rcout << "whole reads\n";
                // for(int m = 0; m < exists.size(); ++m)
                //   Rcout << exists[m] << " ";
                // Rcout << "\n";
                samples_hap.push_back(exists);
              }
            }
            iter++;
          }
  }
  hap_pool = IntegerMatrix(samples_hap.size(), num);
  for(i = 0; i < samples_hap.size(); ++i)
    for(j = 0; j < num; ++j)
      hap_pool(i, j) = samples_hap[i][j];
  return hap_pool;
}


IntegerMatrix determine_hap(IntegerMatrix hap_pool, double lower_ab, IntegerMatrix linkage) {
  // hash the hap_pool get at most top 4 haplotypes
  int i, j, k;
  IntegerMatrix out_hap;
  List final_info = hash_mat(hap_pool);
  IntegerVector hap_id = final_info["idx"];
  List hap_ids = final_info["all_id"];
  IntegerVector abundance(hap_id.size());
  vector<int> hap_abun;

  for(i = 0; i < hap_ids.size(); ++i) {
    IntegerVector tp = hap_ids[i];
    abundance[i] = tp.size();
  }
  // Rcout << "in_hap\n";
  // IntegerMatrix tmp = ss(hap_pool, hap_id);
  // print_intmat(tmp);
  // Rcout << "all hap abundance: " << abundance << "\n";
  if(hap_id.size() == 2) {
    // Rcout << "abundance: ";
    // Rcout << abundance << "\n";
    out_hap = ss(hap_pool, hap_id);
    // print_intmat(out_hap);
  } else {
    vector<int> ab_id;
    int prop = 0;
    IntegerVector sorted = clone(abundance).sort(true);
    // Rcout << "all hap abundance: " << sorted << "\n";
    vector<int> order(sorted.size(), 0);
    for (i = 0 ; i < order.size(); ++i)
      order[i] = i;
    sort(order.begin(), order.end(),
         [&](const int& a, const int& b) {
           return (abundance[a] > abundance[b]);
         });
    int end;
    if(sorted.size() <= 4) { // if less than 4 haplotypes
      end = sorted.size();
    } else
      end = sorted[4] == sorted[3] ? 5:4;

    vector<int> new_top_id;
    double total = 0;
    for(i = 0; i < end; ++i)
      total += sorted[i];
    int new_end = 0;
    for(i = 0; i < end; ++i) {
      new_top_id.push_back(hap_id[order[i]]);
      hap_abun.push_back(sorted[i]);
      prop += sorted[i];
      new_end++;
      if(prop/total >= lower_ab)
        break;
    }
    for(i = 0; i < new_end; ++i)
      ab_id.push_back(i);
    // should check if the kept haplotypes have all the nucs
    out_hap = ss(hap_pool, wrap(new_top_id));
    // Rcout << "candidate\n";
    // print_intmat(out_hap);
    for(j = 0; j < hap_pool.ncol(); ++j) {
      IntegerVector col_nuc = out_hap(_, j);
      IntegerVector nucs = unique(col_nuc);
      // Rcout << "unique at " << j << ": " << nucs << "\n";
      if(nucs.size() == 1)
        // search from end until find another nuc, add it to top_id
        for(i = new_end; i < order.size(); ++i)
          if(hap_pool(hap_id[order[i]], j) != nucs[0]) {
            // Rcout << "need to add more hap\n";
            // Rcout << "nuc? " << order[i] << ": "<< hap_pool(hap_id[order[i]], j) << "\n";
            int flag = 0;
            for(k = ab_id.size() - 1; k >= new_end; --k)
              if(ab_id[k] == i) {
                flag = 1;
                break;
              }
              if(!flag)
                ab_id.push_back(i);
              break;
          }
    }
    if(ab_id.size() > new_end)// meaning add more haplotypes
      for(i = ab_id.size() - 1; i >= new_end; --i) {
        hap_abun.push_back(sorted[ab_id[i]]);
        new_top_id.push_back(hap_id[order[ab_id[i]]]);
      }
      out_hap = ss(hap_pool, wrap(new_top_id));
    // check if hap include all appeared read

    for(i = 0; i < linkage.nrow(); ++i) {
      int cnt = 0;
      for(j = 0; j < linkage.ncol(); ++j)
        if(linkage(i, j) == -1)
          cnt++;
        if(cnt == linkage.ncol() - 1)
          continue;
        int found = 0;
        for(k = 0; k < out_hap.ncol(); ++k) {
          int flag = 1;
          for(j = 0; j < linkage.ncol(); ++j) {
            if(linkage(i, j) != -1 && linkage(i, j) != out_hap(k, j)) {
              flag = 0;
              break;
            }
          }
          if(flag) {
            found = 1;
            break;
          }
        }
        if(!found) {
          // check the original percent of such read if lower than 10% then exclude
          // Rcout << "missed one hap\n";
          IntegerVector tmp = linkage(i, _);
          // Rcout << tmp << "\n";
          for(k = new_end; k < order.size(); ++k) {
            int flag = 1;
            for(j = 0; j < linkage.ncol(); ++j) {
              if(linkage(i, j) != -1 && linkage(i, j) != hap_pool(hap_id[order[k]], j)) {
                flag = 0;
                break;
              }
            }
            if(flag) {
              IntegerVector tmp1 = hap_pool(hap_id[order[k]], _);
              // Rcout << tmp1 << "\n";
              if(find(new_top_id.begin(), new_top_id.end(), hap_id[order[k]]) == new_top_id.end()) {
                hap_abun.push_back(sorted[k]);
                new_top_id.push_back(hap_id[order[k]]);
                ab_id.push_back(k);
              }
              break;
            }
          }
        }
    }
    // Rcout << "abundance: ";
    // print_c_intvec(hap_abun);
    // check if should include more
    if(hap_abun.size() < order.size()) {
      IntegerVector tmp = wrap(hap_abun);
      int min_id = which_min(tmp);
      // Rcout << "min abundance id: " << min_id << "," << hap_abun[min_id] << "\n";
      if(ab_id[min_id] < order.size() - 1) { // check sorted abundance
        // Rcout << sorted[ab_id[min_id] + 1] << "," << sorted[ab_id[min_id]] << "\n";
        if(double(sorted[ab_id[min_id] + 1])/sorted[ab_id[min_id]] >= 0.6 ||
           (sorted[ab_id[min_id] + 1] < 10 && sorted[ab_id[min_id]] < 10)) {
          new_top_id.push_back(hap_id[order[ab_id[min_id] + 1]]);
          hap_abun.push_back(sorted[ab_id[min_id] + 1]);
        }
      }
    }
    out_hap = ss(hap_pool, wrap(new_top_id));

    // Rcout << "abundance: ";
    // print_c_intvec(hap_abun);
    // Rcout << "\n";
    // print_intmat(out_hap);
  }
  return out_hap;
}

List determine_snp_type(IntegerMatrix out_hap) {
  List nuc_unique(out_hap.ncol());
  List nuc_count(out_hap.ncol());
  for(int j = 0; j < out_hap.ncol(); ++j) {
    unordered_map<int, int> snp_type;
    for(int k = 0; k < out_hap.nrow(); ++k) {
      int tmp = out_hap(k, j);
      if(tmp != -1)
        snp_type[tmp]++;
    }
    IntegerVector counts;
    IntegerVector nucs;
    for(auto &snp: snp_type) {
      counts.push_back(snp.second);
      int nuc = snp.first;
      if(nuc == 4)
        nuc = -1;
      nucs.push_back(nuc);
    }
    nuc_unique[j] = nucs;
    nuc_count[j] = counts;
  }
  List ls = List::create(Named("nuc_unique") = nuc_unique,
                         Named("nuc_count") = nuc_count);
  return ls;
}

// get all the possible hap block based on the overlapped info, at least two reads covers two connected nuc
// [[Rcpp::export]]
List determine_hidden(List overlap_info, IntegerMatrix link_in, List opt, List hmm_info,
                      CharacterVector uni_alignment, int hap_min_pos)  {
  int sample = opt["n_sample"];
  IntegerVector start_t = overlap_info["start_t"];
  List loci = overlap_info["location"];
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector n_row = hmm_info["n_row"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  int t_max = hmm_info["t_max"];
  int given_prior = 0;
  int i, t, j, k;
  unordered_map<int, int> all_pos; // INDEX OF VARIANT SITES

  for(t = 0; t < undecided_pos.size(); ++t)
    all_pos[undecided_pos[t]] = t;
  start_t.push_back(t_max);

  IntegerVector pos_possibility;
  List hap_block(start_t.size() - 1); // start_t.size() may change due to only 1 read linked situation
  List block_sites(start_t.size() - 1);
  opt["sampling"] = 1;
  if(given_prior)
    pos_possibility = hmm_info["pos_possibility"];;
  // double lower_ab = opt["lower_ab"];
  // get sites in each hap block
  //avoid cases like A:1, 2; B: 3, C:2, in 2 blocks
  vector<std::set<int>> index_sets;
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
    index_sets.push_back(index_set1);
  }
  unordered_set<int> targets;
  for(t = 0; t < index_sets.size() - 1; ++t) {
    std::set<int> index_set1 = index_sets[t];
    std::set<int> index_set2 = index_sets[t + 1];
    for(auto it = index_set2.begin(); it != index_set2.end(); )
      if(index_set1.find(*it) != index_set1.end()) {// merge this one
        it = index_set2.erase(it);
      } else
        ++it;
      index_sets[t + 1] = index_set2;
  }
  for(t = 0; t < start_t.size() - 1; ++t) {
    std::set<int> index_set = index_sets[t];
    // need to loop though since the last e-1 does not necessary contain the last site in the block
    int fst_site = *index_set.begin();
    int lst_site = *index_set.rbegin();
    int start = all_pos[fst_site];
    int end = all_pos[lst_site];
    block_sites[t] = index_set;
    // Rcout << "t:" << t <<" "<< start  <<" "<< end <<"\n";
    if(start == end && !given_prior) {
      unsigned int ref_j = lst_site + hap_min_pos;
      // Rcout << "j: " << fst_site << "\n";
      IntegerVector sub_link = link_in(_, start);
      List tmp = unique_map(sub_link);
      IntegerVector counts = tmp["lengths"];
      IntegerVector nucs = tmp["values"];
      if(nucs[0] == -1) {
        counts.erase(counts.begin());
        nucs.erase(nucs.begin());
      }
      // Rcout << "counts" << counts << "|" << nucs << "\n"; // limit the repeated?
      opt["single"] = 1;
      List out = sbs_state(2, ref_j, nucs, counts, uni_alignment, opt);
      opt["single"] = 0;
      List hap_temp = out["haplotype"];
      n_row[fst_site] = out["n_row"];
      pos_possibility.push_back(out["n_row"]);
      hidden_states(fst_site) = hap_temp[0];
    } else {
      IntegerMatrix sub_link = link_in(_, Range(start, end));
      int count = 0;
      int num = sub_link.ncol();
      IntegerMatrix link_pre(sub_link.nrow(), num);
      for(i = 0; i < sub_link.nrow(); ++i) {
        IntegerVector read = sub_link(i, _);
        int rowsum = sum(read);
        if(rowsum == -num)
          continue;
        link_pre(count++, _) = sub_link(i, _); // contain replicates
      }
      IntegerMatrix link = link_pre(Range(0, count - 1), _);

      std::unordered_map<string, vector<int>> map = hash_mat2(link);
      int nres = map.size();
      vector<string> original_reads;
      IntegerVector idx(nres);
      List all_id(nres);
      int total = 0;
      for (auto itr = map.begin(); itr != map.end(); ++itr) {
        idx[total] = itr->second[0];
        all_id[total++] = wrap(itr->second);
        original_reads.push_back(itr->first);
      }
      IntegerMatrix uni_read = ss(link, idx);
      // Rcout << "uni_read:\n";
      // print_intmat(uni_read);

      IntegerMatrix hap_pool = sample_hap(original_reads, uni_read, all_id, sample);
      // Rcout << "sampling hap:\n";
      // print_intmat(hap_pool);
      hap_block(t) = hap_pool;

      if(!given_prior) {
        List hap_pool_info = determine_snp_type(link);
        List nucs_in = hap_pool_info["nuc_unique"];
        List count_in = hap_pool_info["nuc_count"];

        for(j = 0; j < uni_read.ncol(); ++j) {
          IntegerVector counts = count_in[j];
          IntegerVector nucs = nucs_in[j];
          unsigned int ref_j = undecided_pos[start + j] + hap_min_pos;
          // Rcout << "j: " << undecided_pos[start + j] << "\n";
          // Rcout << "counts" << counts << "|" << nucs << "\n";
          int num = nucs.size();
          // need to change this?
          List out = sbs_state(num, ref_j, nucs, counts, uni_alignment, opt);
          pos_possibility.push_back(out["n_row"]);
          n_row[undecided_pos[start + j]] = out["n_row"];
          List hap_temp = out["haplotype"];
          hidden_states(undecided_pos[start + j]) = hap_temp[0];
        }
      }
    }
  }
  List ls;
  unordered_map<int, int> update_state;
  IntegerVector num_states(t_max, 1);
  for (int j = 0; j < undecided_pos.size(); ++j)
    update_state[undecided_pos[j]] = pos_possibility[j];
  for(t = 0; t < t_max; ++t) {
    IntegerVector location_t = loci[t];
    if(location_t[0] != -1)
      for(i = 0; i < location_t.size(); ++i)
        num_states[t] *= update_state[location_t[i]];
    else
      num_states[t] = 1;
  }
  if(!given_prior) {
    ls = List::create(Named("block_sites") = block_sites,
                      Named("hap_block") = hap_block, // need to be adjusted
                      Named("n_row") = n_row,
                      Named("hidden_states") = hidden_states,
                      Named("pos_possibility") = pos_possibility,
                      Named("num_states") = num_states); // hidden_states at j with 1 unique nuc
  } else {
    ls = List::create(Named("block_sites") = block_sites,
                      Named("hap_block") = hap_block,
                      Named("num_states") = num_states);
  }
  return ls;
}

// make hap and transition info:
List limit_comb_samp(IntegerMatrix combination, List hidden_states, IntegerVector location,
                     List hap_block, List block_sites, int t, unsigned int num_states,
                     double lower_ab, IntegerMatrix linkage) {
  int i, j;
  // find the covered block sites in this hidden state
  IntegerVector block_loci = block_sites[t];
  IntegerMatrix block_hap = hap_block[t];
  IntegerVector index;
  int s = 0;
  for(i = 0; i < location.size(); ++i)
    for(j = s; j < block_loci.size(); ++j)
      if(location[i] == block_loci[j]) {
        index.push_back(j);
        s++;
        break;
      }

      IntegerMatrix overlapped_hap = ss(block_hap, index, 0); // then choose the top few
      IntegerMatrix kept_hap = determine_hap(overlapped_hap, lower_ab, linkage);
      // Rcout << "hap pool\n";
      // print_intmat(kept_hap);
      List ls = filter_combination(kept_hap, combination, hidden_states, location,
                                   kept_hap.ncol(), num_states);
      return(ls);
}

// [[Rcpp::export]]
List full_hap_samp2 (List hmm_info, List overlap_info, IntegerMatrix linkage,
                     unsigned int hap_length, int hap_min_pos, double lower_ab) {
  List hap_block = hmm_info["hap_block"];
  List block_sites = hmm_info["block_sites"];
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector n_row = hmm_info["n_row"];
  IntegerVector pos_possibility = hmm_info["pos_possibility"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int t, m, i;
  List loci = overlap_info["location"];
  List overlapped = overlap_info["overlapped"];
  IntegerVector overlapped_id = overlap_info["overlapped_id"];
  IntegerVector start_t = overlap_info["start_t"];

  List full_hap(t_max);
  List comb(t_max);
  IntegerMatrix hap = fill_all_hap(hidden_states, hap_length, n_row);
  IntegerVector new_num_states(t_max);
  // map for block_sites
  unordered_map<int, int> block_map;
  IntegerVector n_sites(block_sites.size());
  for(t = 0; t < block_sites.size(); ++t) {
    IntegerVector index = block_sites[t];
    n_sites[t] = index.size();
    for(m = 0; m < index.size(); ++m)
      block_map[index[m]] = t;
  }
  int n_start = start_t.size();
  List ori_comb(t_max);
  unordered_map<int, int> snp_mp;
  for(m = 0; m < undecided_pos.size(); ++m)
    snp_mp[undecided_pos[m]] = m;

  // find the states for the first state in each block (do this separately since same hidden states can be removed)
  for(t = 0; t < n_start; ++t) {
    IntegerVector loci_currt = loci[start_t[t]];
    int which_t = block_map[loci_currt[0]];
    // Rcout << start_t[t] << ", which t: " << which_t<< "\n";
    int no_sites = n_sites[which_t];
    int count = 0;
    List comb_info = find_combination(loci_currt, undecided_pos, pos_possibility);
    IntegerMatrix combination = comb_info["combination"];

    ori_comb[start_t[t]] = combination;
    IntegerMatrix final_comb;
    if(no_sites > 1) {
      int start = snp_mp[loci_currt[0]];
      int end = snp_mp[loci_currt[loci_currt.size() - 1]];
      IntegerMatrix sub_link = linkage(_, Range(start, end));
      int cnt = 0;
      int num = sub_link.ncol();
      IntegerMatrix link_pre(sub_link.nrow(), num);
      for(i = 0; i < sub_link.nrow(); ++i) {
        IntegerVector read = sub_link(i, _);
        int rowsum = sum(read);
        if(rowsum == -num)
          continue;
        link_pre(cnt++, _) = sub_link(i, _); // contain replicates
      }
      IntegerMatrix link = link_pre(Range(0, cnt - 1), _);
      List mat_info = hash_mat(link);
      IntegerVector idx = mat_info["idx"];
      IntegerMatrix uni_read = ss(link, idx);

      List exclude_info = limit_comb_samp(combination, hidden_states, loci_currt, hap_block,
                                          block_sites, which_t, num_states[start_t[t]], lower_ab, uni_read);
      IntegerVector exclude = exclude_info["exclude"];
      new_num_states[start_t[t]] = exclude_info["num_states"];
      IntegerMatrix new_comb(new_num_states[start_t[t]], combination.ncol());

      for(m = 0; m < num_states[start_t[t]]; ++m)
        if(!exclude[m])
          new_comb(count++, _) = combination(m, _);
        // dereplicate
        final_comb = dereplicate_states(new_comb, hidden_states, loci_currt,
                                        new_comb.ncol(), new_num_states[start_t[t]], 0);
    } else {
      final_comb = combination;
    }
    comb[start_t[t]] = final_comb;
    new_num_states[start_t[t]] = final_comb.nrow();
    // Rcout << "FINAL comb\n";
    // print_intmat(final_comb);
  }
  // get the unique, longest states which covers the whole hap block

  for(t = 0; t < t_max; ++t) {
    if(num_states[t] != 1 && overlapped_id[t] != -1) {
      int count = 0;
      IntegerVector loci_currt = loci[t];
      int which_t = block_map[loci_currt[0]];
      int last_t = overlapped_id[t];
      IntegerVector overlapped_t = overlapped[t];
      IntegerVector loci_lastt = loci[last_t];

      // skip the identical states
      if(loci_lastt[0] <= loci_currt[0] && loci_lastt[loci_lastt.size() - 1] >= loci_currt[loci_currt.size() - 1])
        continue;
      // first make the combination then loop to the previous states, add or remove some states
      // Rcout << t << "\n";
      // Rcout << "loci_currt: " << loci_currt << "\n";

      List comb_info = find_combination(loci_currt, undecided_pos, pos_possibility);
      IntegerMatrix combination = comb_info["combination"];
      ori_comb[t] = combination;

      int start = snp_mp[loci_currt[0]];
      int end = snp_mp[loci_currt[loci_currt.size() - 1]];
      IntegerMatrix sub_link = linkage(_, Range(start, end));
      int cnt = 0;
      int num = sub_link.ncol();
      IntegerMatrix link_pre(sub_link.nrow(), num);
      for(i = 0; i < sub_link.nrow(); ++i) {
        IntegerVector read = sub_link(i, _);
        int rowsum = sum(read);
        if(rowsum == -num)
          continue;
        link_pre(cnt++, _) = sub_link(i, _); // contain replicates
      }
      IntegerMatrix link = link_pre(Range(0, cnt - 1), _);
      List mat_info = hash_mat(link);
      IntegerVector idx = mat_info["idx"];
      IntegerMatrix uni_read = ss(link, idx);

      List exclude_info = limit_comb_samp(combination, hidden_states, loci_currt, hap_block,
                                          block_sites, which_t, num_states[t], lower_ab, uni_read);
      IntegerVector exclude = exclude_info["exclude"];
      new_num_states[t] = exclude_info["num_states"];
      IntegerMatrix new_comb(new_num_states[t], combination.ncol());

      for(m = 0; m < num_states[t]; ++m)
        if(!exclude[m])
          new_comb(count++, _) = combination(m, _);
        comb[t] = new_comb;
        // Rcout << new_num_states[t] << ",FINAL comb\n";
        // print_intmat(new_comb);
        // depth search and modify the former combination
        vector<int> index;
        // Rcout << "need to adjust\n";
        for(m = 0; m <= t; ++m)
          if(new_num_states[m] != 0) {
            IntegerVector loci_past = loci[m];
            if(which_t == block_map[loci_past[0]])
              index.push_back(m);
          }

          // print_c_intvec(index);
          // backward adding states

          for(m = index.size() - 1; m > 0; --m) {
            int first = index[m - 1];
            int second = index[m];
            IntegerVector loci_past = loci[first];
            int l_id = 0, r_id = 0;
            IntegerVector loci_cur = loci[second];
            IntegerVector mutual = intersect(loci_past, loci_cur);
            std::sort(mutual.begin(), mutual.end());
            // Rcout << "mutual loci: " << mutual << "\n";
            IntegerMatrix comb_cur = comb[second];
            IntegerMatrix comb_last = comb[first];
            //
            // Rcout<< first  << ",comb_last: \n";
            // print_intmat(comb_last);
            //
            for(int i = 0; i < loci_cur.size(); ++i)
              if(mutual[0] == loci_cur[i]) {
                r_id = i;
                break;
              }
              for(int i = 0; i < loci_past.size(); ++i)
                if(mutual[0] == loci_past[i]) {
                  l_id = i;
                  break;
                }

                IntegerMatrix comb_cur_ss = comb_cur(_, Range(r_id, r_id + mutual.size() - 1));
                IntegerMatrix comb_last_ss = comb_last(_, Range(l_id, l_id + mutual.size() - 1));
                unordered_set<string> cur_map;
                unordered_map<string, vector<int>> lst_map;
                //
                // Rcout << "subset comb_cur: ";
                // Rcout << r_id << " " << r_id + mutual.size() - 1 << "\n";
                //
                // print_intmat(comb_cur_ss);
                // Rcout << "subset comb_last: ";
                // Rcout << l_id << " " << l_id + mutual.size() - 1 << "\n";

                // print_intmat(comb_last_ss);

                for(int i = 0; i < comb_cur_ss.nrow(); ++i) {
                  string s;
                  for(int j = 0; j < mutual.size(); ++j)
                    s += to_string(comb_cur_ss(i, j));
                  if(cur_map.find(s) == cur_map.end())
                    cur_map.insert(s);
                }

                for(int i = 0; i < comb_last_ss.nrow(); ++i) {
                  string s;
                  for(int j = 0; j < mutual.size(); ++j)
                    s += to_string(comb_last_ss(i, j));
                  lst_map[s].push_back(i);
                }

                // check what is missing/redundant in last t
                unordered_set<string> missing;

                for(auto &s:cur_map)
                  if(lst_map.find(s) == lst_map.end())
                    missing.insert(s);

                  // print_c_intvec(redundant);
                  // Rcout << "missing: \t";
                  // for(auto &s:missing)
                  //   Rcout << s << "\t";
                  // Rcout << "\n";
                  // add missing
                  IntegerMatrix needed;
                  // if(!redundant.empty() && redundant.size() < comb_last.size()) // redundant meaning need to be added
                  //   comb_last = ss(comb_last, wrap(redundant));
                  if(!missing.empty()) {
                    // Rcout << "old comb:\n";
                    IntegerMatrix old_comb = ori_comb[first];
                    IntegerMatrix sub_comb = old_comb(_, Range(l_id, l_id + mutual.size() - 1));
                    vector<int> index;
                    unordered_set<string> record;
                    for(int i = 0; i < sub_comb.nrow(); ++i) {
                      string s;
                      for(int j = 0; j < mutual.size(); ++j)
                        s += to_string(sub_comb(i, j));
                      if(missing.find(s) != missing.end()) {
                        if(record.find(s) == record.end())
                          index.push_back(i);
                        record.insert(s);
                      }
                    }
                    // Rcout << "need to add:\n";
                    needed = ss(old_comb, wrap(index));
                    arma::Mat<int> m2 = as<arma::Mat<int>>(needed);
                    arma::Mat<int> m1 = as<arma::Mat<int>>(comb_last);
                    m1.insert_rows(0, m2);
                    comb[first] = wrap(m1);
                    new_num_states[first] = needed.nrow() + comb_last.nrow();
                    // Rcout << "adjusted comb: " << first << ", " << new_num_states[first] << "\n";
                    // print_intmat(comb[first]);
                  }

          }
          // froward pass
          for(m = 1; m < index.size(); ++m) {
            int first = index[m - 1];
            int second = index[m];
            IntegerVector loci_past = loci[first];
            // Rcout << "forward, " << first << ": " << second << "\n";
            int l_id = 0, r_id = 0;
            IntegerVector loci_cur = loci[second];
            IntegerVector mutual = intersect(loci_past, loci_cur);
            // Rcout << "loci_cur: " << loci_cur << "\n";
            std::sort(mutual.begin(), mutual.end());
            // Rcout << "mutual loci: " << mutual << "\n";
            IntegerMatrix comb_cur = comb[second];
            IntegerMatrix comb_last = comb[first];

            for(int i = 0; i < loci_cur.size(); ++i)
              if(mutual[0] == loci_cur[i]) {
                r_id = i;
                break;
              }
              for(int i = 0; i < loci_past.size(); ++i)
                if(mutual[0] == loci_past[i]) {
                  l_id = i;
                  break;
                }

                IntegerMatrix comb_cur_ss = comb_cur(_, Range(r_id, r_id + mutual.size() - 1));
                IntegerMatrix comb_last_ss = comb_last(_, Range(l_id, l_id + mutual.size() - 1));
                unordered_set<string> lst_map;
                unordered_set<string> cur_map;
                //
                //         Rcout << "subset comb_cur: ";
                //         Rcout << r_id << " " << r_id + mutual.size() - 1 << "\n";
                //
                //         print_intmat(comb_cur_ss);
                //         Rcout << "subset comb_last: ";
                //         Rcout << l_id << " " << l_id + mutual.size() - 1 << "\n";

                // print_intmat(comb_last_ss);

                for(int i = 0; i < comb_last_ss.nrow(); ++i) {
                  string s;
                  for(int j = 0; j < mutual.size(); ++j)
                    s += to_string(comb_last_ss(i, j));
                  lst_map.insert(s);
                }

                for(int i = 0; i < comb_cur_ss.nrow(); ++i) {
                  string s;
                  for(int j = 0; j < mutual.size(); ++j)
                    s += to_string(comb_cur_ss(i, j));
                  cur_map.insert(s);
                }

                // check what is missing in t
                unordered_set<string> missing;
                for(auto &s:lst_map)
                  if(cur_map.find(s) == cur_map.end())
                    missing.insert(s);
                  // Rcout << "missing: \t";
                  // for(auto &s:missing)
                  //   Rcout << s << "\t";
                  // Rcout << "\n";
                  IntegerMatrix needed;
                  if(!missing.empty()) {
                    IntegerMatrix old_comb = ori_comb[index[m]];
                    IntegerMatrix sub_comb = old_comb(_, Range(r_id, r_id + mutual.size() - 1));
                    vector<int> index;
                    unordered_set<string> record;
                    for(int i = 0; i < sub_comb.nrow(); ++i) {
                      string s;
                      for(int j = 0; j < mutual.size(); ++j)
                        s += to_string(sub_comb(i, j));
                      if(missing.find(s) != missing.end()) {
                        if(record.find(s) == record.end())
                          index.push_back(i);
                        record.insert(s);
                      }
                    }
                    needed = ss(old_comb, wrap(index));
                    arma::Mat<int> m2 = as<arma::Mat<int>>(needed);
                    arma::Mat<int> m1 = as<arma::Mat<int>>(comb_cur);
                    m1.insert_rows(0, m2);
                    comb[second] = wrap(m1);
                    new_num_states[second] = needed.nrow() + comb_cur.nrow();
                  }
                  // Rcout << "adjusted comb: " << new_num_states[second] << "\n";
                  // print_intmat(comb[second]);
          }
    }
  }
  // Rcout << "end for comb\n";
  for(t = 0; t < t_max; ++t) {
    if(num_states[t] != 1 && overlapped_id[t] != -1) {
      IntegerVector loci_currt = loci[t];
      int last_t = overlapped_id[t];
      IntegerVector overlapped_t = overlapped[t];
      IntegerVector loci_lastt = loci[last_t];
      //  Rcout << t <<"\n";
      // Rcout << loci_currt << "\n";
      // Rcout << loci_lastt << "\n";
      if(loci_lastt[0] <= loci_currt[0] && loci_lastt[loci_lastt.size() - 1] >= loci_currt[loci_currt.size() - 1]) {
        if(loci_lastt.size() > loci_currt.size()) {
          IntegerMatrix comb_in = comb[last_t];
          // Rcout << "contained in\n";
          // print_intmat(comb_in);
          IntegerMatrix new_comb = unique_overlap(overlapped_t, comb_in, loci_lastt, new_num_states[last_t]);
          new_num_states[t] = new_comb.nrow();
          comb[t] = new_comb;
        }
        else if(loci_lastt.size() == loci_currt.size()) {
          // Rcout << "same sites as before\n";
          comb[t] = comb[last_t];
          new_num_states[t] = new_num_states[last_t];
        }
      }
    } else if(num_states[t] == 1) {
      new_num_states[t] = 1;
      comb[t] = -1;
    }
  }
  // Rcout << "end for all comb\n";
  List full_hap_t;
  for(t = 0; t < t_max; ++t) {
    int count = 0;
    // Rcout << t << "\n";
    // should be num_states since those are non-variant sites
    if(num_states[t] == 1) {
      full_hap_t = List(1);
      full_hap_t[0] = hap(_, Range(time_pos[t] - hap_min_pos, time_pos[t] + p_tmax[t] - hap_min_pos - 1));
    } else {
      // Rcout << "adjusted comb: \n";
      IntegerMatrix new_comb = comb[t];
      IntegerVector loci_currt = loci[t];
      full_hap_t = List(new_num_states[t]);

      // print_intmat(new_comb);
      for(m = 0; m < new_num_states[t]; ++m) {
        IntegerMatrix haplotype = make_hap(hidden_states, hap, loci_currt, p_tmax[t], new_comb(m, _),
                                           time_pos[t], loci_currt.size(), hap_min_pos);
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

