#include <RcppArmadillo.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "utils.h"
#include "hmm_state.h"
using namespace Rcpp;
using namespace std;
#define NUM_CLASS 4

// [[Rcpp::depends(RcppArmadillo)]]

IntegerMatrix doubleheter(IntegerVector a) {
  IntegerVector temp = {a[0], a[0], a[1], a[1], a[1], a[1], a[0], a[0],
                        a[1], a[0], a[0], a[1], a[0], a[1], a[1], a[0]};
  temp.attr("dim") = Dimension(4, 4);
  return(as<IntegerMatrix>(temp));
}

IntegerMatrix Twopossible(IntegerVector a) {
  IntegerVector temp = {a[0], a[1], a[0], a[1], a[1], a[0], a[1], a[0]};
  temp.attr("dim") = Dimension(2, 4);
  return(as<IntegerMatrix>(temp));
}

IntegerMatrix Fourpossible(IntegerVector small, int big) {
  IntegerVector temp = {small[0], small[1], big, big, small[1], small[0], big, big,
                        big, big, small[0], small[1], big, big, small[1], small[0]};
  temp.attr("dim") = Dimension(4, 4);
  return(as<IntegerMatrix>(temp));
}

IntegerMatrix call_permute_N(vector<int> a, unsigned int genome_A) {
  Permutation per;
  vector<vector<int> > res =  per.permuteUnique(a);
  IntegerMatrix permutation(res.size(), res[0].size());
  for (int i = 0; i < res.size(); i++)
    for (int j = 0; j < res[i].size(); j++)
      permutation(i, j) = res[i][j];
  IntegerVector del(res.size() * 2, -1);
  del.attr("dim") = Dimension(res.size(), 2);
  arma::mat m1 = as<arma::mat>(permutation);
  arma::mat m2 = as<arma::mat>(del);
  arma::mat out;
  if(genome_A)
    out = join_rows(m2, m1);
  else
    out = join_rows(m1, m2);

  return(wrap(out));
}

List aux_noN_S2(IntegerVector sum_site, IntegerVector hap_site,
                List opt) {
  int n_row;
  double sum;
  int max_id = which_max(sum_site);

  int left = !max_id;
  sum = sum_site[1] + sum_site[0];
  List ls(2);
  double majority = opt["majority"];
  int emprical = opt["emprical"];
  int sampling = opt["sampling"];
  int single = opt["single"];
  int dbheter = opt["dbheter"];
  // if one nuc appears 90% or one nuc appears only once
  if(!sampling && !single) {
    if (((sum - sum_site[max_id]) == 1 && sum > 4) || (sum_site[max_id]/sum > majority)) {
      n_row = 1;
      IntegerVector temp = {hap_site[max_id], hap_site[max_id], hap_site[max_id], hap_site[max_id]};
      ls["temp"] = temp;
    } else {
      IntegerMatrix temp;
      IntegerMatrix temp2;
      double three_hap = opt["three_hap"];
      if (sum_site[max_id]/sum >= three_hap) {//largest one appears 3 times as much as the rest
        n_row = 6;
        if(emprical) {
          temp = call_permute({hap_site[max_id], hap_site[max_id], hap_site[max_id], hap_site[left]});
          temp2 = Twopossible(hap_site);
          arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
          arma::Mat<int> m2 = as<arma::Mat<int>>(temp2);
          arma::Mat<int> m3 = join_cols(m1, m2);
          temp = wrap(m3);
        }
      } else {
        if(emprical)
          temp = Twopossible(hap_site);
        n_row = 2;
        if(dbheter) {
          temp2 = doubleheter(hap_site);
          arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
          arma::Mat<int> m2 = as<arma::Mat<int>>(temp2);
          arma::Mat<int> m3 = join_cols(m1, m2);
          temp = wrap(m3);
          n_row = 6;
        }
      }
      ls["temp"] = temp;
    }
  } else {
    IntegerMatrix temp;
    IntegerMatrix temp2;
    double three_hap = opt["three_hap"];
    if (sum_site[max_id]/sum >= three_hap) {//largest one appears 3 times as much as the rest
      n_row = 6;
      temp = call_permute({hap_site[max_id], hap_site[max_id], hap_site[max_id], hap_site[left]});
      temp2 = Twopossible(hap_site);
      arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
      arma::Mat<int> m2 = as<arma::Mat<int>>(temp2);
      arma::Mat<int> m3 = join_cols(m1, m2);
      temp = wrap(m3);
    } else {
      temp = Twopossible(hap_site);
      n_row = 2;
    }
    if(single) {
      if(n_row == 6) {
        IntegerMatrix temp3 = ss(temp, {0, 2, 4}); // reduce the replicated hidden states
        n_row = 2;
        ls["temp"] = temp3;
      } else {
        ls["temp"] = temp;
      }
    } else if(sampling) {
      ls["temp"] = temp;
    }
  }

  ls["n_row"] = n_row;
  return ls;
}

List aux_noN_S3(IntegerVector sum_site, IntegerVector hap_site, List opt) {
  int n_row;
  double sum;
  IntegerMatrix temp;
  List from_2(2);
  List ls(2);

  int min_id = which_min(sum_site);
  IntegerVector sum_site2(2);
  IntegerVector hap_site2(2);
  int num = 0;
  for(unsigned int i = 0; i < 3; ++i) {
    if(i != min_id) {
      sum_site2[num] = sum_site[i];
      hap_site2[num] = hap_site[i];
      num++;
    }
  }
  int third_nuc = opt["third_nuc"];
  int emprical = opt["emprical"];
  if(third_nuc == 0) {
    Rcout << "exclude 3rd nuc\n";
    ls = aux_noN_S2(sum_site2, hap_site2, opt);
  } else {
    sum = sum_site[2] + sum_site[1] + sum_site[0];
    // double p0 = sum_site[0]/sum;
    // double p1 = sum_site[1]/sum;
    // double p2 = sum_site[2]/sum;
    IntegerVector small(2);
    IntegerVector small_count(2);
    int max_id = which_max(sum_site);
    double two_hap = opt["two_hap"];
    double left_range = opt["left_range"];
    double right_range = 2 - left_range;
    // different rules to reduce the possibilities, basically rule out things like ATAG
    if (((sum - sum_site[max_id])/sum_site[max_id] <= left_range
           && (sum - sum_site[max_id])/sum_site[max_id] >= right_range) ||
             (sum_site[max_id]/sum >= two_hap)) {
      if(emprical) {
        IntegerVector wset = {0, 1, 2};
        IntegerVector bset = {max_id};
        IntegerVector diff = setdiff(wset, bset);
        small = {hap_site[diff[0]], hap_site[diff[1]]};
        temp = Fourpossible(small, hap_site[max_id]);
      }
      n_row = 4;
    } else { //hopefully this won't happen, when three roughly equal happens
      if(emprical) {
        small = {hap_site[1], hap_site[2]};
        temp = Fourpossible(small, hap_site[0]);
        arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
        small = {hap_site[1], hap_site[0]};
        temp = Fourpossible(small, hap_site[2]);
        arma::Mat<int> m2 = as<arma::Mat<int>>(temp);
        small = {hap_site[2], hap_site[0]};
        temp = Fourpossible(small, hap_site[1]);
        arma::Mat<int> m3 = as<arma::Mat<int>>(temp);
        arma::Mat<int> out = join_cols(m1, m2, m3);
        temp = wrap(out);
      }
      n_row = 12;
    }
    // also include the situation of 2 possibles
    List list_2pos = aux_noN_S2(sum_site2, hap_site2, opt);
    int more_row = list_2pos["n_row"];
    n_row += more_row;
    if(emprical) {
      IntegerMatrix temp2 = list_2pos["temp"];
      arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
      arma::Mat<int> m2 = as<arma::Mat<int>>(temp2);
      arma::Mat<int> out = join_cols(m1, m2);
      temp = wrap(out);
    }
    ls["n_row"] = n_row;
    ls["temp"] = temp;
  }

  return ls;
}

//3 different nuc with gaps in the alignment
List N3_gap(IntegerVector hap_site, IntegerVector sum_site, unsigned int ref_j,
            CharacterVector uni_alignment, List opt) {
  double r = double(sum_site[1])/(sum_site[2] + sum_site[1]);
  int n_row;
  List ls(2);
  int emprical = opt["emprical"];
  if(r >= 0.45 && r <= 0.55) {
    if(emprical) {
      IntegerVector temp(4 * NUM_CLASS);
      if (uni_alignment[ref_j] != "M")
        temp = {-1, -1, hap_site[1], hap_site[2], -1, -1, //J in universal alignment
                hap_site[2], hap_site[1], hap_site[1], hap_site[2], -1,  -1,
                hap_site[2],hap_site[1], -1, -1};
      temp.attr("dim") = Dimension(4, 4); // by column
      IntegerMatrix final = as<IntegerMatrix>(temp);
      ls["temp"] = final;
    }
    n_row = 4;
  } else {
    IntegerVector tp = {sum_site[1], sum_site[2]};
    int id = which_max(tp) + 1;
    if(uni_alignment[ref_j] != "M") {
      if(emprical) {
        IntegerVector temp(2 * NUM_CLASS);
        temp = {-1, hap_site[id], -1, hap_site[id], hap_site[id], -1, hap_site[id], -1};
        temp.attr("dim") = Dimension(2, 4);
        IntegerMatrix final = as<IntegerMatrix>(temp);
        ls["temp"] = final;
      }
      n_row = 2;
    } else {
      IntegerVector temp(NUM_CLASS);
      temp = {hap_site[id], hap_site[id], hap_site[id], hap_site[id]};
      ls["temp"] = temp;
      n_row = 1;
    }
  }
  ls["n_row"] = n_row;
  return ls;
}
/*
 * determine the number of hidden states site by site
 */
List sbs_state(unsigned int num, unsigned int ref_j, IntegerVector hap_site, IntegerVector sum_site,
               CharacterVector uni_alignment, List opt) {
  double sum;
  unsigned int n_row;
  List haplotype(1);
  //record possible hidden states[only suitable for ployploids]
  if (num == 2) {
    sum = sum_site[0] + sum_site[1];
    if(hap_site[0] == -1) {
      IntegerVector temp;
      if(uni_alignment[ref_j] != "M" && sum_site[0]/sum >= 0.4) {
        temp = IntegerVector(2 * NUM_CLASS);
        temp = {-1, hap_site[1], -1, hap_site[1],hap_site[1], -1, hap_site[1], -1}; //gaps
        temp.attr("dim") = Dimension(2, 4);
        n_row = 2;
      } else {
        temp = IntegerVector(NUM_CLASS);
        temp = {hap_site[1], hap_site[1], hap_site[1], hap_site[1]};
        n_row = 1;
      }
      haplotype[0] = temp;
    } else { // if N not appears
      List out = aux_noN_S2(sum_site, hap_site, opt);
      haplotype[0] = out["temp"];
      n_row = out["n_row"];
    }
  }
  else if (num == 3) {
    List out;
    sum = sum_site[2] + sum_site[1] + sum_site[0];
    if (hap_site[0] == -1) {
      // if N appears
      if(sum_site[0]/sum >= 0.4) {
        out = N3_gap(hap_site,sum_site, ref_j, uni_alignment, opt);
      } else {
        IntegerVector sum_site_cleaned = {sum_site[1], sum_site[2]};
        IntegerVector hap_site_cleaned = {hap_site[1], hap_site[2]};
        out = aux_noN_S2(sum_site_cleaned, hap_site_cleaned, opt);
      }
    } else {
      out = aux_noN_S3(sum_site, hap_site, opt);
    }
    haplotype[0] = out["temp"];
    n_row = out["n_row"];
  }
  else if (num == 4) {
    List out;
    sum = sum_site[3] + sum_site[2] + sum_site[1] + sum_site[0];
    if(hap_site[0] == -1) {
      if(sum_site[0]/sum >= 0.4) {
        int min_id = which_min(sum_site);
        IntegerVector sum_site2(3);
        IntegerVector hap_site2(3);
        int num = 0;
        for (unsigned int i = 0; i < 4; ++i)
          if (i != min_id) {
            sum_site2[num] = sum_site[i];
            hap_site2[num] = hap_site[i];
            num++;
          }
          out = N3_gap(hap_site2, sum_site2, ref_j, uni_alignment, opt);
      } else {
        IntegerVector sum_site_cleaned = {sum_site[1], sum_site[2], sum_site[3]};
        IntegerVector hap_site_cleaned = {hap_site[1], hap_site[2], hap_site[3]};
        out = aux_noN_S3(sum_site_cleaned, hap_site_cleaned, opt);
      }
    } else {
      int min_id = which_min(sum_site);
      IntegerVector sum_site2(3);
      IntegerVector hap_site2(3);
      int num = 0;
      for (unsigned int i = 0; i < 4; ++i)
        if (i != min_id) {
          sum_site2[num] = sum_site[i];
          hap_site2[num] = hap_site[i];
          num++;
        }
        out = aux_noN_S3(sum_site2, hap_site2, opt);
    }
    haplotype[0] = out["temp"];
    n_row = out["n_row"];
  }
  else if(num == 5) { // include -1 0 1 2 3
    List out;
    int min_id = which_min(sum_site);
    IntegerVector sum_site2(4);
    IntegerVector hap_site2(4);
    int num = 0;
    for (unsigned int i = 0; i < 5; ++i)
      if (i != min_id) {
        sum_site2[num] = sum_site[i];
        hap_site2[num] = hap_site[i];
        num++;
      }
      min_id = which_min(sum_site2);
      hap_site2.erase(min_id);
      sum_site2.erase(min_id);
      sum = sum_site2[2] + sum_site2[1] + sum_site2[0];
      if(hap_site2[0] == -1 && sum_site2[0]/sum >= 0.4) {
        out = N3_gap(hap_site2, sum_site2, ref_j, uni_alignment, opt);
      } else {
        out = aux_noN_S3(sum_site2, hap_site2, opt);
      }
      haplotype[0] = out["temp"];
      n_row = out["n_row"];
  }
  List ls = List::create(
    Named("n_row") = n_row,
    Named("haplotype") = haplotype);

  return(ls);
}

// remake linkage (based on observed data)
IntegerMatrix remake_linkage(IntegerMatrix sub_link, unsigned int num) {
  unsigned int i, j, k, i1;
  arma::mat sub_uni = unique_rows(as<arma::mat>(sub_link));
  IntegerMatrix link_uni = wrap(sub_uni);
  List new_link(link_uni.nrow());
  unsigned int total_row = 0;
  for (i = 0; i < link_uni.nrow(); ++i) {
    List nuc_info = unique_map(link_uni(i, _));
    IntegerVector nuc = nuc_info["values"];
    IntegerVector nuc_count = nuc_info["lengths"];
    if (nuc.size() == 1 && nuc[0] == -1) // skip the read does not cover any site
      continue;
    // if(nuc_count[0] == num - 1 || nuc_count[0] == 1)
    //   continue;
    if(nuc[0] != -1) { // read covers all site
      total_row++;
      new_link[i] = link_uni(i, _);
      continue;
    }
    IntegerVector idx(num);
    for(j = 0; j < num; ++j)
      if(link_uni(i, j) == -1)
        idx(j) = 1; // indicate -1 is here

    int move_out = 0;
    // if this read is contained in others
    for (i1 = 0; i1 < link_uni.nrow(); ++i1) {
      if (i1 == i)
        continue;
      int count = 0;
      // List nuc_info = unique_map(link_uni(i1, _));
      // IntegerVector nuc1 = nuc_info["values"];
      // IntegerVector nuc_count1 = nuc_info["lengths"];
      // if(nuc_count1[0] >= nuc_count[0])
      //   continue;
      for(j = 0; j < num; ++j) {
        if(!idx(j))
          if(link_uni(i, j) == link_uni(i1, j))
            count++;
      }
      if(count == num - nuc_count[0]) {
        // Rcout << i << "move" << "\n";
        move_out = 1;
        break;
      }
    }
    if(move_out)
      continue;
    List missing(num);
    IntegerVector flag(num);
    int in_row_num = 0;
    for (j = 0; j < num; ++j) {
      if (link_uni(i, j) == -1) {
        List nuc_col = unique_map(link_uni(_, j));
        IntegerVector nuc_unique = nuc_col["values"];
        // Rcout << nuc_unique << "\n";
        if(nuc_unique[0] == -1)
          missing[j] = nuc_unique[Range(1, nuc_unique.size() - 1)];
        else
          missing[j] = nuc_unique;
        in_row_num++;
      } else
        flag[j] = 1;
    }

    int add_row = 0;
    if(in_row_num != 1) {
      IntegerMatrix missing_rows = comb_element(missing, flag, in_row_num);
      add_row = missing_rows.nrow();
      IntegerMatrix new_link_i(add_row, num);
      int count = 0;
      for (j = 0; j < num; ++j) { // make fake reads with missing linkage info
        if(flag[j] != 1)
          new_link_i(_, j) = missing_rows(_, count++);
        else
          for (k = 0; k < add_row; ++k)
            new_link_i(k, j) = link_uni(i, j); // repeat the non-missing ones
      }
      new_link[i] = new_link_i;
    } else {
      IntegerMatrix new_link_i;
      IntegerVector tmp;
      for (j = 0; j < num; ++j)
        if(flag[j] != 1) {
          tmp = missing[j];
          add_row = tmp.size();
          new_link_i = IntegerMatrix(add_row, num);
        }
        for (j = 0; j < num; ++j) { // make fake reads with missing linkage info
          if(flag[j] != 1) {
            new_link_i(_, j) = tmp;
          } else {
            for (k = 0; k < add_row; ++k)
              new_link_i(k, j) = link_uni(i, j);
          }
        }
        new_link[i] = new_link_i;
    }
    total_row += add_row;
  }

  IntegerMatrix new_link_out(total_row, num);
  total_row = 0;
  for(i = 0; i < link_uni.nrow(); ++i) {
    if(new_link[i] == R_NilValue)
      continue;
    IntegerVector tmp = new_link[i];
    tmp.attr("dim") = Dimension(tmp.size()/num, num);
    IntegerMatrix new_link_i = as<IntegerMatrix>(tmp);
    for (k = 0; k < tmp.size()/num; ++k)
      new_link_out(total_row++, _) = new_link_i(k, _);
  }
  // finally, remove duplicated rows
  arma::mat new_linkage = unique_rows(as<arma::mat>(new_link_out));
  IntegerMatrix out = wrap(new_linkage);
  return(out);
}

IntegerVector best_branch(IntegerMatrix link, List transition, NumericVector initial,
                          List possi_nuc, int i) {
  unsigned int j, k, l, m , w;
  List comb_in(link.ncol());
  List llk_in(link.ncol());
  int id = 0;
  // possible or determined nuc at each position
  for(j = 0; j < link.ncol(); ++j) {
    IntegerVector nuc = possi_nuc[j];
    if(link(i, j) != -1) {
      for(l = 0 ; l < nuc.size(); ++l)
        if(link(i, j) == nuc[l]) {
          id = l;
          break;
        }
        comb_in(j) = link(i, j);
    } else {
      comb_in(j) = nuc;
    }
  }
  IntegerVector flag(link.ncol());
  // IntegerMatrix poss_reads = comb_element(comb_in, flag, link.ncol());
  // get state likelihood
  for(j = 0; j < link.ncol(); ++j) {
    // Rcout << j << "\t";
    IntegerVector nuc = possi_nuc[j];
    if(link(i, j) != -1) {// if nuc at j is known
      for(l = 0 ; l < nuc.size(); ++l)
        if(link(i, j) == nuc[l]) {
          id = l;
          // Rcout << "nuc " << nuc[l] << "\t";
          break;
        }
        if(j == 0) {
          llk_in(j) = initial[id];
          // Rcout << "ini " << initial[id] << "\t";
        } else {
          NumericMatrix trans = transition[j - 1];
          IntegerVector nuc2 = possi_nuc[j - 1];
          int id1 = 0;
          if(link(i, j - 1) != -1) {// if nuc at j-1 is known, trans prob is known
            for(l = 0 ; l < nuc2.size(); ++l)
              if(link(i, j - 1) == nuc2[l]) {
                // Rcout << "nuc(j-1) " << nuc2[l] << "\t";
                id1 = l;
                break;
              }
              llk_in(j) = trans(id1, id);
              // Rcout << "trans " << trans(id1, id) << "\n";
          } else{
            // Rcout << "trans: all\n";
            llk_in(j) = trans(_, id);
          }
        }
    } else {
      if(j == 0)
        llk_in[j] = initial;
      else {
        NumericMatrix trans = transition[j - 1];
        IntegerVector nuc2 = possi_nuc[j - 1];
        if(link(i, j - 1) != -1) {
          for(l = 0 ; l < nuc2.size(); ++l)
            if(link(i, j) == nuc2[l]) {
              // Rcout << "nuc(j-1) " << nuc2[l] << "\n";
              id = l;
              break;
            }
            llk_in(j) = trans(id, _);
        } else {
          llk_in(j) = trans;
        }
      }
    }
  }
  IntegerVector hidden_state(llk_in.size());
  List path(llk_in.size());
  List backptr(llk_in.size() - 1);
  int b_next = 0;
  for(k = 0; k < llk_in.size(); ++k) {
    // Rcout << k << "\n";
    IntegerVector nuc = comb_in(k);
    NumericVector path_t(nuc.size());
    IntegerVector backptr_t(nuc.size());
    if(k == 0) {
      NumericVector trans = llk_in(k);
      for(l = 0; l < nuc.size(); ++l) {
        path_t(l) = trans(l);
      }
      // Rcout << path_t << "\n";
    } else {
      NumericVector tran = llk_in(k);
      int len = tran.size();
      int nrow = len/nuc.size();
      tran.attr("dim") = Dimension(nrow, nuc.size());
      NumericMatrix trans = as<NumericMatrix>(tran);
      NumericVector path_last = path[k - 1];
      for(m = 0; m < trans.ncol(); ++m) {
        double max = -INFINITY;
        int max_id = 0;
        for(w = 0; w < trans.nrow(); ++w) {
          double max_prob = path_last(w) + trans(w, m);
          if (max_prob > max) {
            max = max_prob;
            max_id = w;
          }
        }
        path_t(m) = max;
        backptr_t[m] = max_id;
      }
      backptr(k - 1) = backptr_t;
      // Rcout << path_t << "\n";
      // Rcout << backptr_t << "\n";
    }
    path(k) = path_t;
  }
  // Rcout << "decode\n";
  k = llk_in.size() - 1;
  double max = -INFINITY;
  IntegerVector nuc = comb_in(k);
  NumericVector path_t = path(k);
  for(m = 0; m < nuc.size(); ++m) {
    if (path_t(m) > max) {
      b_next = m;
      max = path_t(m);
    }
  }
  hidden_state(k) = nuc[b_next];

  while (k--) {
    // Rcout << k << "\n";
    IntegerVector nuc = comb_in(k);
    IntegerVector backptr_t = backptr(k);
    // Rcout << backptr_t << "\n";
    b_next = backptr_t[b_next];
    hidden_state(k) = nuc[b_next];
  }

  return(hidden_state);
}


// use markov chain to make the linkage
IntegerMatrix mc_linkage(IntegerMatrix sub_link, int num) {
  unsigned int i, j, k, l;
  //remove non-covered reads
  NumericVector initial;
  IntegerMatrix uni;
  IntegerMatrix link_pre(sub_link.nrow(), sub_link.ncol());
  int count = 0;
  for(i = 0; i < sub_link.nrow(); ++i) {
    IntegerVector read = sub_link(i, _);
    int rowsum = sum(read);
    if(rowsum == -num)
      continue;
    link_pre(count++, _) = sub_link(i, _);
  }
  IntegerMatrix link = link_pre(Range(0, count - 1), _);
  // Rcout << "cleaned link\n";
  // print_intmat(link);
  // if there is only one read link them, then include all the possible combinations
  List compression = hash_mat(link);
  IntegerVector id = compression["idx"];
  IntegerMatrix comp_link = ss(link, id);
  int complete_seq = 0;
  for(k = 0; k < comp_link.nrow(); ++k) {
    int flg = 0;
    for(j = 0 ; j < comp_link.ncol(); ++j)
      if(comp_link(k, j) == -1) {
        flg = 1;
        break;
      }
      if(!flg)
        complete_seq++;
  }
  if(complete_seq == 1) {
    // Rcout << "only 1 linked read\n";
    uni = remake_linkage(link, num);
  } else {
    List transition(link.ncol() - 1);
    List possi_nuc(link.ncol());
    for(j = 0 ; j < link.ncol() - 1; ++j) {
      List nuc_info = unique_map(link(_, j));
      IntegerVector nuc = nuc_info["values"];
      IntegerVector nuc_count = nuc_info["lengths"];
      // Rcout << j << ":\t" << nuc << "\t" << nuc_count << "\n";
      int state1 = nuc.size();
      int start = 0;
      possi_nuc[j] = nuc;
      if(nuc[0] == -1) {
        state1 = nuc.size() - 1;
        start = 1;
        possi_nuc[j] = nuc[Range(1, nuc.size() - 1)];
      }

      if(j == 0) {
        double total = count;
        initial = IntegerVector(state1);
        if(nuc[0] == -1) {
          for(k = 1; k < nuc.size(); ++k)
            initial[k - 1] = log(nuc_count[k]/(total - nuc_count[0])); // in the ascending order
        }
        else {
          for(k = 0; k < nuc.size(); ++k)
            initial[k] = log(nuc_count[k]/total);
        }
      }
      // unique rows and the count
      List nuc1_info = unique_map(link(_, j + 1));
      IntegerVector nuc1 = nuc1_info["values"];
      int state2 = nuc1.size();
      possi_nuc[j + 1] = nuc1;
      if(nuc1[0] == -1) {
        possi_nuc[j + 1] = nuc1[Range(1, nuc1.size() - 1)];
        state2 = nuc1.size() - 1;
      }

      IntegerMatrix link_unique(link.nrow(), 2);
      int unique_ct = 0;
      for(k = 0; k < link.nrow(); ++k)
        if(link(k, j + 1) != -1 && link(k, j) != -1) {
          link_unique(unique_ct, 0) = link(k, j);
          link_unique(unique_ct++, 1) = link(k, j + 1);
        }
      // for(k = 0; k < link.nrow(); ++k)
      //   for(l = 0; l < 2; ++l)
      //     Rcout << link_unique(k, j) << "\t";
      IntegerMatrix in_link = link_unique(Range(0, unique_ct - 1), _);
      // Rcout << "in_link:\n";
      // print_intmat(in_link);
      List unique_row = hash_mat(in_link);
      List all_id = unique_row["all_id"];
      IntegerVector idx = unique_row["idx"];
      NumericMatrix trans(state1, state2);
      for(k = 0; k < state1; ++k)
        for(i = 0; i < state2; ++i)
          trans(k, i) = R_NegInf;
      // Rcout << state1 << "\t" << state2 << "\n";
      List nuc_info_uni = unique_map(in_link(_, 0));
      IntegerVector nuc_count_uni = nuc_info_uni["lengths"];
      IntegerVector nuc_uni = nuc_info_uni["values"];

      for(k = start; k < nuc.size(); ++k) {
        // Rcout<< "\n" << nuc[k] << "\n";
        int nuc_id = 0;
        for(l = 0; l < nuc_uni.size(); ++l)
          if(nuc_uni[l] == nuc[k]) {
            nuc_id = l;
            break;
          }
        for(i = 0; i < idx.size(); ++i) {
          double de = nuc_count_uni[nuc_id];
          IntegerVector sub_read = in_link(idx[i], _);
          // Rcout << "unique sub " << sub_read << "\n";
          // IntegerVector sub_read = ordered_read(i, _);
          // if(sub_read[1] == -1) {
          //   de--;
          //   continue;
          // }
          if(nuc[k] != sub_read[0])
            continue;
          IntegerVector sub_id = all_id[i];
          // Rcout << "uniques " << sub_id << "\n";
          double nu = sub_id.size();
          int col_id = 0;
          int flag = 0;
          if(state2 < nuc1.size())
            flag = 1;
          for(l = 0; l < nuc1.size(); ++l)
            if(nuc1[l] != -1 && nuc1[l] == sub_read[1])
              col_id = l - flag;
            // Rcout <<"col_id " << col_id << "\n" ;
            // if(nuc[k] == sub_read[0]) {
            // Rcout <<"de " << de << " nu " << nu<< "\n" ;
          trans(k - start, col_id) = log(nu/de); // log likelihood
            // }
        }
      }
      transition[j] = trans;
    }
    //
    // // now use MC to impute the missing nuc
    arma::mat uniqu_link = unique_rows(as<arma::mat>(link));
    IntegerMatrix sub_uni_link = wrap(uniqu_link);
    int len = sub_uni_link.ncol();
    int row_num = sub_uni_link.nrow();
    IntegerMatrix mc_reads(row_num, len);
    // NumericVector llk(row_num);
    count = 0;
    for(i = 0; i < row_num; ++i) {
      int flag = 0; // record the no. of -1 in read
      for(j = 0; j < len; ++j)
        if(sub_uni_link(i, j) == -1) {
          flag = 1;
          break;
        }

      if(flag == 0) {
        mc_reads(i, _) = sub_uni_link(i, _);
        // llk[i] = 0;
        continue;
      }
      IntegerVector tmp = sub_uni_link(i, _);
      // Rcout << i << ", " << tmp << "\n";
      mc_reads(i, _) = best_branch(sub_uni_link, transition, initial, possi_nuc, i);
      // llk[i] = inferred_read["llk"];
      // mc_reads(i, _) = inferred_read["read"];
    }
    arma::mat uniqu = unique_rows(as<arma::mat>(mc_reads));
    uni = wrap(uniqu);
  }
  return(uni);
}

List filter_combination(IntegerMatrix sub_link, IntegerMatrix combination,
                        List hidden_states, IntegerVector location,
                        unsigned int num, unsigned int num_states) {
  unsigned int i, j, k, idx, m;
  IntegerMatrix sub_hap(NUM_CLASS, num);
  IntegerVector exclude(num_states);
  int count, linkage_len, all_excluded;
  linkage_len = num - 1;
  unsigned int n_observation = sub_link.nrow();
  all_excluded = 0;
  for (m = 0; m < num_states; ++m) {
    exclude(m) = 0;
    IntegerVector comb = combination(m, _);
    // Rcout << "comb: " << comb << "\n";
    count = 0;
    for (k = 0; k < NUM_CLASS; ++k) {
      // Rcout << "k" << k << ": ";
      for (j = 0; j < num; ++j) {
        IntegerMatrix hidden = hidden_states[location[j]];
        idx = comb[j];
        sub_hap(k, j) = hidden(idx, k);
        // Rcout << sub_hap(k, j) << "\t";
      }
      // Rcout << "\n read " << "\n";
      for (i = 0; i < n_observation; i++) {
        int flag = 0;
        for (j = 0; j < num - 1; ++j) {
          // Rcout << sub_link(i, j)  << sub_link(i, j + 1) <<  "\t";
          // if sub_link(i, j) == 4, meanning gap in genome, for comparison, convert the coding to -1
          int a = sub_link(i, j), b = sub_link(i, j + 1);
          if(a == 4)
            a = -1;
          if(b == 4)
            b = -1;
          if (sub_hap(k, j) == a && sub_hap(k, j + 1) == b)
            flag++;
        }
        if (flag >= linkage_len) {
          // Rcout << "kept\n" ;
          count++;
          break;
        }
      }
    }
    // Rcout << "count "<< count << "\n";
    if (count != NUM_CLASS) {
      exclude(m) = 1;
      all_excluded++;
    }
  }
  if(num_states == all_excluded) {
    Rcout << "not enough linkage information, so include all\n";
    all_excluded = 0;
    exclude = IntegerVector(num_states);
  }

  List out = List::create(
    Named("num_states") = num_states - all_excluded,
    Named("exclude") = exclude);
  return(out);
}
// limit combination based on linkage information. For the states not start a new sequence(overlap,
// keep the last overlapped combination and then make the rest)
List limit_comb_t0(IntegerMatrix combination, List hidden_states, IntegerVector location,
                   IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx,
                   unsigned int num_states, unsigned int use_MC) {
  IntegerMatrix old_sub_link = linkage_info(_, Range(start_idx, start_idx + num - 1));
  // int cut_off;
  // all_excluded = num_states;
  //remake the linkage
  IntegerMatrix sub_link;
  if(!use_MC)
    sub_link = remake_linkage(old_sub_link, num);
  else
    sub_link = mc_linkage(old_sub_link, num);
  // Rcout << "linkage:\n";
  // print_intmat(sub_link);
  List ls = filter_combination(sub_link, combination, hidden_states, location,
                               num, num_states);
  return(ls);
}

// only applicable for
//  0    4    0    0    2 = 0    4    0    0    3
// if max = 6, 1=2, 3=4, 5=6; else 0=1, 2=3, 4=5
IntegerMatrix dereplicate_states(IntegerMatrix new_combination, List hidden_states, IntegerVector location,
                                 unsigned int num, unsigned int num_states, unsigned int db_heter) {
  unsigned int i, j, k;
  IntegerMatrix sub_hap(NUM_CLASS, num);
  IntegerMatrix long_hap(num_states, NUM_CLASS * num);
  for(i = 0; i < num_states; ++i) {
    IntegerVector comb = new_combination(i, _);
    for (k = 0; k < NUM_CLASS; ++k) {
      for (j = 0; j < num; ++j) {
        IntegerMatrix hidden = hidden_states[location[j]];
        int id = comb[j];
        sub_hap(k, j) = hidden(id, k);
      }
    }
    // pick out the db heter hidden states and then
    IntegerMatrix ordered_hap = sort_mat(sub_hap, NUM_CLASS, num, db_heter);
    // Rcout << "sorted" << "\n";
    // print_intmat(ordered_hap);
    IntegerVector tmp = matrix2vec(ordered_hap);
    long_hap(i, _) = tmp;
  }
  // hash long_hap
  List info = hash_mat(long_hap);
  IntegerVector idx = info["idx"];
  IntegerMatrix new_comb = ss(new_combination, idx);
  // print_intmat(new_comb);
  return(new_comb);
}

// trans_indicator: indicate which state can transfer to which, further_limit indicate some states should not be considered
List prepare_ini_hmm (unsigned int t_max, IntegerVector num_states, List trans_indicator, List further_limit) {
  List trans_new_ind(t_max - 1);
  unsigned int w, m, t;

  for(t = 0; t < t_max - 1; ++t) {
    // Rcout << t << "\n";
    IntegerVector more_limit_last = further_limit[t];
    IntegerVector more_limit = further_limit[t + 1];
    IntegerMatrix trans_ind = trans_indicator[t];
    IntegerMatrix trans_ind_new(trans_ind.nrow(), num_states[t + 1]);
    int count2 = 0;
    int count = 0;
    for (w = 0; w < trans_ind.ncol(); ++w)
      if(!more_limit[w])
        trans_ind_new(_, count++) = trans_ind(_, w);
      IntegerMatrix trans_ind_new2(num_states[t], num_states[t + 1]);
      for (m = 0; m < trans_ind.nrow(); ++m)
        if(!more_limit_last[m])
          trans_ind_new2(count2++, _) = trans_ind_new(m, _);
        trans_new_ind[t] = trans_ind_new2;
  }

  return(trans_new_ind);
}

// get the unique rows for overlapped region
IntegerMatrix unique_overlap(IntegerVector overlapped, IntegerMatrix overlap_comb, IntegerVector overlap_loci,
                             unsigned int overlap_new_states) {
  unsigned int i, m;
  // unsigned int n_observation = linkage_info.nrow();
  //decide the first t
  unsigned int overlap_len = overlapped.size();
  // limit the space based on the last limition first
  int start_overlap = 0;
  for(i = 0; i < overlap_loci.size(); ++i)
    if (overlap_loci[i] == overlapped[0]) {
      start_overlap = i;
      break;
    }

    //find the unique states
    IntegerMatrix comb_last(overlap_new_states, overlap_len);
    unsigned int count = 0;
    for (m = 0; m < overlap_new_states; ++m) {
      // if(!exclude_last[m]) {
      IntegerVector tmp = overlap_comb(m, _);
      comb_last(count++, _) = tmp[Range(start_overlap, start_overlap + overlap_len - 1)];
    }
    arma::mat out = unique_rows(as<arma::mat>(comb_last));
    return(wrap(out));
}

//fill haplotype at the rest sites (with variation) at time t
IntegerMatrix make_hap(List hidden_states, IntegerMatrix haplotype, IntegerVector location, unsigned int p_tmax,
                       IntegerVector combination, unsigned int time_pos, unsigned int num, int hap_min_pos) {
  unsigned int j, k, idx;

  for(j = 0; j < num; ++j) {
    IntegerMatrix hidden = hidden_states[location[j]];
    // Rcout << "hidden " << location[j] << "\n";
    // print_intmat(hidden);
    idx = combination[j];
    // Rcout << idx << "\n";
    for(k = 0; k < NUM_CLASS; ++k)
      haplotype(k, location[j]) = hidden(idx, k);
  }
  // Rcout << time_pos - hap_min_pos << " " << time_pos + p_tmax - hap_min_pos - 1 << "\n";
  return(haplotype(_, Range(time_pos - hap_min_pos, time_pos + p_tmax - hap_min_pos - 1)));
}

/*
 * Return the combination of sites with variation within one time t[NOTE:time_pos might not start from 0, but undecided_pos starts from 0]
 */
List find_combination(IntegerVector location, IntegerVector undecided_pos, IntegerVector pos_possibility) {
  //possible combination of the rest non-unique loci
  IntegerVector location_len(location.size());
  unsigned int num = 0;
  unordered_map<int, int> mp;
  for(unsigned int i = 0; i < undecided_pos.size(); ++i)
    mp[undecided_pos[i]] = i;

  for(unsigned int i = 0; i < location.size(); ++i) {
    int id = mp[location[i]];
    location_len(num++) = pos_possibility[id];
  }

  IntegerMatrix combination = call_cart_product(location_len);
  List ls = List::create(
    Named("combination") = combination, // possible comb
    Named("num") = num);
  return(ls);
}

IntegerMatrix new_combination(List hmm_info, IntegerVector location, IntegerVector overlapped, IntegerMatrix overlap_comb,
                              IntegerVector overlap_loci, IntegerMatrix linkage_info, unsigned int overlap_new_states,
                              unsigned int use_MC) {

  IntegerMatrix first_comb = unique_overlap(overlapped, overlap_comb, overlap_loci, overlap_new_states);
  //find combination of the rest position(include 1 overlap to make sure the linkage)
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector pos_possibility = hmm_info["pos_possibility"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  unsigned int i, j, m, w, k;
  int overlap_len = overlapped.size();

  IntegerVector left_loci = location[Range(overlap_len - 1, location.size() - 1)];
  int count = 0;
  int len = left_loci.size();
  IntegerVector left_possible(len);
  for(i = 0; i < undecided_pos.size(); ++i) {
    if (undecided_pos[i] >= left_loci[0] && undecided_pos[i] <= left_loci[len - 1]) {
      left_possible[count++] = pos_possibility[i];
    }
  }
  // get the combination of the non-overlapped location (include last overlapped)
  IntegerMatrix combination = call_cart_product(left_possible[Range(0, count - 1)]);

  // get the appeared possibilities at the overlapped position
  IntegerVector last_col = first_comb(_, first_comb.ncol() - 1);
  List first_uni = unique_map(last_col); // start might not from 0 (e.g. ailgnment starts from 2)
  IntegerVector n1 = first_uni["lengths"];
  IntegerVector allowed = first_uni["values"];
  // IntegerVector allowed = unique(last_col);
  IntegerVector first_col = combination(_, 0);
  IntegerVector exist = unique(first_col);
  List exclude_info(2);
  int flag = 0;
  int num = 0;
  IntegerMatrix new_combination(combination.nrow(), combination.ncol());
  unsigned int start_idx = 0;
  for(i = 0; i < undecided_pos.size(); ++i)
    if (undecided_pos[i] == left_loci[0]) {
      start_idx = i;
      break;
    }

  if(!setequal(exist, allowed))
    flag = 1;

  if(flag) {  // if the exist contains more possibility than allowed, only keeps the allowed
    for(m = 0; m < combination.nrow(); ++m)
      for(w = 0; w < allowed.size(); ++w)
        if(allowed(w) == combination(m, 0)) {
          new_combination(num++, _) = combination(m, _);
          break;
        }
        exclude_info = limit_comb_t0(new_combination(Range(0, num - 1), _), hidden_states,
                                     left_loci, linkage_info, combination.ncol(), start_idx, num, use_MC);
  } else {
    exclude_info = limit_comb_t0(combination, hidden_states, left_loci, linkage_info,
                                 combination.ncol(), start_idx, combination.nrow(), use_MC);
  }

  // Now give the limited combination
  int num_states = exclude_info["num_states"];
  IntegerVector exclude = exclude_info["exclude"];
  // Rcout << "exclude "<< exclude << "\n";
  IntegerMatrix next_comb(num_states, combination.ncol());
  count = 0;
  // Now combine first and second part[make sure the connection states appear]
  if(flag) {
    for(m = 0; m < num; ++m)
      if(!exclude[m])
        next_comb(count++, _) = new_combination(m, _);
  } else {
    for(m = 0; m < combination.nrow(); ++m)
      if(!exclude[m])
        next_comb(count++, _) = combination(m, _);
  }
  // this will introduce more states not shown in the reads linkage, but to keep the trans works, have to...
  IntegerVector new_exist = next_comb(_, 0);
  // Rcout << "new_exist: " << new_exist << "\n";
  if(!setequal(new_exist, allowed)) {
    IntegerVector diff = setdiff(allowed, new_exist);
    Rcout << "add more possibilities " << diff << "\n";;
    IntegerMatrix extra(diff.size(), combination.ncol());
    arma::Mat<int> m1 = as<arma::Mat<int>>(next_comb);
    next_comb = IntegerMatrix(num_states + diff.size(), combination.ncol());
    count = 0;
    if(flag) {
      for(w = 0; w < diff.size(); ++w)
        for(m = 0; m < new_combination.nrow(); ++m) {
          IntegerVector tmp = new_combination(m, _);
          if(new_combination(m, 0) == diff[w]) {
            extra(count++, _) = new_combination(m, _);
            break;
          }
        }
    } else {
      for(w = 0; w < diff.size(); ++w)
        for(m = 0; m < combination.nrow(); ++m) {
          IntegerVector tmp = combination(m, _);
          if(combination(m, 0) == diff[w]) {
            extra(count++, _) = combination(m, _);
            break;
          }
        }
    }
    // print_intmat(extra);
    arma::Mat<int> m2 = as<arma::Mat<int>>(extra);
    m1.insert_rows(1, m2);
    next_comb = wrap(m1);
  }

  IntegerVector new_col = next_comb(_, 0);
  List second_uni = unique_map(new_col); // start might not from 0 (e.g. ailgnment starts from 2)
  IntegerVector n2 = second_uni["lengths"];
  // IntegerVector possible = second_uni["values"];
  int all = 0;
  for(m = 0; m < n2.size(); ++m)
    all += n2[m] * n1[m];
  IntegerMatrix final_comb(all, location.size());

  all = 0;
  for(k = 0; k < last_col.size(); ++k)
    for(w = 0; w < new_col.size(); ++w) {
      if(new_col(w) == last_col(k)) {
        for(j = 0; j < overlap_len; ++j)
          final_comb(all, j) = first_comb(k, j);
        for(i = 1; i < len; ++i)
          final_comb(all, i + overlap_len - 1) = next_comb(w, i);
        all++;
      }
    }

  return(final_comb);
}

/*
 * Fill the positions without variation
 */
IntegerMatrix fill_all_hap(List hidden_states, unsigned int hap_length, IntegerVector n_row) {
  unsigned int j, k;
  IntegerMatrix haplotype(NUM_CLASS, hap_length);
  IntegerVector nuc_j(NUM_CLASS);
  for(j = 0; j < hap_length; ++j)
    if(n_row(j) == 1) { // fill the loci with only 1 possibility, can be optimized by using the part has been filled and add/trim
      nuc_j = hidden_states[j];
      for(k = 0; k < NUM_CLASS; ++k)
        haplotype(k, j) = nuc_j[k];
    }
    return haplotype;
}
