#include <RcppArmadillo.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "hmm_state.h"
#include "utils.h"

using namespace Rcpp;
using namespace std;

#define MLOGIT_CLASS 4
#define NUM_CLASS 4
#define QUA_MAX 106
#define QUA_ERR 10

// [[Rcpp::depends(RcppArmadillo)]]

//https://stackoverflow.com/questions/51110244/in-rcpp-how-to-get-a-user-defined-structure-from-c-into-r
// namespace Rcpp {
//   template <>
//   SEXP wrap(const par& x);
// }
//
//
//
// namespace Rcpp {
//   template <>
//   SEXP wrap(const par& x) {
//     Rcpp::NumericVector beta;
//     return Rcpp::wrap(Rcpp::Named("beta") = Rcpp::wrap(x.beta));
//   };
// }

typedef unsigned char xy_t;
/**
 * Convert char to xy.
 *
 * @param c	ASCII char
 * @return	xy_t
 */
inline xy_t char_to_xy(char c)
{
  if (c == 'A' || c == 'a')
    return 'A' >> 1 & 3L;
  else if (c == 'C' || c == 'c')
    return 'C' >> 1 & 3L;
  else if (c == 'G' || c == 'g')
    return 'G' >> 1 & 3L;
  else if (c == 'T' || c == 't' || c == 'U' || c == 'u')
    return 'T' >> 1 & 3L;
  else
    return 1 << 7;	/* non-nuc */
} /* char_to_xy */

CharacterVector r_xy_to_char = {"A", "C", "T", "G"};


/**
 * Types of nucleotide encodings.
 */
enum {
  IUPAC_ENCODING,		/*!< iupac_t */
XY_ENCODING,		/*!< xy_t */
};

int uni_sum(List unique_map, unsigned int cut_off)
{
  NumericVector result_vals = unique_map["lengths"];
  NumericVector result_keys = unique_map["values"];
  int n = result_keys.size();
  int cumsum = 0;
  int key = 0;
  for (unsigned int i = 0; i < n; ++i) {
    cumsum += result_vals(i);
    if(cumsum >= cut_off) {
      key = result_keys(i);
      break;
    }
  }
  return key;
}

int top_n_map(List unique_map)
{
  NumericVector result_vals = unique_map["lengths"];
  NumericVector result_keys = unique_map["values"];
  int n = result_keys.size();
  int key = 0;
  for (unsigned int i = n; i --> 0;) {
    if(result_vals(i) >= NUM_CLASS) {
      key = result_keys(i);
      break;
    }
  }
  return key;
}

// [[Rcpp::export]]
List read_data(std::string path, unsigned int old_v) {
  int j, k, l;
  unsigned int i, m, count = 0, count_del = 0;
  char c;
  char str[100];
  int n_observation = 0;

  FILE* fp = fopen(path.c_str(), "r");
  if(!fp) {
    Rcpp::stop("File opening failed");
  }
  int temp_n = -1;
  while (fgets(str, sizeof(str), fp)) {
    std::sscanf(str, "%d %d %d %d %c", &i, &j, &k, &l, &c);
    /* exclude -1 in ref_pos and qua */
    if (l == -1)
      count_del++;
    else if(k != -1) {
      if (temp_n != i) {
        temp_n = i;
        n_observation++;
      }
      count++;
    }
  }

  unsigned int del_num = 0;
  IntegerVector del_obs_index;
  IntegerVector del_ref_pos;
  IntegerVector del_read_pos;
  IntegerVector del_id;
  IntegerVector del_flag;

  if(count_del) {
    del_obs_index = IntegerVector(count_del);
    del_ref_pos = IntegerVector(count_del);
    del_read_pos = IntegerVector(count_del);
    del_id = IntegerVector(n_observation);
    del_flag = IntegerVector(n_observation);
  }

  IntegerVector qua(count);
  IntegerVector obs(count);
  IntegerVector obs_index(count);
  IntegerVector ref_pos(count);
  IntegerVector read_pos(count);
  int total;

  rewind(fp);

  count = 0;
  count_del = 0;
  // count_ins = 0;
  while (fgets(str, sizeof(str), fp)) {
    std::sscanf(str, "%d %d %d %d %c", &i, &j, &k, &l, &c);
    /* exclude -1 in ref_pos and qua */
    if (l == -1) {
      del_obs_index[count_del] = i;
      del_read_pos[count_del] = j;
      del_ref_pos[count_del] = k;
      del_flag[i-1] = 1;
      count_del++;
    }
    else if(k != -1) {
      obs_index[count] = i;
      read_pos[count] = j;
      ref_pos[count] = k;
      qua[count] = l;
      obs[count] = c;
      count++;
    }
  }
  fclose(fp);
  fp = NULL;

  total = count;
  IntegerVector length(n_observation);

  // deletion
  IntegerVector del_count;
  IntegerVector del_length_all(n_observation);
  IntegerVector del_strat_id(n_observation);
  if (count_del) {
    List tmp = unique_map(del_obs_index);
    del_id = tmp["values"];
    del_num = del_id.size();
    del_count = tmp["lengths"];
    for (m = 0; m < del_num; ++m)
      del_length_all[del_id[m] - 1] = del_count[m];

    for (i = 1; i < n_observation; ++i) {
      if(del_flag[i] == 1)
        for (j = 0; j < i; ++j)
          del_strat_id[i] += del_length_all[j];
      else
        del_strat_id[i] = -1;
    }
  }

  // non indel
  for (m = 0; m < count; ++m) {
    i = obs_index[m];
    ++length[i-1];
    obs[m] = char_to_xy(obs[m]);
  }

  /* index of read in non-indel loci */
  IntegerVector index(n_observation);
  for (m = 1; m < n_observation; ++m)
    index[m] = index[m - 1] + length[m - 1];

  // true length (include insertion) and fake length (include deletion) (some reads' alignment does not start from 0, other one also count that in)
  IntegerVector true_length(n_observation);
  IntegerVector fake_length(n_observation);
  for (i = 0; i < n_observation; ++i) {
    true_length[i] = length[i];
    fake_length[i] = length[i] + del_length_all[i];
    //fake_length[i] = length[i] + del_length_all[i] + ref_pos[index[i]];
  }

  int max_len = 0; // largest position of reference
  int min_len = 0; // smallest position of reference
  min_len = min(ref_pos);
  int over_hapmax = 0; // indicate if the length of reads is more than the hap_max
  if(old_v == 1) { // use the old version for genotyping targeted region not WGS [the starting aligned position has to be 0]
    /* find the longest reference position && appears more than no. of classes (4) */
    List uni_map = unique_map(ref_pos);
    max_len = top_n_map(uni_map);
    for(i = 0; i < total; ++i)
      if (max_len < ref_pos[i]) {
        over_hapmax = 1;
        break;
      }
  } else { // for WGS
    max_len = max(ref_pos);
  }
  max_len = max_len + 1;

  IntegerMatrix ref_index(n_observation, max_len - min_len); // start from the first aligned pos
  for(i = 0; i < n_observation; ++i)
    for(m = 0; m < max_len - min_len; ++m)
      ref_index(i, m) = -1;

  /* Find the index of ref position under different read of every j, which ref pos aligned to which read pos */
  for(i = 0; i < n_observation; ++i)
    for(m = 0; m < max_len - min_len; ++m)
      for(j = 0; j < length[i]; ++j)
        if(ref_pos[index[i] + j] == m + min_len) {
          ref_index(i, m) = j;
          break;
        }

        IntegerVector non_covered_site(max_len - min_len);
        unsigned int num;
        for (m = 0; m < max_len - min_len; ++m) {
          num = 0;
          for(i = 0; i < n_observation; ++i)
            if (ref_index(i, m) == -1)
              num++;
            if (num == n_observation)
              non_covered_site[m] = 1;
        }
        List del(10);
        if(count_del) {
          del = List::create(
            Named("del_id") = del_id[Range(0, del_num - 1)],
                                    Named("del_id_all") = del_obs_index,
                                    Named("del_flag") = del_flag,
                                    Named("del_read_pos") = del_read_pos,
                                    Named("del_ref_pos") = del_ref_pos,
                                    Named("del_length") = del_count, //no. of deletion in each read
                                    Named("del_num") = del_num,  //no. of reads have deletion
                                    Named("del_total") = count_del,
                                    Named("del_length_all") = del_length_all,
                                    Named("del_strat_id") = del_strat_id);}

        List ls = List::create(
          Named("id") = obs_index,
          Named("read_pos") = read_pos,
          Named("ref_pos") = ref_pos,
          Named("nuc") = obs,
          Named("qua") = qua,
          Named("n_observation") = n_observation,
          Named("length") = length,
          Named("true_length") = true_length,
          Named("fake_length") = fake_length,
          Named("total") = total,
          Named("start_id") = index,
          Named("ref_start") = min_len,
          Named("ref_length_max") = max_len,
          Named("ref_idx") = ref_index,
          Named("over_hapmax") = over_hapmax,
          Named("non_covered_site") = non_covered_site,
          Named("deletion") = del);

        return ls;
}

IntegerMatrix format_data_simple(List dat_info, IntegerMatrix haplotype, int time_pos = -1) {
  unsigned int i, k;
  int total = dat_info["total"];
  int hap_max_pos = dat_info["ref_length_max"];
  // int hap_length = hap_max_pos - hap_min_pos;
  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector obs_index = dat_info["id"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];

  unsigned int len = total * NUM_CLASS;
  IntegerMatrix df_new(len, 6);
  /* pick out the haplotypes according to ref_pos */
  for (i = 0; i < total; ++i) {
    for (k = 0; k < NUM_CLASS; ++k)
      if(time_pos != -1) {
        //Rcout << haplotype(k, ref_pos[i] - tp) << "\t"; //Both ref_pos[i] and time_pos has the starting aligned location
        df_new(i * NUM_CLASS + k, 5) = haplotype(k, ref_pos[i] - time_pos);
      } else {
        int hap_min_pos = dat_info["ref_start"];
        // the reference position might be longer than the sampled haplotypes
        if(ref_pos[i] > hap_max_pos - 1)
          df_new(i * NUM_CLASS + k, 5) = -1;
        else
          df_new(i * NUM_CLASS + k, 5) = haplotype(k, ref_pos[i] - hap_min_pos);
      }
  }

  /* repeat the data for mnlogit */
  for (k = 0; k < total; ++k)
    for (i = 0; i < NUM_CLASS; ++i) {
      df_new(k * NUM_CLASS + i, 1) = ref_pos[k];
      df_new(k * NUM_CLASS + i, 2) = read_pos[k];
      df_new(k * NUM_CLASS + i, 4) = qua[k];
      df_new(k * NUM_CLASS + i, 0) = obs_index[k];
      df_new(k * NUM_CLASS + i, 3) = obs[k];
    }
    return(df_new);
}

//prepare data for mnlogit, mnlogit only takes data without indels in read or in haplotypes
//time_pos is used when calling hmm method
// [[Rcpp::export]]
DataFrame format_data(List dat_info, IntegerMatrix haplotype, int time_pos = -1) {
  unsigned int i, k, l;
  int input_arr[] = {0, 1, 2, 3};
  size_t input_arr_sz = sizeof input_arr / sizeof *input_arr;

  int total = dat_info["total"];
  int hap_max_pos = dat_info["ref_length_max"];
  // int hap_length = hap_max_pos - hap_min_pos;
  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector obs_index = dat_info["id"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];

  unsigned int len = total * MLOGIT_CLASS * NUM_CLASS;
  IntegerVector r_ref_pos(len);
  IntegerVector r_read_pos(len);
  IntegerVector r_qua(len);
  IntegerVector r_obs(len);
  IntegerVector r_hap_nuc(len);
  IntegerVector mode(len);
  IntegerVector id(len);

  /* pick out the haplotypes according to ref_pos */
  for (i = 0; i < total; ++i) {
    for (k = 0; k < NUM_CLASS; ++k)
      for (l = 0; l < MLOGIT_CLASS; ++l) {
        if(time_pos != -1) {
          //Rcout << haplotype(k, ref_pos[i] - tp) << "\t"; //Both ref_pos[i] and time_pos has the starting aligned location
          r_hap_nuc[i * MLOGIT_CLASS * NUM_CLASS + MLOGIT_CLASS * k + l] = haplotype(k, ref_pos[i] - time_pos);
        } else {
          int hap_min_pos = dat_info["ref_start"];
          // the reference position might be longer than the sampled haplotypes
          if(ref_pos[i] > hap_max_pos - 1)
            r_hap_nuc[i * MLOGIT_CLASS * NUM_CLASS + MLOGIT_CLASS * k + l] = -1;
          else
            r_hap_nuc[i * MLOGIT_CLASS * NUM_CLASS + MLOGIT_CLASS * k + l] = haplotype(k, ref_pos[i] - hap_min_pos);
        }
      }
  }

  /* repeat the data for mnlogit */
  for (k = 0; k < total; ++k)
    for (i = 0; i < MLOGIT_CLASS * NUM_CLASS; ++i) {
      r_ref_pos[k * MLOGIT_CLASS * NUM_CLASS + i] = ref_pos[k];
      r_read_pos[k * MLOGIT_CLASS * NUM_CLASS + i] = read_pos[k];
      r_qua[k * MLOGIT_CLASS * NUM_CLASS + i] = qua[k];
      id[k * MLOGIT_CLASS * NUM_CLASS + i] = obs_index[k];
    }

    for (i = 0; i < total * MLOGIT_CLASS * NUM_CLASS; ++i)
      r_obs[i] = input_arr[i % input_arr_sz];

  for (i = 0; i < total; ++i)
    for (k = 0; k < MLOGIT_CLASS * NUM_CLASS; ++k)
      if (r_obs[k + MLOGIT_CLASS * NUM_CLASS * i] == obs[i])
        mode[k + MLOGIT_CLASS * NUM_CLASS * i] = 1;

  DataFrame df_new = DataFrame::create(
    Named("id") = id,
    Named("mode") = mode,
    Named("read_pos") = r_read_pos,
    Named("ref_pos") = r_ref_pos,
    Named("qua") = r_qua,
    Named("nuc") = r_obs,
    Named("hap_nuc") = r_hap_nuc);

  return(df_new);
}

// [[Rcpp::export]]
CharacterVector to_char(IntegerVector nuc) {
  unsigned int i;
  CharacterVector nuc_char(nuc.size());
  for(i = 0; i < nuc.size(); ++i) {
    nuc_char[i] = r_xy_to_char[nuc[i]];
  }
  return(nuc_char);
}

// [[Rcpp::export]]
List hmm_info(List dat_info, CharacterVector uni_alignment,
              List opt, unsigned int sbs = 1) {
  unsigned int i, j, t;
  unsigned int total = dat_info["total"];
  int n_observation = dat_info["n_observation"];
  int hap_max_pos = dat_info["ref_length_max"];
  int hap_min_pos = dat_info["ref_start"];
  int hap_length = hap_max_pos - hap_min_pos;
  IntegerMatrix ref_index = dat_info["ref_idx"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector index = dat_info["start_id"];
  IntegerVector fake_length = dat_info["fake_length"];
  IntegerVector length = dat_info["length"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector qua = dat_info["qua"];
  IntegerVector non_covered = dat_info["non_covered_site"];
  double cut_off = opt["cut_off"];

  /* Find the number of reads with alignment start from each t (hash) */
  IntegerVector read_start(n_observation);
  for(i = 0; i < n_observation; ++i)
    read_start[i] = ref_pos[index[i]];

  List start_t = unique_map(read_start); // start might not from 0 (e.g. ailgnment starts from 2)
  IntegerVector n_t = start_t["lengths"];
  IntegerVector time_pos = start_t["values"];
  IntegerVector pos;
  NumericVector prop(hap_length);
  // get the coverage of each site (including -)
  if(cut_off != 0) {
    List deletion = dat_info["deletion"];
    IntegerVector del_ref_pos;
    unsigned int del_total = 0;
    if(deletion[0] != R_NilValue) {
      del_ref_pos = deletion["del_ref_pos"];
      del_total = deletion["del_total"];
    }
    unsigned int enuma = total + del_total;
    pos = IntegerVector(enuma);
    for(i = 0; i < total; ++i)
      pos(i) = ref_pos(i);
    for(i = 0; i < del_total; ++i)
      pos(i + total) = del_ref_pos(i);
    List coverage_stat = unique_map(pos);
    IntegerVector coverage = coverage_stat["lengths"];
    for(i = 0; i < hap_length; ++i)
      prop[i] = coverage[i] * cut_off;
  }

  /* Find the max p in each t and record the observations in each set t*/
  unsigned int t_max = time_pos.length();
  List n_in_t(t_max);
  IntegerVector p_tmax(t_max);
  unsigned int max, count_nt;
  for(t = 0; t < t_max; ++t) {
    count_nt = 0;
    IntegerVector nt(n_observation);
    max = 0;
    for(i = 0; i < n_observation; ++i)
      if(ref_pos[index[i]] == time_pos[t]) {
        nt[count_nt++] = i;
        if(max < fake_length[i])
          max = fake_length[i];
      }
      p_tmax[t] = max;
      n_in_t(t) = nt[Range(0, count_nt - 1)];
  }

  unsigned int count, num;
  IntegerVector nuc_j(n_observation);
  IntegerVector qua_ij(n_observation);
  List nuc_unique(hap_length);
  List nuc_count(hap_length);
  List nuc(2);
  NumericVector keys(MLOGIT_CLASS + 1);
  NumericVector vals(MLOGIT_CLASS + 1);
  // List haplotype(hap_length);
  // IntegerVector n_row(hap_length);
  // IntegerVector pos_possibility(hap_length);
  // IntegerVector undecided_pos(hap_length);
  IntegerVector record_cov(hap_length);
  int cov_count = 0;
  for (j = 0; j < hap_length; ++j) {
    if(!non_covered[j]) {
      unsigned int ref_j = j + hap_min_pos;
      int qua_min = QUA_MAX;
      count = 0;
      int non_del = 0;

      for (i = 0; i < n_observation; ++i)
        if((ref_j >= ref_pos[index[i]]) && (ref_j <= ref_pos[index[i] + length[i] - 1])) {
          if(ref_index(i, j) == -1) { // deletion appears
            nuc_j(count++) = -1;
          } else {
            qua_ij(count) = qua[index[i] + ref_index(i, j)];
            if(qua_ij[count] < qua_min)
              qua_min = qua_ij[count];
            nuc_j(count++) = obs[index[i] + ref_index(i, j)];
            non_del++;
          }
        }

      if(non_del == 1) {
        record_cov[cov_count++] = ref_j;
        //record the min and the max index that covered by one nuc, exclude them when compare with other algorithms
        Rcout << "site " << ref_j;
        Rcout << " is covered by 1 read \n";
        // n_row[j] = 1;
        // IntegerVector tmp = {nuc_j[0], nuc_j[0], nuc_j[0], nuc_j[0]};
        // haplotype(j) = tmp;
        nuc_unique[j] = nuc_j[0];
        nuc_count[j] = 1;
        continue;
      }

      // Rcout << j << ":\t" << non_del<< ":" << count<<"\t" << qua_min<<  "\n";
      IntegerVector tmp(count);
      IntegerVector tmp_nuc;
      // if there are only two reads covered and only keep the one with the highest prob
      if(non_del == 2) {
        record_cov[cov_count++] = ref_j;
        int max_qua = 0;
        int keep_id = 0;
        if(nuc_j[0] == nuc_j[1])
          keep_id = 0;
        else {
          for(i = 0; i < count; ++i)
            if(nuc_j[i] != -1)
              if(qua_ij[i] >= max_qua) {
                max_qua = qua_ij[i];
                keep_id = i;
              }
        }
        Rcout << "site " << ref_j;
        Rcout << " is covered by 2 reads \n";
        // n_row[j] = 1;
        // IntegerVector tmp = {nuc_j[keep_id], nuc_j[keep_id], nuc_j[keep_id], nuc_j[keep_id]};
        // haplotype(j) = tmp;
        nuc_unique[j] = nuc_j[keep_id];
        nuc_count[j] = 1;
        continue;
      }
        // remove the ones with the lowest quality score
      if (qua_min < QUA_ERR && non_del > 1) {
        IntegerVector ind_que(count);
        for(i = 0; i < count; ++i)
          if(qua_ij[i] == qua_min)
            ind_que[i] = 1;
          // Rcout << ind_que << "\n";
          int enu = 0;
          for(i = 0; i < count; ++i)
            if(!ind_que[i])
              tmp[enu++] = nuc_j[i];
            if(enu == 0)
              tmp_nuc = nuc_j[Range(0, count - 1)];
            else
              tmp_nuc = tmp[Range(0, enu - 1)];
      } else
        tmp_nuc = nuc_j[Range(0, count - 1)];
      // Rcout << tmp_nuc << "\n";
      // /* skip the noncovered site */
      // if(count == 0)
      //   continue;
      /* get the nuc table at each site */
      // Rcout << j << "\n";
      nuc = unique_map(tmp_nuc);
      IntegerVector key = nuc["values"];
      IntegerVector val = nuc["lengths"];
      num = 0;

      /* remove some unlikely occurred nuc (notice the situation that only - appears) */
      if(cut_off != 0) {
        for(t = 0; t < val.length(); ++t)
          if(val[t] >= prop(j)) {
            keys(num) = key(t);
            vals(num++) = val(t);
          }
          nuc_unique[j] = keys[Range(0, num - 1)];
          nuc_count[j] = vals[Range(0, num - 1)];
      } else {
        num = key.size();
        nuc_unique[j] = key;
        nuc_count[j] = val;
      }

      IntegerVector hap_site = nuc_unique[j];
      IntegerVector sum_site = nuc_count[j]; // if there are three and minor are the same amount, order the minor two by the quality score, and eliminate the bad one
    }
  }
  IntegerVector record;
  if(cov_count != 0)
    record = record_cov[Range(0, cov_count - 1)];

  // compute the number of hidden states
  List ls = List::create(
    Named("nuc_unique") = nuc_unique, // unique nuc
    Named("nuc_count") = nuc_count,
    Named("p_tmax") = p_tmax, // max length at each t
    Named("n_t") = n_t,// no. emission reads at each time t
    Named("time_pos") = time_pos, // emission start position
    Named("n_in_t") = n_in_t, // reads index in each t
    Named("t_max") = t_max, //biggest t
    Named("cov_record") = record // record small cov loci
  );

  return ls;
}

// [[Rcpp::export]]
List get_hidden_state(List hmm_info, List dat_info, CharacterVector uni_alignment, List opt) {
  IntegerVector non_covered = dat_info["non_covered_site"];
  int hap_max_pos = dat_info["ref_length_max"];
  int hap_min_pos = dat_info["ref_start"];
  int hap_length = hap_max_pos - hap_min_pos;
  List nuc_unique = hmm_info["nuc_unique"];
  List nuc_count = hmm_info["nuc_count"];
  int t_max = hmm_info["t_max"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];

  unsigned int j, t;
  List haplotype(hap_length);
  IntegerVector n_row(hap_length);
  IntegerVector pos_possibility(hap_length);
  IntegerVector undecided_pos(hap_length);
  int emprical = opt["emprical"];
  int dbheter = opt["dbheter"];

  unsigned int more_than1 = 0;
  for (j = 0; j < hap_length; ++j) {
    if(!non_covered[j]) {
      unsigned int ref_j = j + hap_min_pos;
      IntegerVector hap_site = nuc_unique[j];
      IntegerVector sum_site = nuc_count[j];
      // Rcout << "j: " << j << ": " << hap_site << "|" << sum_site << "\n";
      unsigned int num = hap_site.size();
      if(num == 1) {
        n_row[j] = 1;
        IntegerVector tmp = {hap_site[0], hap_site[0], hap_site[0], hap_site[0]};
        haplotype(j) = tmp;
      } else {
        List out = sbs_state(num, ref_j, hap_site, sum_site, uni_alignment, opt);
        n_row[j] = out["n_row"];
        List hap_temp = out["haplotype"];
        if(n_row[j] == 1 || emprical)
          haplotype(j) = hap_temp[0];
        if(n_row[j] > 1) {
          undecided_pos(more_than1) = j; // relative to the haplotype currently try to infer (alignment start from 0)
          if(emprical)
            pos_possibility(more_than1) = n_row[j];
          more_than1++;
        }
      }
    } else {
      n_row[j] = 1;
      Rcout << "site " << j;
      Rcout << " is non_covered, cannot genotype \n";
      IntegerVector tmp = {-1, -1, -1, -1};
      haplotype(j) = tmp;
    }
  }
  if(!more_than1)
    stop("Filtering did not find any SNPs (neither homeologous or allelic)\n");
  // find the number of hidden states at each t
  List ls;
  IntegerVector num_states(t_max, 1);
  for(t = 0; t < t_max; ++t)
    for(j = time_pos[t] - hap_min_pos; j < time_pos[t] + p_tmax[t] - hap_min_pos; ++j)
      num_states[t] *= n_row[j];

  vector<int> db_loci;
  if(dbheter) { // loop through all the states
    for(j = 0; j < more_than1; ++j) {
      int idx = undecided_pos[j];
      IntegerMatrix haps = haplotype[idx];
      if((haps(5, 0) != haps(5, 1)) && (haps(5, 2) != haps(5, 3)))
        db_loci.push_back(idx);
    }
  }

  if(emprical) {
    ls = List::create(Named("undecided_pos") = undecided_pos[Range(0, more_than1 - 1)],
                      Named("pos_possibility") = pos_possibility[Range(0, more_than1 - 1)], //number of possible snp site
                                                                Named("n_row") = n_row, // number of possible genome at each site
                                                                Named("num_states") = num_states,// number of states at each t
                                                                Named("hidden_states") = haplotype,  // hidden_states at j
                                                                Named("db_loci") = db_loci);
  }
  else
    ls = List::create(Named("undecided_pos") = undecided_pos[Range(0, more_than1 - 1)],
                      Named("num_states") = num_states, // need to be adjusted
                      Named("n_row") = n_row,
                      Named("hidden_states") = haplotype,
                      Named("db_loci") = db_loci);
  return ls;
}

// [[Rcpp::export]]
List linkage_info(List dat_info, IntegerVector undecided_pos, CharacterVector uni_alignment,
                  List nuc_count, List nuc_unique) {
  IntegerMatrix ref_index = dat_info["ref_idx"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector index = dat_info["start_id"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector length = dat_info["length"];
  int hap_min_pos = dat_info["ref_start"];
  int n_observation = dat_info["n_observation"];

  List coverage;
  List deletion_snp; // record read with deletion at snp site
  IntegerMatrix link(n_observation, undecided_pos.size());
  IntegerVector new_nuc = clone(obs);
  unsigned int i, j;
  int idx;

  for (j = 0; j < undecided_pos.size(); ++j) {
    // Rcout<< "\n" << undecided_pos[j] << "\n";
    unsigned int ref_j = undecided_pos[j] + hap_min_pos;
    IntegerVector count = nuc_count[undecided_pos[j]];
    IntegerVector nucs = nuc_unique[undecided_pos[j]];

    if(nucs[0] == -1) {
      nucs.erase(0);
      count.erase(0);
    }
    vector<int> kept_nuc;
    int flag = 0;
    if(count.size() > 2) {
      flag = 1;
      int id = which_max(count);
      kept_nuc.push_back(nucs[id]);
      nucs.erase(id);
      count.erase(id);
      id = which_max(count);
      kept_nuc.push_back(nucs[id]);
    }
    IntegerVector tmp;
    vector<int> n_reads;
    for (i = 0; i < n_observation; i++) {
      if (ref_pos[index[i]] <= ref_j && ref_j <= ref_pos[index[i] + length[i] - 1]) {
        idx = ref_index(i, ref_j - hap_min_pos); // read pos, start from 0
        n_reads.push_back(i);
        if (idx != -1) {
          int in_nuc = obs[index[i] + idx];
          // Rcout << i << ":nuc:" << in_nuc << "|";
          if((flag && in_nuc != kept_nuc[0] && in_nuc != kept_nuc[1]) ||
             (!flag && in_nuc != nucs[0] && in_nuc != nucs[1])) // reads with third nuc
            in_nuc = -1;
          if(in_nuc == -1) {
            new_nuc[index[i] + idx] = -1;
            tmp.push_back(i);
            n_reads.pop_back();
          }
          link(i, j) = in_nuc;
        } else if (uni_alignment[ref_j] == "M" && idx == -1) {
          link(i, j) = -1; // meaning deletion
        } else
          link(i, j) = 4; // meaning gaps in the genomes
      } else
        link(i, j) = -1; // meaning not covered
    }
    if(tmp.size() != 0){
      tmp.push_back(undecided_pos[j]);
      deletion_snp.push_back(tmp);
    }
    coverage.push_back(n_reads);
  }
  List ls = List::create(Named("link") = link,
                         Named("deletion_snp") = deletion_snp,
                         Named("coverage") = coverage,
                         Named("obs") = new_nuc);
  return(ls);
}












