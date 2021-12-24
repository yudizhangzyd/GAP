#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>

#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector sort_ind(CharacterVector input, int n) {
  int i;
  NumericVector X(n);
  NumericVector re(n);
  
  for (i = 0; i < n; ++i) {
    if (input[i] == "A") {
      X[i] = 1;
    }
    if (input[i] == "C") {
      X[i] = 2;
    }
    if (input[i] == "G") {
      X[i] = 3;
    }
    if (input[i] == "T") {
      X[i] = 4;
    }
  }
  // Although using sort_index is faster, it will reorder the orginal sequence
  vec x = as<vec>(X);
  uvec indices = arma::stable_sort_index(x);
  re = as<NumericVector>(wrap(indices)); 
  return(re);
}

// [[Rcpp::export]]
NumericMatrix formDesignMat_c(DataFrame dat, int N, int PD_LENGTH) {
    IntegerVector qua = dat["qua"];
    IntegerVector ref_pos = dat["ref_pos"];
    IntegerVector read_pos = dat["read_pos"];
    CharacterVector hap_nuc = dat["hap_nuc"];

    NumericMatrix X(N, PD_LENGTH);

    for (int i = 0; i < N; ++i) {
    	X(i, 0) = 1;
    	X(i, 1) = read_pos[i];
    	X(i, 2) = ref_pos[i];
    	X(i, 3) = qua[i];
    	if (hap_nuc[i] == "C") {
    		X(i, 4) = 1;
    		X(i, 7) = qua[i];
    	}
    	if (hap_nuc[i] == "G") {
    		X(i, 5) = 1;
    		X(i, 8) = qua[i];
    	}
    	if (hap_nuc[i] == "T") {
    		X(i, 6) = 1;
    		X(i, 9) = qua[i];
    	}
    }
    return(X);
}

// 
// DataFrame sort_df(DataFrame dat, vec index) {
//   IntegerVector qua = dat["qua"];
//   IntegerVector obs = dat["nuc"];
//   IntegerVector ref_pos = dat["ref_pos"];
//   IntegerVector read_pos = dat["read_pos"];
//   IntegerVector mode = dat["mode"];
//   IntegerVector hap_nuc = dat["hap_nuc"];
// 
// IntegerVector q(N);
// IntegerVector ref_p(N);
// IntegerVector read_p(N);
// IntegerVector hap_n(N);
// IntegerVector q_h(N);
// //   
//   vec q = as<vec>(qua);
//   q = q(index);
//   o = o(index);
//   
//   
//   return DataFrame::create(
//     
//   );
// }
// 
// int compare_matrix(NumericMatrix x, NumericMatrix y, 
//                    unsigned int ncol, unsigned int nrow) {
//   int flag = 0;
//   for(unsigned int i = 0; i < nrow; ++i)
//     for(unsigned int j = 0; j < ncol; ++j)
//       if(x(i, j) != y(i, j))
//         flag = 1;
//       return(flag);
// }




