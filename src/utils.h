#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

class Permutation {
public:
  vector<vector<int> > permuteUnique(vector<int>& nums) {
    vector<vector<int> > res;
    sort(nums.begin(), nums.end());
    res.push_back(nums);
    while (next_permutation(nums.begin(), nums.end())) {
      res.push_back(nums);
    }
    return res;
  }
};

vector<vector<double> > cart_product_dbl (const vector<vector<double> > &v);
vector<vector<int> > cart_product (const vector<vector<int> > &v);
IntegerMatrix call_cart_product(IntegerVector len);
IntegerMatrix comb_element(List len, IntegerVector flag, unsigned int row);
List unique_map(const Rcpp::IntegerVector & v);
List hash_mat(IntegerMatrix x);
std::unordered_map<string, vector<int>> hash_mat2(IntegerMatrix x);
IntegerMatrix ss(IntegerMatrix X_, IntegerVector ind_, int is_row = 1);
arma::mat unique_rows(const arma::mat& m);
IntegerVector find_max(List ls, int n_obs);
IntegerMatrix call_permute(vector<int> a);
IntegerMatrix sort_mat(IntegerMatrix mat, int nrow, int ncol, unsigned int db_heter);
IntegerVector matrix2vec(IntegerMatrix m, const bool byrow = true);
void print_intmat(IntegerMatrix m);
void print_c_intvec(vector<int> a);
List hash_intvec(IntegerVector x);
IntegerVector sample_int(IntegerVector x, int size, bool replace,
                         NumericVector prob = NumericVector::create());
#endif
