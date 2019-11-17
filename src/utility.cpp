
#include <Rcpp.h>
using namespace Rcpp;

template<int RTYPE>
NumericMatrix normalize_cols_cpp_imp(Matrix<RTYPE> M) {
  NumericMatrix result(M.nrow(), M.ncol());
  for (int j = 0; j < M.ncol(); j++) {
    double mean = 0;
    for (int i = 0; i < M.nrow(); i++) { mean += M(i,j); }
    mean /= M.nrow();
    for (int i = 0; i < M.nrow(); i++) { result(i,j) = M(i,j) - mean; }
    double norm = 0;
    for (int i = 0; i < M.nrow(); i++) { norm += result(i,j) * result(i,j); }
    norm = 1 / sqrt(norm);
    for (int i = 0; i < M.nrow(); i++) { result(i,j) *= norm; }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix normalize_cols_cpp(SEXP M) {
  switch (TYPEOF(M)) {
    case INTSXP: return normalize_cols_cpp_imp<INTSXP>(M);
    case REALSXP: return normalize_cols_cpp_imp<REALSXP>(M);
  }
}
