// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// create_Q_grid
arma::mat create_Q_grid(int n);
RcppExport SEXP _VoA_create_Q_grid(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(create_Q_grid(n));
    return rcpp_result_gen;
END_RCPP
}
// Q
Rcpp::DataFrame Q(const arma::mat& Q_grid, const arma::mat& C_grid);
RcppExport SEXP _VoA_Q(SEXP Q_gridSEXP, SEXP C_gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q_grid(Q_gridSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C_grid(C_gridSEXP);
    rcpp_result_gen = Rcpp::wrap(Q(Q_grid, C_grid));
    return rcpp_result_gen;
END_RCPP
}
// seq_int
arma::uvec seq_int(long int a, long int b);
RcppExport SEXP _VoA_seq_int(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< long int >::type a(aSEXP);
    Rcpp::traits::input_parameter< long int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_int(a, b));
    return rcpp_result_gen;
END_RCPP
}
// rnd_idx
arma::uvec rnd_idx(int n);
RcppExport SEXP _VoA_rnd_idx(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rnd_idx(n));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_samples
std::vector<arma::uvec> bootstrap_samples(int n, int MC);
RcppExport SEXP _VoA_bootstrap_samples(SEXP nSEXP, SEXP MCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type MC(MCSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_samples(n, MC));
    return rcpp_result_gen;
END_RCPP
}
// arma_copula_mc
arma::mat arma_copula_mc(arma::vec Rx, arma::vec Ry, int MC);
RcppExport SEXP _VoA_arma_copula_mc(SEXP RxSEXP, SEXP RySEXP, SEXP MCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Rx(RxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ry(RySEXP);
    Rcpp::traits::input_parameter< int >::type MC(MCSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_copula_mc(Rx, Ry, MC));
    return rcpp_result_gen;
END_RCPP
}
// arma_copula
arma::mat arma_copula(arma::vec Rx, arma::vec Ry);
RcppExport SEXP _VoA_arma_copula(SEXP RxSEXP, SEXP RySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Rx(RxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ry(RySEXP);
    rcpp_result_gen = Rcpp::wrap(arma_copula(Rx, Ry));
    return rcpp_result_gen;
END_RCPP
}
// calculate_copula_grid
arma::mat calculate_copula_grid(const arma::vec X, const arma::vec Y);
RcppExport SEXP _VoA_calculate_copula_grid(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_copula_grid(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// calculate_copula_mc_grid
arma::mat calculate_copula_mc_grid(const arma::vec& X, const arma::vec& Y, int MC);
RcppExport SEXP _VoA_calculate_copula_mc_grid(SEXP XSEXP, SEXP YSEXP, SEXP MCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type MC(MCSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_copula_mc_grid(X, Y, MC));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VoA_create_Q_grid", (DL_FUNC) &_VoA_create_Q_grid, 1},
    {"_VoA_Q", (DL_FUNC) &_VoA_Q, 2},
    {"_VoA_seq_int", (DL_FUNC) &_VoA_seq_int, 2},
    {"_VoA_rnd_idx", (DL_FUNC) &_VoA_rnd_idx, 1},
    {"_VoA_bootstrap_samples", (DL_FUNC) &_VoA_bootstrap_samples, 2},
    {"_VoA_arma_copula_mc", (DL_FUNC) &_VoA_arma_copula_mc, 3},
    {"_VoA_arma_copula", (DL_FUNC) &_VoA_arma_copula, 2},
    {"_VoA_calculate_copula_grid", (DL_FUNC) &_VoA_calculate_copula_grid, 2},
    {"_VoA_calculate_copula_mc_grid", (DL_FUNC) &_VoA_calculate_copula_mc_grid, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_VoA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
