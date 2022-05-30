// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// loglik_cpp_R6
double loglik_cpp_R6(arma::vec V, const arma::mat& betahat, const arma::mat& shat2, const arma::vec& prior_weight, const int nancestry, arma::uvec diag_index);
RcppExport SEXP _meSuSie_loglik_cpp_R6(SEXP VSEXP, SEXP betahatSEXP, SEXP shat2SEXP, SEXP prior_weightSEXP, SEXP nancestrySEXP, SEXP diag_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type shat2(shat2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior_weight(prior_weightSEXP);
    Rcpp::traits::input_parameter< const int >::type nancestry(nancestrySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type diag_index(diag_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_cpp_R6(V, betahat, shat2, prior_weight, nancestry, diag_index));
    return rcpp_result_gen;
END_RCPP
}
// mvlmm_reg
SEXP mvlmm_reg(arma::mat betahat, arma::mat shat2, arma::mat V_mat);
RcppExport SEXP _meSuSie_mvlmm_reg(SEXP betahatSEXP, SEXP shat2SEXP, SEXP V_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type shat2(shat2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V_mat(V_matSEXP);
    rcpp_result_gen = Rcpp::wrap(mvlmm_reg(betahat, shat2, V_mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_meSuSie_loglik_cpp_R6", (DL_FUNC) &_meSuSie_loglik_cpp_R6, 6},
    {"_meSuSie_mvlmm_reg", (DL_FUNC) &_meSuSie_mvlmm_reg, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_meSuSie(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
