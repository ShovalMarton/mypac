// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// kernel
Rcpp::NumericVector kernel(Rcpp::NumericVector x, int krnl_type);
RcppExport SEXP _mypac_kernel(SEXP xSEXP, SEXP krnl_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type krnl_type(krnl_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel(x, krnl_type));
    return rcpp_result_gen;
END_RCPP
}
// getBeta0
double getBeta0(const arma::vec X, const arma::vec Y, const arma::vec w);
RcppExport SEXP _mypac_getBeta0(SEXP XSEXP, SEXP YSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(getBeta0(X, Y, w));
    return rcpp_result_gen;
END_RCPP
}
// ifelseYt
arma::vec ifelseYt(arma::vec x, double num);
RcppExport SEXP _mypac_ifelseYt(SEXP xSEXP, SEXP numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type num(numSEXP);
    rcpp_result_gen = Rcpp::wrap(ifelseYt(x, num));
    return rcpp_result_gen;
END_RCPP
}
// ifelseNt
arma::vec ifelseNt(arma::vec x, arma::vec delta1, double num);
RcppExport SEXP _mypac_ifelseNt(SEXP xSEXP, SEXP delta1SEXP, SEXP numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< double >::type num(numSEXP);
    rcpp_result_gen = Rcpp::wrap(ifelseNt(x, delta1, num));
    return rcpp_result_gen;
END_RCPP
}
// bw_RT_bool
double bw_RT_bool(arma::vec x, arma::vec y, int krnl, double deg);
RcppExport SEXP _mypac_bw_RT_bool(SEXP xSEXP, SEXP ySEXP, SEXP krnlSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type krnl(krnlSEXP);
    Rcpp::traits::input_parameter< double >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(bw_RT_bool(x, y, krnl, deg));
    return rcpp_result_gen;
END_RCPP
}
// bw_RT
double bw_RT(arma::vec x, arma::vec y, int krnl, double deg);
RcppExport SEXP _mypac_bw_RT(SEXP xSEXP, SEXP ySEXP, SEXP krnlSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type krnl(krnlSEXP);
    Rcpp::traits::input_parameter< double >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(bw_RT(x, y, krnl, deg));
    return rcpp_result_gen;
END_RCPP
}
// est_h
arma::mat est_h(arma::mat data, int krnl, arma::vec S_s, arma::vec T1_s);
RcppExport SEXP _mypac_est_h(SEXP dataSEXP, SEXP krnlSEXP, SEXP S_sSEXP, SEXP T1_sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type krnl(krnlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_s(S_sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type T1_s(T1_sSEXP);
    rcpp_result_gen = Rcpp::wrap(est_h(data, krnl, S_s, T1_s));
    return rcpp_result_gen;
END_RCPP
}
// Lambda_given_s
arma::mat Lambda_given_s(arma::vec S_s, arma::vec T1_s, arma::mat est_int1);
RcppExport SEXP _mypac_Lambda_given_s(SEXP S_sSEXP, SEXP T1_sSEXP, SEXP est_int1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type S_s(S_sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type T1_s(T1_sSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type est_int1(est_int1SEXP);
    rcpp_result_gen = Rcpp::wrap(Lambda_given_s(S_s, T1_s, est_int1));
    return rcpp_result_gen;
END_RCPP
}
// final_CIF
arma::mat final_CIF(arma::mat Surv_given_s, arma::vec T1_s, arma::vec dval);
RcppExport SEXP _mypac_final_CIF(SEXP Surv_given_sSEXP, SEXP T1_sSEXP, SEXP dvalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Surv_given_s(Surv_given_sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type T1_s(T1_sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dval(dvalSEXP);
    rcpp_result_gen = Rcpp::wrap(final_CIF(Surv_given_s, T1_s, dval));
    return rcpp_result_gen;
END_RCPP
}
// new_CIF1
arma::mat new_CIF1(arma::mat data, int krnl, arma::vec S_s, arma::vec T1_s);
RcppExport SEXP _mypac_new_CIF1(SEXP dataSEXP, SEXP krnlSEXP, SEXP S_sSEXP, SEXP T1_sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type krnl(krnlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_s(S_sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type T1_s(T1_sSEXP);
    rcpp_result_gen = Rcpp::wrap(new_CIF1(data, krnl, S_s, T1_s));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _mypac_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _mypac_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _mypac_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _mypac_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mypac_kernel", (DL_FUNC) &_mypac_kernel, 2},
    {"_mypac_getBeta0", (DL_FUNC) &_mypac_getBeta0, 3},
    {"_mypac_ifelseYt", (DL_FUNC) &_mypac_ifelseYt, 2},
    {"_mypac_ifelseNt", (DL_FUNC) &_mypac_ifelseNt, 3},
    {"_mypac_bw_RT_bool", (DL_FUNC) &_mypac_bw_RT_bool, 4},
    {"_mypac_bw_RT", (DL_FUNC) &_mypac_bw_RT, 4},
    {"_mypac_est_h", (DL_FUNC) &_mypac_est_h, 4},
    {"_mypac_Lambda_given_s", (DL_FUNC) &_mypac_Lambda_given_s, 3},
    {"_mypac_final_CIF", (DL_FUNC) &_mypac_final_CIF, 3},
    {"_mypac_new_CIF1", (DL_FUNC) &_mypac_new_CIF1, 4},
    {"_mypac_rcpparma_hello_world", (DL_FUNC) &_mypac_rcpparma_hello_world, 0},
    {"_mypac_rcpparma_outerproduct", (DL_FUNC) &_mypac_rcpparma_outerproduct, 1},
    {"_mypac_rcpparma_innerproduct", (DL_FUNC) &_mypac_rcpparma_innerproduct, 1},
    {"_mypac_rcpparma_bothproducts", (DL_FUNC) &_mypac_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_mypac(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
