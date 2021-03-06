// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calculateDeltaHat
NumericVector calculateDeltaHat(NumericMatrix sData, List batches);
RcppExport SEXP sva_calculateDeltaHat(SEXP sDataSEXP, SEXP batchesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type sData(sDataSEXP);
    Rcpp::traits::input_parameter< List >::type batches(batchesSEXP);
    __result = Rcpp::wrap(calculateDeltaHat(sData, batches));
    return __result;
END_RCPP
}
// parametricAdjust
NumericMatrix parametricAdjust(NumericMatrix sDat, NumericMatrix gammaHat, NumericMatrix deltaHat, NumericVector gammaBar, NumericVector t2, NumericVector aPrior, NumericVector bPrior, List batches, int numBatches, double conv);
RcppExport SEXP sva_parametricAdjust(SEXP sDatSEXP, SEXP gammaHatSEXP, SEXP deltaHatSEXP, SEXP gammaBarSEXP, SEXP t2SEXP, SEXP aPriorSEXP, SEXP bPriorSEXP, SEXP batchesSEXP, SEXP numBatchesSEXP, SEXP convSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type sDat(sDatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gammaHat(gammaHatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type deltaHat(deltaHatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gammaBar(gammaBarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type aPrior(aPriorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bPrior(bPriorSEXP);
    Rcpp::traits::input_parameter< List >::type batches(batchesSEXP);
    Rcpp::traits::input_parameter< int >::type numBatches(numBatchesSEXP);
    Rcpp::traits::input_parameter< double >::type conv(convSEXP);
    __result = Rcpp::wrap(parametricAdjust(sDat, gammaHat, deltaHat, gammaBar, t2, aPrior, bPrior, batches, numBatches, conv));
    return __result;
END_RCPP
}
// nonparametricAdjust
NumericMatrix nonparametricAdjust(NumericMatrix sDat, NumericMatrix gammaHatMatrix, NumericMatrix deltaHatMatrix, List batches, int numBatches);
RcppExport SEXP sva_nonparametricAdjust(SEXP sDatSEXP, SEXP gammaHatMatrixSEXP, SEXP deltaHatMatrixSEXP, SEXP batchesSEXP, SEXP numBatchesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type sDat(sDatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gammaHatMatrix(gammaHatMatrixSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type deltaHatMatrix(deltaHatMatrixSEXP);
    Rcpp::traits::input_parameter< List >::type batches(batchesSEXP);
    Rcpp::traits::input_parameter< int >::type numBatches(numBatchesSEXP);
    __result = Rcpp::wrap(nonparametricAdjust(sDat, gammaHatMatrix, deltaHatMatrix, batches, numBatches));
    return __result;
END_RCPP
}
