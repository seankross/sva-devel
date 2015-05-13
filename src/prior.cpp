#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List paraStar(NumericMatrix sData, NumericMatrix gammaHat, NumericMatrix deltaHat,
            NumericVector gammaBar, NumericVector t2, NumericVector aPrior, 
            NumericVector bPrior, Function itSol, List batches) {
  int numBatches = batches.size();
  int sDataRows = sData.nrow();
  NumericMatrix gammaStar(numBatches, sDataRows);
  NumericMatrix deltaStar(numBatches, sDataRows);
  NumericMatrix temp;
  for(int i = 0; i < numBatches; i++){
    NumericVector currentBatches = as<NumericVector>(batches[i]) - 1;
    int currentBatchLength = currentBatches.size();
    
    NumericMatrix sliced(sDataRows, currentBatchLength);
    for(int j = 0; j < currentBatchLength; j++){
      sliced(_,j) = sData(_,currentBatches[j]);
    }
    temp = itSol(sliced, gammaHat(i,_), 
              deltaHat(i,_) , gammaBar[i] , t2[i], aPrior[i], bPrior[i]);
    gammaStar(i,_) = temp(0,_);
    deltaStar(i,_) = temp(1,_);
  }
  
  return Rcpp::List::create(Rcpp::Named("gammaStar") = gammaStar,
                          Rcpp::Named("deltaStar") = deltaStar);
}

// [[Rcpp::export]]
List nonParaStar(NumericMatrix sData, NumericMatrix gammaHat, NumericMatrix deltaHat,
            Function intEPrior, List batches) {
  int numBatches = batches.size();
  int sDataRows = sData.nrow();
  NumericMatrix gammaStar(numBatches, sDataRows);
  NumericMatrix deltaStar(numBatches, sDataRows);
  NumericMatrix temp;
  for(int i = 0; i < numBatches; i++){
    NumericVector currentBatches = as<NumericVector>(batches[i]) - 1;
    int currentBatchLength = currentBatches.size();
    
    NumericMatrix sliced(sDataRows, currentBatchLength);
    for(int j = 0; j < currentBatchLength; j++){
      sliced(_,j) = sData(_,currentBatches[j]);
    }
    temp = intEPrior(sliced, gammaHat(i,_), deltaHat(i,_));
    gammaStar(i,_) = temp(0,_);
    deltaStar(i,_) = temp(1,_);
  }
  
  return Rcpp::List::create(Rcpp::Named("gammaStar") = gammaStar,
                          Rcpp::Named("deltaStar") = deltaStar);
}