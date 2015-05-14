#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector sva_calculateDeltaHat(NumericMatrix sData, List batches) {
  int numBatches = batches.size();
  int sDataRows = sData.nrow();
  NumericMatrix result(numBatches, sDataRows);
  
  for(int i = 0; i < numBatches; i++){
    NumericVector currentBatches = as<NumericVector>(batches[i]) - 1;
    int currentBatchLength = currentBatches.size();
    
    NumericMatrix sliced(sDataRows, currentBatchLength);
    for(int j = 0; j < currentBatchLength; j++){
      sliced(_,j) = sData(_,currentBatches[j]);
    }
    
    for(int j = 0; j < sDataRows; j++){
      //NumericVector to_var = as<NumericVector>(sliced(j,_));
      result(i, j) = var(sliced(j,_));
    }
    
  }
  
  return result;
}
