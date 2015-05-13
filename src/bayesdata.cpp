#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int bayesData(NumericMatrix sData, List batches, NumericMatrix batchDesign,
                  NumericMatrix gammaStar, NumericMatrix deltaStar) {
  int numBatches = batches.size();
  int sDataRows = sData.nrow();
  int sDataCols = sData.ncol();
  NumericMatrix result(sDataRows, sDataCols);
  
  for(int i = 0; i < numBatches; i++){
    NumericVector currentBatches = as<NumericVector>(batches[i]) - 1;
    int currentBatchLength = currentBatches.size();
    
    NumericMatrix sliced(sDataRows, currentBatchLength);
    for(int j = 0; j < currentBatchLength; j++){
      sliced(_,j) = sData(_,currentBatches[j]);
    }
    (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))
    / (sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
    for(int j = 0; j < sDataRows; j++){
      //NumericVector to_var = as<NumericVector>(sliced(j,_));
      result(i, j) = var(sliced(j,_));
    }
    
  }
}
