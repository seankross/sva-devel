#include <Rcpp.h>
using namespace Rcpp;

//postmean <- function(g.hat,g.bar,n,d.star,t2){
//    (t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)

// [[Rcpp::export]]
NumericVector postmeanCPP(NumericMatrix gammaHat, NumericVector gammaBar,
    NumericVector n, NumericMatrix dStar, NumericVector t2, int batchNum) {

    int length = gammaHat(batchNum-1, _).size();
    NumericVector gammaNew(length);

    for (int i = 0; i < length; i++) {
        gammaNew[i] =
            ((t2(batchNum-1) * n(i) * gammaHat(batchNum-1, i)) +
                (dStar(batchNum-1, i) * gammaBar(batchNum-1))
            ) / ((t2(batchNum-1) * n(i)) + dStar(batchNum-1, i));
    }

    return gammaNew;
}

//postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}

// [[Rcpp::export]]
NumericVector postvarCPP(NumericVector sum2, NumericVector n, 
    NumericVector a, NumericVector b, int batchNum) {
  
    int length = sum2.size();
    NumericVector deltaNew(length);
    
    for (int i = 0; i < length; i++) {
        deltaNew[i] = (0.5 * sum2[i] + b[batchNum-1]) / 
                      ((n[i]/2) + a[batchNum-1] -1);
    }
    
    return deltaNew;
  
}
