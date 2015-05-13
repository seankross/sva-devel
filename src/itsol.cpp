#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector postmeanCPP(NumericVector gammaHat, NumericVector gammaBar,
    NumericVector n, NumericVector dStar, NumericVector t2, int batchNum) {

    int length = gammaHat.size();
    NumericVector gammaNew(length);

    for (int i = 0; i < length; i++) {
        gammaNew[i] =
            ((t2(batchNum-1) * n(i) * gammaHat(i)) +
                (dStar(i) * gammaBar(batchNum-1))
            ) / ((t2(batchNum-1) * n(i)) + dStar(i));
    }

    return gammaNew;
}

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

// [[Rcpp::export]]
NumericMatrix submat(NumericMatrix m, NumericVector cols) {
    NumericMatrix msub(m.nrow(), cols.size());

    for (int i=0; i < cols.size(); i++) {
        msub(_,i) = m(_, cols(i));
    }

    return msub;
}

// [[Rcpp::export]]
NumericMatrix calculateItSol(NumericMatrix sDat, NumericMatrix gammaHat, NumericMatrix deltaHat, NumericVector gammaBar, NumericVector t2,
    NumericVector aPrior, NumericVector bPrior, List batches, int numBatches, double conv = 0.0001) {

    int length = gammaHat(0, _).size();
    NumericMatrix gammaDeltaStar(numBatches*2, length);

    for (int batchNum = 1; batchNum <= numBatches; batchNum++) {
        NumericVector currentBatches = as<NumericVector>(batches[batchNum-1]) - 1;
        NumericMatrix currentDat = submat(sDat, currentBatches);

        NumericVector n(currentDat.nrow());
        for (int i = 0; i < currentDat.nrow(); i++) {
            n[i] = sum(!is_na(currentDat(i,_)));
        }

        int count = 0;
        double change = 1;
        NumericVector gammaNew, deltaNew;
        NumericVector gammaOld = gammaHat(batchNum-1,_);
        NumericVector deltaOld = deltaHat(batchNum-1,_);

        while(change > conv) {

            gammaNew = postmeanCPP(gammaHat(batchNum-1,_), gammaBar, n, deltaOld, t2, batchNum);

            // TODO need to check for NA?
            NumericVector sum2(currentDat.nrow());
            NumericVector diff;
            for (int i = 0; i < currentDat.nrow(); i++) {
                diff = (currentDat(i,_) - gammaNew[i]);
                sum2[i] = sum(diff*diff);
            }

            deltaNew = postvarCPP(sum2, n, aPrior, bPrior, batchNum);

            double gammaMax = max(abs(gammaNew-gammaOld) / gammaOld);
            double deltaMax = max(abs(deltaNew-deltaOld) / deltaOld);
            change = std::max(gammaMax,deltaMax);

            gammaOld = gammaNew;
            deltaOld = deltaNew;
            count++;
        }

        printf("Batch %d took %d iterations until convergence\n", batchNum, count);

        gammaDeltaStar(batchNum-1,_) = clone(gammaNew);
        gammaDeltaStar(((batchNum-1)+numBatches),_) = clone(deltaNew);
    }

    // the first numBatches rows are gamma.star and
    // the second numBtaches rows are delta.star
    return gammaDeltaStar;
}
