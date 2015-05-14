#include <tgmath.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Move src into dst, skipping the ith element
void icomp(NumericVector src, NumericVector dst, int skip) {
    for (int i=0; i < skip; i++) {
        dst[i] = src[i];
    }
    for (int i=skip+1; i < src.size(); i++) {
        dst[i-1] = src[i];
    }
}

// Return a subset of matrix m containing only columns
// with indexes in cols
NumericMatrix submat(NumericMatrix m, NumericVector cols) {
    NumericMatrix msub(m.nrow(), cols.size());

    for (int i=0; i < cols.size(); i++) {
        msub(_,i) = m(_, cols(i));
    }

    return msub;
}

NumericVector postmean(NumericVector gammaHat, NumericVector gammaBar,
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

NumericVector postvar(NumericVector sum2, NumericVector n,
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
NumericMatrix sva_parametricAdjust(NumericMatrix sDat, NumericMatrix gammaHat, NumericMatrix deltaHat, NumericVector gammaBar, NumericVector t2,
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

            gammaNew = postmean(gammaHat(batchNum-1,_), gammaBar, n, deltaOld, t2, batchNum);

            NumericVector sum2(currentDat.nrow());
            NumericVector diff;
            for (int i = 0; i < currentDat.nrow(); i++) {
                diff = (currentDat(i,_) - gammaNew[i]);
                sum2[i] = sum(diff*diff);
            }

            deltaNew = postvar(sum2, n, aPrior, bPrior, batchNum);

            double gammaMax = max(abs(gammaNew-gammaOld) / gammaOld);
            double deltaMax = max(abs(deltaNew-deltaOld) / deltaOld);
            change = std::max(gammaMax,deltaMax);

            gammaOld = gammaNew;
            deltaOld = deltaNew;
            count++;
        }

        //printf("Batch %d took %d iterations until convergence\n", batchNum, count);

        gammaDeltaStar(batchNum-1,_) = clone(gammaNew);
        gammaDeltaStar(((batchNum-1)+numBatches),_) = clone(deltaNew);
    }

    // the first numBatches rows are gamma.star and
    // the second numBtaches rows are delta.star
    return gammaDeltaStar;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix sva_nonparametricAdjust(NumericMatrix sDat, NumericMatrix gammaHatMatrix,
    NumericMatrix deltaHatMatrix, List batches, int numBatches) {

    int length = gammaHatMatrix(0, _).size();
    NumericMatrix gammaDeltaStar(numBatches*2, length);

    for (int batchNum = 1; batchNum <= numBatches; batchNum++) {
        NumericVector currentBatches = as<NumericVector>(batches[batchNum-1]) - 1;
        NumericMatrix currentDat = submat(sDat, currentBatches);
        NumericVector gammaHat = gammaHatMatrix(batchNum-1,_);
        NumericVector deltaHat = deltaHatMatrix(batchNum-1,_);

        int numRows = gammaHat.size();
        NumericVector gamma(numRows-1);
        NumericVector delta(numRows-1);
        NumericMatrix dat(numRows-1,currentDat.ncol());

        NumericVector gammaStar(numRows), deltaStar(numRows);
        NumericVector ones(currentDat.ncol(), 1.0);
        vec j = as<NumericVector>(wrap(ones));

        for (int i = 0; i < sDat.nrow(); i++) {
            icomp(gammaHat, gamma, i);
            icomp(deltaHat, delta, i);

            NumericVector x = currentDat(i,_);
            double n = x.size();

            for (int k = 0; k < n; k++) {
                for (int l = 0; l < dat.nrow(); l++) {
                    dat(l,k) = x[k];
                }
            }

            for (int k = 0; k < n; k++) {
                NumericMatrix::Column col = dat(_,k);
                col = col - gamma;
            }

            mat resid2 (dat.begin(), dat.nrow(), dat.ncol(), false);
            resid2 = resid2 % resid2;
            vec sum2 = resid2 * j;
            NumericVector s2 = as<NumericVector>(wrap(sum2));

            NumericVector LH = (1/pow((2*M_PI*delta),(n/2)))*exp(-s2/(2*delta));
            LH[is_na(LH)] = 0;

            gammaStar[i] = sum(gamma * LH) / sum(LH);
            deltaStar[i] = sum(delta * LH) / sum(LH);
            //printf("gammaStar[%d]: %f\n", i, gammaStar[i]);
        }

        gammaDeltaStar(batchNum-1,_) = clone(gammaStar);
        gammaDeltaStar(((batchNum-1)+numBatches),_) = clone(deltaStar);
    }

    // the first numBatches rows are gamma.star and
    // the second numBtaches rows are delta.star
    return gammaDeltaStar;
}
