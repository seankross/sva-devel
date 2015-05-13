#include <tgmath.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

NumericVector icomp(NumericVector v, int skip) {
    NumericVector vv(v.size()-1);
    for (int i=0; i < skip; i++) {
        vv[i] = v[i];
    }
    for (int i=skip+1; i < v.size(); i++) {
        vv[i-1] = v[i];
    }
    return vv;
}

NumericMatrix submat(NumericMatrix m, NumericVector cols) {
    NumericMatrix msub(m.nrow(), cols.size());

    for (int i=0; i < cols.size(); i++) {
        msub(_,i) = m(_, cols(i));
    }

    return msub;
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericMatrix calculateIntEprior(NumericMatrix sDat, NumericMatrix gammaHatMatrix,
    NumericMatrix deltaHatMatrix, List batches, int numBatches) {

    int length = gammaHatMatrix(0, _).size();
    NumericMatrix gammaDeltaStar(numBatches*2, length);

    for (int batchNum = 1; batchNum <= numBatches; batchNum++) {
        NumericVector currentBatches = as<NumericVector>(batches[batchNum-1]) - 1;
        NumericMatrix currentDat = submat(sDat, currentBatches);
        NumericVector gammaHat = gammaHatMatrix(batchNum-1,_);
        NumericVector deltaHat = deltaHatMatrix(batchNum-1,_);

        int numRows = gammaHat.size();
        NumericVector gammaStar(numRows), deltaStar(numRows);

        for (int i = 0; i < sDat.nrow(); i++) {
            // TODO make icomp more efficient by not rebuilding each time
            NumericVector gamma = icomp(gammaHat, i);
            NumericVector delta = icomp(deltaHat, i);

            NumericVector x = currentDat(i,_);
            double n = x.size();

            // Create matrix of all 1s
            vec j(n);
            j = j+1;

            NumericMatrix dat(numRows-1,n);
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

            long double zzz = sum(LH);
            printf("sum(LH): %Le\n", zzz);

            gammaStar[i] = sum(gamma * LH) / sum(LH);
            deltaStar[i] = sum(delta * LH) / sum(LH);
            printf("gammaStar: %f\n", gammaStar[i]);
        }

        gammaDeltaStar(batchNum-1,_) = clone(gammaStar);
        gammaDeltaStar(((batchNum-1)+numBatches),_) = clone(deltaStar);
    }

    // the first numBatches rows are gamma.star and
    // the second numBtaches rows are delta.star
    return gammaDeltaStar;
}
