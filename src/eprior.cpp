#include <tgmath.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// # Monte Carlo integration functients
// int.eprior <- function(sdat,g.hat,d.hat){
//     g.star <- d.star <- NULL
//     r <- nrow(sdat)
//     for(i in 1:r){
//     # numeric vector, all elements but i
//         g <- g.hat[-i]
//     # numeric vector, all elements but i
//         d <- d.hat[-i]
//     # numeric vector, one row of the table, no nas
//         x <- sdat[i,!is.na(sdat[i,])]
//     # length of x (usually 11)
//         n <- length(x)
//     # vector of eleven 1s [1 1 1 1 1 1 1 1 1 1 ]
//         j <- numeric(n)+1
//     # matrix with 22282 rows that are each exactly x
//         dat <- matrix(as.numeric(x),length(g),n,byrow=T)
//     # subtract g from each column, now everythings different!
//         resid2 <- (dat-g)^2
//     # matrix multiplication by j
//         sum2 <- resid2%*%j
//     # Vector of size 22282
//         LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
//     # remove nan
//         LH[LH=="NaN"]=0
//     # add single value to this vector at position i
//         g.star <- c(g.star,sum(g*LH)/sum(LH))
//         d.star <- c(d.star,sum(d*LH)/sum(LH))
//         if(i%%1000==0){cat(i,'\n')}
//         }
//     adjust <- rbind(g.star,d.star)
//     rownames(adjust) <- c("g.star","d.star")
//     adjust
//     }
    // MAIN
    // for (i in 1:n.batch){
    //   temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
    //   gamma.star <- rbind(gamma.star,temp[1,])
    //   delta.star <- rbind(delta.star,temp[2,])
    // }

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

            gammaStar[i] = sum(gamma * LH) / sum(LH);
            deltaStar[i] = sum(delta * LH) / sum(LH);
            printf("gamma[%d]: %f, delta: %f\n", i, gammaStar[i], deltaStar[i]);
        }

        gammaDeltaStar(batchNum-1,_) = clone(gammaStar);
        gammaDeltaStar(((batchNum-1)+numBatches),_) = clone(deltaStar);
    }

    // the first numBatches rows are gamma.star and
    // the second numBtaches rows are delta.star
    return gammaDeltaStar;
}
