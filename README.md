# sva + Rcpp

We re-implemented the ComBat algorithm using Rcpp and increased the speed of
processing the canonical dataset by over 10x. In our tests ComBat version
3.14.0 completed the batch correction in 5.86 seconds, where ComBat 3.15.1
(which includes Rcpp) did the same batch correction in 0.447 seconds. Our
major contributions are the files `src/calculateDeltaHapp.cpp` and
`src/helpers.cpp`. You can reproduce our speedup analysis with the script below.

```r
install.packages("devtools")
install.packages("microbenchmark")
source("http://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("bladderbatch")

library(devtools)
library(microbenchmark)
library(sva)
library(bladderbatch)
data(bladderdata)

pheno <- pData(bladderEset)

edata <- exprs(bladderEset)

batch <- pheno$batch

modcombat <- model.matrix(~1, data=pheno)

microbenchmark(
  combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE), times = 5)

detach("package:sva", unload=TRUE)
remove.packages("sva")

install_github("seankross/sva-devel", build_vignettes = FALSE, quick = TRUE)
library(sva)

microbenchmark(
  combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plot=FALSE)
  , times = 5)
```