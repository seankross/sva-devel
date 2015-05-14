
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

install_github("seankross/sva-devel")
library(sva)

microbenchmark(
  combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plot=FALSE), times = 5)
```