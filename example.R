library("foreach")
library("doParallel")
library('MASS')
#source('uHDSet_package.R')
library('uHDSet')

n=100      # number of indviduals
p=200      # number of SNPs
X=matrix(rnorm(n*p),n,p)
y=matrix(rnorm(n),n,1)
Kinship=0.5*diag(n)             # kinship matrix

T1=proc.time()
res=uHDSet(y,X,Kin,p_norm=0.3,option="typical",k=2,permute=100000)     # non-parallelized 
T_lapse=proc.time()-T1
T_lapse
P.uFineMap=res$P.uFineMap       # uFineMape test (marker wise test)
P.uHDSet=res$P.uHDSet           # uHDSet test (regional test)

T1=proc.time()
res_para=uHDSet_parallel(y,X,Kin,p_norm=0.3,option="typical",k=2,permute=100000)   # parallelized version
T_lapse=proc.time()-T1
T_lapse