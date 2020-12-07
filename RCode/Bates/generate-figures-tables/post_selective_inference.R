  library(glmnet)
#
#script for post-sel inference testing H0: sum(beta_A)=0
# i modified fixedLassoInf from selectiveInference pkg, adding argument sumtest
# if sumtest=TRUE, then first pvalue returned is pvalue for this H0 (a hack)
# the fixedLassonf pkg also returns CIs, but I don't look at  them here

options(error=dump.frames)

source("si.R")  #extracted R code from selectiveInference package


set.seed(344)
n=100
p=10
sigma=.8
x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

                                        # two scenarios

#beta=rep(0,p)  #full  null
 beta=c(2,-2,rep(0,p-2))  #some non-null coefs
 
 tr=which(beta!=0)
  lambda = 1
nsim=1000
pv=size=NULL


for(ii in 1:nsim){
   cat(ii)
   y=x%*%beta+sigma*rnorm(n)
   y=y-mean(y)
   gfit=glmnet(x,y,standardize=F)

     b = coef(gfit, s=lambda/n, exact=TRUE)[-1]
     
   if(sum(b!=0)>0){
     out = fixedLassoInf(x,y,b,lambda,sigma=sigma,sumtest=T)
     
     screened=T

      if(length(tr)>0) screened=sum(is.na(match(tr,out$vars)))==0
      if(screened) pv=c(pv,out$pv[1])
      size=c(size,length(out$vars))
}}


np=length(pv)

pdf(file="pvalues.pdf",width=5,height=5)
plot( (1:np)/np,sort(pv),xlab="Expected p-value",ylab="Observed p-value")
abline(0,1)
dev.off()

