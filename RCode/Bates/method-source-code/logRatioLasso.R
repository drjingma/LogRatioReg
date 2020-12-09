#NOTE; this function and the predict function should be made more
#  efficient by not keeping and working with the full pairwise
#  matrix, which is mostly zero

# version that does just FS, with proper CV

logRatioLasso=function(x,y,family=c("gaussian","binomial"),delta=1,initial.fit=c("constrainedLasso","lasso"),en.alpha=1,nlambda=50,
    lambda=NULL,lambda.min.ratio=.05,maxvar=max(50,trunc(ncol(x)/4)),second.fit=c("FS","constrainedFS"),maxsteps=min(50,n),verbose=1){
    #fit  logs at first stage and then log ratios at second stage
    # fits over a lambda path at first stage
    # family is "gaussian" or "binomial"
    # initial.fit can be lasso, or constrained lasso (sum to 0), both with parameters lam  and  en par en.alpha
    #lambda.min.ratio= depth of path for call to glmnet
    # maxvar-  max # vars allowed from first stage model
     # maxsteps-  max # terms allowed in 2nd stage model
    #nlambda=# of lambdas considered in first stage
    # optional list of lambda values; overrides nlambda
    #second fit can be  forw stepwise (usual or constrained)
 
 
    #delta= amount added to x before taking logs
    #verbose= reporting flag 0=quiet, 1=path only, 2=more info
    #returns:
    # beta=coef from first stage, one model per column (nlam columns)
    # ind0= indices of predictors in each first stage model, in terms of original pred numbering
    # glmnet.obj-  object returned from first stage fit
    #  fsout= list of  objects returned by each 2nd stage fs model
    this.call=match.call()
    family = match.arg(family)
    initial.fit=match.arg(initial.fit)
    second.fit=match.arg(second.fit)
    n=length(y)
    meany=NULL
    if(family=="gaussian") {
        meany=mean(y)
        y=y-mean(y)
    }
    x=x+delta
    n=nrow(x)
    p=ncol(x)
    z=log(x)
 if(!is.null(lambda)){ nlambda=NULL}
    #stage 1 fit
    if(initial.fit=="lasso"){
        ggg=glmnet(z,y,family=family,intercept=F,standardize=F,alpha=en.alpha,lambda=lambda,nlambda=nlambda,lambda.min.ratio=lambda.min.ratio)
        
   # b=as.vector(coef(ggg,s=lam/n))[-1]
    #    yhat=z%*%b
   }
     
     if(initial.fit=="constrainedLasso"){
         ggg=glmnet.constr(z,y,family=family,alpha=en.alpha,nlambda=nlambda,lambda=lambda,lambda.min.ratio=lambda.min.ratio)
        
            }
 if(verbose==2) cat("initial fit done",fill=TRUE)
    #stage 2 fit

    nz=colSums(ggg$beta!=0)
   
    minv=2
    #remove steps with fewer than 2 variables 
  ggg$beta=ggg$beta[,nz>=minv & nz<=maxvar]
    lambda=ggg$lambda[nz>=minv& nz<=maxvar]
    
     nlam=length(lambda)
     th2=array(NA,c(nlam,maxsteps,p,p)) #first phase model x FS x variable pairs
    th0=array(NA,c(nlam,maxsteps))
    numFS=rep(NA,nlam)
     prss.all=rss0=rep(NA,nlam)
     beta=ggg$beta[,1:nlam]
     #loop over solutions from first stage
    fsout=ind0=vector("list",nlam)

   
    for(ii in 1:nlam){
        
        lam=lambda[ii]
        if(verbose>0) cat(c("lam=",lam),fill=T)
        b=as.vector(ggg$beta[,ii])
        
        yhat=z%*%b
       rss0[ii]=sum( (yhat-mean(y))^2)
       nonzero=which(b!=0)
      ind0[[ii]]=nonzero
      zz=makexx(x[,ind0[[ii]],drop=F])
   
 
     resp=y  #NOTE! in this version, I use y not yhat
    
       if(second.fit=="FS") fsout[[ii]]=myfs.new(zz,resp, family=family,nsteps=maxsteps,verbose=FALSE)
       
        if(second.fit=="constrainedFS") fsout[[ii]]=myfs.constrained(zz,resp,verbose=(verbose>1))
       fsout[[ii]]$ind.orig= matrix(nonzero[fsout[[ii]]$ind],ncol=2)
ggg2=fsout[[ii]]
mod=ggg2$pred
        numFS[ii]=length(mod)
      
   }
   
    return(list(family=family,
                delta=delta,
                initial.fit=initial.fit,
                en.alpha=en.alpha,
                lambda.min.ratio=lambda.min.ratio,
                maxvar=maxvar,
                maxsteps=maxsteps,
                second.fit=second.fit,
               verbose=verbose,
                nlambda=nlambda,
                lambda=lambda,
                meany=meany,glmnet.obj=ggg,fsout=fsout,numFS=numFS,ind0=ind0,call=this.call))
}


predict.logRatioLasso=function(fit,x,ilist=NULL,jlist=NULL){
    #return fitted values from "fit" object returned by logRatioLasso
    # ilist is indices of models in first stage path desired. Default is all models
    # jlist is indices of models in second stage path desired. Default is all models

    # ilist, jlist can either be null or a single value
    #returns linear predictor (for both Gaussian and binomial models)
    
    if(is.null(ilist)) {
                 ilist2=1:length(fit$fsout)
                 ilist.act=ilist2[!is.na(fit$numFS)]
                 
       }
    
     if(!is.null(ilist)) {ilist.act=ilist}
    
     maxterms=max(fit$numFS,na.rm=T)
     yhat=array(NA,c(nrow(x),length(fit$fsout),maxterms))
   

    for(i in ilist.act){
        
         fsfit=fit$fsout[[i]]
        if(!is.null(jlist)) jlist2=min(jlist,nrow(fsfit$ind))
         if(is.null(jlist)) jlist2=1:min(fit$numFS[[i]], length(fit$fsout[[i]]$beta))
        for(j in jlist2){
            pairs=fit$fsout[[i]]$ind[1:j,] #pairs in reduced numbering of FS
            if(!is.matrix(pairs)) pairs=matrix(pairs,nrow=1)
            opairs=matrix(fit$ind0[[i]][pairs],ncol=ncol(pairs)) #pairs in original numbering
             zz= makexx2(x+fit$delta,opairs)
            beta=fit$fsout[[i]]$beta[[j]]
            out=cbind(opairs,beta)
        yhat[,i,j]=fit$fsout[[i]]$b0[j]+zz%*%beta
    }}
    #if making prediction at single grid point, reduce the output matrix
  if(length(ilist.act)==1 & length(jlist==1)) yhat=yhat[,ilist.act,jlist]
return(list(yhat=yhat,model=out))
}



cv.logRatioLasso=function(fit,x,y,nfolds=10,foldid=NULL,keep=F,alpha=1){
    # cv for log ratio lasso
    # "fit" is output from logRatioLasso- original version
    #keep-  should it return prevalidated fits?
    p=ncol(x)
   n=length(y)
     
    nz=colSums(fit$glmnet.obj$beta!=0)

   yhat=array(NA,c(n,length(fit$lambda),max(fit$numFS,na.rm=T)))
      if(is.null(foldid))  foldid = sample(rep(seq(nfolds), length = n))
    nfolds=length(table(foldid))
  
   for(ii in 1:nfolds){
       cat(c("FOLD=",ii),fill=T)

         oo=foldid==ii
         ggg=logRatioLasso(x[!oo,],y[!oo],
             family=fit$family,
                delta=fit$delta,
                initial.fit=fit$initial.fit,
                en.alpha=fit$en.alpha,
                lambda.min.ratio=fit$lambda.min.ratio,
                maxvar=fit$maxvar,
             maxsteps=fit$maxsteps,
                  second.fit=fit$second.fit,
             verbose=fit$verbose,
               lambda=fit$lambda)
                
             junk=predict.logRatioLasso(ggg,x[oo,])$yhat
       nzz=colSums(ggg$glmnet.obj$beta!=0)
       #  mat=match(nzz,nz)
      # o=!is.na(mat)
       #mat=mat[o]
       temp=as.numeric(knn1(matrix(nz,ncol=1),matrix(nzz,ncol=1),1:length(nz)))
       temp=temp[temp < 1+dim(junk)[2]]
      
       kkk=min(dim(junk)[3],dim(yhat)[3])
       if(length(temp)>0){  yhat[oo,1:length(temp),1:kkk]=junk[,temp,1:kkk]}

}

    if(fit$family=="gaussian") errfun=function(yhat,y){  ( (y-yhat)^2) }
  if(fit$family=="binomial") {errfun=function(yhat,y){
      eps=1e-4
     pr=1/(1+exp(-yhat))
     o=( y==1 & (pr>1-eps) ) | (y==0 & (pr<eps))
      val=y*log(pr)+(1-y)*log(1-pr)
   #-2*(val[!o])
      -2*val
   }}

    yhat=myimpute(yhat)  #note: I do a 1NN imputation here in the stage 1 x stage 2 space of models, to fill in the model predictions
            # there's lot of NAs, because the models have diff #of predictors for each fold
    
   ym=array(y,dim(yhat))
    err=errfun(yhat,ym)
   cvm=apply(err,c(2,3),mean,na.rm=T)
    nn=apply(!is.na(err),c(2,3),sum,na.rm=T)
    cvse=sqrt(apply(err,c(2,3),var,na.rm=T)/nn)
    yhat.preval=NULL
    if(keep & fit$family=="gaussian") yhat.preval=yhat
    if(keep & fit$family=="binomial") yhat.preval=1/(1+exp(-yhat))

    ijmin=which(cvm == min(cvm,na.rm=T), arr.ind = TRUE)[1,]
    
    o=which(cvm[ijmin[1],]< cvm[ijmin[1],ijmin[2]]+cvse[ijmin[1],ijmin[2]])[1]
    ijmin.1se=c(ijmin[1],o)
    
return(list(lambda=fit$lambda,cvm=cvm,cvse=cvse,yhat.preval=yhat.preval,ijmin=ijmin,ijmin.1se=ijmin.1se))
}


glmnet.constr=function(x,y,family = c("gaussian", "binomial"),alpha=1,nlambda=50,lambda.min.ratio=.01,lambda=NULL){
 #glmnet with constraint sum beta=0
    #Gaussian case: make sure y is centered!!!!
     family = match.arg(family)
   p=ncol(x)
   n=length(y)
   R=c(rep(1,p),rep(-1,p))
if(family=="gaussian"){
    xx=rbind(cbind(x,-x),R)
    yy=c(y,0)
     w=c(rep(1,n),100*n)
    w=w/sum(w)
    if(!is.null(lambda)) {nlambda=NULL}
     ggg=glmnet(xx,yy,standardize=F,intercept=F,family=family,weights=w,lower.limits=rep(0,ncol(xx)+1),upper.limits=rep(Inf,ncol(xx)+1),alpha=alpha,nlambda=nlambda,lambda=lambda,lambda.min.ratio=lambda.min.ratio)
}

if(family=="binomial"){
    #here we add two fake obs, one at y=0, one at y=1
    xx=rbind(cbind(x,-x),R,R)
    yy=as.integer(c(y,0,1))
     w=c(rep(1,n),n*100,n*100)
    w=w/sum(w)
    # I may want to include an intercept here, but having trouble gettin sum b==0
   # browser()
    ggg=glmnet(xx,yy,standardize=F,intercept=F,family=family,weights=w,lower.limits=rep(0,ncol(xx)+1),upper.limits=rep(Inf,ncol(xx)+1),alpha=alpha,nlambda=nlambda,lambda=lambda,lambda.min.ratio=lambda.min.ratio)
}

   
    b=ggg$beta[1:p,]-ggg$beta[-(1:p),]
return(list(beta=b,a0=ggg$a0,percvar=1-ggg$dev.ratio,lambda=ggg$lambda,glmnet.obj=ggg,family=family))
}

cv.glmnet.constr=function(fit,x,y,nfolds=10,foldid=NULL,keep=F,alpha=1){
    # "fit" is output from glmnet.constr
    p=ncol(x)
   n=length(y)
    family=fit$family
   yhat=matrix(NA,n,length(fit$lambda))
      if(is.null(foldid))  foldid = sample(rep(seq(nfolds), length = n))
    nfolds=length(table(foldid))
  
   for(ii in 1:nfolds){
         oo=foldid==ii
         ggg= glmnet.constr(x[!oo,],y[!oo],family=family,alpha=alpha)
         b=matrix(as.numeric(coef(ggg$glmnet.obj,s=fit$lambda)[-1,]),ncol=length(fit$lambda))
         a0=ggg$a0
         bb=b[(1:p),]-b[-(1:p),]
         yhat[oo,]=a0+x[oo,]%*%bb
}

    if(family=="gaussian") errfun=function(yhat,y){  
       mean( (y-yhat)^2)
   }
  if(family=="binomial") {errfun=function(yhat,y){
      eps=1e-4
     pr=1/(1+exp(-yhat))
     o=( y==1 & (pr>1-eps) ) | (y==0 & (pr<eps))
      val=y*log(pr)+(1-y)*log(1-pr)
   -2*mean(val[!o])
   }}
    
   ym=matrix(y,nrow=n,ncol=ncol(yhat))
   cvm=apply(yhat,2,errfun,y)
    yhat.preval=NULL
    if(keep & family=="binomial") yhat.preval=yhat
    if(keep & family=="binomial") yhat.preval=1/(1+exp(-yhat))
return(list(lambda=fit$lambda,cvm=cvm,yhat.preval=yhat.preval))
}

predict.glmnet.constr=function(fit,x){
    yhat=matrix(fit$a0,nrow=nrow(x),ncol=ncol(fit$b))+x%*%fit$b
    if(fit$family=="binomial") {
        yhat=1/(1+exp(-yhat))
    }
    return(yhat)
}

#NOTE: these next two functions should be written more efficiently
#  to avoid the double loop
    makexx2=
function(x,ind){
     #makes matrix of all pairwise differences for pairs in "ind"
    #ind is a k by 2 matrix of  indices
p=ncol(x)
n=nrow(x)
np=nrow(ind)
xx=matrix(NA,n,np)
ii=0
for(i in 1:np){
ii=ii+1
j=ind[i,1]
k=ind[i,2]
    xx[,ii]=log(x[,j])-log(x[,k])
#cat(c(ii,j,k),fill=T)
}
return(xx)
}

makexx=function(x){
    #makes matrix of all pairwise differences
p=ncol(x)
n=nrow(x)
xx=matrix(NA,n,p*(p-1)/2)
ii=0
for(j in 1:(p-1)){
  for(k in (j+1):p){
    ii=ii+1
    xx[,ii]=log(x[,j])-log(x[,k])
#cat(c(ii,j,k),fill=T)
}}
return(xx)
}



thToMat=function(coef){
    #return upp tri matrix rwith rows  filled with elements of coef
    m=length(coef)
    p=(1+sqrt(1+8*m))/2
#map theta to  matrix form
th=matrix(0,p,p)
    th[row(th)>col(th)]=coef
    return(t(th))
}





myfs=function(x,y,nsteps=min(nrow(x),ncol(x)),center=T, family="gaussian"){
p=ncol(x)
n=length(y)
# fs by minimizing scaled ip;
# for use in logRatioLasso; assumes predictors correspond to all possible pairs
# of some raw variables
# first center x and y
# note at this point, it always uses Gaussian model for the stepwise calculation.
#  but with family="binomial", it returns logistic reg coefs
#returns:
#  pred,s =predictors and their signs, in order entered 
#  scor, bhat: scaled ip and ls coef for each predictor entered
# sigmahat:  est of error variance;
# rss of each model
#pss- %var unexplained by each model
# ind- indices of variable pairs in order entered

y.orig=y
if(center){
x=scale(x,T,F)
y=y-mean(y)
}

nv=.5*(1+sqrt(1+4*2*p));
#construct matrix indices
thmat=matrix(NA,nv,nv)
thmat[row(thmat)>col(thmat)]=1:p
thmat=t(thmat)
yhat=matrix(NA,n,nsteps)
pred=s=scor=bhat=sigmahat=rss=rep(NA,nsteps)
   ip=t(x)%*%y/sqrt(diag(t(x)%*%x))
  pred[1]=which.max(abs(ip))
 s[1]=sign(sum(x[,pred[1]]*y))
scor[1]=ip[pred[1]]
bhat[1]=ip[pred[1]]/sqrt(sum(x[,pred[1]]^2))

r=lsfit(x[,pred[1]],y)$res
rss[1]=sum(r^2)
 sigmahat[1]= sqrt(sum(r^2)/(n-1))
  ind=which(thmat ==pred[1], arr.ind = T)

if(nsteps>1){
for(j in 2:nsteps){
  mod=pred[1:(j-1)]
  r= lsfit(x[,mod],r,int=center)$res
  sigmahat[j]= sqrt(sum(r^2)/n)
  xr= lsfit(x[,mod],x,int=center)$res
  ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
 ip[mod]=0
  pred[j]=which.max(abs(ip))
  scor[j]=ip[pred[j]]
  s[j]=sign(sum(xr[,pred[j]]*r))
 bhat[j]=ip[pred[j]]/sqrt(sum(xr[,pred[j]]^2))
  rss[j]=sum( lsfit(x[,pred[1:j]],y)$res^2)
 
     new=which(thmat ==pred[j], arr.ind = T)
  ind=rbind(ind,new)

}}
#compute ls coefs for all models
beta=vector("list",nsteps)
b0=rep(NA,nsteps)
for(j in 1:nsteps){
    if(family=="gaussian"){
        junk=lsfit(x[,pred[1:j],drop=F],y)
        beta[[j]]=junk$coef[-1]
          b0[j]=junk$coef[j]
}
   #  if(family=="binomial") junk=glm(y.orig~x[,pred[1:j],drop=F],family="binomial")
    if(family=="binomial") {
         junk=glmnet(x[,pred[1:j],drop=F],y.orig,alpha=.05,standardize=FALSE,family="binomial")
        nlam=ncol(junk$beta)
       beta[[j]]=junk$beta[,nlam]
         if(sum(is.na(beta[[j]]))>0) {browser()}
          b0[j]=junk$a0[nlam]
                        }
    
}
prss=rss/sum( (y-mean(y))^2)
return(list(pred=pred,ind=ind,s=s,scor=scor,b0=b0,beta=beta,sigmahat=sigmahat,rss=rss,prss=prss))
}

myfs.new=function(x,y,nsteps=min(nrow(x),ncol(x)),center=T, family="gaussian",verbose=FALSE){
p=ncol(x)
n=length(y)
# fs by minimizing scaled ip;
# for use in logRatioLasso; assumes predictors correspond to all possible pairs
# of some raw variables
# first center x and y
# note at this point, it always uses Gaussian model for the stepwise calculation.
#  but with family="binomial", it returns logistic reg coefs
#returns:
#  pred,s =predictors and their signs, in order entered 
#  scor, bhat: scaled ip and ls coef for each predictor entered
# sigmahat:  est of error variance;
# rss of each model
#pss- %var unexplained by each model
# ind- indices of variable pairs in order entered

y.orig=y
if(center){
x=scale(x,T,F)
y=y-mean(y)
}

out=fs(x, y, maxsteps =nsteps, intercept = FALSE, 
    normalize = FALSE, verbose = verbose)



pred=out$act

nv=.5*(1+sqrt(1+4*2*p));
#construct matrix indices
thmat=matrix(NA,nv,nv)
thmat[row(thmat)>col(thmat)]=1:p
thmat=t(thmat)
ind=matrix(NA,length(pred),2)
for(j in 1:nrow(ind)){
  ind[j,]=which(thmat ==pred[j], arr.ind = T)

}


#compute ls coefs for all models
nsteps.act=min(nsteps,length(pred))
beta=vector("list",nsteps.act)
b0=rep(NA,nsteps.act)
for(j in 1:nsteps.act){

if(family=="gaussian"){
        junk=lsfit(x[,pred[1:j],drop=F],y)
        beta[[j]]=junk$coef[-1]
          b0[j]=junk$coef[j]
}
   #  if(family=="binomial") junk=glm(y.orig~x[,pred[1:j],drop=F],family="binomial")
    if(family=="binomial") {
        if(j==1) {junk=glm(y.orig~x[,pred[1:j],drop=F],family="binomial")
                   beta[[j]]=junk$coef[-1]
                      b0[j]=junk$coef[j]
              }
         if(j>1){
              junk=glmnet(x[,pred[1:j],drop=F],y.orig,alpha=.01,standardize=FALSE,family="binomial")
        nlam=ncol(junk$beta)
       beta[[j]]=junk$beta[,nlam]
        # if(sum(is.na(beta[[j]]))>0) {browser()}
          b0[j]=junk$a0[nlam]
 }}}

return(list(pred=pred,ind=ind,b0=b0,beta=beta))
}


myfs.constrained=function(x,y,nsteps=min(nrow(x),ncol(x)),enter=T,verbose=F,center=T){
p=ncol(x)

# fs by minimizing scaled ip; restricted to a variable appearing in only one pair
# for use in logRatioLasso; assumes predictors correspond to all possible pairs
# of some raw variables
# first center x and y
#returns:
#  pred,s =predictors and their signs, in order entered 
#  scor, bhat: scaled ip and ls coef for each predictor entered
# sigmahat:  est of error variance;
# rss of each model
#pss- %var unexplained by each model
# ind- indices of variable pairs in order entered

if(center){
x=scale(x,T,F)
y=y-mean(y)
}

 nv=.5*(1+sqrt(1+4*2*p));

#construct matrix indices

thmat=matrix(NA,nv,nv)
thmat[row(thmat)>col(thmat)]=1:p
thmat-t(thmat)

 thmat0=thmat
pred=s=scor=bhat=sigmahat=rss=rep(NA,nsteps)
if(verbose) cat("step 1",fill=T)
   ip=t(x)%*%y/sqrt(diag(t(x)%*%x))
  pred[1]=which.max(abs(ip))
 s[1]=sign(sum(x[,pred[1]]*y))
scor[1]=ip[pred[1]]
bhat[1]=ip[pred[1]]/sqrt(sum(x[,pred[1]]^2))

r=lsfit(x[,pred[1]],y)$res
rss[1]=sum(r^2)
 sigmahat[1]= sqrt(sum(r^2)/(n-1))
       rem=which(thmat ==pred[1], arr.ind = T)
       thmat[rem[1],]=thmat[,rem[1]]= thmat[rem[2],]=thmat[,rem[2]]=NA
avail=thmat[!is.na(thmat)]
ind=rem
if(nsteps>1){
    j=1
    while(j<nsteps &  sum(!is.na(thmat>0))) {
        j=j+1
        if(verbose) cat(c("step 1, #avail=", length(avail)),fill=T)
      mod=pred[1:(j-1)]
  
    r= lsfit(x[,mod],r,int=center)$res
  sigmahat[j]= sqrt(sum(r^2)/n)
  xr= lsfit(x[,mod],x,int=center)$res
  ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
    ip[-avail]=0
  pred[j]=which.max(abs(ip))
  scor[j]=ip[pred[j]]
  s[j]=sign(sum(xr[,pred[j]]*r))
    bhat[j]=ip[pred[j]]/sqrt(sum(xr[,pred[j]]^2))
   rem=which(thmat ==pred[j], arr.ind = T)
        ind=rbind(ind,rem)
       thmat[rem[1],]=thmat[,rem[1]]= thmat[rem[2],]=thmat[,rem[2]]=NA
     avail=thmat[!is.na(thmat)]
        rss[j]=sum(lsfit(x[,pred[1:j]],y)$res^2)
 
}}


pred=pred[!is.na(pred)]
s=s[!is.na(s)]
scor=scor[!is.na(scor)]
bhat=bhat[!is.na(bhat)]
sigmahat=sigmahat[!is.na(sigmahat)]
rss=rss[!is.na(rss)]
prss=rss/sum( (y-mean(y))^2)


return(list(pred=pred,ind=ind,s=s,scor=scor,bhat=bhat,sigmahat=sigmahat,rss=rss,prss=prss))
}

upTri=function(x){
    #extract upp tri of matrix in proper row order
    xx=t(x)
    xx[row(xx)>col(xx)]
}
findbestpair=function(z,y){
    ip=t(z)%*%y
    o1=which.max(ip)
    o2=which.max(-ip)
    if(o2<o1){temp=o2;o2=o1;o1=temp}
    return(c(o1,o2))
}
    



myfs2=function(z,y,nsteps=min(nrow(z),ncol(z)),center=T,eps=1e-6){

    # for log ratio model, with z=log(x)


# fs by minimizing **unscaled** ip; (can't do scaled efficiently)

# here x are the  raw variables 
# first center x and y- uses Stephens trick to find best pair each time
#returns:
#  pred,s =predictors and their signs, in order entered
#  scor, bhat: scaled ip and ls coef for each predictor entered
# sigmahat:  est of error variance;
# rss of each model
#pss- %var unexplained by each model
# ind- indices of variable pairs in order entered
p=ncol(z)
n=nrow(z)

if(center){
z=scale(z,T,F)
y=y-mean(y)
}

yhat=zz=matrix(NA,n,nsteps)
pred=matrix(NA,nsteps,2)
scor=bhat=sigmahat=rss=rep(NA,nsteps)
  pred[1,]=findbestpair(z,y)
  zz[,1]=z[,pred[1,1]]- z[,pred[1,2]]
scor[1]=sum(y*zz[,1])
temp=lsfit(zz[,1],y)

r=temp$res
rss[1]=sum(r^2)

 

if(nsteps>1){
    go=T
    j=1
while(j<nsteps  &go){
    j=j+1
go=F
   temp2=findbestpair(z,r)
  zz[,j]=z[,temp2[1]]- z[,temp2[2]]
scor[j]=sum(r*zz[,j])
if(abs(scor[j])>eps){
    go=T
    pred[j,]=temp2
temp=lsfit(zz[,j],r)

 
r=lsfit(zz[,1:j],y)$res
   rss[j]=sum(r^2)
}}}

nstepsact=j-1

zz=zz[,1:nstepsact]

pred=pred[1:nstepsact,]
scor=scor[1:nstepsact]
bhat=bhat[1:nstepsact]
rss=rss[1:nstepsact]
bhat=lsfit(zz,y)$coef
return(list(pred=pred,scor=scor,bhat=bhat,rss=rss,nsteps=nstepsact))
}


predict.myfs2=function(fit,z,center=T){
    nsteps=fit$nsteps
    zz=matrix(NA,nrow(z),fit$steps)
         for(ii in 1:nsteps){
             zz[,ii]=z[,fit$pred[j,1]]-z[,fit$pred[j,2]]
         }
   yhat=cbind(1,zz)%*%fit$bhat
    return(yhat)
}
    


myimpute=function(yhat){
 
d1=dim(yhat)[2]
d2=dim(yhat)[3]
xx=cbind(rep(1:d1,d2),sort(rep(1:d2,d1)))
for(i in 1:dim(yhat)[1]){
   
 yy=as.vector(yhat[i,,])
 cl=1:length(yy)
 o=is.na(yy)
if(sum(o)<length(yy)){
 yc=yy[!o]
 yy[o]=yc[knn1(xx[!o,],xx[o,],cl[!o])]
 yhat[i,,]=matrix(yy,d1,d2)
}
}
return(yhat)
}
