findbestcut=function(y,yhat,patid){
    # find best cutp for patient classification
    #: y, yhat are vectors of pixels
    cutplist=seq(.1,.9,length=100)
    ypat=unique(patid)
    npat=length(unique(patid))
        for(i in 1:npat){
            o=which(patid==i)[1]
            ypat[i]=y[o]
        }
 yhatc.pat=matrix(NA,length(cutplist),npat)
        for(k in 1:length(cutplist)){
           yhatc=1*(yhat>cutplist[k])
           for(kk in 1:npat){
                yhatc.pat[k,kk]=1*(mean(yhatc[patid==kk])>.5)
            }}

err=rowMeans(matrix(ypat,nrow=length(cutplist),ncol=npat,byrow=T)-1 != yhatc.pat)
    ihat=which.min(err)
  return(cutplist[ihat])
}
 logratio.QP.path=function(x,y,lamlist=NULL,nlam=NULL){
     if(is.null(lamlist)){
     lammax=max(abs(t(x)%*%y))
     if(is.null(nlam)) nlam=50
      lamlist=seq(lammax,lammax/100,length=nlam) # should instead be equally spaced on log scale
 }
        b=matrix(NA,ncol(x),length(lamlist))
    for(k in 1:length(lamlist)){
        cat("lambda=",lamlist[k],fill=T)
        b[,k]=logratio.QP(x,y,lamlist[k])
    }
        return(list(b=b,lamlist=lamlist))
    }

logratio.QP=function(x,y,gam){
    #do lasso on x with constraint sum b=0
 p=ncol(x)
 SMALL=1e-7
eps=.0001 #needed- o/w complains that dmat is singular
xx=cbind(x,-x)
dmat=t(xx)%*%xx+eps*diag(2*p)
dvec=t(xx)%*%y-gam*rep(1,2*p)

Amat=c(rep(1,p),rep(-1,p))
Amat=rbind(Amat,diag(2*p))
bvec=rep(0,2*p+1)

a=solve.QP(dmat,dvec,t(Amat),bvec,meq=1)
bbb=a$sol[1:p]-a$sol[-(1:p)]
bbb[abs(bbb)<SMALL]=0
return(bbb)
}

 soft=
function(x,lam){sign(x)*(abs(x)-lam)*(abs(x)>lam)}

 psoft=
function(x,lam){(x-lam)*(x>lam)}

soft2=function(beta,tt,w){
tt=rep(tt,length(beta))
val=c(beta[1],soft(beta[-1],tt[-1]/w[-1]))
return(val)
}
psoft2=function(beta,tt,w){
tt=rep(tt,length(beta))
val=c(beta[1],psoft(beta[-1],tt[-1]/w[-1]))
return(val)
}


critf=function(beta,x,y,lam,s=0,sum.to.zero=F){
    p=ncol(x)/2
yhat=x%*%beta
    extra=0
    
    if(sum.to.zero) extra=s*sum(beta[1:p]-beta[-(1:p)])
val=.5*sum( (y-yhat)^2)+lam*sum(abs(beta))+extra
return(val)
}



critf.logistic=function(beta,x,y,lam,s=0,sum.to.zero=F){
    #note- assumes that  beta includes an  intercept
       p=(ncol(x)-1)/2
       beta2=beta[-1]
eta=x%*%beta
pr=as.vector(1/(1+exp(-eta)))
         extra=0
    if(sum.to.zero) extra=s*sum(beta2[1:p]-beta2[-(1:p)])
val=-sum(y*log(pr)+(1-y)*log(1-pr)) +lam*sum(abs(beta2))+extra
   
return(val)
}

ggrad.path=function(x,y,nlambda=50,tt=.01,eps=.001,maxiter=500,trace=F,backtrack.flag=F,sum.to.zero=F){
     #do lasso on x w/o or with constraint sum b=0, over a path am lambda values
    #- eps= conv threshold for all loops
    # sum.to.zero=  enforce constraint?
    # maxiter= max #iterations in outer and inner loops
    #backtrack.flag- backtrack in gen grad?
    #Output: beta- matrix of solutions
    #   lambda=lam values used
    #  ierr.inner, ierr- error flags- 1 means not converged
    #  niter= number of outer iterations used
 lammax=max(abs(t(x)%*%y))
  ierr=ierr.inner=niter=rep(NA,nlambda)
 
 lamlist=seq(lammax,lammax/100,length=nlambda) # should instead be equally spaced on log scale
 betaall=matrix(NA,ncol(x),nlambda)
 b=rep(0,ncol(x))
 beta=rep(0,2*length(b))
  slow=-50;shi=50
 s=0
 for(j in 1:nlambda){
     cat(c("lambda=",lamlist[j]),fill=T)
 #   for(j in 1:15){
     #if(s==0) { slow=-50;shi=50}
    # if(s!=0) {slow=s/5;shi=5*s}
     if(sum.to.zero)   out=ggrad.constr(x,y,lamlist[j],tt,beta=beta,slow=slow,shi=shi,maxiter=maxiter,s=s,eps=eps,trace=trace,backtrack.flag=backtrack.flag)
     if(!sum.to.zero)   out=ggrad(x,y,lamlist[j],b,tt,maxiter=maxiter,eps=eps,trace=trace,backtrack.flag=backtrack.flag)
     betaall[,j]=out$b
     b=betaall[,j]
    if(sum.to.zero)  s=out$s
   #  beta=c(b*(b>0),-b*(b<0))
   if(sum.to.zero)  ierr[j]=out$ierr
      ierr.inner[j]=out$ierr.inner
     niter[j]=out$niter
 }
 return(list(beta=betaall,lambda=lamlist,ierr.inner=ierr.inner,ierr=ierr,niter=niter))
}

ggrad.constr=function(x,y,lam,tt,beta=rep(0,2*ncol(x)),maxiter=500,s=0,slow=-50,shi=50,maxiter.outer=500,eps=.01,backtrack.flag=F,trace=F){
    p=ncol(x)
xx=cbind(x,-x)


a0=ggrad(xx,y,lam,beta,tt=tt,maxiter=maxiter,sum.to.zero=T,s=slow,trace=trace)
b0=a0$beta[1:p]-a0$beta[-(1:p)]

a1=ggrad(xx,y,lam,beta,tt=tt,maxiter=maxiter,sum.to.zero=T,s=shi,trace=trace)
b1=a1$beta[1:p]-a1$beta[-(1:p)]
if(sum(b0)<0| sum(b1)>0) stop("initial interval not wide enough")

ii=0
go=TRUE
while(ii < maxiter.outer & go){
    ii=ii+1
go=FALSE
    
a2=ggrad(xx,y,lam,beta,tt=tt,maxiter=maxiter,sum.to.zero=T,eps=eps,backtrack.flag=F,s=s,trace=trace)
beta=a2$b
b=a2$b[1:p]-a2$b[-(1:p)]
    ierr.inner=a2$ierr.inner
   

#cat(c(s,sum(b)),fill=T)
if(sum(b)>0){
    slow=s
    s=(slow+shi)/2
}
if(sum(b)<0){
    shi=s
    s=(slow+shi)/2
}
   go=abs(sum(b))>eps
 if(trace)   cat(c(slow,shi,sum(b)),fill=T)
}
    ierr=0
if(ii==maxiter.outer){
    cat("outer loop failed to converge",fill=T);ierr=1
}
return(list(b=b,s=s,slow=slow,shi=shi,ierr.inner=ierr.inner,ierr=ierr,niter=ii))
}


ggrad=function(x,y,lam,beta,maxiter=200,tt=.01,eps=.01,trace=T,backtrack.flag=F,sum.to.zero=F,s=0){
 #beta is starting parameter vector
    #betat is true pop'n value (for simulation study)

       #note- does not fit an intercept!
rat=rep(NA,maxiter)

R=rep(0,length(beta))
shrinker=soft
if(sum.to.zero){
    R=c(rep(1,length(beta)/2),rep(-1,length(beta)/2))
     shrinker=psoft
 }
crit=10e7
ii=0
go=TRUE
while(ii<maxiter & go){
    go=FALSE
 ii=ii+1
 g=beta+tt*(t(x)%*%(y-x%*%beta)-s*R)
 beta.old=beta
 beta=shrinker(g,tt*lam)
 
if(backtrack.flag){
 beta=backtrack(beta.old,lam,x,y,s=s,tt=tt)$beta
}
critold=crit
crit=critf(beta,x,y,lam,s,sum.to.zero)
    del=abs(critold-crit)
    go=del>eps
if(trace){cat(crit,fill=T)}

}
ierr.inner=0
if(ii==maxiter) {cat("inner loop failed to converge",fill=T);ierr.inner=1
   cat(c(maxiter,eps,del),fill=T)
             }
return(list(b=beta,niter=ii,ierr.inner=ierr.inner))
}

#NOT YET IMPLEMENTED for some to zreo case
nest=function(x,y,beta,betat,lam,niter=100,tt=.01,gam=.8,trace=T,traceb=T,backtrack.flag=F){
# nesterovs momentum method for lasso
#  beta is starting vector, betat is final solution from QP, for reporting purposes
rat=rep(NA,niter)
critt=critf(betat,x,y,lam)
theta=beta
for(ii in 1:niter){
 g=-t(x)%*%(y-x%*%theta)
 beta.old=beta
 r=theta-tt*g
 beta=soft(r,tt*lam)
if(backtrack.flag){
 tt=backtrack.nest(theta,lam,x,y,tt=tt,gam=gam,trace=traceb)
 beta=soft(r,tt*lam)
}
theta=beta+(ii/(ii+3))*(beta-beta.old)
crit=critf(beta,x,y,lam)
if(trace){cat(crit,fill=T)}
rat[ii]=(crit-critt)/critt
}
return(list(beta=beta,rat=rat,niter=ii))
}


backtrack=function(b,lam,x,y,tt=1,gam=.8,s=0,sum.to.zero=F,trace=F){
# for non-smooth functions
g0=critf(b,x,y,0,s,sum.to.zero)
dg=-t(x)%*%(y-x%*%b)
go=T
while(go){
 go=F
 r=b-tt*dg
 G=(1/tt)*(b-soft(r,tt*lam))
 val=g0-tt*sum(dg*G)+(tt/2)*sum(G^2)
 beta1=b-tt*G
 g1=critf(beta1,x,y,0)
if(g1>val){go=T;tt=gam*tt}
if(trace){cat(c(tt,val,g1),fill=T)}
}
return(list(beta=beta1,tt=tt))
}

#NOT YET IMPLEMENTED for some to zreo case
backtrack.nest=function(theta,lam,x,y,tt=1,gam=.8,trace=F){
g0=critf(theta,x,y,0)
g=-t(x)%*%(y-x%*%theta)
go=T
while(go){
 go=F
 r=theta-tt*g
 G=(1/tt)*(theta-soft(r,tt*lam))
 val=g0-tt*sum(g*G)+(tt/2)*sum(G^2)
 theta1=theta-tt*G
 g1=critf(theta1,x,y,0)
if(g1>val){go=T;tt=gam*tt}
if(trace){cat(c(tt,val,g1),fill=T)}
}
return(tt)
}


ggrad.logistic=function(x,y,beta,lam,maxiter=100,tt=.01,eps=.0001,trace=T,backtrack.flag=F,hess=F,sum.to.zero=F,s=0){
    #note- fits an intercept
xx=cbind(1,x)
R=rep(0,length(beta))
shrinker=soft2
if(sum.to.zero){
    pp=(length(beta)-1)/2
    R=c(0,rep(1,pp),rep(-1,pp))
     shrinker=psoft2
 }
crit=10e7


go=TRUE
while(ii<maxiter & go){
    go=FALSE
 ii=ii+1
eta=xx%*%beta
pr=as.vector(1/(1+exp(-eta)))
w=rep(1,ncol(xx))
if(hess){w=diag(t(xx)%*%diag(pr*(1-pr))%*%xx)}
    temp=t(xx)%*%(y-pr)-s*R
g=beta+tt*as.vector(temp)/w
beta.old=beta
beta=soft2(g,tt*lam,w)
if(backtrack.flag){
beta=backtrack.logistic(beta.old,xx,y,lam,tt=tt,w=w,s=s,sum.to.zero=sum.to.zero)$beta

}
    critold=crit
crit=critf.logistic(beta,xx,y,lam,s,sum.to.zero)
      del=abs(critold-crit)
    go=del>eps
if(trace){cat(crit,fill=T)}
}
return(list(beta=beta,niter=ii))
}

genx=function(n=20, rr=0, p=8){

#    generate x's multivariate normal
        inds <- 1:p
        Sigma <- rr^abs(outer(inds, inds, "-")) #covariance of the x's
        bb <- svd(Sigma)
        hh <- bb$u %*% (sqrt(diag(bb$d))) %*% t(bb$v)
        x <- matrix(rnorm(p * n), n, p) %*% hh
        return(list(x=x, Sigma=Sigma))
}



 
backtrack.logistic=function(b,x,y,lam,tt=1,alpha=.5,gam=.8,s=0,w=rep(1,length(b))){
g0=critf.logistic(b,x,y,0,s)
eta=x%*%b
pr=as.vector(1/(1+exp(-eta)))
dg=-t(x)%*%(y-pr)/w
go=T
while(go){
go=F
r=as.vector(b-tt*dg)
G=(1/tt)*(b-soft2(r,tt*lam,w))
val=g0-tt*sum(dg*G)+(tt/2)*sum(G^2)
beta1=b-tt*G
g1=critf.logistic(beta1,x,y,0)
if(g1>val){go=T;tt=gam*tt}
cat(c(tt,val,g1),fill=T)
}
return(list(beta=beta1,tt=tt))
}




ggrad.logistic.constr=function(x,y,lam,tt,maxiter=100,slow=-1,shi=1,maxiter.outer=50,eps=.001){
xx=cbind(x,-x)
beta=rep(0,ncol(xx)+1)
p=ncol(x)
s=0
a0=ggrad.logistic(xx,y,beta,lam,tt=tt,maxiter=maxiter,sum.to.zero=T,s=slow)
beta2=beta[-1]
b0=a0$beta2[1:p]-a0$beta2[-(1:p)]

a1=ggrad.logistic(xx,y,beta,lam,tt=tt,maxiter=maxiter,sum.to.zero=T,s=shi)
beta2=beta[-1]
b1=a1$beta2[1:p]-a1$beta2[-(1:p)]
if(sum(b0<0) | sum(b1)>0) stop("initial interval not wide enough")

ii=0
go=TRUE
while(ii < maxiter.outer & go){
    ii=ii+1
go=FALSE
    
a2=ggrad.logistic(xx,y,beta,lam,tt=tt,maxiter=maxiter,sum.to.zero=T,s=s)
beta=a2$beta
    beta2=beta[-1]
b=a2$beta2[1:p]-a2$beta2[-(1:p)]

cat(c(s,sum(b)),fill=T)
if(sum(b)>0){
    slow=s
    s=(slow+shi)/2
}
if(sum(b)<0){
    shi=s
    s=(slow+shi)/2
}
   go=abs(sum(b))>eps
}
if(ii==maxiter.outer) cat("outer loop failed to converge",fill=T)
return(b)
}


betaToTheta=function(beta){
    p=length(beta)
    theta=theta2=matrix(NA,p,p)
    ord=order(-beta)
    bb=beta[ord]
    pplus=sum(beta>=0)
    if(sum(beta!=0)>0){
    for(i in 1:pplus){
        for(j in (pplus+1):p){
            theta2[i,j]=2*abs(bb[i])*abs(bb[j])/sum(abs(bb))
        }}
   theta[ord,ord]=theta2
  #  theta=theta2[ord,ord]
}
    theta[is.na(theta)]=0
    return(theta)
}

 thetaToBeta=function(theta){
p=nrow(theta)
s1=s2=rep(0,p)
for(k in 1:p){
if(k>1) s1[k]=s1[k]- sum(theta[1:(k-1), k])
if(k<p) s2[k]=s2[k]+sum(theta[k,(k+1):p])
}
bhat=s1+s2
    return(bhat)
}
makexx=function(x){
p=ncol(x)
n=nrow(x)
xx=matrix(NA,n,p*(p-1)/2)
ii=0
for(j in 1:(p-1)){
  for(k in (j+1):p){
ii=ii+1
    xx[,ii]=log(x[,j])-log(x[,k])

}}
return(xx)
}

thToB=function(coef,p){
#map theta to beta
    #first fill in theta matrix from vector solution
th=matrix(0,p,p)
ii=0
for(j in 1:(p-1)){
  for(k in (j+1):p){
   ii=ii+1
 th[j,k]=coef[ii]
}}
#now transform to betas
s1=s2=rep(0,p)
for(k in 1:p){
if(k>1) s1[k]=s1[k]- sum(th[1:(k-1), k])
if(k<p) s2[k]=s2[k]+sum(th[k,(k+1):p])
}
bhat=s1+s2
    return(list(b=bhat,th=th))
}


thToMat=function(coef){
    m=length(coef)
    p=(1+sqrt(1+8*m))/2
#map theta to  matrix form
th=matrix(0,p,p)
ii=0
for(j in 1:(p-1)){
  for(k in (j+1):p){
   ii=ii+1
 th[j,k]=coef[ii]
}}

    return(th)
}
fit1=function(x,y,lam,nvar.bss=NULL,initial.fit="lasso",en.alpha=1,second.fit="lasso",lam.second=NULL,en.alpha.constr=.1,noisesd=.001, second.resp="y"){
    #fit  logs at first stage and then log ratios at second stage
    # initial.fit can be lasso, or constrained lasso (sum to 0), both with paramaters lam  and  en par en.alpha
    #second fit can be lasso, sparsenet (par lam.second), or bestsub (#terms =nvar.bss)
    # second.resp = response at second stage, either "y" or "yhat"
    n=length(y)
    p=ncol(x)
    z=log(x)
    if(initial.fit=="lasso"){
        ggg=glmnet(z,y,intercept=F,standardize=F,alpha=en.alpha)
    b=as.vector(coef(ggg,s=lam/n))[-1]
        yhat=z%*%b
   }
     if(initial.fit=="constrained.lasso"){
    z3=rbind(z,sqrt( (1-en.alpha.constr)*lam)*diag(p))
    yy=c(y,rep(0,p))
    b=logratio.QP(z3,yy,en.alpha.constr*(lam/2))
    yhat=(z3%*%b)[1:n]
     }
    
    ind0=which(b!=0)
    zz=makexx(x[,ind0])
    zzn=zz+sqrt(en.alpha.constr)*matrix(rnorm(n*ncol(zz)),n,ncol(zz))
    rss0=NULL
    rss=NULL
    norm=NULL
    if(second.fit=="lasso"){
    ggg=glmnet(zzn,resp,intercept=F,standardize=F)
    th=as.numeric(coef(ggg,s=lam.second/n))[-1]
  
    th2=matrix(0,p,p)
        th2[ind0,ind0]=thToMat(th)
        }
    resp=y
    if(second.resp=="yhat") resp=yhat
    if(second.fit=="sparsenet"){
        zzc=scale(zz,T,F)
        ggg2=sparsenet(zzc,resp)
        th=as.vector(coef(ggg2,s=.01*lam.second/n,which.gamma=5))[-1]
        yhat=zz%*%th
        rss=sum( (y-yhat)^2)
        th2=matrix(0,p,p)
        th2[ind0,ind0]=thToMat(th)
         }
    if(second.fit=="bestsub"){
        rss=norm=rep(NA,nvar.bss)
        a=regsubsets(zzn,resp,nbest=nvar.bss,intercept=F)
        aa=summary(a)
        o=which(as.numeric(rownames(aa$which))==nvar.bss)
        
        th=matrix(0,length(o),ncol(zz))
        for(j in 1:length(o)){
            ind=which(aa$which[o[j],])
        aaa=lsfit(zzn[,ind],yhat,int=F)$coef
        th[j,ind]=aaa
            junk=lsfit(zz[,ind],y,int=F)
        rss[j]=sum(junk$res^2)
            norm[j]=sum(abs(junk$coef))
        }
    
     th2=array(0,c(p,p,nvar.bss))
    for(j in 1:nvar.bss){
     th2[ind0,ind0,j]=thToMat(th[j,])
     }
    }

        rss0=sum( (y-mean(y))^2)
    return(list(b=b,th=th2,rss0=rss0,rss=rss, norm=norm))
}


fit2=function(x,y,lam,pfrat=1,method="lasso",noisesd=.001){
    #fit log main+log ratio 
    n=length(y)
    p=ncol(x)
    z=log(x)
    gg=glmnet(z,y,intercept=F,standardize=F)
   
    a=as.numeric(coef(gg,s=lam/n))[-1]
    ind=which(a!=0)
    yhat=y-lsfit(z[,ind],y)$res
    zz=makexx(x[,ind])
    zz=zz+noisesd*matrix(rnorm(nrow(zz)*ncol(zz)),nrow(zz),ncol(zz))
    zall=cbind(z[,ind],zz)
    pf=rep(1,ncol(zall))
    pf[-(1:length(ind))]=pfrat
    if(method=="lasso"){
    ggg=glmnet(zall,yhat,intercept=F,penalty.factor=pf,standardize=F)
    aa=as.numeric(coef(ggg,s=lam/n))[-1]
       }
    if(method=="sparsenet"){
        ggg2=sparsenet(zall,yhat)
        aa=coef(ggg2,s=.1*lam/n,which.gamma=5)
    }
    b=aa[1:length(ind)]
    bb=rep(0,p)
    bb[ind]=b
    th=aa2[-(1:length(ind))]
    th2=matrix(0,p,p)
    ii=0
    for(j in ind){
        ind2=ind[ind>j]
        for(k in ind2){
            ii=ii+1
            th2[j,k]=th[ii]
        }}
    return(list(b=bb,th=th2,ind=ind))
}


glmnet.constr=function(x,y,family="gaussian",alpha=1){
 #glmnet with constraint sum beta=0
    #Gaussian case: make sure y is centered!!!!
p=ncol(x)
n=length(y)
R=c(rep(1,p),rep(-1,p))
if(family=="gaussian"){
    xx=rbind(cbind(x,-x),R)
    yy=c(y,0)
     w=c(rep(1,n),100*n)
    w=w/sum(w)
     ggg=glmnet(xx,yy,standardize=F,intercept=F,family="gaussian",weights=w,lower.limits=rep(0,ncol(xx)+1),upper.limits=rep(Inf,ncol(xx)+1),alpha=alpha)
}

if(family=="binomial"){
    #here we add two fake obs, one at y=0, one at y=1
    xx=rbind(cbind(x,-x),R,R)
    yy=c(y,0,1)
     w=c(rep(1,n),n*100,n*100)
    w=w/sum(w)
    # I want to include an intercept here, but having trouble gettin sum b==0
    ggg=glmnet(xx,yy,standardize=F,intercept=F,family="binomial",weights=w,lower.limits=rep(0,ncol(xx)+1),upper.limits=rep(Inf,ncol(xx)+1),alpha=alpha)
}

   
    b=ggg$beta[1:p,]-ggg$beta[-(1:p),]
return(list(b=b,a0=ggg$a0,lambda=ggg$lambda,glmnet.obj=ggg,family=family))
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
    
