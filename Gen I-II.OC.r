


Gen.OC=function(ncohort1,ncohort2,cohortsize=3,intera=3,ndose=4,toxm,effm,pfsm,toxc,effc,pfsc,ttox,teff,c,beta0s,ints,betaes,betats,gammads,score,a.p,lambda,rho,t1,t2,n2,ntrial){
  
  ## ncohort1: number of cohorts for the dose-finding period in stage I
  ## ncohort2: number of cohorts for the adaptive randomization period in stage I
  ## cohortsize: number of patients in each cohort
  ## intera: inter-arrival time for each cohort of patients  
  ## ndose: number of dose level
  ## toxm: maximumal acceptable toxicity rate
  ## effm: minimal acceptable efficacy rate  
  ## pfsm: minimal acceptable pfs rate
  ## toxc: cut-off posterior probability for toxicity of the admissible set
  ## effc: cut-off posterior probability for efficacy of the admissible set
  ## pfsc: cut-off posterior probability for pfs      of the admissible set
  ## ttox: vector of DLT rates at all the dose levels
  ## teff: matrixs of efficacy rates. First, second and thrid rows represnt the PR/CR, SD and PD rates at all the dose levels
  ## c: co-variance of the latent bi-variate normal to generate ttox and teff outcomes
  ## beta0s: the baseline log-hazards for piecewise exponential survival
  ## ints: the length of each piece
  ## betaes: proportional log-hazard for the efficacy outcome
  ## betats: proportional log-hazard for the toxicity outcome
  ## gammads: proportional log-hazard for dose levels; row for interval; column for dose level
  ## score: desirability score for utility function
  ## a.p: alpha prior for the dirichlet distribution
  ## lambda: tunning parameter adjusting the adaptive randomization rates  
  ## rho: percentage of utility to determine candidate set for stage II
  ## t1: short-term outcome follow-up time
  ## t2: long-term outcome follow-up time
  ## n2: pre-determined sample size of each dose in the candidate set
  ## ntrial: number of simulated trials
  
  
  
  
  
  library(rjags)
  library(mvtnorm)
  library(gtools)
  
  biv=function(ttox,teff,c){
    ## ttox: vector of DLT rates at all the dose levels
    ## ttex: matrixs of efficacy rates. First, second and thrid rows represnt the PR/CR, SD and PD rates at all the dose levels
    at=qnorm(ttox)
    ae=qnorm(teff[1,])
    lambda=ae-qnorm(teff[1,]+teff[2,])
    re=matrix(0,nrow=6,ncol=length(ttox))
    cor=matrix( c(1,c,c,1),nrow=2     )
    for(i in 1: length(ttox)){
      re[1,i]=pmvnorm(  upper=c(0,lambda[i]), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## t=0, e=0 PD
      re[2,i]=pmvnorm( lower=c( -Inf, lambda[i]  ),  upper=c(0,0), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## t=0, e=1 SD
      re[3,i]=pmvnorm( lower=c( -Inf, 0  ),  upper=c(0,Inf), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## t=0, e=2 CR/PR
      
      re[4,i]=pmvnorm( lower=c(0, -Inf),  upper=c(Inf,lambda[i]), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## t=1, e=0 PD
      re[5,i]=pmvnorm( lower=c( 0, lambda[i]  ),  upper=c(Inf,0), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## t=1, e=1 SD
      re[6,i]=pmvnorm( lower=c( 0, 0  ), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## t=0, e=2 CR/PR
      
    }
    return(re)
  }
  
  
  
  
  
  uti=function(ttox,teff,c,score){
    pi=biv(ttox,teff,c)
    re=NULL
    ndose=dim(pi)[2]
    for (i in 1:ndose){
      re[i]=sum(pi[,i]*score)
    }
    return(re)
  }  
  
  
  
  
  
  
  rpwexp <- function(n, rate, intervals=NULL, cumulative=FALSE){
    if(is.null(intervals)){
      if (cumulative){return(cumsum(rexp(n,rate[1])))}else
        return(rexp(n,rate[1]))}
    k <- length(rate)
    if (k==1){
      if(cumulative){return(cumsum(rexp(n,rate)))}else
        return(rexp(n,rate))
    }
    if (length(intervals) < k-1) stop("length(intervals) must be at least length(rate) - 1")
    tx <- 0
    j <- 1
    times <- array(0,n)
    timex <- cumsum(intervals)
    indx <- array(TRUE,n)
    for(i in 1:k){
      nindx <- sum(indx)
      if (nindx==0) break
      increment <- rexp(nindx,rate[i])
      if (cumulative) times[indx] <- tx + cumsum(increment)
      else times[indx] <- tx + increment
      if (i<k){
        tx <- timex[i]
        indx <- (times > timex[i])
      }
    }
    return(times)
  }
  
  
  
  pfs.gen=function(beta0,int,betae,betat,gammad,dose,ye,yt){
    ## gammad: row for interval; column for dose level
    
    ha=c(exp(beta0+betae*ye+betat*yt+gammad[,dose]) ,1)
    return(rpwexp(1,ha,int))
  }
  
  
  
  upfs=function(ttox,teff,c,beta0,int,betae,betat,gammad){
    pi=biv(ttox,teff,c)
    ndose=dim(pi)[2]
    re=NULL
    for(j in 1: ndose){
      re[j]=pi[2,j]*exp(-sum(   exp(beta0+gammad[,j])*int        )       )+pi[3,j]*exp(-sum(   exp(beta0+betae+gammad[,j])*int        )       )+pi[5,j]*exp(-sum(   exp(beta0+betat+gammad[,j])*int        )       )+pi[6,j]*exp(-sum(   exp(beta0+betat+betae+gammad[,j])*int        )       )
      
    }
    return(re)
  }
  
  post=function( data,n,toxm,effm,score,a.p    ){
    ndose=dim(data)[2]
    data.tox=data[4,]+data[5,]+data[6,]
    data.eff=data[3,]+data[6,]
    ut=NULL
    post.tox=NULL
    post.eff=NULL
    for (i in 1: ndose){
      ut[i]=sum( score* ( data[,i]+rep(a.p,6)     )      )/(n[i]+6*a.p)
      post.tox[i]=pbeta(toxm,data.tox[i]+3*a.p, n[i]-data.tox[i]+3*a.p    )
      post.eff[i]=1-pbeta(effm,data.eff[i]+2*a.p, n[i]-data.eff[i]+4*a.p)
    }
    return(list( ut,post.tox,post.eff    ))
    
  }
  
  
  
  
  
  post.pfs=function(data,a.p,beta0,betae,betat,alpha,gamma,t){
    nmcmc=dim(gamma)[1]
    ndose=dim(gamma)[2]
    re=matrix(0,nrow=nmcmc,ncol=ndose)
    for (i in 1:nmcmc){
      
      for (j in 1:ndose){
        
        pi= as.vector(   rdirichlet(1, rep(1/6,6)+data[,j]      )  )
        
        re[i,j]=1-( pi[1]+pi[2]*pweibull(t, alpha[i],  (1/exp(beta0[i]+gamma[i,j])      )  )+pi[3]*pweibull(t, alpha[i],  (1/exp(beta0[i]+betae[i]+gamma[i,j])      )  )+pi[4]+pi[5]*pweibull(t, alpha[i],  (1/exp(beta0[i]+betat[i]+gamma[i,j])      )  )+pi[6]*pweibull(t, alpha[i],  (1/exp(beta0[i]+betat[i]+betae[i]+gamma[i,j])      )  )    )
        
        
      }
      
    }
    
    return(re)
    
    
  }
  
  
  model<-"model{
  
  
  for (i in 1:n.obs){
  t.obs[i]~dweib(alpha,lambda.obs[i])
  lambda.obs[i]<-(exp(betaw0+betawe*eff.obs[i]+betawt*tox.obs[i]+gammawd[dose.obs[i]]     ))^alpha
  }
  
  for (i in 1:n.cen){
  status.cen[i]~dbern(S[i])
  S[i]<-pweib(t.cen[i],alpha,lambda.cen[i])
  lambda.cen[i]<-(exp(betaw0+betawe*eff.cen[i]+betawt*tox.cen[i]+gammawd[dose.cen[i]]       ))^alpha
  }
  
  betaw0~dnorm(0.0,0.01)
  betawe~dnorm(0.0,0.01)
  betawt~dnorm(0.0,0.01)
  alpha~dgamma(0.01,0.01)
  for(i in 1: ndose){
  gammawd[i]~dnorm(0.0,0.01)
  }
  
}"

  
  stageIa=function(ncohort,cohortsize=3,intera=3,ndose=4,toxm,effm,toxc,effc,ttox,teff,c,beta0s,ints=c(1,1),betaes,betats,gammads,score,a.p){
    ## intera: inter-arrival time for each cohort
    ## a.p: alpha prior for the dirichlet distribution
    pi=biv(ttox,teff,c)
    d=1 ## start from the lowest dose
    dur=0 ## trial duration
    datate=NULL
    datate2=matrix(0, nrow=6,ncol=ndose)
    tox=NULL
    sd=NULL 
    eff=NULL
    pfs=NULL
    dose=NULL
    
    for (i in 1:ncohort){
      dur=dur+intera
      dose.temp=rep(d,cohortsize)
      dose=c(dose,dose.temp)
      datate.temp=rmultinom(cohortsize,1,pi[,d])
      datate=cbind(datate,datate.temp)
      datate2[,d]=datate2[,d]+apply(datate.temp,1,sum)
      tox.temp=datate.temp[4,]+datate.temp[5,]+datate.temp[6,]
      tox=c(tox,tox.temp)
      sd.temp=datate.temp[2,]+datate.temp[5,]
      sd=c(sd,sd.temp)
      eff.temp=datate.temp[3,]+datate.temp[6,]
      eff=c(eff,eff.temp)
      pfs.temp=rep(0,cohortsize)
      for(l in 1:cohortsize){
        if(eff.temp[l]==1){pfs.temp[l]=pfs.gen(beta0s,ints,betaes,betats,gammads,d,1,tox.temp[l])}
        else if(sd.temp[l]==1){pfs.temp[l]=pfs.gen(beta0s,ints,betaes,betats,gammads,d,0,tox.temp[l])}
        else {pfs.temp[l]=0}
      }
      pfs=c(pfs,pfs.temp)
      n=apply(datate2,2,sum)
      d.high=max(dose) ## the highest dose has been tried 
      post.est=post( as.matrix(datate2[, 1: d.high]), n[1:d.high],toxm,effm,score,a.p    )
      ut.est=post.est[[1]]
      prob.tox=post.est[[2]]
      prob.eff=post.est[[3]]
      if(length(which(  prob.tox<=toxc     ))==0){
        adm.tox.high=d.high
      } else { adm.tox.high=min(  which(  prob.tox<=toxc   )     )-1          }
      if ( adm.tox.high==0    ){
        adm=0
        break
      } else {
        adm=intersect(   which(prob.eff>effc), 1:adm.tox.high       )
        if (   ( d==d.high   )&( d.high==adm.tox.high    )&( d<ndose  )    ){
          d=d+1
        } else if (    length(adm)==0      ){
          adm=0
          break
        } else {
          d=which(ut.est==max(ut.est[adm]))[1]
        }
        
        
      }  
    }
    return(list( dur,datate,datate2,tox,sd,eff,pfs,dose,adm,n,ut.est[adm]))
    
  }
  
  
  
  
  stageI=function(ncohort1,ncohort2,cohortsize=3,intera=3,ndose=4,toxm,effm,toxc,effc,ttox,teff,c,beta0s,ints=c(1,1),betaes,betats,gammads, score,a.p,lambda,rho){
    pi=biv(ttox,teff,c)
    dataIa=stageIa(ncohort1,cohortsize,intera,ndose,toxm,effm,toxc,effc,ttox,teff,c,beta0s,ints,betaes,betats,gammads,score,a.p)
    dur=dataIa[[1]]
    datate=dataIa[[2]]
    datate2=dataIa[[3]]
    tox=dataIa[[4]]
    sd=dataIa[[5]]
    eff=dataIa[[6]]
    pfs=dataIa[[7]]
    dose=dataIa[[8]]
    adm=dataIa[[9]]
    n=dataIa[[10]]
    post.ut.adm=dataIa[[11]]
    for (i in 1: ncohort2){
      if (sum(adm)==0){
        break
      } else{
        dur=dur+intera
        if (length(post.ut.adm)==1){
          dose.temp=rep(adm,cohortsize)
        }  else {
          ar.prob=post.ut.adm^lambda/sum(post.ut.adm^lambda) 
          dose.temp=sample(adm,cohortsize,replace=T,prob=ar.prob)
        }
        dose=c(dose,dose.temp)
        datate.temp=matrix(0,nrow=6,ncol=cohortsize)
        for(l in 1: cohortsize){
          datate.temp[,l]=rmultinom(1,1,pi[,dose.temp[l]])
          datate2[,dose.temp[l]]=datate2[,dose.temp[l]]+datate.temp[,l]
        }
        datate=cbind(datate,datate.temp)
        tox.temp=datate.temp[4,]+datate.temp[5,]+datate.temp[6,]
        tox=c(tox,tox.temp)
        sd.temp=datate.temp[2,]+datate.temp[5,]
        sd=c(sd,sd.temp)
        eff.temp=datate.temp[3,]+datate.temp[6,]
        eff=c(eff,eff.temp)
        pfs.temp=rep(0,cohortsize)
        for(j in 1:cohortsize){       
          if(eff.temp[j]==1){pfs.temp[j]=pfs.gen(beta0s,ints,betaes,betats,gammads,dose.temp[j],1,tox.temp[j])}
          else if(sd.temp[j]==1){pfs.temp[j]=pfs.gen(beta0s,ints,betaes,betats,gammads,dose.temp[j],0,tox.temp[j])}
          else {pfs.temp[j]=0}        
        }
        pfs=c(pfs,pfs.temp)
        n=apply(datate2,2,sum)
        
        d.high=max(dose) ## the highest dose has been tried 
        post.est=post( as.matrix(datate2[, 1: d.high]), n[1:d.high],toxm,effm,score,a.p    )
        ut.est=post.est[[1]]
        prob.tox=post.est[[2]]
        prob.eff=post.est[[3]]
        if(length(which(  prob.tox<=toxc     ))==0){
          adm.tox.high=d.high
        } else { adm.tox.high=min(  which(  prob.tox<=toxc   )     )-1          }
        if ( adm.tox.high==0    ){
          adm=0
          break
        }  else {
          adm=intersect(   which(prob.eff>effc), 1:adm.tox.high       )
          if (length(adm)==0){
            adm=0
            break
          } else {
            post.ut.adm=ut.est[adm]
          }
          
        }
        
        
      }
    }
    
    if((sum(adm)==0)|(length(adm)==1)    ){
      can=adm
    } else {
      can=intersect(    which( ut.est>=rho*max( post.ut.adm    )[1]             )      ,   adm )
    } 
    return(list( dur,datate,datate2,tox,sd,eff,pfs,dose,n,can))
  }
  
  
  
  Gen=function(ncohort1,ncohort2,cohortsize=3,intera=3,ndose=4,toxm,effm,toxc,effc,ttox,teff,c,beta0s,ints,betaes,betats,gammads,score,a.p,lambda,rho,t1,t2,n2,pfsm,pfsc){ 
    
#     ut=uti(ttox,teff,c,score) ## true utility 
#     uts=upfs(ttox,teff,c,beta0s,ints,betaes,betats,gammads) ## true utility for PFS
    pi=biv(ttox,teff,c)
    dataI=stageI(ncohort1,ncohort2,cohortsize,intera,ndose,toxm,effm,toxc,effc,ttox,teff,c,beta0s,ints,betaes,betats,gammads,score,a.p,lambda,rho)
    dur=dataI[[1]]
    datate=dataI[[2]]
    datate2=dataI[[3]]
    tox=dataI[[4]]
    sd=dataI[[5]]
    eff=dataI[[6]]
    pfs=dataI[[7]]
    dose=dataI[[8]]
    n=dataI[[9]]
    ss.stageI=sum(n)
    can=dataI[[10]]
    pfs[pfs>(t2-t1)]=(t2-t1)
    if (  sum(can)==0  ){
      dur=dur-intera+t2
      obd=0
    }else {
      n.stageII=n2-n[can]
      n.stageII[n.stageII<0]=0
      dur=dur+   (    ceiling(sum(n.stageII)/cohortsize)        -1)*intera+t2
      if (sum(n.stageII)==0){
        dose=dose
        datate=datate
        datate2=datate2
        tox=tox
        sd=sd
        eff=eff
        n=n
        pfs=pfs
      } else {
        
        dose.temp=rep(  can, n.stageII      )
        dose=c(dose,dose.temp)
        datate.temp=matrix(0,nrow=6,ncol=sum(n.stageII))
        for(l in 1: sum(n.stageII)){
          datate.temp[,l]=rmultinom(1,1,pi[,dose.temp[l]])
          datate2[,dose.temp[l]]=datate2[,dose.temp[l]]+datate.temp[,l]
        }
        datate=cbind(datate,datate.temp)
        tox.temp=datate.temp[4,]+datate.temp[5,]+datate.temp[6,]
        tox=c(tox,tox.temp)
        sd.temp=datate.temp[2,]+datate.temp[5,]
        sd=c(sd,sd.temp)
        eff.temp=datate.temp[3,]+datate.temp[6,]
        eff=c(eff,eff.temp)
        n=apply(datate2,2,sum)
        pfs.temp=rep(0,sum(n.stageII))
        for(j in 1:sum(n.stageII)){
          
          if (  eff.temp[j]==1 ){
            pfs.temp[j]=pfs.gen(beta0s,ints,betaes,betats,gammads,dose.temp[j],1,tox.temp[j])
          } else if (sd.temp[j]==1){
            pfs.temp[j]=pfs.gen(beta0s,ints,betaes,betats,gammads,dose.temp[j],0,tox.temp[j])
          } else {pfs.temp[j]=0}
          
        }
        pfs=c(pfs,pfs.temp) 
        pfs[pfs>(t2-t1)]=(t2-t1)        
      }
      status=rep(1,length(pfs))
      status[pfs==(t2-t1)]=0
      data.final=cbind(pfs[sd==1|eff==1],status[sd==1|eff==1],dose[sd==1|eff==1],tox[sd==1|eff==1],eff[sd==1|eff==1])
      data.obs<-subset(data.final,status[sd==1|eff==1]==1)
      data.cen<-subset(data.final,status[sd==1|eff==1]==0)
      n.obs<-length(data.obs[,1])
      n.cen<-length(data.cen[,1])
      data.mcmc<-list(n.obs=n.obs,n.cen=n.cen,t.obs=data.obs[,1],t.cen=data.cen[,1],
                      dose.obs=data.obs[,3], dose.cen=data.cen[,3],tox.obs=data.obs[,4], tox.cen=data.cen[,4],
                      eff.obs=data.obs[,5], eff.cen=data.cen[,5],status.cen=data.cen[,2],ndose=ndose)
      mcmc.final<- jags.model(textConnection(model), data=data.mcmc,  quiet=T)
      
      
      update(mcmc.final, 1000, progress.bar="none")
      post<- jags.samples(mcmc.final, n.iter=1000,variable.names=c("betaw0","betawe","betawt","alpha","gammawd"), progress.bar="none")
      post.betaw0.sample=post$betaw0[1:    length(post$betaw0)   ]
      post.betawe.sample=post$betawe[1:    length(post$betawe)   ]
      post.betawt.sample=post$betawt[1:    length(post$betawt)   ]
      post.alpha.sample=post$alpha[1:    length(post$alpha)   ]
      post.gammawd.sample=matrix(post$gammawd,ncol=ndose,byrow=T)
      d.high=max(dose) ## the highest dose has been tried 
      post.est=post( as.matrix(datate2[, 1: d.high]), n[1:d.high],toxm,effm,score,a.p    )
      prob.tox=post.est[[2]]
      post.pfs.sample=post.pfs( datate2, a.p, post.betaw0.sample, post.betawe.sample, post.betawt.sample, post.alpha.sample, post.gammawd.sample, t2-t1     )
      post.pfs.est=apply(post.pfs.sample, 2,mean)
      prob.pfs=apply( post.pfs.sample>pfsm, 2, mean      )
      if(length(which(  prob.tox<=toxc     ))==0){
        adm.tox.high=d.high
      } else { adm.tox.high=min(  which(  prob.tox<=toxc   )     )-1          }
      
      if ( adm.tox.high==0    ){
        obd=0
        
      }  else{
        
        adm= intersect( intersect(   which(prob.pfs>pfsc), 1:adm.tox.high       ), can)
        
        if (length(adm)==0){
          obd=0
        } else{
          obd=which(post.pfs.est==max(post.pfs.est[adm]))[1]
        }
        
      }
      
      
    } 
    
    return(list(dur,sum(tox)/sum(n),sum(eff)/sum(n),  sum(sd)/sum(n) ,n,can,obd,ss.stageI))
  }
  
  
  
  
  ut=uti(ttox,teff,c,score) ## true utility 
  uts=upfs(ttox,teff,c,beta0s,ints,betaes,betats,gammads) ## true utility for PFS
  
  
  
  t.adm=which(  (ttox<=toxm)&(uts>=pfsm)    )
  if (  length(t.adm)==0    ){
    t.obd=0
  } else {    t.obd=which(  uts==max(    uts[t.adm] )      )   }
  dur.trial=NULL
  toxp.trial=NULL
  effp.trial=NULL
  sdp.trial=NULL
  n.trial=matrix(0,nrow=ntrial,ncol=ndose)
  can.trial=NULL
  obd.trial=NULL
  ss.stageI.trial=NULL
  ss.stageII.trial=NULL
  ss.trial=NULL 
  
for (i in 1: ntrial){
  data=Gen(ncohort1,ncohort2,cohortsize=3,intera=3,ndose=4,toxm,effm,toxc,effc,ttox,teff,c,beta0s,ints,betaes,betats,gammads,score,a.p,lambda,rho,t1,t2,n2,pfsm,pfsc) 
  dur.trial=c( dur.trial ,     data[[1]])
  toxp.trial=c(toxp.trial, data[[2]])
  effp.trial=c(effp.trial,    data[[3]])
  sdp.trial=c(sdp.trial,    data[[4]])
  n.trial[i,]=data[[5]]
  can.trial=c(can.trial,data[[6]])
  obd.trial=c(obd.trial,data[[7]])
  ss.stageI.trial=c(ss.stageI.trial, data[[8]])
  ss.stageII.trial=c(ss.stageII.trial, sum(data[[5]])-data[[8]]   )
  ss.trial=c(ss.trial,sum(data[[5]]))
} 

  
if ( t.obd==0   ){
  R=NA
} else {
  obd2.trial=obd.trial[obd.trial!=0]
  temp=0
  for (j in 1: length(obd2.trial   )){
    temp=temp+uts[  obd2.trial[j]   ]/length(obd2.trial   )
  }
  R=temp/uts[t.obd]
}
  
  
return( list("true toxicity rate"=ttox, "true tumor response rate"=teff[1,], "true utility"=ut, "true PFS utility"=uts,"true optimal dose"=t.obd,
             "dose-finding percentage"=table(obd.trial)/ntrial, "patients allocation"=apply(n.trial,2,mean), "trial duration"=mean(dur.trial), "candidate doses percentage"=table(can.trial)/ntrial, 
             "sample size for stage I"=mean(ss.stageI.trial), "sample size for stage II"=mean(ss.stageII.trial), "sample size for the entire trial"=mean(ss.trial), "R"=R       )    )
 
  
} 


date()
# Sce 1.

set.seed(123456)
Gen.OC(ncohort1=5,ncohort2=9,cohortsize=3,intera=3,ndose=4,toxm=0.3,effm=0.5,pfsm=0.4,toxc=0.1,effc=0.1,pfsc=0.1,ttox=c(0.1,0.2,0.4,0.5), teff=matrix(   c( 0.3,0.4,0.5,0.55,0.6,0.55,0.45,0.4,0.1,0.05,0.05,0.05       ), byrow=T, nrow=3    ),c=0.2,beta0s=c(log(-log(0.3)/2)-0.1, log(-log(0.3)/2)+0.1     ),
       ints=c(1,1),betaes=rep(  -0.5    ,2)         ,betats=rep( 0.1     ,2),  gammads=matrix(   c( log( log(0.029)/log(0.3)    )    ,log( log(0.029)/log(0.3)    ),log( log(0.085)/log(0.35)    )    ,log( log(0.085)/log(0.35)    ),log( log(0.135)/log(0.36)    )    ,log( log(0.135)/log(0.36)    ),log( log(0.298)/log(0.37)    )    ,log( log(0.298)/log(0.37)    )    )      , nrow=2  ),score=c(20,50,100,0,30,60),a.p=1/6,lambda=0.5,rho=0.7,t1=1,
       t2=3,n2=21,ntrial=5000)



## Sce 2.



set.seed(123456)
Gen.OC(ncohort1=5,ncohort2=9,cohortsize=3,intera=3,ndose=4,toxm=0.3,effm=0.5,pfsm=0.4,toxc=0.1,effc=0.1,pfsc=0.1,ttox=c(0.05,0.2,0.4,0.5), teff=matrix(   c(0.6,0.65,0.55,0.4,0.35,0.1,0.1,0.1,0.05,0.25,0.35,0.5), byrow=T, nrow=3    ),c=0.2,beta0s=c(log(-log(0.3)/2)-0.1, log(-log(0.3)/2)+0.1     ),
       ints=c(1,1),betaes=rep(  -0.5    ,2)         ,betats=rep( 0.1     ,2),
       gammads=matrix(   c(log(  log(0.085)/log(0.39) ),log(  log(0.085)/log(0.39) ),log(  log(0.112)/log(0.34) ),log(  log(0.112)/log(0.34) ),log(  log(0.103)/log(0.28) ),log(  log(0.103)/log(0.28) ),log(  log(0.110)/log(0.21) ),log(  log(0.110)/log(0.21) )   )      , nrow=2  )
       ,score=c(20,50,100,0,30,60),a.p=1/6,lambda=0.5,rho=0.7,t1=1,
       t2=3,n2=21,ntrial=5000)



## Sce 3.



set.seed(123456)
Gen.OC(ncohort1=5,ncohort2=9,cohortsize=3,intera=3,ndose=4,toxm=0.3,effm=0.5,pfsm=0.4,toxc=0.1,effc=0.1,pfsc=0.1,ttox=c(0.15,0.2,0.25,0.4), teff=matrix(   c( 0.6,0.7,0.6,0.5,0.3,0.3,0.3,0.3,0.1,0,0.1,0.2), byrow=T, nrow=3    ),c=0.2,beta0s=c(log(-log(0.3)/2)-0.1, log(-log(0.3)/2)+0.1     ),
       ints=c(1,1),betaes=rep(  -0.5    ,2)         ,betats=rep( 0.1     ,2),
       gammads=matrix(   c(log(   log(0.46)/log(0.37)  ), log(  log(0.46)/log(0.37)     )  ,log(  log(0.605)/log(0.42)  ),log(  log(0.605)/log(0.42) ),log( log(0.463)/log(0.37)    ),log( log(0.463)/log(0.37)   ),log(  log(0.43)/log(0.32)  ),log( log(0.43)/log(0.32)    )   )      , nrow=2  )
       ,score=c(20,50,100,0,30,60),a.p=1/6,lambda=0.5,rho=0.7,t1=1,
       t2=3,n2=21,ntrial=5000)



## Sce 4.


set.seed(123456)
Gen.OC(ncohort1=5,ncohort2=9,cohortsize=3,intera=3,ndose=4,toxm=0.3,effm=0.5,pfsm=0.4,toxc=0.1,effc=0.1,pfsc=0.1,ttox=c(0.01,0.03,0.06,0.08), teff=matrix(   c(0.6,0.7,0.7,0.7,0.15,0.2,0.25,0.25,0.25,0.1,0.05,0.05), byrow=T, nrow=3    ),c=0.2,beta0s=c(log(-log(0.3)/2)-0.1, log(-log(0.3)/2)+0.1     ),
       ints=c(1,1),betaes=rep(  -0.5    ,2)         ,betats=rep( 0.1     ,2),
       gammads=matrix(   c(log(log(0.425)/log(0.33)    ),log(  log(0.425)/log(0.33) ),log( log(0.455)/log(0.39)  ),log(   log(0.455)/log(0.39) ),log(  log(0.675)/log(0.41) ) ,log(  log(0.675)/log(0.41) ),log(  log(0.51)/log(0.41) ),log(  log(0.51)/log(0.41) )   )      , nrow=2  )
       ,score=c(20,50,100,0,30,60),a.p=1/6,lambda=0.5,rho=0.7,t1=1,
       t2=3,n2=21,ntrial=5000)



## Sce 5.


set.seed(123456)
Gen.OC(ncohort1=5,ncohort2=9,cohortsize=3,intera=3,ndose=4,toxm=0.3,effm=0.5,pfsm=0.4,toxc=0.1,effc=0.1,pfsc=0.1,ttox=c(0.03,0.05,0.1,0.15), teff=matrix(   c(0.4,0.6,0.75,0.65,0.4,0.3,0.2,0.2,0.2,0.1,0.05,0.15), byrow=T, nrow=3    ),c=0.2,beta0s=c(log(-log(0.3)/2)-0.1, log(-log(0.3)/2)+0.1     ),
       ints=c(1,1),betaes=rep(  -0.5    ,2)         ,betats=rep( 0.1     ,2),
       gammads=matrix(   c(log( log(0.297)/log(0.31)   ),log( log(0.297)/log(0.31)   ),log( log(0.466)/log(0.38)   ),log( log(0.466)/log(0.38)   ),log( log(0.51)/log(0.42)   ),log( log(0.51)/log(0.42)   ),log( log(0.735)/log(0.37)   ),log( log(0.735)/log(0.37)   )   )      , nrow=2  )
       ,score=c(20,50,100,0,30,60),a.p=1/6,lambda=0.5,rho=0.7,t1=1,
       t2=3,n2=21,ntrial=5000)



## sce 6.

 set.seed(123456)
 Gen.OC(ncohort1=5,ncohort2=9,cohortsize=3,intera=3,ndose=4,toxm=0.3,effm=0.5,pfsm=0.4,toxc=0.1,effc=0.1,pfsc=0.1,ttox=c(0.03,0.05,0.1,0.15), teff=matrix(   c(0.4,0.5,0.6,0.5,0.5,0.4,0.3,0.45,0.1,0.1,0.1,0.05), byrow=T, nrow=3    ),c=0.2,beta0s=c(log(-log(0.3)/2)-0.1, log(-log(0.3)/2)+0.1     ),
         ints=c(1,1),betaes=rep(  -0.5    ,2),betats=rep( 0.1     ,2),  
        gammads=matrix(   c(0.12,0.12,-0.294,-0.294,-1.01,-1.01,-0.095,-0.095)      , nrow=2  )   ,score=c(20,50,100,0,30,60),a.p=1/6,lambda=0.5,rho=0.7,t1=1,
         t2=3,n2=21,ntrial=5000)




## sce 7.

 set.seed(123456)
 Gen.OC(ncohort1=5,ncohort2=9,cohortsize=3,intera=3,ndose=4,toxm=0.3,effm=0.5,pfsm=0.4,toxc=0.1,effc=0.1,pfsc=0.1,ttox=c(0.04,0.06,0.08,0.1), teff=matrix(   c(0.4,0.5,0.6,0.6,0.35,0.3,0.3,0.35,0.25,0.2,0.1,0.05), byrow=T, nrow=3    ),c=0.2,beta0s=c(log(-log(0.3)/2)-0.1, log(-log(0.3)/2)+0.1     ),
        ints=c(1,1),betaes=rep(  -0.5    ,2),betats=rep( 0.1     ,2),  
        gammads=matrix(   c(0.36,0.36,-0.26,-0.26,-0.41,-0.41,-1.095,-1.095)      , nrow=2  )   ,score=c(20,50,100,0,30,60),a.p=1/6,lambda=0.5,rho=0.7,t1=1,
        t2=3,n2=21,ntrial=5000)


## sce 8.

set.seed(123456)
Gen.OC(ncohort1=5,ncohort2=9,cohortsize=3,intera=3,ndose=4,toxm=0.3,effm=0.5,pfsm=0.4,toxc=0.1,effc=0.1,pfsc=0.1,ttox=c(0.02,0.04,0.07,0.1), teff=matrix(   c(0.4,0.6,0.65,0.75,0.1,0.2,0.1,0.05,0.5,0.2,0.25,0.2), byrow=T, nrow=3    ),c=0.2,beta0s=c(log(-log(0.3)/2)-0.1, log(-log(0.3)/2)+0.1     ),
       ints=c(1,1),betaes=rep(  -0.5    ,2),betats=rep( 0.1     ,2),  
       gammads=matrix(   c(-2.059,-2.059,-1.86,-1.86,-0.68,-0.68,-0.289,-0.289)      , nrow=2  )   ,score=c(20,50,100,0,30,60),a.p=1/6,lambda=0.5,rho=0.7,t1=1,
       t2=3,n2=21,ntrial=5000)




date()








  
  