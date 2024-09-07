
# number generation for Zero-inflated bell distribution
rzibell<- function(theta, zerop){
           library(bellreg)
           ifelse(rbinom(1, 1, zerop) == 0, rbell(1, theta), 0)
          }


#Log likelihood function of Zero-inflated bell model
loglik<- function(theta){
  library(bellreg)
  library(LambertW) 
  beta1<- theta[1:p]
  gam1<- theta[(p+1):k]
  eta1<- X%*%beta1
  eta2<- S%*%gam1
  
  mu1<- exp(eta1)
  pi<- exp(eta2)/(1+exp(eta2))
  
  loglik0<- 0
  loglik1<- 0
  for (i in 1:length(y)) {
    if(y[i]==0)
    {loglik0<- loglik0+log(pi[i]+(1-pi[i])*dbell(0, W(mu1[i])))}
    else
    {loglik1<- loglik1+ log(1-pi[i])+(1-exp(W(mu1[i])))+dbell(y[i], W(mu1[i]), log=T)}
  }
  
  loglik0+loglik1
}


# MSE of estimators
MSE<- function(R,TR){
  p<- length(TR)
  
  ABS<- matrix(nrow = nrow(R), ncol = ncol(R))
  
  for (i in 1:p) {
    ABS[,i]<- (R[,i]-TR[i])^2
  }
  
  msec<- apply(ABS, 2, mean)
  tmse<- sum(msec)
  tmse
}




