#Shrinkage estimators in zero-inflated Bell regression model with application
#This code is associated with the article referenced as arXiv:2403.00749. 
#Please ensure to cite it if you utilize this code.

rm(list = ls())

# Recall the functions
source("bellloglikfun.R")

#-----------------------experiment design----------------------#
m=1000     # the number of replications
p = 7      # the number of count paraneters 
q = 3      # the number of zero parameters 
k = p + q  # the number of total parameters in the model
nobs= 200   # the number of observations

# The true value of parameters
betav = c(0.5, 1,-1.5, rep(0, times=p-3))
zetav= c(0.5, -1, 0)

# Restrictions on the model
r<- p-2  # the number of restrictions
H<- matrix(0, nrow=r, ncol=(p+q))

H[1, q]<- 1
if (r>1){for(i in 1:(r-1)){ H[i+1, 6+i]<- 1}}else{H}


# Define matrixes for storing estimators
# When delta=0
BMLE<- matrix(nrow = m, ncol = k)    # For UNZIBE
BRE<- matrix(nrow = m, ncol = k)     # For REZIBE
BJS<- matrix(nrow = m, ncol = k)     # For JSZIBE
BPJS<- matrix(nrow = m, ncol = k)    # For PJSZIBE
BPT<- matrix(nrow = m, ncol = k)     # For PTZIBE

# When delta=0.2
BRE2<- matrix(nrow = m, ncol = k)    # For REZIBE
BJS2<- matrix(nrow = m, ncol = k)    # For JSZIBE
BPJS2<- matrix(nrow = m, ncol = k)   # For PJSZIBE
BPT2<- matrix(nrow = m, ncol = k)    # For PTZIBE

# When delta=0.4
BRE4<- matrix(nrow = m, ncol = k)    # For REZIBE
BJS4<- matrix(nrow = m, ncol = k)    # For JSZIBE 
BPJS4<- matrix(nrow = m, ncol = k)   # For PJSZIBE
BPT4<- matrix(nrow = m, ncol = k)    # For PTZIBE

# When delta=0.6
BRE6<- matrix(nrow = m, ncol = k)    # For UNZIBE
BJS6<- matrix(nrow = m, ncol = k)    # For JSZIBE
BPJS6<- matrix(nrow = m, ncol = k)   # For PJSZIBE
BPT6<- matrix(nrow = m, ncol = k)    # For PTZIBE

# When delta=0.8
BRE8<- matrix(nrow = m, ncol = k)    # For UNZIBE
BJS8<- matrix(nrow = m, ncol = k)    # For JSZIBE
BPJS8<- matrix(nrow = m, ncol = k)   # For PJSZIBE
BPT8<- matrix(nrow = m, ncol = k)    # For PTZIBE

# When delta=1.0
BRE1<- matrix(nrow = m, ncol = k)    # For UNZIBE
BJS1<- matrix(nrow = m, ncol = k)    # For JSZIBE
BPJS1<- matrix(nrow = m, ncol = k)   # For PJSZIBE
BPT1<- matrix(nrow = m, ncol = k)    # For PTZIBE

# When delta=1.2
BRE12<- matrix(nrow = m, ncol = k)   # For UNZIBE
BJS12<- matrix(nrow = m, ncol = k)   # For JSZIBE
BPJS12<- matrix(nrow = m, ncol = k)  # For PJSZIBE
BPT12<- matrix(nrow = m, ncol = k)   # For PTZIBE

# When delta=1.4
BRE14<- matrix(nrow = m, ncol = k)   # For UNZIBE
BJS14<- matrix(nrow = m, ncol = k)   # For JSZIBE
BPJS14<- matrix(nrow = m, ncol = k)  # For PJSZIBE
BPT14<- matrix(nrow = m, ncol = k)   # For PTZIBE

# When delta=1.6
BRE16<- matrix(nrow = m, ncol = k)   # For UNZIBE
BJS16<- matrix(nrow = m, ncol = k)   # For JSZIBE
BPJS16<- matrix(nrow = m, ncol = k)  # For PJSZIBE
BPT16<- matrix(nrow = m, ncol = k)   # For PTZIBE

# When delta=1.8
BRE18<- matrix(nrow = m, ncol = k)   # For UNZIBE
BJS18<- matrix(nrow = m, ncol = k)   # For JSZIBE
BPJS18<- matrix(nrow = m, ncol = k)  # For PJSZIBE
BPT18<- matrix(nrow = m, ncol = k)   # For PTZIBE

# When delta=2.0
BRE20<- matrix(nrow = m, ncol = k)   # For UNZIBE
BJS20<- matrix(nrow = m, ncol = k)   # For JSZIBE
BPJS20<- matrix(nrow = m, ncol = k)  # For PJSZIBE
BPT20<- matrix(nrow = m, ncol = k)   # For PTZIBE

# loop for simulation
for (j in 1:m) {
  
  #generate covariates observations for mean
  X<- cbind(rep(1, times=nobs), matrix(rnorm((p-1)*nobs), nrow = nobs, ncol = (p-1)))
  Xn<- X[,-1]
  
  #generate covariates observations for zero part   
  S<- cbind(rep(1, times=nobs),matrix(rexp((q-1)*nobs,1), nrow = nobs, ncol = (q-1)))  
  Sn<- S[,-1]
  
  #calculate mean and pi
  meant<- exp(X%*%betav)
  pit<- exp(S%*%zetav)/(1+exp(S%*%zetav))
  
  # calculate parameter of bell parameter
  library(LambertW)
  theta<- W(meant)
  
  # generate observation of zibell distribution
  library(bellreg)
  
  y<- c()
  for (i in 1:nobs) {
    y[i]<- rzibell(theta[i], pit[i])
  }
  
  # data frame 
  ndata<- data.frame(y, Xn, Sn)
  library(bellreg)
  
  # derive UNZIBE
  md1 <- zibellreg(y~ Sn|Xn, data=ndata, approach = "mle")
  BMLE[j,]<- c(coef(md1)[[1]], coef(md1)[[2]])
  
  # obtain information matrix
  Vb<- vcov(md1)    
  
  # inverse of the information matrix
  XWX<- ginv(Vb,tol = sqrt(.Machine$double.eps))
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=0
  
  # derive REZIBE  
  h<- rep(0, times=r)
  DTV<- H%*%BMLE[j,]-h
  D1<- ginv(H%*%Vb%*%t(H),tol = sqrt(.Machine$double.eps))
  BRE[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV
  
  # Calculate test statistic
  Tstat<- t(DTV)%*%D1%*%DTV
  df<- r
  achit<- qchisq(0.975, df)
  
  # derive JSZIBE
  WJS<- 1-((r-2)/as.numeric(Tstat))
  BJS[j,]<- BRE[j,]+WJS*(BMLE[j,]-BRE[j,])
  
  # derive PJSZIBE
  ZZ<- ifelse(WJS>0, WJS, 0)
  BPJS[j,]<- BRE[j,]+ZZ*(BMLE[j,]-BRE[j,])
  
  # derive PTZIBE
  BPT[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat<= achit, 1, 0))*(BMLE[j,]-BRE[j,])
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=0.2
  
  # derive REZIBE
  h2<- c(sqrt(0.2),rep(0, times=r-1))
  DTV2<- H%*%BMLE[j,]-h2
  BRE2[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV2
  
  # calculate test statistic
  Tstat2<- t(DTV2)%*%D1%*%DTV2
  
  # derive JSZIBE
  WJS2<- 1-((r-2)/as.numeric(Tstat2))
  BJS2[j,]<- BRE2[j,]+WJS2*(BMLE[j,]-BRE2[j,])
  
  # derive PJSZIBE
  ZZ2<- ifelse(WJS2>0, WJS2, 0)
  BPJS2[j,]<- BRE2[j,]+ZZ2*(BMLE[j,]-BRE2[j,])
  
  # derive PTZIBE
  BPT2[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat2<= achit, 1, 0))*(BMLE[j,]-BRE2[j,]) 
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=0.4
  
  # derive REZIBE
  h4<- c(sqrt(0.4),rep(0, times=r-1))
  DTV4<- H%*%BMLE[j,]-h4
  BRE4[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV4
  
  # calculate test statistic
  Tstat4<- t(DTV4)%*%D1%*%DTV4
  
  # derive JSZIBE
  WJS4<- 1-((r-2)/as.numeric(Tstat4))
  BJS4[j,]<- BRE4[j,]+WJS4*(BMLE[j,]-BRE4[j,])
  
  # derive PJSZIBE
  ZZ4<- ifelse(WJS4>0, WJS4, 0)
  BPJS4[j,]<- BRE4[j,]+ZZ4*(BMLE[j,]-BRE4[j,])
  
  # derive PTZIBE
  BPT4[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat4<= achit, 1, 0))*(BMLE[j,]-BRE4[j,]) 
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=0.6
  
  # derive REZIBE
  h6<- c(sqrt(0.6),rep(0, times=r-1))
  DTV6<- H%*%BMLE[j,]-h6
  BRE6[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV6
  
  # calculate test statistic
  Tstat6<- t(DTV6)%*%D1%*%DTV6
  
  # derive JSZIBE
  WJS6<- 1-((r-2)/as.numeric(Tstat6))
  BJS6[j,]<- BRE6[j,]+WJS6*(BMLE[j,]-BRE6[j,])
  
  # derive PJSZIBE
  ZZ6<- ifelse(WJS6>0, WJS6, 0)
  BPJS6[j,]<- BRE6[j,]+ZZ6*(BMLE[j,]-BRE6[j,])
  
  # derive PTZIBE
  BPT6[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat6<= achit, 1, 0))*(BMLE[j,]-BRE6[j,]) 
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=0.8
  
  # derive REZIBE
  h8<- c(sqrt(0.8),rep(0, times=r-1))
  DTV8<- H%*%BMLE[j,]-h8
  BRE8[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV8
  
  # calculate test statistic
  Tstat8<- t(DTV8)%*%D1%*%DTV8
  
  # derive JSZIBE
  WJS8<- 1-((r-2)/as.numeric(Tstat8))
  BJS8[j,]<- BRE8[j,]+WJS8*(BMLE[j,]-BRE8[j,])
  
  # derive PJSZIBE
  ZZ8<- ifelse(WJS8>0, WJS8, 0)
  BPJS8[j,]<- BRE8[j,]+ZZ8*(BMLE[j,]-BRE8[j,])
  
  # derive PTZIBE
  BPT8[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat8<= achit, 1, 0))*(BMLE[j,]-BRE8[j,]) 
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=1.0
  
  # derive REZIBE
  h1<- c(1,rep(0, times=r-1))
  DTV1<- H%*%BMLE[j,]-h1
  BRE1[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV1
  
  # calculate test statistic
  Tstat1<- t(DTV1)%*%D1%*%DTV1
  
  # derive JSZIBE
  WJS1<- 1-((r-2)/as.numeric(Tstat1))
  BJS1[j,]<- BRE1[j,]+WJS1*(BMLE[j,]-BRE1[j,])
  
  # derive PJSZIBE
  ZZ1<- ifelse(WJS1>0, WJS1, 0)
  BPJS1[j,]<- BRE1[j,]+ZZ1*(BMLE[j,]-BRE1[j,])
  
  # derive PTZIBE
  BPT1[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat1<= achit, 1, 0))*(BMLE[j,]-BRE1[j,])
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=1.2
  
  # derive REZIBE
  h12<- c(sqrt(1.2),rep(0, times=r-1))
  DTV12<- H%*%BMLE[j,]-h12
  BRE12[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV12
  
  # calculate test statistic
  Tstat12<- t(DTV12)%*%D1%*%DTV12
  
  # derive JSZIBE
  WJS12<- 1-((r-2)/as.numeric(Tstat12))
  BJS12[j,]<- BRE12[j,]+WJS12*(BMLE[j,]-BRE12[j,])
  
  # derive PJSZIBE
  ZZ12<- ifelse(WJS12>0, WJS12, 0)
  BPJS12[j,]<- BRE12[j,]+ZZ12*(BMLE[j,]-BRE12[j,])
  
  # derive PTZIBE
  BPT12[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat12<= achit, 1, 0))*(BMLE[j,]-BRE12[j,])
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=1.4
  
  # derive REZIBE
  h14<- c(sqrt(1.4),rep(0, times=r-1))
  DTV14<- H%*%BMLE[j,]-h14
  BRE14[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV14
  
  # calculate test statistic
  Tstat14<- t(DTV14)%*%D1%*%DTV14
  
  # derive JSZIBE
  WJS14<- 1-((r-2)/as.numeric(Tstat14))
  BJS14[j,]<- BRE14[j,]+WJS14*(BMLE[j,]-BRE14[j,])
  
  # derive PJSZIBE
  ZZ14<- ifelse(WJS14>0, WJS14, 0)
  BPJS14[j,]<- BRE14[j,]+ZZ14*(BMLE[j,]-BRE14[j,])
  
  # derive PTZIBE
  BPT14[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat14<= achit, 1, 0))*(BMLE[j,]-BRE14[j,])
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=1.6
  
  # derive REZIBE
  h16<- c(sqrt(1.6),rep(0, times=r-1))
  DTV16<- H%*%BMLE[j,]-h16
  BRE16[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV16
  
  # calculate test statistic
  Tstat16<- t(DTV16)%*%D1%*%DTV16
  
  # derive JSZIBE
  WJS16<- 1-((r-2)/as.numeric(Tstat16))
  BJS16[j,]<- BRE16[j,]+WJS16*(BMLE[j,]-BRE16[j,])
  
  # derive PJSZIBE
  ZZ16<- ifelse(WJS16>0, WJS16, 0)
  BPJS16[j,]<- BRE16[j,]+ZZ16*(BMLE[j,]-BRE16[j,])
  
  # derive PTZIBE
  BPT16[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat16<= achit, 1, 0))*(BMLE[j,]-BRE16[j,])
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=1.8
  
  # derive REZIBE
  h18<- c(sqrt(1.8),rep(0, times=r-1))
  DTV18<- H%*%BMLE[j,]-h18
  BRE18[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV18
  
  # calculate test statistic
  Tstat18<- t(DTV18)%*%D1%*%DTV18
  
  # derive JSZIBE
  WJS18<- 1-((r-2)/as.numeric(Tstat18))
  BJS18[j,]<- BRE18[j,]+WJS18*(BMLE[j,]-BRE18[j,])
  
  # derive PJSZIBE
  ZZ18<- ifelse(WJS18>0, WJS18, 0)
  BPJS18[j,]<- BRE18[j,]+ZZ18*(BMLE[j,]-BRE18[j,])
  
  # derive PTZIBE
  BPT18[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat18<= achit, 1, 0))*(BMLE[j,]-BRE18[j,])
  
  #---------------------------------------------------------------------#
  # Caculate estimators when delta=2.0
  
  # derive REZIBE
  h20<- c(sqrt(2),rep(0, times=r-1))
  DTV20<- H%*%BMLE[j,]-h20
  BRE20[j,]<- BMLE[j,]- Vb%*%t(H)%*%D1%*%DTV20
  
  # calculate test statistic
  Tstat20<- t(DTV20)%*%D1%*%DTV20
  
  # derive JSZIBE
  WJS20<- 1-((r-2)/as.numeric(Tstat20))
  BJS20[j,]<- BRE20[j,]+WJS20*(BMLE[j,]-BRE20[j,])
  
  # derive PJSZIBE
  ZZ20<- ifelse(WJS20>0, WJS20, 0)
  BPJS20[j,]<- BRE20[j,]+ZZ20*(BMLE[j,]-BRE20[j,])
  
  # derive PTZIBE
  BPT20[j,]<-BMLE[j,]-as.numeric(ifelse(Tstat20<= achit, 1, 0))*(BMLE[j,]-BRE20[j,])
  
}

# true value of the parameters
b<- betav

# The output of the simulation
MSE_MLE<- MSE(BMLE[,-c(1:3)],b)
MSE_BRE<- MSE(BRE[,-c(1:3)],b)
MSE_BJS<- MSE(BJS[,-c(1:3)],b)
MSE_BPJS<- MSE(BPJS[,-c(1:3)],b)
MSE_BPT<- MSE(BPT[,-c(1:3)],b)

MSE_BRE2<- MSE(BRE2[,-c(1:3)],b)
MSE_BJS2<- MSE(BJS2[,-c(1:3)],b)
MSE_BPJS2<- MSE(BPJS2[,-c(1:3)],b)
MSE_BPT2<- MSE(BPT2[,-c(1:3)],b)

MSE_BRE4<- MSE(BRE4[,-c(1:3)],b)
MSE_BJS4<- MSE(BJS4[,-c(1:3)],b)
MSE_BPJS4<- MSE(BPJS4[,-c(1:3)],b)
MSE_BPT4<- MSE(BPT4[,-c(1:3)],b)

MSE_BRE6<- MSE(BRE6[,-c(1:3)],b)
MSE_BJS6<- MSE(BJS6[,-c(1:3)],b)
MSE_BPJS6<- MSE(BPJS6[,-c(1:3)],b)
MSE_BPT6<- MSE(BPT6[,-c(1:3)],b)

MSE_BRE8<- MSE(BRE8[,-c(1:3)],b)
MSE_BJS8<- MSE(BJS8[,-c(1:3)],b)
MSE_BPJS8<- MSE(BPJS8[,-c(1:3)],b)
MSE_BPT8<- MSE(BPT8[,-c(1:3)],b)

MSE_BRE1<- MSE(BRE1[,-c(1:3)],b)
MSE_BJS1<- MSE(BJS1[,-c(1:3)],b)
MSE_BPJS1<- MSE(BPJS1[,-c(1:3)],b)
MSE_BPT1<- MSE(BPT1[,-c(1:3)],b)

MSE_BRE12<- MSE(BRE12[,-c(1:3)],b)
MSE_BJS12<- MSE(BJS12[,-c(1:3)],b)
MSE_BPJS12<- MSE(BPJS12[,-c(1:3)],b)
MSE_BPT12<- MSE(BPT12[,-c(1:3)],b)

MSE_BRE14<- MSE(BRE14[,-c(1:3)],b)
MSE_BJS14<- MSE(BJS14[,-c(1:3)],b)
MSE_BPJS14<- MSE(BPJS14[,-c(1:3)],b)
MSE_BPT14<- MSE(BPT14[,-c(1:3)],b)

MSE_BRE16<- MSE(BRE16[,-c(1:3)],b)
MSE_BJS16<- MSE(BJS16[,-c(1:3)],b)
MSE_BPJS16<- MSE(BPJS16[,-c(1:3)],b)
MSE_BPT16<- MSE(BPT16[,-c(1:3)],b)

MSE_BRE18<- MSE(BRE18[,-c(1:3)],b)
MSE_BJS18<- MSE(BJS18[,-c(1:3)],b)
MSE_BPJS18<- MSE(BPJS18[,-c(1:3)],b)
MSE_BPT18<- MSE(BPT18[,-c(1:3)],b)

MSE_BRE20<- MSE(BRE20[,-c(1:3)],b)
MSE_BJS20<- MSE(BJS20[,-c(1:3)],b)
MSE_BPJS20<- MSE(BPJS20[,-c(1:3)],b)
MSE_BPT20<- MSE(BPT20[,-c(1:3)],b)


# Derive the relative efficiency 
RE<- MSE_MLE/c(MSE_BRE,MSE_BJS,MSE_BPJS,MSE_BPT)
RE2<- MSE_MLE/c(MSE_BRE2,MSE_BJS2,MSE_BPJS2,MSE_BPT2)
RE4<- MSE_MLE/c(MSE_BRE4,MSE_BJS4,MSE_BPJS4,MSE_BPT4)
RE6<- MSE_MLE/c(MSE_BRE6,MSE_BJS6,MSE_BPJS6,MSE_BPT6)
RE8<- MSE_MLE/c(MSE_BRE8,MSE_BJS8,MSE_BPJS8,MSE_BPT8)
RE1<- MSE_MLE/c(MSE_BRE1,MSE_BJS1,MSE_BPJS1,MSE_BPT1)
RE12<- MSE_MLE/c(MSE_BRE12,MSE_BJS12,MSE_BPJS12,MSE_BPT12)
RE14<- MSE_MLE/c(MSE_BRE14,MSE_BJS14,MSE_BPJS14,MSE_BPT14)
RE16<- MSE_MLE/c(MSE_BRE16,MSE_BJS16,MSE_BPJS16,MSE_BPT16)
RE18<- MSE_MLE/c(MSE_BRE18,MSE_BJS18,MSE_BPJS18,MSE_BPT18)
RE20<- MSE_MLE/c(MSE_BRE20,MSE_BJS20,MSE_BPJS20,MSE_BPT20)

Total_RE<- matrix(ncol = 5, nrow = 11)
Total_RE[,1]<- seq(0,2,0.2)
colnames(Total_RE)<- c("delta","RE","JS", "PJS", "PT")
Total_RE[1,-1]<- RE
Total_RE[2,-1]<- RE2
Total_RE[3,-1]<- RE4
Total_RE[4,-1]<- RE6
Total_RE[5,-1]<- RE8
Total_RE[6,-1]<- RE1
Total_RE[7,-1]<- RE12
Total_RE[8,-1]<- RE14
Total_RE[9,-1]<- RE16
Total_RE[10,-1]<- RE18
Total_RE[11,-1]<- RE20

print(Total_RE)

