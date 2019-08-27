library(e1071)
library(boot)
library(spTDyn)
library(spNNGP)
library(mvtnorm)
library(truncnorm)
library(ROCR)
library(polspline)
library(spatstat)
library(MCMCpack) # inverse wishart
library(MASS) # also for ordered logit / probit
library(parallel)
library(Matrix)
library(sparseMVN)
library(fields)
library(RcppTN)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(tmvtnorm)
library(class) #for knn
library(nnet) #for multinomial logistic regression
library(randomForest)
# library(mgcv) # GAM multinomial logistic regression
library(spatstat)
library(data.table) # for rbindlist: stack fillin into DATA
library(doSNOW)
library(parallel)
S=100

n.cores=17
index=1:34
resolution4=4
resolution9=9
Z=2 # number of different segmentations considered
source("~/Google Drive/mpMRI/17.1210/basic_functions.R")
source("sl_functions.R")

base_functions = c("rf","knn","qda")
base_function = "glm" #the chosen base function 
knn.k=250 # giving the smallest cv.error.rate
n.tree = 136
n.X = 3 # number of X in SL
slcovnames = c("Global_pred","Region4_pred","Region9_pred")
covariates = c("Global_pred","Region4_pred","Region9_pred")
### Define MRIpars:
MRIpars = c("ADC","KTRANS","KEP","AUGC")
rf.baseformula = as.formula(paste("cancer ~ ", paste(MRIpars, collapse= "+")))


rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}


########################################################################
########################################################################
#### Using complete data to get region specific MRI distributions: #####
load("~/Google Drive/mpMRI/server/mpmri/msi_fillrecenter.RData")
fillin=fillrecenter
rm(fillrecenter)
index=1:46
for (i in 1:46)
{
  fillin[[i]][,"ADC"]=fillin[[i]][,"ADC"]/100
  fillin[[i]][,c("KTRANS","KEP","AUGC")]=log(fillin[[i]][,c("KTRANS","KEP","AUGC")])
}
########## Generate data for jags: 
data=cbind(as.data.frame(fillin[[1]][,c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone")]),1)
colnames(data)=c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone","subject")
for (k in c(index[-1]))
{
  temp=cbind(as.data.frame(fillin[[k]][,c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone")]),k)
  colnames(temp)=c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone","subject")
  data=rbind(data,temp)
}
dataNC=data[data[,"cancer"]==0,]
dataC=data[data[,"cancer"]==1,]
DATA=data
rm(data)
####### average pcancer per region:
ppz=nrow(DATA[(DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ"),])/nrow(DATA[(DATA[,"zone"]=="PZ"),])
pcg=nrow(DATA[(DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG"),])/nrow(DATA[(DATA[,"zone"]=="CG"),])

############# mean: binary response
meancpz=colMeans(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
meanccg=colMeans(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
meanncpz=colMeans(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
meannccg=colMeans(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
meannc=colMeans(DATA[DATA[,"cancer"]==0,c("ADC","KTRANS","KEP","AUGC")])
meanc=colMeans(DATA[DATA[,"cancer"]==1,c("ADC","KTRANS","KEP","AUGC")])
means=list(meancpz,meanccg,meanncpz,meannccg)
covcpz=var(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])/1.5
covccg=var(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])/1.5
covncpz=var(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])/1.5
covnccg=var(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])/1.5
covnc=var(DATA[(DATA[,"cancer"]==0),c("ADC","KTRANS","KEP","AUGC")])
covc=var(DATA[(DATA[,"cancer"]==1),c("ADC","KTRANS","KEP","AUGC")])
covs=list(covcpz,covccg,covncpz,covnccg)

Gammas=list()
Gammas[[1]]=diag(c(1,0.05,0.05,0.05))
Gammas[[2]]=diag(c(5,0.18,0.18,0.18))
Gammas[[3]]=diag(c(10,0.36,0.36,0.36))

pcg.pool=c(0.25,0.2,0.15)#small diff, big diff
ppz.pool=c(0.3,0.4,0.55)#small diff, big diff
tausq=1
the1=c(10,4,1) #sigmasq for c_ij
the2=c(2,5,10) #phi for c_ij
the3=c(1.5,0.8,0.5) #nu for c_ij
omega.delta=diag(rep(1,4))
Settings=expand.grid(1:3,1:3)


folds=list()
folds[[1]]=1:7
folds[[2]]=8:14
folds[[3]]=15:21
folds[[4]]=22:28
folds[[5]]=29:34

#############################################################
#############################################################
###### detect the cancer prevalence in each subregion: ######
#############################################################
#############################################################
index=1:34
b.posthoc.glm = b.posthoc.qda = b.posthoc.rf = matrix(NA,nrow(Settings),3)
#b.posthoc[1]=0.14
for (set in c(5,7,8))
{
  ############## simulate data:
  sigmasq_c = the1[Settings[set,1]]
  phi_c = the2[Settings[set,1]]
  nu_c = the3[Settings[set,1]]
  
  Gamma = Gammas[[Settings[set,2]]]
  pcg = pcg.pool[Settings[set,2]]
  ppz = ppz.pool[Settings[set,2]]
  
  ################################################################################
  ################## simulation 1: choosing smoothing bandwidth: #################
  ################################################################################
  load("~/Google Drive/mpMRI/server/mpmri/msi_fillrecenter_missing.RData") # including the locations/voxels that have missing mpMRI parameters
  fillin=fillrecenter
  rm(fillrecenter)
  for (i in 1:length(fillin))
  {
    tem = fillin[[i]]
    Region4=1*((tem$x<=0)&(tem$y>0)) + 2*((tem$x>0)&(tem$y>0)) + 3*((tem$x<=0)&(tem$y<=0))+ 4*((tem$x>0)&(tem$y<=0))
    Region9=1*((tem$x<=2/3-1)&(tem$y>1-2/3)) + 2*((tem$x>2/3-1)&(tem$x<=2*2/3-1)&(tem$y>1-2/3)) + 3*((tem$x>2*2/3-1)&(tem$y>1-2/3))+ 4*((tem$x<=2/3-1)&(tem$y<=1-2/3)&(tem$y>1-2*2/3)) + 5*((tem$x>2/3-1)&(tem$x<=2*2/3-1)&(tem$y<=1-2/3)&(tem$y>1-2*2/3)) + 6*((tem$x>2*2/3-1)&(tem$y<=1-2/3)&(tem$y>1-2*2/3))+ 7*((tem$x<=2/3-1)&(tem$y<=1-2*2/3)) + 8*((tem$x>2/3-1)&(tem$x<=2*2/3-1)&(tem$y<=1-2*2/3)) + 9*((tem$x>2*2/3-1)&(tem$y<=1-2*2/3))
    tem$Region4=Region4
    tem$Region9=Region9
    fillin[[i]] = tem
  }
  ### select 34 out of 46 shapes
  select=sample(1:46,size=34,replace=TRUE)
  ##### simulate data:
  sse4=sse9=list()
  for (i in 1:4)
  {
    sse4[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma)
  }
  for (i in 1:9)
  {
    sse9[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma)
  }
  fillnew = lapply(1:length(select),FUN=function(x){simu_binary_data(fillin,select[x],x,sse4,sse9)})
  DATA=rbindlist(lapply(fillnew,as.data.frame))
  DATA=as.data.frame(DATA)
  Region4=1*((DATA$x<=0)&(DATA$y>0)) + 2*((DATA$x>0)&(DATA$y>0)) + 3*((DATA$x<=0)&(DATA$y<=0))+ 4*((DATA$x>0)&(DATA$y<=0))
  Region9=1*((DATA$x<=2/3-1)&(DATA$y>1-2/3)) + 2*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y>1-2/3)) + 3*((DATA$x>2*2/3-1)&(DATA$y>1-2/3))+ 4*((DATA$x<=2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 5*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 6*((DATA$x>2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3))+ 7*((DATA$x<=2/3-1)&(DATA$y<=1-2*2/3)) + 8*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2*2/3)) + 9*((DATA$x>2*2/3-1)&(DATA$y<=1-2*2/3))
  DATA$Region4=Region4
  DATA$Region9=Region9
  DATA$cancer=as.factor(DATA$cancer)
  DATA=as.data.frame(DATA)
  ni=sapply(1:34,FUN=function(x){nrow(DATA[DATA[,"subject"]==x,])})
  
  #####
  b.posthoc.glm[set,1:3] = b.posthoc.tuning.multiresolution('glm')
  b.posthoc.qda[set,1:3] = b.posthoc.tuning.multiresolution('qda')
  b.posthoc.rf[set,1:3] = b.posthoc.tuning.multiresolution('rf')

  print(paste0("Select ",set,"-th post-hoc smoothing bandwidth: finished"))
  save(b.posthoc.glm,b.posthoc.qda,b.posthoc.rf,file=paste0("~/Google Drive/mpMRI/Paper 3 Super Learner/new.b.posthoc_simu_binary_permethod_multiresolution.RData"))
}
