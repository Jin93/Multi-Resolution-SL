set = 1
SS = 1:100
load(paste0("/home/zhangl4/jinxx493/mpmri/new.b.posthoc_simu_binary_permethod_multiresolution.RData")) 
# for binary, also use this for ordinal classification
setwd("/home/zhangl4/jinxx493/mpmri/")
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
#library(doSNOW)
library(parallel)
S=100

n.cores=20
index=1:40
resolution4=4
resolution9=9
Z=3 # number of different segmentations considered
source("basic_functions.R")
source("sl_functions.R")

base_functions = c("glm","qda","rf")
#base_function = "rf" #the chosen base function 
knn.k=250 # giving the smallest cv.error.rate
n.tree = 136
covsl = "prob" 
n.X= 3 * (covsl == "category") + 6 * (covsl == "prob") #4 # number of spatially varying coefficients # 10: redundant: rank-deficient
### Define MRIpars:
MRIpars = c("ADC","KTRANS","KEP","AUGC")
rf.baseformula = as.formula(paste("category ~ ", paste(MRIpars, collapse= "+")))

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
load("DATA9.RData")
####### average pcancer per region:
ppz=nrow(DATA[(DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ"),])/nrow(DATA[(DATA[,"zone"]=="PZ"),])
pcg=nrow(DATA[(DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG"),])/nrow(DATA[(DATA[,"zone"]=="CG"),])
####### average p0, p1, p2 per region:? necessary?
#fit=polymars(DATA[,c("ADC","KTRANS","KEP","AUGC")], DATA[,c("x","y")], classify = F,maxsize=6,knots=10)

############# mean: binary response
m0pz=colMeans(DATA[((DATA[,"category"]==0)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
m0cg=colMeans(DATA[((DATA[,"category"]==0)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
m1pz=colMeans(DATA[((DATA[,"category"]==1)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
m1cg=colMeans(DATA[((DATA[,"category"]==1)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
m2pz=colMeans(DATA[((DATA[,"category"]==2)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
m2cg=colMeans(DATA[((DATA[,"category"]==2)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
means=list(m0pz,m0cg,m1pz,m1cg,m2pz,m2cg)
c0pz=var(DATA[((DATA[,"category"]==0)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])/1.5
c0cg=var(DATA[((DATA[,"category"]==0)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])/1.5
c1pz=var(DATA[((DATA[,"category"]==1)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])/1.5
c1cg=var(DATA[((DATA[,"category"]==1)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])/1.5
c2pz=var(DATA[((DATA[,"category"]==2)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])/1.5
c2cg=var(DATA[((DATA[,"category"]==2)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])/1.5
covs=list(c0pz,c0cg,c1pz,c1cg,c2pz,c2cg)

Means = list(means)#small diff, big diff
Covs = list(covs)#small diff, big diff
Gammas=list()
Gammas[[1]]=diag(c(1,0.05,0.05,0.05))
Gammas[[2]]=diag(c(5,0.18,0.18,0.18))
Gammas[[3]]=diag(c(10,0.36,0.36,0.36))
pcg.pool=c(0.25,0.2,0.15)#small diff, big diff
ppz.pool=c(0.3,0.4,0.6)#small diff, big diff
tausq=1
the1=c(10,4,1) #sigmasq for c_ij
the2=c(2,5,10) #phi for c_ij
the3=c(1.5,0.8,0.5) #nu for c_ij
omega.delta=diag(rep(1,4))
Settings=expand.grid(1:3,1:3)
t1 = 0.5
t2 = 0.7

folds=list()
folds[[1]]=1:10
folds[[2]]=11:20
folds[[3]]=21:30
folds[[4]]=31:40

#############################################################
#############################################################
###### detect the cancer prevalence in each subregion: ######
#############################################################
#############################################################
index=1:40
#for (set in 1:nrow(Settings))
#{
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

b.posthoc.glm = b.posthoc.glm[set,]
b.posthoc.qda = b.posthoc.qda[set,]
b.posthoc.rf = b.posthoc.rf[set,]

b.posthoc = c(b.posthoc.glm,b.posthoc.qda,b.posthoc.rf)
b.posthoc = rep(b.posthoc,each=2)

simuresults = mclapply(SS,FUN=function(x){fullsimu_category(x,b.posthoc)},mc.cores = S/4)
save(simuresults,file=paste0("/home/zhangl4/jinxx493/mpmri/results-paper3/new.simuresults-category-setting",set,"-",SS[1],"-",SS[length(SS)],".RData"))

