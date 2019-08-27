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
library(hablar) #for s(): remove NaN

n.cores=4
index=1:34
resolution4=4
resolution9=9
Z=2 # number of different segmentations considered
source("basic_functions.R")
sourceCpp("sl_functions.R")

#### Select tuning parameters:
slf= "Unordered" # ordered or Unordered
base_functions = c("glm","qda","rf")
base_function = "glm" #the chosen base function 
knn.k=250 # 143 for categorical outcome # giving the smallest cv.error.rate
n.tree = 136
n.X = 3 # number of X in SL
slcovnames = covariates = c("Global_pred","Region4_pred","Region9_pred")

### Define MRIpars:
presmooth =0 #1: yes, 90:no
if (presmooth==1){MRIpars = c("sadc","sktrans","skep","saugc")}
if (presmooth==0){MRIpars = c("ADC","KTRANS","KEP","AUGC")}
rf.baseformula = as.formula(paste("cancer ~ ", paste(MRIpars, collapse= "+")))

#############################################################
##################### Cross validation: #####################
#############################################################
folds=list()
folds[[1]]=c(2:6,8,34)
folds[[2]]=c(7,9:14)
folds[[3]]=15:21
folds[[4]]=c(22:24,26:28,33)
folds[[5]]=c(1,25,29:32)

tp.spsl=tn.spsl=err.spsl=auc.spsl=s90.spsl=s80.spsl=numeric()
table.spsl=list()
tp.sl=tn.sl=err.sl=auc.sl=s90.sl=s80.sl=numeric()
table.sl=list()
tp.base=tn.base=err.base=auc.base=s90.base=s80.base=numeric()
table.base=list()
pred.mat=matrix(NA,nrow(DATA),3)
colnames(pred.mat) = c("pred.spsl","pred.sl","pred.base")


#################### delete duplicated voxels:
load("DATA9.RData")
DATA$cancer=as.factor(DATA$cancer)
DATA$xnew = DATA$x
DATA$ynew = DATA$y
DATA$x = DATA$x_orig
DATA$y = DATA$y_orig

fill = list()
for (k in index) 
{ 
  fill[[k]] = DATA[DATA$subject == k,]
  temp = fill[[k]]
  #temp = unique(fill[[k]])
  fill[[k]] = temp
}

####### Choosing Smoothing Bandwidth:
b.posthoc.glm = b.posthoc.tuning.multiresolution('glm')
b.posthoc.qda = b.posthoc.tuning.multiresolution('qda')
b.posthoc.rf = b.posthoc.tuning.multiresolution('rf')

B.posthoc = c(b.posthoc.glm,b.posthoc.qda,b.posthoc.rf)
save(B.posthoc,file="~/Google Drive/mpMRI/Paper 3 Super Learner/b.posthoc.realdata.RData")
#########################################
select = index = 1:34
fillnew = fill
ni=sapply(1:34,FUN=function(x){nrow(DATA[DATA[,"subject"]==x,])})
#############################################################
##################### Model Evaluation: #####################
#############################################################
realresults.combine = simu_binary(DATA,"combine",b.posthoc)
save(realresults.combine, file = "~/Google Drive/mpMRI/Paper 3 Super Learner/application_binary_results.RData")

