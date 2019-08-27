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
library(data.table)
n.cores=4
index=1:34
resolution4=4
resolution9=9
Z=3 # number of different segmentations considered
source("basic_functions.R")
source("sl_functions.R")
slf= "Unordered" # ordered or Unordered  
base_function = "rf" #the chosen base function 
base_functions = c("glm","qda","rf")
knn.k=143 # 143 for categorical outcome # giving the smallest cv.error.rate
n.tree = 136

covsl = "prob" # or prob or category
n.X= 3 * (covsl == "category") + 6 * (covsl == "prob") #4 # number of spatially varying coefficients # 10: redundant: rank-deficient
if (covsl == "category") {slcovnames = covariates = c("Global_pred","Region4_pred","Region9_pred")}
if (covsl == "prob") {slcovnames = covariates = c("Global_pred1","Global_pred2","Region4_pred1","Region4_pred2","Region9_pred1","Region9_pred2")}

### Define MRIpars:
MRIpars = c("ADC","KTRANS","KEP","AUGC")
rf.baseformula = as.formula(paste("category ~ ", paste(MRIpars, collapse= "+")))

################################################################################
########################## Cross validation ###################################
################################################################################
load("DATA9.RData")
DATA$cancer=as.factor(DATA$cancer)
DATA$xnew = DATA$x
DATA$ynew = DATA$y
DATA$x = DATA$x_orig
DATA$y = DATA$y_orig
folds=list()
folds[[1]]=c(1:8)
folds[[2]]=c(10:17)
folds[[3]]=c(9,18:24)
folds[[4]]=25:34
ni=sapply(1:34,FUN=function(x){nrow(DATA[DATA[,"subject"]==x,])})
presmooth=0

fill = list()
for (k in index) 
{ 
  fill[[k]] = DATA[DATA$subject == k,]
  temp = fill[[k]]
  #temp = unique(fill[[k]])
  fill[[k]] = temp
}

load("~/Google Drive/mpMRI/Paper 3 Super Learner/b.posthoc.realdata.RData")
base_function = "glm"

if (base_function == "glm")
{
  if (covsl == "prob") {b.posthoc = rep(B.posthoc[1:3],each=2)}
  if (covsl == "category") {b.posthoc = B.posthoc[1:3]}
}
if (base_function == "qda")
{
  if (covsl == "prob") {b.posthoc = rep(B.posthoc[4:6],each=2)}
  if (covsl == "category") {b.posthoc = B.posthoc[4:6]}
}
if (base_function == "rf")
{
  if (covsl == "prob") {b.posthoc = rep(B.posthoc[7:9],each=2)}
  if (covsl == "category") {b.posthoc = B.posthoc[7:9]}
}
if (base_function == "combine")
{
  if (covsl == "prob") {b.posthoc = rep(B.posthoc,each=2)}
  if (covsl == "category") {b.posthoc = B.posthoc}
}
Levels = sort(unique(DATA$category))

if (base_function != "combine")
{
  if (covsl == "category") {slcovnames = covariates = c("Global_pred","Region4_pred","Region9_pred")}
  if (covsl == "prob") {slcovnames = covariates = c("Global_pred1","Global_pred2","Region4_pred1","Region4_pred2","Region9_pred1","Region9_pred2")}
}
if (base_function == "combine")
{
  if (covsl == "category") {slcovnames = covariates =   c(sapply(1:length(base_functions),FUN=function(x){c(paste0("Global_pred_",base_functions[x]),paste0("Region4_pred_",base_functions[x]),paste0("Region9_pred_",base_functions[x]))}))}
  if (covsl == "prob") {slcovnames = covariates =   c(sapply(1:length(base_functions),FUN=function(x){c(paste0("Global_pred1_",base_functions[x]),paste0("Global_pred2_",base_functions[x]),paste0("Region4_pred1_",base_functions[x]),paste0("Region4_pred2_",base_functions[x]),paste0("Region9_pred1_",base_functions[x]),paste0("Region9_pred2_",base_functions[x]))}))}
}
pred.mat=matrix(NA,nrow(DATA),5)
colnames(pred.mat) = paste0("pred.",c("base","sl10","sl1","sl20","sl2"))
for (split in 1:length(folds))
{
  TRAIN=DATA[!(DATA$subject%in%folds[[split]]),]
  tr.subj = unique(TRAIN$subject)
  holdout = DATA[(DATA$subject%in%folds[[split]]),]
  holdout.ind = which(DATA[,"subject"]%in%folds[[split]])
  n.holdout = nrow(holdout)
  n.whole=nrow(TRAIN)
  p=3
 
  ########### Weights for GLM:
  weights.train=list()
  weights.train[[1]] = rep(1,n.whole)
  weights.train[[2]] = 1/sum(TRAIN$category==0) * (TRAIN$category==0)+1/sum(TRAIN$category==1) * (TRAIN$category==1)+1/sum(TRAIN$category==2) * (TRAIN$category==2)
  weights.train[[2]]=weights.train[[2]]/(sum(weights.train[[2]]))
  
  ############## Obtain the new X for super learner:
  if (base_function != "combine")
  {
    if (covsl == "category") 
    {
      whole=cbind(TRAIN,Global_pred=0,Region4_pred=0,Region9_pred=0)
    }
    if (covsl == "prob") 
    {
      whole=cbind(TRAIN,Global_pred1=0,Global_pred2=0,Region4_pred1=0,Region4_pred2=0,Region9_pred1=0,Region9_pred2=0)
    }
  }
  if (base_function == "combine")
  {
    tem = matrix(NA,nrow(TRAIN),n.X*length(base_functions))
    colnames(tem) = slcovnames
    whole = data.frame(TRAIN, tem)
  }  

  TRAIN.folds = c(1:length(folds))[-split]
  for (fold in 1:length(TRAIN.folds))
  {
    k=folds[[TRAIN.folds[fold]]]
    ind=index[-k]
    train=TRAIN[!TRAIN$subject%in%k,]
    test=TRAIN[TRAIN$subject%in%k,]
    test.rownames=rownames(test)
    
    if (base_function != "combine")
    {
      regional_prediction_test=matrix(0,nrow(test),n.X)
    }
    if (base_function == "combine")
    {
      regional_prediction_test=matrix(0,nrow(test),n.X*length(base_functions))
    }
    
    MRIpars_whole <- TRAIN[,MRIpars]
    MRIpars_train <- train[,MRIpars]
    MRIpars_test <- test[,MRIpars]
    
    ##### Global model:
    ###### train the candidate model on the training blocks:
    if (base_function == "rf")
    {
      global_fit <- randomForest(rf.baseformula, data = train, ntree=n.tree, mtry=2,
                                 nodesize=1, importance=T)
      if (covsl == "prob")
      {
        global_pred <- predict(global_fit, newdata=test,type="prob")
        pred.test=global_pred
        regional_prediction_test[,1:2]<- pred.test[,c(1,2)]
      } 
      if (covsl == "category")
      {
        global_pred <- predict(global_fit, newdata=test,type="response")
        pred.test=as.numeric(as.character(global_pred)) #factor -> numerical
        regional_prediction_test[,1]<- pred.test
      }
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                  nodesize=1, importance=T)
        if (covsl == "prob")
        {
          local_pred=predict(local_fit, newdata=test,type="prob")
          pred.test=local_pred
          regional_prediction_test[,3:4]<- pred.test[,c(1,2)] * (test[,"Region4"]==reg) + regional_prediction_test[,3:4]
        }
        if (covsl == "category")
        {
          local_pred=predict(local_fit, newdata=test,type="response")
          pred.test=as.numeric(as.character(local_pred)) #factor -> numerical
          regional_prediction_test[,2]<- pred.test * (test[,"Region4"]==reg) + regional_prediction_test[,2]
        }
      }
      
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                  nodesize=1, importance=T)
        if (covsl == "prob")
        {
          local_pred=predict(local_fit, newdata=test,type="prob")
          pred.test=local_pred
          regional_prediction_test[,5:6]<- pred.test[,c(1,2)] * (test[,"Region9"]==reg) + regional_prediction_test[,5:6]
        }
        if (covsl == "category")
        {
          local_pred=predict(local_fit, newdata=test,type="response")
          pred.test=as.numeric(as.character(local_pred)) #factor -> numerical
          regional_prediction_test[,3]<- pred.test * (test[,"Region9"]==reg) + regional_prediction_test[,3]
        }
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    
    
    if (base_function == "knn")
    {
      global_fit=knn(train=train[,MRIpars],test=test[,MRIpars],cl=train$category,k=knn.k,prob=FALSE)
      regional_prediction_test[,1]<- as.numeric(as.character(global_fit))
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- knn(train=regional_train[,MRIpars],test=test[,MRIpars],cl=regional_train$category,k=knn.k,prob=FALSE)
        local_pred = as.numeric(as.character(local_fit))
        regional_prediction_test[,2] <- local_pred * (test[,"Region4"]==reg) + regional_prediction_test[,2]
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- knn(train=regional_train[,MRIpars],test=test[,MRIpars],cl=regional_train$category,k=knn.k,prob=FALSE)
        local_pred = as.numeric(as.character(local_fit))
        regional_prediction_test[,3] <- local_pred * (test[,"Region9"]==reg) + regional_prediction_test[,3]
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    
    if (base_function == "qda")
    {
      global_fit<-qda(train[,MRIpars], train$category)
      global_pred<-predict(global_fit, test[,MRIpars])
      if (covsl == "prob")
      {
        pred.test=global_pred$posterior
        regional_prediction_test[,1:2]<- pred.test[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.test=as.numeric(as.character(global_pred$class)) #factor -> numerical
        regional_prediction_test[,1]<- pred.test
      }
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- qda(regional_train[,MRIpars], regional_train$category)
        local_pred=predict(local_fit,test[,MRIpars])
        if (covsl == "prob")
        {
          pred.test=local_pred$posterior
          regional_prediction_test[,3:4]<- pred.test[,c(1,2)] * (test[,"Region4"]==reg) + regional_prediction_test[,3:4]
        }
        if (covsl == "category")
        {
          pred.test=as.numeric(as.character(local_pred$class))
          regional_prediction_test[,2]<- pred.test * (test[,"Region4"]==reg) + regional_prediction_test[,2]
        }
      }
      
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- qda(regional_train[,MRIpars], regional_train$category) 
        local_pred=predict(local_fit,test[,MRIpars])
        if (covsl == "prob")
        {
          pred.test=local_pred$posterior
          regional_prediction_test[,5:6]<- pred.test[,c(1,2)] * (test[,"Region9"]==reg) + regional_prediction_test[,5:6]
        }
        if (covsl == "category")
        {
          pred.test=as.numeric(as.character(local_pred$class))
          regional_prediction_test[,3]<- pred.test * (test[,"Region9"]==reg) + regional_prediction_test[,3]
        }
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    if (base_function == "glm")
    {
      global_fit <- polr(rf.baseformula, data = train, method = "probit",Hess=TRUE) 
      if (covsl == "prob")
      {
        global_pred <- predict(global_fit, newdata=test,type="probs")
        pred.test=global_pred
        regional_prediction_test[,1:2]<- pred.test[,c(1,2)]
      } 
      if (covsl == "category")
      {
        global_pred <- predict(global_fit, newdata=test,type="class")
        pred.test=as.numeric(as.character(global_pred)) 
        regional_prediction_test[,1]<- pred.test
      }
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- polr(rf.baseformula, data = regional_train, method = "probit",Hess=TRUE) 
        if (covsl == "prob")
        {
          local_pred=predict(local_fit, newdata=test,type="prob")
          pred.test=local_pred
          regional_prediction_test[,3:4]<- pred.test[,c(1,2)] * (test[,"Region4"]==reg) + regional_prediction_test[,3:4]
        }
        if (covsl == "category")
        {
          local_pred=predict(local_fit, newdata=test,type="class")
          pred.test=as.numeric(as.character(local_pred)) 
          regional_prediction_test[,2]<- pred.test * (test[,"Region4"]==reg) + regional_prediction_test[,2]
        }
      }
      
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- polr(rf.baseformula, data = regional_train, method = "probit",Hess=TRUE) 
        if (covsl == "prob")
        {
          local_pred=predict(local_fit, newdata=test,type="prob")
          pred.test=local_pred
          regional_prediction_test[,5:6]<- pred.test[,c(1,2)] * (test[,"Region9"]==reg) + regional_prediction_test[,5:6]
        }
        if (covsl == "category")
        {
          local_pred=predict(local_fit, newdata=test,type="class")
          pred.test=as.numeric(as.character(local_pred)) 
          regional_prediction_test[,3]<- pred.test * (test[,"Region9"]==reg) + regional_prediction_test[,3]
        }
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    
    if (base_function == "combine")
    {
      #####################
      ######## GLM: #######
      #####################
      global_fit <- polr(rf.baseformula, data = train, method = "probit",Hess=TRUE) 
      if (covsl == "prob")
      {
        global_pred <- predict(global_fit, newdata=test,type="probs")
        pred.test=global_pred
        regional_prediction_test[,1:2]<- pred.test[,c(1,2)]
      } 
      if (covsl == "category")
      {
        global_pred <- predict(global_fit, newdata=test,type="class")
        pred.test=as.numeric(as.character(global_pred)) 
        regional_prediction_test[,1]<- pred.test
      }
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- polr(rf.baseformula, data = regional_train, method = "probit",Hess=TRUE) 
        if (covsl == "prob")
        {
          local_pred=predict(local_fit, newdata=test,type="prob")
          pred.test=local_pred
          regional_prediction_test[,3:4]<- pred.test[,c(1,2)] * (test[,"Region4"]==reg) + regional_prediction_test[,3:4]
        }
        if (covsl == "category")
        {
          local_pred=predict(local_fit, newdata=test,type="class")
          pred.test=as.numeric(as.character(local_pred)) 
          regional_prediction_test[,2]<- pred.test * (test[,"Region4"]==reg) + regional_prediction_test[,2]
        }
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- polr(rf.baseformula, data = regional_train, method = "probit",Hess=TRUE) 
        if (covsl == "prob")
        {
          local_pred=predict(local_fit, newdata=test,type="prob")
          pred.test=local_pred
          regional_prediction_test[,5:6]<- pred.test[,c(1,2)] * (test[,"Region9"]==reg) + regional_prediction_test[,5:6]
        }
        if (covsl == "category")
        {
          local_pred=predict(local_fit, newdata=test,type="class")
          pred.test=as.numeric(as.character(local_pred)) 
          regional_prediction_test[,3]<- pred.test * (test[,"Region9"]==reg) + regional_prediction_test[,3]
        }
      }
      #####################
      ######## QDA: #######
      #####################
      global_fit<-qda(train[,MRIpars], train$category)
      global_pred<-predict(global_fit, test[,MRIpars])
      if (covsl == "prob")
      {
        pred.test=global_pred$posterior
        regional_prediction_test[,7:8]<- pred.test[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.test=as.numeric(as.character(global_pred$class)) #factor -> numerical
        regional_prediction_test[,4]<- pred.test
      }
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- qda(regional_train[,MRIpars], regional_train$category)
        local_pred=predict(local_fit,test[,MRIpars])
        if (covsl == "prob")
        {
          pred.test=local_pred$posterior
          regional_prediction_test[,9:10]<- pred.test[,c(1,2)] * (test[,"Region4"]==reg) + regional_prediction_test[,9:10]
        }
        if (covsl == "category")
        {
          pred.test=as.numeric(as.character(local_pred$class))
          regional_prediction_test[,5]<- pred.test * (test[,"Region4"]==reg) + regional_prediction_test[,5]
        }
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- qda(regional_train[,MRIpars], regional_train$category) 
        local_pred=predict(local_fit,test[,MRIpars])
        if (covsl == "prob")
        {
          pred.test=local_pred$posterior
          regional_prediction_test[,11:12]<- pred.test[,c(1,2)] * (test[,"Region9"]==reg) + regional_prediction_test[,11:12]
        }
        if (covsl == "category")
        {
          pred.test=as.numeric(as.character(local_pred$class))
          regional_prediction_test[,6]<- pred.test * (test[,"Region9"]==reg) + regional_prediction_test[,6]
        }
      }
      #####################
      ######## RF：########
      #####################
      global_fit <- randomForest(rf.baseformula, data = train, ntree=n.tree, mtry=2,
                                 nodesize=1, importance=T)
      if (covsl == "prob")
      {
        global_pred <- predict(global_fit, newdata=test,type="prob")
        pred.test=global_pred
        regional_prediction_test[,13:14]<- pred.test[,c(1,2)]
      } 
      if (covsl == "category")
      {
        global_pred <- predict(global_fit, newdata=test,type="response")
        pred.test=as.numeric(as.character(global_pred)) #factor -> numerical
        regional_prediction_test[,7]<- pred.test
      }
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                  nodesize=1, importance=T)
        if (covsl == "prob")
        {
          local_pred=predict(local_fit, newdata=test,type="prob")
          pred.test=local_pred
          regional_prediction_test[,15:16]<- pred.test[,c(1,2)] * (test[,"Region4"]==reg) + regional_prediction_test[,15:16]
        }
        if (covsl == "category")
        {
          local_pred=predict(local_fit, newdata=test,type="response")
          pred.test=as.numeric(as.character(local_pred)) #factor -> numerical
          regional_prediction_test[,8]<- pred.test * (test[,"Region4"]==reg) + regional_prediction_test[,8]
        }
      }
      
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                  nodesize=1, importance=T)
        if (covsl == "prob")
        {
          local_pred=predict(local_fit, newdata=test,type="prob")
          pred.test=local_pred
          regional_prediction_test[,17:18]<- pred.test[,c(1,2)] * (test[,"Region9"]==reg) + regional_prediction_test[,17:18]
        }
        if (covsl == "category")
        {
          local_pred=predict(local_fit, newdata=test,type="response")
          pred.test=as.numeric(as.character(local_pred)) 
          regional_prediction_test[,9]<- pred.test * (test[,"Region9"]==reg) + regional_prediction_test[,9]
        }
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    print(fold)
  }

  ##### CV training:
  ####################################################################################
  ####################################################################################
  ###################  Train candidate models on the entire dataset: #################
  ####################################################################################
  ############ Evaluate the trained super learner on the holdout dataset: ############
  ####################################################################################
  if (base_function != "combine")
  {
    tem = matrix(NA,nrow(holdout),n.X)
    colnames(tem) = slcovnames
    holdout = data.frame(holdout, tem)
  }  
  if (base_function == "combine")
  {
    tem = matrix(NA,nrow(holdout),n.X*length(base_functions))
    colnames(tem) = slcovnames
    holdout = data.frame(holdout, tem)
  }  
  
  if (base_function == "rf")
  {
    baselearner.global <- randomForest(rf.baseformula, data = TRAIN, ntree=n.tree, mtry=2,
                                       nodesize=1, importance=T)
    ##### Local model with 4 regions (2 by 2 split):
    baselearner.local4=list()
    for (reg in 1:resolution4)
    {
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      baselearner.local4[[reg]] <- randomForest(rf.baseformula, data = regional_whole, ntree=n.tree, mtry=2,
                                                nodesize=1, importance=T)
    }
    ##### Local model with 9 regions (3 by 3 split):
    baselearner.local9=list()
    for (reg in 1:resolution9)
    {
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      baselearner.local9[[reg]] <- randomForest(rf.baseformula, data = regional_whole, ntree=n.tree, mtry=2,
                                                nodesize=1, importance=T)
    }
    if (covsl == "prob")
    {
      pred.global <- predict(baselearner.global, newdata=holdout,type="prob")
      holdout[,slcovnames[1:2]]<- pred.global[,c(1,2)]
    }
    if (covsl == "category")
    {
      pred.global <- predict(baselearner.global, newdata=holdout,type="response") #default type: response
      holdout[,slcovnames[1]]<- as.numeric(as.character(pred.global))
    }
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      if (covsl == "prob")
      {
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="prob")
        holdout[region4.ind,slcovnames[3:4]]<- pred.local4[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="response")
        holdout[region4.ind,slcovnames[2]]<- as.numeric(as.character(pred.local4))
      }
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      if (covsl == "prob")
      {
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="prob")
        holdout[region9.ind,slcovnames[5:6]]<- pred.local9[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="response")
        holdout[region9.ind,slcovnames[3]]<- as.numeric(as.character(pred.local9))
      }
    }
  }
  
  if (base_function == "knn")
  {
    pred.global <- knn(train=TRAIN[,MRIpars],test=holdout[,MRIpars],cl=TRAIN$category,k=knn.k,prob=FALSE)
    holdout[,slcovnames[1]]<- as.numeric(as.character(pred.global))
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      pred.local4 <- knn(train=regional_whole[,MRIpars],test=regional_holdout[,MRIpars],cl=regional_whole$category,k=knn.k,prob=FALSE)
      holdout[region4.ind,slcovnames[2]]<- as.numeric(as.character(pred.local4))
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      pred.local9 <- knn(train=regional_whole[,MRIpars],test=regional_holdout[,MRIpars],cl=regional_whole$category,k=knn.k,prob=FALSE)
      holdout[region9.ind,slcovnames[3]]<- as.numeric(as.character(pred.local9))
    }
  }
  
  if (base_function == "qda")
  {
    global.fit = qda(TRAIN[,MRIpars],TRAIN$category)
    pred.global = predict(global.fit, holdout[,MRIpars])
    
    if (covsl == "prob")
    {
      holdout[,slcovnames[1:2]]<- pred.global$posterior[,c(1,2)]
    }
    if (covsl == "category")
    {
      holdout[,slcovnames[1]]<- as.numeric(as.character(pred.global$class))
    }
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      local.fit4 = qda(regional_whole[,MRIpars],regional_whole$category)
      pred.local4 = predict(local.fit4, regional_holdout[,MRIpars])
      
      if (covsl == "prob")
      {
        holdout[region4.ind,slcovnames[3:4]]<- pred.local4$posterior[,c(1,2)]
      }
      if (covsl == "category")
      {
        holdout[region4.ind,slcovnames[2]]<- as.numeric(as.character(pred.local4$class))
      }
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      local.fit9 = qda(regional_whole[,MRIpars],regional_whole$category)
      pred.local9 = predict(local.fit9, regional_holdout[,MRIpars])
      if (covsl == "prob")
      {
        holdout[region9.ind,slcovnames[5:6]]<- pred.local9$posterior[,c(1,2)]
      }
      if (covsl == "category")
      {
        holdout[region9.ind,slcovnames[3]]<- as.numeric(as.character(pred.local9$class))
      }
    }
  }
  #################### glm ######################
  if (base_function == "glm")
  {
    baselearner.global <- polr(rf.baseformula, data = TRAIN, method = "probit",Hess=TRUE) 
    ##### Local model with 4 regions (2 by 2 split):
    baselearner.local4=list()
    for (reg in 1:resolution4)
    {
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      baselearner.local4[[reg]] <- polr(rf.baseformula, data = regional_whole, method = "probit",Hess=TRUE) 
    }
    ##### Local model with 9 regions (3 by 3 split):
    baselearner.local9=list()
    for (reg in 1:resolution9)
    {
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      baselearner.local9[[reg]] <- polr(rf.baseformula, data = regional_whole, method = "probit",Hess=TRUE) 
    }
    ############# predictions:
    if (covsl == "prob")
    {
      pred.global <- predict(baselearner.global, newdata=holdout,type="prob")
      holdout[,slcovnames[1:2]]<- pred.global[,c(1,2)]
    }
    if (covsl == "category")
    {
      pred.global <- predict(baselearner.global, newdata=holdout,type="class") #default type: response
      holdout[,slcovnames[1]]<- as.numeric(as.character(pred.global))
    }
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      if (covsl == "prob")
      {
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="prob")
        holdout[region4.ind,slcovnames[3:4]]<- pred.local4[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="class")
        holdout[region4.ind,slcovnames[2]]<- as.numeric(as.character(pred.local4))
      }
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      if (covsl == "prob")
      {
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="prob")
        holdout[region9.ind,slcovnames[5:6]]<- pred.local9[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="class")
        holdout[region9.ind,slcovnames[3]]<- as.numeric(as.character(pred.local9))
      }
    }
  }
  #################### combine ###################
  if(base_function == "combine")
  {
    ##########################
    ########### GLM ##########
    ##########################
    baselearner.global <- polr(rf.baseformula, data = TRAIN, method = "probit",Hess=TRUE) 
    ##### Local model with 4 regions (2 by 2 split):
    baselearner.local4=list()
    for (reg in 1:resolution4)
    {
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      baselearner.local4[[reg]] <- polr(rf.baseformula, data = regional_whole, method = "probit",Hess=TRUE) 
    }
    ##### Local model with 9 regions (3 by 3 split):
    baselearner.local9=list()
    for (reg in 1:resolution9)
    {
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      baselearner.local9[[reg]] <- polr(rf.baseformula, data = regional_whole, method = "probit",Hess=TRUE) 
    }
    ############# predictions:
    if (covsl == "prob")
    {
      pred.global <- predict(baselearner.global, newdata=holdout,type="prob")
      holdout[,slcovnames[1:2]]<- pred.global[,c(1,2)]
    }
    if (covsl == "category")
    {
      pred.global <- predict(baselearner.global, newdata=holdout,type="class") #default type: response
      holdout[,slcovnames[1]]<- as.numeric(as.character(pred.global))
    }
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      if (covsl == "prob")
      {
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="prob")
        holdout[region4.ind,slcovnames[3:4]]<- pred.local4[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="class")
        holdout[region4.ind,slcovnames[2]]<- as.numeric(as.character(pred.local4))
      }
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      if (covsl == "prob")
      {
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="prob")
        holdout[region9.ind,slcovnames[5:6]]<- pred.local9[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="class")
        holdout[region9.ind,slcovnames[3]]<- as.numeric(as.character(pred.local9))
      }
    }
    ##########################
    ########### QDA ##########
    ##########################
    global.fit = qda(TRAIN[,MRIpars],TRAIN$category)
    pred.global = predict(global.fit, holdout[,MRIpars])
    if (covsl == "prob")
    {
      holdout[,slcovnames[7:8]]<- pred.global$posterior[,c(1,2)]
    }
    if (covsl == "category")
    {
      holdout[,slcovnames[4]]<- as.numeric(as.character(pred.global$class))
    }
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      local.fit4 = qda(regional_whole[,MRIpars],regional_whole$category)
      pred.local4 = predict(local.fit4, regional_holdout[,MRIpars])
      if (covsl == "prob")
      {
        holdout[region4.ind,slcovnames[9:10]]<- pred.local4$posterior[,c(1,2)]
      }
      if (covsl == "category")
      {
        holdout[region4.ind,slcovnames[5]]<- as.numeric(as.character(pred.local4$class))
      }
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      local.fit9 = qda(regional_whole[,MRIpars],regional_whole$category)
      pred.local9 = predict(local.fit9, regional_holdout[,MRIpars])
      if (covsl == "prob")
      {
        holdout[region9.ind,slcovnames[11:12]]<- pred.local9$posterior[,c(1,2)]
      }
      if (covsl == "category")
      {
        holdout[region9.ind,slcovnames[6]]<- as.numeric(as.character(pred.local9$class))
      }
    }
    ##########################
    ########### RF ###########
    ##########################
    baselearner.global <- randomForest(rf.baseformula, data = TRAIN, ntree=n.tree, mtry=2,
                                       nodesize=1, importance=T)
    ##### Local model with 4 regions (2 by 2 split):
    baselearner.local4=list()
    for (reg in 1:resolution4)
    {
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      baselearner.local4[[reg]] <- randomForest(rf.baseformula, data = regional_whole, ntree=n.tree, mtry=2,
                                                nodesize=1, importance=T)
    }
    ##### Local model with 9 regions (3 by 3 split):
    baselearner.local9=list()
    for (reg in 1:resolution9)
    {
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      baselearner.local9[[reg]] <- randomForest(rf.baseformula, data = regional_whole, ntree=n.tree, mtry=2,
                                                nodesize=1, importance=T)
    }
    if (covsl == "prob")
    {
      pred.global <- predict(baselearner.global, newdata=holdout,type="prob")
      holdout[,slcovnames[13:14]]<- pred.global[,c(1,2)]
    }
    if (covsl == "category")
    {
      pred.global <- predict(baselearner.global, newdata=holdout,type="response") #default type: response
      holdout[,slcovnames[7]]<- as.numeric(as.character(pred.global))
    }
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      if (covsl == "prob")
      {
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="prob")
        holdout[region4.ind,slcovnames[15:16]]<- pred.local4[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="response")
        holdout[region4.ind,slcovnames[8]]<- as.numeric(as.character(pred.local4))
      }
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      if (covsl == "prob")
      {
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="prob")
        holdout[region9.ind,slcovnames[17:18]]<- pred.local9[,c(1,2)]
      }
      if (covsl == "category")
      {
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="response")
        holdout[region9.ind,slcovnames[9]]<- as.numeric(as.character(pred.local9))
      }
    }
    
  }
  
  ############################################################################
  ########################## SL without Smoothing ############################
  ############################################################################
  ####### Model fitting using the new predictors:
  nospatial_polr_fomula = as.formula(paste("category ~ ", paste(slcovnames, collapse= "+")))
  
  sl.model <- polr(nospatial_polr_fomula, data = whole, method = "probit",Hess=TRUE,weights=weights.train[[1]]) 
  ####### Prediction:
  prediction.superlearner<-predict(sl.model, holdout, type="class")
  pred.sl10<- as.numeric(as.character(prediction.superlearner))
  holdout = data.frame(holdout, pred.sl10 = pred.sl10)
  
  sl.model <- polr(nospatial_polr_fomula, data = whole, method = "probit",Hess=TRUE,weights=weights.train[[2]]) 
  ####### Prediction:
  prediction.superlearner<-predict(sl.model, holdout, type="class")
  pred.sl20<- as.numeric(as.character(prediction.superlearner))
  holdout = data.frame(holdout, pred.sl20 = pred.sl20)
  
  ################################################################
  ########################## Smoothing: ##########################
  ################################################################
  wholesmooth = whole
  for (ii in tr.subj)
  {
    tem = whole[whole$subject == ii,]
    row.ind = which(whole$subject == ii)
    Xcoord=tem[,"x_orig"];Ycoord=tem[,"y_orig"];
    for (jj in 1:length(slcovnames))
    {
      ######## smooth pij:
      X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[jj]])
      #b=0.2
      #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
      Spre=as.numeric(Smooth(X, sigma=b.posthoc[jj],at="points",edge=TRUE, diggle=F))
      tem[,covariates[jj]] = Spre
      wholesmooth[row.ind,covariates[jj]] = Spre
      ######## somooth qij:
    }
  }
  whole=wholesmooth
  rm(wholesmooth)
  holdoutsmooth = holdout
  for (ii in folds[[split]])
  {
    tem = holdoutsmooth[holdoutsmooth$subject == ii,]
    row.ind = which(holdoutsmooth$subject == ii)
    Xcoord=round(tem[,"x_orig"],2);Ycoord=round(tem[,"y_orig"],2);
    for (jj in 1:length(slcovnames))
    {
      ######## smooth pij:
      X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[jj]])
      #b=0.2
      #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
      Spre=as.numeric(Smooth(X, sigma=b.posthoc[jj],at="points",edge=TRUE, diggle=F))
      tem[,covariates[jj]] = Spre
      holdoutsmooth[row.ind,covariates[jj]] = Spre
      ######## somooth qij:
    }
  }
  holdout = holdoutsmooth
  rm(holdoutsmooth)
  categorytrue = holdout[,"category"]  ##### notice that in each fold: should update "category" for the test slices
  
  ############################################################################
  ########################### SL with Smoothing ##############################
  ############################################################################
  ####### Model fitting using the new predictors:
  nospatial_polr_fomula = as.formula(paste("category ~ ", paste(slcovnames, collapse= "+")))
  
  sl.model <- polr(nospatial_polr_fomula, data = whole, method = "probit",Hess=TRUE,weights=weights.train[[1]]) 
  ####### Prediction:
  prediction.superlearner<-predict(sl.model, holdout,type="class")
  pred.sl1<- as.numeric(as.character(prediction.superlearner))
  holdout = data.frame(holdout, pred.sl1 = pred.sl1)
  
  sl.model <- polr(nospatial_polr_fomula, data = whole, method = "probit",Hess=TRUE,weights=weights.train[[2]]) 
  ####### Prediction:
  prediction.superlearner<-predict(sl.model, holdout,type="class")
  pred.sl2<- as.numeric(as.character(prediction.superlearner))
  holdout = data.frame(holdout, pred.sl2 = pred.sl2)

  ###########################################################
  ######################## Baseline: ########################
  ###########################################################
  if (base_function != "combine")
  {
    if (base_function == "knn")
    {
      baseline.model <- knn(train=TRAIN[,MRIpars],test=holdout[,MRIpars],cl=TRAIN$category,k=knn.k,prob=FALSE)
      pred.base <- as.numeric(as.character(baseline.model))
    }
    if (base_function == "rf")
    {
      baseline.model <- randomForest(rf.baseformula, data = TRAIN, ntree=n.tree, mtry=2,
                                     nodesize=1, importance=T)
      baseline.pred <- predict(baseline.model, newdata=holdout, type="response")
      pred.base <- as.numeric(as.character(baseline.pred))
    }
    if (base_function == "qda")
    {
      baseline.model <- qda(TRAIN[,MRIpars], TRAIN$category)
      base_pred <-predict(baseline.model, holdout[,MRIpars])
      pred.base <- as.numeric(as.character(base_pred$class))
    }
    if (base_function == "spline")
    {
      global_fit<-polymars(train$category, MRIpars_train, classify = T, maxsize=t[tn,1], knots=t[tn,2])
      global_pred<-predict(global_fit, x = MRIpars_test, classify = T)
      base_prediction[test.rows]<- as.numeric(as.character(global_pred))-1
    }
    if (base_function == "glm")
    {
      baseline.model <- polr(rf.baseformula, data = TRAIN, method = "probit",Hess=TRUE)
      baseline.pred <- predict(baseline.model, newdata=holdout, type="class")
      pred.base <- as.numeric(as.character(baseline.pred))
    }
    #if (base_function == "kqr") # kernel quantile regression or quantile regression?
    holdout = data.frame(holdout, pred.base = pred.base)
    pred.mat[holdout.ind,paste0("pred.",c("base","sl10","sl1","sl20","sl2"))] = as.matrix(holdout[,paste0("pred.",c("base","sl10","sl1","sl20","sl2"))])
  }
  pred.mat[holdout.ind,paste0("pred.",c("sl10","sl1","sl20","sl2"))] = as.matrix(holdout[,paste0("pred.",c("sl10","sl1","sl20","sl2"))])
}
###### overall tables:
otable.base = evaluation_table(DATA$category,pred.mat[,"pred.base"],Levels)
otable.sl10 = evaluation_table(DATA$category,pred.mat[,"pred.sl10"],Levels)
otable.sl20 = evaluation_table(DATA$category,pred.mat[,"pred.sl20"],Levels)
otable.sl1 = evaluation_table(DATA$category,pred.mat[,"pred.sl1"],Levels)
otable.sl2 = evaluation_table(DATA$category,pred.mat[,"pred.sl2"],Levels)
###### overall error rates:
oerr.base = 1 - sum(DATA$category == pred.mat[,"pred.base"])/nrow(DATA)
oerr.sl10 = 1 - sum(DATA$category == pred.mat[,"pred.sl10"])/nrow(DATA)
oerr.sl20 = 1 - sum(DATA$category == pred.mat[,"pred.sl20"])/nrow(DATA)
oerr.sl1 = 1 - sum(DATA$category == pred.mat[,"pred.sl1"])/nrow(DATA)
oerr.sl2 = 1 - sum(DATA$category == pred.mat[,"pred.sl2"])/nrow(DATA)
oerr = round(c(oerr.base,oerr.sl10,oerr.sl1,oerr.sl20,oerr.sl2),3)
names(oerr) = c("base","sl10","sl1","sl20","sl2")


###### overall fpr:
ofpr=round(matrix(c(sapply(1:Z,FUN=function(z){1-sum(otable.base[z,z])/sum(DATA$category==(z-1))}),
                    sapply(1:Z,FUN=function(z){1-sum(otable.sl10[z,z])/sum(DATA$category==(z-1))}),
                    sapply(1:Z,FUN=function(z){1-sum(otable.sl1[z,z])/sum(DATA$category==(z-1))}),
                    sapply(1:Z,FUN=function(z){1-sum(otable.sl20[z,z])/sum(DATA$category==(z-1))}),
                    sapply(1:Z,FUN=function(z){1-sum(otable.sl2[z,z])/sum(DATA$category==(z-1))})),
                           nrow=5,ncol=Z,byrow=TRUE),3)
rownames(ofpr) = c("base","sl10","sl1","sl20","sl2")
ofdr=round(matrix(c(sapply(1:Z,FUN=function(z){1-sum(otable.base[z,z])/sum(otable.base[,z])}),
                    sapply(1:Z,FUN=function(z){1-sum(otable.sl10[z,z])/sum(otable.sl10[,z])}),
                    sapply(1:Z,FUN=function(z){1-sum(otable.sl1[z,z])/sum(otable.sl1[,z])}),
                    sapply(1:Z,FUN=function(z){1-sum(otable.sl20[z,z])/sum(otable.sl20[,z])}),
                    sapply(1:Z,FUN=function(z){1-sum(otable.sl2[z,z])/sum(otable.sl2[,z])})),
                  nrow=5,ncol=Z,byrow=TRUE),3)
rownames(ofdr) = c("base","sl10","sl1","sl20","sl2")
save(otable.base, otable.sl10,otable.sl20, otable.sl1,otable.sl2, oerr, ofpr,ofdr,
    file=paste0("results/category/",base_function,"/results_sl=",covsl,"_",base_function,".RData"))

  
