############## Calculate the weight of each resolution, i.e. the coefficients of the stage-two model, for binary classification:
coefs0=coefs = matrix(NA,5,4)
base_function = "rf"
if (base_function == "glm") {b.posthoc = B.posthoc[1:3]}
if (base_function == "qda") {b.posthoc = B.posthoc[4:6]}
if (base_function == "rf") {b.posthoc = B.posthoc[7:9]}
for (split in 1:5)
{
  TRAIN=DATA[!(DATA$subject%in%folds[[split]]),]
  tr.subj = unique(TRAIN$subject)
  n.tr.subj = length(unlist(folds)) - length(folds[[split]])
  holdout = DATA[(DATA$subject%in%folds[[split]]),]
  holdout.ind = which(DATA[,"subject"]%in%folds[[split]])
  n.holdout = nrow(holdout)
  n.whole=nrow(TRAIN)
  p_r=p=3
  
  ############## Obtain the new X for super learner:
  whole=cbind(TRAIN,Global_pred=0,Region4_pred=0,Region9_pred=0)
  TRAIN.folds = c(1:length(folds))[-split]
  for (fold in 1:length(TRAIN.folds))
  {
    #load("DATAraw.RData")
    k=folds[[TRAIN.folds[fold]]]
    #u=j=slice=k
    ind=index[-k]
    ###### Spline X for beta:
    train=TRAIN[!TRAIN$subject%in%k,]
    test=TRAIN[TRAIN$subject%in%k,]
    test.rownames=rownames(test)
    
    regional_prediction_whole=matrix(0,nrow(TRAIN),n.X)
    regional_prediction_train=matrix(0,nrow(train),n.X)
    regional_prediction_test=matrix(0,nrow(test),n.X)
    
    MRIpars_whole <- TRAIN[,MRIpars]
    MRIpars_train <- train[,MRIpars]
    MRIpars_test <- test[,MRIpars]
    
    ##### Global model:
    ###### train the candidate model on the training blocks:
    if (base_function == "rf")
    {
      global_fit <- randomForest(rf.baseformula, data = train, ntree=n.tree, mtry=2,nodesize=1, importance=T)
      global_pred <- predict(global_fit, newdata=test,type="prob")
      regional_prediction_test[,1]<- global_pred[,2]
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                  nodesize=1, importance=T)
        local_pred=predict(local_fit, newdata=test,type="prob")
        regional_prediction_test[,2]<- local_pred[,2] * (test[,"Region4"]==reg) + regional_prediction_test[,2]
      }
      
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                  nodesize=1, importance=T)
        local_pred=predict(local_fit, newdata=test,type="prob")
        regional_prediction_test[,3]<- local_pred[,2] * (test[,"Region9"]==reg) + regional_prediction_test[,3]
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    
    if (base_function == "knn")
    {
      global_fit=knn(train=train[,MRIpars],test=test[,MRIpars],cl=train$cancer,k=knn.k,prob=TRUE)
      regional_prediction_test[,1]<- 1-attr(global_fit,"prob")
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- knn(train=regional_train[,MRIpars],test=test[,MRIpars],cl=regional_train$cancer,k=knn.k,prob=TRUE)
        local_pred = 1-attr(local_fit,"prob")
        regional_prediction_test[,2] <- local_pred * (test[,"Region4"]==reg) + regional_prediction_test[,2]
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- knn(train=regional_train[,MRIpars],test=test[,MRIpars],cl=regional_train$cancer,k=knn.k,prob=TRUE)
        local_pred = 1-attr(local_fit,"prob")
        regional_prediction_test[,3] <- local_pred * (test[,"Region9"]==reg) + regional_prediction_test[,3]
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    
    if (base_function == "qda")
    {
      global_fit<-qda(train[,MRIpars], train$cancer)
      global_pred<-predict(global_fit, test[,MRIpars])
      regional_prediction_test[,1]<- global_pred$posterior[,2]
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- qda(regional_train[,MRIpars], regional_train$cancer)
        local_pred=predict(local_fit,test[,MRIpars])
        regional_prediction_test[,2]<- local_pred$posterior[,2] * (test[,"Region4"]==reg) + regional_prediction_test[,2]
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- qda(regional_train[,MRIpars], regional_train$cancer)
        local_pred=predict(local_fit,test[,MRIpars])
        regional_prediction_test[,3]<- local_pred$posterior[,2] * (test[,"Region9"]==reg) + regional_prediction_test[,3]
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    if (base_function == "glm")
    {
      global_fit <- glm(rf.baseformula, data = train, family = binomial(link = "probit")) 
      global_pred <- predict(global_fit, newdata=test,type="response")
      regional_prediction_test[,1]<- global_pred
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_train=train[train[,"Region4"]==reg,]
        local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "probit")) 
        local_pred=predict(local_fit, newdata=test,type="response")
        regional_prediction_test[,2]<- local_pred * (test[,"Region4"]==reg) + regional_prediction_test[,2]
      }
      
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_train=train[train[,"Region9"]==reg,]
        local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "probit")) 
        local_pred=predict(local_fit, newdata=test,type="response")
        regional_prediction_test[,3]<- local_pred * (test[,"Region9"]==reg) + regional_prediction_test[,3]
      }
      colnames(regional_prediction_test)=slcovnames
      whole[test.rownames,slcovnames]=regional_prediction_test
    }
    print(fold)
  }
  
  ########################## Train the super learner using the new X:
  cancer = whole[,"cancer"]  ##### notice that in each fold: should update "cancer" for the test slices
  cancer_num = as.numeric(as.character(cancer))
  
  
  ########################################################################################
  ############# Train each candidate learner on the whole training data set: #############
  ########################################################################################
  ##### Global model:
  ###### train the candidate models on the training blocks:
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
  }
  
  ##################################################################################
  ############ Evaluate the trained super learner on the holdout dataset: ############
  ##################################################################################
  ######### Obtain the new X on the holdout data set:
  holdout = cbind(holdout,Global_pred=0,Region4_pred=0,Region9_pred=0)
  if (base_function == "rf")
  {
    pred.global <- predict(baselearner.global, newdata=holdout,type="prob") #default type: response
    holdout[,slcovnames[1]]<- pred.global[,2]
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      #local4_fit <- randomForest(cancer ~ ADC + KTRANS + KEP + AUGC, data = regional_whole, ntree=n.tree, mtry=2,
      #                           nodesize=1, importance=T)
      pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="prob")
      holdout[region4.ind,slcovnames[2]]<- pred.local4[,2] 
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      #local9_fit <- randomForest(cancer ~ ADC + KTRANS + KEP + AUGC, data = regional_whole, ntree=n.tree, mtry=2,
      #                           nodesize=1, importance=T)
      pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="prob")
      holdout[region9.ind,slcovnames[3]]<- pred.local9[,2]
    }
  }
  if (base_function == "knn")
  {
    pred.global <- knn(train=TRAIN[,MRIpars],test=holdout[,MRIpars],cl=TRAIN$cancer,k=knn.k,prob=TRUE)
    holdout[,slcovnames[1]]<- 1-attr(pred.global,"prob")
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      pred.local4 <- knn(train=regional_whole[,MRIpars],test=regional_holdout[,MRIpars],cl=regional_whole$cancer,k=knn.k,prob=TRUE)
      holdout[region4.ind,slcovnames[2]]<- 1-attr(pred.local4,"prob")
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      pred.local9 <- knn(train=regional_whole[,MRIpars],test=regional_holdout[,MRIpars],cl=regional_whole$cancer,k=knn.k,prob=TRUE)
      holdout[region9.ind,slcovnames[3]]<- 1-attr(pred.local9,"prob")
    }
  }
  if (base_function == "qda")
  {
    global.fit = qda(TRAIN[,MRIpars],TRAIN$cancer)
    pred.global = predict(global.fit, holdout[,MRIpars])
    holdout[,slcovnames[1]]<- pred.global$posterior[,2]
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      local.fit4 = qda(regional_whole[,MRIpars],regional_whole$cancer)
      pred.local4 = predict(local.fit4, regional_holdout[,MRIpars])
      holdout[region4.ind,slcovnames[2]]<- pred.local4$posterior[,2]
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      local.fit9 = qda(regional_whole[,MRIpars],regional_whole$cancer)
      pred.local9 = predict(local.fit9, regional_holdout[,MRIpars])
      holdout[region9.ind,slcovnames[3]]<- pred.local9$posterior[,2]
    }
  }
  if (base_function == "glm")
  {
    global.fit <- glm(rf.baseformula, data = TRAIN, family = binomial(link = "probit")) 
    pred.global <- predict(global.fit, newdata=holdout,type="response")
    holdout[,slcovnames[1]]<- pred.global
    
    ##### Local model with 4 regions (2 by 2 split):
    for (reg in 1:resolution4)
    {
      region4.ind = which(holdout[,"Region4"] == reg)
      regional_holdout=holdout[holdout[,"Region4"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
      local.fit4 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "probit"))
      pred.local4 <- predict(local.fit4, newdata=regional_holdout,type="response")
      holdout[region4.ind,slcovnames[2]]<- pred.local4
    }
    ##### Local model with 9 regions (3 by 3 split):
    for (reg in 1:resolution9)
    {
      region9.ind = which(holdout[,"Region9"] == reg)
      regional_holdout=holdout[holdout[,"Region9"]==reg,]
      regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
      local.fit9 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "probit"))
      pred.local9 <- predict(local.fit9, newdata=regional_holdout,type="response")
      holdout[region9.ind,slcovnames[3]]<- pred.local9
    }
  }
  cancertrue = as.numeric(as.character(holdout[,"cancer"]))  ##### notice that in each fold: should update "category" for the test slices
  
  ############################################################################
  ############################ SL without smoothing ##########################
  ############################################################################
  ####### Model fitting using the new predictors:
  nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames, collapse= "+")))
  ####################### use numeric cancer in the training gave the same result:
  sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "probit")) 
  coefs0[split,]=sl.model$coefficients
  pred.sl0<-predict(sl.model, newdata=holdout,type="response")
  holdout$pred.sl0 = pred.sl0
  bpred <- prediction(predictions=pred.sl0, labels=cancertrue)
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  s90.sl0[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  s80.sl0[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
  auc.sl0[split] <- unlist(slot(bau, "y.values"))
  
  
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
  
  ############################################################################
  ############################## Non-spatial SL ##############################
  ############################################################################
  ####### Model fitting using the new predictors:
  nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames, collapse= "+")))
  ####################### use numeric cancer in the training gave the same result:
  sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "probit")) 
  coefs[split,]=sl.model$coefficients
  pred.sl<-predict(sl.model, newdata=holdout,type="response")
  holdout$pred.sl = pred.sl
  sl_prediction[holdout.ind] = pred.sl
  bpred <- prediction(predictions=pred.sl, labels=cancertrue)
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  s90.sl[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  s80.sl[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
  auc.sl[split] <- unlist(slot(bau, "y.values"))
  auc.sl[split]
  
  
  ###########################################################
  ######################## Baseline: ########################
  ###########################################################
  base_prediction = rep(NA,n.holdout)
  k=folds[[split]]
  ind=index[-k]
  
  ##### Global model:
  ###### train the candidate model on the training blocks:
  if (base_function == "knn")
  {
    baseline.model <- knn(train=TRAIN[,MRIpars],test=holdout[,MRIpars],cl=TRAIN$cancer,k=knn.k,prob=TRUE)
    pred.base <- 1-attr(baseline.model,"prob")
  }
  
  if (base_function == "rf")
  {
    baseline.model <- randomForest(rf.baseformula, data = TRAIN, ntree=n.tree, mtry=2,
                                   nodesize=1, importance=T)
    baseline.pred <- predict(baseline.model, newdata=holdout,type="prob")
    pred.base <- baseline.pred[,2]
  }
  if (base_function == "qda")
  {
    global_fit<-qda(TRAIN[,MRIpars], TRAIN$cancer)
    global_pred<-predict(global_fit, holdout[,MRIpars])
    pred.base <- global_pred$posterior[,2]
  }
  if (base_function == "spline")
  {
    global_fit<-polymars(train$category, MRIpars_train, classify = T, maxsize=t[tn,1], knots=t[tn,2])
    global_pred<-predict(global_fit, x = MRIpars_test, classify = T)
    base_prediction[test.rows]<- as.numeric(as.character(global_pred))-1
  }
  if (base_function == "glm")
  {
    baseline.model <- glm(rf.baseformula, data = TRAIN, family = binomial(link = "probit")) 
    pred.base <- predict(baseline.model, newdata=holdout,type="response")
  }
  #if (base_function == "kqr") # kernel quantile regression or quantile regression?
  
  holdout$pred.base = pred.base
  #b_prediction[holdout.ind] = pred.base
  bpred <- prediction(pred.base, cancertrue)
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  s90.base[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  s80.base[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
  auc.base[split] <- unlist(slot(bau, "y.values"))
  auc.base[split] #0.6235
  
  ###################### overall AUC:
  #bpred <- prediction(spsl_prediction, as.numeric(DATA$cancer)-1)
  #bperf<- performance(bpred,"tpr","fpr")
  #bau <- performance(bpred,"auc")
  #os90.spsl=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  #os80.spsl=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
  #oauc.spsl <- unlist(slot(bau, "y.values"))
  #oauc.spsl #0.62  
  
  
  #bpred <- prediction(sl_prediction, as.numeric(DATA$cancer)-1)
  #bperf<- performance(bpred,"tpr","fpr")
  #bau <- performance(bpred,"auc")
  #os90.sl=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  #os80.sl=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
  #oauc.sl <- unlist(slot(bau, "y.values"))
  #oauc.sl #0.62    
  
  #bpred <- prediction(b_prediction, as.numeric(DATA$cancer)-1)
  #bperf<- performance(bpred,"tpr","fpr")
  #bau <- performance(bpred,"auc")
  #os90.base=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  #os80.base=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
  #oauc.base <- unlist(slot(bau, "y.values"))
  #oauc.base #0.62  
  save(coefs0,coefs,file=paste0("results/binary/",base_function,"/coefs_",base_function,".RData"))
}

round(colMeans(coefs0),3) # coefficients when not using the intermediate spatial smoothing step
round(colMeans(coefs),3) # coefficients when using the intermediate spatial smoothing step
