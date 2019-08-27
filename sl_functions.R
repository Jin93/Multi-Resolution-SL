newloc.results = function(coord.new,coords,m,phi,nu,ws)
{
  coord.new = as.matrix(coord.new)
  ##### step 1: get the indexes of the neighboring voxels
  dists = sqrt((coord.new[1,"x"]-coords[,"x"])^2 + (coord.new[1,"y"]-coords[,"y"])^2)
  nn.voxel = order(dists)[1:m]
  coord.Ns = as.matrix(coords[nn.voxel,])
  ##### calculate distance matrix:
  #### D: length 849265: distances for N(ij)'s:
  dis.s.Ns=distmatcpp(coord.new, coord.Ns)
  dis.Ns=distmatcpp(coord.Ns, coord.Ns)
  
  ##### step 2: calculate B and F:
  BFresults = newlocBFcpp(dis.s.Ns, dis.Ns, phi, nu, 1)
  wnewhat = kronecker(BFresults$B,diag(1,n.X)) %*% matrix(t(ws[nn.voxel,]),m*n.X,1) #ws[nn.voxel,]
  return(wnewhat)
}


######## simulate binary data:
simu_binary_data = function(fillin,i,indx,sse4,sse9)
{
  temp=fillin[[i]][,c("ADC","KTRANS","KEP","AUGC","x","y","zone","cancer","Region4","Region9")]
  zone=1*(temp[,"zone"]=="PZ")+0*(temp[,"zone"]=="CG") # 1: PZ
  coords=temp[,c("x","y")]
  D <- as.matrix(dist(coords))
  #R <- exp(-phi*D)
  R <- materncovcpp(D,sigmasq_c,phi_c,nu_c)
  ##### 
  qs0 <- rmvn(1, rep(0,nrow(temp)), R) + sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
  qs <- qs0 + rnorm(n=nrow(temp),mean=0,sd=sqrt(tausq))
  while (sum(qs>0)==0) {
    qs0 <- rmvn(1, rep(0,nrow(temp)), R)+sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
    qs <- qs0 + rnorm(n=nrow(temp),mean=0,sd=sqrt(tausq))
  }
  ps <- sapply(1:nrow(temp),FUN=function(x){pnorm(q=qs[x],mean=0,sd=1)})
  
  cancer=1*(qs>0)+0*(qs<0)
  temp$cancer=cancer
  temp=cbind(temp,qs0,qs,ps)
  #B=mvrnorm(n=1,mu=rep(0,4),Sigma=Sigma)
  for (j in 1:nrow(temp))
  {
    indicator=1*((cancer[j]==1)&(zone[j]==1))+2*((cancer[j]==1)&(zone[j]==0))+3*((cancer[j]==0)&(zone[j]==1))+4*((cancer[j]==0)&(zone[j]==0))
    A=mvrnorm(n=1,mu=means[[indicator]],Sigma=covs[[indicator]])
    #temp[j,c("ADC","KTRANS","KEP","AUGC")]=A + B + sse4[[temp[j,"Region4"]]] + sse9[[temp[j,"Region9"]]]
    temp[j,c("ADC","KTRANS","KEP","AUGC")]=A + sse4[[temp[j,"Region4"]]] + sse9[[temp[j,"Region9"]]]
  }
  temp$subject = indx
  #fillnew[[new]]=temp
  #print(new)
  #new=new+1
  #print(i)
  return(temp)
}

######## simulate categorical data:
simu_category_data = function(fillin,i,indx)
{
  temp=fillin[[i]][,c("ADC","KTRANS","KEP","AUGC","x","y","zone","cancer","Region4","Region9")]
  zone=1*(temp[,"zone"]=="PZ")+0*(temp[,"zone"]=="CG") # 1: PZ
  coords=temp[,c("x","y")]
  D <- as.matrix(dist(coords))
  R <- materncovcpp(D,sigmasq_c,phi_c,nu_c)
  ##### 
  qs0 <- rmvn(1, rep(0,nrow(temp)), R) + sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
  qs <- qs0 + rnorm(n=nrow(temp),mean=0,sd=sqrt(tausq))
  while (sum(qs>0)==0) {
    qs0 <- rmvn(1, rep(0,nrow(temp)), R)+sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
    qs <- qs0 + rnorm(n=nrow(temp),mean=0,sd=sqrt(tausq))
  }
  ps <- sapply(1:nrow(temp),FUN=function(x){pnorm(q=qs[x],mean=0,sd=1)})
  
  temp=cbind(temp,qs0,qs,ps)
  temp$subject = indx
  #qs.total = c(qs.total,qs)
  return(list(temp=temp,qs = qs))
  print(indx)
}

simu_category_data.addition = function(fillnewi,sse4,sse9,qs.t1,qs.t2)
{
  temp = fillnewi
  zone=1*(temp[,"zone"]=="PZ")+0*(temp[,"zone"]=="CG") # 1: PZ
  cat.temp = 0*(temp$qs<qs.t1) + 1*((temp$qs>qs.t1)&(temp$qs<qs.t2)) + 2*(temp$qs>qs.t2)
  temp$category=cat.temp
  #B=mvrnorm(n=1,mu=rep(0,4),Sigma=Sigma)
  for (j in 1:nrow(temp))
  {
    indicator=1*((cat.temp[j]==0)&(zone[j]==1)) + 2*((cat.temp[j]==0)&(zone[j]==0)) + 3*((cat.temp[j]==1)&(zone[j]==1)) + 4*((cat.temp[j]==1)&(zone[j]==0)) + 5*((cat.temp[j]==2)&(zone[j]==1)) + 6*((cat.temp[j]==2)&(zone[j]==0))
    A=mvrnorm(n=1,mu=means[[indicator]],Sigma=covs[[indicator]])
    temp[j,c("ADC","KTRANS","KEP","AUGC")]= A + sse4[[temp[j,"Region4"]]] + sse9[[temp[j,"Region9"]]]
  }
  #fillnew[[i]]=temp
  return(temp)
}



############################################################################
####### Function for Tuning Smoothness Parameter for Binary Outcome: #######
############################################################################
##### binary
b.posthoc.tuning = function(t,base_function)
{
  if (base_function != "combine")
  {
    slcovnames = covariates = c("Global_pred","Region4_pred","Region9_pred")
  }
  if (base_function == "combine")
  {
    slcovnames = covariates = c(sapply(1:length(base_functions),FUN=function(x){c(paste0("Global_pred_",base_functions[x]),paste0("Region4_pred_",base_functions[x]),paste0("Region9_pred_",base_functions[x]))}))
  }
  b.posthoc = B[t]
  auc.sl=s90.sl=s80.sl=numeric()
  auc.sl0=s90.sl0=s80.sl0=numeric()
  auc.base=s90.base=s80.base=numeric()
  for (split in c(1:length(folds)))
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
    if (base_function != "combine")
    {
      whole=cbind(TRAIN,Global_pred=0,Region4_pred=0,Region9_pred=0)
    }
    if (base_function == "combine")
    {
      whole=cbind(TRAIN,Global_pred_rf=0,Region4_pred_rf=0,Region9_pred_rf=0,Global_pred_knn=0,Region4_pred_knn=0,Region9_pred_knn=0,Global_pred_qda=0,Region4_pred_qda=0,Region9_pred_qda=0)
    }
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
        global_fit <- glm(rf.baseformula, data = train, family = binomial(link = "logit")) 
        global_pred <- predict(global_fit, newdata=test,type="response")
        regional_prediction_test[,1]<- global_pred
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_train=train[train[,"Region4"]==reg,]
          local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "logit")) 
          local_pred=predict(local_fit, newdata=test,type="response")
          regional_prediction_test[,2]<- local_pred * (test[,"Region4"]==reg) + regional_prediction_test[,2]
        }
        
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_train=train[train[,"Region9"]==reg,]
          local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "logit")) 
          local_pred=predict(local_fit, newdata=test,type="response")
          regional_prediction_test[,3]<- local_pred * (test[,"Region9"]==reg) + regional_prediction_test[,3]
        }
        colnames(regional_prediction_test)=slcovnames
        whole[test.rownames,slcovnames]=regional_prediction_test
      }
      if (base_function == "combine")
      {
        regional_prediction_test=matrix(0,nrow(test),n.X*length(base_functions))
        
        #####################
        ######## GLM: #######
        #####################
        global_fit=glm(rf.baseformula, data = train, family = binomial(link = "logit")) 
        regional_prediction_test[,1]<- predict(global_fit, newdata=test,type="response")
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_train=train[train[,"Region4"]==reg,]
          local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "logit")) 
          local_pred = predict(local_fit, newdata=test,type="response")
          regional_prediction_test[,2] <- local_pred * (test[,"Region4"]==reg) + regional_prediction_test[,2]
        }
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_train=train[train[,"Region9"]==reg,]
          local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "logit")) 
          local_pred = predict(local_fit, newdata=test,type="response")
          regional_prediction_test[,3] <- local_pred * (test[,"Region9"]==reg) + regional_prediction_test[,3]
        }
        
        #####################
        ######## QDA: #######
        #####################        
        global_fit<-qda(train[,MRIpars], train$cancer)
        global_pred<-predict(global_fit, test[,MRIpars])
        regional_prediction_test[,4]<- global_pred$posterior[,2]
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_train=train[train[,"Region4"]==reg,]
          local_fit <- qda(regional_train[,MRIpars], regional_train$cancer)
          local_pred=predict(local_fit,test[,MRIpars])
          regional_prediction_test[,5]<- local_pred$posterior[,2] * (test[,"Region4"]==reg) + regional_prediction_test[,5]
        }
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_train=train[train[,"Region9"]==reg,]
          local_fit <- qda(regional_train[,MRIpars], regional_train$cancer)
          local_pred=predict(local_fit,test[,MRIpars])
          regional_prediction_test[,6]<- local_pred$posterior[,2] * (test[,"Region9"]==reg) + regional_prediction_test[,6]
        }
        
        #####################
        ######## RF：########
        #####################
        global_fit <- randomForest(rf.baseformula, data = train, ntree=n.tree, mtry=2,nodesize=1, importance=T)
        global_pred <- predict(global_fit, newdata=test,type="prob")
        regional_prediction_test[,7]<- global_pred[,2]
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_train=train[train[,"Region4"]==reg,]
          local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                    nodesize=1, importance=T)
          local_pred=predict(local_fit, newdata=test,type="prob")
          regional_prediction_test[,8]<- local_pred[,2] * (test[,"Region4"]==reg) + regional_prediction_test[,8]
        }
        
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_train=train[train[,"Region9"]==reg,]
          local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                    nodesize=1, importance=T)
          local_pred=predict(local_fit, newdata=test,type="prob")
          regional_prediction_test[,9]<- local_pred[,2] * (test[,"Region9"]==reg) + regional_prediction_test[,9]
        }
        
        colnames(regional_prediction_test)=slcovnames
        whole[test.rownames,slcovnames]=regional_prediction_test
      }
    }
    
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
    ############ Evaluate the trained super learner on the holdout dataset: ##########
    ##################################################################################
    ######### Obtain the new X on the holdout data set:
    if (base_function != "combine")
    {
      holdout = cbind(holdout,Global_pred=0,Region4_pred=0,Region9_pred=0)
    }
    if (base_function == "combine")
    {
      holdout = cbind(holdout,Global_pred_glm=0,Region4_pred_glm=0,Region9_pred_glm=0,Global_pred_qda=0,Region4_pred_qda=0,Region9_pred_qda=0,Global_pred_rf=0,Region4_pred_rf=0,Region9_pred_rf=0)
    }
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
      global.fit <- glm(rf.baseformula, data = TRAIN, family = binomial(link = "logit")) 
      pred.global <- predict(global.fit, newdata=holdout,type="response")
      holdout[,slcovnames[1]]<- pred.global
      
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        region4.ind = which(holdout[,"Region4"] == reg)
        regional_holdout=holdout[holdout[,"Region4"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
        local.fit4 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "logit"))
        pred.local4 <- predict(local.fit4, newdata=regional_holdout,type="response")
        holdout[region4.ind,slcovnames[2]]<- pred.local4
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        region9.ind = which(holdout[,"Region9"] == reg)
        regional_holdout=holdout[holdout[,"Region9"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
        local.fit9 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "logit"))
        pred.local9 <- predict(local.fit9, newdata=regional_holdout,type="response")
        holdout[region9.ind,slcovnames[3]]<- pred.local9
      }
    }
    if (base_function == "combine")
    {
      ##########################
      ########### GLM ##########
      ##########################
      global.fit <- glm(rf.baseformula, data = TRAIN, family = binomial(link = "logit")) 
      holdout[,slcovnames[1]]<- predict(global.fit, newdata=holdout,type="response")
      
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        region4.ind = which(holdout[,"Region4"] == reg)
        regional_holdout=holdout[holdout[,"Region4"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
        local.fit4 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "logit"))
        pred.local4 <- predict(local.fit4, newdata=regional_holdout,type="response")
        holdout[region4.ind,slcovnames[2]]<- pred.local4
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        region9.ind = which(holdout[,"Region9"] == reg)
        regional_holdout=holdout[holdout[,"Region9"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
        local.fit9 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "logit"))
        pred.local9 <- predict(local.fit9, newdata=regional_holdout,type="response")
        holdout[region9.ind,slcovnames[3]]<- pred.local9
      }
      ##########################
      ########### QDA ##########
      ##########################
      global.fit = qda(TRAIN[,MRIpars],TRAIN$cancer)
      pred.global = predict(global.fit, holdout[,MRIpars])
      holdout[,slcovnames[4]]<- pred.global$posterior[,2]
      
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        region4.ind = which(holdout[,"Region4"] == reg)
        regional_holdout=holdout[holdout[,"Region4"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
        local.fit4 = qda(regional_whole[,MRIpars],regional_whole$cancer)
        pred.local4 = predict(local.fit4, regional_holdout[,MRIpars])
        holdout[region4.ind,slcovnames[5]]<- pred.local4$posterior[,2]
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        region9.ind = which(holdout[,"Region9"] == reg)
        regional_holdout=holdout[holdout[,"Region9"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
        local.fit9 = qda(regional_whole[,MRIpars],regional_whole$cancer)
        pred.local9 = predict(local.fit9, regional_holdout[,MRIpars])
        holdout[region9.ind,slcovnames[6]]<- pred.local9$posterior[,2]
      }
      
      ##########################
      ########### RF ###########
      ##########################
      ###### train the candidate models on the training blocks:
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
      ############## Prediction:
      pred.global <- predict(baselearner.global, newdata=holdout,type="prob") #default type: response
      holdout[,slcovnames[7]]<- pred.global[,2]
      
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        region4.ind = which(holdout[,"Region4"] == reg)
        regional_holdout=holdout[holdout[,"Region4"]==reg,]
        #local4_fit <- randomForest(cancer ~ ADC + KTRANS + KEP + AUGC, data = regional_whole, ntree=n.tree, mtry=2,
        #                           nodesize=1, importance=T)
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="prob")
        holdout[region4.ind,slcovnames[8]]<- pred.local4[,2] 
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        region9.ind = which(holdout[,"Region9"] == reg)
        regional_holdout=holdout[holdout[,"Region9"]==reg,]
        #local9_fit <- randomForest(cancer ~ ADC + KTRANS + KEP + AUGC, data = regional_whole, ntree=n.tree, mtry=2,
        #                           nodesize=1, importance=T)
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="prob")
        holdout[region9.ind,slcovnames[9]]<- pred.local9[,2]
      }
    }
    cancertrue = as.numeric(as.character(holdout[,"cancer"]))  ##### notice that in each fold: should update "category" for the test slices
    
    ################################################################
    ########################## Smoothing: ##########################
    ################################################################
    wholesmooth = whole
    for (ii in tr.subj)
    {
      tem = whole[whole$subject == ii,]
      row.ind = which(whole$subject == ii)
      Xcoord=tem[,"x"];Ycoord=tem[,"y"];
      for (jj in 1:length(slcovnames))
      {
        ######## smooth pij:
        X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[jj]])
        #b=0.2
        #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
        Spre=as.numeric(Smooth(X, sigma=b.posthoc,at="points",edge=TRUE, diggle=F))
        tem[,covariates[jj]] = Spre
        wholesmooth[row.ind,covariates[jj]] = Spre
        ######## somooth qij:
      }
    }
    whole=wholesmooth
    holdoutsmooth = holdout
    for (ii in folds[[split]])
    {
      tem = holdoutsmooth[holdoutsmooth$subject == ii,]
      row.ind = which(holdoutsmooth$subject == ii)
      Xcoord=round(tem[,"x"],2);Ycoord=round(tem[,"y"],2);
      for (jj in 1:length(slcovnames))
      {
        ######## smooth pij:
        X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[jj]])
        #b=0.2
        #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
        Spre=as.numeric(Smooth(X, sigma=b.posthoc,at="points",edge=TRUE, diggle=F))
        tem[,covariates[jj]] = Spre
        holdoutsmooth[row.ind,covariates[jj]] = Spre
        ######## somooth qij:
      }
    }
    holdout = holdoutsmooth
    
    ############################################################################
    ############################ SL with smoothing #############################
    ############################################################################
    ####### Model fitting using the new predictors:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames, collapse= "+")))
    ####################### use numeric cancer in the training gave the same result:
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    auc.sl[split] <- unlist(slot(bau, "y.values"))
    print(paste0("Bandwidth ",t," - Split ",split,": finished"))
  }
  cvauc=mean(auc.sl)
  print(paste0("Bandwidth ",t,": finished "))
  return(cvauc)
}




############################################################################
#### Function for Tuning Smoothness Parameter for categorical outcome: #####
############################################################################
##### criteria: error rate using equal weights
b.posthoc.tuning.category = function(t)
{
  b.posthoc = B[t]
  err.sl=numeric()
  pred.mat=matrix(NA,nrow(DATA),1)
  colnames(pred.mat) = "pred.sl"
  for (split in c(1:length(folds)))
  {
    TRAIN=DATA[!(DATA$subject%in%folds[[split]]),]
    tr.subj = unique(TRAIN$subject)
    n.tr.subj = length(unlist(folds)) - length(folds[[split]])
    holdout = DATA[(DATA$subject%in%folds[[split]]),]
    holdout.ind = which(DATA[,"subject"]%in%folds[[split]])
    n.holdout = nrow(holdout)
    n.whole=nrow(TRAIN)
    p_r=p=3
    
    ########### Weights for GLM:
    weights.train=list()
    weights.train[[1]] = rep(1,n.whole)
    
    ############## Obtain the new X for super learner:
    whole=cbind(TRAIN,Global_pred1=0,Global_pred2=0,Region4_pred1=0,Region4_pred2=0,Region9_pred1=0,Region9_pred2=0)
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
      #print(fold)
    }
    
    ########################################################################################
    ############# Train each candidate learner on the whole training data set: #############
    ########################################################################################
    ##### Global model:
    ###### train the candidate models on the training blocks:
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
    ################################################################
    ########################## Smoothing: ##########################
    ################################################################
    wholesmooth = whole
    for (ii in tr.subj)
    {
      tem = whole[whole$subject == ii,]
      row.ind = which(whole$subject == ii)
      Xcoord=tem[,"x"];Ycoord=tem[,"y"];
      for (jj in 1:length(slcovnames))
      {
        ######## smooth pij:
        X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[jj]])
        #b=0.2
        #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
        Spre=as.numeric(Smooth(X, sigma=b.posthoc,at="points",edge=TRUE, diggle=F))
        tem[,covariates[jj]] = Spre
        wholesmooth[row.ind,covariates[jj]] = Spre
        ######## somooth qij:
      }
    }
    whole=wholesmooth
    holdoutsmooth = holdout
    for (ii in folds[[split]])
    {
      tem = holdoutsmooth[holdoutsmooth$subject == ii,]
      row.ind = which(holdoutsmooth$subject == ii)
      Xcoord=round(tem[,"x"],2);Ycoord=round(tem[,"y"],2);
      for (jj in 1:length(slcovnames))
      {
        ######## smooth pij:
        X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[jj]])
        #b=0.2
        #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
        Spre=as.numeric(Smooth(X, sigma=b.posthoc,at="points",edge=TRUE, diggle=F))
        tem[,covariates[jj]] = Spre
        holdoutsmooth[row.ind,covariates[jj]] = Spre
        ######## somooth qij:
      }
    }
    holdout = holdoutsmooth
    
    ############################################################################
    ############################ SL with smoothing #############################
    ############################################################################
    ####### Model fitting using the new predictors:
    nospatial_polr_fomula = as.formula(paste("category ~ ", paste(slcovnames, collapse= "+")))
    
    sl.model <- polr(nospatial_polr_fomula, data = whole, method = "probit",Hess=TRUE,weights=weights.train[[1]]) 
    ####### Prediction:
    prediction.superlearner<-predict(sl.model, holdout,type="class")
    pred.sl<- as.numeric(as.character(prediction.superlearner))
    holdout = data.frame(holdout, pred.sl = pred.sl)
    
    pred.mat[holdout.ind,"pred.sl"] = holdout[,"pred.sl"]
    
    print(paste0("Bandwidth ",t," - Split ",split,": finished"))
  }
  err = 1 - sum(DATA$category == pred.mat[,"pred.sl"])/nrow(DATA)
  print(paste0("Bandwidth ",t,": finished "))
  return(err)
}




############################################################################
################# Function for Binary Simulation per Methd: ################
############################################################################
simu_binary = function(DATA,base_function,b.posthoc)
{
  if (base_function != "combine")
  {
    slcovnames = covariates = c("Global_pred","Region4_pred","Region9_pred")
  }
  if (base_function == "combine")
  {
    slcovnames = covariates = c(sapply(1:length(base_functions),FUN=function(x){c(paste0("Global_pred_",base_functions[x]),paste0("Region4_pred_",base_functions[x]),paste0("Region9_pred_",base_functions[x]))}))
  }
  
  auc.sl.glm=s90.sl.glm=s80.sl.glm=numeric()
  auc.sl.qda=s90.sl.qda=s80.sl.qda=numeric()
  auc.sl.rf=s90.sl.rf=s80.sl.rf=numeric()
  auc.sl.combine=s90.sl.combine=s80.sl.combine=numeric()
  auc.sl0.glm=s90.sl0.glm=s80.sl0.glm=numeric()
  auc.sl0.qda=s90.sl0.qda=s80.sl0.qda=numeric()
  auc.sl0.rf=s90.sl0.rf=s80.sl0.rf=numeric()
  auc.sl0.combine=s90.sl0.combine=s80.sl0.combine=numeric()
  auc.base.glm=s90.base.glm=s80.base.glm=numeric()
  auc.base.qda=s90.base.qda=s80.base.qda=numeric()
  auc.base.rf=s90.base.rf=s80.base.rf=numeric()
  ################################################################################
  ################################################################################
  sl_prediction=rep(NA,nrow(DATA))
  for (split in c(1:length(folds)))
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
    if (base_function != "combine")
    {
      whole=cbind(TRAIN,Global_pred=0,Region4_pred=0,Region9_pred=0)
    }
    if (base_function == "combine")
    {
      whole=cbind(TRAIN,Global_pred_glm=0,Region4_pred_glm=0,Region9_pred_glm=0,Global_pred_qda=0,Region4_pred_qda=0,Region9_pred_qda=0,Global_pred_rf=0,Region4_pred_rf=0,Region9_pred_rf=0)
    }
    TRAIN.folds = c(1:length(folds))[-split]
    for (fold in 1:length(TRAIN.folds))
    {
      k=folds[[TRAIN.folds[fold]]]
      ind=index[-k]
      train=TRAIN[!TRAIN$subject%in%k,]
      test=TRAIN[TRAIN$subject%in%k,]
      test.rownames=rownames(test)
      
      MRIpars_whole <- TRAIN[,MRIpars]
      MRIpars_train <- train[,MRIpars]
      MRIpars_test <- test[,MRIpars]
      
      ##### Global model:
      ###### train the candidate model on the training blocks:
      if (base_function == "rf")
      {
        regional_prediction_test=matrix(0,nrow(test),n.X)
        
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
        regional_prediction_test=matrix(0,nrow(test),n.X)
        
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
      
      if (base_function == "glm")
      {
        regional_prediction_test=matrix(0,nrow(test),n.X)
        
        global_fit <- glm(rf.baseformula, data = train, family = binomial(link = "logit")) 
        global_pred <- predict(global_fit, newdata=test,type="response")
        regional_prediction_test[,1]<- global_pred
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_train=train[train[,"Region4"]==reg,]
          local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "logit")) 
          local_pred=predict(local_fit, newdata=test,type="response")
          regional_prediction_test[,2]<- local_pred * (test[,"Region4"]==reg) + regional_prediction_test[,2]
        }
        
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_train=train[train[,"Region9"]==reg,]
          local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "logit")) 
          local_pred=predict(local_fit, newdata=test,type="response")
          regional_prediction_test[,3]<- local_pred * (test[,"Region9"]==reg) + regional_prediction_test[,3]
        }
        colnames(regional_prediction_test)=slcovnames
        whole[test.rownames,slcovnames]=regional_prediction_test
      }
      
      if (base_function == "qda")
      {
        regional_prediction_test=matrix(0,nrow(test),n.X)
        
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
      
      if (base_function == "combine")
      {
        regional_prediction_test=matrix(0,nrow(test),n.X*length(base_functions))
        
        #####################
        ######## GLM: #######
        #####################
        global_fit=glm(rf.baseformula, data = train, family = binomial(link = "logit")) 
        regional_prediction_test[,1]<- predict(global_fit, newdata=test,type="response")
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_train=train[train[,"Region4"]==reg,]
          local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "logit")) 
          local_pred = predict(local_fit, newdata=test,type="response")
          regional_prediction_test[,2] <- local_pred * (test[,"Region4"]==reg) + regional_prediction_test[,2]
        }
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_train=train[train[,"Region9"]==reg,]
          local_fit <- glm(rf.baseformula, data = regional_train, family = binomial(link = "logit")) 
          local_pred = predict(local_fit, newdata=test,type="response")
          regional_prediction_test[,3] <- local_pred * (test[,"Region9"]==reg) + regional_prediction_test[,3]
        }
        
        #####################
        ######## QDA: #######
        #####################        
        global_fit<-qda(train[,MRIpars], train$cancer)
        global_pred<-predict(global_fit, test[,MRIpars])
        regional_prediction_test[,4]<- global_pred$posterior[,2]
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_train=train[train[,"Region4"]==reg,]
          local_fit <- qda(regional_train[,MRIpars], regional_train$cancer)
          local_pred=predict(local_fit,test[,MRIpars])
          regional_prediction_test[,5]<- local_pred$posterior[,2] * (test[,"Region4"]==reg) + regional_prediction_test[,5]
        }
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_train=train[train[,"Region9"]==reg,]
          local_fit <- qda(regional_train[,MRIpars], regional_train$cancer)
          local_pred=predict(local_fit,test[,MRIpars])
          regional_prediction_test[,6]<- local_pred$posterior[,2] * (test[,"Region9"]==reg) + regional_prediction_test[,6]
        }
        
        #####################
        ######## RF：########
        #####################
        global_fit <- randomForest(rf.baseformula, data = train, ntree=n.tree, mtry=2,nodesize=1, importance=T)
        global_pred <- predict(global_fit, newdata=test,type="prob")
        regional_prediction_test[,7]<- global_pred[,2]
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_train=train[train[,"Region4"]==reg,]
          local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                    nodesize=1, importance=T)
          local_pred=predict(local_fit, newdata=test,type="prob")
          regional_prediction_test[,8]<- local_pred[,2] * (test[,"Region4"]==reg) + regional_prediction_test[,8]
        }
        
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_train=train[train[,"Region9"]==reg,]
          local_fit <- randomForest(rf.baseformula, data = regional_train, ntree=n.tree, mtry=2,
                                    nodesize=1, importance=T)
          local_pred=predict(local_fit, newdata=test,type="prob")
          regional_prediction_test[,9]<- local_pred[,2] * (test[,"Region9"]==reg) + regional_prediction_test[,9]
        }
        
        colnames(regional_prediction_test)=slcovnames
        whole[test.rownames,slcovnames]=regional_prediction_test
      }
    }
    
    ########################## Train the super learner using the new X:
    cancer = whole[,"cancer"]  ##### notice that in each fold: should update "cancer" for the test slices
    cancer_num = as.numeric(as.character(cancer))
    
    ########################################################################################
    ############# Train each candidate learner on the whole training data set: #############
    ########################################################################################
    ##### Global model:
    
    
    ##################################################################################
    ############ Evaluate the trained super learner on the holdout dataset: ############
    ##################################################################################
    ######### Obtain the new X on the holdout data set:
    if (base_function != "combine")
    {
      holdout = cbind(holdout,Global_pred=0,Region4_pred=0,Region9_pred=0)
    }
    if (base_function == "combine")
    {
      holdout = cbind(holdout,Global_pred_glm=0,Region4_pred_glm=0,Region9_pred_glm=0,Global_pred_qda=0,Region4_pred_qda=0,Region9_pred_qda=0,Global_pred_rf=0,Region4_pred_rf=0,Region9_pred_rf=0)
    }
    
    if (base_function == "rf")
    {
      ###### train the candidate models on the training blocks:
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
      ############## Prediction:
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
      global.fit <- glm(rf.baseformula, data = TRAIN, family = binomial(link = "logit")) 
      pred.global <- predict(global.fit, newdata=holdout,type="response")
      holdout[,slcovnames[1]]<- pred.global
      
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        region4.ind = which(holdout[,"Region4"] == reg)
        regional_holdout=holdout[holdout[,"Region4"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
        local.fit4 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "logit"))
        pred.local4 <- predict(local.fit4, newdata=regional_holdout,type="response")
        holdout[region4.ind,slcovnames[2]]<- pred.local4
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        region9.ind = which(holdout[,"Region9"] == reg)
        regional_holdout=holdout[holdout[,"Region9"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
        local.fit9 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "logit"))
        pred.local9 <- predict(local.fit9, newdata=regional_holdout,type="response")
        holdout[region9.ind,slcovnames[3]]<- pred.local9
      }
    }
    
    if (base_function == "combine")
    {
      ##########################
      ########### GLM ##########
      ##########################
      global.fit <- glm(rf.baseformula, data = TRAIN, family = binomial(link = "logit")) 
      holdout[,slcovnames[1]]<- predict(global.fit, newdata=holdout,type="response")
      
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        region4.ind = which(holdout[,"Region4"] == reg)
        regional_holdout=holdout[holdout[,"Region4"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
        local.fit4 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "logit"))
        pred.local4 <- predict(local.fit4, newdata=regional_holdout,type="response")
        holdout[region4.ind,slcovnames[2]]<- pred.local4
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        region9.ind = which(holdout[,"Region9"] == reg)
        regional_holdout=holdout[holdout[,"Region9"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
        local.fit9 <- glm(rf.baseformula, data = regional_whole, family = binomial(link = "logit"))
        pred.local9 <- predict(local.fit9, newdata=regional_holdout,type="response")
        holdout[region9.ind,slcovnames[3]]<- pred.local9
      }
      ##########################
      ########### QDA ##########
      ##########################
      global.fit = qda(TRAIN[,MRIpars],TRAIN$cancer)
      pred.global = predict(global.fit, holdout[,MRIpars])
      holdout[,slcovnames[4]]<- pred.global$posterior[,2]
      
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        region4.ind = which(holdout[,"Region4"] == reg)
        regional_holdout=holdout[holdout[,"Region4"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region4"]==reg,]
        local.fit4 = qda(regional_whole[,MRIpars],regional_whole$cancer)
        pred.local4 = predict(local.fit4, regional_holdout[,MRIpars])
        holdout[region4.ind,slcovnames[5]]<- pred.local4$posterior[,2]
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        region9.ind = which(holdout[,"Region9"] == reg)
        regional_holdout=holdout[holdout[,"Region9"]==reg,]
        regional_whole=TRAIN[TRAIN[,"Region9"]==reg,]
        local.fit9 = qda(regional_whole[,MRIpars],regional_whole$cancer)
        pred.local9 = predict(local.fit9, regional_holdout[,MRIpars])
        holdout[region9.ind,slcovnames[6]]<- pred.local9$posterior[,2]
      }
      
      ##########################
      ########### RF ###########
      ##########################
      ###### train the candidate models on the training blocks:
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
      ############## Prediction:
      pred.global <- predict(baselearner.global, newdata=holdout,type="prob") #default type: response
      holdout[,slcovnames[7]]<- pred.global[,2]
      
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        region4.ind = which(holdout[,"Region4"] == reg)
        regional_holdout=holdout[holdout[,"Region4"]==reg,]
        #local4_fit <- randomForest(cancer ~ ADC + KTRANS + KEP + AUGC, data = regional_whole, ntree=n.tree, mtry=2,
        #                           nodesize=1, importance=T)
        pred.local4 <- predict(baselearner.local4[[reg]], newdata=regional_holdout,type="prob")
        holdout[region4.ind,slcovnames[8]]<- pred.local4[,2] 
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        region9.ind = which(holdout[,"Region9"] == reg)
        regional_holdout=holdout[holdout[,"Region9"]==reg,]
        #local9_fit <- randomForest(cancer ~ ADC + KTRANS + KEP + AUGC, data = regional_whole, ntree=n.tree, mtry=2,
        #                           nodesize=1, importance=T)
        pred.local9 <- predict(baselearner.local9[[reg]], newdata=regional_holdout,type="prob")
        holdout[region9.ind,slcovnames[9]]<- pred.local9[,2]
      }
    }
    cancertrue = as.numeric(as.character(holdout[,"cancer"]))  ##### notice that in each fold: should update "category" for the test slices
    
    ############################################################################
    ############################ SL without smoothing ##########################
    ############################################################################
    ####### GLM:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames[1:3], collapse= "+")))
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    s90.sl0.glm[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    s80.sl0.glm[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
    auc.sl0.glm[split] <- unlist(slot(bau, "y.values"))
    ####### QDA:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames[4:6], collapse= "+")))
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    s90.sl0.qda[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    s80.sl0.qda[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
    auc.sl0.qda[split] <- unlist(slot(bau, "y.values"))
    ####### RF:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames[7:9], collapse= "+")))
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    s90.sl0.rf[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    s80.sl0.rf[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
    auc.sl0.rf[split] <- unlist(slot(bau, "y.values"))
    ####### combine:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames, collapse= "+")))
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    s90.sl0.combine[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    s80.sl0.combine[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
    auc.sl0.combine[split] <- unlist(slot(bau, "y.values"))   
    
    ################################################################
    ########################## Smoothing: ##########################
    ################################################################
    wholesmooth = whole
    for (ii in tr.subj)
    {
      tem = whole[whole$subject == ii,]
      row.ind = which(whole$subject == ii)
      Xcoord=tem[,"x"];Ycoord=tem[,"y"];
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
      Xcoord=round(tem[,"x"],2);Ycoord=round(tem[,"y"],2);
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
    ############################ SL with smoothing #############################
    ############################################################################
    ####### GLM:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames[1:3], collapse= "+")))
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    s90.sl.glm[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    s80.sl.glm[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
    auc.sl.glm[split] <- unlist(slot(bau, "y.values"))
    ####### QDA:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames[4:6], collapse= "+")))
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    s90.sl.qda[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    s80.sl.qda[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
    auc.sl.qda[split] <- unlist(slot(bau, "y.values"))
    ####### RF:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames[7:9], collapse= "+")))
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    s90.sl.rf[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    s80.sl.rf[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
    auc.sl.rf[split] <- unlist(slot(bau, "y.values"))
    ####### combine:
    nospatial_polr_fomula = as.formula(paste("cancer ~ ", paste(slcovnames, collapse= "+")))
    sl.model <- glm(nospatial_polr_fomula, data = whole, family = binomial(link = "logit")) 
    pred.sl<-predict(sl.model, newdata=holdout,type="response")
    holdout$pred.sl = pred.sl
    bpred <- prediction(predictions=pred.sl, labels=cancertrue)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    s90.sl.combine[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    s80.sl.combine[split]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
    auc.sl.combine[split] <- unlist(slot(bau, "y.values"))    
    
    
    ###########################################################
    ######################## Baseline: ########################
    ###########################################################
    #if (base_function != "combine")
    #{
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
      ### rf 
        baseline.model <- randomForest(rf.baseformula, data = TRAIN, ntree=n.tree, mtry=2,
                                       nodesize=1, importance=T)
        baseline.pred <- predict(baseline.model, newdata=holdout,type="prob")
        pred.base <- baseline.pred[,2]
        bpred <- prediction(pred.base, cancertrue)
        bperf<- performance(bpred,"tpr","fpr")
        bau <- performance(bpred,"auc")
        s90.base.rf = max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
        s80.base.rf = max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
        auc.base.rf = unlist(slot(bau, "y.values"))
       #### qda:
        global_fit<-qda(TRAIN[,MRIpars], TRAIN$cancer)
        global_pred<-predict(global_fit, holdout[,MRIpars])
        pred.base <- global_pred$posterior[,2]
        bpred <- prediction(pred.base, cancertrue)
        bperf<- performance(bpred,"tpr","fpr")
        bau <- performance(bpred,"auc")
        s90.base.qda = max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
        s80.base.qda = max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
        auc.base.qda = unlist(slot(bau, "y.values"))
        
      if (base_function == "spline")
      {
        global_fit<-polymars(train$category, MRIpars_train, classify = T, maxsize=t[tn,1], knots=t[tn,2])
        global_pred<-predict(global_fit, x = MRIpars_test, classify = T)
        base_prediction[test.rows]<- as.numeric(as.character(global_pred))-1
      }
    ### glm:
        baseline.model <- glm(rf.baseformula, data = TRAIN, family = binomial(link = "logit")) 
        pred.base <- predict(baseline.model, newdata=holdout,type="response")
        bpred <- prediction(pred.base, cancertrue)
        bperf<- performance(bpred,"tpr","fpr")
        bau <- performance(bpred,"auc")
        s90.base.glm = max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
        s80.base.glm = max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.2])
        auc.base.glm = unlist(slot(bau, "y.values"))
    #}
  }
  AUC.base.glm = mean(auc.base.glm)
  S80.base.glm = mean(s80.base.glm)
  S90.base.glm = mean(s90.base.glm)
  AUC.sl0.glm = mean(auc.sl0.glm)
  S80.sl0.glm = mean(s80.sl0.glm)
  S90.sl0.glm = mean(s90.sl0.glm)
  AUC.sl.glm = mean(auc.sl.glm)
  S80.sl.glm = mean(s80.sl.glm)
  S90.sl.glm = mean(s90.sl.glm)
  AUC.base.qda = mean(auc.base.qda)
  S80.base.qda = mean(s80.base.qda)
  S90.base.qda = mean(s90.base.qda)
  AUC.sl0.qda = mean(auc.sl0.qda)
  S80.sl0.qda = mean(s80.sl0.qda)
  S90.sl0.qda = mean(s90.sl0.qda)
  AUC.sl.qda = mean(auc.sl.qda)
  S80.sl.qda = mean(s80.sl.qda)
  S90.sl.qda = mean(s90.sl.qda)  
  AUC.base.rf = mean(auc.base.rf)
  S80.base.rf = mean(s80.base.rf)
  S90.base.rf = mean(s90.base.rf)
  AUC.sl0.rf = mean(auc.sl0.rf)
  S80.sl0.rf = mean(s80.sl0.rf)
  S90.sl0.rf = mean(s90.sl0.rf)
  AUC.sl.rf = mean(auc.sl.rf)
  S80.sl.rf = mean(s80.sl.rf)
  S90.sl.rf = mean(s90.sl.rf) 
  AUC.sl0.combine = mean(auc.sl0.combine)
  S80.sl0.combine = mean(s80.sl0.combine)
  S90.sl0.combine = mean(s90.sl0.combine)
  AUC.sl.combine = mean(auc.sl.combine)
  S80.sl.combine = mean(s80.sl.combine)
  S90.sl.combine = mean(s90.sl.combine)
  #if (base_function != "combine")
  #{
    outputs = c(AUC.base.glm,S80.base.glm,S90.base.glm,AUC.sl0.glm,S80.sl0.glm,S90.sl0.glm,AUC.sl.glm,S80.sl.glm,S90.sl.glm,
                AUC.base.qda,S80.base.qda,S90.base.qda,AUC.sl0.qda,S80.sl0.qda,S90.sl0.qda,AUC.sl.qda,S80.sl.qda,S90.sl.qda,
                AUC.base.rf,S80.base.rf,S90.base.rf,AUC.sl0.rf,S80.sl0.rf,S90.sl0.rf,AUC.sl.rf,S80.sl.rf,S90.sl.rf,
                AUC.sl0.combine,S80.sl0.combine,S90.sl0.combine,AUC.sl.combine,S80.sl.combine,S90.sl.combine)
    names(outputs) = c(sapply(1:length(base_functions),FUN=function(x){paste0(c("AUC.base","S80.base","S90.base","AUC.sl0","S80.sl0","S90.sl0","AUC.sl","S80.sl","S90.sl"),".",base_functions[x])}),
                       paste0(c("AUC.sl0","S80.sl0","S90.sl0","AUC.sl","S80.sl","S90.sl"),".combine"))
  #}
  #if (base_function == "combine")
  #{
  #  outputs = c(AUC.sl0,S80.sl0,S90.sl0,AUC.sl,S80.sl,S90.sl)
  #  names(outputs) = c("AUC.combine0","S80.combine0","S90.combine0","AUC.combine","S80.combine","S90.combine")
  #}
  return(outputs)
}






############################################################################
################# Function for the Full Binary Simulation: #################
############################################################################
fullsimu_binary = function(si,b.posthoc.glm,b.posthoc.qda,b.posthoc.rf)
{
  load("msi_fillrecenter_missing.RData")
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
    sse4[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma/1.5)
  }
  for (i in 1:9)
  {
    sse9[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma/1.5)
  }
  fillnew = mclapply(1:length(select),FUN=function(x){simu_binary_data(fillin,select[x],x,sse4,sse9)},mc.cores=n.cores)
  # fillnew = lapply(1:length(select),FUN=function(x){simu_binary_data(fillin,select[x],x,sse4,sse9)})
  print(paste0("Data Set ",si," Simulated"))
  DATA=rbindlist(lapply(fillnew,as.data.frame))
  DATA=as.data.frame(DATA)
  Region4=1*((DATA$x<=0)&(DATA$y>0)) + 2*((DATA$x>0)&(DATA$y>0)) + 3*((DATA$x<=0)&(DATA$y<=0))+ 4*((DATA$x>0)&(DATA$y<=0))
  Region9=1*((DATA$x<=2/3-1)&(DATA$y>1-2/3)) + 2*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y>1-2/3)) + 3*((DATA$x>2*2/3-1)&(DATA$y>1-2/3))+ 4*((DATA$x<=2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 5*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 6*((DATA$x>2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3))+ 7*((DATA$x<=2/3-1)&(DATA$y<=1-2*2/3)) + 8*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2*2/3)) + 9*((DATA$x>2*2/3-1)&(DATA$y<=1-2*2/3))
  DATA$Region4=Region4
  DATA$Region9=Region9
  DATA$cancer=as.factor(DATA$cancer)
  DATA=as.data.frame(DATA)
  ni=sapply(1:34,FUN=function(x){nrow(DATA[DATA[,"subject"]==x,])})
  #############################################################
  ##################### Model Evaluation: #####################
  #############################################################
  simuresults.glm = simu_binary(DATA,base_functions[1],b.posthoc.glm)
  print(paste0(base_functions[1]," ",si," Finished"))
  simuresults.qda = simu_binary(DATA,base_functions[2],b.posthoc.qda)
  print(paste0(base_functions[2]," ",si," Finished"))
  simuresults.rf = simu_binary(DATA,base_functions[3],b.posthoc.rf)
  print(paste0(base_functions[3]," ",si," Finished"))
  b.posthoc.combine = c(b.posthoc.glm,b.posthoc.qda,b.posthoc.rf)
  simuresults.combine = simu_binary(DATA,"combine",b.posthoc.combine)
  print(paste0("Combined ",si," Finished"))
  Outputs = c(simuresults.glm,simuresults.qda,simuresults.rf,simuresults.combine)
  return(Outputs)
}

############################################################################
############### Function for the Combined Binary Simulation: ###############
############################################################################
combinesimu_binary = function(si)
{
  load("msi_fillrecenter_missing.RData")
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
    sse4[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma/1.5)
  }
  for (i in 1:9)
  {
    sse9[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma/1.5)
  }
  fillnew = mclapply(1:length(select),FUN=function(x){simu_binary_data(fillin,select[x],x,sse4,sse9)},mc.cores=n.cores)
  print(paste0("Data Set ",si," Simulated"))
  DATA=rbindlist(lapply(fillnew,as.data.frame))
  DATA=as.data.frame(DATA)
  Region4=1*((DATA$x<=0)&(DATA$y>0)) + 2*((DATA$x>0)&(DATA$y>0)) + 3*((DATA$x<=0)&(DATA$y<=0))+ 4*((DATA$x>0)&(DATA$y<=0))
  Region9=1*((DATA$x<=2/3-1)&(DATA$y>1-2/3)) + 2*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y>1-2/3)) + 3*((DATA$x>2*2/3-1)&(DATA$y>1-2/3))+ 4*((DATA$x<=2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 5*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 6*((DATA$x>2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3))+ 7*((DATA$x<=2/3-1)&(DATA$y<=1-2*2/3)) + 8*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2*2/3)) + 9*((DATA$x>2*2/3-1)&(DATA$y<=1-2*2/3))
  DATA$Region4=Region4
  DATA$Region9=Region9
  DATA$cancer=as.factor(DATA$cancer)
  DATA=as.data.frame(DATA)
  ni=sapply(1:34,FUN=function(x){nrow(DATA[DATA[,"subject"]==x,])})
  #############################################################
  ##################### Model Evaluation: #####################
  #############################################################
  simuresults.combine = simu_binary(DATA,"combine")
  print(paste0("Combined ",si," Finished"))
  Outputs = simuresults.combine
  return(Outputs)
}


############################################################################
############### Function for the Combined Binary Simulation: ###############
############################################################################
newsimu_binary = function(si,b.posthoc)
{
  load("msi_fillrecenter_missing.RData")
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
    sse4[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma/1.5)
  }
  for (i in 1:9)
  {
    sse9[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma/1.5)
  }
  fillnew = mclapply(1:length(select),FUN=function(x){simu_binary_data(fillin,select[x],x,sse4,sse9)},mc.cores=n.cores)
  print(paste0("Data Set ",si," Simulated"))
  DATA=rbindlist(lapply(fillnew,as.data.frame))
  DATA=as.data.frame(DATA)
  Region4=1*((DATA$x<=0)&(DATA$y>0)) + 2*((DATA$x>0)&(DATA$y>0)) + 3*((DATA$x<=0)&(DATA$y<=0))+ 4*((DATA$x>0)&(DATA$y<=0))
  Region9=1*((DATA$x<=2/3-1)&(DATA$y>1-2/3)) + 2*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y>1-2/3)) + 3*((DATA$x>2*2/3-1)&(DATA$y>1-2/3))+ 4*((DATA$x<=2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 5*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 6*((DATA$x>2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3))+ 7*((DATA$x<=2/3-1)&(DATA$y<=1-2*2/3)) + 8*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2*2/3)) + 9*((DATA$x>2*2/3-1)&(DATA$y<=1-2*2/3))
  DATA$Region4=Region4
  DATA$Region9=Region9
  DATA$cancer=as.factor(DATA$cancer)
  DATA=as.data.frame(DATA)
  ni=sapply(1:34,FUN=function(x){nrow(DATA[DATA[,"subject"]==x,])})
  #############################################################
  ##################### Model Evaluation: #####################
  #############################################################
  simuresults.combine = simu_binary(DATA,"combine",b.posthoc)
  print(paste0("Combined ",si," Finished"))
  Outputs = simuresults.combine
  return(Outputs)
}


#########################################################################################
#########################################################################################


############################################################################
############## Function for Categorical Simulation per Methd: ##############
############################################################################
simu_category = function(DATA,base_function,Levels,b.po)
{
  if (base_function != "combine")
  {
    if (covsl == "category") {slcovnames = covariates = c("Global_pred","Region4_pred","Region9_pred")}
    if (covsl == "prob") {slcovnames = covariates = c("Global_pred1","Global_pred2","Region4_pred1","Region4_pred2","Region9_pred1","Region9_pred2")}
    pred.mat=matrix(NA,nrow(DATA),5)
    colnames(pred.mat) = paste0("pred.",c("base","sl10","sl1","sl20","sl2"))
  }
  if (base_function == "combine")
  {
    if (covsl == "category") {slcovnames = covariates =   c(sapply(1:length(base_functions),FUN=function(x){c(paste0("Global_pred_",base_functions[x]),paste0("Region4_pred_",base_functions[x]),paste0("Region9_pred_",base_functions[x]))}))}
    if (covsl == "prob") {slcovnames = covariates =   c(sapply(1:length(base_functions),FUN=function(x){c(paste0("Global_pred1_",base_functions[x]),paste0("Global_pred2_",base_functions[x]),paste0("Region4_pred1_",base_functions[x]),paste0("Region4_pred2_",base_functions[x]),paste0("Region9_pred1_",base_functions[x]),paste0("Region9_pred2_",base_functions[x]))}))}
    pred.mat=matrix(NA,nrow(DATA),4)
    colnames(pred.mat) = paste0("pred.",c("sl10","sl1","sl20","sl2"))
  }
  #######################################################
  for (split in c(1:length(folds)))
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
        # if (length(unique(train$category)) == 3)
        # {
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
          #}
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
      #print(fold)
    }
    
    ##################################################################################
    ############ Evaluate the trained super learner on the holdout dataset: ############
    ##################################################################################
    ######### Obtain the new X on the holdout data set:
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
    ############################ SL without smoothing ##########################
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
      Xcoord=tem[,"x"];Ycoord=tem[,"y"];
      for (jj in 1:length(slcovnames))
      {
        ######## smooth pij:
        X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[jj]])
        #b=0.2
        #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
        Spre=as.numeric(Smooth(X, sigma=b.po[jj],at="points",edge=TRUE, diggle=F))
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
      Xcoord=round(tem[,"x"],2);Ycoord=round(tem[,"y"],2);
      for (jj in 1:length(slcovnames))
      {
        ######## smooth pij:
        X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[jj]])
        #b=0.2
        #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
        Spre=as.numeric(Smooth(X, sigma=b.po[jj],at="points",edge=TRUE, diggle=F))
        tem[,covariates[jj]] = Spre
        holdoutsmooth[row.ind,covariates[jj]] = Spre
        ######## somooth qij:
      }
    }
    holdout = holdoutsmooth
    rm(holdoutsmooth)
    categorytrue = holdout[,"category"]  ##### notice that in each fold: should update "category" for the test slices
    
    ############################################################################
    ############################ SL with smoothing #############################
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
  
  if (base_function != "combine")
  {
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
    
    outputs = list(otable.base=otable.base, otable.sl10=otable.sl10, otable.sl20=otable.sl20, otable.sl1=otable.sl1, otable.sl2=otable.sl2, oerr=oerr, ofpr=ofpr, ofdr=ofdr)
    #names(outputs) = c("AUC.combine0","S80.combine0","S90.combine0","AUC.combine","S80.combine","S90.combine")
  }
  if (base_function == "combine")
  {
    ###### overall tables:
    otable.sl10 = evaluation_table(DATA$category,pred.mat[,"pred.sl10"],Levels)
    otable.sl20 = evaluation_table(DATA$category,pred.mat[,"pred.sl20"],Levels)
    otable.sl1 = evaluation_table(DATA$category,pred.mat[,"pred.sl1"],Levels)
    otable.sl2 = evaluation_table(DATA$category,pred.mat[,"pred.sl2"],Levels)
    ###### overall error rates:
    oerr.sl10 = 1 - sum(DATA$category == pred.mat[,"pred.sl10"])/nrow(DATA)
    oerr.sl20 = 1 - sum(DATA$category == pred.mat[,"pred.sl20"])/nrow(DATA)
    oerr.sl1 = 1 - sum(DATA$category == pred.mat[,"pred.sl1"])/nrow(DATA)
    oerr.sl2 = 1 - sum(DATA$category == pred.mat[,"pred.sl2"])/nrow(DATA)
    oerr = round(c(oerr.sl10,oerr.sl1,oerr.sl20,oerr.sl2),3)
    names(oerr) = c("sl10","sl1","sl20","sl2")
    ###### overall fpr:
    ofpr=round(matrix(c(sapply(1:Z,FUN=function(z){1-sum(otable.sl10[z,z])/sum(DATA$category==(z-1))}),
                        sapply(1:Z,FUN=function(z){1-sum(otable.sl1[z,z])/sum(DATA$category==(z-1))}),
                        sapply(1:Z,FUN=function(z){1-sum(otable.sl20[z,z])/sum(DATA$category==(z-1))}),
                        sapply(1:Z,FUN=function(z){1-sum(otable.sl2[z,z])/sum(DATA$category==(z-1))})),
                      nrow=4,ncol=Z,byrow=TRUE),3)
    rownames(ofpr) = c("sl10","sl1","sl20","sl2")
    ofdr=round(matrix(c(sapply(1:Z,FUN=function(z){1-sum(otable.sl10[z,z])/sum(otable.sl10[,z])}),
                        sapply(1:Z,FUN=function(z){1-sum(otable.sl1[z,z])/sum(otable.sl1[,z])}),
                        sapply(1:Z,FUN=function(z){1-sum(otable.sl20[z,z])/sum(otable.sl20[,z])}),
                        sapply(1:Z,FUN=function(z){1-sum(otable.sl2[z,z])/sum(otable.sl2[,z])})),
                      nrow=4,ncol=Z,byrow=TRUE),3)
    rownames(ofdr) = c("sl10","sl1","sl20","sl2")
    
    outputs = list(otable.sl10=otable.sl10, otable.sl20=otable.sl20, otable.sl1=otable.sl1, otable.sl2=otable.sl2, oerr=oerr, ofpr=ofpr, ofdr=ofdr)
  }
  return(outputs)
}


############################################################################
################# Function for the Full Binary Simulation: #################
############################################################################
fullsimu_category = function(si,bb)
{
  load("msi_fillrecenter_missing.RData")
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
  select=sample(1:46,size=length(index),replace=TRUE)
  ##### simulate data:
  sse4=sse9=list()
  for (i in 1:4)
  {
    sse4[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma/1.5)
  }
  for (i in 1:9)
  {
    sse9[[i]] = mvrnorm(n=1,mu=rep(0,4),Sigma=Gamma/1.5)
  }
  
  fillnewresults = mclapply(1:length(select),FUN=function(x){simu_category_data(fillin,select[x],x)},mc.cores=n.cores)
  fillnew = lapply(1:length(select),FUN=function(x){fillnewresults[[x]]$temp})
  qs.total = sapply(1:length(select),FUN=function(x){fillnewresults[[x]]$qs})
  qs.total = unlist(qs.total)
  qs.t1 = quantile(qs.total,t1)
  qs.t2 = quantile(qs.total,t2)
  
  Temp = mclapply(1:length(select),FUN=function(x){simu_category_data.addition(fillnew[[x]],sse4,sse9,qs.t1,qs.t2)},mc.cores=n.cores)
  fillnew = Temp
  rm(Temp)
  DATA = fillnew[[1]]
  for (i in index[-1])
  {
    DATA = rbind(DATA,fillnew[[i]])
  }
  #DATA=rbindlist(lapply(fillnew,as.data.frame))
  #DATA=as.data.frame(DATA)
  #Region4=1*((DATA$x<=0)&(DATA$y>0)) + 2*((DATA$x>0)&(DATA$y>0)) + 3*((DATA$x<=0)&(DATA$y<=0))+ 4*((DATA$x>0)&(DATA$y<=0))
  #Region9=1*((DATA$x<=2/3-1)&(DATA$y>1-2/3)) + 2*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y>1-2/3)) + 3*((DATA$x>2*2/3-1)&(DATA$y>1-2/3))+ 4*((DATA$x<=2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 5*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3)) + 6*((DATA$x>2*2/3-1)&(DATA$y<=1-2/3)&(DATA$y>1-2*2/3))+ 7*((DATA$x<=2/3-1)&(DATA$y<=1-2*2/3)) + 8*((DATA$x>2/3-1)&(DATA$x<=2*2/3-1)&(DATA$y<=1-2*2/3)) + 9*((DATA$x>2*2/3-1)&(DATA$y<=1-2*2/3))
  #DATA$Region4=Region4
  #DATA$Region9=Region9
  DATA$category=factor(DATA$category, ordered=TRUE, levels = c("0", "1", "2"))
  #DATA$category=as.factor(DATA$category)
  #DATA=as.data.frame(DATA)
  ni=sapply(1:length(select),FUN=function(x){nrow(DATA[DATA[,"subject"]==x,])})
  Levels = sort(unique(DATA$category))
  
  print(paste0("Data Set ",si," Simulated"))
  #############################################################
  ##################### Model Evaluation: #####################
  #############################################################
  b = bb[1:6]
  simuresults.glm = simu_category(DATA,base_functions[1],Levels,b)
  print(paste0(base_functions[1]," ",si," Finished"))
  b = bb[7:12]
  simuresults.qda = simu_category(DATA,base_functions[2],Levels,b)
  print(paste0(base_functions[2]," ",si," Finished"))
  b = bb[13:18]
  simuresults.rf = simu_category(DATA,base_functions[3],Levels,b)
  print(paste0(base_functions[3]," ",si," Finished"))
  simuresults.combine = simu_category(DATA,"combine",Levels,bb)
  print(paste0("Combined ",si," Finished"))
  Outputs = c(simuresults.glm,simuresults.qda,simuresults.rf,simuresults.combine)
  return(Outputs)
}


######### 08.14.2019
###### modified tuning functions:
b.posthoc.tuning.multiresolution = function(base_function) 
{
  #### output: 3 b values, for the 3 different resolutions
  if (base_function != "combine")
  {
    slcovnames = covariates = c("Global_pred","Region4_pred","Region9_pred")
  }
  if (base_function == "combine")
  {
    slcovnames = covariates = c(sapply(1:length(base_functions),FUN=function(x){c(paste0("Global_pred_",base_functions[x]),paste0("Region4_pred_",base_functions[x]),paste0("Region9_pred_",base_functions[x]))}))
  }
  b.posthoc = numeric()
  auc.sl=s90.sl=s80.sl=numeric()
  auc.sl0=s90.sl0=s80.sl0=numeric()
  auc.base=s90.base=s80.base=numeric()
  
  whole=cbind(DATA,Global_pred=0,Region4_pred=0,Region9_pred=0)
  
  for (split in c(1:length(folds)))
  {
    TRAIN=DATA[!(DATA$subject%in%folds[[split]]),]
    tr.subj = unique(TRAIN$subject)
    n.tr.subj = length(unlist(folds)) - length(folds[[split]])
    holdout = DATA[(DATA$subject%in%folds[[split]]),]
    holdout.ind = which(DATA[,"subject"]%in%folds[[split]])
    n.holdout = nrow(holdout)
    holdout.rownames=rownames(holdout)
    
    n.whole=nrow(TRAIN)
    p_r=p=3
    
    regional_prediction_holdout=matrix(0,nrow(holdout),n.X)
    
    MRIpars_TRAIN <- TRAIN[,MRIpars]
    MRIpars_holdout <- holdout[,MRIpars]
    
    ##### Global model:
    ###### train the candidate model on the training blocks:
    if (base_function == "rf")
    {
      global_fit <- randomForest(rf.baseformula, data = TRAIN, ntree=n.tree, mtry=2,nodesize=1, importance=T)
      global_pred <- predict(global_fit, newdata=holdout,type="prob")
      regional_prediction_holdout[,1]<- global_pred[,2]
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_TRAIN=TRAIN[TRAIN[,"Region4"]==reg,]
        local_fit <- randomForest(rf.baseformula, data = regional_TRAIN, ntree=n.tree, mtry=2,
                                  nodesize=1, importance=T)
        local_pred=predict(local_fit, newdata=holdout,type="prob")
        regional_prediction_holdout[,2]<- local_pred[,2] * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,2]
      }
      
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_TRAIN=TRAIN[TRAIN[,"Region9"]==reg,]
        local_fit <- randomForest(rf.baseformula, data = regional_TRAIN, ntree=n.tree, mtry=2,
                                  nodesize=1, importance=T)
        local_pred=predict(local_fit, newdata=holdout,type="prob")
        regional_prediction_holdout[,3]<- local_pred[,2] * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,3]
      }
      colnames(regional_prediction_holdout)=slcovnames
      whole[holdout.rownames,slcovnames]=regional_prediction_holdout
    }
    
    if (base_function == "knn")
    {
      global_fit=knn(train=TRAIN[,MRIpars],test=holdout[,MRIpars],cl=TRAIN$cancer,k=knn.k,prob=TRUE)
      regional_prediction_holdout[,1]<- 1-attr(global_fit,"prob")
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_TRAIN=TRAIN[TRAIN[,"Region4"]==reg,]
        local_fit <- knn(train=regional_TRAIN[,MRIpars],test=holdout[,MRIpars],cl=regional_TRAIN$cancer,k=knn.k,prob=TRUE)
        local_pred = 1-attr(local_fit,"prob")
        regional_prediction_holdout[,2] <- local_pred * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,2]
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_TRAIN=TRAIN[TRAIN[,"Region9"]==reg,]
        local_fit <- knn(train=regional_TRAIN[,MRIpars],test=holdout[,MRIpars],cl=regional_TRAIN$cancer,k=knn.k,prob=TRUE)
        local_pred = 1-attr(local_fit,"prob")
        regional_prediction_holdout[,3] <- local_pred * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,3]
      }
      colnames(regional_prediction_holdout)=slcovnames
      whole[holdout.rownames,slcovnames]=regional_prediction_holdout
    }
    
    if (base_function == "qda")
    {
      global_fit<-qda(TRAIN[,MRIpars], TRAIN$cancer)
      global_pred<-predict(global_fit, holdout[,MRIpars])
      regional_prediction_holdout[,1]<- global_pred$posterior[,2]
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_TRAIN=TRAIN[TRAIN[,"Region4"]==reg,]
        local_fit <- qda(regional_TRAIN[,MRIpars], regional_TRAIN$cancer)
        local_pred=predict(local_fit,holdout[,MRIpars])
        regional_prediction_holdout[,2]<- local_pred$posterior[,2] * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,2]
      }
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_TRAIN=TRAIN[TRAIN[,"Region9"]==reg,]
        local_fit <- qda(regional_TRAIN[,MRIpars], regional_TRAIN$cancer)
        local_pred=predict(local_fit,holdout[,MRIpars])
        regional_prediction_holdout[,3]<- local_pred$posterior[,2] * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,3]
      }
      colnames(regional_prediction_holdout)=slcovnames
      whole[holdout.rownames,slcovnames]=regional_prediction_holdout
    }
    if (base_function == "glm")
    {
      global_fit <- glm(rf.baseformula, data = TRAIN, family = binomial(link = "logit")) 
      global_pred <- predict(global_fit, newdata=holdout,type="response")
      regional_prediction_holdout[,1]<- global_pred
      ##### Local model with 4 regions (2 by 2 split):
      for (reg in 1:resolution4)
      {
        regional_TRAIN=TRAIN[TRAIN[,"Region4"]==reg,]
        local_fit <- glm(rf.baseformula, data = regional_TRAIN, family = binomial(link = "logit")) 
        local_pred=predict(local_fit, newdata=holdout,type="response")
        regional_prediction_holdout[,2]<- local_pred * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,2]
      }
      
      ##### Local model with 9 regions (3 by 3 split):
      for (reg in 1:resolution9)
      {
        regional_TRAIN=TRAIN[TRAIN[,"Region9"]==reg,]
        local_fit <- glm(rf.baseformula, data = regional_TRAIN, family = binomial(link = "logit")) 
        local_pred=predict(local_fit, newdata=holdout,type="response")
        regional_prediction_holdout[,3]<- local_pred * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,3]
      }
      colnames(regional_prediction_holdout)=slcovnames
      whole[holdout.rownames,slcovnames]=regional_prediction_holdout
    }
  }
  cancertrue = as.numeric(as.character(whole[,"cancer"]))  ##### notice that in each fold: should update "category" for the test slices
  ################################################################
  ########################## Smoothing: ##########################
  ################################################################
  #### tuning:
  B = seq(0.01,0.3,by=0.01)
  for (reso in 1:3)
  {
    auc=numeric()
    for (t in 1:length(B))
    {
      bb = B[t]
      wholesmooth = whole
      for (ii in index)
      {
        tem = whole[whole$subject == ii,]
        row.ind = which(whole$subject == ii)
        Xcoord=tem[,"x"];Ycoord=tem[,"y"];
        X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[reso]])
        Spre=as.numeric(Smooth(X, sigma=bb,at="points",edge=TRUE, diggle=F))
        tem[,covariates[reso]] = Spre
        wholesmooth[row.ind,covariates[reso]] = Spre
      }
      bpred <- prediction(predictions=wholesmooth[,covariates[reso]], labels=cancertrue)
      bperf<- performance(bpred,"tpr","fpr")
      bau <- performance(bpred,"auc")
      auc[t] <- unlist(slot(bau, "y.values"))
    }
    b.posthoc[reso] = B[which.max(auc)]
  }
  return(b.posthoc)
}


###################### new tuning function for b (categorical)
b.posthoc.tuning.category.multiresolution = function(base_function)
{
  whole=cbind(DATA,Global_pred1=0,Global_pred2=0,Region4_pred1=0,Region4_pred2=0,Region9_pred1=0,Region9_pred2=0)
  pred.mat=matrix(NA,nrow(DATA),1)
  colnames(pred.mat) = "pred.sl"
  for (split in c(1:length(folds)))
  {
    TRAIN=DATA[!(DATA$subject%in%folds[[split]]),]
    tr.subj = unique(TRAIN$subject)
    n.tr.subj = length(unlist(folds)) - length(folds[[split]])
    holdout = DATA[(DATA$subject%in%folds[[split]]),]
    holdout.ind = which(DATA[,"subject"]%in%folds[[split]])
    n.holdout = nrow(holdout)
    holdout.rownames=rownames(holdout)
    p_r=p=3
    n.whole=nrow(TRAIN)
    ############## Obtain the new X for super learner:
    regional_prediction_test=matrix(0,nrow(test),n.X)
    MRIpars_TRAIN <- TRAIN[,MRIpars]
    MRIpars_holdout <- holdout[,MRIpars]

      ##### Global model:
      ###### train the candidate model on the training blocks:
      if (base_function == "rf")
      {
        global_fit <- randomForest(rf.baseformula, data = TRAIN, ntree=n.tree, mtry=2,
                                   nodesize=1, importance=T)
        if (covsl == "prob")
        {
          global_pred <- predict(global_fit, newdata=holdout,type="prob")
          pred.holdout=global_pred
          regional_prediction_holdout[,1:2]<- pred.holdout[,c(1,2)]
        } 
        if (covsl == "category")
        {
          global_pred <- predict(global_fit, newdata=holdout,type="response")
          pred.holdout=as.numeric(as.character(global_pred)) #factor -> numerical
          regional_prediction_test[,1]<- pred.holdout
        }
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_TRAIN=TRAIN[TRAIN[,"Region4"]==reg,]
          local_fit <- randomForest(rf.baseformula, data = regional_TRAIN, ntree=n.tree, mtry=2,
                                    nodesize=1, importance=T)
          if (covsl == "prob")
          {
            local_pred=predict(local_fit, newdata=holdout,type="prob")
            pred.holdout=local_pred
            regional_prediction_holdout[,3:4]<- pred.holdout[,c(1,2)] * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,3:4]
          }
          if (covsl == "category")
          {
            local_pred=predict(local_fit, newdata=holdout,type="response")
            pred.holdout=as.numeric(as.character(local_pred)) #factor -> numerical
            regional_prediction_holdout[,2]<- pred.holdout * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,2]
          }
        }
        
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_TRAIN=TRAIN[TRAIN[,"Region9"]==reg,]
          local_fit <- randomForest(rf.baseformula, data = regional_TRAIN, ntree=n.tree, mtry=2,
                                    nodesize=1, importance=T)
          if (covsl == "prob")
          {
            local_pred=predict(local_fit, newdata=holdout,type="prob")
            pred.holdout=local_pred
            regional_prediction_holdout[,5:6]<- pred.holdout[,c(1,2)] * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,5:6]
          }
          if (covsl == "category")
          {
            local_pred=predict(local_fit, newdata=holdout,type="response")
            pred.holdout=as.numeric(as.character(local_pred)) #factor -> numerical
            regional_prediction_holdout[,3]<- pred.holdout * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,3]
          }
        }
        colnames(regional_prediction_holdout)=slcovnames
        whole[holdout.rownames,slcovnames]=regional_prediction_holdout
      }
      
      
      if (base_function == "knn")
      {
        global_fit=knn(train=TRAIN[,MRIpars],test=holdout[,MRIpars],cl=TRAIN$category,k=knn.k,prob=FALSE)
        regional_prediction_holdout[,1]<- as.numeric(as.character(global_fit))
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_TRAIN=TRAIN[TRAIN[,"Region4"]==reg,]
          local_fit <- knn(train=regional_TRAIN[,MRIpars],test=holdout[,MRIpars],cl=regional_TRAIN$category,k=knn.k,prob=FALSE)
          local_pred = as.numeric(as.character(local_fit))
          regional_prediction_holdout[,2] <- local_pred * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,2]
        }
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_TRAIN=TRAIN[TRAIN[,"Region9"]==reg,]
          local_fit <- knn(train=regional_TRAIN[,MRIpars],test=holdout[,MRIpars],cl=regional_TRAIN$category,k=knn.k,prob=FALSE)
          local_pred = as.numeric(as.character(local_fit))
          regional_prediction_holdout[,3] <- local_pred * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,3]
        }
        colnames(regional_prediction_holdout)=slcovnames
        whole[holdout.rownames,slcovnames]=regional_prediction_holdout
      }
      
      if (base_function == "qda")
      {
        global_fit<-qda(TRAIN[,MRIpars], TRAIN$category)
        global_pred<-predict(global_fit, holdout[,MRIpars])
        if (covsl == "prob")
        {
          pred.holdout=global_pred$posterior
          regional_prediction_holdout[,1:2]<- pred.holdout[,c(1,2)]
        }
        if (covsl == "category")
        {
          pred.holdout=as.numeric(as.character(global_pred$class)) #factor -> numerical
          regional_prediction_holdout[,1]<- pred.holdout
        }
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_TRAIN=TRAIN[TRAIN[,"Region4"]==reg,]
          local_fit <- qda(regional_TRAIN[,MRIpars], regional_TRAIN$category)
          local_pred=predict(local_fit,holdout[,MRIpars])
          if (covsl == "prob")
          {
            pred.holdout=local_pred$posterior
            regional_prediction_holdout[,3:4]<- pred.holdout[,c(1,2)] * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,3:4]
          }
          if (covsl == "category")
          {
            pred.holdout=as.numeric(as.character(local_pred$class))
            regional_prediction_holdout[,2]<- pred.holdout * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,2]
          }
        }
        
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_TRAIN=TRAIN[TRAIN[,"Region9"]==reg,]
          local_fit <- qda(regional_TRAIN[,MRIpars], regional_TRAIN$category) 
          local_pred=predict(local_fit,holdout[,MRIpars])
          if (covsl == "prob")
          {
            pred.holdout=local_pred$posterior
            regional_prediction_holdout[,5:6]<- pred.holdout[,c(1,2)] * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,5:6]
          }
          if (covsl == "category")
          {
            pred.holdout=as.numeric(as.character(local_pred$class))
            regional_prediction_holdout[,3]<- pred.holdout * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,3]
          }
        }
        colnames(regional_prediction_holdout)=slcovnames
        whole[holdout.rownames,slcovnames]=regional_prediction_holdout
      }
      if (base_function == "glm")
      {
        global_fit <- polr(rf.baseformula, data = TRAIN, method = "probit",Hess=TRUE) 
        if (covsl == "prob")
        {
          global_pred <- predict(global_fit, newdata=holdout,type="probs")
          pred.holdout=global_pred
          regional_prediction_holdout[,1:2]<- pred.holdout[,c(1,2)]
        } 
        if (covsl == "category")
        {
          global_pred <- predict(global_fit, newdata=holdout,type="class")
          pred.holdout=as.numeric(as.character(global_pred)) 
          regional_prediction_holdout[,1]<- pred.holdout
        }
        ##### Local model with 4 regions (2 by 2 split):
        for (reg in 1:resolution4)
        {
          regional_TRAIN=TRAIN[TRAIN[,"Region4"]==reg,]
          local_fit <- polr(rf.baseformula, data = regional_TRAIN, method = "probit",Hess=TRUE) 
          if (covsl == "prob")
          {
            local_pred=predict(local_fit, newdata=holdout,type="prob")
            pred.holdout=local_pred
            regional_prediction_holdout[,3:4]<- pred.holdout[,c(1,2)] * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,3:4]
          }
          if (covsl == "category")
          {
            local_pred=predict(local_fit, newdata=holdout,type="class")
            pred.holdout=as.numeric(as.character(local_pred)) 
            regional_prediction_holdout[,2]<- pred.holdout * (holdout[,"Region4"]==reg) + regional_prediction_holdout[,2]
          }
        }
        
        ##### Local model with 9 regions (3 by 3 split):
        for (reg in 1:resolution9)
        {
          regional_TRAIN=TRAIN[TRAIN[,"Region9"]==reg,]
          local_fit <- polr(rf.baseformula, data = regional_TRAIN, method = "probit",Hess=TRUE) 
          if (covsl == "prob")
          {
            local_pred=predict(local_fit, newdata=holdout,type="prob")
            pred.holdout=local_pred
            regional_prediction_holdout[,5:6]<- pred.holdout[,c(1,2)] * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,5:6]
          }
          if (covsl == "category")
          {
            local_pred=predict(local_fit, newdata=holdout,type="class")
            pred.holdout=as.numeric(as.character(local_pred)) 
            regional_prediction_holdout[,3]<- pred.holdout * (holdout[,"Region9"]==reg) + regional_prediction_holdout[,3]
          }
        }
        colnames(regional_prediction_holdout)=slcovnames
        whole[holdout.rownames,slcovnames]=regional_prediction_holdout
      }
      #print(fold)
    }

    ################################################################
    ########################## Smoothing: ##########################
    ################################################################
  B = seq(0.01,0.3,by=0.01)
  for (reso in 1:length(slcovnames))
  {
    auc=numeric()
    for (t in 1:length(B))
    {
      bb = B[t]
      wholesmooth = whole
      for (ii in index)
      {
        tem = whole[whole$subject == ii,]
        row.ind = which(whole$subject == ii)
        Xcoord=tem[,"x"];Ycoord=tem[,"y"];
        X <- ppp(Xcoord, Ycoord, c(min(Xcoord),max(Xcoord)),c(min(Ycoord),max(Ycoord)),marks=tem[,covariates[reso]])
        Spre=as.numeric(Smooth(X, sigma=bb,at="points",edge=TRUE, diggle=F))
        tem[,covariates[reso]] = Spre
        wholesmooth[row.ind,covariates[reso]] = Spre
      }
      bpred <- prediction(predictions=wholesmooth[,covariates[reso]], labels=cancertrue)
      bperf<- performance(bpred,"tpr","fpr")
      bau <- performance(bpred,"auc")
      auc[t] <- unlist(slot(bau, "y.values"))
    }
    b.posthoc[reso] = B[which.max(auc)]
  }
  return(b.posthoc)
}



