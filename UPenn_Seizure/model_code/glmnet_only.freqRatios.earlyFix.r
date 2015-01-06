#GLMnet models

setwd("/home/dan/Kaggle/UPenn_Seizure/clips")

library(caret)
library(data.table)
library(stringr)
library(glmnet)

library(doMC)
registerDoMC()

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]

#lists to hold models and predictions per subject
allModels <- list()
allPreds <- list()

#Loop along each subject and generate models and predictions for ictal/non-ictal and latency <= 15
for(i in seq_along(dirs)) {
  start = proc.time()
  load(file=paste0(dirs[[i]],"/train+test.freqRatios.wide.rbin"))

  currSub <- str_replace(dirs[[i]],"[/ .][/ /]","")
  clips = paste0(currSub,"_test_segment_",test.wide$trialNum,".mat")
  
  #Create test file
  test.x <- as.matrix(test.wide[,c("trialType","trialNum") := NULL])
  
  #**************
  #ictal or not 
  #**************
  train.copy <- copy(train.wide)
  #A few changes to the data for use with glmnet
  train.y <- as.factor(ifelse(train.wide$trialType=="ictal","X1","X0"))
  train.x <- as.matrix(train.copy[,c("trialType","trialNum","latency_15") := NULL])
  rm(train.copy)
  
  #upsample seizures to deal with severe class imbalance
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y
  
  
  #elastic net
  glmnet.fit.ictal <- cv.glmnet(x = as.matrix(up.x), 
                                y=up.y, 
                                family="binomial", 
                                nfolds=10, 
                                nlambda=300, 
                                type.measure = "auc", 
                                parallel=TRUE)   


  glmnet.fit.ictal <- allModels$ictal[[i]]
  #Make predictions
  ictalProb <- round(predict(glmnet.fit.ictal,test.x,type="response"),2)
  ictalClass <- ifelse(ictalProb < 0.5, "interictal","ictal")
  
  #add to preds data.table
  preds <- data.table(clip=clips, seizure=as.numeric(ictalProb))
  
  
  #---------------------------------------------------------------------------------------------------
  #*************
  #Early Latency
  #************
  train.copy <- copy(train.wide)
  train.y <- ifelse(train.copy$latency_15 == 1,"X1","X0")
  train.y <- as.factor(ifelse(is.na(train.y),"X0",train.y))
  train.x <- as.matrix(train.copy[,c("trialType","trialNum","latency_15") := NULL])
  
  
  #upsampling
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y
  
  
  glmnet.fit.early <- cv.glmnet(x = as.matrix(up.x), 
                                y=up.y, 
                                family="binomial", 
                                nfolds=10, 
                                nlambda=300, 
                                type.measure = "auc",
                                parallel=TRUE)  
  
  #Get predictions
  earlyProb <- round(predict(glmnet.fit.early,test.x,type="response"),2)
  preds[,early := earlyProb]
  
  #add models and predictions to lists holding 
  allModels$ictal[[i]] <- glmnet.fit.ictal
  allModels$early[[i]] <- glmnet.fit.early
  allPreds[[i]] <- preds
  

  print(paste0("Done with subject ", dirs[[i]]," in: ", (proc.time()-start)[[3]]))
  save(list=c("allPreds","allModels"), file=paste0("../models/glmnetOnly.freqRatios.earlyFix.rbin"))
  
}

#Combine predictions across subjects into a single .csv file
allPredComb <- rbindlist(allPreds)
write.csv(allPredComb,file=paste0("../models/glmnetOnly.freqRatios.preds.earlyFix.csv"),row.names=FALSE)
