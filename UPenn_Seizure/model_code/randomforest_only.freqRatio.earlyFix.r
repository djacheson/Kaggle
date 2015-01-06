#Random forest models using caret

setwd("/home/dan/Kaggle/UPenn_Seizure/clips")

library(caret)
library(data.table)
library(stringr)

library(doMC)
registerDoMC()

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]

load(file=paste0("../models/randomForest.freqRatios.rbin"))

#list object to hold all predictions
allModels <- list()
allPreds <- list()

for(i in seq_along(dirs)) {
  start = proc.time()
  load(file=paste0(dirs[[i]],"/train+test.freqRatios.wide.rbin"))
  
  currSub <- str_replace(dirs[[i]],"[/ .][/ /]","")
  
  clips = paste0(currSub,"_test_segment_",test.wide$trialNum,".mat")
  
  #Create test file
  test.x <- test.wide[,c("trialType","trialNum") := NULL]
  
  
  #train control for all models
  trControl <- trainControl(## 10-fold CV
    method = "repeatedcv",
    number = 10,
    ## repeated two
    repeats = 2,
    classProbs = TRUE)
  
  #random forest tuning grid for all models
  rf_Grid <-  expand.grid(mtry = seq(1,21,by=2))
  
  
  #**************
  #ictal or not 
  #**************
  train.copy <- copy(train.wide)
  train.y <- as.factor(ifelse(train.wide$trialType=="ictal","X1","X0"))
  train.x <- as.matrix(train.copy[,c("trialType","trialNum","latency_15") := NULL])
  rm(train.copy)
  
  #upsample to deal with class imbalance
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y
  
  
  #Random forest
  ictal.rf <- train(train.y ~ . , data=upSamp,
                  method = "parRF",
                  trControl = trControl,
                  verbose = TRUE,
                  tuneGrid = rf_Grid)
  
  
  ictal.rf <- allModels$ictal[[i]]
  #Make predictions
  ictalProb <- round(predict(ictal.rf,test.x,type = "prob"),2)
  #ictalClass <- ifelse(ictalProb < 0.5, "interictal","ictal")
  
  #add to preds data.table
  preds <- data.table(clip=clips, seizure=as.numeric(ictalProb$X1))
  
  
  #---------------------------------------------------------------------------------------------------
  #*************
  #Early Latency
  #************
  train.copy <- copy(train.wide)
  train.y <- ifelse(train.copy$latency_15 == 1,"X1","X0")
  train.y <- as.factor(ifelse(is.na(train.y),"X0",train.y))
  train.x <- as.matrix(train.copy[,c("trialType","trialNum","latency_15") := NULL])
  
  
  #upsample to deal witih class imbalance
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y

  
  early.rf <- train(train.y ~ . , data=upSamp,
                    method = "parRF",
                    trControl = trControl,
                    verbose = TRUE,
                    tuneGrid = rf_Grid)
  
  
  #Make predictions
  earlyProb <- round(predict(early.rf,test.x,type = "prob"),2)

  preds[,early := earlyProb$X1]
  
  allModels$subs[[i]] <- currSub
  allModels$ictal[[i]] <- ictal.rf
  allModels$early[[i]] <- early.rf
  allPreds[[i]] <- preds
  
  rm(early.rf)
  rm(ictal.rf)
  rm(preds)
  gc()
  print(paste0("Done with subject ", dirs[[i]]," in: ", (proc.time()-start)[[3]]))
  save(allPreds, file=paste0("../models/randomForest.freqRatios.allPreds.earlyFix.rbin"))
}


#Create table of model predictions
allPredComb <- rbindlist(allPreds)
write.csv(allPredComb,file=paste0("../models/randomForest.freqRatios.preds.earlyFix.csv"),row.names=FALSE)
