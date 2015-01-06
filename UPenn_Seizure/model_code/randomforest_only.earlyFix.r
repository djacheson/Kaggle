#Random forest on original feature space (i.e., no frequency ratios included)

setwd("/home/dan/Kaggle/UPenn_Seizure/clips")

library(caret)
library(data.table)

library(doMC)
registerDoMC()

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]

#list objects to hold all models and predictions
allModels <- list()
allPreds <- list()

#loop along directories to generate models and predictions
for(i in seq_along(dirs)) {
  start = Sys.time()
  load(file=paste0(dirs[[i]],"/train+test.wide.rbin"))
 
  currSub <- str_replace(dirs[[i]],"[/ .][/ /]","")
  
  clips = paste0(currSub,"_test_segment_",test.wide$trialNum,".mat")
  #Create test file
  test.x <- test.wide[,c("trialType","trialNum") := NULL]
  
  
  #train control for all models
  trControl <- trainControl(## 10-fold CV
    method = "repeatedcv",
    number = 10,
    ## repeated five
    repeats = 5,
    classProbs = TRUE)
  
  #random forest tuning grid for all models
  rf_Grid <-  expand.grid(mtry = seq(5,21,by=2))
  
  
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
  
  
  #Parallel Random forest
  ictal.rf <- train(train.y ~ . , data=upSamp,
                  method = "parRF",
                  trControl = trControl,
                  verbose = TRUE,
                  tuneGrid = rf_Grid
                  )
  
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
  
  
  #upsampling
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y

  
  early.rf <- train(train.y ~ . , data=upSamp,
                    method = "parRF",
                    trControl = trControl,
                    verbose = TRUE,
                    tuneGrid = rf_Grid
                    )
  
  #Make predictions
  earlyProb <- round(predict(early.rf,test.x,type = "prob"),2)
  preds[,early := earlyProb$X1]
  
  #add predictions and models to list
  allModels$ictal[[i]] <- ictal.rf
  allModels$early[[i]] <- early.rf
  allPreds[[i]] <- preds
  
  save(list=c("allPreds","allModels"), file=paste0("../models/randomForest.earlyFix.rbin"))
  print(paste0("Done with subject ", dirs[[i]]," in: ", Sys.time()-start))
}

#Combine predictions into a data.table
allPredComb <- rbindlist(allPreds)
write.csv(allPredComb,file=paste0("../models/randomForest.preds.earlyFix.csv"),row.names = FALSE)

