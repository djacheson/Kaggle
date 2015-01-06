setwd("/home/dan/Kaggle/UPenn_Seizure/clips")

library(caret)
library(data.table)
library(stringr)
library(kernlab)

library(doMC)
registerDoMC()

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]


#Calculate optimal sigma parameter for each subject using sigest from kernlab
sigma <- list()
for(i in seq_along(dirs)) {

    load(file=paste0(dirs[[i]],"/train+test.freqRatios.wide.rbin"))
    train.x <- as.matrix(train.wide[,c("trialType","trialNum","latency_15") := NULL])
    sigma[[i]] <- sigest(train.x, frac=0.75)[[2]]
}


#list object to hold all predictions
allModels <- list()
allPreds <- list()

#Loop through all subjects
for(i in seq_along(dirs)) {
  start = proc.time()
  load(file=paste0(dirs[[i]],"/train+test.freqRatios.wide.rbin"))

  
  currSub <- str_replace(dirs[[i]],"[/ .][/ /]","")
  
  clips = paste0(currSub,"_test_segment_",test.wide$trialNum,".mat")
  
  #Create test file
  test.x <- test.wide[,c("trialType","trialNum") := NULL]
  
  
  #train control for all models
  trControl <- trainControl(## 10-fold CV
    method = "cv",
    number = 10,
    ## repeated two
    #repeats = 2,
    classProbs = TRUE)
  
  #svm tuning grid for all models
  #using sigma estimate from above
  svm_Grid <-  expand.grid(sigma = sigma[[i]],
                           C = 2^seq(-6,1,length=15))
  
#   
#   **************
#   ictal or not 
#   **************
  train.copy <- copy(train.wide)
  train.y <- as.factor(ifelse(train.wide$trialType=="ictal","X1","X0"))
  train.x <- as.matrix(train.copy[,c("trialType","trialNum","latency_15") := NULL])
  rm(train.copy)
  
  #upsample to deal with class imbalance
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y
  
  
  #SVM
  
  ictal.svm <- train(train.y ~ . , data=upSamp,
                  method = "svmRadial",
                  trControl = trControl,
                  verbose = TRUE,
                  tuneGrid = svm_Grid)
  
  
  ictal.svm <- allModels$ictal[[i]]
  
  #Make predictions
  ictalProb <- round(predict(ictal.svm,test.x,type = "prob"),2)
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
  
  
  #upsampling to deal with class imbalance
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y

  
  early.svm <- train(train.y ~ . , data=upSamp,
                    method = "svmRadial",
                    trControl = trControl,
                    verbose = TRUE,
                    tuneGrid = svm_Grid)
  
  
  #Make predictions
  earlyProb <- round(predict(early.svm,test.x,type = "prob"),2)

  preds[,early := earlyProb$X1]
  
  #Add models and predictions to list objects, clean workspace and save
  allModels$subs[[i]] <- currSub
  allModels$ictal[[i]] <- ictal.svm
  allModels$early[[i]] <- early.svm
  allPreds[[i]] <- preds
  
  rm(early.svm)
  rm(ictal.svm)
  rm(preds)
  gc()

  print(paste0("Done with subject ", dirs[[i]]," in: ", (proc.time()-start)[[3]]))
  save(list=c("allPreds","allModels"), file=paste0("../models/SVM.freqRatios.earlyFix.rbin"))
} # end big directory loop

#Combine all predictions into a single table
allPredComb <- rbindlist(allPreds)
write.csv(allPredComb,file=paste0("../models/SVM.freqRatios.preds.earlyFix.csv"), row.names=FALSE)
