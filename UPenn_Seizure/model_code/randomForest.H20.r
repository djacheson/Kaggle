#Random forest using H2O

setwd("/home/dan/Kaggle/UPenn_Seizure/clips")

library(caret)
library(data.table)
library(stringr)
library(h2o)

#initialize an h2o instance
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE,
                    Xmx = '2g')

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]


#list object to hold all predictions
allPerformance <- list()
allPreds <- list()
trainPreds <- list()

for(i in #seq_along(dirs)) {
  start = proc.time()
  load(file=paste0(dirs[[i]],"/train+test.freqRatios.wide.rbin"))
  
  currSub <- str_replace(dirs[[i]],"[/ .][/ /]","")
  
  clips = paste0(currSub,"_test_segment_",test.wide$trialNum,".mat")
  
  #Create test file
  test.x <- test.wide[,c("trialType","trialNum") := NULL]
  test.h2o <- as.h2o(localH2O, test.x, key = 'test.h20')
  
  
  #**************
  #ictal or not 
  #**************
  train.copy <- data.frame(copy(train.wide))
  part <- createDataPartition(train.copy$trialType,p = 0.8, times = 1,list=FALSE)
  train.copy$trialNum <- NULL
  train.copy$latency_15 <- NULL
  
  validation <- train.copy[-part,]
  train <- train.copy[part,]
  #train.y <- as.factor(ifelse(train.wide$trialType=="ictal","X1","X0"))
  
  #train.trialNum <- train.copy$trialNum
  
  #train.x <- as.matrix(train.copy[,c("trialType","trialNum","latency_15") := NULL])
  rm(train.copy)
  

  #load training and validation data into h2o
  ictal.validation.h2o <- as.h2o(localH2O, validation, key = 'ictal.validation.h20')
  ictal.train.h2o <- as.h2o(localH2O, train, key = 'ictal.train.h20')
  
  colTrain <- ncol(train)

  ictal.rf <- h2o.randomForest(x = 2:colTrain, # column numbers for predictors
                               y = 1, # column number for label
                               data = ictal.train.h2o,
                               validation = ictal.validation.h2o,
                               classification=TRUE,
                               balance.classes = TRUE,
                               ntree=400,
                               depth = 2
                              )
  
  #save(list=c("ictal.deepBelief"), file=paste0("../models/DeepBelief/",currSub,"_ictalDeepBelief.rbin"))
  
  trainPredict <- h2o.predict(ictal.rf)
  performance_ictal <- h2o.performance(trainPredict$interictal, ictal.train.h2o$trialType,measure = "F1")
  cutoff <- performance_ictal@model$best_cutoff
  interictal <- as.numeric(as.character(data.frame(as.matrix(h2o.predict(ictal.rf, test.h2o)))$interictal))
  
  #get performance using cutoff, scaled by max value
  ictal <- (1 - interictal) / cutoff
  ictal <- ictal / max(ictal)
  print(performance_ictal)

  preds <- data.table(clip=clips, seizure=ictal)

  
  #Some code to save model performance
  #allPerformance$subject[[i]] <- currSub
  #allPerformance$ictal[[i]] <- print(performance_ictal)
  
  
  #predict training
#   trainProb <- predict(ictal.deepBelief, test.h20)$ictal
#   trainProb <- Reduce(cbind,trainProbs)
#   trainProb <- trainProb - 1
#   trainProb <- rowSums(trainProb) / ncol(trainProb)
#   
#   trainPreds$ictal[[i]] <- data.table(subject = currSub, 
#                                       trialNum = train.trialNum, 
#                                       elm.freq = trainProb, 
#                                       target=train.y)
#   
  
  
  #---------------------------------------------------------------------------------------------------
  #*************
  #Early Latency
  #************
  train.copy <- data.frame(copy(train.wide))
  train.copy$trialNum <- NULL
  train.copy$trialType <- NULL
  train.copy$latency_15 <- ifelse(is.na(train.copy$latency_15) |train.copy$latency_15==0 ,0,1)
  train.copy$latency_15 <- as.factor(ifelse(train.copy$latency_15 == 1, "early","not_early"))
  part <- createDataPartition(train.copy$latency_15,p = 0.8, times = 1,list=FALSE)
  
  validation <- train.copy[-part,]
  train <- train.copy[part,]
  rm(train.copy)
  
  #load data into h2o
  early.validation.h2o <- as.h2o(localH2O, validation, key = 'early.train.h20')
  early.train.h2o <- as.h2o(localH2O, train, key = 'early.train.h20')
  
  colTrain <- ncol(train)
  early.rf <- h2o.randomForest(x = 1:(colTrain-1), # column numbers for predictors
                               y = colTrain, # column number for label
                               data = early.train.h2o,
                               validation = early.validation.h2o,
                               classification=TRUE,
                               balance.classes = TRUE,
                               ntree=400,
                               depth = 2
                              )
  

  trainPredict <- h2o.predict(early.rf)
  performance_early <- h2o.performance(trainPredict$not_early, early.train.h2o$latency_15,measure = "F1")
  cutoff <- performance_early@model$best_cutoff
  not_earlyProb <- as.numeric(as.character(data.frame(as.matrix(h2o.predict(early.rf, test.h2o)))$not_early))
  early <- (1 - not_earlyProb) / cutoff
  early <- early / max(early)
  print(performance_early)  

  preds[,early := early]

#allPerformance$early[[i]] <- performance_early
  
  #predict training
#   trainProb <- Reduce(cbind,trainProbs)
#   trainProb <- trainProb - 1
#   trainProb <- rowSums(trainProb) / ncol(trainProb)
#   
#   trainPreds$early[[i]] <- data.table(subject = currSub, 
#                                       trialNum = train.trialNum, 
#                                       elm.freq = trainProb, 
#                                       target=train.y)
#   
  
  allPreds[[i]] <- preds
  
  rm(preds)
# 
#   file2rm <- h2o.ls(localH2O)
#   h2o.rm(localH2O,file2rm$Key)  

  print(paste0("Done with subject ", dirs[[i]]," in: ", (proc.time()-start)[[3]]))
  save(list=c("allPreds"), file=paste0("../models/h2o_randomForest.depth2.freqRatios.allPreds.rbin"))

}

#Combine predictions
allPredComb <- rbindlist(allPreds)
write.csv(allPredComb,file=paste0("../models/h2o_randomForest.depth2.freqRatios.preds.csv"),row.names=FALSE)
