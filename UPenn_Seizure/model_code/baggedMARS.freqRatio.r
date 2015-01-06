setwd("/home/dan/Kaggle/UPenn_Seizure/clips")

library(data.table)
library(caret)
#library(doMC)
#registerDoMC()
library(stringr)

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]

#list object to hold all predictions
allPreds <- list()
trainPreds <- list()


for(i in seq_along(dirs))  {
  start = proc.time()
  load(file=paste0(dirs[[i]],"/train+test.freqRatios.wide.rbin"))
  
  currSub <- str_replace(dirs[[i]],"[/ .][/ /]","")
  
  clips = paste0(currSub,"_test_segment_",test.wide$trialNum,".mat")
  #Create test file
  test.x <- test.wide[,c("trialType","trialNum") := NULL]

  
  #**************
  #ictal or not 
  #**************
  train.copy <- copy(train.wide)
  train.y <- as.factor(ifelse(train.wide$trialType=="ictal","X1","X0"))

  train.trialNum <- train.copy$trialNum

  train.x <- as.matrix(train.copy[,c("trialType","trialNum","latency_15") := NULL])
  rm(train.copy)
  
  #upsampling
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  
  shuffle <- sample(1:nrow(upSamp), nrow(upSamp))
  upSamp <- upSamp[shuffle,]
  
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y
  
  
#   #Mars
  
  ictal.mars <- bagEarth(x=up.x, 
                         y = up.y, 
                         B = 25, 
                         glm=list(family=binomial), 
                         degree = 5, 
                         keepX = FALSE, 
                         minspan=1) #- patient 8 only
  
  #Make predictions
  ictalProb <- round(predict(ictal.mars,as.matrix(test.x), type = "response"),2)
  trainProb <- round(predict(ictal.mars,as.matrix(train.x), type = "response"),2)
  #ictalClass <- ifelse(ictalProb < 0.5, "interictal","ictal")
  
  #save it and save memory :)
  save(list=c("ictal.mars"), file=paste0("../model_save/bagMARS/",currSub,"_ictalMARS.freqRatio.rbin"))
  rm(ictal.mars)
  gc()
  #add to preds data.table

  preds <- data.table(clip=clips, seizure=as.numeric(ictalProb))
  
  train.preds.ictal <- data.table(subject = currSub, 
                          trialNum = train.trialNum, 
                          bagMARS.freq = trainProb, 
                          target=train.y)
  
  #---------------------------------------------------------------------------------------------------
  #*************
  #Early Latency
  #************
  train.copy <- copy(train.wide)
  #lat15 <- train.copy[latency_15 >= 0]
  train.trialNum <- train.copy$trialNum

  train.y <- ifelse(train.copy$latency_15 == 1,"X1","X0")
  train.y <- as.factor(ifelse(is.na(train.y),"X0",train.y))

  train.x <- as.matrix(train.copy[,c("trialType","trialNum","latency_15") := NULL])
  
  #upsampling
  upSamp <- upSample(train.x,train.y, yname = "train.y")
  shuffle <- sample(1:nrow(upSamp), nrow(upSamp))
  upSamp <- upSamp[shuffle,]
  
  up.x <- upSamp
  up.x$train.y <- NULL
  up.y <- upSamp$train.y

  #MODEL
  early.mars <- bagEarth(x=up.x, 
                         y = up.y, 
                         B = 25, 
                         glm=list(family=binomial), 
                         degree = 5, 
                         keepX = FALSE, 
                         minspan=1) 
  
  #Make predictions
  earlyProb <- round(predict(early.mars,as.matrix(test.x), type = "response"),2)
  trainProb <- round(predict(early.mars,as.matrix(train.x), type = "response"),2)
  
  save(list=c("early.mars"), file=paste0("../model_save/bagMARS/",currSub,"_earlyMARS.freqRatio.rbin"))
  rm(early.mars)
  gc()

  preds[,early := earlyProb]

  train.preds.early <- data.table(subject = currSub, 
                                trialNum = train.trialNum, 
                                bagMARS.freq = trainProb, 
                                target=train.y)

  
  allPreds[[i]] <- preds
  
  trainPreds$ictal[[i]] <- train.preds.ictal
  trainPreds$early[[i]] <- train.preds.early
  #save(list=c("allPreds","allmodel_save"), file=paste0("../model_save/bagMARS1.rbin"))
  print(paste0("Done with subject ", dirs[[i]]," in: ", (proc.time()-start)[[3]]))
  gc()
  save(allPreds, file=paste0("../model_save/bagMARS/allPreds.freqRatios.rbin"))
  save(trainPreds, file=paste0("../model_save/bagMARS/trainPreds.freqRatios.rbin"))
  
}


allPredComb <- rbindlist(allPreds)
write.csv(allPredComb,file=paste0("../model_save/bagMARS.preds.freqRatios.csv"),row.names=FALSE)

trainPreds$ictal <- rbindlist(trainPreds$ictal)
trainPreds$early <- rbindlist(trainPreds$early)
write.csv(trainPreds$ictal,file=paste0("../model_save/trainPreds/bagMARS.freqRatio.TRAINpreds.ictal.csv"),row.names=FALSE)
write.csv(trainPreds$early,file=paste0("../model_save/trainPreds/bagMARS.freqRatio.TRAINpreds.early.csv"),row.names=FALSE)
