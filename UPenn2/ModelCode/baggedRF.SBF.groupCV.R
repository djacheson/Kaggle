#Bagged random forests using selection by filtering (SBF) to reduce the features space

#NOTES:
# 1. A big feature space is build by combining features generated during preprocessing 
#     (see ../preProcessing/createFilesX.r files)
# 2. CV folds are determined using CV_splits.csv file
# 3. SBF filters variables according to variable importance in CV folds
# 4. Individual models are saved per subject. The files can get large depending on the number of bagged models


setwd('/home/dan/Kaggle/UPenn2/')
library(caret)
library(data.table)
library(stringr)



library(doMC)
registerDoMC(cores = 4)


dirs <- c("./Dog_1", "./Dog_2", "./Dog_3", "./Dog_4", "./Dog_5", "./Patient_1",  "./Patient_2")
# dirs <- list.dirs()
# dirs <- dirs[2:length(dirs)]
# dirs <- c(dirs[1:5],dirs[9:10]) #get rid of model

makeClipNames <- function(dir,testSet) {
  name <- str_replace(dir,"[//.][//]","")
  clipNums <- as.character(testSet$trialNum)
  clipNums <- sapply(clipNums, function(x) ifelse(nchar(x)==1,paste0("000",x),
                                                  ifelse(nchar(x)==2,paste0("00",x),
                                                         ifelse(nchar(x)==3, paste0("0",x),x))))
  paste0(name,"_test_segment_",clipNums,".mat")
}


#load CV splits file
CV_splits <- fread("./CV_splits.csv")
#create big list of all files
all_models <- list()
all_preds <- list()
for (i in seq_along(dirs)) { #loop to go through each subject
  curr_subj <- str_replace(dirs[[i]],"[//.][//]","")
  
  #Get CV information
  curr_CV <- CV_splits[File == curr_subj,]
  #fix some ordering issues so the CV and train set line up
  curr_CV$trialNum <- as.numeric(curr_CV$trialNum)
  setorderv(curr_CV, c("trialType","trialNum"),c(-1,1))
  curr_CV$rowNums <- 1:nrow(curr_CV)
  
  #list with CV folds to be used in fit control below
  index <- list()
  for(fold in seq_along(unique(curr_CV$CV))) {
    index[[fold]] = curr_CV[CV == fold]$rowNums
  }
  
  #specify fit control using the CV folds defined above
  # Selection by filter here uses cross validation performance to select best variables
  
  fitControl <- sbfControl(
    functions = treebagSBF,
    method = "cv",
    number = max(curr_CV$CV),
    index = index,
    saveDetails = TRUE,
    #     classProbs=TRUE,
    #     summaryFunction = twoClassSummary,
    allowParallel = TRUE
  )
  
  #------------------------------------------------------------------------
  #LOAD FEATURE SPACES AND COMBINE THEM TOGETHER... LATER FILTERED USING SBF
  #NOTE: For details about different features used, see ../preProcessing/createFilesX.r files
  
  #Load original train / test file 
  load(file=paste0(dirs[[i]],"/train.rbin"))
  load(file=paste0(dirs[[i]],"/test.rbin"))
  test_orig <- as.data.frame(test[[1]])
  train_orig <- as.data.frame(rbindlist(train))
  train_orig <- train_orig[,!(str_detect(names(train_orig),"cor"))]
  test_orig <- test_orig[,!(str_detect(names(test_orig),"cor"))]
  #         train_orig <- train_orig[,!(str_detect(names(train_orig),"autoCor"))] #get rid of autoCor variable
  #          test_orig <- test_orig[,!(str_detect(names(test_orig),"autoCor"))]
  
  #first and second derivative features
  load(file = paste0(dirs[[i]],"/train_derivs2.rbin"))
  load(file = paste0(dirs[[i]],"/test_derivs2.rbin"))
  test_der <- as.data.frame(test[[1]])
  train_der <- as.data.frame(rbindlist(train))
  
  #correlation and mutual information
  features = "cor+MI"
  load(file=paste0(dirs[[i]],"/train_corr+MI.rbin"))
  load(file=paste0(dirs[[i]],"/test_corr+MI.rbin"))
  test_cor <- as.data.frame(test[[1]])
  train_cor <- as.data.frame(rbindlist(train))
  
  #train <- cbind(train_orig,train_der, train_cor)
  test  <- cbind(test_orig,test_der, test_cor) 
  
  
  #------------------------------------------------------------
  #PREPROCESSING AND TRAINING
  
  #upSample, center scale and spatial sign trainsform the data
  upSamp <- as.data.frame(train)
  scale <- preProcess(upSamp[,3:ncol(upSamp)],
                      method = c("center","scale","spatialSign")
  )
  upSamp_scale <- predict(scale, upSamp[,3:ncol(upSamp)])
  test_scale <- predict(scale, test[,3:ncol(test),with=F])
  upSamp <- cbind(upSamp[,1:2],upSamp_scale)
  
  rm(train_orig, test_orig, test2, train2, train_coh, train_cor, train_freq, test_coh, test_cor, test_freq, train_der, test_der, train_der2, test_der2)
  gc()
  
  start = proc.time()
  
  model <-  sbf(x = upSamp[,3:ncol(upSamp)],
                y = upSamp$trialType,
                sbfControl = fitControl#,
                ntree = 1000
                metric= "Kappa"
  )
    
  
  #model
  #plot(all_models[[i]])
  
  #   
  all_preds[[i]] <- data.table(clip = makeClipNames(dirs[[i]],test),
                               preictal = predict(model,test[,3:ncol(test)],type="prob")[[2]]
  )
  save(model,file=paste0("./ModelSave/bagTree/",curr_subj,"treebagSBF_groupCV.rbin"))
  rm(model)
  gc()
  
  print(paste0(dirs[[i]],": ",(proc.time()-start)[[3]]))
}# end of inner loop

#save(all_models,file=paste0("./ModelSave/treebagSBF_groupCV.rbin"))

#Combine predictions
all_preds_orig = all_preds

all_preds <- rbindlist(all_preds)
write.csv(all_preds,file=paste0("./ModelPredictions/treebagSBF_groupCV.csv"),row.names=F)

# scale all preds so highest prob = 1
all_preds = all_preds_orig
for(i in seq_along(length(all_preds))) {
  all_preds[[i]]$preictal = all_preds[[i]]$preictal / max(all_preds[[i]]$preictal)
}

all_preds <- rbindlist(all_preds)
write.csv(all_preds,file=paste0("./ModelPredictions/treebagSBF_groupCV.maxScale.csv"),row.names=F)

)