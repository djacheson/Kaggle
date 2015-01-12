#Script to run svmRadial models across subjects using different combinations of feature spaces tuning using ROC
#NOTES ON TRAINING:
# 1. Data is upsampled to deal with class imbalance
# 2. Different combinations of features can be added using 'feature_set' variable (see ../preProcessing/createFilesX.r files)
# 3. Given large feature space, variable importance of a previously trained model is used for dimensionality reduction
# 4. Data is centered, scaled and a spatial sign transform is applied prior to training
# 5. Different types of models can be run (e.g., random forest, glmnet, nnet, etc. simply by changing the model type in the call to train in caret)


setwd('/home/dan/Kaggle/UPenn2/')
library(caret)
library(data.table)
library(stringr)
library(kernlab)
#library(mda)


library(doMC)
registerDoMC(cores = 4)


dirs <- c("./Dog_1", "./Dog_2", "./Dog_3", "./Dog_4", "./Dog_5", "./Patient_1",  "./Patient_2")
# dirs <- list.dirs()
# dirs <- dirs[2:length(dirs)]
# dirs <- c(dirs[1:5],dirs[9:10]) #get rid of model

makeClipNames <- function(dir,testSet) {
  #Helper function to generate clip names for submission
  #
  #INPUT:
  #   dir  : path to directory
  #   testSet  : data object for test items
  
  name <- str_replace(dir,"[//.][//]","")
  
  #Get clip numbers from trial numbers in testSet$trialNum vector
  clipNums <- as.character(testSet$trialNum)
  clipNums <- sapply(clipNums, function(x) ifelse(nchar(x)==1,paste0("000",x),
                                                  ifelse(nchar(x)==2,paste0("00",x),
                                                         ifelse(nchar(x)==3, paste0("0",x),x))))
  paste0(name,"_test_segment_",clipNums,".mat")
}


classWeights <- function(factor_column) {
  #function to generate class weights, used for weighting cost in SVMs (but not used here)
  
  #INPUT:
  #   factor_column  :  vector with factor (i.e., categorical) data
  #
  lngth = length(factor_column)
  class_weights <- data.frame(class = unique(factor_column), 
                              total = sapply(unique(factor_column), function(x) lngth / sum(factor_column == x)))
  weights <- sapply(factor_column, function(x) class_weights[class_weights[[1]] == x,2])
  weights
}

#caret fit control used by all models
fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 1,
  classProbs=TRUE,
  summaryFunction = twoClassSummary,
  allowParallel = TRUE
)




#Loop through different feature sets and subjects generating models and predictions

#outer loop to go through different feature sets
#NOTE: If feature_set is set to an individual number, only one feature space is used, 
#       otherwise, different features spaces are combined

for(feature_set in 1:1) { 
  #create big list of all models and predictions
  all_models2 <- list()
  all_preds <- list()
  for (i in seq_along(dirs)) { #inner loop to go through each subject
    
    #--------------------------------------------
    #LOAD THE FEATURE SET TO BE USED FOR TRAINING
    #original feature set
    load(file=paste0(dirs[[i]],"/train.rbin"))
    load(file=paste0(dirs[[i]],"/test.rbin"))
    test_orig <- as.data.frame(test[[1]])
    train_orig <- as.data.frame(rbindlist(train))
    #     train <- train[,!(str_detect(names(train),"cor")),with=F]
    #     test <- test[,!(str_detect(names(test),"cor")),with=F]
    
    if(feaature_set==1) {
      #orig + frequency correlations
      features = "Orig+freqCors"
      #load previous model
      if(i==1) {
        load(file=paste0("./ModelSave/svmRadial_",features,"_upSamp_scale+spatialSign_ROC.rbin"))
      }
      
      #SELECT SUBSET OF FEATURES FROM PREVIOUS MODEL USING VARIABLE IMPORTANCE
      print(paste0("calculating variable importance for ", dirs[[i]]))
      var_imps <- varImp(all_models[[i]])
      #here... choose best using a cutoff of importance >30 (out of 100)
      vars2keep <- row.names(var_imps$importance)[var_imps$importance$preictal > 30]
      vars2keep <- c("trialType","trialNum",vars2keep)
      
      load(file=paste0(dirs[[i]],"/train_freqCorrs.rbin"))
      load(file=paste0(dirs[[i]],"/test_freqCorrs.rbin"))
      test_freq <- as.data.frame(test[[1]])
      train_freq <- as.data.frame(rbindlist(train))
      
      #here, filter features based on variable importance above
      test_freq <- test_freq[,names(test_freq) %in% vars2keep]
      train_freq <- train_freq[,names(train_freq) %in% vars2keep]
      test_orig <- test_orig[,names(test_orig) %in% vars2keep]
      train_orig <- train_orig[,names(train_orig) %in% vars2keep]
      
      train <- cbind(train_orig,train_freq) #get rid of repeat names
      test <- cbind(test_orig,test_freq) 
    }
    #Add other feature spaces to the original set of features
    if(j > 1) {
      #only keep trialType and trialNum
      train_orig <- as.data.frame(train_orig[,1:2,with=F])
      test_orig <- as.data.frame(test_orig[,1:2, with=F])
      if(feaature_set== 2) {
        #add coherence between channels 
        features = "coherence"
        #load previous model
        if(i==1) {
          load(file=paste0("./ModelSave/svmRadial_",features,"_scale+spatialSign_upSamp_ROC.rbin"))
        }
        
        load(file = paste0(dirs[[i]],"/train_coherence.rbin"))
        load(file = paste0(dirs[[i]],"/test_coherence.rbin"))
        test_coh <- test[[1]]
        train_coh <- rbindlist(train)
        train_coh <- train_coh[,unique(names(train_coh)),with=F] #get rid of repeat names
        test_coh <- test_coh[,unique(names(train_coh)),with=F] #get rid of repeat names
        train <- cbind(train_orig,train_coh) 
        test <- cbind(test_orig,test_coh) }
      if(feaature_set==3) {
        #add correlations and mutual information between channels
        features = "cor+MI"
        #load previous model
        if(i==1) {
          load(file=paste0("./ModelSave/svmRadial_",features,"_scale+spatialSign_upSamp_ROC.rbin"))
        }
        load(file=paste0(dirs[[i]],"/train_corr+MI.rbin"))
        load(file=paste0(dirs[[i]],"/test_corr+MI.rbin"))
        test_cor <- test[[1]]
        train_cor <- rbindlist(train)
        train <- cbind(train_orig,train_cor) 
        test <- cbind(test_orig,test_cor) }
      if(feaature_set==4) {
        #add first derivative features
        features = "derivs"
        #load previous model
        if(i==1) {
          load(file=paste0("./ModelSave/svmRadial_",features,"_scale+spatialSign_upSamp_ROC.rbin"))
        }
        load(file = paste0(dirs[[i]],"/train_derivs.rbin"))
        load(file = paste0(dirs[[i]],"/test_derivs.rbin"))
        test_der <- as.data.frame(test[[1]])
        train_der <- as.data.frame(rbindlist(train))
        train <- cbind(train_orig,train_der) 
        test <- cbind(test_orig,test_der) }
      if(feaature_set== 5) {
        #expanded derivative features, keeping only second derivative featuresS
        features = "derivs2_2D"
        #load previous model
        if(i==1) {
          load(file=paste0("./ModelSave/svmRadial_",features,"_scale+spatialSign_upSamp_ROC.rbin"))
        }
        load(file = paste0(dirs[[i]],"/train_derivs2.rbin"))
        load(file = paste0(dirs[[i]],"/test_derivs2.rbin"))
        test_der <- as.data.frame(test[[1]])
        train_der <- as.data.frame(rbindlist(train))
        
        #only keep the second derivative features
        var2keep <- str_detect(names(train_der),"2D")
        test_der <- test_der[,var2keep]
        train_der <- train_der[,var2keep]
        
        train <- cbind(train_orig,train_der)
        test  <- cbind(test_orig,test_der) }
    }
    
    
    
    #--------------------------------------------
    #SET UP DATA FOR TRAINING, TRAIN AND PREDICT
    
    #Upsample to deal with class imbalance
    upSamp <- upSample(train,as.factor(train$trialType))
    #upSamp <- as.data.frame(train)
    
    #randomize the upSamp
    upSamp <- upSamp[sample(1:nrow(upSamp)),]
    upSamp$Class <- NULL
    
    #Scale, center an spatial sign transformations
    scale <- preProcess(upSamp[,3:ncol(upSamp)],
                        method = c("center","scale","spatialSign")
    )
    upSamp_scale <- predict(scale, upSamp[,3:ncol(upSamp)])
    test_scale <- predict(scale, test[,3:ncol(test)])
    upSamp <- cbind(upSamp[,1:2],upSamp_scale)
    
    #clean up the space
    rm(train_orig, test_origtest2, train2, train_coh, train_cor, train_freq, test_coh, test_cor, test_freq, train_der, test_der)
    
    #Estimate sigma and set tuning control for all models
    sigma <- sigest(as.matrix(upSamp[,3:ncol(upSamp)]), frac = 0.75)
    tuneControl <-  expand.grid(sigma = sigma[[3]],
                                C = 10000 / 10^(0:7)#,
                                #Weight = c(seq(1,61,by=10))
    )
    
    
    start = proc.time()
    all_models2[[i]] <- train(x = upSamp[,3:ncol(upSamp)],
                             y = upSamp$trialType,#y = as.factor(ifelse(upSamp$trialType == "preictal",1,0)), #this might need to be changed depending on model type
                             method = 'svmRadial', #change method if you want to train using a different model type
                             trControl = fitControl,
                             tuneGrid = tuneControl,
                             metric = "ROC"
    )
    
    print(paste0(dirs[[i]],": ",(proc.time()-start)[[3]]))
    print(all_models2[[i]])
    plot(all_models2[[i]])
    
    all_preds[[i]] <- data.table(clip = makeClipNames(dirs[[i]],test),
                                 preictal = predict(all_models2[[i]],test_scale,type="prob")[[1]]
    )
  }#end of inner loop across subjects
  
  
  all_models <- all_models2
  save(all_models,file=paste0("./ModelSave/svmRadial_",features,"_scale+spatialSign_upSamp_ROC.varImp.rbin"))
  
  all_preds_orig = all_preds
  
  #Combine all predictions and save them in original format
  all_preds <- rbindlist(all_preds)
  write.csv(all_preds,file=paste0("./ModelPredictions/svmRadial_",features,"_scale+spatialSign_upSamp_ROC.varImp.csv"),row.names=F)
  
  #scale all preds so highest prob = 1
  all_preds = all_preds_orig
  for(i in seq_along(length(all_preds))) {
    all_preds[[i]]$preictal = all_preds[[i]]$preictal / max(all_preds[[i]]$preictal)
  }
  all_preds <- rbindlist(all_preds)
  write.csv(all_preds,file=paste0("./ModelPredictions/svmRadial_",features,"_scale+spatialSign_upSamp_ROC.maxScale.varImp.csv"),row.names=F)

}#end of outer loop across features

