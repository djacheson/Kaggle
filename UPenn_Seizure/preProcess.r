# This file creates the data to be used for training and testing in caret
#   It generates various features using the functions in 'featureGeneration.R'
#   and ultimately transforms the data from long to wide format 
#   resulting in a tidy data.table with 1 row per clip

setwd("/home/dan/Kaggle/UPenn_Seizure/clips")

source("../featureGeneration.r")
library(data.table)

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]

#Loop through all subjects, training and testing files and create features from raw data
#Features include;
#  overall variance 
#  minMax difference 
#  extreme values (difference in 5th and 95th percentile for rectified, mean subtracted data)
#  mean global spectral power
#  mean spectral power at different frequency bands

for (i in seq_along(dirs)) {
  
  print("Starting train file")
  start = proc.time()
  
  #Training file
  load(file=paste0(dirs[[i]],"/train.rbin"))
  train <- rbindlist(train) 
  
  #code variance of each trial
  train[,variance := var(data), by=list(trialType,trialNum,channel)]
  
  #minMax and extreme values (currently 5% and 95%)
  train[,minMax := minMax(data), by=list(trialType,trialNum,channel)]
  train[,extremeVals := extremeVals(data), by=list(trialType,trialNum,channel)]
  
  #Add frequency info
  sampleRate = as.numeric(unique(train$freq))
  
  #Global Power Across Frequencies
  train[,globalPower := getGlobalPower(data, sampleRate),by=list(trialType,trialNum,channel)]
  
  #Average Power at Individual Frequencies
  train[,c("highGamma","lowGamma","beta","alpha","theta","delta") := getSpecPower_all(data,sampleRate), 
        by=list(trialType,trialNum,channel)]
  
  #code latency_15 variable, coding for if latency <=15
  train[,latency_15 := ifelse(latency <= 15,1,0)]
  
  #save it
  save(train,file=paste0(dirs[[i]],"/train+features.rbin"))
  print(paste0("Done with trainfile. That took", (proc.time()-start)[[3]]))
  rm(train)
  
  #--------------------------------
  #Test file
  print("Starting test file")
  start = Sys.time()
  
  load(file=paste0(dirs[[i]],"/test.rbin"))
  test <- rbindlist(test)
  #code variance of each trial
  test[,variance := var(data), by=list(trialNum,channel)]
  #minMax and extreme values (currently 5% and 95%)
  test[,minMax := minMax(data), by=list(trialType,trialNum,channel)]
  test[,extremeVals := extremeVals(data), by=list(trialType,trialNum,channel)]
  
  # Mean Spectral Power
  sampleRate = as.numeric(unique(test$freq))
  test[,globalPower := getGlobalPower(data, sampleRate),by=list(trialType,trialNum,channel)]
  test[,c("highGamma","lowGamma","beta","alpha","theta","delta") := getSpecPower_all(data,sampleRate), 
       by=list(trialType,trialNum,channel)]
  print(paste0("Done with test file. That took", Sys.time()-start))
  save(test,file=paste0(dirs[[i]],"/test+features.rbin"))
  rm(test)
}




# Create table to display some summary stats for different frequency bands
# train[,
#       j = list(highGamma = mean(highGamma),
#                lowGamma = mean(lowGamma),
#                beta = mean(beta),
#                alpha = mean(alpha),
#                theta=mean(theta)
#       ),
#       by=latency_15]



combineFeaturesTrain <- function(DT) {
  #Function that combines features from a data table to make them useable with caret.
  # Reshapes data from long to wide format for a bunch of variables, 
  # then combining them into one data.table at the end
  
  
  #Mean alpha power
  # Long to wide
  alpha = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'alpha', mean)
  
  # get channel names to column names can be modified accordingly
  colNames <- names(alpha)
  colNames <- colNames[!(colNames %in% c("trialType","trialNum","latency_15"))]
  
  # add alpha_ in front of each channel name
  setnames(alpha,colNames,paste0("alpha_",colNames))
  setkeyv(alpha,c("trialType","trialNum","latency_15"))
  
  #same for other frequencies, and variance
  delta = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'delta', mean)
  setnames(delta,colNames,paste0("delta_",colNames))
  setkeyv(delta,c("trialType","trialNum","latency_15"))
  
  theta = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'theta', mean)
  setnames(theta,colNames,paste0("theta_",colNames))
  setkeyv(theta,c("trialType","trialNum","latency_15"))
  
  beta = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'beta', mean)
  setnames(beta,colNames,paste0("beta_",colNames))
  setkeyv(beta,c("trialType","trialNum","latency_15"))
  
  lowGamma = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'lowGamma', mean)
  setnames(lowGamma,colNames,paste0("lowGamma_",colNames))
  setkeyv(lowGamma,c("trialType","trialNum","latency_15"))
  
  highGamma = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'highGamma', mean)
  setnames(highGamma,colNames,paste0("highGamma_",colNames))
  setkeyv(highGamma,c("trialType","trialNum","latency_15"))
  
  variance = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'variance', mean)
  setnames(variance,colNames,paste0("variance_",colNames))
  setkeyv(variance,c("trialType","trialNum","latency_15"))
  
  minMax = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'minMax', mean)
  setnames(minMax,colNames,paste0("minMax_",colNames))
  setkeyv(minMax,c("trialType","trialNum","latency_15"))
  
  extremeVals = dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = 'extremeVals', mean)
  setnames(extremeVals,colNames,paste0("extremeVals_",colNames))
  setkeyv(extremeVals,c("trialType","trialNum","latency_15"))
  
  # Combine all columns into one data.table
  wideTable <- Reduce(merge,list(delta,theta,alpha,beta,lowGamma,highGamma,variance,minMax,extremeVals))
  wideTable
}

combineFeaturesTest<- function(DT) {
  #Function that combines features from a data table to make them useable with caret/other ML functions.
  # Reshapes data from long to wide format for a bunch of variables, 
  # then combining them into one data.table at the end
  
  
  #make lists for training
  alpha = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'alpha', mean)
  
  #get channel names so this can be changed for each
  colNames <- names(alpha)
  colNames <- colNames[!(colNames %in% c("trialType","trialNum"))]
  
  setnames(alpha,colNames,paste0("alpha_",colNames))
  setkeyv(alpha,c("trialType","trialNum"))
  
  #same for other frequencies, and variance
  delta = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'delta', mean)
  setnames(delta,colNames,paste0("delta_",colNames))
  setkeyv(delta,c("trialType","trialNum"))
  
  theta = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'theta', mean)
  setnames(theta,colNames,paste0("theta_",colNames))
  setkeyv(theta,c("trialType","trialNum"))
  
  beta = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'beta', mean)
  setnames(beta,colNames,paste0("beta_",colNames))
  setkeyv(beta,c("trialType","trialNum"))
  
  lowGamma = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'lowGamma', mean)
  setnames(lowGamma,colNames,paste0("lowGamma_",colNames))
  setkeyv(lowGamma,c("trialType","trialNum"))
  
  highGamma = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'highGamma', mean)
  setnames(highGamma,colNames,paste0("highGamma_",colNames))
  setkeyv(highGamma,c("trialType","trialNum"))
  
  variance = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'variance', mean)
  setnames(variance,colNames,paste0("variance_",colNames))
  setkeyv(variance,c("trialType","trialNum"))
  
  minMax = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'minMax', mean)
  setnames(minMax,colNames,paste0("minMax_",colNames))
  setkeyv(minMax,c("trialType","trialNum"))
  
  extremeVals = dcast.data.table(DT, trialType + trialNum  ~ channel, value.var = 'extremeVals', mean)
  setnames(extremeVals,colNames,paste0("extremeVals_",colNames))
  setkeyv(extremeVals,c("trialType","trialNum"))
  
  wideTable <- Reduce(merge,list(delta,theta,alpha,beta,lowGamma,highGamma,variance,minMax,extremeVals))
  wideTable
}





#Big Loop to make feature sets to be used with caret
for (i in seq_along(dirs)) {
  #train
  print(paste0("Starting subject:", dirs[[i]]))
  start = proc.time()
  load(file=paste0(dirs[[i]],"/train+features.rbin"))
  train.wide <- combineFeaturesTrain(train)
  rm(train)
  
  #test
  load(paste0(dirs[[i]],"/test+features.rbin"))
  test.wide <- combineFeaturesTest(test)
  rm(test)
  
  save(list=c("train.wide","test.wide"),file=paste0(dirs[[i]],"/train+test.wide.rbin"))
  print(paste0("Done: ", (proc.time()-start)[[3]]))
}
