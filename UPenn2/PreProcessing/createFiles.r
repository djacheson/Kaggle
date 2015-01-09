#This file reads in raw data in matlab format, produces a raw data object and does some initial
#  processing of the data to create features

#The data object is a list composed of the following elements:
#  $trialType : what type of trial (ictal / interictal)
#  $trialNum  : the clip number
#  $data      : chan X time matrix of ecog data (uV / time)
#  $length    : number of samples in the data
#  $freq      : sampling rate of the data
#  $channels  : vector with channel numbers
#  $time      : vector of time in seconds


setwd('/home/dan/Kaggle/UPenn2/')
source("../featureGeneration.R")
library(R.matlab)
library(data.table)
library(stringr)


DownsampleRaw <- function(data,newRate) {
  #Downsample raw data
  #INPUT:
  #  data    : data object
  #  newRate : new sampling rate of the data
  #
  #OUTPUT: downsampled data object
  #  data   : raw data object
  #
  
  if(round(data$freq) %% newRate > 0 ) {
    stop(cat(c("\nNew rate must be even multiple of old rate","\nOld rate is:",data$freq)))
  }
  
  dataDim <- dim(data$data)
  #Specify number of samples to keep
  N = data$freq / newRate
  seed <- c(TRUE,rep(FALSE,N-1))
  #create logical filter the length of second data dimension (time)
  filt <- rep(seed,ceiling(dataDim[[2]]/N))[1:dataDim[[2]]]
  #filter raw data matrix and time
  data$data <- data$data[,which(filt)]
  data$time <- data$time[which(filt)] #1:dim(data$data)[[2]] / newRate
  #log the new sampling rate
  data$freq <- as.numeric(newRate)
  data
}



readFile <- function(filename) {
  #function to read in a raw matlab file and create a data list object
  
  #INPUT:
  #  filename : path to file
  
  #OUTPUT:
  #  data object
  
  
  #File organization
  #[[1]] data
  #[[2]] data.length.sec
  #[[3]] sampling frequency
  #[[4]] channels
  #[[5]] sequence
  
  library(data.table)
  library(stringr)
  file <- readMat(filename)
  trialType <- ifelse(str_detect(names(file)[[1]],"interictal"),"interictal",ifelse(str_detect(names(file)[[1]],"preictal"),"preictal","test"))
  trialNum <- str_extract(names(file)[[1]],"\\d+")
  
  #print where we are
  print(paste(trialType,":",trialNum))
        

  file$trialType <- trialType
  file$trialNum <- trialNum
  #Extract data from the embedded list...
  file$data <- file[[1]][[1]]
  file$length <- file[[1]][[2]]
  file$freq<- round(file[[1]][[3]])
  file$channels <- seq(1:length(file[[1]][[4]]))
  file$time <- 1:dim(file[[1]][[1]])[[2]] / file$freq
  file[[1]] <- NULL #get rid of the embedded list

  return(file)
  
}


makeFeatures <- function(filename) {
  #Function to make an initial set of features using functions in featureGeneration.R
  #  reads in a raw file, downsamples the data, then creates the features
  
  #INPUT:
  #  filename  :  path to a raw matlab file
  
  file <- readFile(filename)
  fileInfo <- data.frame(trialType = file$trialType, trialNum = file$trialNum)
  
  #downsample if sampling rate >500
  if(file$freq > 500) {
      file <- DownsampleRaw(file,500)
  }
  
  #generate various features sets
  power <- getSpecPower_all(file) #spectral power
  SD <- getSD(file) #standard deviation
  cors <- getCors(file) #correlation across channels
  skew <- get_skew(file) #skew in data across time
  kurt <- get_kurtosis(file) #kurtosis of data across time
  autoCorVar <- get_autoCorVar(file,maxLag = 20) #variance in autocorrelation of signal
  canCorMean <- get_canCorMean(file,maxLag = 20) #mean autocorrelation in the signal
  
  
  cbind(fileInfo,power,SD,cors,skew,kurt,canCorVar,canCorMean)
}

#benchmarking
#start = proc.time()
#g <- makeFeatures(filename)
#print(proc.time() - start)

#----------------------------------------------
#Process all files
dirs <- list.dirs("../Data")
dirs <- dirs[2:length(dirs)]

all_files <- list()

#create big list of all files in the different directories
for (i in seq_along(dirs)) {
    all_files[[dirs[[i]]]]$preictal <- list.files(dirs[[i]], pattern = "*_preictal_*")
    all_files[[dirs[[i]]]]$interictal <- list.files(dirs[[i]], pattern = "*_interictal_*")
    all_files[[dirs[[i]]]]$test <- list.files(dirs[[i]], pattern = "*_test_*")
    
}

#read in files to a list data.tables for each trial type (ictal, interictal, test)
#make features for each and and save as an .rbin file for later loading
for (i in seq_along(dirs)) {
    start=proc.time()
    train <- list()
    test <- list()
    print(dirs[[i]])
    train$preictal <- lapply(all_files[[dirs[[i]]]]$preictal, function(x) makeFeatures(paste0(dirs[[i]],"/",x)))
    train$interictal  <- lapply(all_files[[dirs[[i]]]]$interictal, function(x) makeFeatures(paste0(dirs[[i]],"/",x)))

    train$preictal <- rbindlist(train$preictal)
    train$interictal <- rbindlist(train$interictal)
    save(train,file=paste0(dirs[[i]],"/train.rbin"))
    rm(train)
    
    test$test  <- lapply(all_files[[dirs[[i]]]]$test, function(x) makeFeatures(paste0(dirs[[i]],"/",x)))
    test$test <- rbindlist(test$test)
    save(test,file=paste0(dirs[[i]],"/test.rbin"))
    rm(test)
    print((proc.time()-start)[[3]])
}
