#Script to load raw matlab data and convert it into R list objects for later processing
# Downsamples the raw data along the way

setwd('/home/dan/Kaggle/UPenn_Seizure/clips')
library(R.matlab)
library(data.table)
library(stringr)


DownsampleRaw <- function(data,newRate) {
  #Downsample raw data
  #Assumes input of a data object (i.e., list) including:
  #  freq :    sampling rate
  #  data :    channel X time matrix
  
  #INPUT:
  #  data    : data list
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
  data$time <- seq(0,1 - 1 / newRate, by = 1 / newRate)
  #log the new sampling rate
  data$freq <- as.numeric(newRate)
  data
}



readFile <- function(filename) {
  #Function to read in raw Matlab data and output in tabular format
  # Extracts relevant data values and downsamples the data
  # Designed to be used within a loop / apply call over directories
  
  #INPUT:
  #  filename :   string of filename to process
  
  #OUTPUT
  #  flatList :   table where rows represent channels and columns are variables / data of interest (see below)
  
  library(data.table)
  library(stringr)
  #read matlab data in using R.matlab
  file <- readMat(filename)
  #determine trialtype and number using filename
  trialType <- ifelse(str_detect(filename,"interictal"),"interictal",ifelse(str_detect(filename,"ictal"),"ictal","test"))
  trialNum <- str_extract(filename,"\\d+\\.mat")
  trialNum <- str_replace(trialNum,"\\.mat","")
  
  #print where we are
  print(paste(trialType,":",trialNum))
  
  #If it's an ictal trial, we want to include the latency
  if(trialType=="ictal") {
    latency = file$latency
  } else {
    latency = NA
  }
  
  
  file$freq<- round(file$freq)
  #downsample the data to 400 Hz or 500 Hz
  if(is.character(try(file <- DownsampleRaw(file, 400),silent=TRUE))) {
    file <- DownsampleRaw(file, 500)
  } 
  #rename channels
  file$channels <- seq(1:length(file$channels))
  
  #create a big list of data
  flatList <- lapply(seq(1:length(file$channels)), function(i) data.table(trialType=trialType, # type of trial: ictal, intei
                                                                          trialNum = trialNum, # trial number
                                                                          freq=file$freq, # sampling rate
                                                                          channel=file$channels[[i]], # channel number
                                                                          latency=latency, # latency in seconds
                                                                          time=file$time, # vector of time
                                                                          data=file$data[i,])) # raw data for a channel (voltage over time)
  fileTable <- rbindlist(flatList)
  #output
  fileTable
  
}



#----------------------------------------------
dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]

all_files <- list()

#create big list of all the files for each directory (Dog / Patient)
for (i in seq_along(dirs)) {
    all_files[[dirs[[i]]]]$ictal <- list.files(dirs[[i]], pattern = "*_ictal_*")
    all_files[[dirs[[i]]]]$interictal <- list.files(dirs[[i]], pattern = "*_interictal_*")
    all_files[[dirs[[i]]]]$test <- list.files(dirs[[i]], pattern = "*_test_*")
    
}

#read in files to a list of data.tables for each trial type and save as an .rbin file for later processing
#Note: This drastically reduces the file size and loading time 
for (i in seq_along(dirs)) {
    train <- list()
    test <- list()
    print(dirs[[i]])
    train$ictal <- lapply(all_files[[dirs[[i]]]]$ictal, function(x) readFile(paste0(dirs[[i]],"/",x)))
    train$interictal  <- lapply(all_files[[dirs[[i]]]]$interictal, function(x) readFile(paste0(dirs[[i]],"/",x)))

    train$ictal <- rbindlist(train$ictal)
    train$interictal <- rbindlist(train$interictal)
    save(train,file=paste0(dirs[[i]],"/train.rbin"))
    rm(train)
    
    test$test  <- lapply(all_files[[dirs[[i]]]]$test, function(x) readFile(paste0(dirs[[i]],"/",x)))
    test$test <- rbindlist(test$test)
    save(test,file=paste0(dirs[[i]],"/test.rbin"))
    rm(test)
}
