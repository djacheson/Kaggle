#SPECTRAL COHERENCE

#This function was added to expand the feature space to include spectral coherence
#  It's similar to the steps performed in createFiles.r
#  but it does it for a single directory, so is thus easily parallelized across cores / nodes in a cluster

# It's ultimately called by the createFilesBatch.R scripts

createFiles_oneDir <- function(dirs) {
    #setwd('/home/dan/Kaggle/UPenn2/')
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
        
        #   if(trialType=="preictal") {
        #     latency = file$latency
        #   } else {
        #     latency = NA
        #   }
        # 
        file$trialType <- trialType
        file$trialNum <- trialNum
        #Extract data from the embedded list...
        file$data <- file[[1]][[1]]
        file$length <- file[[1]][[2]]
        file$freq<- round(file[[1]][[3]])
        file$channels <- seq(1:length(file[[1]][[4]]))
        file$time <- 1:dim(file[[1]][[1]])[[2]] / file$freq
        file[[1]] <- NULL #get rid of the embedded list
        
        file
        
    }
    
    
    makeFeatures <- function(filename) {
        #here, generating spectral coherence instead of a big set of features
        file <- readFile(filename)
        fileInfo <- data.frame(trialType = file$trialType, trialNum = file$trialNum)
        
        if(file$freq > 500) {
            file <- DownsampleRaw(file,500)
        }
        
        coherence <- getCoherence_all(file)
        
        
        coherence
    }
    
    # start = proc.time()
    # g <- makeFeatures(filename)
    # print(proc.time() - start)
    
    #----------------------------------------------
    # dirs <- list.dirs()
    # dirs <- dirs[2:length(dirs)]
    
    all_files <- list()
    
    #create big list of all files
    for (i in seq_along(dirs)) {
        all_files[[dirs]]$preictal <- list.files(dirs, pattern = "*_preictal_*")
        all_files[[dirs]]$interictal <- list.files(dirs, pattern = "*_interictal_*")
        all_files[[dirs]]$test <- list.files(dirs, pattern = "*_test_*")
        
    }
    
    #read in files to a list data.tables for each data type and save as an .rbin file for later loading
    for (i in seq_along(dirs)) {
        #for (i in seq(from=8, to=12, by=1)) {
        start=proc.time()
        train <- list()
        test <- list()
        print(dirs)
        train$preictal <- lapply(all_files[[dirs]]$preictal, function(x) makeFeatures(paste0(dirs,"/",x)))
        train$interictal  <- lapply(all_files[[dirs]]$interictal, function(x) makeFeatures(paste0(dirs,"/",x)))
        
        train$preictal <- rbindlist(train$preictal)
        train$interictal <- rbindlist(train$interictal)
        save(train,file=paste0(dirs,"/train_coherence.rbin"))
        rm(train)
        
        test$test  <- lapply(all_files[[dirs]]$test, function(x) makeFeatures(paste0(dirs,"/",x)))
        test$test <- rbindlist(test$test)
        save(test,file=paste0(dirs,"/test_coherence.rbin"))
        rm(test)
        print((proc.time()-start)[[3]])
    }
}
