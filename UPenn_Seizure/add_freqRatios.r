#This file added variables for the ratio of spectral power of different frequency bands
#  all possible combinations were added
#  These features were then added to the feature space including mean spectral power, variance, min-max and extreme values


library(data.table)

freqRatio <- function(delta,theta,alpha,beta,lowGamma,highGamma, combo="dt") {
    #function to calculate the ratio of mean spectral powers for different frequency bands
  
    #INPUT:
    #  mean spectral power for delta, theta, alpha, beta, lowGamm and highGamma
    #  combo : character string of ratio desired
  
    #OUTPUT:
    #  ratio of mean spectral power
  
    if(combo=="dt"){ ratio <- mean(delta/theta)} #dt = delta / theta
    if(combo=="da"){ ratio <- mean(delta/alpha)} # da = delta / alpha
    if(combo=="db"){ ratio <- mean(delta/beta)}  # etc. 
    if(combo=="dLg"){ ratio <- mean(delta/lowGamma)}
    if(combo=="dHg"){ ratio <- mean(delta/highGamma)}
    if(combo=="ta"){ ratio <- mean(theta/alpha)}
    if(combo=="tb"){ ratio <- mean(theta/beta)}
    if(combo=="tLg"){ ratio <- mean(theta/lowGamma)}
    if(combo=="tHg"){ ratio <- mean(theta/highGamma)}
    if(combo=="ab"){ ratio <- mean(alpha/beta)}
    if(combo=="aLg"){ ratio <- mean(alpha/lowGamma)}
    if(combo=="aHg"){ ratio <- mean(alpha/highGamma)}
    if(combo=="bLg"){ ratio <- mean(beta/lowGamma)}
    if(combo=="bHg"){ ratio <- mean(beta/highGamma)}
    if(combo=="LgHg"){ ratio <- mean(lowGamma/highGamma)}
    return(ratio)
}

combineFeaturesTrain <- function(DT,wideTable) {
  #Function that combines new features with a previously wide data table
  # in this instance, combining the different frequency ratios generated using freqRatio function
  
  #INPUT:
  #  DT : long data table with different frequency ratios
  #  wideTable : previously generated, wide data.table with features such as max-min, extreme values, etc. 
  
  #OUTPUT:
  #  wide table with all features added

  setkeyv(DT,c("trialType","trialNum","latency_15"))
  vars <- c("dt","da","db","dLg","dHg","ta","tb","tLg","tHg","ab","aLg","aHg","bLg","bHg","LgHg")
  #make lists for training
  
  #loop along all the different frequency ratios, going long to wide then combining with wideTable
  for (i in seq_along(vars)) {
    #current data.table corresponds to one of the frequency ratios in vars
    # long to wide
    currDT <- dcast.data.table(DT, trialType + trialNum + latency_15 ~ channel, value.var = vars[[i]], mean)
  
    #add column names corresponding to frequency ratio and channel number
    if(i==1) {
      colNames <- names(currDT)
      colNames <- colNames[!(colNames %in% c("trialType","trialNum","latency_15"))]
    }
    
    setnames(currDT, colNames, paste0(vars[[i]],"_",colNames))
    setkeyv(currDT, c("trialType","trialNum","latency_15"))
    
    #combine with wideTable
    wideTable <- Reduce(merge,list(wideTable,currDT))
    rm(currDT)
  
  }
  
  wideTable

}

combineFeaturesTest<- function(DT,wideTable) {
  #Same as combineFeaturesTrain but for test files
  setkeyv(DT,c("trialType","trialNum"))
  vars <- c("dt","da","db","dLg","dHg","ta","tb","tLg","tHg","ab","aLg","aHg","bLg","bHg","LgHg")
  #make lists for training
  
  for (i in seq_along(vars)) {
    currDT <- dcast.data.table(DT, trialType + trialNum ~ channel, value.var = vars[[i]], mean)
    
    #get channel names so this can be changed for each
    if(i==1) {
      colNames <- names(currDT)
      colNames <- colNames[!(colNames %in% c("trialType","trialNum"))]
    }
    
    setnames(currDT,colNames,paste0(vars[[i]],"_",colNames))
    setkeyv(currDT,c("trialType","trialNum"))
    
    wideTable <- Reduce(merge,list(wideTable,currDT))
    rm(currDT)
    
  }
  
  wideTable  

}



setwd("/home/dan/Kaggle/UPenn_Seizure/clips")

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)]

#start = proc.time()

#loop through all directories and add frequency ratios to each
for (i in seq_along(dirs)) {
  ratios <- c("dt","da","db","dLg","dHg","ta","tb","tLg","tHg","ab","aLg","aHg","bLg","bHg","LgHg")
  #long original, long train and test files
  print(dirs[[i]])
  
  #Load the earlier derived features in long format
  load(file=paste0(dirs[[i]],"/train+features.rbin"))
  
  #calculate ratio of frequencies by looping along the ratio names
  for(j in seq_along(ratios)) {
    train[,ratios[[j]]:= freqRatio(delta,theta,alpha,beta,lowGamma,highGamma,ratios[[j]]), 
            by = list(trialType,trialNum)]
  }
  
  #load wide train and test files, and combine with ratio features
  load(file=paste0(dirs[[i]],"/train+test.wide.rbin"))
  train.wide <- combineFeaturesTrain(train,train.wide)
  rm(train)
  gc()
  
#----------------------------
  #TEST
  load(file=paste0(dirs[[i]],"/test+features.rbin"))
  #calculate new variables
  for(j in seq_along(ratios)) {
    test[,ratios[[j]]:= freqRatio(delta,theta,alpha,beta,lowGamma,highGamma,ratios[[j]]), 
          by = list(trialType,trialNum)]
  }
  
  test.wide <- combineFeaturesTest(test,test.wide)
  rm(test)
  gc()

  save(list=c("train.wide","test.wide"),file=paste0(dirs[[i]],"/train+test.freqRatios.wide.rbin"))
  rm(list=c("train.wide","test.wide"))
  gc()

}

