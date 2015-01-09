#Functions to generate features from raw data.
#Features include:
#  spectral power at various frequency bands
#  range - abs(min - max)
#  extreme values - abs(5th - 95th percentile)
#  ratio of spectral powers to each other (see add_freqRatios.r)

#INCLUDED BUT NOT USED IN FINAL MODELS:
#  Average and Variance of Canonical Correlation (i.e., autocorrelation of lagged signal)
#  Global spectral power (used to normalize spectral power per trial)

getSpecPower <- function(data,sampleRate=400,loFreq=100,hiFreq='max') {
  #function to return the mean spectral power within a given frequency range for a given timeseries
  #NOTE: later.... if want to do this for a whole matrix (e.g., chanXtime), there are more efficient ways, e.g., via biwavelet package
  
  #INPUT:
  # data  : a timeseries of voltage changes (if it's just numeric vector, it will be coerced below)
  # sampleRate  : the sampling rate of the timeseries
  # loFreq :  low frequency cutoff
  # hiFreq : high frequency cutoff, default is to use maximum frequency 
  
  #OUTPUT:
  # mean spectral power for given timeseries
  
  require(multitaper)
  #print(str(ts(data,deltat=1/sampleRate)))
  spec <- spec.mtm(ts(data,deltat=1/sampleRate),plot=FALSE)
  #print(str(spec))
  if(hiFreq=='max') {
    hiFreq = max(spec$freq) 
  }
  return(mean(spec$spec[spec$freq>=loFreq & spec$freq <= hiFreq]))
}

getSpecPower_all <- function(data,sampleRate=400) {
  #function to return the mean spectral power of all major frequency bands
  # intended for use with data.table to create multiple frequency bands at once
  # e.g., dt[,c("highGamma","lowGamma","beta","alpha","theta") := getSpecPower_all(data,400)] 
  
  #INPUT:
  # data  : a timeseries of voltage changes (if it's just numeric vector, it will be coerced below)
  # sampleRate  : the sampling rate of the timeseries

  #OUTPUT
  # a list containing mean power for each frequency band
  
  require(multitaper)
  #print(str(ts(data,deltat=1/sampleRate)))
  spec <- spec.mtm(ts(data,deltat=1/sampleRate),plot=FALSE)
  #print(str(spec))
  
  hiFreq = max(spec$freq) 

  highGamma <- mean(spec$spec[spec$freq>=100 & spec$freq <= hiFreq])
  lowGamma <- mean(spec$spec[spec$freq>=32 & spec$freq <= 99.99])
  beta  <- mean(spec$spec[spec$freq>=15 & spec$freq <= 31.99])
  alpha <-  mean(spec$spec[spec$freq>=8 & spec$freq <= 14.99])
  theta <- mean(spec$spec[spec$freq>=4 & spec$freq <= 7])
  delta <- mean(spec$spec[spec$freq < 4])

  return(list(highGamma,lowGamma,beta,alpha,theta,delta))
}


getGlobalPower <- function(data,sampleRate=400) {
  #function to return the mean spectral power across all frequencies 
  #INPUT:
  # data  : a timeseries of voltage changes (if it's just numeric vector, it will be coerced below)
  # sampleRate  : the sampling rate of the timeseries
  
  #OUTPUT
  # mean global power
  
  #NOTE: This didn't produce useful features beyond variance estimates
  
  require(multitaper)
  #print(str(ts(data,deltat=1/sampleRate)))
  spec <- spec.mtm(ts(data,deltat=1/sampleRate),plot=FALSE)
  #print(str(spec))
  
  hiFreq = max(spec$freq) 
  
  globPower <- mean(spec$spec)
  return(list(globPower))
}

minMax <- function(data) {
  #function to return the absolute value of the difference between minimum and maximum values of a timeseries
  #INPUT:
  # data  : vector of voltage changes
  
  #OUTPUT
  # abs(min-max)
  
  return(abs(min(data)-max(data)))
}

extremeVals <- function(data,probs=c(0.05,0.95)) {
  #function to return mean voltage difference of 5th and 95th percentile
  # Note: data is rectified and overall mean is removed
  
  #INPUT:
  # data: vector of voltage changes
  
  #OUTPUT:
  # abs(5th - 95th percentile) 
  
  #rectify and remove the mean.. this will hopefully deal with drift issues in the signal
  data <- abs(data) - mean(data)

  ends <- quantile(data,probs)
  return(as.numeric(abs(ends[1]-ends[2])))
}

autoCor <- function(data,lag) {
  #function to calculate the autoorrelation of a signal 
  #  (i.e., correlation of signal with a lagged version of itself)
  
  # INPUT:
  #   data : a vector of voltage changes
  #   lag : number of data points to correlate signal with
  
  #OUTPUT:
  #  correlation of signal with lagged version of itself

  
  autoCor <- cor(data[(lag+1):length(data)],data[1:(length(data)-lag)])
  autoCor
}

autoCorVar <- function(data,maxLag) {
  #function to calculate the variance in autocorrelation of a signal across a range of lags
  #INPUT:
  #  data :  a vector of voltage changes
  #  maxLag : the maximum lag (in sampling points) over which the canonical correlation is calculated
  
  #OUTPUT:
  #  variance of canonical correlations
  
  #NOTE: This didn't produce a useful set of features for this data
  
  #Create a list of canonical correlations starting at lag 1 to lag maxLag by steps of 2
  autoCorList <- list()
  autoCorList <- lapply(seq(1,maxLag,2), function(x) autoCor(data,x))
  autoCorVar <- var(unlist(autoCorList))
  autoCorVar
}

autoCorMean <- function(data,maxLag) {
  #function to calculate the mean autocorrelation of a signal across a range of lags
  #INPUT:
  #  data :  a vector of voltage changes
  #  maxLag : the maximum lag (in sampling points) over which the canonical correlation is calculated
  
  #OUTPUT:
  #  mean of canonical correlations
  
  #NOTE: This didn't produce a useful set of features for this data
  
  autoCorList <- list()
  autoCorList <- lapply(seq(1,maxLag,2), function(x) autoCor(data,x))
  autoCorMean<- mean(unlist(autoCorList))
  autoCorMean
}



#Helper functions to flatten list of matrices into a data frame organized as trial X flattened_Matrix
toDataMat <- function(dataList) {
  #convert a list of matrices into a big data matrix by flattening each matrix
  require(plyr)
  dataMat <- laply(dataList, function(x) as.vector(t(x)))
}

toDataList <- function(dataMat, ncol) {
  #convert a data matrix (trial X flattened_Matrix, into a list of matrices (chan X time)
  require(plyr)
  dataList <- alply(dataMat, 1, function(x) matrix(x,ncol=ncol, byrow=TRUE))
}