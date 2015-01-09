#Functions to process raw ECog data to produce various features

#All functions are designed to be used with a data list created in 'createFiles.r'
#These functions are applied for every clip for every subject resulting in 1 row per clip
# and columns the features for each channel (e.g., theta_1_mean = mean spectral power for theta band at channel 1)

#The data list has various elements such as:
# data$data : chanXsamples matrix of ECoG data
# data$sampleRate : the sampling rate of the timeseries

# Example use of getSpecPower_allwith data.table :  
# dt[,c("highGamma","lowGamma","beta","alpha","theta") := getSpecPower_all(data,400)] 

#Functions with _Half perform the function on the first and second half of a clip
#  separately and calculate the difference.
#  The idea was that there might be differences in the signal the closer to a seizure you are



getSpecPower_all <- function(data) {
  #function to return the mean and standard deviation of spectral power for
  # all major frequency bands using the multitaper method

  #INPUT:
  # data  :  raw data object 

  #OUT: 
  # dataframe with columns representing mean / sd of log10 spectral power for each channel
  
  sampleRate <- data$freq
  library(multitaper)
  #start = proc.time()
  #get full spectral power across all frequency bands using multitaper method
  spec <- apply(data$data,1,function(x) spec.mtm(ts(x,deltat=1/sampleRate),k = 3,plot=FALSE))
  #print(proc.time()-start)[[3]]
  
  
  hiFreq = max(spec[[1]]$freq) 

  #extract various frequency bands from the matrix of spectral power and rename 
  # columns accordingly
  highGamma_mean <- data.frame(lapply(spec,function(x) mean(log10(x$spec[x$freq>100 & x$freq <= hiFreq]))))
  names(highGamma_mean) <- paste0("hg_m_",data$channels)
  
  highGamma_sd <- data.frame(lapply(spec,function(x) sd(log10(x$spec[x$freq>100 & x$freq <= hiFreq]))))
  names(highGamma_sd) <- paste0("hg_sd_",data$channels)
  
  lowGamma_mean <- data.frame(lapply(spec,function(x) mean(log10(x$spec[x$freq>32 & x$freq <= 100]))))
  names(lowGamma_mean) <- paste0("lg_m_",data$channels)
  
  lowGamma_sd <- data.frame(lapply(spec,function(x) sd(log10(x$spec[x$freq>32 & x$freq <= 100]))))
  names(lowGamma_sd) <- paste0("lg_sd_",data$channels)
  
  beta_mean <- data.frame(lapply(spec,function(x) mean(log10(x$spec[x$freq>15 & x$freq <= 32]))))
  names(beta_mean) <- paste0("beta_m_",data$channels)
  
  beta_sd <- data.frame(lapply(spec,function(x) sd(log10(x$spec[x$freq>15 & x$freq <= 32]))))
  names(beta_sd) <- paste0("beta_sd_",data$channels)
  
  alpha_mean <- data.frame(lapply(spec,function(x) mean(log10(x$spec[x$freq>8 & x$freq <= 15]))))
  names(alpha_mean) <- paste0("alpha_m_",data$channels)
  
  alpha_sd <- data.frame(lapply(spec,function(x) sd(log10(x$spec[x$freq>8 & x$freq <= 15]))))
  names(alpha_sd) <- paste0("alpha_sd_",data$channels)
  
  theta_mean <- data.frame(lapply(spec,function(x) mean(log10(x$spec[x$freq>=4 & x$freq <= 8]))))
  names(theta_mean) <- paste0("theta_m_",data$channels)
  
  theta_sd <- data.frame(lapply(spec,function(x) sd(log10(x$spec[x$freq>=4 & x$freq <= 8]))))
  names(theta_sd) <- paste0("theta_sd_",data$channels)
  
  delta_mean <- data.frame(lapply(spec,function(x) mean(log10(x$spec[x$freq < 4]))))
  names(delta_mean) <- paste0("delta_m_",data$channels)
  
  delta_sd <- data.frame(lapply(spec,function(x) sd(log10(x$spec[x$freq < 4]))))
  names(delta_sd) <- paste0("delta_sd_",data$channels)
  

  return(cbind(highGamma_mean,lowGamma_mean,beta_mean,alpha_mean,theta_mean,delta_mean,
               highGamma_sd,lowGamma_sd,beta_sd,alpha_sd,theta_sd,delta_sd))
}

getSpecPower_all_Half <- function(data) {
  #function to get spectral power for the 1st and 2nd half of a clip as well
  # as well as the difference between the two 
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing mean / sd spectral power for each channel
  #  for each half and the difference between halves
  
  #Split the clip in half
  allSamps <- dim(data$data)[2]
  first_half <- data
  first_half$data <- data$data[,1:round(allSamps)/2]
  second_half <- data
  second_half$data <- data$data[,round(allSamps)/2:allSamps]
  
  #get spectral power (mean / sd) for each half
  first <- getSpecPower_all(first_half)
  second <- getSpecPower_all(second_half)
  
  #difference between halves
  dat <- second-first
  
  #rename the columns
  names(dat) <- paste0(names(dat),"_HalfDiff")
  names(first) <- paste0(names(first),"_firstHalf")
  names(second) <- paste0(names(second),"_secondHalf")
  dat <- cbind(dat,first,second)
  dat
}


getCoherence_all <- function(data) {
  #function to return the mean and sd of the spectral coherence between channels 
  #  at different frequency bands using the multitaper method
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing spectral coherence at different frequency
  #  bands between pairs of channels
  
  sampleRate <- data$freq
  library(multitaper)
  #start = proc.time()
  #calculate power spectrum using multitaper
  spec <- apply(data$data,1,function(x) spec.mtm(ts(x,deltat=1/sampleRate),k = 8,nw=20,plot=FALSE, returnInternals = TRUE))
  #print(proc.time()-start)[[3]]
  
  #loop through channels to get the coherence of all combinations of channels
  all_coherence <- list()
  counter = 1
  for(i in data$channels) {
    for(j in data$channels) {
      if(i==j){} else {
        all_coherence[[counter]] <- mtm.coh(spec[[i]],spec[[j]],plot=F)
        all_coherence[[counter]]$channels <- paste0("c",i,"_c",j)
        counter = counter +1
      }
    }
  }

  
  freqs <- all_coherence[[1]]$freq
  
  #extract mean coherence in different frequency bands and rename columns
  delta_mean <- data.frame(lapply(1:length(all_coherence), function(x) mean(log10(all_coherence[[x]]$msc[freqs<=4]))))
  names(delta_mean) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_mean_delta_coh")
  
  theta_mean <- data.frame(lapply(1:length(all_coherence), function(x) mean(log10(all_coherence[[x]]$msc[freqs>4 & freqs<=8]))))
  names(theta_mean) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_mean_theta_coh")
  
  alpha_mean <- data.frame(lapply(1:length(all_coherence), function(x) mean(log10(all_coherence[[x]]$msc[freqs>=8 & freqs<14]))))
  names(alpha_mean) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_mean_alpha_coh")
  
  beta_mean <- data.frame(lapply(1:length(all_coherence), function(x) mean(log10(all_coherence[[x]]$msc[freqs>14 & freqs<=32]))))
  names(beta_mean) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_mean_beta_coh")
  
  lowGamma_mean <- data.frame(lapply(1:length(all_coherence), function(x) mean(log10(all_coherence[[x]]$msc[freqs>32 & freqs<=100]))))
  names(lowGamma_mean) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_mean_lowGamma_coh")
  
  highGamma_mean <- data.frame(lapply(1:length(all_coherence), function(x) mean(log10(all_coherence[[x]]$msc[freqs>100 & freqs<=200]))))
  names(highGamma_mean) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_mean_highGamma_coh")
  
  #SD of coherence in different bands
  delta_sd <- data.frame(lapply(1:length(all_coherence), function(x) sd(log10(all_coherence[[x]]$msc[freqs<=4]))))
  names(delta_sd) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_sd_delta_coh")
  
  theta_sd <- data.frame(lapply(1:length(all_coherence), function(x) sd(log10(all_coherence[[x]]$msc[freqs>4 & freqs<=8]))))
  names(theta_sd) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_sd_theta_coh")
  
  alpha_sd <- data.frame(lapply(1:length(all_coherence), function(x) sd(log10(all_coherence[[x]]$msc[freqs>=8 & freqs<14]))))
  names(alpha_sd) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_sd_alpha_coh")
  
  beta_sd <- data.frame(lapply(1:length(all_coherence), function(x) sd(log10(all_coherence[[x]]$msc[freqs>14 & freqs<=32]))))
  names(beta_sd) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_sd_beta_coh")
  
  lowGamma_sd <- data.frame(lapply(1:length(all_coherence), function(x) sd(log10(all_coherence[[x]]$msc[freqs>32 & freqs<=100]))))
  names(lowGamma_sd) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_sd_lowGamma_coh")
  
  highGamma_sd <- data.frame(lapply(1:length(all_coherence), function(x) sd(log10(all_coherence[[x]]$msc[freqs>100 & freqs<=200]))))
  names(highGamma_sd) <- paste0(sapply(1:length(all_coherence), function(x) all_coherence[[x]]$channels),"_sd_highGamma_coh")
  
  all <- cbind(delta_mean,theta_mean,alpha_mean,beta_mean,lowGamma_mean,highGamma_mean,delta_sd,theta_sd,alpha_sd,beta_sd,lowGamma_sd,highGamma_sd)
  
  return(all)
}

getCoherence_all_Half <- function(data) {
  #Function to get spectral coherence at different bands for each half of clip,
  #  as well as the difference between halves
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing spectral coherence at different frequency
  #  bands between pairs of channels for each half, and the difference between halves
  
  
  allSamps <- dim(data$data)[2]
  first_half <- data
  first_half$data <- data$data[,1:round(allSamps)/2]
  second_half <- data
  second_half$data <- data$data[,round(allSamps)/2:allSamps]
  
  first <- getCoherence_all(first_half)
  second <- getCoherence_all(second_half)
  
  dat <- second-first
  names(dat) <- paste0(names(dat),"_HalfDiff")
  names(first) <- paste0(names(first),"_firstHalf")
  names(second) <- paste0(names(second),"_secondHalf")
  dat <- cbind(dat,first,second)
  dat
}



getCors <- function(data, scale = "none") {
  #Function to calculate the correlations amongst channels and eigenvalue
  #  decomposition of the correlation matrix 
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing correlation in signal between channels
  # and the eigenvalues of the correlation matrix 
    
    #Normalize the data either across a channel or across time (samples)
    if(scale == "samples") {
        data$data <- scale(data$data)
    } else if(scale == "channels") {
        data$data = t(scale(t(data$data)))
    }
    
    #correlation matrix
    allCors <- cor(t(data$data))
    
    #eigenvalue decomposition
    allCors_eigen <- data.frame(t(eigen(allCors)[[1]]))
    
    #rename according to the type of scaling done
    if(scale == "samples") {
        names(allCors_eigen) <- paste0("Cor_eigen_",1:ncol(allCors_eigen),"_sampleScale")
    } else if(scale == "channels") {
        names(allCors_eigen) <- paste0("Cor_eigen_",1:ncol(allCors_eigen),"_channelScale")
    } else {
        names(allCors_eigen) <- paste0("Cor_eigen_",1:ncol(allCors_eigen))
    }
    
    #here, get the names of combinations of channels
    temp <- upper.tri(allCors)
    rownames(temp) <- paste0("c",1:nrow(temp))
    colnames(temp) <- paste0("c",1:ncol(temp))
    
    #provide channel names to the columns
    allNames <- c()
    for(i in 1:(nrow(temp)-1)) {
        if(scale == "samples") {
            allNames <- c(allNames,paste0("Cor_",rownames(temp)[[i]],"_", names(temp[i,][temp[i,]]),"_sampleScale"))
        } else if(scale == "channels") {
            allNames <- c(allNames,paste0("Cor_",rownames(temp)[[i]],"_", names(temp[i,][temp[i,]]),"_channelScale"))
        } else {
            allNames <- c(allNames,paste0("Cor_",rownames(temp)[[i]],"_", names(temp[i,][temp[i,]])))
        }
        
    }
    
    #only get upper triangle of correlation matrix, then assign appropriate names
    allCors <- data.frame(t(allCors[upper.tri(allCors)]))
    
    names(allCors) <- allNames
    allCors <- cbind(allCors, allCors_eigen)
    return(allCors)
}

getCors_Half <- function(data) {
  #Function to calculate the correlations amongst channels and eigenvalue
  #  decomposition of the correlation matrix for each half of the data as 
  #  well as the difference between halves
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing spectral coherence at different frequency
  #  bands between pairs of channels for each half, and the difference between halves
  
  allSamps <- dim(data$data)[2]
  first_half <- data
  first_half$data <- data$data[,1:round(allSamps)/2]
  second_half <- data
  second_half$data <- data$data[,round(allSamps)/2:allSamps]
  
  first <- getCors(first_half)
  second <- getCors(second_half)
  
  dat <- second-first
  names(dat) <- paste0(names(dat),"_HalfDiff")
  names(first) <- paste0(names(first),"_firstHalf")
  names(second) <- paste0(names(second),"_secondHalf")
  dat <- cbind(dat,first,second)
  dat
  
}


getMI_all <- function(data, scale = "none") {
  # get mutual information between channels as well as eigenvalues of MI matrix
  # scale options: none, samples (i.e., across channels), channels (i.e., across samples)
  
  #Basically the same as getCor_all, but using MI instead of correlation
  
  #!!!!!!!!!
  # WARNING: This function takes a very long time to run
  #!!!!!!!!!
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing mutual information between channels as well
  #  as the eigenvalues of the MI matrix
  
    
    library(infotheo)
    
    #Normalize first?
    if(scale == "samples") {
        data$data <- scale(data$data)
    } else if(scale == "channels") {
        data$data = t(scale(t(data$data)))
    }
    
    #Mutual information
    allMI <- mutinformation(discretize(t(data$data),nbins = ncol(data$data)^0.6)) #defalut is 1/3; 0.6 based on Reshef MIC article...see correlation comparsion article
    #eigen values of MI matrix
    allMI_eigen <- data.frame(t(eigen(allMI)[[1]]))
    if(scale == "samples") {
        names(allMI_eigen) <- paste0("MI_eigen_",1:ncol(allMI_eigen),"_sampleScale")
    } else if(scale == "channels"){
        names(allMI_eigen) <- paste0("MI_eigen_",1:ncol(allMI_eigen),"_channelScale")
    } else {
        names(allMI_eigen) <- paste0("MI_eigen_",1:ncol(allMI_eigen))
    }
    
    #here, get the names of combinations of channels
    temp <- upper.tri(allMI)
    rownames(temp) <- paste0("c",1:nrow(temp))
    colnames(temp) <- paste0("c",1:ncol(temp))
    
    allNames <- c()
    for(i in 1:(nrow(temp)-1)) {
        if(scale == "samples") {
            allNames <- c(allNames,paste0("MI_",rownames(temp)[[i]],"_", names(temp[i,][temp[i,]]),"_sampleScale"))
        } else if(scale == "channels") {
            allNames <- c(allNames,paste0("MI_",rownames(temp)[[i]],"_", names(temp[i,][temp[i,]]),"_channelScale"))
        } else {
            allNames <- c(allNames,paste0("MI_",rownames(temp)[[i]],"_", names(temp[i,][temp[i,]])))
        }
    }
    
    
    #only get upper triangle of correlation matrix, then assign appropriate names
    allMI <- data.frame(t(allMI[upper.tri(allMI)]))
    names(allMI) <- allNames
    allMI <- cbind(allMI, allMI_eigen)
    allMI
    
}

  
getSD <- function(data) {
  #function to calculate the standard deviation of each channel
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing sd of each channel
  
  dataSD <- data.frame(t(apply(data$data,1,function(x) sd(x))))
  names(dataSD) <- paste0("sd_",1:ncol(dataSD))
  dataSD
}

getSD_Half <- function(data) {
  #function to calculate the standard deviation of each channel separately for
  #  the first and second half as well as the difference
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing sd of each channel
    allSamps <- dim(data$data)[2]
    first_half <- data
    first_half$data <- data$data[,1:round(allSamps)/2]
    second_half <- data
    second_half$data <- data$data[,round(allSamps)/2:allSamps]
    
    first <- getSD(first_half)
    second <- getSD(second_half)
    
    dat <- second-first
    names(dat) <- paste0(names(dat),"_HalfDiff")
    names(first) <- paste0(names(first),"_firstHalf")
    names(second) <- paste0(names(second),"_secondHalf")
    dat <- cbind(dat,first,second)
    dat
  }


get_autoCor <- function(data,lag) {
  #function to calculate the autocorrelation with a given lag
  #  autcorrelation is just correlation of a signal with a lagged version of itself
  
  #Designed to be called by get_autoCorVar, get_autoCorMean below
  
  #NOTE: This wasn't very predictive of the difference between ictal / interictal
  
  #INPUT:
  # data  :  raw data object
  # lag   :  how many samples back to calculate autoCor on
  
  #OUT: 
  # correlation matrix of signal correlated with lagged version of itself
  
  
  autoCor <- cor(data[(lag+1):length(data)],data[1:(length(data)-lag)])
  autoCor
}

get_autoCorVar <- function(data,maxLag) {
  #function to calculate the variance in canonical correlation of a signal over a range of lags
  
  #INPUT:
  # data  :  raw data object
  # maxLag   :  maximum lag in samples to go backwards
  
  #OUT: 
  # dataframe with columns corresponding to the variance in autocorrelation for each channel
  
  #Function to apply get the variance of the autocorrelation for a given row
  autoCorVarTrial <- function(row) {
    var(sapply(seq(1,maxLag,2), function(x) get_autoCor(row,x)))
  }

  data_autoCorVar <- data.frame(t(apply(data$data,1,function(x) autoCorVarTrial(x))))
  names(data_autoCorVar) <- paste0("autoCor_var_",1:ncol(data_autoCorVar))
  return(data_autoCorVar)
}
  
get_autoCorVar_Half <- function(data,maxLag) {
  #Apply autoCor var to each half separately
    allSamps <- dim(data$data)[2]
    first_half <- data
    first_half$data <- data$data[,1:round(allSamps)/2]
    second_half <- data
    second_half$data <- data$data[,round(allSamps)/2:allSamps]
    
    first <- get_autoCorVar(first_half, maxLag)
    second <- get_autoCorVar(second_half,maxLag)
    
    dat <- second-first
    names(dat) <- paste0(names(dat),"_HalfDiff")
    names(first) <- paste0(names(first),"_firstHalf")
    names(second) <- paste0(names(second),"_secondHalf")
    dat <- cbind(dat,first,second)
    dat
}
  
  

get_autoCorMean <- function(data,maxLag) {
  #function to calculate the mean canonical correlation of a signal over a range of lags
  
  #INPUT:
  # data  :  raw data object
  # maxLag   :  maximum lag in samples to go backwards
  
  #OUT: 
  # dataframe with columns corresponding to the variance in autocorrelation for each channel
  
  #Function to apply get the mean of the autocorrelation for a given row
  autoCorMeanTrial <- function(row) {
    mean(sapply(seq(1,maxLag,2), function(x) get_autoCor(row,x)))
    
  }
  
  data_autoCorMean <- data.frame(t(apply(data$data,1,function(x) autoCorMeanTrial(x))))
  names(data_autoCorMean) <- paste0("autoCor_mean_",1:ncol(data_autoCorMean))
  data_autoCorMean
}
  
get_autoCorMean_Half <- function(data,maxLag) {
    #Function to apply autoCorMean to each half
    allSamps <- dim(data$data)[2]
    first_half <- data
    first_half$data <- data$data[,1:round(allSamps)/2]
    second_half <- data
    second_half$data <- data$data[,round(allSamps)/2:allSamps]
    
    first <- get_autoCorMean(first_half, maxLag)
    second <- get_autoCorMean(second_half,maxLag)
    
    dat <- second-first
    names(dat) <- paste0(names(dat),"_HalfDiff")
    names(first) <- paste0(names(first),"_firstHalf")
    names(second) <- paste0(names(second),"_secondHalf")
    dat <- cbind(dat,first,second)
    dat
  }
  

get_skew <- function(data) {
  #function to calculate the skew of data in each channel
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing skew of data in each channel
  
  library(e1071)
  data_skew <- data.frame(t(apply(data$data,1,function(x) skewness(x))))
  names(data_skew) <- paste0("skew_",1:ncol(data_skew))
  return(data_skew)
  
}
  
get_skew_Half <- function(data) {
  #Function to apply get_skew to each half separately and calculate the difference
  allSamps <- dim(data$data)[2]
  first_half <- data
  first_half$data <- data$data[,1:round(allSamps)/2]
  second_half <- data
  second_half$data <- data$data[,round(allSamps)/2:allSamps]
  
  first <- get_skew(first_half)
  second <- get_skew(second_half)
  
  dat <- second-first
  names(dat) <- paste0(names(dat),"_HalfDiff")
  names(first) <- paste0(names(first),"_firstHalf")
  names(second) <- paste0(names(second),"_secondHalf")
  dat <- cbind(dat,first,second)
  dat
}


  
  
get_kurtosis<- function(data) {
  #function to calculate the kurtosis of data in each channel
  
  #NOTE: This didn't produce useful features in this dataset
  
  #INPUT:
  # data  :  raw data object 
  
  #OUT: 
  # dataframe with columns representing kurtosis of data in each channel
  
  
  
  library(e1071)
  data_kurt <- data.frame(t(apply(data$data,1,function(x) kurtosis(x))))
  names(data_kurt) <- paste0("kurt_",1:ncol(data_kurt))
  data_kurt
  
}

  
get_kurtosis_Half <- function(data) {
  #fucntion to apply get_kurtosis to each half and calculate the difference
  
    allSamps <- dim(data$data)[2]
    first_half <- data
    first_half$data <- data$data[,1:round(allSamps)/2]
    second_half <- data
    second_half$data <- data$data[,round(allSamps)/2:allSamps]
    
    first <- get_kurtosis(first_half)
    second <- get_kurtosis(second_half)
    
    dat <- second-first
    names(dat) <- paste0(names(dat),"_HalfDiff")
    names(first) <- paste0(names(first),"_firstHalf")
    names(second) <- paste0(names(second),"_secondHalf")
    dat <- cbind(dat,first,second)
    dat
}
  
  
applyButterworth <- function(data, band = c(0,4)) {
    #Function to return a bandbass butterworth filtered version of the data
    # Used in various functions below to calculate various features for specific frequency bands

    #INPUT: 
    # data: ecog data object, 
    # band: frequency band of interest
    
    #OUTPUT:
    #  ecog data object with filtered data
    
    library(signal)
    maxFreq = data$freq
    Nyquist = maxFreq / 2
    #express band as a proportion of the Nyquist frequency (0-1)
    band = band / Nyquist
    
    #third order butterworth filter
    butter_filt <-  butter(3,W = band)
    #applied to each row (channel) of the data
    data$data <- t(apply(data$data,1,function(x) filter(butter_filt, x)))
    data
}

#---------------------------
#BELOW- using butterworth filter to get various stas for specific frequency bands

getCors_freq <- function(data, freqBand = "alpha") {
    #Function to get correlations and eigenvalues for pre-specified frequency bands using
    # First performsn bandpass filter on the data, 
    # then gets correlations and eigenvalues using the getCors function
    
    #INPUT:
    #  data  :  raw data object
    #  freqBand  : character string representing a band of interest
  
    #OUTPUT:
    #  dataframe with correlations of bandpassed signal between channels
  
    freqs <- c("delta","theta","alpha","beta","lg","hg")
    if(!(freqBand %in% freqs)) {
        print("Error: Need to specify frequency band:")
        print(paste(freqs,sep = ","))
    } else {
        if(freqBand == "delta") band = c(0,4) 
        if(freqBand == "theta") band = c(4,8)
        if(freqBand == "alpha") band = c(8,14)    
        if(freqBand == "beta") band = c(14,32)
        if(freqBand == "lg") band = c(32,100) #low gamma
        if(freqBand == "hg") band = c(100,200) #high gamma
        
        
        data <- applyButterworth(data,band)
        corrs <- getCors(data)
        names(corrs) <- paste0(freqBand, "_", names(corrs))
        corrs
    }
}

getCors_all_freq <- function(data) {
    #Function to apply getCors_freq to all the different frequency bands
    
    #!!!!!
    #WARNING: This creates a pretty large feature space depending on the number channels
    #!!!!!!
  
    #INPUT:
    #  data  :  raw data object
  
    #OUTPUT:
    #  dataframe with columns representing the correlation across channels and eigenvalues for 
    #    each of the major frequency bands
  

    delta <- getCors_freq(data,"delta")
    theta <- getCors_freq(data,"theta")
    alpha <- getCors_freq(data,"alpha")
    beta <- getCors_freq(data,"beta")
    lg <- getCors_freq(data,"lg")
    hg <- getCors_freq(data,"hg")
    
    data.frame(c(delta,theta,alpha,beta,lg,hg))
}



getCors_all_freq_Half <- function(data) {
    #Function to apply getCors_all_freq to each half separately
    allSamps <- dim(data$data)[2]
    first_half <- data
    first_half$data <- data$data[,1:round(allSamps)/2]
    second_half <- data
    second_half$data <- data$data[,round(allSamps)/2:allSamps]
    
    first <- getCors_all_freq(first_half)
    second <- getCors_all_freq(second_half)
    
    dat <- second-first
    names(dat) <- paste0(names(dat),"_HalfDiff")
    names(first) <- paste0(names(first),"_firstHalf")
    names(second) <- paste0(names(second),"_secondHalf")
    dat <- cbind(dat,first,second)
    dat
    
}

getDerivs_freq <- function(data, freqBand = "alpha") {
  #Function to get the mean, standard deviation, skew and kurtosis of the first
  #  derivative of the data for a given frequency band
  
  #!!!!!
  #WARNING: This creates a pretty large feature space depending on the number channels
  #!!!!!!
  
  #INPUT:
  #  data  :  raw data object
  #  freqBand  : character string representing a band of interest
  
  #OUTPUT:
  #  dataframe with columns representing the  mean, standard deviation, skew and 
  #   kurtosis of the first derivative for each channel
  
  
    library(prospectr)
    library(e1071)
    #get correlations for pre-specified frequency bands
    #first, perform bandpass butterworth filter, then get correlations and eigenvalues using getCors function
    freqs <- c("delta","theta","alpha","beta","lg1","lg2","hg","all")
    if(!(freqBand %in% freqs)) {
        print("Error: Need to specify frequency band:")
        print(paste(freqs,sep = ","))
    } else {
        if(freqBand == "delta") band = c(0,4) 
        if(freqBand == "theta") band = c(4,8)
        if(freqBand == "alpha") band = c(8,14)    
        if(freqBand == "beta") band = c(14,32)
        if(freqBand == "lg1") band = c(32,55)
        if(freqBand == "lg2") band = c(55,100)
        if(freqBand == "hg") band = c(100,200)
        if(freqBand == "all") band = c(0,200)        
        
        #filter then get first derivative using Gap-Segment derivative
        data <- applyButterworth(data,band)
        derivs <- gapDer(data$data)
        #names(derivs) <- paste0(freqBand, "_", names(derivs))
        
        #mean of derivatives (note rectifed by sqrt(deriv^2))
        mean_derivs <- data.frame(t(apply(derivs, 1, function(x) mean(sqrt(x^2)))))
        names(mean_derivs) <- paste0(freqBand,"_1D_mean_",1:ncol(mean_derivs))
        
        #sd of derivatives
        sd_derivs <- data.frame(t(apply(derivs, 1, function(x) sd(x))))
        names(sd_derivs) <- paste0(freqBand,"_1D_sd_",1:ncol(sd_derivs))
        
        #skew of derivatives
        skew_derivs <- data.frame(t(apply(derivs, 1, function(x) skewness(x))))
        names(skew_derivs) <- paste0(freqBand,"_1D_skew_",1:ncol(skew_derivs))
        
        #kurtosis of derivatives
        kurt_derivs <- data.frame(t(apply(derivs, 1, function(x) kurtosis(x))))
        names(kurt_derivs) <- paste0(freqBand,"_1D_kurt_",1:ncol(skew_derivs))
        
        cbind(mean_derivs, sd_derivs, skew_derivs, kurt_derivs)
    }
}

getDerivs_freq2 <- function(data, freqBand = c(0,10)) {
  #Function to get the mean, standard deviation, skew and kurtosis as well
  # as correlations across channels for the given frequency band 
  # for the raw, first and second derivatives of the signal 
  
  
  #!!!!!
  #WARNING: This creates a very large feature space depending on the number channels
  #!!!!!!
  
  #INPUT:
  #  data  :  raw data object
  #  freqBand  : character string representing a band of interest
  
  #OUTPUT:
  #  dataframe with columns representing the  mean, standard deviation, skew and 
  #   kurtosis of the first derivative for each channel
  
    #This does more than just get derivatives for a given frequency band
    #it calculates summary statistics for raw, first and second derivatives
    #as well as correlations across channels for the given frequency band and raw / derivative signals
    library(prospectr)
    library(e1071)
    #first, perform bandpass butterworth filter, then get correlations and eigenvalues using getCors function
    band = freqBand
    bandName = paste0("b",band[[1]],"_",band[[2]])
    
    #filter using bandbass butterworth filter
    data <- applyButterworth(data,band)
    
    #raw data first
    raw <- data$data
        
    #mean (note rectifed by sqrt(deriv^2))
    mean_raw <- data.frame(t(apply(raw, 1, function(x) mean(sqrt(x^2)))))
    names(mean_raw) <- paste0(bandName,"_raw_mean_",1:ncol(mean_raw))
    
    #sd 
    sd_raw <- data.frame(t(apply(raw, 1, function(x) sd(x))))
    names(sd_raw) <- paste0(bandName,"_raw_sd_",1:ncol(sd_raw))
    
    #skew 
    skew_raw <- data.frame(t(apply(raw, 1, function(x) skewness(x))))
    names(skew_raw) <- paste0(bandName,"_raw_skew_",1:ncol(skew_raw))
    
    #kurtosis 
    kurt_raw <- data.frame(t(apply(raw, 1, function(x) kurtosis(x))))
    names(kurt_raw) <- paste0(bandName,"_raw_kurt_",1:ncol(skew_raw))
    
    #correlations
    cor_raw <- getCors(data)
    names(cor_raw) <- paste0(bandName,"_",names(cor_raw),"_raw")
    
    #--------------------------------------
    #first derivs - using gapDerivative with a width of 3
    derivs <- gapDer(data$data, m = 1, w = 3)
    data2 <- data
    data2$data <- derivs
    #names(derivs) <- paste0(freqBand, "_", names(derivs))
    
    #mean of derivatives (note rectifed by sqrt(deriv^2))
    mean_derivs <- data.frame(t(apply(derivs, 1, function(x) mean(sqrt(x^2)))))
    names(mean_derivs) <- paste0(bandName,"_1D_mean_",1:ncol(mean_derivs))
    
    #sd of derivatives
    sd_derivs <- data.frame(t(apply(derivs, 1, function(x) sd(x))))
    names(sd_derivs) <- paste0(bandName,"_1D_sd_",1:ncol(sd_derivs))
    
    #skew of derivatives
    skew_derivs <- data.frame(t(apply(derivs, 1, function(x) skewness(x))))
    names(skew_derivs) <- paste0(bandName,"_1D_skew_",1:ncol(skew_derivs))
    
    #kurtosis of derivatives
    kurt_derivs <- data.frame(t(apply(derivs, 1, function(x) kurtosis(x))))
    names(kurt_derivs) <- paste0(bandName,"_1D_kurt_",1:ncol(skew_derivs))
    
    #correlations
    data2 <- data
    data2$data <- derivs
    cor_deriv <- getCors(data2)
    names(cor_deriv) <- paste0(bandName,"_",names(cor_deriv),"_1D")
    
    #---------------------------------------
    #second derivs
    derivs2 <- gapDer(data$data,m = 2, w= 1)
    #names(derivs) <- paste0(freqBand, "_", names(derivs))
    
    #mean of derivatives (note rectifed by sqrt(deriv^2))
    mean_derivs2 <- data.frame(t(apply(derivs2, 1, function(x) mean(sqrt(x^2)))))
    names(mean_derivs2) <- paste0(bandName,"_2D_mean_",1:ncol(mean_derivs2))
    
    #sd of derivatives
    sd_derivs2 <- data.frame(t(apply(derivs2, 1, function(x) sd(x))))
    names(sd_derivs2) <- paste0(bandName,"_2D_sd_",1:ncol(sd_derivs2))
    
    #skew of derivatives
    skew_derivs2 <- data.frame(t(apply(derivs2, 1, function(x) skewness(x))))
    names(skew_derivs2) <- paste0(bandName,"_2D_skew_",1:ncol(skew_derivs2))
    
    #kurtosis of derivatives
    kurt_derivs2 <- data.frame(t(apply(derivs2, 1, function(x) kurtosis(x))))
    names(kurt_derivs2) <- paste0(bandName,"_2D_kurt_",1:ncol(skew_derivs2))
        
    #correlations
    data2 <- data
    data2$data <- derivs2
    cor_deriv2 <- getCors(data2)
    names(cor_deriv2) <- paste0(bandName,"_",names(cor_deriv2),"_2D")
    
    cbind(mean_raw, sd_raw, skew_raw, kurt_raw,cor_raw,
          mean_derivs, sd_derivs, skew_derivs, kurt_derivs,cor_deriv,
          mean_derivs2, sd_derivs2, skew_derivs2, kurt_derivs2,cor_deriv2
          )
}


getDerivs_all_freq <- function(data) {
   #Function to apply getDerivs_freq to all the major frequency bands
  
    delta <- getDerivs_freq(data,"delta")
    theta <- getDerivs_freq(data,"theta")
    alpha <- getDerivs_freq(data,"alpha")
    beta <- getDerivs_freq(data,"beta")
    lg1 <- getDerivs_freq(data,"lg1")
    lg2 <- getDerivs_freq(data,"lg2")
    hg <- getDerivs_freq(data,"hg")
    all <- getDerivs_freq(data,"all")
    
    
    data.frame(c(delta,theta,alpha,beta,lg1,lg2,hg,all))
}

getDerivs_all_freq2 <- function(data) {
  #Function to apply getDerivs_freq2 to various frequency bans, here, 10 hz bands from 0.1-150 hz
    
    
    freqs <- c(0.1,seq(10,150,by = 10))
    data_list <- list()
    for(i in 1:(length(freqs)-1)) {
        data_list[[i]] <- getDerivs_freq2(data, freqBand = c(freqs[i],freqs[i+1]))
    }
    
    
    data.frame(Reduce(cbind, data_list))
}

#!!!!!!!!!!!!!!!
#WARNING: THE FUNCTIONS BELOW TAKE AN EXCESSIVELY LONG TIME TO RUN ON A SINGLE MACHINE
#!!!!!!!!!!!!!!!

getMI_freq <- function(data, freqBand = "alpha") {
  #Function to get mutual information for pre-specified frequency bands 
  # First performs bandpass filter on the data, 
  
  #INPUT:
  #  data  :  raw data object
  #  freqBand  : character string representing a band of interest
  
  #OUTPUT:
  #  dataframe with correlations of bandpassed signal between channels
  
    freqs <- c("delta","theta","alpha","beta","lg","hg")
    if(!(freqBand %in% freqs)) {
        print("Error: Need to specify frequency band:")
        print(paste(freqs,sep = ","))
    } else {
        if(freqBand == "delta") band = c(0,4) 
        if(freqBand == "theta") band = c(4,8)
        if(freqBand == "alpha") band = c(8,14)    
        if(freqBand == "beta") band = c(14,32)
        if(freqBand == "lg") band = c(32,100)
        if(freqBand == "hg") band = c(100,200)
        
        
        data <- applyButterworth(data,band)
        MI <- getMI_all(data)
        names(MI) <- paste0(freqBand, "_", names(MI))
        MI
    }
}

getMI_all_freq <- function(data) {
    #Function to apply getMI_freq to all the major frequency bands
    #get MI for all frequency bands using getMI_freq 
    delta <- getMI_freq(data,"delta")
    theta <- getMI_freq(data,"theta")
    alpha <- getMI_freq(data,"alpha")
    beta <- getMI_freq(data,"beta")
    lg <- getMI_freq(data,"lg")
    hg <- getMI_freq(data,"hg")
    
    data.frame(c(delta,theta,alpha,beta,lg,hg))
}

getMI_all_freq_Half <- function(data) {
    #Function to apply getMI_all_freq to each half separately
    allSamps <- dim(data$data)[2]
    first_half <- data
    first_half$data <- data$data[,1:round(allSamps)/2]
    second_half <- data
    second_half$data <- data$data[,round(allSamps)/2:allSamps]
    
    first <- getMI_all_freq(first_half)
    second <- getMI_all_freq(second_half)
    
    dat <- second-first
    names(dat) <- paste0(names(dat),"_HalfDiff")
    names(first) <- paste0(names(first),"_firstHalf")
    names(second) <- paste0(names(second),"_secondHalf")
    dat <- cbind(dat,first,second)
    dat
    
}


