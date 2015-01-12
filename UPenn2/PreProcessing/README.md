##Files for preprocessing the data and generating cross-validation folds
---

###Two sets of files are included here:
1. creatFiles[X].r
    Reads in raw data, and generates various features by processing data.
    E.G.:
    - Standard deviation, skewness and kurtosis of data
    - Mean and standard deviation of spectral power
    - Correlations and mutual information between signal in channels
2. extractSequences.R - cross-validation files
    -  generates cross-validation folds based on how the data was collected 
        (i.e., sequences closer in time are in the same fold)
