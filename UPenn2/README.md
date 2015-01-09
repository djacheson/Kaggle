##American Epilepsy Society Seizure Prediction Challenge
---

This repo contains files for my solution to the [Second Seizure Detection Challenge](http://www.kaggle.com/c/seizure-prediction)


###The competition: Generate a model to predict
    1. Ictal / Interictal : is a clip 10 minutes before a seizure or not?
    
###The data:
    -Direct EEG brain recordings from Dogs and Humans suffering from epilepsy.
    -Thus, a very high-dimensional, multivariate time-series 

###The challenges:
    1. Feature engineering:  Generating an appropriate feature space to differentiate seizure from non-seizure and early vs. late seizures
    2. Severe class imbalance: there are very few clips preceding seizures
    3. Data size - relative to the first [Seizure Detection Challenge](https://www.kaggle.com/c/seizure-detection)
        Here, the recorded data is 10 minutes vs 1 minute. 

---

###How to use the files in this repo

With raw data included in the Data/[Subject]_ directories (not included here), the files are intended to be run in the following order:

```
#To read in data and create an initial features space
PreProcessing/createFiles.r
#To create additional features
createFilesBatch.R (which calls on createFiles_oneDir.X.r files)

```

This will create the feature space that is then used by various models located in the _ModelCode_ directory, with predictions in _ModelPredictions_ and saved models in _ModelSave_

__NOTE:__

    1. Raw data isn't included here, nor is processed data due to size limitations
    2. As with all Kaggle competitions, many more models were run than are included in this repo using minor variations of the models located in _ModelCode_
