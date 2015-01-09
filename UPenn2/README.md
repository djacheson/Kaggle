##UPenn and Mayo Clinic Seizure Detection Challenge
---

This repo contains files for my solution to the [Seizure Detection Challenge](https://www.kaggle.com/c/seizure-detection)


###The competition: Generate two types of models
    1. Ictal / Non-Ictal: Predict whether a clip containes a seizure or not
    2. Early / Not-Early: If a clip does contain a seizure, predict whether it's within 15 seconds of the start of the clip
    
###The data:
    -Direct EEG brain recordings from Dogs and Humans suffering from epilepsy.
    -Thus, a very high-dimensional, multivariate time-series 

###The challenges:
    1. Feature engineering:  Generating an appropriate feature space to differentiate seizure from non-seizure and early vs. late seizures
    2. Severe class imbalance: there are very few clips with seizures

---

###How to use the files in this repo

With raw data included in the _clips/[Subject]_ directories (not included here), the files are intended to be run in the following order:

```
create_R_files.r
preProcess.r
add_freqRatios.r
```

This will create the feature space that is then used by various models located in the _model_code_ directory, with predictions and models saved in _model_save_

__NOTE:__

    1. Raw data isn't included here, but processe / derived data is.
    2. As with all Kaggle competitions, many more models were run than are included in this repo using minor variations of the models located in _model_code_
