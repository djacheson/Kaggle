### Code for generating different models and ultimately combining the predictions

Two models are generated for each subject:
1. Ictal / Non-Ictal - whether a segment contains a seizure or not
2. Latency <= 15 seconds - whether the seizure is within 15 seconds of the start of the clip

The general flow for each is as follows:
For each subject:
1. Load tidy data (1 row per clip, columns are features)
2. Set tuning parameters for different models
3. Upsample seizure / early segments to deal with severe class imbalance
4. Run models with parameters tuned using cross-validation
5. Generate probabilistic predictions

After this is run for each subject and each , predictions are combined across subjects

---
__NOTE__: In addition to the general model code, there is code for combining (i.e., stacking) model predictions located in _combined_model_predictions.r_
