#Code to combine (i.e., stack) model predictions

#NOTE:
#  1. many more models were generated than were used in the final stack
#  2. The weighting of models in the stack is decidedly a hack, reflecting my 
#     confidence in the models based on their leaderboard performance

library(data.table)

setwd("/home/dan/Kaggle/UPenn_Seizure/clips")


rf <- fread("../models/randomForest.preds.earlyFix.csv")
rfFreq <-fread("../models/randomForest.freqRatios.preds.earlyFix.csv")
obliqueRF <- fread("../models//oblique_randomForest.freqRatios.preds.csv")
bagMARS <- fread("../models/bagMARS.preds.earlyFix.csv") #NOTE: This hasn't been re-run with early Fix
bagMARS.freq <- fread("../models//bagMARS.preds.freqRatios.csv")
glmNet <- fread("../models//glmnetOnly.preds.earlyFix.csv") #(0.85)
glmNet.freq <-  fread("../models//glmnetOnly.freqRatios.preds.earlyFix.csv") #(0.89)
gbm <- fread("../models//GBM.freqRatios.preds.earlyFix.csv") #(0.90)
svm <- fread("../models/SVM.freqRatios.preds.earlyFix.csv") #0.87

deepBelief.freq <- fread("../models/deepBelief.freqRatios.preds.csv")
rf.h2o <- fread("../models/h2o_randomForest.freqRatios.preds.csv")
elm <- fread("../models/extremeLearningMachine.freqRatios.preds.csv")

rf.h2o.freqRatOnly <- fread(file=paste0("../models/h2o_randomForest.freqRatiosOnly.preds.csv"))
rf.h2o.globalPower <- fread(file=paste0("../models/h2o_randomForest.GlobalPower.preds.csv"))
rf.h2o.globalPowerNorm <- fread(file=paste0("../models/h2o_randomForest.GlobalPowerNormed.preds.csv"))

#Combine model predictions
seizure <- data.table(clip=rf$clip,seizure = (2*rf$seizure + 2*rfFreq$seizure + 2*bagMARS$seizure + 2*bagMARS.freq$seizure + 1*gbm$seizure + 1*glmNet.freq$seizure + 1*svm$seizure) / 11)
early <- data.table(clip=rf$clip,early = (2*rf$early + 2*rfFreq$early + 2*bagMARS$early + 2*bagMARS$early + 1*gbm$seizure + 1*glmNet.freq$early + 1*svm$seizure) / 11)

combined <- merge(seizure,early,by="clip")
write.csv(combined,file="../models/rf_rfFreq_bagMars_bagMars.freq_gbm_glmNetFreq_svm.combined.preds.earlyfix.All.csv",row.names=FALSE)


