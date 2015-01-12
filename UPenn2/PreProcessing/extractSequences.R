#This file is designed to manually create cross-validation (CV) splits based on whether clips were part of the same sequence or not

#The data was generated where combinations of clips occurred closer in time (i.e., in the same sequence). 
# Taking this grouping into account during CV led to better model tuning
# and more accurate leaderboard performance

setwd('/home/dan/Kaggle/UPenn2/')

library(data.table)
library(stringr)

#Generate csv with different sequences coded by clip, trialtype and trial number
seqs <- fread("./Segment_Sequence.csv") 
g <- lapply(as.list(seqs$Segment), function(x) str_split(x,"[//]")[[1]][[2]])
g <- unlist(g)

g2 <- lapply(g, function(x) t(data.frame(str_split(x[[1]],"_"))))
g2 <- data.frame(Reduce(rbind, g2),row.names=NULL)
g2[[5]] <- as.numeric(str_replace(g2[[5]],".mat",""))

File <- paste0(g2[[1]],"_",g2[[2]])
trialType <- paste0(g2[[3]])
trialNum <- paste0(g2[[5]])

allSeqs <- data.table(File = File,
                      trialType = trialType,
                      trialNum = trialNum,
                      Sequence = seqs$Sequence)

write.csv(allSeqs, file = "Sequence_Processed.csv", row.names = F)


allSeqs[,length(unique(Sequence)),by = c("File","trialType")]

files <- unique(allSeqs$File)

#Create data.table with which CV fold according to clip ("File") and trialType
#Done separately for Dogs and Patients due to differing numbers of clips (more CV folds for dogs)
#Later merged with all sequences

Dog_CV <- list()
count = 1
for(i in files[1:5]) {
  interict <- allSeqs[File == i & trialType == "interictal", unique(Sequence), by = c("File","trialType")]
  setnames(interict,"V1","Sequence")
  interict$CV <- rep(1:4,length.out = nrow(interict))
  preict <- allSeqs[File == i & trialType == "preictal", unique(Sequence), by = c("File","trialType")]
  setnames(preict,"V1","Sequence")
  preict$CV <- rep(1:4,length.out = nrow(preict))
  Dog_CV[[i]] = rbind(preict,interict)
}


Patient_CV <- list()
for(i in files[6:7]) {
  interict <- allSeqs[File == i & trialType == "interictal", unique(Sequence), by = c("File","trialType")]
  setnames(interict,"V1","Sequence")
  interict$CV <- rep(1:3,length.out = nrow(interict))
  preict <- allSeqs[File == i & trialType == "preictal", unique(Sequence), by = c("File","trialType")]
  setnames(preict,"V1","Sequence")
  preict$CV <- rep(1:3,length.out = nrow(preict))
  Patient_CV[[i]] = rbind(preict,interict)
}

allCV <- rbind(rbindlist(Dog_CV),rbindlist(Patient_CV))

#Merge data.tables with sequences and CV fold number
allSeqs2 <- merge(allSeqs, allCV, by = c("File","trialType","Sequence"))

write.csv(allSeqs, file = "./CV_splits.csv", row.names = F)


