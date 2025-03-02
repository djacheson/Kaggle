##Files for running various models and generating predictions
---

###NOTES ABOUT MODELS
1. Two different types of models are presented here:
    - Support vector machines with a radial basis function
    - Bagged random forests
2. In both instances, a very large feature space is generated using features generated in ../PreProcessing/ createFile[X].r
3. To reduce the number of features:
    -The SVM model uses variable importance based on earlier instances of running the models
    -The bagged decision trees uses carets selection by filter (SBF) functionality
4. The SVM code allows you to select different feature spaces
5. The Bagged Decision Tree code uses CV folds determined using the extractSequences.R file in ../PreProcessing
6. As usual, many more models were run, but this basic code can easily be modified to specify different types of models (e.g., bagged neural networks, random forests, bagged partial least squares, etc.)
