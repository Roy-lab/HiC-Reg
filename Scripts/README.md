## Scripts
### aggregateSignalInRegion/
Aggregate region-level features to prepare input for genPairFeatures

### genPairFeatures/
Generate PAIR-CONCAT or WINDOW features as input for HiC-Reg

### genMULTICELLfeats/
Generate MULTI-CELL features as input for HiC-Reg

### LearnModelsPred/
Script for training and predicting in bash script

### Evaluation/
- CorrelationPlot.m (Run RFdistbinEval.m to make distance-stratified correlation plot and compuete AUC)
- RFdistbinEval.m (Compute distance-stratified correlation, global correlation and MSE)

### FeatureAnalysis/
Tree path analysis and rules counting of how many times a feature or a pair of features are used to predict a test pair
- featureanalysis.sh (bash script to run steps in one)
