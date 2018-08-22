#!/bin/sh

## Step1: Generate tree path for predicting a test pair 
echo ./regForest_logcnt -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt -s out/model/regtree_node
./regForest_logcnt -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt -s out/model/regtree_node

## Step2: Use the output testset_featureanalysis.txt from Step 1 to count how many times a feature or a pair of features are used to predict a test pair
filter=0  # 0:no fitering or 1:fitering based on $pairs
level=1 #1: individual feature analysis 2:pair-wise feature analysis
Error=all  #all or bottom5high or top5low
output=0 ## output testset_featureanalysis_${Error}.txt or not
pairs=pairs_5lowerror_hard.txt
echo ./analyzePathv3 testset_featureanalysis.txt $level $feature testset_featureanalysis_path_f${level}_${Error}err.txt $filter $pairs $output testset_featureanalysis_${Error}.txt
./analyzePathv3 testset_featureanalysis.txt $level feature.txt testset_featureanalysis_path_f${level}_${Error}err.txt $filter $pairs $output testset_featureanalysis_${Error}.txt

