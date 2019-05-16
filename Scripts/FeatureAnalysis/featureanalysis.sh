#!/bin/sh

## Step1: Generate tree path for predicting a test pair 
./regForest_logcnt -t Examples/Data/Gm12878_chr17_WINDOW_train0.txt -o out/ -k1 -l10 -n20 -b Examples/Data/prior_window.txt -d Examples/Data/Gm12878_chr17_WINDOW_test0.txt -s Examples/out/regtree_node

## Step2: Use the output testset_featureanalysis.txt from Step 1 to count how many times a feature is used to predict a test pair
level=1
./analyzePathv3 testset_featureanalysis.txt $level feature.txt testset_featureanalysis_f${level}.txt

## Use the output testset_featureanalysis.txt from Step 1 to count how many times a pair of features are used to predict a test pair
level=2
./analyzePathv3 testset_featureanalysis.txt $level feature.txt testset_featureanalysis_f${level}.txt





