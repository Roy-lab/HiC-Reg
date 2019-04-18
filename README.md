# HiC-Reg: In silico prediction of high-resolution Hi-C interaction matrices

## Generate pair features as input for HiC-Reg:
Program in Scripts/genPairFeatures/

foldcv=5

splittype=regionwise

prerandomize_pairs=yes

featype=Window

./genDatasetsRH Gm12878_chr17_5kb_SQRTVC_counts_pairs.tab 1000000 $foldcv $splittype Gm12878_norm_featurefiles_PIQ_5kb.txt no out/ $prerandomize_pairs $featype



## Run HiC-Reg:
### 1. Train on train.txt and Predict on hard test set:
./regForest -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt

### 2. Predict on test sets with Models trained:
./regForest -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt -s out/model/regtree_node

### Parameters:
-t is the training data

-o is the output directory

-k is maxfactorsize (1 should be OK)

-l is leaf size

-n is the number of trees

-d is the test data

-s is the prefix of saved models

-b is the prior structure (potential regulators)

prior_merge.txt is a file to list the regulators, that is like this:

CTCF_E    Count

CTCF_P    Count

CTCF_W    Count

### 3. Make predictions in a new cell line:
#### Train models:
./regForest -t Data/CrossCell/train0_Gm12878.txt -o out/ -k1 -l10 -n20 -b prior_merge_Gm12878toK562.txt -d Data/CrossCell/test0_Gm12878.txt 

#### Make predictions:
./regForest -t Data/CrossCell/train0_Gm12878.txt -o out/ -k1 -l10 -n20 -b prior_merge_Gm12878toK562.txt -d Data/CrossCell/test0_K562.txt -s out/model/regtree_node



## Memory Usage
![alt text](https://github.com/Roy-lab/HiC-Reg/blob/master/Images/Memory_usage_Window_Merge_DepNormPerCell_Gm12878_summary.png "Gm12878")
