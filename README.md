# HiC-Reg: In silico prediction of high-resolution Hi-C interaction matrices

## Step1: Generate pair features as input for HiC-Reg:
Program in Scripts/genPairFeatures/

#### Usage:
```
./genDatasetsRH HiCsparseMatrix maxdistance foldcv splittype[regionwise|pairwise] featureinputfile correlation[yes|no] outputdir prerandomize_pairs[yes|no] featype[Window|Pconcat|Outerprod]
```

- region
- pair

#### Arguments:

- HiCsparseMatrix: sparse hic matrix
- maxdistance: max genomic distance for pair of regions
- foldcv: number of CV folds
- splittype: 1) regionwise--split the total regions in the sparse hic matrix into N folds. 2) pairwise--split the total pairs in the sparse hic matrix into N folds.
- featureinputfile: input file with path for each feature signal, see Gm12878_norm_featurefiles_PIQ_5kb.txt for example.
- correlation: calculate the correlation of features in region1 and features in region2 or not.
- outputdir: path for output directory.
- prerandomize_pairs: pre-randomize the pairs in the sparse hic matrix or not.
- featype: 1)Pconcat: generate feature signal for region1 and region2. 2) Window: generate feature signal for region1 and region2, and average feature signal for the window between these two regions. 3) Outerprod: calculate the output product of the features in region1 and features in region2.

#### Example: 
./genDatasetsRH Gm12878_chr17_5kb_SQRTVC_counts_pairs.tab 1000000 5 regionwise Gm12878_norm_featurefiles_PIQ_5kb.txt no out/ yes Window



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





