# HiC-Reg: In silico prediction of high-resolution Hi-C interaction matrices

## Usage:
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
