# HiC-Reg: In silico prediction of high-resolution Hi-C interaction matrices

## Usage:
### 1. Train on train.txt and Predict on hard test set:
./regForest -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt

### 2. Predict on test sets with Models trained:
./regForest -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt -s out/model/regtree_node


-t is the training data
-o is the output directory
-k is maxfactorsize (1 should be OK)
-l is leaf size
-n is the number of trees
-d is the test set
-s is the prefix to saved trees
-b is the prior structure (potential regulators)

prior_merge.txt is a file to list the regulators, that is like this:
CTCF_E    Count
CTCF_P    Count
CTCF_W    Count

