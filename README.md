# HiC-Reg: In silico prediction of high-resolution Hi-C interaction matrices
![GitHub Logo](/Images/HiC-Reg.png| width=100)


## Step 1: Generate pair features as input for HiC-Reg:
### Step 1.1: Generate PAIR-CONCAT or WINDOW features:
Program in Scripts/genPairFeatures/

#### Usage:
```
./genDatasetsRH HiCsparseMatrix maxdistance foldcv splittype[regionwise|pairwise] featureinputfile correlation[yes|no] outputdir prerandomize_pairs[yes|no] featype[Window|Pconcat]
```

#### Arguments:

- HiCsparseMatrix: sparse hic matrix
- maxdistance: max genomic distance for pair of regions
- foldcv: number of CV folds
- splittype: 
  - regionwise--split the total regions in the sparse hic matrix into N folds. 
  - pairwise--split the total pairs in the sparse hic matrix into N folds.
- featureinputfile: input file with path for each feature signal, see Gm12878_norm_featurefiles_test.txt for example.
- correlation: calculate the correlation of features in region1 and features in region2 or not.
- outputdir: path for output directory.
- prerandomize_pairs: pre-randomize the pairs in the sparse hic matrix or not.
- featype: 
  - Pconcat: generate feature signal for region1 and region2. 
  - Window: generate feature signal for region1 and region2, and average feature signal for the window between these two regions. 

#### Input Files:
1.HiCsparseMatrix: Gm12878_chr17_5kb_SQRTVC_counts_pairs_100.tab

Column1: region1 Column2: region2 Column3: HiC count (tab deliminated)
```
chr17_0_5000	chr17_5000_10000	485.3207854051
chr17_0_5000	chr17_10000_15000	212.4988640932
chr17_0_5000	chr17_15000_20000	130.1059543547
```
2.featureinputfile: Gm12878_norm_featurefiles_test.txt

Column1: mark name Column2: feature signal file Column3: source data type, C for Chip-seq data, M for motif data (tab deliminated)
```
H3k4me1	Gm12878_RawData_5000bp_seqdepth_norm_H3k4me1.txt	C
H3k4me2	Gm12878_RawData_5000bp_seqdepth_norm_H3k4me2.txt	C
```

#### Example: 
```
./genDatasetsRH Gm12878_chr17_5kb_SQRTVC_counts_pairs_100.tab 1000000 5 regionwise Gm12878_norm_featurefiles_test.txt no out/ yes Window
```

### Step 1.2: Generate MULTI-CELL features:



## Step2: Train HiC-Reg models and make predictions
### 1. Training mode:
#### Train on training data and Predict on test set:
```
./regForest -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt
```
### 2. Prediction mode:
#### Predict on test sets with Models trained:
```
./regForest -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt -s out/model/regtree_node
```
#### Arguments:
- -t is the training data
- -o is the output directory
- -k is maxfactorsize (1 should be OK)
- -l is leaf size
- -n is the number of trees
- -b is the prior structure (potential regulators), see prior_merge.txt.
- -d is the test data
- -s is the prefix of saved models

#### Input Files:
1. prior_merge.txt is a file to list the regulators, that is like this:
```
H3k4me1_E    Count
H3k4me1_P    Count
H3k4me1_W    Count
H3k4me2_E    Count
H3k4me2_P    Count
H3k4me2_W    Count
```
2. training and test data:
```
Pair	H3k4me1_E	H3k4me2_E	H3k4me1_P	H3k4me2_P	H3k4me1_W	H3k4me2_W	Distance	Count
chr17_0_5000-chr17_10000_15000	0.714904	0.435131	3.76434	1.3964	0.83982	0.418721	5000	5.36363
chr17_0_5000-chr17_15000_20000	0.714904	0.435131	1.73117	0.869143	2.30208	0.90756	10000	4.87601
chr17_0_5000-chr17_25000_30000	0.714904	0.435131	0.881374	0.708816	1.73773	0.78015	20000	4.76591
chr17_0_5000-chr17_30000_35000	0.714904	0.435131	1.49854	0.760584	1.56646	0.765883	25000	4.5534
```


## Application: Make predictions in a new cell line:
#### Train models:
./regForest -t Data/CrossCell/train0_Gm12878.txt -o out/ -k1 -l10 -n20 -b prior_merge_Gm12878toK562.txt -d Data/CrossCell/test0_Gm12878.txt 

#### Make predictions:
./regForest -t Data/CrossCell/train0_Gm12878.txt -o out/ -k1 -l10 -n20 -b prior_merge_Gm12878toK562.txt -d Data/CrossCell/test0_K562.txt -s out/model/regtree_node





