### Tree path analysis:
#### Usage:
```
./regForest_logcnt -t Examples/Data/Gm12878_chr17_WINDOW_train0.txt -o out/ -k1 -l10 -n20 -b Examples/Data/prior_window.txt -d Examples/Data/Gm12878_chr17_WINDOW_test0.txt -s Examples/out/regtree_node
```
#### Arguments (sames as regForest): 
- -t is the training data
- -o is the output directory
- -k is maxfactorsize (1 should be OK)
- -l is leaf size (e.g. 10)
- -n is the number of trees (e.g. 20)
- -b is the prior structure (potential regulators), see Examples/Data/prior_window.txt 
- -d is the test data
- -s is the prefix of saved models (e.g. Examples/out/regtree_node)


#### Input Files (sames as regForest)
1. training and test data genreated by Step 1: (tab deliminated)
```
Pair	H3k4me1_E	H3k4me2_E	H3k4me1_P	H3k4me2_P	H3k4me1_W	H3k4me2_W	Distance	Count
chr17_0_5000-chr17_10000_15000	0.714904	0.435131	3.76434	1.3964	0.83982	0.418721	5000	5.36363
chr17_0_5000-chr17_15000_20000	0.714904	0.435131	1.73117	0.869143	2.30208	0.90756	10000	4.87601
chr17_0_5000-chr17_25000_30000	0.714904	0.435131	0.881374	0.708816	1.73773	0.78015	20000	4.76591
chr17_0_5000-chr17_30000_35000	0.714904	0.435131	1.49854	0.760584	1.56646	0.765883	25000	4.5534
```
header of training and test data should be the same.

2. prior_window.txt is a file to list the feature set used for training models, that is like this: (tab deliminated)
```
H3k4me1_E    Count
H3k4me1_P    Count
H3k4me1_W    Count
H3k4me2_E    Count
H3k4me2_P    Count
H3k4me2_W    Count
Distance    Count
```
Column1: feature set (header of column 2 to second last column) in training/test data, like H3k4me1_E, H3k4me1_P, H3k4me1_W, Distance.

Column2: Count (i.e. the header of the last column in training/test data)

#### Output Files:
1. testset_featureanalysis.txt (tree path rules)
```
chr17_10035000_10040000-chr17_10070000_10075000	30000		H3k27ac_W>	Distance<=	TBP_P>	3.1667	4.28671
chr17_10035000_10040000-chr17_10070000_10075000	30000		Dnase_W>	Ctcf_W<=	Dnase_P<=	H3k9ac_P<=	H3k9me3_P>	H3k4me1_P>	3.52945	4.28671
chr17_10035000_10040000-chr17_10070000_10075000	30000		H3k79me2_W>	H3k9ac_P>	H3k4me2_W>	H3k4me1_P>	3.27295	4.28671
```
- Column1: Pair of two regions
- Column2: genomic distance of the two regions
- Column 3 to n-2: feature rules for spliting the node, e.g. H3k27ac_W> (H3k27ac_W is used to split the node)
- Second last column (n-1): Predicted HiC Count for this pair for tree i in the log scale 
- Last Column (n): True HiC Count for this pair in the log scale 	


### Rules counting:
#### Usage:
```
./analyzePathv3 testset_featureanalysis.txt $level feature.txt $outputfile 0 pairs_5LowHigherror_hard.txt 0 ${featurepath}/out.txt
```