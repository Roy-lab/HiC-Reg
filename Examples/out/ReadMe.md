## Exammples of output directory with training models and predictions from HiC-Reg
#### Output Files:
1. testset_error.txt
(Predictions on the test data)
- Column1: Pair of two regions
- Column2: True HiC Count for this pair in the log scale 	
- Column3: Predicted HiC Count for this pair in the log scale 	
- Column4: squared prediction error
- Column5: genomic distance of the two regions

2. trainset_error.txt
(Predictions on the training data)
- Column1: Pair of two regions
- Column2: True HiC Count for this pair in the log scale 	
- Column3: Predicted HiC Count for this pair in the log scale 	
- Column4: squared prediction error
- Column5: genomic distance of the two regions

3. regtree_node_0.txt
(Saved regression tree models for tree 0, examples in Examples/out/)




