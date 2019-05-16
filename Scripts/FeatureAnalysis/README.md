### Tree path analysis:
#### Usage:
```
./regForest_logcnt -t Examples/Data/Gm12878_chr17_WINDOW_train0.txt -o out/ -k1 -l10 -n20 -b Examples/Data/prior_window.txt -d Examples/Data/Gm12878_chr17_WINDOW_test0.txt -s Examples/out/regtree_node
```
#### Arguments:
- -t is the training data
- -o is the output directory
- -k is maxfactorsize (1 should be OK)
- -l is leaf size (e.g. 10)
- -n is the number of trees (e.g. 20)
- -b is the prior structure (potential regulators), see Examples/Data/prior_window.txt 
- -d is the test data
- -s is the prefix of saved models (e.g. Examples/out/regtree_node)
