#!/bin/sh

### 2.1 Training mode:
#### Train on training data and Predict on test set:
./regForest -t Examples/Data/Gm12878_chr17_WINDOW_train0.txt -o Examples/out/ -k1 -l10 -n20 -b Examples/Data/prior_window.txt -d Examples/Data/Gm12878_chr17_WINDOW_test0.txt

### 2.2 Prediction mode:
#### Predict on test sets with Models trained:

./regForest -t Examples/Data/Gm12878_chr17_WINDOW_train0.txt -o out/ -k1 -l10 -n20 -b Examples/Data/prior_window.txt -d Examples/Data/Gm12878_chr17_WINDOW_test0.txt -s Examples/out/regtree_node
