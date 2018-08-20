#!/bin/sh

## 1. Train on train.txt and Predict on hard test set:
./regForest -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt

## 2. Predict on test sets with Models trained:
./regForest -t train0.txt -o out/ -k1 -l10 -n20 -b prior_merge.txt -d test0.txt -s out/model/regtree_node

