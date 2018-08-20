#!/bin/sh

CELL=Gm12878
sparseMatrix=${CELL}_normalized/Gm12878_chr17_5kb_SQRTVC_counts_pairs.tab
foldcv=5
splittype=regionwise
prerandomize_pairs=yes
featuretype=Window

#Usage: getFeatures sparseMatrix dist foldcv splittype[regionwise|pairwise] featurefile correlation[yes|no] outputdir prerandomize_pairs[yes|no] featuretype[Window]

./genDatasetsRH $sparseMatrix 1000000 $foldcv $splittype ${CELL}_norm_featurefiles_PIQ_5kb.txt no out/ $prerandomize_pairs $featype

