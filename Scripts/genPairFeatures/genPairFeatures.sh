sparseMatrix=${CELL}_normalized/SQRTVC/${CHR}_${res}kb_SQRTVC_counts_pairs.tab
#Usage: getFeatures sparseMatrix dist foldcv splittype[regionwise|pairwise] featurefile correlation[yes|no] outputdir prerandomize_pairs[yes|no] featuretype[Window]
./genDatasetsRH $sparseMatrix 1000000 $foldcv $splittype ${CELL}_norm_featurefiles_PIQ_${res}kb.txt no $out/ $prerandomize_pairs $featype

