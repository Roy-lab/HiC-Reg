# Aggregate region-level features

#### Usage:
```
./aggregateSignal REGION hg19.fa.fai countfile outputfile
```
#### Examples:
```
./aggregateSignal hg19_5kbp_chr17.txt hg19.fa.fai wgEncodeBroadHistoneGm12878CtcfStdRawDataRep1_chr17.counts wgEncodeBroadHistoneGm12878CtcfStdRawDataRep1_chr17.txt
```

#### Arguments:
- REGION: file of 5kb bin with format like hg19_5kbp_chr17.txt
- countfile: four columns count file generated from bam files
  
#### Input Files:
1. hg19_5000bp_chr17.txt
```
chr17	CONVERT	gene	0	4999	.	+	.	chr17_0_4999
chr17	CONVERT	gene	5000	9999	.	+	.	chr17_5000_9999
chr17	CONVERT	gene	10000	14999	.	+	.	chr17_10000_14999
```
2. count file: wgEncodeBroadHistoneGm12878CtcfStdRawDataRep1_chr17.counts 
```
chr17	342	389	1
chr17	1005	1056	1
chr17	1060	1111	1
chr17	1346	1397	1
```
Column1: chromosome Column2: start position Column3: end position Column4: reads (tab deliminated)



  
