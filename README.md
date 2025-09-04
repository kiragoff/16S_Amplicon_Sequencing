# 16S_Amplicon_Sequencing
R scripts for processing and analyzing 16S amplicon data with DADA2, and workflows for QIIME2 processing.

This repository contains scripts and workflows for processing 16S amplicon sequencing data. 
There are different scripts and workflows for handling data from a single sequencing run vs doing combined analysis of data from multiple sequencing runs. 
Each sequencing run has distinct error patterns, so combining results runs results in error correction profiles that introduce more errors. 
For ASV output, use the R-based DADA2 scripts. For OTU output, use the QIIME2 workflows. If you plan to collapse the genus (for environmental data or other reasons), choose whichever you'd like.

The R scripts contain several places where I say something has to be formatted a certain way. If you're comfortable changing the relevant R code, feel free to play around.

For variable-length regions (18S, fungal ITS), use the scripts and workflows from the variable-length region repository.
