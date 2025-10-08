#!/bin/bash
#This script Kira Goff, 2025. Based on https://amplicon-docs.qiime2.org/en/latest/ (retrieved August 2025)
#Uses QIIME2 version 2025.4
set -e

#Every sequencing run MUST be run separately through the denoising step
#You can merge it after
#This script is for a single sequencing run, qiime2-multi-run.sh combines output from multiple runs

#For demultiplexed samples, you will need a sample manifest for each sequencing run
#Example format for a sample-manifest.tsv
#sample-id forward-absolute-filepath   reverse-absolute-filepath #header
#sample-1  /scratch/microbiome/sample1_R1.fastq.gz  /scratch/microbiome/sample1_R2.fastq.gz
#sample-2  /scratch/microbiome/sample2_R1.fastq.gz  /scratch/microbiome/sample2_R2.fastq.gz

conda activate qiime2

qiime metadata tabulate \
    --m-input-file meta1.tsv \
    --o-visualization meta1-viz.qzv

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path sample-manifest.tsv \
    --output-path demux.qza \
    --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
    --i-data demux.qza \
    --o-visualization demux.qzv

#adjust trim parameters based on your primer lengths
#adjust trunc parameters based on the length of your region, leaving enough overlap to merge
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux.qza \
    --p-trim-left-f 19 --p-trim-left-r 20 \
    --p-trunc-len-f 240 --p-trunc-len-r 200 \
    --p-pooling-method pseudo \
    --p-n-threads 12 \
    --o-table table1.qza --o-representative-sequences rep-seq1.qza --o-denoising-stats stats1.qza

qiime metadata tabulate \
    --m-input-file stats1.qza \
    --o-visualization stats1.qzv

qiime feature-table summarize \
    --i-table table1.qza \
    --m-sample-metadata-file meta1.tsv \
    --o-visualization table1.qzv

qiime feature-table tabulate-seqs \
    --i-data rep-seq1.qza \
    --o-visualization rep-seq1.qzv

qiime metadata tabulate \
    --m-input-file metadata.tsv \
    --o-visualization metadata-viz.qzv

#qiime2 subsamples data: whatever you use as your sampling depth will be the # of reads it randomly subsamples
#You want to keep as many reads per sample while least #of samples possible
#Example output from my data use using a variety of cutoff thresholds
#   cutoff  
#   22,000  Retained 1,958,000 (29.48%) observations in 89 (97.80%) samples  
#   39,800  Retained 3,502,400 (52.74%) observations in 88 (96.70%) samples
#   43,500  Retained 3,784,500 (56.98%) observations in 87 (95.60%) samples
#   49,500  Retained 4,059,000 (61.12%) observations in 82 (90.11%) samples  

#Below I generate two objects using two different cutoff thresholds for comparison

#Using a 20k read cutoff
qiime feature-table filter-samples \
    --i-table table1.qza \
    --p-min-frequency 20000 \
    --o-filtered-table table1-filtered.20k.qza

qiime feature-table summarize \
    --i-table table1-filtered.20k.qza \
    --o-visualization table1-filtered.20k.qzv \
    --m-sample-metadata-file metadata.tsv

qiime feature-table filter-seqs \
    --i-data rep-seq1.qza \
    --i-table table1-filtered.20k.qza \
    --o-filtered-data rep-seq1.20k.qza 

qiime feature-table tabulate-seqs \
    --i-data rep-seq1.20k.qza \
    --o-visualization rep-seq1.20k.qzv

#Using a 39.8k read cutoff
qiime feature-table filter-samples \
    --i-table table1.qza \
    --p-min-frequency 39800 \
    --o-filtered-table table1-filtered.39k.qza

qiime feature-table summarize \
    --i-table table1-filtered.39k.qza \
    --o-visualization table1-filtered.39k.qzv \
    --m-sample-metadata-file metadata.tsv

qiime feature-table filter-seqs \
    --i-data rep-seq1.qza \
    --i-table table1-filtered.39k.qza \
    --o-filtered-data rep-seq1.39k.qza 

qiime feature-table tabulate-seqs \
    --i-data rep-seq1.39k.qza \
    --o-visualization rep-seq1.39k.qzv

#Train the silva classifier
qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza

#get rid of ugly sequences 
qiime rescript cull-seqs \
    --i-sequences silva-138.2-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza

#get rid of sequences outside expected threshold
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza 

#Extract your 16S region using your primers
qiime feature-classifier extract-reads \
    --i-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer GTGYCAGCMGCCGCGGTAA \
    --p-r-primer CCGYCWATTYMTTTRAGTTT \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.2-ssu-nr99-seqs-515f-926r.qza

#Dereplicate sequences
qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-515f-926r.qza \
    --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-515f-926r-uniq.qza \
    --o-dereplicated-taxa  silva-138.2-ssu-nr99-tax-515f-926r-derep-uniq.qza

#Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.2-ssu-nr99-seqs-515f-926r-uniq.qza \
    --i-reference-taxonomy silva-138.2-ssu-nr99-tax-515f-926r-derep-uniq.qza \
    --o-classifier silva-138.2-ssu-nr99-515f-926r-classifier.qza

#Classify your sequences
#Here I'm just doing the 39k filter
qiime feature-classifier classify-sklearn \
    --i-classifier classifier/silva-138.2-ssu-nr99-515f-926r-classifier.qza \
    --i-reads rep-seq1-filtered.39k.qza \
    --o-classification taxonomy.qza
#Add your metadata    
qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv

#Make your bar plots, then open them in the qiime2 viewer and play around!
qiime taxa barplot \
    --i-table table1-filtered.39k.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization taxa-bar-plots.qzv

#let's get the format to match dada2 output - only necessary if you are going to be comparing with DADA2 output.
qiime tools export --input-path table1-filtered.39k.qza --output-path exported
qiime tools export --input-path taxonomy.qza --output-path exported