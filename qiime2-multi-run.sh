#!/bin/bash
#This script Kira Goff, 2025. Based on https://amplicon-docs.qiime2.org/en/latest/ (retrieved August 2025)
#Uses QIIME2 version 2025.4
set -e

#Every sequencing run MUST be run separately through the denoising step
#You can merge it after
#This script uses an array to process and then combine data from sequencing runs
#Modify as needed for your data

#For demultiplexed samples, you will need a sample manifest for each sequencing run
#Example format for a sample-manifest.tsv
#sample-id forward-absolute-filepath   reverse-absolute-filepath #header
#sample-1  /scratch/microbiome/sample1_R1.fastq.gz  /scratch/microbiome/sample1_R2.fastq.gz
#sample-2  /scratch/microbiome/sample2_R1.fastq.gz  /scratch/microbiome/sample2_R2.fastq.gz

conda activate qiime2

array=(01 02 03 04 05 06 07 08 09 10 11)
for i in "${array[@]}";
do 
qiime tools import \
    --type "SampleData[PairedEndSequencesWithQuality]" \
    --input-path manifests/manifest-"$i".tsv \
    --output-path run-"$i"/demux-"$i".qza \
    --input-format PairedEndFastqManifestPhred33V2 ;   

qiime demux summarize \
    --i-data run-"$i"/demux-"$i".qza \
    --o-visualization run-"$i"/demux-"$i".qzv ;

#adjust trim parameters based on your primer lengths
#adjust trunc parameters based on the length of your region, leaving enough overlap to merge
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs run-"$i"/demux-"$i".qza \
    --p-trim-left-f 19 --p-trim-left-r 20 \
    --p-trunc-len-f 240 --p-trunc-len-r 190 \
    --p-pooling-method pseudo \
    --p-n-threads 14 \
    --o-table run-"$i"/table-"$i".qza \
    --o-representative-sequences run-"$i"/rep-seq-"$i".qza \
    --o-denoising-stats run-"$i"/stats-"$i".qza ;

qiime metadata tabulate \
    --m-input-file run-"$i"/stats-"$i".qza \
    --o-visualization run-"$i"/stats-"$i".qzv    
done    

#merge seqs after denoising

qiime feature-table merge \
    --i-tables run-01/table-01.qza  run-02/table-02.qza  run-03/table-03.qza  run-04/table-04.qza  run-05/table-05.qza  run-06/table-06.qza  run-07/table-07.qza  run-08/table-08.qza  run-09/table-09.qza run-10/table-10.qza run-11/table-11.qza --o-merged-table merged-table.qza
qiime feature-table merge-seqs \
    --i-data run-01/rep-seq-01.qza run-02/rep-seq-02.qza run-03/rep-seq-03.qza run-04/rep-seq-04.qza run-05/rep-seq-05.qza run-06/rep-seq-06.qza run-07/rep-seq-07.qza run-08/rep-seq-08.qza run-09/rep-seq-09.qza run-10/rep-seq-10.qza run-11/rep-seq-11.qza --o-merged-data merged-rep-seqs.qza

qiime feature-table summarize \
    --i-table merged-table.qza \
    --o-visualization merged-table.qzv \
    --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
--i-data merged-rep-seqs.qza \
--o-visualization merged-rep-seqs.qzv

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
    --i-table merged-table.qza \
    --p-min-frequency 20000 \
    --o-filtered-table merged-table-filtered.20k.qza

qiime feature-table summarize \
    --i-table merged-table-filtered.20k.qza \
    --o-visualization merged-table-filtered.20k.qzv \
    --m-sample-metadata-file metadata.tsv

qiime feature-table filter-seqs \
    --i-data merged-rep-seqs.qza \
    --i-table merged-table-filtered.20k.qza \
    --o-filtered-data merged-rep-seqs.20k.qza 

qiime feature-table tabulate-seqs \
    --i-data merged-rep-seqs.20k.qza \
    --o-visualization merged-rep-seqs.20k.qzv

#Using a 39.8k read cutoff
qiime feature-table filter-samples \
    --i-table merged-table.qza \
    --p-min-frequency 39800 \
    --o-filtered-table merged-table-filtered.39k.qza

qiime feature-table summarize \
    --i-table merged-table-filtered.39k.qza \
    --o-visualization merged-table-filtered.39k.qzv \
    --m-sample-metadata-file metadata.tsv

qiime feature-table filter-seqs \
    --i-data merged-rep-seqs.qza \
    --i-table merged-table-filtered.39k.qza \
    --o-filtered-data merged-rep-seqs.39k.qza 

qiime feature-table tabulate-seqs \
    --i-data merged-rep-seqs.39k.qza \
    --o-visualization merged-rep-seqs.39k.qzv

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
#Here I'm just doing the 39k object
qiime feature-classifier classify-sklearn \
    --i-classifier classifier/silva-138.2-ssu-nr99-515f-926r-classifier.qza \
    --i-reads merged-rep-seqs.39k.qza \
    --o-classification taxonomy.qza
#Add your metadata    
qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv

#Make your bar plots, then open them in the qiime2 viewer and play around!
qiime taxa barplot \
    --i-table merged-table-filtered.39k.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization taxa-bar-plots.qzv

#let's get the format to match dada2 output - only necessary if you are going to be comparing with DADA2 output.
qiime tools export --input-path merged-table-filtered.39k.qza --output-path exported
qiime tools export --input-path taxonomy.qza --output-path exported