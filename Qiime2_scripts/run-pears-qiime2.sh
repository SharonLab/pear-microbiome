#!/bin/bash

conda activate qiime2-2022.8
## the next line is for the bug
export TMPDIR=$HOME/tmp/
##  Import the data - this is done in directory ../qiime2.itai-backup-of-originl-folder-from-server/ by itai. 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --output-path 01.pear-flowers.qza \
  --input-format PairedEndFastqManifestPhred33V2&


qiime demux summarize \
  --i-data 01.pear-flowers.qza \
  --o-visualization 01.pear-flowers.qzv

# Denoise - the trimming parameters are based on the fastqc analysis. 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ~/Fire-Blight/qiime-pseudo/01.pear-flowers.qza \
  --p-trim-left-f 17 \
  --p-trim-left-r 7 \
  --p-trunc-len-f 256 \
  --p-trunc-len-r 256 \
  --p-pooling-method 'pseudo' \
  --o-table 02.table-pseudo.qza \
  --p-n-threads 20 \
  --o-representative-sequences 02.rep-seqs-pseudo.qza \
  --o-denoising-stats 02.denoising-stats-pseudo.qza \
  --verbose&

qiime feature-table summarize \
  --i-table 02.table-pseudo.qza \
  --o-visualization 02.table-pseudo.qzv \
  --m-sample-metadata-file metadata-filtered-copy.tsv


qiime feature-table tabulate-seqs \
  --i-data 02.rep-seqs-pseudo.qza \
  --o-visualization 02.rep-seqs-pseudo.qzv

qiime metadata tabulate \
  --m-input-file 02.denoising-stats-pseudo.qza \
  --o-visualization 02.denoising-stats-pseudo.qzv


# Generate a phylogennetic tree.  This is required for later analyses (e.g. unifrac-based)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 02.rep-seqs-pseudo.qza \
  --o-alignment 03.aligned-rep-seqs.qza \
  --o-masked-alignment 03.masked-aligned-rep-seqs.qza \
  --o-tree 03.unrooted-tree.qza \
  --o-rooted-tree 03.rooted-tree.qza


#####
# Assign taxonomy
# Get the silva classifier
#wget https://data.qiime2.org/2020.8/common/silva-138-99-nb-classifier.qza
# classifier downloaded on 7.12.2022
qiime feature-classifier classify-sklearn \
  --i-classifier feature-classifiers/silva-138-99-nb-classifier.qza \
  --i-reads 02.rep-seqs-pseudo.qza \
  --o-classification 04.taxonomy-silva.qza \
 --p-n-jobs 3


qiime metadata tabulate \
  --m-input-file 04.taxonomy-silva.qza \
  --o-visualization 04.taxonomy-silva.qzv

## filter out samples with low feature count. 
qiime feature-table filter-samples \
  --i-table 02.table-pseudo.qza \
  --p-min-frequency 1500 \
  --o-filtered-table 04.sample-frequency-filtered-table.qza
qiime tools export --input-path 04.sample-frequency-filtered-table.qza --output-path exported
## filter out the mitochondria and chloroplasts
qiime taxa filter-table \
  --i-table 04.sample-frequency-filtered-table.qza \
  --i-taxonomy 04.taxonomy-silva.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table 06.table-no-mitochondria-no-chloroplast-silva.qza  
