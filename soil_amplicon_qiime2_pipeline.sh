#!/bin/bash
# QIIME2 soil amplicon analysis pipeline

###@email yanhuizhen@westlake.edu.cn

# Import raw sequences into QIIME2
cd /home/yanhuizhen/soil_amplicon
time qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path metadata.txt \
  --output-path paired-seqs.qza \
  --input-format PairedEndFastqManifestPhred33

# Denoising with DADA2
mkdir asv
time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-seqs.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 212 \
  --o-representative-sequences asv/rep-seqs-dada2.qza \
  --o-table asv/table_dada2.qza \
  --o-denoising-stats asv/denoising_stats_dada2.qza \
  --p-n-threads 20 \
  --verbose

# Download SILVA database and prepare classifier
mkdir silva138db
cd silva138db
wget https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza
wget https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza

# Extract V4-V5 region
time qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer CCGYCAATTYMTTTRAGTTT \
  --p-n-jobs 20 \
  --o-reads silva-138-99-seqs_v4v5.qza

# Train naive Bayes classifier
time qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-seqs_v4v5.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99_v4v5_classifier.qza

# Assign taxonomy to ASVs
cd ..
time qiime feature-classifier classify-sklearn \
  --i-classifier silva138db/silva-138-99_v4v5_classifier.qza \
  --i-reads asv/rep-seqs-dada2.qza \
  --o-classification asv/rep-seqs-dada2_taxonomy.qza \
  --p-n-jobs 0

# OTU clustering at 97% similarity
mkdir otu
time qiime vsearch cluster-features-de-novo \
  --i-table asv/table_dada2.qza \
  --i-sequences asv/rep-seqs-dada2.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table otu/table-97.qza \
  --o-clustered-sequences otu/rep-seqs-97.qza \
  --p-threads 40 \
  --verbose

# Filter OTU table and sequences by abundance
time qiime feature-table filter-features \
  --i-table otu/table-97.qza \
  --p-min-frequency 2 \
  --p-min-samples 2 \
  --o-filtered-table otu/filtered-table-97.qza \
  --verbose

time qiime feature-table filter-seqs \
  --i-data otu/rep-seqs-97.qza \
  --i-table otu/filtered-table-97.qza \
  --o-filtered-data otu/filtered-rep-seqs-97.qza

# Remove mitochondria, chloroplast, unassigned, and eukaryota
mkdir filter_UECM
time qiime taxa filter-seqs \
  --i-sequences otu/filtered-rep-seqs-97.qza \
  --i-taxonomy asv/rep-seqs-dada2_taxonomy.qza \
  --p-exclude Mitochondria,Chloroplast,Unassigned,Eukaryota \
  --o-filtered-sequences filter_UECM/sequences-with-phyla-no-mitochondria-no-chloroplast.qza

qiime taxa filter-table \
  --i-table otu/filtered-table-97.qza \
  --i-taxonomy asv/rep-seqs-dada2_taxonomy.qza \
  --p-exclude Mitochondria,Chloroplast,Unassigned,Eukaryota \
  --o-filtered-table filter_UECM/table-with-phyla-no-mitochondria-no-chloroplast.qza

# Summarize table
qiime feature-table summarize \
  --i-table filter_UECM/table-with-phyla-no-mitochondria-no-chloroplast.qza \
  --o-visualization filter_UECM/table-with-phyla-no-mitochondria-no-chloroplast.qzv

# Build phylogenetic tree
mkdir phylogeny
time qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads auto \
  --i-sequences filter_UECM/sequences-with-phyla-no-mitochondria-no-chloroplast.qza \
  --o-alignment phylogeny/sequences_aligned.qza \
  --o-masked-alignment phylogeny/sequences_mask.qza \
  --o-tree phylogeny/sequences_unrooted-tree.qza \
  --o-rooted-tree phylogeny/sequences_rooted-tree.qza

# Calculate alpha and beta diversity metrics
mkdir diversity
time qiime diversity core-metrics-phylogenetic \
  --i-table filter_UECM/table-with-phyla-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny phylogeny/sequences_rooted-tree.qza \
  --p-sampling-depth 33130 \
  --m-metadata-file metadata-forq2divsity.tsv \
  --p-n-jobs-or-threads 20 \
  --o-rarefied-table diversity/table_noUECM_rare33130.qza \
  --o-faith-pd-vector diversity/rare33130_adiv_pd.qza \
  --o-observed-features-vector diversity/rare33130_adiv_richness.qza \
  --o-shannon-vector diversity/rare33130_adiv_shannon.qza \
  --o-evenness-vector diversity/rare33130_adiv_evenness.qza \
  --o-unweighted-unifrac-distance-matrix diversity/rare33130_bdiv_unweightedunifracdm.qza \
  --o-weighted-unifrac-distance-matrix diversity/rare33130_bdiv_weightedunifracdm.qza \
  --o-jaccard-distance-matrix diversity/rare33130_bdiv_jaccarddm.qza \
  --o-bray-curtis-distance-matrix diversity/rare33130_bdiv_bcdm.qza \
  --o-unweighted-unifrac-pcoa-results diversity/rare33130_bdiv_unweightedunifracpcoa.qza \
  --o-weighted-unifrac-pcoa-results diversity/rare33130_bdiv_weightedunifracpcoa.qza \
  --o-jaccard-pcoa-results diversity/rare33130_bdiv_jaccardpcoa.qza \
  --o-bray-curtis-pcoa-results diversity/rare33130_bdiv_bcpcoa.qza \
  --o-unweighted-unifrac-emperor diversity/rare33130_bdiv_unweightedunifracemp.qzv \
  --o-weighted-unifrac-emperor diversity/rare33130_bdiv_weightedunifracemp.qzv \
  --o-jaccard-emperor diversity/rare33130_bdiv_jaccardemp.qzv \
  --o-bray-curtis-emperor diversity/rare33130_bdiv_bcemp.qzv

# Collapse taxonomy at the phylum level
mkdir taxa_collapse
qiime taxa collapse \
  --i-table diversity/table_noUECM_rare33130.qza \
  --i-taxonomy asv/rep-seqs-dada2_taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table taxa_collapse/table_noUECM_rare33130-phylum.qza
