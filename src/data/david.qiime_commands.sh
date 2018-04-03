#!/usr/bin/env bash

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../all_raw_gzip_manifest.csv \
  --output-path david_raw_seqs_from_gzip.qza \
  --source-format PairedEndFastqManifestPhred33

qiime quality-filter q-score \
  --i-demux david_raw_seqs_from_gzip.qza \
  --o-filtered-sequences david.demux-filtered.qza \
  --o-filter-stats david.demux-filter-stats.qza \
  --verbose    

qiime deblur denoise-16S \
  --i-demultiplexed-seqs david.demux-filtered.qza \
  --p-trim-length 95 \
  --o-representative-sequences david.rep-seqs-deblur.qza \
  --o-table david.table-deblur.qza \
  --p-sample-stats \
  --o-stats david.deblur-stats.qza \
  --verbose

qiime tools export \
  david.table-deblur.qza 
  --output-dir david_export

biom convert \
  -i david_export/feature-table.biom \
  -o david_export/david.feature_table.txt \
  --to-tsv

qiime tools export \
  david.rep-seqs-deblur.qza \
  --output-dir david_export/
