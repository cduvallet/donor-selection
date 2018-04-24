#!/usr/bin/env bash

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.hsiao.csv \
  --output-path hsiao.raw-data.qza \
  --source-format PairedEndFastqManifestPhred33

qiime quality-filter q-score \
  --i-demux hsiao.raw-data.qza \
  --o-filtered-sequences hsiao.demux-filtered.qza \
  --o-filter-stats hsiao.demux-filter-stats.qza \
  --verbose    

qiime demux summarize \
  --i-data hsiao.raw-data.qza \
  --o-visualization hsiao.data_qual.qzv

qiime deblur denoise-16S \
  --i-demultiplexed-seqs hsiao.demux-filtered.qza \
  --p-trim-length 200 \
  --o-representative-sequences hsiao.rep-seqs-deblur.qza \
  --o-table hsiao.table-deblur.qza \
  --p-sample-stats \
  --o-stats hsiao.deblur-stats.qza \
  --verbose

qiime tools export \
  hsiao.table-deblur.qza \ 
  --output-dir hsiao_export

biom convert \
  -i hsiao_export/feature-table.biom \
  -o hsiao_export/hsiao.feature_table.txt \
  --to-tsv

qiime tools export \
  hsiao.rep-seqs-deblur.qza \
  --output-dir hsiao_export/

mv hsiao_export/dna-sequences.fasta hsiao_export/hsiao.dna-sequences.fasta
mv hsiao_export/feature-table.biom hsiao_export/hsiao.feature-table.biom
