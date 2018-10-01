## Make sure you're in a QIIME 2 environment when
## running these commands!
# source activate qiime2-2018.6

# Import data into qiime on aws:
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path jacob2017.manifest.csv \
  --output-path jacob2017.paired-end-demux.qza \
  --source-format PairedEndFastqManifestPhred33

# Summarize data to visualize quality, etc
qiime demux summarize \
  --i-data jacob2017.paired-end-demux.qza \
  --o-visualization jacob2017.demux.qzv

## Optional: visualize quality plots
# qiime tools view demux.qzv

# Quality filter data (with default values)
qiime quality-filter q-score \
  --i-demux jacob2017.paired-end-demux.qza \
  --o-filtered-sequences jacob2017.demux-quality-filtered.qza \
  --o-filter-stats jacob2017.demux-quality-filter-stats.qza

# Convert quality filtered data to visualization
qiime metadata tabulate \
  --m-input-file jacob2017.demux-quality-filter-stats.qza \
  --o-visualization jacob2017.demux-quality-filter-stats.qzv

# denoise sequences with deblur
qiime deblur denoise-16S \
  --i-demultiplexed-seqs jacob2017.demux-quality-filtered.qza \
  --p-trim-length 220 \
  --o-representative-sequences jacob2017.rep-seqs-deblur.qza \
  --o-table jacob2017.table-deblur.qza \
  --p-sample-stats \
  --o-stats jacob2017.deblur-stats.qza

# Convert deblur stats to visualization
qiime deblur visualize-stats \
  --i-deblur-stats jacob2017.deblur-stats.qza \
  --o-visualization jacob2017.deblur-stats.qzv

# Assign taxonomy

# First, you'll need to prep the classifier
# See prep_qiime_classifier.sh script for eventual makefile

# TODO: change path to classifier to something smarter/better
qiime feature-classifier classify-sklearn \
  --i-reads jacob2017.rep-seqs-deblur.qza \
  --i-classifier ../../silva-132-99-nb-classifier.qza \
  --o-classification jacob2017.taxonomy.silva-132-99.qza

# Classify with GG
qiime feature-classifier classify-sklearn \
  --i-reads jacob2017.rep-seqs-deblur.qza \
  --i-classifier ../../gg-13-8-99-nb-classifier.qza \
  --o-classification jacob2017.taxonomy.gg-13-8-99.qza

# export data
qiime tools export \
  jacob2017.taxonomy.silva-132-99.qza \
  --output-dir exported_data

mv exported_data/taxonomy.tsv exported_data/jacob2017.taxonomy.silva-132-99.tsv

qiime tools export \
  jacob2017.taxonomy.gg-13-8-99.qza \
  --output-dir exported_data

mv exported_data/taxonomy.tsv exported_data/jacob2017.taxonomy.gg-13-8-99.tsv

qiime tools export \
  jacob2017.table-deblur.qza \
  --output-dir exported_data

mv exported_data/feature-table.biom exported_data/jacob2017.feature-table.biom

biom convert \
  -i exported_data/jacob2017.feature-table.biom \
  -o exported_data/jacob2017.feature-table.txt \
  --to-tsv
