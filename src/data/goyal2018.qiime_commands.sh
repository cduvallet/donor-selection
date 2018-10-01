## Make sure you're in a QIIME 2 environment when
## running these commands!
# source activate qiime2-2018.6

dataset=goyal2018
# Update according to sequence qualities (after demux summarize step)
length=220

# Import data into qiime on aws:
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ${dataset}.manifest.csv \
  --output-path ${dataset}.demux.qza \
  --source-format SingleEndFastqManifestPhred33

# Summarize data to visualize quality, etc
qiime demux summarize \
  --i-data ${dataset}.demux.qza \
  --o-visualization ${dataset}.demux.qzv

## Optional: visualize quality plots
# qiime tools view demux.qzv

# Quality filter data (with default values)
qiime quality-filter q-score \
  --i-demux ${dataset}.demux.qza \
  --o-filtered-sequences ${dataset}.demux-quality-filtered.qza \
  --o-filter-stats ${dataset}.demux-quality-filter-stats.qza

# Convert quality filtered data to visualization
qiime metadata tabulate \
  --m-input-file ${dataset}.demux-quality-filter-stats.qza \
  --o-visualization ${dataset}.demux-quality-filter-stats.qzv

# denoise sequences with deblur
qiime deblur denoise-16S \
  --i-demultiplexed-seqs ${dataset}.demux-quality-filtered.qza \
  --p-trim-length ${length} \
  --o-representative-sequences ${dataset}.rep-seqs-deblur.qza \
  --o-table ${dataset}.table-deblur.qza \
  --p-sample-stats \
  --o-stats ${dataset}.deblur-stats.qza

# Convert deblur stats to visualization
qiime deblur visualize-stats \
  --i-deblur-stats ${dataset}.deblur-stats.qza \
  --o-visualization ${dataset}.deblur-stats.qzv

# Assign taxonomy

# First, you'll need to prep the classifier
# See prep_qiime_classifier.sh script for eventual makefile

# TODO: change path to classifier to something smarter/better
qiime feature-classifier classify-sklearn \
  --i-reads ${dataset}.rep-seqs-deblur.qza \
  --i-classifier ../../silva-132-99-nb-classifier.qza \
  --o-classification ${dataset}.taxonomy.silva-132-99.qza

# Classify with GG
qiime feature-classifier classify-sklearn \
  --i-reads ${dataset}.rep-seqs-deblur.qza \
  --i-classifier ../../gg-13-8-99-nb-classifier.qza \
  --o-classification ${dataset}.taxonomy.gg-13-8-99.qza

# export data
qiime tools export \
  ${dataset}.taxonomy.silva-132-99.qza \
  --output-dir exported_data

mv exported_data/taxonomy.tsv exported_data/${dataset}.taxonomy.silva-132-99.tsv

qiime tools export \
  ${dataset}.taxonomy.gg-13-8-99.qza \
  --output-dir exported_data

mv exported_data/taxonomy.tsv exported_data/${dataset}.taxonomy.gg-13-8-99.tsv

qiime tools export \
  ${dataset}.table-deblur.qza \
  --output-dir exported_data

mv exported_data/feature-table.biom exported_data/${dataset}.feature-table.biom

biom convert \
  -i exported_data/${dataset}.feature-table.biom \
  -o exported_data/${dataset}.feature-table.txt \
  --to-tsv

