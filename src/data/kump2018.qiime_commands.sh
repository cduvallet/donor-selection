## Make sure you're in a QIIME 2 environment
## before running these!
# source activate qiime2-2018.6

# Import fastq files into QIIME 2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../fastq \
  --source-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path kump2018.demux-paired-end.qza

## Optional: visualize data quality
#qiime tools view kump2018.demux.qzv

# Quality filter reads
qiime quality-filter q-score \
  --i-demux kump2018.demux-paired-end.qza \
  --o-filtered-sequences kump2018.demux-quality-filtered.qza \
  --o-filter-stats kump2018.demux-quality-filtered-stats.qza

# Visualize quality filtered data
qiime metadata tabulate \
  --m-input-file kump2018.demux-quality-filtered-stats.qza \
  --o-visualization kump2018.demux-quality-filtered-stats.qzv

# Use deblur to denoise sequences
qiime deblur denoise-16S \
  --i-demultiplexed-seqs kump2018.demux-quality-filtered.qza \
  --p-trim-length 210 \
  --o-representative-sequences kump2018.rep-seqs-deblur.qza \
  --o-table kump2018.table-deblur.qza \
  --p-sample-stats \
  --o-stats kump2018.deblur-stats.qza

# Visualize deblur output stats
qiime deblur visualize-stats \
  --i-deblur-stats kump2018.deblur-stats.qza \
  --o-visualization kump2018.deblur-stats.qzv

# Assign taxonomy

# First, you'll need to prep the classifier
# See prep_qiime_classifier.sh script for eventual makefile

# TODO: change path to classifier to something smarter/better
qiime feature-classifier classify-sklearn \
  --i-reads kump2018.rep-seqs-deblur.qza \
  --i-classifier ../../silva-132-99-nb-classifier.qza \
  --o-classification kump2018.taxonomy.qza
