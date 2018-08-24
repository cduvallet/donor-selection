# script to prepare to assign taxonomy with qiime2

## download pre-trained classifier from qiime2 data resources page
# https://docs.qiime2.org/2018.6/data-resources/
wget https://data.qiime2.org/2018.6/common/silva-132-99-nb-classifier.qza

## check md5sum
# user should replace this with the expected md5sum from the data resources page
expected="a02c3f7473fa4369bbc66158c799d39a (expected md5)"

# calculate the actual md5sum of the downloaded file
actual="$(md5sum silva-132-99-nb-classifier.qza)"

echo $expected
echo $actual
