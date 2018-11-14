## Makefile for donor selection paper

all: data

## Process data
data: tidy_ibd metadata_ibd bn10


################################################
################### IBD DATA ###################
################################################


# Need to figure out/whether to include the qiime2 steps in the Makefile.
# Currently, I'm keeping things separate per dataset. For the most
# part, I'm trying to make datasetID.qiime_commands.sh that have all
# the QIIME 2 commands.


## TO DO: Figure out how to parallelize this code. For now, I need results
## so I'll just define and make everything manually. Womp.

# DATASETS = kump2018 goyal2018 jacob2017
#
# # Define input files (QIIME 2 processing outputs)
# QIIME_STUB = $(foreach D, $(DATASETS), data/qiime-proc/$(D)/exported_data/$(D))
# TABLES = $(addsuffix feature-table.txt, $(QIIME_STUB))
# SILVA = $(addsuffix .taxonomy.silva-132-99.tsv, $(QIIME_STUB))
# GG = $(addsuffix .taxonomy.gg-13-8-99.tsv, $(QIIME_STUB))
#
# # Define output clean data files
# TIDY_TABLES = $(foreach D, $(DATASETS), data/clean/$(D).tidy_otu_w_taxonomy.txt)

## Add taxonomy and tidy OTU table
D = kump2018
kump_table := data/qiime-proc/$(D)/exported_data/$(D).feature-table.txt
kump_silva := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.silva-132-99.tsv
kump_gg := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.gg-13-8-99.tsv
kump_tidy := data/clean/$(D).tidy_otu_w_taxonomy.txt

$(kump_tidy): src/data/add_taxonomy_and_tidy_table.py # To add: OTU table and taxonomy files
	python $< --otu-table $(kump_table) \
		      --silva $(kump_silva) \
		      --gg $(kump_gg) \
		      --tidy-file $(kump_tidy)

D = jacob2017
jacob_table := data/qiime-proc/$(D)/exported_data/$(D).feature-table.txt
jacob_silva := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.silva-132-99.tsv
jacob_gg := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.gg-13-8-99.tsv
jacob_tidy := data/clean/$(D).tidy_otu_w_taxonomy.txt

$(jacob_tidy): src/data/add_taxonomy_and_tidy_table.py # To add: OTU table and taxonomy files
	python $< --otu-table $(jacob_table) \
		      --silva $(jacob_silva) \
		      --gg $(jacob_gg) \
		      --tidy-file $(jacob_tidy)

D = goyal2018
goyal_table := data/qiime-proc/$(D)/exported_data/$(D).feature-table.txt
goyal_silva := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.silva-132-99.tsv
goyal_gg := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.gg-13-8-99.tsv
goyal_tidy := data/clean/$(D).tidy_otu_w_taxonomy.txt

$(goyal_tidy): src/data/add_taxonomy_and_tidy_table.py # To add: OTU table and taxonomy files
	python $< --otu-table $(goyal_table) \
		      --silva $(goyal_silva) \
		      --gg $(goyal_gg) \
		      --tidy-file $(goyal_tidy)


tidy_ibd: $(kump_tidy) $(jacob_tidy) $(goyal_tidy)

## Metadata files

D = kump2018
kump_meta := data/clean/$(D).metadata.txt
raw_kump_meta := data/raw/kump2018/mapping_file.txt

$(kump_meta): src/data/clean_metadata.kump2018.py $(raw_kump_meta)
	python $<

D = goyal2018
goyal_meta := data/clean/$(D).metadata.txt
raw_goyal_meta := data/raw/goyal2018/FMT_study_log_23Sept2016_edited.xlsx \
			data/raw/goyal2018/goyal2018.PRJNA380944.txt

$(goyal_meta): src/data/clean_metadata.goyal2018.py $(raw_goyal_meta)
	python $<

D = jacob2017
jacob_meta := data/clean/$(D).metadata.txt
raw_jacob_meta := data/raw/jacob2017/FMT\ Clinical\ Response\ Remission.xlsx \
			data/raw/jacob2017/jacob2017.PRJNA388210.txt

$(jacob_meta): src/data/clean_metadata.jacob2017.py $(raw_jacob_meta)
	python $<

metadata_ibd: $(kump_meta) $(jacob_meta) $(goyal_meta)

################################################
################## BN 10 DATA ##################
################################################

raw_bn10 := data/raw/bn10/BN10_newIDs.csv
tidy_bn10 := data/clean/bn10.tidy_metabolomics.feather

$(raw_bn10):
	echo -e "Need to download raw BN10 data from internet"

$(tidy_bn10): src/data/bn10.tidy_data.py $(raw_bn10)
	python $< $(raw_bn10) $@

bn10: $(tidy_bn10)

################################################
################### FIGURES ####################
################################################
