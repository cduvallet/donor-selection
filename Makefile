## Makefile for donor selection paper

all: data

## Process data
data: tidy

# Need to figure out/whether to include these steps in the Makefile.
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


tidy: $(kump_tidy) $(jacob_tidy) $(goyal_tidy)
