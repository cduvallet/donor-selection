## Makefile for donor selection paper

## Some Makefile notes
# Automatic variables: https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html#Automatic-Variables
# $@ is the target file
# $* is the stem that files have in common
# $< is the first file listed in the dependencies (usually a .py or .ipynb file)
# Multiple targets for one rule: https://www.gnu.org/software/automake/manual/html_node/Multiple-Outputs.html


all: data figures

## Process data
data: tidy_ibd metadata_ibd bn10


################################################
############### IBD FMT DATASETS ###############
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
############ IBD BUTYRATE ANALYSIS #############
################################################

# These commands grab the butyrate producers from the IBD datasets
# Note: as of Nov 14, I haven't checked or debugged these!

but_jacob := data/analysis/butyrate_producers.jacob2017.txt
but_goyal := data/analysis/butyrate_producers.goyal2018.txt
but_kump := data/analysis/butyrate_producers.kump2018.txt

but_script := src/analysis/calculate_butyrate_producer_abundance.py

$(but_jacob): $(but_script) $(jacob_meta) $(jacob_tidy)
	python $< jacob2017

$(but_goyal): $(but_script) $(goyal_meta) $(goyal_tidy)
	python $< goyal2018

$(but_kump): $(but_script) $(kump_meta) $(kump_tidy)
	python $< kump2018

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
################ BN10 ANALYSIS #################
################################################

donor_ranks := data/analysis/bn10_metabolomics.donor_ranks.txt

$(donor_ranks): src/analysis/calculate_donor_ranks.py $(tidy_bn10)
	python $< $(tidy_bn10) $@

################################################
############# POWER SIM DATASETS  ##############
################################################

download_src := src/data/download_case_control.sh
clean_src := src/data/clean_case_control_datasets.py

D = cdi_schubert
cdi_raw_otu := data/raw/$(D)_results/RDP/$(D).otu_table.100.denovo.rdp_assigned
cdi_raw_meta := data/raw/$(D)_results/$(D).metadata.txt
cdi_clean_otu := data/clean/$(D).otu_table.genus.feather
cdi_clean_meta := data/clean/$(D).metadata.feather

# This also downloads the metadata file, just FYI
$(cdi_raw_otu): $(download_src)
	$< cdi_schubert
	touch $(cdi_raw_otu)

# This also cleans the metadata file, just FYI
$(cdi_clean_otu): $(clean_src) $(cdi_raw_otu)
	$< cdi_schubert

D = crc_baxter
crc_raw_otu := data/raw/$(D)_results/RDP/$(D).otu_table.100.denovo.rdp_assigned
crc_raw_meta := data/raw/$(D)_results/$(D).metadata.txt
crc_clean_otu := data/clean/$(D).otu_table.genus.feather
crc_clean_meta := data/clean/$(D).metadata.feather

$(crc_raw_otu): $(download_src)
	$< crc_baxter
	touch $(crc_raw_otu)

$(crc_clean_otu): $(clean_src) $(crc_raw_otu)
	$< crc_baxter

D = ob_goodrich
ob_raw_otu := data/raw/$(D)_results/RDP/$(D).otu_table.100.denovo.rdp_assigned
ob_raw_meta := data/raw/$(D)_results/$(D).metadata.txt
ob_clean_otu := data/clean/$(D).otu_table.genus.feather
ob_clean_meta := data/clean/$(D).metadata.feather

$(ob_raw_otu): $(download_src)
	$< ob_goodrich
	touch $(ob_raw_otu)

$(ob_clean_otu): $(clean_src) $(ob_raw_otu)
	$< ob_goodrich

power_sim_data: $(cdi_clean_otu) $(cdi_clean_meta) $(crc_clean_otu) $(crc_clean_meta) $(ob_clean_otu) $(ob_clean_meta)

################################################
############## POWER SIMULATION  ###############
################################################

nreps := 50
nsig := data/analysis/power_simulation.n_sig.$(nreps)_reps.txt
tophits := data/analysis/power_simulation.top_hits_sig.$(nreps)_reps.txt
powersim_src := src/analysis/power_simulation.py

# This script also makes the data/analysis/population_effects.dataset.txt
# files, which have the population-level signal-to-noise ratios,
# but I won't put these in the makefile for simplicity.
$(nsig): $(powersim_src) $(power_sim_data)
	python $< --nreps $(nreps) --fout-nsig $(nsig) --fout-tophits $(tophits)

$(tophits): $(nsig)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

# No, it doesn't look like the latest version of power_simulation.py makes the population-level effects. That's in calculate_snr.py, I think...
# python calculate_snr.py -- looks like it's all hard-coded in there

################################################
################### FIGURES ####################
################################################

figures: fig2 fig3 fig4

## IBD case study
# Figures
fig_donor_butyrate_abun := figures/final/fig2.donor_butyrate_abun.png
fig_butyrate_response := figures/final/fig2.butyrate_vs_response.png
# Notebooks
ibd_fig_notebook_final := src/figures/figures.ibd_donors_butyrate_producers.ipynb
ibd_fig_notebook_src := src/exploration/2018-10-26.figures.ibd_donors_butyrate_producers.ipynb

fig2: $(fig_donor_butyrate_abun) $(fig_butyrate_response)

$(ibd_fig_notebook_final): $(ibd_fig_notebook_src)
	cp $< $@

$(fig_donor_butyrate_abun): $(ibd_fig_notebook_final) $(but_jacob) $(but_goyal) $(but_kump)
	jupyter nbconvert --execute $<

# Same notebook makes both figures
$(fig_butyrate_response): $(fig_donor_butyrate_abun)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

## Liver cirrhosis case study (BN10 metabolomics data)
# Figures
fig_distribution_mtabs := figures/final/fig3.scfas_bile_acid_conversion_bn10_donors.distribution.png
fig_scfa_bile_acid := figures/final/fig3.scfas_bile_acid_conversion_bn10_donors.ranked.png
fig_donor_ranks := figures/final/fig3.scfas_bile_acid_conversion_bn10_donors.bile_acid_vs_scfa.png

# Notebooks
bn10_fig_notebook_src := src/exploration/2018-11-15.figures.liver_cirrhosis_bn10_metabolomics.ipynb
bn10_fig_notebook_final := src/figures/figures.liver_cirrhosis_bn10_metabolomics_donors.ipynb

fig3: $(fig_distribution_mtabs) $(fig_scfa_bile_acid) $(fig_donor_ranks)

$(bn10_fig_notebook_final): $(bn10_fig_notebook_src)
	cp $< $@

$(fig_distribution_mtabs): $(bn10_fig_notebook_final) $(tidy_bn10) $(donor_ranks)
	jupyter nbconvert --execute $<

# Same notebook makes all panels
$(fig_scfa_bile_acid): $(fig_distribution_mtabs)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

# Same notebook makes all panels
$(fig_donor_ranks): $(fig_distribution_mtabs)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi


### Power simulation figure
fig_power_sim := figures/final/fig4.power_simulation.png
power_sim_fig_notebook_src := src/exploration/2018-11-19.figures.power_simulation.ipynb
power_sim_fig_notebook_final := src/figures/figures.power_simulation.ipynb
fig4: $(fig_power_sim)

$(power_sim_fig_notebook_final): $(power_sim_fig_notebook_src)
	cp $< $@

$(fig_power_sim): $(power_sim_fig_notebook_final) $(tophits)
	jupyter nbconvert --execute $<
