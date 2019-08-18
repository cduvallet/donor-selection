## Makefile for donor selection paper

## Some Makefile notes
# Automatic variables: https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html#Automatic-Variables
# $@ is the target file
# $* is the stem that files have in common
# $< is the first file listed in the dependencies (usually a .py or .ipynb file)
# Multiple targets for one rule: https://www.gnu.org/software/automake/manual/html_node/Multiple-Outputs.html

# If you've never seen a Makefile before, the way that you read it is:
# file_that_you_want_to_make: files you need to make it
#    thing that you run to make the file from the dependency files

# i.e. it can look like:
# results.txt: script.py data.csv
#    python script.py --input data.csv --output results.txt

all: data figures supp_figs

## Process data
data: tidy_ibd metadata_ibd bn10 power_sim_data power_sim_res

################################################
############### IBD FMT DATASETS ###############
################################################

# Note: downloading and processing these datasets is
# unforunately not fully reproducible with make. I will
# write the steps I took here, but don't have time to
# make it fully reproducible via make.

######### DOWNLOAD RAW DATA #########

### Jacob 2017
## Notes available in data/raw/jacob2017/jacob2017.readme.txt
## Download data via FTP (URLs gotten from ENA)
## From the data/raw/jacob2017 directory:
# wget -P data/ -i jacob2017.fastq_ftp.txt

### Kump 2018
## Notes available in data/raw/kump2018/kump2018.readme.txt
## Download data via FTP (URLs gotten from ENA)
## From the data/raw/kump2018 directory:
# wget -P data/ -i kump2018.fastq_ftp.txt

### Goyal 2018
## No readme, but I think I did the same steps as
## for the others... :(
## Download data via FTP (URLs gotten from ENA)
## From the data/raw/goyal2018 directory:
# wget -P data/ -i goyal2018.fastq_ftp.txt

######### PREP METADATA #########

# Note: there are also recipes for this process below,
# under the "add taxonomy section".

### Jacob 2017
# Run src/data/jacob2017.make_metadata.sh
# (This just calls the iPython notebook that does it)

### Kump 2018
# From the data/raw/kump2018/ directory:
# wget https://raw.githubusercontent.com/LGPW/FMT_in_TRUC_APT2017/master/mapping_file.txt

### Goyal 2018
# # Define raw file names
#data/raw/goyal2018/FMT_study_log_23Sept2016_edited.xlsx : acquired via email communication
#data/raw/goyal2018/goyal2018.PRJNA380944.txt: downloaded from ENA
# src/data/goyal2018.make_metadata.sh (which just calls the iPython notebook that the .py file contains a copy of)
# creates data/clean/goyal2018.metadata.txt

# All of these datasets also have a clean_metadata.{dataset}.py file,
# which has these commands in them (and/or a copy of the respective iPython nb)

############ PROCESS WITH QIIME 2 ############

## For all datasets, I tracked all of the qiime2 commands I ran in
## src/data/{dataset}.qiime_commands.sh

## Importing data into qiime 2
# Kump 2018: import with format CasavaOneEightSingleLanePerSampleDirFmt
# Goyal 2018: make manifest file with code pasted in goyal2018.processing_notes.txt
# Jacob 2017: looks like I have a jacob2017.manifest.csv in data/qiime-proc/jacob2017, but no code to make it. I probably made it via a quick python session, like in the Goyal study.

## The only other relevant file in src/data is prep_qiime_classifier.sh,
## which downloads the taxonomy classifier from the QIIME 2 website

###### PROCESS AND CLEAN DATA #######

## Add taxonomy and tidy OTU table
# Define file names
D = kump2018
# Inputs
kump_table := data/qiime-proc/$(D)/exported_data/$(D).feature-table.txt
kump_silva := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.silva-132-99.tsv
kump_gg := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.gg-13-8-99.tsv
# Output
kump_tidy := data/clean/$(D).tidy_otu_w_taxonomy.txt

# Run script to tidyfy and add latin names to OTU table
$(kump_tidy): src/data/add_taxonomy_and_tidy_table.py $(kump_table) $(kump_gg) $(kump_silva)
	python $< --otu-table $(kump_table) \
		      --silva $(kump_silva) \
		      --gg $(kump_gg) \
		      --tidy-file $(kump_tidy)

# Define file names
D = jacob2017
# Inputs
jacob_table := data/qiime-proc/$(D)/exported_data/$(D).feature-table.txt
jacob_silva := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.silva-132-99.tsv
jacob_gg := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.gg-13-8-99.tsv
# Output
jacob_tidy := data/clean/$(D).tidy_otu_w_taxonomy.txt

# Run script to tidyfy and add latin names to OTU table
$(jacob_tidy): src/data/add_taxonomy_and_tidy_table.py  $(jacob_table) $(jacob_silva) $(jacob_gg)
	python $< --otu-table $(jacob_table) \
		      --silva $(jacob_silva) \
		      --gg $(jacob_gg) \
		      --tidy-file $(jacob_tidy)

# Define file names
D = goyal2018
# Inputs
goyal_table := data/qiime-proc/$(D)/exported_data/$(D).feature-table.txt
goyal_silva := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.silva-132-99.tsv
goyal_gg := data/qiime-proc/$(D)/exported_data/$(D).taxonomy.gg-13-8-99.tsv
# Outputs
goyal_tidy := data/clean/$(D).tidy_otu_w_taxonomy.txt

# Run script to tidyfy and add latin names to OTU table
$(goyal_tidy): src/data/add_taxonomy_and_tidy_table.py $(goyal_table) $(goyal_silva) $(goyal_gg)
	python $< --otu-table $(goyal_table) \
		      --silva $(goyal_silva) \
		      --gg $(goyal_gg) \
		      --tidy-file $(goyal_tidy)

# Gather all tidy OTU tables as one target
tidy_ibd: $(kump_tidy) $(jacob_tidy) $(goyal_tidy)

## Make associated clean metadata files
# This is also described above, in the "Prep metadata" section

D = kump2018
# Output
kump_meta := data/clean/$(D).metadata.txt
# Input
raw_kump_meta := data/raw/kump2018/mapping_file.txt

# Script to clean metadata (after you've downloaded it, see section above)
$(kump_meta): src/data/clean_metadata.kump2018.py $(raw_kump_meta)
	python $<

D = goyal2018
#Output
goyal_meta := data/clean/$(D).metadata.txt
# Inputs, acquired as described in the section above
raw_goyal_meta := data/raw/goyal2018/FMT_study_log_23Sept2016_edited.xlsx \
			data/raw/goyal2018/goyal2018.PRJNA380944.txt

# Script to clean and combine inputs
$(goyal_meta): src/data/clean_metadata.goyal2018.py $(raw_goyal_meta)
	python $<

D = jacob2017
# Output
jacob_meta := data/clean/$(D).metadata.txt
# Inputs, acquired as described in the section above
raw_jacob_meta := data/raw/jacob2017/FMT\ Clinical\ Response\ Remission.xlsx \
			data/raw/jacob2017/jacob2017.PRJNA388210.txt

# Script to clean and combine inputs
$(jacob_meta): src/data/clean_metadata.jacob2017.py $(raw_jacob_meta)
	python $<

# Gather all clean metadata files as one target
metadata_ibd: $(kump_meta) $(jacob_meta) $(goyal_meta)

################################################
############ IBD BUTYRATE ANALYSIS #############
################################################

# These commands grab the butyrate producers from the IBD datasets

# Define output files
but_jacob := data/analysis/butyrate_producers.jacob2017.txt
but_goyal := data/analysis/butyrate_producers.goyal2018.txt
but_kump := data/analysis/butyrate_producers.kump2018.txt

# Script that does the butyrate abundance calculation
but_script := src/analysis/calculate_butyrate_producer_abundance.py

# Run for all datasets, with the metadata and tidy OTU tables as inputs
$(but_jacob): $(but_script) $(jacob_meta) $(jacob_tidy)
	python $< jacob2017

$(but_goyal): $(but_script) $(goyal_meta) $(goyal_tidy)
	python $< goyal2018

$(but_kump): $(but_script) $(kump_meta) $(kump_tidy)
	python $< kump2018

################################################
################## BN 10 DATA ##################
################################################

# Note: BN10 referes to the Alm lab project that this
# metabolomics data was generated from. So anything labeled
# BN10 refers to the metabolomics / liver cirrhosis case study.

# Raw data file acquired from communication with the authors
# Note: the version they upload may be slightly different than
# the version they gave me, so apologies in advance if this code doesn't work

# Input, acquired from authors
raw_bn10 := data/raw/bn10/BN10_newIDs.csv
# Output
tidy_bn10 := data/clean/bn10.tidy_metabolomics.feather

# Just a placeholder recipe that tells you that you need to get the raw data yourself
$(raw_bn10):
	echo -e "Need to download raw BN10 data from the internet or the authors"

# Run script to tidy the raw data
$(tidy_bn10): src/data/bn10.tidy_data.py $(raw_bn10)
	python $< $(raw_bn10) $@

bn10: $(tidy_bn10)

################################################
################ BN10 ANALYSIS #################
################################################

# Output file
donor_ranks := data/analysis/bn10_metabolomics.donor_ranks.txt

# Make output ranks file from the tidyfied metabolomics data
$(donor_ranks): src/analysis/calculate_donor_ranks.py $(tidy_bn10)
	python $< $(tidy_bn10) $@

################################################
############# POWER SIM DATASETS  ##############
################################################

# Recipes to download and clean the case-control datasets
# used in the power simulation

# Input scripts
download_src := src/data/download_case_control.sh
clean_src := src/data/clean_case_control_datasets.py

# Define file names
D = cdi_schubert
# Files downloaded from Zenodo
cdi_raw_otu := data/raw/$(D)_results/RDP/$(D).otu_table.100.denovo.rdp_assigned
cdi_raw_meta := data/raw/$(D)_results/$(D).metadata.txt
# Clean output files
cdi_clean_otu := data/clean/$(D).otu_table.genus.feather
cdi_clean_meta := data/clean/$(D).metadata.feather

# Download the data
# This also downloads the metadata file, just FYI
$(cdi_raw_otu): $(download_src)
	$< cdi_schubert
	touch $(cdi_raw_otu)

# Clean up the data
# This also cleans the metadata file, just FYI
$(cdi_clean_otu): $(clean_src) $(cdi_raw_otu)
	$< cdi_schubert

# Define file names
D = crc_baxter
# Files downloaded from Zenodo
crc_raw_otu := data/raw/$(D)_results/RDP/$(D).otu_table.100.denovo.rdp_assigned
crc_raw_meta := data/raw/$(D)_results/$(D).metadata.txt
# Clean output files
crc_clean_otu := data/clean/$(D).otu_table.genus.feather
crc_clean_meta := data/clean/$(D).metadata.feather

# Download the OTU table and metadata
$(crc_raw_otu): $(download_src)
	$< crc_baxter
	touch $(crc_raw_otu)

# Clean the OTU table and metadata
$(crc_clean_otu): $(clean_src) $(crc_raw_otu)
	$< crc_baxter

# Define file names
D = ob_goodrich
# Files downloaded from Zenodo
ob_raw_otu := data/raw/$(D)_results/RDP/$(D).otu_table.100.denovo.rdp_assigned
ob_raw_meta := data/raw/$(D)_results/$(D).metadata.txt
# Clean output files
ob_clean_otu := data/clean/$(D).otu_table.genus.feather
ob_clean_meta := data/clean/$(D).metadata.feather

# Download the OTU table and metadata
$(ob_raw_otu): $(download_src)
	$< ob_goodrich
	touch $(ob_raw_otu)

# Clean the OTU table and metadata
$(ob_clean_otu): $(clean_src) $(ob_raw_otu)
	$< ob_goodrich

# Gather all the data into one target
power_sim_data: $(cdi_clean_otu) $(cdi_clean_meta) $(crc_clean_otu) $(crc_clean_meta) $(ob_clean_otu) $(ob_clean_meta)

################################################
############## POWER SIMULATION  ###############
################################################

# Do the actual power simulation

# Set the number of reps to do (one rep = one subsetting of the samples, for
# each of the responder proportions)
nreps := 50
# Define output files
# nsig has the number of significant bugs total
nsig := data/analysis/power_simulation.n_sig.$(nreps)_reps.txt
# tophits has the number of significant top hits
tophits := data/analysis/power_simulation.top_hits_sig.$(nreps)_reps.txt
# Script that performs the power simulation
powersim_src := src/analysis/power_simulation.py

## Make the data/analysis/population_effects.dataset.txt, which define
## the "true" differential bugs based on signal-to-noise
# The calculate_snr.py script has all the datasets hard-coded in it,
# and automatically makes the other population effects files.
pop_effects_script := src/analysis/calculate_snr.py

# Only do this for one dataset, since all the other datasets are
# processed at the same time
D = ob_goodrich
pop_goodrich := data/analysis/population_effects.{}.txt
$(pop_goodrich): $(pop_effects_script) $(ob_clean_otu) $(ob_clean_meta)
	python $<

## Do the power simulation
# This script reads in the data/analysis/population_effects.{dataset}.txt
# files for all datasets (which are hard-coded).
# It runs the power simulation and returns the number of significant
# bugs (in the *.n_sig.* file) and the number of top hits which were
# significant (*.top_hits_sig.*). The *.top_hits_sig.* is the file
# used to make the final figure in the paper.
$(nsig): $(powersim_src) $(power_sim_data)
	python $< --nreps $(nreps) --fout-nsig $(nsig) --fout-tophits $(tophits)

# This is just a make recipe that checks for both files which are
# made by the script
$(tophits): $(nsig)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

power_sim_res: $(tophits)

################################################
################### FIGURES ####################
################################################

# Make the main figures

figures: fig2 fig3 fig4

### IBD case study
## Figures
# Figure 2A
fig_donor_butyrate_abun := figures/final/fig2.donor_butyrate_abun.png
# Figure 2B
fig_butyrate_response := figures/final/fig2.butyrate_vs_response.png

## Notebooks
# Final src/figures notebook
ibd_fig_notebook_final := src/figures/figures.ibd_donors_butyrate_producers.ipynb
# Input notebook in src/exploration
ibd_fig_notebook_src := src/exploration/2018-10-26.figures.ibd_donors_butyrate_producers.ipynb

fig2: $(fig_donor_butyrate_abun) $(fig_butyrate_response)

# Copy notebook in src/exploration to src/figures
$(ibd_fig_notebook_final): $(ibd_fig_notebook_src)
	cp $< $@

# Execute the notebook, which makes an html and saves the figures
$(fig_donor_butyrate_abun): $(ibd_fig_notebook_final) $(but_jacob) $(but_goyal) $(but_kump)
	jupyter nbconvert --execute $<

# Same notebook makes both figures
$(fig_butyrate_response): $(fig_donor_butyrate_abun)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

### Liver cirrhosis case study (BN10 metabolomics data)
## Figures
fig_distribution_mtabs := figures/final/fig3.scfas_bile_acid_conversion_bn10_donors.distribution.png
fig_scfa_bile_acid := figures/final/fig3.scfas_bile_acid_conversion_bn10_donors.ranked.png
fig_donor_ranks := figures/final/fig3.scfas_bile_acid_conversion_bn10_donors.bile_acid_vs_scfa.png

# Notebooks
# Input notebook in src/exploration
bn10_fig_notebook_src := src/exploration/2018-11-15.figures.liver_cirrhosis_bn10_metabolomics.ipynb
# Final src/figures notebook
bn10_fig_notebook_final := src/figures/figures.liver_cirrhosis_bn10_metabolomics_donors.ipynb

fig3: $(fig_distribution_mtabs) $(fig_scfa_bile_acid) $(fig_donor_ranks)

# Copy notebook in src/exploration to src/figures
$(bn10_fig_notebook_final): $(bn10_fig_notebook_src)
	cp $< $@

# Execute the notebook, which makes an html and saves the figures
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
# Input notebook in src/exploration
power_sim_fig_notebook_src := src/exploration/2018-11-19.figures.power_simulation.ipynb
# Final src/figures notebook
power_sim_fig_notebook_final := src/figures/figures.power_simulation.ipynb
fig4: $(fig_power_sim)

# Copy notebook in src/exploration to src/figures
$(power_sim_fig_notebook_final): $(power_sim_fig_notebook_src)
	cp $< $@

# Execute the notebook, which makes an html and saves the figures
$(fig_power_sim): $(power_sim_fig_notebook_final) $(tophits)
	jupyter nbconvert --execute $<

## Supplementary figures

supp_figs: supp13 supp4

# Supp Fig 1-3
fig_donor_genera_goyal := figures/final/suppfig.butyrate_producer_genera.goyal.png
fig_donor_genera_kump := figures/final/suppfig.butyrate_producer_genera.kump.png
fig_donor_genera_jacob := figures/final/suppfig.butyrate_producer_genera.jacob.png

## Notebooks
# Final src/figures notebook
genera_suppfig_notebook_final := src/figures/suppfigures.ibd_donors_butyrate_all_genera.ipynb
# Input notebook in src/exploration
genera_suppfig_notebook_src := src/exploration/2018-11-30.detailed_butyrate_producer_abundances.ipynb

supp13: $(fig_donor_genera_goyal) $(fig_donor_genera_kump) $(fig_donor_genera_jacob)

# Copy notebook in src/exploration to src/figures
$(genera_suppfig_notebook_final): $(genera_suppfig_notebook_src)
	cp $< $@

# Execute the notebook, which makes an html and saves the figures
$(fig_donor_genera_goyal): $(genera_suppfig_notebook_final) $(but_jacob) $(but_goyal) $(but_kump)
	jupyter nbconvert --execute $<

# Same notebook makes all three Supp Figs
$(fig_donor_genera_kump): $(fig_donor_genera_goyal)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

# Same notebook makes all three Supp Figs
$(fig_donor_genera_jacob): $(fig_donor_genera_goyal)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

## Fig S4: donor-patient delta butyrate abundance vs. response
fig_delta_butyrate := figures/final/suppfig.delta_butyrate_vs_response.png

## Notebooks
# Final src/figures notebook
delta_but_suppfig_notebook_final := src/figures/suppfigure.ibd_butyrate_match_patient_donors.ipynb
# Input notebook in src/exploration
delta_but_suppfig_notebook_src := src/exploration/2018-10-29.figures.ibd_butyrate_match_patient_donors.ipynb

supp4: $(fig_delta_butyrate)

# Copy notebook in src/exploration to src/figures
$(delta_but_suppfig_notebook_final): $(delta_but_suppfig_notebook_src)
	cp $< $@

# Execute the notebook, which makes an html and saves the figures
$(fig_delta_butyrate): $(delta_but_suppfig_notebook_final) $(but_jacob) $(but_goyal) $(but_kump)
	jupyter nbconvert --execute $<
