#!/usr/bin/env bash

# Script to download the case-control studies from MicrobiomeHD
# to do the "power/sample size" simulations

# ${1} is the name of the target tar file, in format:
# data/raw_otu_tables/nash_chan_results.tar.gz
## Usage: download_tar_folders.sh data/raw_otu_tables/target_folder.tar.gz

# If ${1} is data/raw_otu_tables/nash_chan_results.tar.gz, then orig is
# nash_chan_results.tar.gz

download_dataset() {

    datadir=data/raw

    file=${1}_results.tar.gz
    target_path=${datadir}/${file}

    ## Download *.tar.gz from Zenodo to data/raw/*.tar.gz
    # --no-check-certificate in case you're on a secure connection, like MIT Secure...
    wget -O ${target_path} https://zenodo.org/record/1146764/files/${file} --no-check-certificate

    # Need to update timestamp on file, otherwise it keeps the timestamp from
    # upload day to Zenodo
    touch ${target_path}

    ## Extract file, in the data/raw/ directory as well
    targetdir=${target_path%/*}
    tar -C ${targetdir} -xvf ${target_path}

}

download_dataset $1
