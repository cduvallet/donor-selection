#!/usr/bin/env python
"""
This script calculates the abundance of butyrate producers for each of the
three IBD datasets.
"""

import pandas as pd


def get_genera():
    """
    Returns the list of genera we are considering butyrate producers.

    The taxa from the "[Colonic Butyrate-Producing Communities in Humans: an Overview Using Omics Data](https://msystems.asm.org/content/2/6/e00130-17)" paper are:

    - Odoribacter
        - Table S2: 5/5 genomes have Acetyl-CoA
    - E. ventriosum
        - this is part of Lachnospiraceae_incertae_sedis
        - 1/1 genomes have Acetyl-CoA
    - E. hallii
        - 3/3 E. hallii genomes have Acetyl-CoA
        - this is one of the Lachno incertae sedis species
    - E. rectale
        - 7/7 genomes have Acetyl-CoA
        - this is one of the Lachno incertae sedis species
    - B. crossotus
        - this one is not present in >70% of people, it was just added to Fig 1 bc it's high abundance in a few people
    - S. variabile
        - this was the only subdoligranulum genome (and does have Acetyl-CoA)
    - ClostridiumXIVa
        - 42% of the 75 genomes have Acetyl-CoA
        - Paper doesn't address this low percentage, and in fact includes Clostridium XIVa as part of the "global core community"
    - Coprococcus
        - 92% of the 13 genomes have Acetyl-CoA
    - Butyricicoccus
        - 3/3 genomes have Acetyl-CoA
    - Pseudoflavonifractor
        - 66% of 3 genomes have Acetyl-Coa
    - Flavonifractor
        - 10/10 genomes have Acetyl-CoA
    - Anaerostipes
        - 14/14 genomes have Acetyl-CoA
    - R. intestinalis, R. inulinivorans, and R. faecis
        - 12/12 Roseburia genomes have Acetyl-CoA
        - "Manual inspections of major, abundant  genera led us to resolve Roseburia and Lachnospiraceae incertae sedis at the species level, as sequences  of individual species displayed high phylogenetic distances for all pathway genes"
    - Oscillibacter
        - 83% of 6 genomes have Acetyl-CoA
    - F. prausnitzii
        - 7/7 Faecalibacterium genomes have Acetyl-CoA
        - "For genera that encompassed only one species, such as F. prausnitzii and S. variabile, the species name is displayed."

    We'll use the easy/cleanest ones (genus-level):

    Odoribacter
    ClostridiumXIVa
    Coprococcus
    Butyricicoccus
    Pseudoflavonifractor
    Flavonifractor
    Anaerostipes
    Roseburia
    Oscillibacter
    Faecalibacterium
    Subdoligranulum

    """
    genera = ['Odoribacter',
              'Clostridium_XIVa',
              'Coprococcus',
              'Butyricicoccus',
              'Pseudoflavonifractor',
              'Flavonifractor',
              'Anaerostipes',
              'Roseburia',
              'Oscillibacter',
              'Faecalibacterium',
              'Subdoligranulum'
             ]
    return genera

def get_butyrate_genera(tidydf):
    """
    Get the OTU IDs for butyrate-producing genera in a given dataset.

    Note that the list of genera is not an input, it is acquired from
    the function get_genera().

    Parameters
    ----------
    tidydf : pandas DataFrame
        tidy OTU table, with column labeled 'otu_id_gg' that has the
        GreenGenes OTU ID

    Returns
    -------
    but_gg : list
        list of OTU IDs which contain a butyrate-producing genus
    """
    genera = get_genera()
    alltaxa_gg = tidydf['otu_id_gg'].unique().tolist()

    but_gg = []
    for g in genera:
        but_gg += [i for i in alltaxa_gg if g in i]
    but_gg = list(set(but_gg))

    return but_gg

def calculate_butyrate_abundance(tidydf, but_gg, dataset):
    """
    Calculate the total abundance of butyrate producers for each dataset.

    This is wrapped in a function because the different datasets have
    different query calls (because they have different metadata)

    Parameters
    ----------
    tidydf : pandas DataFrame
        tidy OTU table merged with respective dataset's metadata (on sample ID)
    but_gg : list
        list of OTU IDs that are butyrate producers
    dataset : str

    Returns
    -------
    all_but : pandas DataFrame
        dataframe with the dataset's respective metadata columns plus
        a column called "butyrate_abun" with the total abundance of
        OTUs in but_gg list per sample (including all samples, not just
        donors)
    """
    # The groupby is mostly to keep track of other metadata of interest
    if dataset == "goyal2018":
        all_but = (
            tidydf
                .query('otu_id_gg == @but_gg')
                #.query('(sample_type == "D") | (sample_type == "D2")')
                .groupby(['sample_id', 'patient_id',
                          'sample_type', 'time_point',
                          'remission_m1', 'response_m1',
                          'remission_m6', 'response_m6'])
                .sum()
                ['rel_abun']
            ).reset_index()

    elif dataset == "jacob2017":
        all_but = (
            tidydf
                .query('otu_id_gg == @but_gg')
                #.query('donor_patient == "donor"')
                .groupby(['sample_id', 'sample_type',
                          'donor_patient', 'remission_w4',
                          'response_w4', 'patient_id'])
                .sum()
                ['rel_abun']
            ).reset_index()

    elif dataset == "kump2018":
        all_but = (
            tidydf
                .query('otu_id_gg == @but_gg')
                #.query('Sampletype == "Donorstool"')
                .groupby(['sample_id', 'DonorID', 'PatientID',
                          'Sampling_day', 'Sampletype', 'prepost',
                          'Response'])
                .sum()
                ['rel_abun']
            ).reset_index()
        # Kump 2018 has multiple days per donor, so also get the
        # average abundance across all donor samples
        avg_but = (
            all_but
                .query('Sampletype == "Donorstool"')
                .groupby(['DonorID', 'Response', 'PatientID'])
                .mean()
                ['rel_abun']
            ).reset_index()
        avg_but['Sampling_day'] = 'average'
        all_but = pd.concat((all_but, avg_but), sort=False)

    else:
        raise ValueError('Unrecognized dataset')

    all_but = all_but.rename(columns={'rel_abun': 'butyrate_abun'})

    return all_but

###

datasets = ['jacob2017', 'kump2018', 'goyal2018']

for d in datasets:
    print(d)
    fmeta = 'data/clean/{}.metadata.txt'.format(d)
    ftidy = 'data/clean/{}.tidy_otu_w_taxonomy.txt'.format(d)

    # Define out file
    fbutyrate = 'data/analysis/butyrate_producers.{}.txt'.format(d)

    # Read in metadata
    meta = pd.read_csv(fmeta, sep='\t', index_col=0)

    # Read in tidy OTU table
    tidydf = pd.read_csv(ftidy, sep='\t')
    # Add relative abundance column
    tidydf['rel_abun'] = tidydf['reads'] / tidydf['total_reads']

    # Combine tidy OTU table with metadata
    if d == "goyal2018":
        tidydf = pd.merge(tidydf, meta, left_on='sample_id',
                          right_on='sample_alias')
    else:
        tidydf = pd.merge(tidydf, meta, left_on='sample_id',
                          right_index=True)

    # Get butyrate-producing OTUs
    but_gg = get_butyrate_genera(tidydf)

    # Calculate total abundance of butyrate producers per donor stool
    donor_but = calculate_butyrate_abundance(tidydf, but_gg, d)

    # Write file
    donor_but.to_csv(fbutyrate, sep='\t', index=False)
