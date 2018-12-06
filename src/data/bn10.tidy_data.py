#!/usr/bin/env python
"""
This data cleans up the BN10 metabolomics data and writes
it to a tidy feather file.
"""
import argparse

import pandas as pd
import feather

p = argparse.ArgumentParser()
p.add_argument('raw_data')
p.add_argument('clean_data')
args = p.parse_args()

# Raw data file
fname = args.raw_data

# Out file
ftidy = args.clean_data

# Read data
df = pd.read_csv(fname)
# Convert to tidy data, with metabolite (+info), sample, intensity in each row
tidydf = df.melt(id_vars=['Compound', 'mz', 'RT', 'Method',
                          'Compound.ID', 'HMDB.ID...representative.',
                          'Metabolite', 'Number.of.Members',
                          'Cluster.Name', 'Cluster.Number',
                          'Cluster.major.ion', 'potential.adduct'],
                var_name='sample',
                value_name='intensity')
# Rename columns
tidydf = tidydf.rename(
    columns={'Compound': 'compound',
             'RT': 'rt',
             'Method': 'method',
             'Compound.ID': 'compound_id',
             'HMDB.ID...representative.': 'hmdb_id',
             'Metabolite': 'metabolite',
             'Number.of.Members': 'n_members',
             'Cluster.Name': 'cluster_name',
             'Cluster.Number': 'cluster_number',
             'Cluster.major.ion': 'cluster_major_ion',
             'potential.adduct': 'potential_adduct'})

# Split sample ID into donor ID and sample number
tidydf[['donor', 'sample_number']] = (
    tidydf['sample'].str.split('.', expand=True)
    )

# Clean up lithocholate, which was measured in C18neg and HILIC neg
# All the other bile acids were in C18, so we'll use that one here too
# See 2018-11-14.bn10_metabolomics.tidy_julians_new_data.ipynb for more work
litho_hiln_idx = (
    tidydf
        .query('metabolite == "lithocholate"')
        .query('method == "HILn"')
        .index
    )
tidydf.loc[litho_hiln_idx, 'metabolite'] = 'lithocholate-HILn'

# Fill NaN's with zeros
tidydf['intensity'] = tidydf['intensity'].fillna(0.0)

# Write to file
feather.write_dataframe(tidydf, ftidy)
