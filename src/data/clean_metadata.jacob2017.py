#!/usr/bin/env python
"""
Clean the Jacob 2017 metadata
"""

import pandas as pd

# This is copied from 2018-09-05.jacob2017.metadata.ipynb

# This file was acquired via email communication, and should not be public
foutcomes = 'data/raw/jacob2017/FMT Clinical Response Remission.xlsx'
# This file was downloaded from the ENA.
fmeta = 'data/raw/jacob2017/jacob2017.PRJNA388210.txt'
fout = 'data/clean/jacob2017.metadata.txt'

outcomes = pd.read_excel(foutcomes)
meta = pd.read_csv(fmeta, sep='\t')

# Expand the sample ID to different columns
sampledf = meta['sample_alias'].str.split('.', expand=True)

# Replace 'Donor' with 'donor'
sampledf = sampledf.replace('Donor', 'donor')

# Combine back with full metadata
sampledf = sampledf.rename(
    columns={0: 'study_id',
             1: 'FMT',
             2:'patient_id',
             3: 'sample_type'})
allmeta = pd.merge(
    sampledf[['patient_id', 'sample_type']], meta,
    left_index=True, right_index=True)

# Change the patient ID to an integer in allmeta, so that it
# matches the outcome data (which is an Excel file, so was
# automatically converted to an int)
allmeta['patient_id'] = allmeta['patient_id'].astype(int)
outcomes = outcomes.rename(
    columns={'Remission W4 ': 'remission_w4',
             'Response W4 ': 'response_w4'})
allmeta = pd.merge(allmeta, outcomes, left_on='patient_id', right_on='Patient')

# Add sample type indicator
allmeta['donor_patient'] = (
    allmeta['sample_type']
        .apply(lambda x: 'donor' if x == 'donor' else 'patient'))

# Reorder columns for my own sanity
first_cols = ['run_accession', 'patient_id',
              'sample_type', 'donor_patient',
              'remission_w4', 'response_w4',
              'sample_title']
col_order = first_cols + [i for i in allmeta.columns if i not in first_cols]
allmeta = allmeta[col_order]

# Write to file
allmeta.to_csv(fout, sep='\t', index=False)
