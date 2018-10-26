#!/usr/bin/env python
"""
Clean up Goyal 2018 metadata
"""

import pandas as pd
import numpy as np

# This is copied from the 2018-10-01.goyal2018.metadata.ipynb notebook

# Define raw file names
fexcel = 'data/raw/goyal2018/FMT_study_log_23Sept2016_edited.xlsx'
fprjna = 'data/raw/goyal2018/goyal2018.PRJNA380944.txt'

# Out file
fout = 'data/clean/goyal2018.metadata.txt'

# Read in metadata files
outcomes = pd.read_excel(fexcel, sheet_name="PUCAI-PCDAI", skiprows=2)
outcomes = outcomes.rename(columns={'Unnamed: 5': 'Notes'})
meta = pd.read_csv(fprjna, sep='\t')

## Clean up Sample IDs
# Sample IDs contain metadata info, in the following order
# (from Andrew's email)
smpls = meta['sample_alias'].str.split('.', expand=True)
smpls = smpls.rename(columns={0: 'FMT',
                              1: 'sequencing_run',
                              2: 'patient_id',
                              3: 'sample_type'})
meta = pd.merge(meta, smpls, left_index=True, right_index=True)

## Parse sample types
# Let's parse the sample types. From Andrew's email:
#
# > D refers to the donor sample. When there is a 1 or 2 that means we sent the donor sample for repeat sequencing (most likely because it failed). P is the patient sample before FMT. 4M, 6M, 9M refer to the 4 month, 6 month, and 9 month sample. W is the 1 week sample. WX, P1, P2, PA, PB, 6MX are just samples that underwent repeat sequencing because they failed.
#
# Not sure what just `M` is, but the paper includes a 1 month sample so let's guess that it's that. (Also note: the paper doesn't actually look at 4 month or 9 month samples, I wonder why...)

time_dict = {'W': '1_week',
             'P': 'pre_fmt',
             'D': 'donor',
             'M': '1_month',
             '4M': '4_month',
             '6M': '6_month',
             'PB': 'pre_fmt',
             'D1': 'donor',
             'D2': 'donor',
             'PX': 'pre_fmt',
             '6MX': '6_month',
             'P1': 'pre_fmt',
             'WX': '1_week',
             'P2': 'pre_fmt',
             'PA': 'pre_fmt',
             '9M': '9_month'}

meta['time_point'] = meta['sample_type'].apply(lambda x: time_dict[x])

## Clean up clinical metadata
outcomes['patient_id'] = (
    outcomes['Patient # dds 1/26/15'].str.split(expand=True)[0])

# Add responder/remission flags, based on the clinical outcomes.
#
# The paper says:
#
# > Response was defined as a decrease of 15 points in PUCAI or 12.5 points in PCDAI at 1 month, as used in previous studies. Remission was defined as normalization of previously elevated fecal biomarkers and a PCDAI/PUCAI of 0 points. If subjects required escalation of medical therapy prior to 1-month evaluation, they were considered to be nonresponders. Subsequently, any escalation of medical therapy was considered a loss of response.
#
# I don't really know how to figure out who required escalation of medical therapy or something, but I'm assuming that would be reflected in the PCDAI/PUCAI scores..? Anyway, for now let's just go from those.
#
# Note: I think that the numbers highlighted in red are the PCDAI scores and the black ones are PUCAI. There aren't very many red ones, so let's just go with everything being PUCAI...
outcomes.loc[23, 'Notes'] = outcomes.loc[23, 'Month 6']
outcomes.loc[23, 'Month 6'] = np.nan

# We'll define response as _either_ remission (value is 0) OR decrease in at least 15 points.
outcomes['remission_m1'] = outcomes['Month 1'] == 0
outcomes['response_m1'] = ((outcomes['Screen'] - outcomes['Month 1']) >= 15) | outcomes['remission_m1']

outcomes['remission_m6'] = outcomes['Month 6'] == 0
outcomes['response_m6'] = ((outcomes['Screen'] - outcomes['Month 6']) >= 15) | outcomes['remission_m6']

## Merge clinical and ENA metadata
fullmeta = pd.merge(meta, outcomes)
fullmeta.to_csv(fout, sep='\t', index=False)
