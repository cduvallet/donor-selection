#!/usr/bin/env python
# Code to clean up OTU tables and metadata
import pandas as pd
import janitor as jn

# Hsiao et al
otu_out = 'data/clean/hsiao.dada2_otu_table.txt'
meta_out = 'data/clean/hsiao.metadata.txt'

## OTU table
## Needs: remove first line (#converted from biom), transpose so that OTUs are in columns

######### OTU TABLES #######

## DADA 2 table
fotu = 'data/raw/hsiao/qiime/hsiao_export/hsiao.feature_table.txt'

# Read in, remove biom header, and transpose
df = pd.read_csv(fotu, sep='\t', skiprows=1, index_col=0).T

# Write to clean file
df.to_csv(otu_out, sep='\t')

######## METADATA ##########

## Clean file
fmeta = 'data/raw/hsiao/TableS3_human16S_tidy.xlsx'
meta = pd.read_excel(fmeta)

# Rename columns
meta = jn.clean_names(meta, strip_underscores='both')
meta = meta.rename(columns={'fecal_sampleid*': 'fecal_sampleid'})

# Grab only samples with OTU data
samples = df.index.tolist()
meta = meta.query('fecal_sampleid == @samples')

# Set index to be sample IDs
meta.index = meta['fecal_sampleid']

# And add "sample_id" column
meta['sample_id'] = meta['fecal_sampleid']

## Make useful columns

# Add disease label: group A children are coded in the sample ID,
# group G children are all healthy
a_group = meta['group_description'].str.startswith('a')
meta.loc[a_group, 'disease'] = (meta
    .loc[a_group, 'sample_id']
    .str.split('.')
    .str[1])
g_group = meta['group_description'].str.startswith('g')
meta.loc[g_group, 'disease'] = 'healthy'

# Add day of infection or post-infection (raw_day)
for subj, subdf in meta.loc[a_group].groupby('subject_id'):
    tmp = subdf['fecal_sampleid'].tolist()
    tmp_days = []
    day = 0
    last_diarrhea = 0
    for i in tmp:
        if i.split('.')[1] == "diarrhea":
            if i.split('.')[-1] != "end":
                day = float(i.split('.')[-1])
            else:
                day += 1
        elif i.split('.')[1] == "recovery":
            # Grab final day
            if last_diarrhea == 0:
                last_diarrhea = day
            day = last_diarrhea + float(i.split('.')[-1][1:])
        tmp_days.append(day)
    meta.loc[subdf.index, 'raw_day'] = tmp_days

# Add progression through infection, code the post-infection days as "post", and healthy as "healthy" (infection_day)
for subj, subdf in meta.loc[a_group].groupby('subject_id'):
    # Re-code diarrhea samples as percent (0-1) through infection
    endday = float(subdf.query('disease == "diarrhea"').max()['raw_day'])
    diarrhea_idx = subdf.query('disease == "diarrhea"').index
    meta.loc[diarrhea_idx, 'infection_progression'] = subdf.loc[diarrhea_idx, 'raw_day'] / endday

    # Also code recovery progression, from 1-2 (will allow diverging scales woo)
    lastrecover = float(
        subdf
        .query('disease == "recovery"')
        .max()['raw_day']
        ) - endday
    recovery_idx = subdf.query('disease == "recovery"').index
    meta.loc[recovery_idx, 'infection_progression'] = (
        subdf.loc[recovery_idx,'raw_day'] - endday) / lastrecover + 1.0

# And code healthy patients as healthy
meta.loc[g_group, 'infection_progression'] = 'healthy'

# Write to clean file
meta.to_csv(meta_out, sep='\t')
