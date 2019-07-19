#!/usr/bin/env python
"""
This script calculates the signal to noise ratio and differential
presence for each genus and OTU.
"""

import pandas as pd
import numpy as np
import feather

def read_dataframe(fn):
    df = feather.read_dataframe(fn)
    df.index = df.iloc[:,0]
    df = df.iloc[:, 1:]
    return df

datasets = ['cdi_schubert', 'crc_baxter', 'ibd_papa', 'ob_goodrich']

labels = {'cdi_schubert': {'control': 'H', 'case': 'CDI'},
          'crc_baxter': {'control': 'H', 'case': 'CRC'},
          'ibd_papa': {'control': 'nonIBD', 'case': 'IBD'},
          'ob_goodrich': {'control': 'H', 'case': 'OB'}
          }

for d in datasets:
    print(d)

    fmeta = 'data/clean/' + d + '.metadata.feather'
    meta = read_dataframe(fmeta)

    ctrl_lbl = labels[d]['control']
    case_lbl = labels[d]['case']

    for taxa in ['otu', 'genus']:
        print(taxa)
        if taxa == 'otu':
            fotu = 'data/clean/' + d + '.otu_table.feather'
        else:
            fotu = 'data/clean/' + d + '.otu_table.genus.feather'

        df = read_dataframe(fotu)

        if d == 'ob_goodrich':
            # Remove some samples from OB Goodrich to reduce size of data
            meta = meta.query('(DiseaseState == @ctrl_lbl) | (DiseaseState == @case_lbl)')
            df = df.loc[meta.index]

        df.index.name = 'sample_id'

        ## Prep data
        print('Preparing data...')
        # Convert OTU table to tidy
        tidydf = (
            df
            .reset_index()
            .melt(id_vars='sample_id', value_name='abun', var_name='otu')
            )
        # Add column with only the last part of the OTU string
        tidydf['denovo'] = tidydf['otu'].str.split(';', 7).str[-1]

        # Add metadata
        tidydf = pd.merge(
            tidydf, meta,
            left_on='sample_id', right_index=True)

        # Convert to log abundance
        min_abun = tidydf['abun'].replace(0, np.nan).min() / 2
        tidydf['log_abun'] = np.log10(tidydf['abun'] + min_abun)

        # Get the mean abundance per group
        meandf = tidydf.groupby(['denovo', 'DiseaseState'])['log_abun'].mean()
        meandf.name = 'mean'
        meandf = (
            meandf
            .reset_index()
            .pivot(index='denovo', columns='DiseaseState', values='mean')
            )
        meandf.columns = [i + '_logabun' for i in meandf.columns]
        meandf['delta_logabun'] = (
            meandf[case_lbl + '_logabun']
            - meandf[ctrl_lbl + '_logabun']
            )

        # Get the std of abundance over all samples
        # TODO: do I need to exclude the non-control or case samples here?
        stds = tidydf.groupby('denovo')['log_abun'].std()
        stds.name = 'std'

        # And differential presence/absence
        tidydf['present'] = tidydf['abun'] > 0
        mean_present = tidydf.groupby(['denovo', 'DiseaseState'])['present'].mean()
        mean_present.name = 'mean'
        mean_present = (
            mean_present
            .reset_index()
            .pivot(index='denovo', columns='DiseaseState', values='mean')
            )
        mean_present.columns = [i + '_presence' for i in mean_present.columns]

        mean_present['delta_presence'] = (
            mean_present[case_lbl + '_presence']
            - mean_present[ctrl_lbl + '_presence']
            )

        # Concatenate into signal-to-noise ratio df
        snr = pd.concat((stds, meandf), axis=1)

        # Calculate signal to noise and absolute value of it
        snr['snr'] = snr['delta_logabun'] / snr['std']
        snr['abs_snr'] = abs(snr['snr'])

        # Also add differential presence
        snr = pd.concat((snr, mean_present), axis=1)
        snr['abs_delta_presence'] = abs(snr['delta_presence'])

        # And rank them based on these differences
        snr = snr.sort_values(by='abs_snr', ascending=False)
        snr['rank_snr'] = range(snr.shape[0])
        snr = snr.sort_values(by='abs_delta_presence', ascending=False)
        snr['rank_delta_presence'] = range(snr.shape[0])

        print('Done.')

        if taxa == 'otu':
            snr.to_csv('data/analysis/population_effects.otu.{}.txt'.format(d),
                       sep='\t', index=True)
        else:
            snr.to_csv('data/analysis/population_effects.{}.txt'.format(d),
                       sep='\t', index=True)
