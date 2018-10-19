#!/usr/bin/env python

"""
Investigate the top-ranked OTUs based on signal-to-noise and differential
abundance.

This script reads in data and makes some plots:
- delta log(abun) vs std and signal-to-noise
- delta log(abun) vs. delta presence (both raw values and ranks)
- delta log(abun), signal-to-noise, and delta presence vs. qvalue
- rank signal-to-noise and rank delta presence vs. rank qvalue
- boxplots of abundances of the top 12 hits in cases and controls

"""

import pandas as pd
import numpy as np
import feather

import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-white')

def read_dataframe(fn):
    df = feather.read_dataframe(fn)
    df.index = df.iloc[:,0]
    df = df.iloc[:, 1:]
    return df

## Read in the qvalues
fqvals = 'data/analysis/power_simulation.otu_qvalues.30_reps.denovo_otu_only.feather'
print('Reading qvalues...')
qvals = feather.read_dataframe(fqvals)
print('Done.')

datasets = ['cdi_schubert', 'crc_baxter', 'ibd_papa', 'ob_goodrich']

labels = {'cdi_schubert': {'control': 'H', 'case': 'CDI'},
          'crc_baxter': {'control': 'H', 'case': 'CRC'},
          'ibd_papa': {'control': 'nonIBD', 'case': 'IBD'},
          'ob_goodrich': {'control': 'H', 'case': 'OB'}
          }

for d in datasets:
    print(d)
    fotu = 'data/clean/' + d + '.otu_table.genus.feather'
    fmeta = 'data/clean/' + d + '.metadata.feather'

    df = read_dataframe(fotu)
    df.index.name = 'sample_id'
    meta = read_dataframe(fmeta)

    ctrl_lbl = labels[d]['control']
    case_lbl = labels[d]['case']

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

    # Read in signal-to-noise and other effects dataframe
    snr = pd.read_csv('data/analysis/population_effects.{}.txt'.format(d),
                      sep='\t', index_col=0)

    # Add q-values to the dataframe
    snr = pd.merge(
        qvals.query('study == @d')[['denovo', 'q_allsamples']].drop_duplicates(),
        snr,
        left_on = 'denovo',
        right_index = True
    )

    # Add rank based on qvalues
    snr = snr.sort_values(by='q_allsamples')
    snr['rank_qvalue'] = range(snr.shape[0])

    print('Done.')

    ## Make plots
    # delta log(abun) vs std and signal-to-noise
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    ax[0].scatter(snr['std'], snr['delta_logabun'])
    ax[0].set_xlabel('std')
    ax[0].set_ylabel('delta log(abun)')

    ax[1].scatter(snr['delta_logabun'], snr['snr'])
    ax[1].set_xlabel('delta log(abun)')
    ax[1].set_ylabel('logabun / std(logabun)')
    fig.tight_layout()
    fig.savefig('figures/analysis/delta_abun_vs_std_snr.genus.{}.png'.format(d))

    # delta log(abun) vs. delta presence (both raw values and ranks)
    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].scatter(snr['delta_logabun'], snr['delta_presence'])
    ax[0].set_xlabel('delta log(abun)')
    ax[0].set_ylabel('delta presence')

    ax[1].scatter(snr['rank_snr'], snr['rank_delta_presence'])
    ax[1].set_xlabel('rank delta log(abun)')
    ax[1].set_ylabel('rank delta presence')
    fig.tight_layout()
    fig.savefig('figures/analysis/delta_abun_vs_delta_presence.genus.{}.png'.format(d))

    # delta log(abun), signal-to-noise, and delta presence vs. qvalue
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    ax[0].scatter(snr['delta_logabun'], -np.log10(snr['q_allsamples']))
    ax[0].set_xlabel('delta log(abun)')
    ax[0].set_ylabel('-log(qval)')

    ax[1].scatter(snr['snr'], -np.log10(snr['q_allsamples']))
    ax[1].set_xlabel('signal to noise')
    ax[1].set_ylabel('-log(qval)')

    ax[2].scatter(snr['delta_presence'], -np.log10(snr['q_allsamples']))
    ax[2].set_xlabel('delta presence')
    ax[2].set_ylabel('-log(qval)')
    fig.tight_layout()
    fig.savefig('figures/analysis/deltas_vs_qvalues.genus.{}.png'.format(d))

    # rank signal-to-noise and rank delta presence vs. rank qvalue
    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].scatter(snr['rank_snr'], snr['rank_qvalue'])
    ax[0].set_xlabel('rank snr')
    ax[0].set_ylabel('rank qval')

    ax[1].scatter(snr['rank_delta_presence'], snr['rank_qvalue'])
    ax[1].set_xlabel('rank delta presence')
    ax[1].set_ylabel('rank qval')
    fig.tight_layout()
    fig.savefig('figures/analysis/rank_snr_vs_rank_delta_presence.genus.{}.png'.format(d))

    ## Boxplot of abundances for the top 12 hits
    # Based on signal-to-noise ratios
    top12 = (
        snr
        .sort_values(by='abs_snr', ascending=False)
        .head(12)
        ['denovo']
        .tolist()
        )

    order = [case_lbl, ctrl_lbl]

    g = sns.FacetGrid(data=tidydf.query('denovo == @top12'),
                      col='denovo', col_wrap=4,
                      sharey=False)
    g.map(sns.boxplot, 'DiseaseState', 'abun', order=order)
    g.map(sns.stripplot, 'DiseaseState', 'abun', order=order)
    plt.tight_layout()
    plt.savefig('figures/analysis/top12_snr.genus.{}.png'.format(d))

    # Based on differential presence
    top12 =(
        snr
        .sort_values(by='rank_delta_presence')
        .head(12)
        ['denovo']
        .tolist()
        )

    g = sns.FacetGrid(data=tidydf.query('denovo == @top12'),
                      col='denovo', col_wrap=4,
                      sharey=False)
    g.map(sns.boxplot, 'DiseaseState', 'abun', order=order)
    g.map(sns.stripplot, 'DiseaseState', 'abun', order=order)

    plt.tight_layout()
    plt.savefig('figures/analysis/top12_diff_presence.genus.{}.png'.format(d))
