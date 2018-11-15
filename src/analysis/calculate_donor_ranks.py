#!/usr/bin/env python
"""
This script calculates the donor ranks by SCFA production and
secondary bile acid conversion.
"""

import pandas as pd
import feather
import numpy as np
import argparse


p = argparse.ArgumentParser()
p.add_argument('ftidy', help='path to tidy metabolomics file')
#p.add_argument('bile_tidy_out',
#    help='path to tidy output file with bile acid conversion ratios')
p.add_argument('wide_df_out',
    help='path to wide dataframe with SCFA and bile acid conversion ' +
         'ratios and ranks')
args = p.parse_args()

DF = feather.read_dataframe(args.ftidy)

# Grab the molecules of interest
scfas = ['propionate', 'butyrate', 'isovalerate']
primary_bile = ['cholate', 'chenodeoxycholate']
secondary_bile = ['deoxycholate', 'lithocholate']

all_mlcls = scfas + primary_bile + secondary_bile

df = DF.query('metabolite == @all_mlcls')

# Some donors have multiple samples; take the mean abundance
df = (
    df.groupby(['donor', 'metabolite'])
      ['intensity']
      .mean()
      .reset_index(name='mean_intensity')
)

## Bile acid conversion calculations
# LCA + DCA / CDCA + CA
tmp = []
for d, subdf in df.groupby('donor'):
    total_primary = (
        subdf
            .query('metabolite == @primary_bile')
            .sum()
            ['mean_intensity']
        )
    total_secondary = (
        subdf
            .query('metabolite == @secondary_bile')
            .sum()
            ['mean_intensity']
        )
    tmp.append([d, total_primary, total_secondary])

bile_comp = pd.DataFrame(
    tmp, columns=['donor', 'total_primary', 'total_secondary'])
bile_comp['secondary_to_primary'] = bile_comp['total_secondary']/bile_comp['total_primary']

# LCA / CDCA
m2 = "lithocholate"
m1 = 'chenodeoxycholate'
col = m2 + '_' + m1
tmp_ratio = (df
     .query('(metabolite == @m2) | metabolite == @m1')
     .pivot(index='donor', columns='metabolite', values='mean_intensity')
)
tmp_ratio[col] = tmp_ratio[m2] / tmp_ratio[m1]
bile_comp = pd.merge(bile_comp, tmp_ratio, left_on='donor', right_index=True)

# DCA / CA
m2 = "deoxycholate"
m1 = 'cholate'
col = m2 + '_' + m1
tmp_ratio = (df
     .query('(metabolite == @m2) | metabolite == @m1')
     .pivot(index='donor', columns='metabolite', values='mean_intensity')
)
tmp_ratio[col] = tmp_ratio[m2] / tmp_ratio[m1]
bile_comp = pd.merge(bile_comp, tmp_ratio, left_on='donor', right_index=True)

# Add donor ranks, based on the total secondary to total primary ratio
bile_comp['secondary_to_primary_rank'] = bile_comp['secondary_to_primary'].rank(method='dense')

# bile_comp is a wide dataframe with 'donor' column and the values
# for each bile acid, bile acid conversion ratio, and the ranks for
# the secondary-to-primary ratio

## Get back into tidy format for plotting and write to file
#bile_tidy = bile_comp.melt(id_vars='donor', value_name='value')
#bile_tidy.to_csv(args.bile_tidy_out, sep='\t', index=False)

## Add SCFA to the wide dataframe

# Convert original df to wide dataframe for easy ranking
df_wide = (
    df.query('metabolite == @scfas')
      .pivot(index='donor', columns='metabolite', values='mean_intensity')
    )
# Convert ion intensities to ranks
# method=dense: like 'min', but rank always increases by 1 between groups
df_ranks = df_wide.fillna(0).rank(axis=0, method='dense')
df_ranks.columns = [i + '_rank' for i in df_ranks.columns]

# Re-melt into tidy form
tidy_ranks = pd.melt(df_ranks.reset_index(),
    id_vars='donor', value_name='rank')

# Rank based on avg SCFA ranks
scfa_ranked = (
    tidy_ranks
#        .query('variable == @keep_vars')
        .groupby('donor')
        .mean()
        .reset_index()
        .sort_values(by='rank')
    )
scfa_ranked = scfa_ranked.rename(columns={'rank': 'avg_scfa_rank'})

# Combine the bile acid wide dataframe with: original SCFA mean intenstity,
# rank in each SCFA, and avg rank across SCFAs
df_wide_final = (
    pd.merge(bile_comp, df_wide.fillna(0), left_on='donor', right_index=True)
      .merge(df_ranks, left_on='donor', right_index=True)
      .merge(scfa_ranked)
    )
df_wide_final.to_csv(args.wide_df_out, sep='\t', index=False)
