#!/usr/bin/env python

"""
This script takes in an OTU table and associated taxonomy file(s)
and writes a tidy OTU table with the taxonomy/ies as another variable.

All of the input files are assumed to be output from QIIME 2.

author: cduvallet@gmail.com
"""
import argparse

import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--otu-table', help='path to OTU table. OTUs in rows and '
        'samples in columns, as in the output feature tables from QIIME 2. '
        'Sample IDs are in the second row (first row has the "#Converted from '
        'BIOM" line and is skipped).')
    parser.add_argument('--silva', help='path to SILVA taxonomy assignments',
        default=None)
    parser.add_argument('--gg', help='path to GreenGenes taxonomy assignments',
        default=None)
    parser.add_argument('--tidy-file', help='path to write tidy OTU table to')
    return parser.parse_args()

def fix_gg_tax_string(t):
    """
    t is a string returned by the GG sklearn qiime2 classifier.
    This function splits t by semicolon and removes trailing spaces,
    and returns a taxonomic string with missing levels filled in.

    e.g.
    input:  'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae'
    output: 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__;s__'
    """

    t = [i.strip() for i in t.split(';')]
    missing_taxa = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

    return ';'.join(t + missing_taxa[len(t):])

def return_gg_otu(row):
    return fix_gg_tax_string(row['Taxon']) + ';d__' + row['Feature ID']

def fix_silva_tax_string(t):
    """
    t is a string returned by the GG sklearn qiime2 classifier.
    This function splits t by semicolon and removes trailing spaces,
    and returns a taxonomic string with missing levels filled in.

    e.g.
    input:  'D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Lachnospiraceae'
    output: 'D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Lachnospiraceae;D_5__;D_6__;D_7__e0a05e5465157a68916b185f11c680fd'
    """

    t = [''.join(i.strip().split()) for i in t.split(';')]
    missing_taxa = ['D_' + str(i) + '__' for i in range(7)]

    return ';'.join(t + missing_taxa[len(t):])

def return_silva_otu(row):
    return fix_silva_tax_string(row['Taxon']) + ';D_7__' + row['Feature ID']

if __name__ == "__main__":
    ## Parse args
    args = parse_args()

    ## OTU Table
    # Read in OTU table and transpose (so that OTUs become columns)
    df = pd.read_csv(args.otu_table, sep='\t',
                     skiprows=1, index_col=0).T
    # Get total reads per sample (to put in tidy table so we can
    # still convert to relative abundance later)
    df['total_reads'] = df.sum(axis=1)
    # Give a meaningful index name for later reshaping
    df.index.name = 'sample_id'

    ## GreenGenes Taxonomy
    if args.gg is not None:
        # Read in GG taxonomy file
        taxgg = pd.read_csv(args.gg, sep='\t')
        # Fix the OTU name - fill in missing taxonomic levels,
        # remove blank spaces, and append the denovo ID
        taxgg['otu_id'] = taxgg.apply(
            lambda row: return_gg_otu(row), axis=1)

    ## SILVA taxonomy
    if args.silva is not None:
        taxsilva = pd.read_csv(args.silva, sep='\t')
        # Fix the OTU name - fill in missing taxonomic levels,
        # remove blank spaces, and append the denovo ID
        taxsilva['otu_id'] = taxsilva.apply(
            lambda row: return_silva_otu(row), axis=1)

    ## Tidy the data
    # Convert OTU table to tidy form
    tidydf = (df
        .reset_index()
        .melt(id_vars=['sample_id', 'total_reads'],
              value_name='reads')
        )
    # Rename taxa df columns so they don't conflict later
    taxgg.columns = [i + '_gg' for i in taxgg.columns]
    taxsilva.columns = [i + '_silva' for i in taxsilva.columns]

    # Add the OTU ID from GG to the tidy OTU table
    tidydf = pd.merge(tidydf, taxgg,
                      how='left', left_on='#OTU ID',
                      right_on='Feature ID_gg')
    # And from SILVA
    tidydf = pd.merge(tidydf, taxsilva,
                      how='left', left_on='#OTU ID',
                      right_on='Feature ID_silva')

    ## Write output
    tidydf.to_csv(args.tidy_file, sep='\t', index=False)
