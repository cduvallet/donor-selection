#!/usr/bin/env python
"""
Script to clean up the OTU tables, using same parameters and steps as I used
in meta-analysis paper (but with a much less intense code-base since I don't
need to do this for 30 studies).

Steps:
    df, meta = clean_up_samples(df, meta, y[dataset_id])
    --> keeps only samples which have both 16S data and metadata
    --> also removes any non-stool and/or duplicate time point samples

    df, meta = clean_up_tables(df, meta, args.n_reads_otu, args.n_reads_sample, args.perc_samples)
    -->     Removes samples with fewer than n_reads_sample reads. (default 100)
            Removes OTUs with fewer than n_reads_otu reads. (default 10)
            Removes OTUs which are present in fewer than perc_samples*100 percent of samples. (default 0.01)
            Removes empty samples/OTUs.

    And writes tables in feather format (to speed up IO in downstream analyses)

Most of this code is copied from my previous code using
this data: https://github.com/cduvallet/microbiomeHD/blob/master/src/data/clean_otu_and_metadata.py

I'll also convert to relative abundance before writing the file.
And I'll collapse to genus level and save that as a clean file as well.

"""
import pandas as pd
import feather
import copy
import os

import argparse

def define_file_paths(dataset, rawdir):
    """
    Returns the file path for the otu table and metadata
    """
    otupath = [rawdir,
               dataset + '_results',
               'RDP',
               dataset + '.otu_table.100.denovo.rdp_assigned']
    otupath = os.path.join(*otupath)

    metapath = [rawdir,
                dataset + '_results',
                dataset + '.metadata.txt']
    metapath = os.path.join(*metapath)

    return otupath, metapath

def read_raw_files(otufile, metafile):

    df = pd.read_csv(otufile, sep='\t', index_col=0)
    meta = pd.read_csv(metafile, sep='\t', index_col=0)

    # If either index wasn't read as a string, explicitly do so
    if meta.index.dtype != 'O':
        meta.index = pd.read_csv(metafile, sep='\t', dtype=str).iloc[:,0]

    if df.index.dtype != 'O':
        df.index = pd.read_csv(clean_otu_file, sep='\t', dtype=str).iloc[:,0]

    return df.T, meta

def clean_up_samples(df, meta, conditions=None):
    """
    Cleans up samples in the OTU table and metadata dataframes.
    Keeps only samples which have both metadata and 16S data. If a condition
    is given, also keeps only samples which have the specified condition.

    Parameters
    ----------
    df, meta : pandas dataframes
        samples in rows, OTUs/metadata labels in columns.
    conditions : dict
        dictionary with {metadata_column: [values to keep]}
        metadata_column should be a str, and there can be multiple
        (i.e. multiple columns)

    Returns
    -------
    df, meta : pandas dataframes
        dataframes containing only samples which have both
        16S data and metadata. If a subset of metadata was
        specified in 'condition' in the yaml file, only returns
        the subset of samples that meet that condition.
    """

    print('Original: {} samples with 16S, '.format(df.shape[0]) +
          '{} samples with metadata.'.format(meta.shape[0]))

    # If a condition is given, keep only samples with that condition
    if conditions is not None:
        for col in conditions:
            print('Keeping only samples with values ' +
                  '{}'.format(', '.join([str(i) for i in conditions[col]]) +
                  ' for metadata column {}'.format(col))
                  )
            meta = meta[meta[col].isin(conditions[col])]
            print('{} samples left in metadata'.format(meta.shape[0]))

    # Remove samples which don't have both 16S and metadata
    keepsmpls = [i for i in df.index if i in meta.index]

    if len(keepsmpls) != len(df.index) or len(keepsmpls) != len(meta.index):
        print('Dataset has {} samples with 16S, '.format(len(df.index)) +
              '{} samples with metadata, '.format(len(meta.index)) +
              'but only {} samples with both'.format(len(keepsmpls))
              )

    df = df.loc[keepsmpls]
    meta = meta.loc[keepsmpls]

    return df, meta

def clean_up_tables(df, meta, n_reads_otu, n_reads_sample, perc_samples):
    """
    Cleans up the OTU table and metadata dataframes in data.
    Removes samples with fewer than n_reads_sample reads.
    Removes OTUs with fewer than n_reads_otu reads.
    Removes OTUs which are present in fewer than perc_samples*100 percent of samples.
    Removes empty samples/OTUs.
    """

    # Remove samples with fewer than n_reads reads.
    df = remove_shallow_smpls(df, n_reads_sample)

    # Remove OTUs with fewer than 10 reads
    old = df.shape[1]
    df = remove_shallow_otus(df, n_reads=n_reads_otu)
    new = df.shape[1]
    if new < old:
        print('Of {} original OTUs, {}'.format(old, new) +
              ' have more than {} reads'.format(n_reads_otu)
              )

    # Remove OTUs which are present in fewer than perc_samples of samples.
    old = df.shape[1]
    df = remove_shallow_otus(df, perc_samples=perc_samples)
    new = df.shape[1]
    if new < old:
        print('Of {} original OTUs, {} are present '.format(old, new) +
              'in more than {}% of samples'.format(perc_samples*100)
              )

    # Remove any samples which now have fewer than n_reads
    df = remove_shallow_smpls(df, n_reads_sample)

    # Double check that both metadata and OTU table have same samples
    # (after we've filtered out samples based on reads, etc)
    # if len(meta.index) > len(df.index):
    #     keepsmpls = [i for i in meta.index if i in df.index]
    # else:
    #     keepsmpls = [i for i in df.index if i in meta.index]
    keepsmpls = [i for i in df.index if i in meta.index]

    if len(keepsmpls) != len(df.index) or len(keepsmpls) != len(meta.index):
        print('After some cleaning, dataset has {} '.format(len(df.index)) +
              ' samples with 16S, {} samples '.format(len(meta.index)) +
              'with metadata, but only {} samples'.format(len(keepsmpls)) +
              ' with both.'
              )

        df = df.loc[keepsmpls]
        meta = meta.loc[keepsmpls]

    # Remove empty metadata columns (just in case...)
    meta = meta.dropna(how='all', axis=1)

    return df, meta

def remove_shallow_smpls(df, n_reads):
    """
    Removes samples with fewer than n_reads from dataframe df.

    Parameters
    -----------
    df : pandas dataframe
        samples are in rows, OTUs are in columns
    n_reads : int
        minimum number of reads per sample for sample to be kept
    """

    total_reads = df.sum(axis=1)
    shallow_smpls = [smpl for smpl in total_reads.index \
                     if total_reads.loc[smpl] <= n_reads]
    df = df.drop(shallow_smpls)

    return df

def remove_shallow_otus(df, perc_samples=None, n_reads=None):
    """
    Removes OTUs which are present in fewer than 100*perc_samples percent of
    samples OR which have fewer than n_reads reads.

    Parameters
    ----------
    df : pandas dataframe
        Samples are in rows. OTUs are in columns.
    perc_samples : float
        min percent of samples that an OTU must be present in to not
        be thrown out.
    n_reads : int
        min number of reads an OTU must have in df to not be thrown
        out.

    Either perc_samples or n_reads must be specified. If both are specified,
    the perc_samples filtering is done first and then OTUs with fewer than
    n_reads total are thrown out.

    """
    if perc_samples is not None:
        presencemap = lambda x: 1 if x else 0
        otus_perc_present = df.applymap(presencemap).sum() / df.shape[0]
        keepotus = list(
            otus_perc_present[otus_perc_present > perc_samples].index)
        df = df[keepotus]

    if n_reads is not None:
        # Removes any OTUs with fewer than n_reads from the raw and abun dfs
        # samples are in rows and OTUs are in columns
        total_reads = df.sum(axis=0)
        shallow_col_indices = [i for i in range(len(total_reads.index)) \
                               if total_reads.iloc[i] < n_reads]
        shallow_otus = df.columns[shallow_col_indices]
        df = df.drop(shallow_otus, axis=1)

    return df

def collapse_to_genus(OTU_table):
    """
    Collapses OTU table to genus level by string-matching.
    This is still based on Thomas Gurry's code (and can
    certainly be improved...)

    Parameters
    ----------
    df : pandas dataframe
        OTUs in columns, samples in rows.
        Taxonomic levels in OTU strings should be semicolon-delimited,
        starting with kingdom level.
        Unannotated taxonomic levels should end with '__' (e.g. ''...;g__Roseburia;s__;d__denovo123')

    Returns
    -------
    newdf : pandas dataframe
        OTUs in columns, samples in rows.
        OTUs are collapsed to the given taxonomic level.
        Matching values (for annotated taxa) are summed for each sample.
        Values corresponding to unannotated taxa are discarded.
    """

    OTU_IDs = list(OTU_table.columns)

    # Collapse to the right level
    OTU_taxa = [';'.join(OTU_ID.split(';')[:6]) for OTU_ID in OTU_IDs]

    # Get indices of each unique taxon
    taxa_indices = {}
    for i in range(len(OTU_taxa)):
        if ((OTU_taxa[i] not in taxa_indices)
                and (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__")):
            taxa_indices[OTU_taxa[i]] = []
        if (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
            taxa_indices[OTU_taxa[i]].append(i)
    # Make new empty df with the same samples as original df and taxa in taxa_indices.keys()
    newdf = pd.DataFrame(
        index=df.index, columns=taxa_indices.keys(), data=0)

    # Get sample contents for each taxa of the chosen level and put into newdf
    for key in taxa_indices:
        indices = taxa_indices[key]
        newcol = OTU_table.iloc[:, indices].sum(axis=1)
        newdf[key] = copy.copy(newcol)

    return newdf

def define_clean_paths(cleandir, dataset):
    """
    Returns the file path for the clean otu table and metadata
    """
    otupath = [cleandir,
            dataset + '.otu_table.feather']
    otupath = os.path.join(*otupath)

    genuspath = [cleandir,
            dataset + '.otu_table.genus.feather']
    genuspath = os.path.join(*genuspath)

    metapath = [cleandir,
             dataset + '.metadata.feather']
    metapath = os.path.join(*metapath)

    return otupath, genuspath, metapath

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('dataset')
    args = p.parse_args()

    # Define data directories, relative to the current path
    rawdir = 'data/raw'
    cleandir = 'data/clean'

    # Define cleaning parameters
    n_reads_otu = 10
    n_reads_sample = 100
    perc_samples = 0.01

    # Condition says which samples to keep based on a metadata column
    conditions = {'cdi_schubert': None,
                  'crc_baxter': None,
                  'ibd_alm': None,
                  'ob_goodrich': {'n_sample': [0]}
                  }

    dataset = args.dataset
    condition = conditions[dataset]
    print(dataset)
    # Read in data and metadata
    df, meta = read_raw_files(*define_file_paths(dataset, rawdir))
    # Clean up the samples
    df, meta = clean_up_samples(df, meta, condition)
    # Clean up the OTU table
    df, meta = clean_up_tables(df, meta,
        n_reads_otu, n_reads_sample, perc_samples)

    # Convert to relative abundance
    df = df.divide(df.sum(axis=1), axis=0)

    # Collapse to genus level
    genusdf = collapse_to_genus(df)

    # Write files
    if dataset == "ibd_alm":
        dataset = "ibd_papa"

        # Also change the DiseaseState to be IBD and nonIBD
        meta['DetailedDiseaseState'] = meta['DiseaseState']
        meta['DiseaseState'] = (
            meta['DiseaseState']
            .apply(lambda x: 'nonIBD' if x == 'nonIBD' else 'IBD')
            )

    dfpath, genuspath, metapath = define_clean_paths(cleandir, dataset)

    feather.write_dataframe(df.reset_index(), dfpath)
    feather.write_dataframe(genusdf.reset_index(), genuspath)
    feather.write_dataframe(meta.reset_index(), metapath)
