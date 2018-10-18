#!/usr/bin/env python
"""
This script runs through the "power simulation," subsampling
different total N patients and changing the proportion of cases
and controls (e.g. FMT responders and non-responders), to show
how different sample sizes are needed to detect a signal in
differentially abundant OTUs.
"""
import pandas as pd
import numpy as np
import feather

import multiprocessing

from scipy.stats.mstats import kruskalwallis
from scipy.stats import ranksums, mannwhitneyu
# FDR correction
from statsmodels.sandbox.stats.multicomp import multipletests

# These function are adapted from util.py in microbiomeHD
def univariate_one_col(xvals, yvals, pfun):
    try:
        h, p = pfun(xvals, yvals)
    except ValueError:
        p = 1
        h = 0
    return pd.Series([p, h])

def compare_otus_teststat(df, Xsmpls, Ysmpls,
                          method='kruskal-wallis',
                          multi_comp=None):
    """
    Compares columns between Xsmpls and Ysmpls, with statistical method=method.
    Returns dataframe with both the qvals ('p') and test statistic ('test-stat')

    parameters
    ----------
    df             dataframe, samples are in rows and OTUs in columns
    X,Ysmpls       list of samples to compare
    method         statistical method to use for comparison
    multi_comp     str, type of multiple comparison test to do.
                   Currently accepts 'fdr' or None

    outputs
    -------
    results        dataframe with OTUs in rows and 'p' and 'test-stat' in columns

    """
    if method == 'kruskal-wallis':
        pfun = kruskalwallis
    elif method == 'wilcoxon' or method == 'ranksums':
        pfun = ranksums
    elif method == 'mann-whitney':
        pfun = mannwhitneyu
        # Note: prob wanna add some kwargs here to say whether 2sided or not

    results = df.apply(
        lambda col: univariate_one_col(
            col.loc[Xsmpls].values,
            col.loc[Ysmpls].values,
            pfun)
            ).T
    results.columns = ['p', 'test_stat']

    if multi_comp == 'fdr':
        _, results['q'], _, _ = multipletests(results['p'], method='fdr_bh')

    return results

def read_dataframe(fn):
    df = feather.read_dataframe(fn)
    df.index = df.iloc[:,0]
    df = df.iloc[:, 1:]

    return df

def make_param_combos(totalNs, perc_success, ctrl, case):
    """
    Make the parameter combinations we'll test.

    Parameters
    ----------
    totalNs : list of ints
        list of total N's to test
    perc_successes : list of floats
        the percent of samples which will be cases (i.e. FMT success rate)
    ctrl, case : lists of str
        the full case and control sample lists

    Returns
    -------
    n, p, n_ctrl, n_case : lists
        all parameter combinations
    """
    N = []
    P = []
    N_CTRL = []
    N_CASE = []

    for n in totalNs:
        for p in perc_success:
            # Define number of cases and controls
            n_ctrl = int(np.ceil((1.0 - p)*n))
            n_case = int(np.ceil(p*n))

            if (n_ctrl > len(ctrl)) or (n_case > len(case)):
                # If you're trying to sample more patients than
                # are available, skip that parameter set
                continue

            N.append(n)
            P.append(p)
            N_CTRL.append(n_ctrl)
            N_CASE.append(n_case)

    return N, P, N_CTRL, N_CASE

def subsample_and_get_pvals(df, genusdf, ctrl, case, n_ctrl, n_case):
    # Subsample our cases and controls
    subctrl = np.random.choice(ctrl, size=n_ctrl)
    subcase = np.random.choice(case, size=n_case)

    # Calculate qvalues, OTU-level
    # For now, let's make it the same as genus df
    psub = compare_otus_teststat(
        df, subctrl, subcase,
        method='kruskal-wallis', multi_comp='fdr')

    # Calculate qvalues, genus-level
    psubgenus = compare_otus_teststat(
        genusdf, subctrl, subcase,
        method='kruskal-wallis', multi_comp='fdr')

    return psub, psubgenus

def get_significant(pvals, alphas=[0.05, 0.1, 0.25]):
    """
    Calculate the number of OTUs which are significant
    a 'q' <= alphas
    """
    res = []
    for a in alphas:
        nsig = sum(pvals['q'] <= a)
        res.append([a, nsig])
    return pd.DataFrame(data=res, columns=['alpha', 'n_sig'])

def add_metadata_to_sigdf(sigdf, n_ctrl, n_case, n, p, taxa_level):
    sigdf['taxa_level'] = taxa_level
    sigdf['n_ctrl'] = n_ctrl
    sigdf['n_case'] = n_case
    sigdf['total_n'] = n
    sigdf['perc_case'] = p
    return sigdf

def format_sig_results(potus, pgenus, n_ctrl, n_case, n, p):
    # Calculate number OTUs/genera significant at different levels
    sigdf = get_significant(potus)
    sigdf = add_metadata_to_sigdf(sigdf, n_ctrl, n_case, n, p, 'otu')

    genussig = get_significant(pgenus)
    genussig = add_metadata_to_sigdf(genussig, n_ctrl, n_case, n, p, 'genus')

    sigdf = pd.concat((sigdf, genussig), ignore_index=True)

    # Add metadata to the full q-value results
    psub = add_metadata_to_sigdf(potus, n_ctrl, n_case, n, p, 'otu')
    psubgenus = add_metadata_to_sigdf(pgenus, n_ctrl, n_case, n, p, 'genus')
    psub = pd.concat((psub, psubgenus))

    return sigdf, psub

def parallel_process((df, genusdf, ctrl, case, n_ctrl, n_case, n, p)):
    """
    Returns:

    sigdf : columns alpha  n_sig taxa_level  n_ctrl  n_case  total_n  perc_case
    psub : 'p', 'test_stat', 'q', 'taxa_level', 'n_ctrl', 'n_case',
           'total_n', 'perc_case'
    """

    print(n, p, n_case, n_ctrl)
    ## Get the pvalues for subsampled OTU and genus-level
    psub, psubgenus = subsample_and_get_pvals(
        df, genusdf, ctrl, case, n_ctrl, n_case)

    ## Format results
    sigdf, psub = format_sig_results(psub, psubgenus, n_ctrl, n_case, n, p)
    return sigdf, psub

## Hard coded values and parameters
n_reps = 30

fout_qvalues = 'power_simulation.otu_qvalues.{}_reps.denovo_otu_only.feather'.format(n_reps)
fout_nsig = 'data/analysis/power_simulation.n_sig.{}_reps.txt'.format(n_reps)

np.random.seed(12345)

# These values should reflect what we'd expect
# from "reasonable" clinical trials
totalNs = [10, 25, 50, 100, 150, 200]
perc_success = [0.1, 0.25, 0.5, 0.75, 0.9]

### DEBUGGING ONLY
#totalNs = [10, 25]#, 50, 100, 150, 200]
#perc_success = [0.1, 0.25]#, 0.5, 0.75, 0.9]

# Note: different datasets represent different "expected effect sizes"
DATASETS = ['cdi_schubert', 'crc_baxter', 'ibd_papa', 'ob_goodrich']
CTRLS = {'cdi_schubert': 'H',
         'crc_baxter': 'H',
         'ibd_papa': 'nonIBD',
         'ob_goodrich': 'H'}
CASES = {'cdi_schubert': 'CDI',
         'crc_baxter': 'CRC',
         'ibd_papa': ['CD', 'UC'],
         'ob_goodrich': 'OB'}

all_qvals = []
all_nsigs = []
for r in range(n_reps):
    print(r)
    for dataset in DATASETS:
        print(dataset)
        ## Read in dataset
        fotu = 'data/clean/' + dataset + '.otu_table.feather'
        fgenus = 'data/clean/' + dataset + '.otu_table.genus.feather'
        fmeta = 'data/clean/' + dataset + '.metadata.feather'

        df, genusdf, meta = (read_dataframe(f) for f in [fotu, fgenus, fmeta])

        ## Set up cases and controls
        ctrl_lbl = CTRLS[dataset]
        case_lbl = CASES[dataset]
        ctrl = meta.query('DiseaseState == @ctrl_lbl').index.tolist()
        case = meta.query('DiseaseState == @case_lbl').index.tolist()

        ## Set up simulation parameters
        N, P, N_CTRL, N_CASE = make_param_combos(
            totalNs, perc_success, ctrl, case)

        ## Calculate pvalues
        p = multiprocessing.Pool()
        allres = p.map(
            parallel_process,
            zip(len(N)*[df],
                len(N)*[genusdf],
                len(N)*[ctrl],
                len(N)*[case],
                N_CTRL, N_CASE,
                N, P))
        p.close()
        # Not sure I need this (not sure what it does...)
        p.join()

        ## Put the parallelized results into two dataframes
        # sig_results and qval_results are both lists of dataframes
        sig_results = pd.concat([i[0] for i in allres], ignore_index=True)
        qval_results = (
            pd.concat([i[1] for i in allres])
            .reset_index()
            .rename(columns={'index': 'otu'})
            )

        ## Calculate pvalues/nsig/etc with original dataset,
        ## and add that to the results
        potus = compare_otus_teststat(
            df, ctrl, case,
            method='kruskal-wallis', multi_comp='fdr')
        pgenus = compare_otus_teststat(
            genusdf, ctrl, case,
            method='kruskal-wallis', multi_comp='fdr')

        # Format
        sigall, pvalsall = format_sig_results(
            potus, pgenus,
            len(ctrl), len(case),
            len(ctrl)+len(case), np.nan)
        pvalsall = pvalsall.reset_index().rename(columns={'index': 'otu'})
        pvalsall = pvalsall.rename(columns={'p': 'p_allsamples',
                                    'q': 'q_allsamples',
                                    'test_stat': 'test_stat_allsamples'})
        sigall = sigall.rename(columns={'n_sig': 'n_sig_allsamples'})

        # Concatenate with the other results full results, so that each
        # combination of parameters has additional columns with the results
        # for all samples

        # This merges on the common columns 'otu' and 'taxa_level'
        qval_results = pd.merge(
            pvalsall[['otu', 'taxa_level', 'p_allsamples', 'q_allsamples']],
            qval_results
        )

        sig_results = pd.merge(
            sigall[['alpha', 'n_sig_allsamples', 'taxa_level']],
            sig_results
        )

        ## Add study label
        qval_results['study'] = dataset
        sig_results['study'] = dataset

        ## Add rep label
        qval_results['rep'] = r
        sig_results['rep'] = r

        all_qvals.append(qval_results)
        all_nsigs.append(sig_results)

## Finally, combine all datasets' results together
all_qvals_df = pd.concat(all_qvals)
all_nsigs_df = pd.concat(all_nsigs)

## Prepare the qvalues dataframe for writing
## Remove the long OTU strings which are too large for feather to handle
all_qvals_df['denovo'] = all_qvals_df['otu'].str.rsplit(';', 1).str[1]
all_qvals_df = all_qvals_df.drop('otu', axis=1)

## Write to files
feather.write_dataframe(all_qvals_df, fout_qvalues)
all_nsigs_df.to_csv(fout_nsig, sep='\t', index=False)
