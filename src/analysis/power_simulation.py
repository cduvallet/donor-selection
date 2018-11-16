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

def get_nsig(potus, pgenus, n_ctrl, n_case, n, p):
    """
    Calculates number of OTUs/genera significant at different
    taxonomic levels. Adds metadata (i.e. simulations settings)
    to the output dataframe.
    """

    # Calculate number OTUs/genera significant at different levels
    sigdf = get_significant(potus)
    sigdf = add_metadata_to_sigdf(sigdf, n_ctrl, n_case, n, p, 'otu')

    genussig = get_significant(pgenus)
    genussig = add_metadata_to_sigdf(genussig, n_ctrl, n_case, n, p, 'genus')

    sigdf = pd.concat((sigdf, genussig), ignore_index=True)

    return sigdf

def get_tophits_rejected(psubgenus, effects):
    """
    Get the number of "top hits" which were rejected. Top hits are defined
    by the 'effects' dataframe. The different numbers of top hits
    iterated over and the alpha cutoff are hard-coded in this function.

    Parameters
    ----------
    psubgenus : pandas DataFrame
        genus names in index, has column 'q'
    effects : pandas DataFrame
        dataframe with the "true" populatin effects. Has column 'rank_snr'
        which ranks genera by their signal to noise ratio. Has column
        'denovo' which contains the genus name (e.g. 'g__Ruminococcus')

    Returns
    -------
    tophits_df : pandas DataFrame
        has columns 'n_top_hits' and 'n_rejected'
    """
    # Hard code these
    alpha = 0.05
    ntophits = [3, 5, 10]

    # Add the short OTU name to dataframe (long OTU IDs are currently index)
    psubgenus['denovo'] = (
        psubgenus
            .reset_index()
            ['index']
            .str.rsplit(';', 1)
            .str[1]
            .values
        )
    n_reject = []
    for n in ntophits:
        # Define the top hits (genera)
        top_genera = effects.query('rank_snr < @n')['denovo'].tolist()

        # Get the qvalues for the top hits
        top_qvals = psubgenus.query('denovo == @top_genera')

        # Count how many were rejected
        n_reject.append(sum(top_qvals['q'] <= alpha))

    ## Put both lists into a dataframe
    tophits_df = pd.DataFrame({'n_top_hits': ntophits,
                              'n_rejected': n_reject})
    return tophits_df

def parallel_process((df, genusdf, effects,
                      ctrl, case, n_ctrl, n_case, total_n, p)):
    """
    Parameters
    ----------
    df : pandas DataFrame
    genusdf : pandas DataFrame
        OTU tables with OTUs/genera in columns and samples in index

    ctrl, case : list
        list of sample IDs, contains all control and case samples

    n_ctrl, n_case : int
        number of control and case samples to subsample to

    total_n : int
        number of total samples in the FMT arm. Only used to label results.

    p : float
        percent of total_n which are case. Only used to label results.

    effects : pandas DataFrame
        dataframe with the "true" effect size of each OTU in the population.
        This dataframe contains only results for genus level (genera are
        labeled in the column 'denovo').

    Returns
    -------
    sigdf : pandas DataFrame
        columns alpha  n_sig taxa_level  n_ctrl  n_case  total_n  perc_case
    tophitsdf : pandas DataFrame
        columns n_top_hits n_rejected (and all the other simulation setting
        columns, as above)
    """

    print(total_n, p, n_case, n_ctrl)
    ## Get the pvalues for subsampled OTU and genus-level
    psub, psubgenus = subsample_and_get_pvals(
        df, genusdf, ctrl, case, n_ctrl, n_case)

    ## Get total number significant
    sigdf = get_nsig(psub, psubgenus, n_ctrl, n_case, total_n, p)

    ## Grab number of significant top hits (genus level only)
    tophitsdf = get_tophits_rejected(psubgenus, effects)
    tophitsdf = add_metadata_to_sigdf(
        tophitsdf, n_ctrl, n_case, total_n, p, 'genus')

    return sigdf, tophitsdf

def run_simulation(n_reps, DATASETS, CTRLS, CASES, totalNs, perc_success):
    """
    Run the entire simulation for each dataset and parameter combo.

    Parameters
    ----------
    n_reps : int
    DATASETS : list of strings
    CTRLS, CASES : dicts
        {dataset_id: str or [list of strs]}
        strs are the case or control labels in the DiseaseState column
    totalNs : list of ints
        list of different total numbers of patients in FMT arm
    perc_successes : list of floats
        list of different "percent cases" (i.e. FMT response rate)

    Returns
    -------
    all_qvals_df : pandas DataFrame
        VERY LARGE dataframe with the qvalue for each rep of each OTU in
        each dataset/simulation setting. Has a 'denovo' column (rather than
        an 'otu' column).
        Columns: taxa_level, p_allsamples, q_allsamples, p, test_stat, q
                 n_ctrl, n_case, total_n, perc_case, study, rep, denovo
    all_nsigs_df : pandas DataFrame
        dataframe with the number of significant hits (q < 0.05) per
        simulation setting.
    """
    all_tophits = []
    all_nsigs = []
    for r in range(n_reps):
        print(r)
        for dataset in DATASETS:
            print(dataset)
            ## Read in dataset
            fotu = 'data/clean/' + dataset + '.otu_table.feather'
            fgenus = 'data/clean/' + dataset + '.otu_table.genus.feather'
            fmeta = 'data/clean/' + dataset + '.metadata.feather'
            feffects = 'data/analysis/population_effects.' + dataset + '.txt'

            df = read_dataframe(fotu)
            genusdf = read_dataframe(fgenus)
            meta = read_dataframe(fmeta)
            effects = pd.read_csv(feffects, sep='\t')

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
                    len(N)*[effects],
                    len(N)*[ctrl],
                    len(N)*[case],
                    N_CTRL, N_CASE,
                    N, P))
            p.close()
            # Not sure I need this (not sure what it does...)
            p.join()

            ## Put the parallelized results into two dataframes
            # sig_results and qval_results are both lists of dataframes
            nsig_results = (
                pd.concat([i[0] for i in allres], ignore_index=True)
                )
            tophits_results = (
                pd.concat([i[1] for i in allres], ignore_index=True)
                )

            ## Calculate nsig with original dataset and add to results
            potus = compare_otus_teststat(
                df, ctrl, case,
                method='kruskal-wallis', multi_comp='fdr')
            pgenus = compare_otus_teststat(
                genusdf, ctrl, case,
                method='kruskal-wallis', multi_comp='fdr')

            # Format
            sigall = get_nsig(
                potus, pgenus,
                len(ctrl), len(case),
                len(ctrl)+len(case),
                np.nan)
            sigall = sigall.rename(columns={'n_sig': 'n_sig_allsamples'})

            # Concatenate with the other results full results, so that each
            # combination of parameters has additional columns with the results
            # for all samples
            nsig_results = pd.merge(
                sigall[['alpha', 'n_sig_allsamples', 'taxa_level']],
                nsig_results
            )

            ## Add study label
            nsig_results['study'] = dataset
            tophits_results['study'] = dataset

            ## Add rep label
            nsig_results['rep'] = r
            tophits_results['rep'] = r

            all_nsigs.append(nsig_results)
            all_tophits.append(tophits_results)

    ## Finally, combine all datasets' results together
    all_tophits_df = pd.concat(all_tophits)
    all_nsigs_df = pd.concat(all_nsigs)

    return all_tophits_df, all_nsigs_df

p = argparse.ArgumentParser()
p.add_argument('--nreps',
    help='number of simulation reps per parameter setting')
p.add_argument('--fout-nsig',
    help='outfile with the number significant per parameter setting per rep')
p.add_argument('--fout-tophits',
    help='outfile with number of significant *top hits* per parameter '
          + 'setting per rep per total N top hits')
args = p.parse_args()

nresp = args.nreps
fout_nsig = args.fout_nsig
fout_tophits = args.fout_tophits
np.random.seed(12345)

# These values should reflect what we'd expect
# from "reasonable" clinical trials
totalNs = [10, 25, 50, 100, 150]
perc_success = [0.1, 0.25, 0.5, 0.75, 0.9]

### DEBUGGING ONLY
#totalNs = [100, 50]#, 50, 100, 150, 200]
#perc_success = [0.1, 0.25]#, 0.5, 0.75, 0.9]

# Note: different datasets represent different "expected effect sizes"
DATASETS = ['cdi_schubert', 'crc_baxter', 'ob_goodrich']

CTRLS = {'cdi_schubert': 'H',
         'crc_baxter': 'H',
         'ibd_papa': 'nonIBD',
         'ob_goodrich': 'H'}
CASES = {'cdi_schubert': 'CDI',
         'crc_baxter': 'CRC',
         'ibd_papa': 'IBD',
         'ob_goodrich': 'OB'}

################# SIMULATION #################

all_tophits_df, all_nsigs_df = run_simulation(
    n_reps, DATASETS, CTRLS, CASES, totalNs, perc_success)

## Write to files
all_nsigs_df.to_csv(fout_nsig, sep='\t', index=False)
all_tophits_df.to_csv(fout_tophits, sep='\t', index=False)
