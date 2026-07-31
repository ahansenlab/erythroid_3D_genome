import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from tqdm import tqdm
import itertools
from hicreppy.hicrep import h_train, genome_scc
import os
import warnings
from multiprocessing import Pool

## these two multiprocessing functions were written by domenic
def process_with_progress(arguments, threads):
    results = []
    with Pool(processes=threads) as pool:
        # tqdm doesn't really work here, but that's okay
        for result in tqdm(pool.imap_unordered(scccompute, arguments), total=len(arguments)):
            results.append(result)
    return results


def scccompute(argument):
    # For RCMC, stops this from printing an error message a million times when it encounters an empty diagonal
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                message="An input array is constant; the correlation coefficient is not defined.")

        mcool1 = argument[1][argument[0][0]]
        mcool2 = argument[1][argument[0][1]]
        # Run comparison and return scc value and which pair it was
        scc = genome_scc(mcool1, mcool2, argument[2], argument[3], argument[4], argument[5])
        pair = argument[0]
    return scc, pair

# Replicate reproducibility, P(s) curves and derivatives, pairtools stats mining

def compute_repro(coolers, whitelist, param_dict = None):
    if param_dict is None:
        param_dict = {'max_dist': 5000000, 'h_max': 10}

    # whitelist = ["chr3", "chr9", "chr14"]  # can pick a random sampling of chromosomes for speed
    h_train_idxs = [0, 1]
    optimal_h = h_train(coolers[h_train_idxs[0]], coolers[h_train_idxs[1]], param_dict['max_dist'],
                        param_dict['h_max'], whitelist)

    # make scc
    cooler_sums = [cooler_file.info["sum"] for cooler_file in coolers]
    downsampling_value = min(cooler_sums)

    k = len(coolers)
    scc_mat = np.zeros((k, k))
    for i, j in tqdm(itertools.product(range(k), range(k))):
        if i == j:
            scc_mat[i, j] = 1.0
        else:
            scc = genome_scc(coolers[i], coolers[j], param_dict['max_dist'], optimal_h, downsampling_value, whitelist)
            scc_mat[i, j] = scc
            scc_mat[j, i] = scc

    return scc_mat

def plot_repro(scc_mat, sd=None):
    sns.clustermap(scc_mat, cmap=sns.color_palette("RdBu_r", as_cmap=True), annot=True)
    if sd is not None:
        plt.savefig(sd, format='svg')
    else:
        plt.show()

# read all the stats files in the folder, grab the metrics
# look for the ones that are the same except for the rep and sum the stats for those in the output table
def summary_metrics(stat_dir, name_delim='.', kw_delim='_', kw_end=-1, kw_names=None, group_kws=None, sd=None, filt_names=None):
    stat_names = [file for file in os.listdir(stat_dir) if 'stats' in file]

    if filt_names is not None:
        stat_names = [file for file in stat_names if filt_names in file]
    if kw_end == None:
        keywords = [file.split(kw_delim)[:] for file in stat_names]
    else:
        keywords = [file.split(kw_delim)[:kw_end] for file in stat_names] # assumes consistent filename splits
    print(keywords)
    cols = ['name']
    if kw_names is not None:
        if len(kw_names) == len(keywords[0]):
            cols.extend(kw_names)
        else:
            print('Wrong number of keyword names provided')
            kw_names = [f"col_{i}" for i in range(len(keywords[0]))]
            cols.extend(kw_names)
    else:
        kw_names = [f"col_{i}" for i in range(len(keywords[0]))]
        cols.extend(kw_names)

    cols.extend(['total', 'total_mapped', 'total_unmapped','unique_mapped', 'dups', 'dup_frac',
                 'cis','trans','cis_frac','trans_frac', 'cis_1kb+', 'cis_20kb+'])

    summary_table = pd.DataFrame(columns=cols)
    summary_table['name'] = [file.split(name_delim)[0] for file in stat_names]
    for i, kw in enumerate(kw_names):
        summary_table[kw] = [kws[i] for kws in keywords]

    metrics_all = dict(zip(stat_names, [grab_metrics(os.path.join(stat_dir, file), 'total', 'total_mapped',
                                                     'total_unmapped','total_nodups', 'total_dups',
                                                     'cis','trans', 'cis_1kb+', 'cis_20kb+') for file in stat_names]))

    summary_table['total'] = [metrics_all[file]['total'] for file in stat_names]
    summary_table['total_mapped'] = [metrics_all[file]['total_mapped'] for file in stat_names]
    summary_table['total_unmapped'] = [metrics_all[file]['total_unmapped'] for file in stat_names]
    summary_table['unique_mapped'] = [metrics_all[file]['total_nodups'] for file in stat_names]
    summary_table['dups'] = [metrics_all[file]['total_dups'] for file in stat_names]
    summary_table['dup_frac'] = np.array(summary_table['dups'])/np.array(summary_table['total_mapped'])
    summary_table['cis'] = [metrics_all[file]['cis'] for file in stat_names]
    summary_table['trans'] = [metrics_all[file]['trans'] for file in stat_names]
    summary_table['cis_frac'] = np.array(summary_table['cis'])/np.array(summary_table['unique_mapped'])
    summary_table['trans_frac'] = np.array(summary_table['trans'])/np.array(summary_table['unique_mapped'])
    summary_table['cis_1kb+'] = [metrics_all[file]['cis_1kb+'] for file in stat_names]
    summary_table['cis_20kb+'] = [metrics_all[file]['cis_20kb+'] for file in stat_names]
    #summary_table['1kb+_fraction'] = np.array(summary_table['cis_1kb+'])/np.array(summary_table['unique_mapped'])
    #summary_table['20kb+_fraction'] = np.array(summary_table['cis_20kb+']) / np.array(summary_table['unique_mapped'])
    print(summary_table.columns)
    if group_kws is not None:
        if kw_names is None: #or ~set(group_kws).issubset(set(kw_names)):
            print('Error in matching keywords and groups')
            grouped_df = None
        else:
            grouped_df = summary_table.groupby(by=group_kws[0]).sum()
    else:
        grouped_df = None

    if sd is not None:
        pd.concat([summary_table, grouped_df], axis=1).to_cvs(sd, sep='\t')

    return summary_table, grouped_df

def grab_metrics(statfile, *args):
    with open(statfile, 'r') as f:
        stats = [line.rsplit('\t', 1) for line in f]
        stats_dict = {key: int(value)/1_000_000_000 for key, value in stats if key in args}

    return stats_dict