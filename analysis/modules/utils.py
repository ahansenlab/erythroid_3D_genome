import re
import pandas as pd
import numpy as np
from setup_logger import logger
from cooltools.lib import numutils
import cooler
import bioframe

import matplotlib.pyplot as plt
from matplotlib import colors
import os
import cooltools
import sys

from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')
from matplotlib.patches import Patch

# given a list of regions in ucsc coordinates, return the viewframes (cooltools compatible)
def vf_from_reg(regions):
    # if it's a dictionary, use the key of each entry as the name
    if isinstance(regions, dict):
        names = regions.keys()
        chroms, starts, ends = zip(*(get_coords(r) for r in regions.values()))
    # else, use the chromosome name
    elif isinstance(regions, list):
        chroms, starts, ends = zip(*(get_coords(r) for r in regions))
        names = chroms
    else:
        logger.error("Invalid regions provided")
        return 1

    df = pd.DataFrame.from_dict({'chrom': chroms, 'start': starts, 'end':ends, 'name': names})
    df['sort'] = df['chrom'].str.extract('(\d+)', expand=False).astype(int)
    return df.sort_values(['sort', 'start', 'end']).drop('sort', axis=1).reset_index(drop=True)

def get_hg38_arms(clr=None, excl_ym=True):
    hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
    hg38_cens = bioframe.fetch_centromeres('hg38')
    hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

    if clr is not None:
        hg38_arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

    if excl_ym:
        hg38_arms = hg38_arms[~hg38_arms["chrom"].isin(['chrY','chrM'])]

    return hg38_arms

def get_locus_from_coord(row, flank):
    return f"{row['chr']}:{int(row['coord'])-flank}-{int(row['coord'])+flank}"

# helper to get the chr and region coords from a ucsc string, with handling for chrX/Y
def get_coords(locus_str):
    m = re.search(r'(chr\w+):(\d+)-(\d+)', locus_str)
    chrom = m.group(1)
    start_coord = int(m.group(2))
    end_coord = int(m.group(3))
    try:
        chrom = int(chrom)
    except:
        return chrom, start_coord, end_coord
    return chrom, start_coord, end_coord

# get a matrix (opt: balanced, normalized) from a cooler
def get_matr(clr, region, balance=True, oe=True):
    # balance
    clr_mat = clr.matrix(balance=balance).fetch(region)
    clr_mat[np.isnan(clr_mat)] = 0

    # normalize
    if oe:
        return numutils.observed_over_expected(clr_mat)[0]
    else:
        return clr_mat

def get_clr(clr_path, res):
    return cooler.Cooler(f"{clr_path}::/resolutions/{res}")

def write_get_ps(zipped_clrs, write=True, proc=8):
    PS_dfs = []

    for p, clr in zipped_clrs:
        if not os.path.exists(p):
            if write:
                print(f"making expected at {p}...")
                PS_df = cooltools.expected_cis(clr=clr, smooth=True, aggregate_smoothed=True, nproc=proc)
                PS_df.to_csv(p, sep='\t', index=False)
            else:
                sys.exit("P(s) file does not exist")
        else:
            PS_df = pd.read_csv(p, sep='\t')

        PS_dfs.append(PS_df)

    return PS_dfs

def bp_to_pix(coords_bp, st, binsize):
    return (coords_bp - st) // binsize

def generate_all_chroms(y=False, m=False):
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    if y:
        chroms.append('chrY')
    if m:
        chroms.append('chrM')
    return chroms

def clean_chroms(df, chrom_col='chr'):
    return df[df[chrom_col].isin(generate_all_chroms())]



# from the fracshift code
def plot_simple(matr, pt=None, ts= None, sd=None):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    im2 = ax.matshow(matr[:], norm=colors.LogNorm(), cmap='fall')
    plt.colorbar(im2)
    if pt is not None:
        plt.scatter(pt[0], pt[1])
        # plt.scatter(0,0)
        # plt.scatter(1, 2)
    if ts is not None:
        plt.title(ts)

    if sd is not None:
        plt.savefig(sd)
    else:
        plt.show()
    pass

#### COOLTOOLS PLOTTING UTILS #####
# Functions to help with plotting
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    import matplotlib.pyplot as plt
    return im


def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


def custom_legend_list(data1, colors):
    # Create custom legend handles
    means = {}
    for x in data1.keys():
        means[x] = np.nanmedian(data1[x])

    legend_handles = [
        Patch(color=colors[x], label=f'{x}: mean={m:.3f}')
        for x, m in means.items()
    ]

    # Add custom legend
    return legend_handles


def get_locus_from_coord(row, flank):
    return f"{row['chr']}:{int(row['coord']) - flank}-{int(row['coord']) + flank}"


def plot_peak(row, clr, flank=200_000, res=2000, fs=5, toplot=None, score_col=None):
    fig, ax = plt.subplots(figsize=(fs, fs))

    locus = get_locus_from_coord(row, flank)

    try:
        mat = get_matr(clr, locus, oe=False)
    except ValueError:
        return fig

    ax.matshow(mat, norm=colors.LogNorm(), cmap='fall')
    ax.axvline((flank // res), ls='--', color='k', alpha=0.4)
    ax.axhline((flank // res), ls='--', color='k', alpha=0.4)

    if toplot is not None:
        toplot_sub = toplot[(toplot['chr'] == row['chr']) & (toplot['start'] > (row['coord'] - flank)) & (
                    toplot['end'] < (row['coord'] + flank))]
        print(toplot_sub)
        print(locus)
        for i, plotrow in toplot_sub.iterrows():
            coord = (plotrow['start'] + plotrow['end']) / 2
            ax.axvline(((coord - (row['coord'] - flank)) // res), ls='--', color='b', alpha=0.3)
            ax.axhline(((coord - (row['coord'] - flank)) // res), ls='--', color='b', alpha=0.3)

    if score_col is not None:
        plt.title(f"{locus}: score {row[score_col]:.5f}")
    return fig


