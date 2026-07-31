from cooltools import insulation
import utils
import matplotlib.pyplot as plt
from matplotlib import colors
import bioframe
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cooltools.lib import plotting
import numpy as np
import pandas as pd
from multiprocessing import Pool

# insulation and boundary functions

# compute insulation as a list of table per cooler.
## TODO multiprocess

def _run_insulation(args):
    clr, windows, view_df = args
    return insulation(clr, windows, view_df, nproc=1)
def compute_insulation(coolers, windows, regions=None):
    if regions is not None:
        view_df = utils.vf_from_reg(regions)
    else:
        # Use bioframe to fetch the genomic features from the UCSC.
        hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
        hg38_cens = bioframe.fetch_centromeres('hg38')
        # create a view with chromosome arms using chromosome sizes and definition of centromeres
        hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

        # select only those chromosomes available in cooler
        # assume all the coolers that you process together are made with the same chromsizers
        view_df = hg38_arms[hg38_arms.chrom.isin(coolers[0].chromnames)].reset_index(drop=True)
    print(view_df)

    with Pool(processes=3) as pool:  # adjust number of processes as needed
        results = pool.map(_run_insulation, [(clr, windows, view_df) for clr in coolers])
    return results

def plot_insulation_base(clr_mat, insulation_table, off_diag, region, res, norm=None, ylim=None):
    if norm is None:
        norm = colors.LogNorm()
    f, ax = plt.subplots(figsize=(20, 10))
    im = utils.pcolormesh_45deg(ax, clr_mat, start=region[1], resolution=res, norm=norm, cmap='fall')

    ax.set_aspect(0.5)
    ax.set_ylim(0, 1000 * off_diag)  # modify later so doesnt have to take windows
    utils.format_ticks(ax, rotate=False)
    ax.xaxis.set_visible(False)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
    plt.colorbar(im, cax=cax)

    insul_region = bioframe.select(insulation_table, region)

    ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)

    ax.set_xlim(region[1], region[2])
    if ylim is not None: # jk, this doesn't work or do what i want it to do
        ax.set_ylim(ylim[0], ylim[1])
    utils.format_ticks(ins_ax, y=False, rotate=False)

    return f, ax, im, insul_region, ins_ax


# plot a single region at single resolution with insulation score track and boundaries
# plotting code modified from https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html
# clr_mat = matrix of clr at region
# region is 3-tuple
def plot_insulation_boundary(clr_mat, insulation_table, windows, region, res,
                             norm=None, off_diag=100, tstr=None, plot_weak_bounds=True, ylim=None, sd=None):
    f, ax, im, insul_region, ins_ax = plot_insulation_base(clr_mat, insulation_table, off_diag, region, res, norm, ylim=ylim)
    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

    boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])]
    weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
    strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]

    if plot_weak_bounds:
        ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
                       weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
    ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
                   strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

    ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)

    if sd is None:
        plt.show()
    else:
        plt.savefig(sd)
    if tstr is not None:
        plt.title(tstr)
    pass


# plot a single region at multiple resolutions with no boundaries
def plot_insulations(clr_mat, insulation_table, windows, region, res, norm=None, off_diag=100):
    f, ax, im, insul_region, ins_ax = plot_insulation_base(clr_mat, insulation_table, off_diag, region, res, norm)

    for res in windows:
        ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), insul_region[f'log2_insulation_score_{res}'],
                    label=f'Window {res} bp')

    plt.show()
    pass


def binary_search_for_match(row, table):
    # find best match for row in table
    return 0


def strengths_per_boundary(bound, matches, w):
    boundary_strengths = []
    for match in matches:
        # ummm binary search for match? merge_asof?
        # subset columns
        best_row = binary_search_for_match(bound, match)
        boundary_strengths.append(best_row[f'boundary_strength_{w}'])

    return boundary_strengths


# boundary strength
def boundary_matching(id_table, match_tables, match_names, window, pad=0, tol=None, merge_method='left'):
    #cols = ['chrom', 'start', 'end']
    #cols.extend([f"strength_{n}_{window}" for n in match_names])
    #all_boundaries = pd.DataFrame(columns=cols)

    # strong boundaries to match with
    all_boundaries = id_table[id_table[f'is_boundary_{window}'] == True].loc[:, 'chrom':'end']
    print(all_boundaries)
    for table, n in zip(match_tables, match_names):
        if merge_method=='outer':
            table_sub = table[table[f'is_boundary_{window}'] == True].loc[:, 'chrom':'end']
        else:
            table_sub = table[
                ['chrom', 'start', 'end', f'boundary_strength_{window}']].dropna()  # only where boundary exists
        all_boundaries = pd.merge(all_boundaries, table_sub, on=['chrom', 'start', 'end'], how=merge_method)
        all_boundaries.rename(columns={f'boundary_strength_{window}': f"strength_{n}_{window}"}, inplace=True)

    #print(all_boundaries[['strength_1_20000', 'strength_2_20000', 'strength_3_20000']])
    print(all_boundaries.columns)
    # out = id_table.apply(lambda x: strengths_per_boundary(x, match_tables, window), axis=1)
    # add out to all_boundaries

    return all_boundaries
