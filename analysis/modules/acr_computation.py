''''
functions for computing attributes about Micro-C around
ATAC-seq peaks of interest
'''

import sys
import pandas as pd
import utils
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.interpolate import make_smoothing_spline
from matplotlib.colors import LogNorm
from cooltools.api import snipping
from cooltools.lib import plotting
from multiprocessing import Pool

"""FUNCTIONS FOR READING AND FORMATTING FILES """

def read_df(f, sep='\t', header=None, columns=None):
    df = pd.read_csv(f, sep=sep, header=header)
    if columns is not None:
        assert len(columns) == len(df.columns), "Column lengths must match"
        df.columns = columns
    return df


def read_peaks(f, sep='\t', header=0):
    # set header if None
    peaks = pd.read_csv(f, sep=sep, header=header)
    if header is None:
        peaks.columns = ['chr', 'start', 'end','p_weight', 'perturbation_name']
    peaks.columns = ['chr', 'start', 'end','p_weight', 'perturbation_name']
    peaks['coord'] = peaks[['start', 'end']].mean(axis=1)
    peaks.sort_values(by=['chr', 'coord', 'perturbation_name'], inplace=True)

    return peaks[['chr', 'coord', 'start', 'end','p_weight', 'perturbation_name']]

def split_loops_into_anchors(loops, format='bedpe'):
    # infer whether the loops are in anchor or coordinate format? for now its a parameter

    if format == 'bed':
        anchor1 = loops.iloc[:, 0:1].copy()
        anchor2 = loops.iloc[:, 2:3].copy()
    elif format == 'bedpe':
        anchor1 = pd.concat([loops.iloc[:,0],loops.iloc[:,1:2].mean(axis=1)], axis=1)
        anchor2 = pd.concat([loops.iloc[:,3],loops.iloc[:,4:5].mean(axis=1)], axis=1)
    else:
        sys.exit("type not supported")

    # reformat both formats to be 'chr', 'coord'
    anchor1.columns= ['chr', 'coord']
    anchor2.columns = ['chr', 'coord']

    # later functions need sorted values (need to assert at some point)
    anchor1.sort_values(by=['chr', 'coord'], inplace=True)
    anchor2.sort_values(by=['chr', 'coord'], inplace=True)

    return anchor1, anchor2


def reformat_loops(loops, format='bedpe', loop_class=None, cols_to_append=None):
    if loop_class is not None:
        # assume class is the last column
        loops = loops[loops.iloc[:, -1]==loop_class]

    if format == 'bed':
        loops_formatted = loops.iloc[:, :4] # chr, coord, chr, coord
    elif format == 'bedpe':
        coord1 = loops.iloc[:, 1:3].mean(axis=1)
        coord2 = loops.iloc[:, 4:6].mean(axis=1)
        chr=loops.iloc[:, 0]
        loops_formatted = pd.concat([chr, coord1, chr, coord2], axis=1)
    elif format == 'classes':
        coord1 = loops.iloc[:, 1:3].mean(axis=1)
        coord2 = loops.iloc[:, 5:7].mean(axis=1)
        chr = loops.iloc[:, 0]
        loops_formatted = pd.concat([chr, coord1, chr, coord2], axis=1)
    else:
        sys.exit("type not supported")

    loops_formatted.columns = ['chr1', 'coord1', 'chr2', 'coord2']

    if cols_to_append is not None:
        print(cols_to_append)
        loops_formatted[cols_to_append] = loops[cols_to_append]
    return loops_formatted.sort_values(by=['chr1', 'coord1', 'coord2'])


def read_bedsites_flex(fname, sort=None):
    bedsites = pd.read_csv(fname, header=None, sep='\t')

    # if sort exists, then sort on that column and return the dict of dfs split by that column along with the split value
    # if sort doesn't exist, then drop all the columns after chr, start, end
    if sort == 'last':
        out = {}
        for val in bedsites[bedsites.columns[-1]].unique():
            bedsub = bedsites[bedsites.iloc[:, -1] == val]
            pname = bedsub.iloc[:, -1]
            bedsub = bedsub.iloc[:, :3]

            bedsub = pd.concat([bedsub, pname], axis=1)
            bedsub.columns = ['chr', 'coord1', 'coord2', 'perturbation_name']
            bedsub['coord'] = bedsub[['coord1', 'coord2']].mean(axis=1)

            out[val] = bedsub.sort_values(by=['chr','coord']).reset_index(drop=True)
            #print(bedsub.head())
    else:
        out = bedsites.iloc[:, :3].reset_index(drop=True)
        out.columns = ['chr', 'coord1', 'coord2']
        out['coord'] = out[['coord1', 'coord2']].mean()
        out = out.sort_values(by=['chr', 'coord']).reset_index(drop=True)
        print(out.head())
    return out

""" COMPUTATION FUNCTIONS """
def find_within_distance(peak, anchor, dist):
    # find all anchors within dist of the peak coord
    anchors_ind = anchor[abs(anchor['coord'] - peak['coord']) <= dist]
    # add
    return anchors_ind

def loops_within_distance(peak, loops, dist):

    loops_dist = loops[(abs(loops['coord1']-peak['coord']) <= dist) &
                       (abs(loops['coord2']-peak['coord']) <= dist)]
    return loops_dist

def collapse_anchors(anchors, pad):
    # align the anchors so that they are probably unique
    return new_anchors

def define_elements(peak, anchors):
    peak_coord = peak['coord']
    labels = anchors.astype('object').copy()

    labels[:] = 'NONE'
    labels[anchors < peak_coord] = 'LEFT'
    labels[anchors > peak_coord] = 'RIGHT'
    # If anchor == peak.coord, remains 'NONE'
    return labels

def define_loops(peak, anchor1, anchor2, loops):
    loc1 = define_elements(peak, anchor1)
    loc2 = define_elements(peak, anchor2)

    loops = loops.copy()
    loops['loc1'] = loc1
    loops['loc2'] = loc2

    def classify(loc1, loc2):
        if 'NONE' in (loc1, loc2):
            return 'NONE'
        elif loc1 == loc2:
            return 'same-side'
        else:
            return 'crosser'

    loops['loop_type'] = [classify(l1, l2) for l1, l2 in zip(loc1, loc2)]

    return loops

def define_loops_all(peaks, loops, dist=100000, count_col=None, agg_fun='mean'):
    all_annotated = []
    ss_counts = []
    cross_counts = []
    ratios = []

    for i, peak in peaks.iterrows():
        loops_chr = loops[loops['chr1'] == peak['chr']]
        loops_sub = loops_within_distance_fast(peak, loops_chr, dist)
        annotated = define_loops(peak, loops_sub['coord1'], loops_sub['coord2'], loops_sub)

        # Add peak info
        annotated['peak_chr'] = peak['chr']
        annotated['peak_coord'] = peak['coord']

        if count_col is None:
            value_counts = annotated['loop_type'].value_counts()
            crosser_count = value_counts.get('crosser', 0)
            same_side_count = value_counts.get('same-side', 0)
        else:
            loop_grouped = annotated.groupby('loop_type')[count_col].agg(agg_fun)
            crosser_count = loop_grouped.get('crosser', default=0)
            same_side_count = loop_grouped.get('same-side', default=0)

        all_annotated.append(annotated)
        ratio = crosser_count / same_side_count if same_side_count > 0 else np.nan
        ss_counts.append(same_side_count)
        cross_counts.append(crosser_count)
        ratios.append(ratio)
    return pd.concat(all_annotated, ignore_index=True), ss_counts, cross_counts, ratios

def loops_within_distance_fast(peak, loops, dist, excl=2000):
    coord = peak['coord']
    coord1 = loops['coord1'].values
    #coord2 = loops['coord2'].values

    # Binary search to limit coord1 search space
    left = np.searchsorted(coord1, coord - dist, side='left')
    right = np.searchsorted(coord1, coord + dist, side='right')

    # Slice the narrowed region
    loops_sub = loops.iloc[left:right]

    # Now filter by coord2
    abs_diff1 = abs(loops_sub['coord1'] - coord)
    abs_diff2 = abs(loops_sub['coord2'] - coord)

    mask = (abs_diff1 >= excl) & (abs_diff1 <= dist) & \
           (abs_diff2 >= excl) & (abs_diff2 <= dist)

    return loops_sub[mask]

def compute_triangle_score(snip, verbose, nan_thresh=0.8):
    n = snip.shape[0]
    if n == 0 or snip.shape[0] != snip.shape[1]:
        return np.nan

    if (np.isnan(snip).sum() / snip.size) > nan_thresh:
        return np.nan

    triu = np.triu(snip)
    mid = n // 2 # adjust for the 2 pix excl in the snipper

    # Top square: rows 0:mid, cols mid:n
    top_square = triu[2:mid, mid:-2]
    top_left_triangle = triu[2:mid, 2:mid]
    bottom_right_triangle = triu[mid:-2, mid:-2]

    numerator = np.nansum(top_square)
    denominator = np.nansum(top_left_triangle) + np.nansum(bottom_right_triangle)

    if verbose:
        plt_regions = snip.copy()
        plt_regions[2:mid, mid:-2] = 3
        plt_regions[2:mid, 2:mid] = 1
        plt_regions[mid:-2, mid:-2] = 2

        top_square_counts = np.count_nonzero(~np.isnan(top_square))
        bottom_right_counts = np.count_nonzero(~np.isnan(bottom_right_triangle))
        top_left_counts = np.count_nonzero(~np.isnan(top_left_triangle))
        print(top_square_counts, bottom_right_counts, top_left_counts)
        f, ax = plt.subplots()
        ax.matshow(snip)
        ax.matshow(plt_regions, alpha=0.2)

        ax.axvline(mid, color='red', linestyle='dotted', linewidth=2)
        ax.axhline(mid, color='red', linestyle='dotted', linewidth=2)
        ax.axvline(2, color='red', linestyle='dotted', linewidth=2)
        ax.axhline(n-2, color='red', linestyle='dotted', linewidth=2)
        plt.show()

    if denominator == 0:
        return np.nan

    return numerator / denominator

"""
Return all the triangle scores for a df of peaks
Redundant to below; there's no reason not to just return the peaks as well 
If verbose, plot the region with the triangles highlighted to check biases 
"""
## TODO test all instances of this and merge with the peak-inclusive version of this function
def snip_region(clr, PS_df, peaks, dist, test_range=None, verbose=False):
    snipper = snipping.ObsExpSnipper(clr, PS_df)
    scores = {p: [] for p in peaks['perturbation_name'].unique()}
    dist = int(dist)

    grouped_peaks = peaks.groupby('chr')
    for chrom, peaks_chrom in grouped_peaks:
        if chrom in clr.chromnames:
            if verbose: print(f"loading chromosome {chrom}")
            mat_chrom = snipper.select(chrom, chrom)
            if verbose: print(f"Finished loading chromosome {chrom}")

            # Precompute constant region values
            region1 = chrom
            region2 = chrom

            for i, peak in peaks_chrom.reset_index(drop=True).iterrows():
                c_lo = int(peak['coord'] - dist)
                c_hi = int(peak['coord'] + dist)
                tup = (c_lo, c_hi, c_lo, c_hi)
                snip = snipper.snip(mat_chrom, region1, region2, tup)

                if i in test_range:
                    test = True

                score = compute_triangle_score(snip, test)
                if verbose and i % 100 == 0: print(f"On peak {i} out of {len(peaks_chrom)}, score {score}...")
                scores[peak['perturbation_name']].append(score) # doesn't save peak information
                test = False
    return scores

"""
Same as above, but also returns the peak coordinate
Return all the triangle scores for a df of peaks
If verbose, plot the region with the triangles highlighted to check biases 
"""
def snip_region_withpeaks(clr, PS_df, peaks, dist, test_range=None, verbose=False):
    snipper = snipping.ObsExpSnipper(clr, PS_df)
    peak_scores = []
    dist = int(dist)

    grouped_peaks = peaks.groupby('chr')
    for chrom, peaks_chrom in grouped_peaks:
        if chrom in clr.chromnames:
            if verbose: print(f"loading chromosome {chrom}")
            mat_chrom = snipper.select(chrom, chrom)
            if verbose: print(f"Finished loading chromosome {chrom}")

            # Precompute constant region values
            region1 = chrom
            region2 = chrom

            for i, peak in peaks_chrom.reset_index(drop=True).iterrows():
                c_lo = int(peak['coord'] - dist)
                c_hi = int(peak['coord'] + dist)
                tup = (c_lo, c_hi, c_lo, c_hi)
                snip = snipper.snip(mat_chrom, region1, region2, tup)

                if i in test_range:
                    test = True

                score = compute_triangle_score(snip, test)
                #curr_score = [chrom, peak['coord'], peak['perturbation_name'], score]
                peak_scores.append({
                        'chr': chrom,
                        'coord': peak['coord'],
                        'perturbation_name': peak['perturbation_name'],
                        'score': score
                    })
                if verbose and i % 100 == 0: print(score)


                test = False
    peak_scores_df = pd.DataFrame(peak_scores)
    return peak_scores_df

"""
Get the loop scores in the crosser and same side regions
Loops are formatted according to reformat_loops() (must be sorted for fast search to work!) 
"""
## TODO could stand to be paralellized. but i dont like the idea of the snipper living in multiple processes
## TODO also attribute partial scores by class
def loop_scores_snip(clr, PS_df, peaks, loops, dist, quantsize=5000, verbose=False):
    snipper = snipping.ObsExpSnipper(clr, PS_df)
    peak_scores = []
    dist = int(dist)

    grouped_peaks = peaks.groupby('chr')
    loops_by_chr = {chrom: df for chrom, df in loops.groupby('chr1')}

    for chrom, peaks_chrom in grouped_peaks:
        if chrom in clr.chromnames:
            if verbose: print(f"loading chromosome {chrom}")
            mat_chrom = snipper.select(chrom, chrom)
            if verbose: print(f"Finished loading chromosome {chrom}")

            # Precompute constant region values
            region1 = chrom
            region2 = chrom

            loops_chr = loops_by_chr.get(chrom, pd.DataFrame())
            for i, peak in peaks_chrom.reset_index(drop=True).iterrows():

                # get the triangle score (can omit i guess and just merge with existing score df later)
                c_lo = int(peak['coord'] - dist)
                c_hi = int(peak['coord'] + dist)
                tup = (c_lo, c_hi, c_lo, c_hi)

                snip = snipper.snip(mat_chrom, region1, region2, tup)
                score = compute_triangle_score(snip, False)

                if len(loops_chr)==0:
                    continue

                loops_sub = loops_within_distance_fast(peak, loops_chr, dist)
                annotated = define_loops(peak, loops_sub['coord1'], loops_sub['coord2'], loops_sub)

                loops_ss = annotated[annotated['loop_type'] == 'same-side']
                loops_crosser = annotated[annotated['loop_type'] == 'crosser']

                crosser_loop_score_list = []
                ss_loop_score_list = []

                crosser_loop_class_list = []
                ss_loop_class_list = []

                crosser_loop_coord_list = []
                ss_loop_coord_list = []


                for loop_type, subset in [('same-side', loops_ss), ('crosser', loops_crosser)]:
                    scores, classes, coords = [], [], []

                    for i, loop_row in subset.iterrows():
                        tup = (int(loop_row.coord1 - quantsize), int(loop_row.coord1 + quantsize),
                               int(loop_row.coord2 - quantsize), int(loop_row.coord2 + quantsize))
                        snip = snipper.snip(mat_chrom, region1, region2, tup)

                        scores.append(np.mean(snip))
                        classes.append(loop_row['pair'])
                        coords.append((loop_row.coord1, loop_row.coord2))

                    if loop_type == 'same-side':
                        ss_loop_score_list, ss_loop_class_list, ss_loop_coord_list = scores, classes, coords
                    elif loop_type == 'crosser':
                        crosser_loop_score_list, crosser_loop_class_list, crosser_loop_coord_list = scores, classes, coords
                    else:
                        print("Somehow, the wrong type of variable was accessed")
                peak_scores.append({
                        'chr': chrom,
                        'coord': peak['coord'],
                        'perturbation_name': peak['perturbation_name'],
                        'triangle_score': score,
                        'crosser_loop_scores': crosser_loop_score_list,
                        'crosser_loop_classes': crosser_loop_class_list,
                        'crosser_loop_coords': crosser_loop_coord_list,
                        'ss_loop_scores': ss_loop_score_list,
                        'ss_loop_classes': ss_loop_class_list,
                        'ss_loop_coords': ss_loop_coord_list,
                        'agg_crosser_score': np.sum(crosser_loop_score_list),
                        'agg_ss_score': np.sum(ss_loop_score_list)
                    })

    peak_scores_df = pd.DataFrame(peak_scores)
    return peak_scores_df

"""
makes the diagonal extension plots
"""
## TODO: i do this first snippet a bunch of times, can i make it all one main fn. that routes to diff. tasks
def intensity_by_distance(clr, PS_df, peaks, dist, plotting=False):
    snipper = snipping.ObsExpSnipper(clr, PS_df)
    dist = int(dist)
    peak_scores = []
    grouped_peaks = peaks.groupby('chr')

    for chrom, peaks_chrom in grouped_peaks:
        if chrom in clr.chromnames:
            print(f"loading {chrom}...")
            mat_chrom = snipper.select(chrom, chrom)
            region1 = chrom
            region2 = chrom

            for i, peak in peaks_chrom.reset_index(drop=True).iterrows():
                c_lo = int(peak['coord'] - dist)
                c_hi = int(peak['coord'] + dist)
                tup = (c_lo, c_hi, c_lo, c_hi)

                snip = snipper.snip(mat_chrom, region1, region2, tup)
                score = compute_triangle_score(snip, False)

                n = snip.shape[0]

                mid = n // 2

                # compute valid range lengths
                steps_up = mid + 1
                steps_right = n - mid
                L = min(steps_up, steps_right)

                # equal-length indices for the 45° anti-diagonal (upper half)
                rows = np.arange(mid, mid - L, -1)
                cols = np.arange(mid, mid + L)

                ortho_stripe = snip[rows, cols]

                # nan-safe max extraction
                try:
                    max_ortho_ind = np.nanargmax(ortho_stripe)

                    peak_scores.append({
                        'chr': chrom,
                        'coord': peak['coord'],
                        'perturbation_name': peak['perturbation_name'],
                        'triangle_score': score,
                        'ortho_max_ind': max_ortho_ind,
                        'ortho_max_val': np.nanmax(ortho_stripe)
                    })

                    if plotting:
                        snip_plot = snip.copy()
                        # snip_plot[rows, cols] = 0

                        f, ax = plt.subplots()
                        ax.matshow(snip_plot)
                        plt.show()

                        print(max_ortho_ind)
                        print(np.nanmax(ortho_stripe))
                        f, ax = plt.subplots()
                        ax.plot(ortho_stripe)
                        ax.scatter(max_ortho_ind, np.nanmax(ortho_stripe))
                        plt.show()

                except ValueError:
                    peak_scores.append({
                        'chr': chrom,
                        'coord': peak['coord'],
                        'perturbation_name': peak['perturbation_name'],
                        'triangle_score': score,
                        'ortho_max_ind': np.nan,
                        'ortho_max_val': np.nan
                    })
                # now get the distance-wise intensity of the snip
                # max_diag = snip.shape[0] - 1
                # for offset in range(-2, max_diag, -1):
                #      diag_vals = np.diagonal(snip, offset=offset)
                #      diag_vals = np.pad(diag_vals, int((len(np.diagonal(snip, 0))-len(diag_vals))/2)) # make all vectors the same length
                # it should be symmetrical, but maybe check

    peak_scores_df = pd.DataFrame(peak_scores)
    return peak_scores_df

""" PLOTTING """

""" 
Make metaplot of loops around a peak 
Basically a heatmap of where the loops are
"""
def make_metaplot(peaks, loops, dist=100000, binsize=2000, pad_pix=0,
                  tstr=None, sd=None, excl=2000, heatmap=None, count_col=None, agg_fun='sum', vmin=None, vmax=None):
    #print(peaks)
    all_loop_pix = []
    avals = []
    # handle extra column across which to plot
    if count_col is not None:
        assert count_col in loops.columns, "Counting column is not in columns"
        if loops[count_col].isnull().any():
            loops = loops.dropna(count_col, axis=1)

    # create df of loops within distance and their values
    for i, peak in peaks.iterrows():
        loops_chr = loops[loops['chr1'] == peak['chr']]
        loops_sub = loops_within_distance_fast(peak, loops_chr, dist, excl=excl)
        loop_pix = utils.bp_to_pix(loops_sub[['coord1', 'coord2']],
                                  peak['coord']-dist-pad_pix*binsize, binsize=binsize)

        # i can compute it directly, but just in case
        peak_pix = utils.bp_to_pix(peak['coord'], peak['coord']-dist-pad_pix*binsize, binsize=binsize)

        if count_col is not None:
            avals.append(loops_sub[count_col])
            loop_pix = pd.concat([loop_pix, loops_sub[count_col]], axis=1)

        all_loop_pix.append(loop_pix)

    all_loop_pix_df = pd.concat(all_loop_pix, ignore_index=True)

    # make it heatmap style where its sum or average loop characteristic
    # if each pixel could only have one value per peak (i.e. no doublecounting loops per pixel so binsize
    # should always be less than loopcall binsize), then sum over all coordinates

    plt.figure()
    if heatmap is None:
        if count_col is not None:
            alphas = pd.concat(avals, ignore_index=True)
            alphas = alphas + np.abs(alphas.min())
            alphas = np.array(alphas.divide(alphas.max()).divide(2))
        else:
            alphas = 0.2

        plt.matshow(np.zeros((2 * dist // binsize + pad_pix, 2 * dist // binsize + pad_pix)), cmap='Grays')
        plt.scatter(all_loop_pix_df.iloc[:, 0], all_loop_pix_df.iloc[:, 1], alpha=alphas, s=1)
        plt.scatter(peak_pix, peak_pix, color='k', s=5)

    else:
        if count_col is not None:
            loop_pix_counts = all_loop_pix_df.groupby(['coord1', 'coord2'], as_index=False)[count_col].agg(agg_fun)
            col_to_agg = count_col
        else:
            loop_pix_counts = all_loop_pix_df.groupby(['coord1', 'coord2'], as_index=False).value_counts()
            col_to_agg = 'count'

        print(loop_pix_counts)
        mat = np.zeros((2 * dist // binsize + pad_pix, 2 * dist // binsize + pad_pix))
        for i, row in loop_pix_counts.iterrows():
            mat[int(row['coord2']), int(row['coord1'])] = row[col_to_agg]
        metaplot = sns.heatmap(mat, cmap='fall', vmin=vmin, vmax=vmax)
        metaplot.tick_params(left=False, bottom=False)  ## other options are right and top

        # plot the intensity profile of the matrix for each diagonal, with a legend describing the diagonal
        # max_diag = mat.shape[0] - 1
        # max_diag = -6
        # plt.figure(figsize=(16,8))
        # for offset in range(-2, max_diag, -1):
        #     diag_vals = np.diagonal(mat, offset=offset)
        #     diag_vals = np.pad(diag_vals, int((len(np.diagonal(mat, 0))-len(diag_vals))/2))
        #     spl = make_smoothing_spline(np.arange(0, len(diag_vals)), diag_vals, lam=1)
        #     plt.scatter(np.linspace(0, len(diag_vals), len(diag_vals)), diag_vals, label=f'diagonal {offset}')
        #     plt.plot(np.arange(0, len(diag_vals)), spl(np.arange(0, len(diag_vals))))
        # plt.xlabel('Index along diagonal')
        # plt.ylabel('Intensity')
        # plt.title('Intensity Profiles by Diagonal')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        # plt.tight_layout()
        # plt.show()

    if tstr is not None:
        plt.title(f"{tstr}, num loops: {len(all_loop_pix_df)}")

    if sd is not None:
        plt.savefig(sd)
    else:
        plt.show()

    if heatmap is not None:
        vmax = metaplot.collections[0].norm.vmax
        return vmax
    else:
        return 0

