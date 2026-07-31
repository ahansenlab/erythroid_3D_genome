## functions for a lot for the matchmaker analyses, to keep them consistent
## a lot of custom formatting so they don't really work for anything else

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from functools import reduce
import seaborn as sns
import ast
import os
from sklearn import cluster, metrics
import utils
from matplotlib.colors import LogNorm
def read_match_file(basedir, phase, nametag, dist, list_cols,
                    score_type="loopscores", do_filter=True, difference_threshold=5000):
    
    if score_type == 'loopscores':
        curr = pd.read_csv(os.path.join(basedir, f"{phase}_{nametag}_loopscores_{dist}_withpeaks.tsv"),
                           sep='\t', index_col=0)
    elif score_type == 'trianglescores':
        curr = pd.read_csv(os.path.join(basedir, f"{phase}_{nametag}_{dist}_withpeaks.tsv"),
                           sep='\t', index_col=0).rename({'score':'triangle_score'}, axis=1)
    else:
        print("Don't have that score type")
        return np.nan

    curr[list_cols] = curr[list_cols].replace('nan', 'None', regex=True).map(ast.literal_eval)
    try:
        curr['mean_crosser_score'] = curr['crosser_loop_scores'].apply(
            lambda x: np.mean(vals) if (vals := [v for v in x if v is not None]) else np.nan)
        curr['mean_ss_score'] = curr['ss_loop_scores'].apply(
            lambda x: np.mean(vals) if (vals := [v for v in x if v is not None]) else np.nan)
    except:
        print("column not found")
    if do_filter:
        curr = filter_by_distance(curr, difference_threshold)
    return curr


def load_matchmaker_scores_phase(basedir, phases, nametag, list_cols,
                                 dist=100000, do_filter=True, difference_threshold=5000, score_type='loopscores',
                                 triangle_score_cols=('triangle_score_P1', 'triangle_score_P2', 'triangle_score_P3')):
    match_scores_list = []
    for i, phase in enumerate(phases):
        curr = read_match_file(basedir, phase, nametag, dist, list_cols,
                               do_filter=False, score_type=score_type)
        match_scores_list.append(curr.rename(columns={c: c + f'_P{i + 1}' for c in curr.columns[3:]}))

    if "perturbation_name" in match_scores_list[0].columns:
        match_scores = reduce(
            lambda left, right: pd.merge(left, right, on=['chr', 'coord', 'perturbation_name'], how='outer'),
            match_scores_list)
    else:
        match_scores = reduce(lambda left, right: pd.merge(left, right, on=['chr', 'coord'], how='outer'),
                              match_scores_list)

    match_scores["locus"] = match_scores.apply(lambda row: f"{row['chr']}:{int(row['coord'])}",
                                               axis=1)
    match_scores.dropna(subset=triangle_score_cols, inplace=True)
    if do_filter:
        match_scores = filter_by_distance(match_scores, difference_threshold, score_to_filt=triangle_score_cols[-1])
    match_scores_dedup = match_scores.drop_duplicates(subset=['chr', 'coord'])
    match_scores_dedup.set_index("locus", inplace=True)

    return match_scores_dedup


def filter_by_distance(df, difference_threshold, score_to_filt='triangle_score_P3'):
    df = df.sort_values(by=['chr', 'coord']).copy()
    df['difference'] = df['coord'].diff().abs()

    # start a new group whenever difference > threshold 
    new_group = (
            (df['difference'] > difference_threshold) |
            (df['chr'] != df['chr'].shift())
    )
    df['group'] = new_group.cumsum()

    df_filtered = (
        df.loc[df.groupby('group')[score_to_filt].idxmax()]
        .drop(columns=['difference', 'group'])
        .sort_index()
    )

    return df_filtered


def generate_score_heatmap(df, ctr=1):
    df_dropna = df.dropna()
    plt.figure(figsize=(10, 10))
    sns.clustermap(df_dropna, cmap='coolwarm', center=ctr)
    plt.show()

    pass

def plot_loops_vs_score(df, loop_col, score_col='triangle_score',
                        plot_color='#CC3311', tstr='Loops vs. Anti-Insulation Score', fax=None,
                        xlim=(0,3), ylim=None):
    if fax is None:
        f, ax = plt.subplots()
    else:
        f, ax = fax

    ax.scatter(df[score_col],
               df[loop_col], alpha=0.1, color=plot_color, rasterized=True)

    ax.axhline(1, color='k', ls='--')
    ax.axvline(1, color='k', ls='--')

    plt.xlabel('Anti-insulation score')
    plt.ylabel('Aggregate loop score')

    plt.title(tstr)
    plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    return f

def do_kmeans(df, nclus, rand=None):
    clus = cluster.KMeans(n_clusters=nclus, random_state=rand)
    cluster_labels = clus.fit_predict(df)

    X_new = clus.transform(df)
    centroids = clus.cluster_centers_

    f = lambda xrow: np.where(xrow == min(xrow))[0][0]
    # find minimum distance for each point and assign it cluster label
    # X_new: ndarray of shape (n_samples, n_clusters)
    grp = np.array(np.apply_along_axis(f, axis=1, arr=X_new))

    return grp

def plot_clusters(df, grp=None, writepath=None):
    if grp is None:
        sns.heatmap(df)
        plt.show()
    else:
        cm = dict(zip(np.unique(grp), sns.color_palette('hls', n_colors=len(np.unique(grp)))))
        df['g'] = grp
        df.sort_values(by='g', inplace=True)
        print(df)
        sns.clustermap(df.drop(columns='g', axis=1), row_colors=[cm[g] for g in df['g']],
                       row_cluster=False, col_cluster=False, cmap='coolwarm', center=1)

        if writepath is not None:
            plt.savefig(f"{writepath}.svg", dpi=300, format='svg')
        else:
            plt.show()

    return df

def merged_attrs(merged_df, test_phase, control_phase):
    merged_df['triangle_fc'] = merged_df[f'triangle_score_{test_phase}'] / merged_df[f'triangle_score_{control_phase}']
    merged_df['ss_fc'] = merged_df[f'agg_ss_score_{test_phase}'] / merged_df[f'agg_ss_score_{control_phase}']
    merged_df['crosser_fc'] = merged_df[f'agg_crosser_score_{test_phase}'] / merged_df[
        f'agg_crosser_score_{control_phase}']
    merged_df['ss_crosser_ratio_control'] = merged_df[f'agg_crosser_score_{control_phase}'] / merged_df[
        f'agg_ss_score_{control_phase}']
    merged_df['ss_crosser_ratio_test'] = merged_df[f'agg_crosser_score_{test_phase}'] / merged_df[
        f'agg_ss_score_{test_phase}']
    merged_df['ratio_fc'] = merged_df['ss_crosser_ratio_test'] / merged_df['ss_crosser_ratio_control']

    merged_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    merged_df = merged_df.sort_values(by='triangle_fc')

    return merged_df

def mean_loop_score(df):
    for loop_type in ['ss', 'crosser']:
        df[f'{loop_type}_num_loops'] = df[f"{loop_type}_loop_scores"].apply(lambda row: len(row))
        df[f'{loop_type}_mean_loop_score'] = df[f"agg_{loop_type}_score"].divide(df[f'{loop_type}_num_loops'])

    df['total_num_loops'] = df['ss_num_loops'] + df['crosser_num_loops']
    df['total_loop_score'] = df['agg_ss_score'] + df['agg_crosser_score']
    df['total_mean_loop_score'] = df['total_loop_score'].divide(df['total_num_loops'])

    return df


def clean_list(lst):
    eltlist = lst.split("[")[1].split("]")[-2].split(',')
    cleaned_list = [re.sub('[^A-Za-z0-9\-\.]+', '', mystring) for mystring in eltlist]
    cleaned_list = [mystr for mystr in cleaned_list if mystr != '']

    return cleaned_list


def get_loops_by_row(row, loop_type, class_list):
    cleaned_list = clean_list(row[f"{loop_type}_loop_classes"])

    # print(cleaned_list)
    # print(list(cleaned_list[i] for i in range(len(cleaned_list))))

    inds = [i for i in range(len(cleaned_list)) if cleaned_list[i] in class_list]

    cleaned2 = clean_list(row[f"{loop_type}_loop_scores"])
    # print(cleaned2)

    loops = [float(cleaned2[i]) for i in inds]
    # print(loops)

    return sum(loops)  # add other agg methods later


def loop_scores_by_class(df, loop_type, class_list):
    df['classwise_score'] = df.apply(lambda row: get_loops_by_row(row, loop_type, class_list), axis=1)

    return df


def read_bedsites_flex(fname, sort=None, header=0):
    bedsites = pd.read_csv(fname, header=header, sep='\t')

    # if sort exists, then sort on that column and return the dict of dfs split by that column along with the split value
    # if sort doesn't exist, then drop all the columns after chr, start, end
    if sort == 'last':
        out = {}
        for val in bedsites[bedsites.columns[-1]].unique():
            bedsub = bedsites[bedsites.iloc[:, -1] == val]
            bedsub = bedsub.iloc[:, :3]
            bedsub.columns = ['chrom', 'start', 'end']
            out[val] = bedsub.reset_index(drop=True)
    else:
        out = bedsites.iloc[:, :3].reset_index(drop=True)
        out.columns = ['chrom', 'start', 'end']

    return out


def search_dist_for_proms(peaks_df, proms_df, dist):
    peaks_df = peaks_df.copy()
    peaks_df['genes'] = [[] for _ in range(len(peaks_df))]
    peaks_df['l2fc'] = [[] for _ in range(len(peaks_df))]
    for peak_idx, peak_row in peaks_df.iterrows():
        peak_chr = peak_row['chr'].split('chr')[1]
        peak_coord = peak_row['coord']

        # Find genes whose TSS is within X of the peak
        nearby_genes = proms_df[
            (proms_df['chromosome_name'] == peak_chr) &
            (abs(proms_df['transcription_start_site'] - peak_coord) <= dist)
            ]

        peaks_df.at[peak_idx, 'genes'] = nearby_genes['external_gene_name'].tolist()
        if 'log2FoldChange' in proms_df.columns:
            peaks_df.at[peak_idx, 'l2fc'] = nearby_genes['log2FoldChange'].mean()

    peaks_df['num_genes'] = peaks_df['genes'].apply(len)

    return peaks_df


def assign_proms_to_mm(peaks_df, proms_df, min_dist, max_dist, score_col='triangle_score_P1',
                       quant_col='quant'):
    peaks = peaks_df.copy()
    proms = proms_df.copy()

    # normalize chromosome naming
    peaks["chromosome_name"] = peaks["chr"].str.replace("chr", "", regex=False)

    # output columns
    proms["mm_quant"] = np.nan
    proms["dist"] = np.nan
    proms["score"] = np.nan
    proms["mm_peak"] = np.nan

    for chrom, proms_sub in proms.groupby("chromosome_name"):
        peaks_sub = peaks[peaks["chromosome_name"] == chrom]
        if peaks_sub.empty:
            continue

        # sort peaks once
        peaks_sub = peaks_sub.sort_values("coord")
        peak_coords = peaks_sub["coord"].values
        peak_quants = peaks_sub[quant_col].values
        tri_scores = peaks_sub[score_col].values
        # promoter coords
        prom_coords = proms_sub["transcription_start_site"].values

        # find insertion positions
        idxs = np.searchsorted(peak_coords, prom_coords)

        # clip to valid range
        idxs_left = np.clip(idxs - 1, 0, len(peak_coords) - 1)
        idxs_right = np.clip(idxs, 0, len(peak_coords) - 1)

        # compute distances to neighbors
        d_left = np.abs(prom_coords - peak_coords[idxs_left])
        d_right = np.abs(prom_coords - peak_coords[idxs_right])

        # choose nearest
        use_left = d_left <= d_right
        nearest_idxs = np.where(use_left, idxs_left, idxs_right)
        d = np.minimum(d_left, d_right)

        valid = (d <= max_dist) & (d >= min_dist)

        proms.loc[proms_sub.index[valid], "mm_quant"] = peak_quants[nearest_idxs[valid]]
        proms.loc[proms_sub.index[valid], "dist"] = d[valid]
        proms.loc[proms_sub.index[valid], "score"] = tri_scores[nearest_idxs[valid]]
        proms.loc[proms_sub.index[valid], "mm_peak"] = peak_coords[nearest_idxs[valid]]

    return proms

def get_locus_from_coord(row, flank):
    return f"{row['chr']}:{int(row['coord'])-flank}-{int(row['coord'])+flank}"

def plot_peak(row, clr, flank=200_000, res=2000, fs=5, score_col=None):
    fig, ax = plt.subplots(figsize=(fs, fs))

    locus = get_locus_from_coord(row, flank)
    try:
        mat = utils.get_matr(clr, locus, oe=False)
    except ValueError:
        return fig

    ax.matshow(mat, norm=LogNorm(), cmap='fall')
    ax.axvline((flank // res), ls='--', color='k', alpha=0.5)
    ax.axhline((flank // res), ls='--', color='k', alpha=0.5)
    if score_col is not None:
        try:
            plt.title(f"{locus}: score {row[score_col]:.5f}")
        except:
            return fig

    return fig

def make_pileup_combo(loops, bed_dfs, tf_method='inner', loop_method='both', vf=None):
    # search loops for either having all of beds on both sides, at least one on both sides, etc...
    loops = loops.rename({'chr1': 'chrom'}, axis=1)
    # print(loops)

    # make sites as the one/two sided intersection/union (inner/outer)
    df_final = reduce(lambda left, right: pd.merge(left, right, on=['chrom', 'start', 'end'], how=tf_method), bed_dfs)

    df_final['coord'] = df_final[['start', 'end']].mean(axis=1)
    df_final.sort_values(by='coord', inplace=True)
    # print(df_final)

    loops_a1 = loops.copy()
    loops_a2 = loops.copy()
    loops_a1['loop_idx'] = loops.index
    loops_a2['loop_idx'] = loops.index
    loops_a1.sort_values(by=['anchor1'], inplace=True)
    loops_a2.sort_values(by=['anchor2'], inplace=True)

    #     print(loops_a1)
    #     print(loops_a2)
    anchor1_matches = pd.merge_asof(
        loops_a1, df_final, right_on='coord', left_on='anchor1', by='chrom', direction='nearest', tolerance=2000
    )
    anchor2_matches = pd.merge_asof(
        loops_a2, df_final, right_on='coord', left_on='anchor2', by='chrom', direction='nearest', tolerance=2000
    )
    # print(anchor1_matches)
    # Mark which anchors matched

    anchor1_matched = anchor1_matches[~anchor1_matches['start'].isna()]
    anchor2_matched = anchor2_matches[~anchor2_matches['start'].isna()]
    idx1 = set(anchor1_matched['loop_idx'])
    idx2 = set(anchor2_matched['loop_idx'])

    if loop_method == 'either':
        matched_idx = idx1 | idx2
    elif loop_method == 'both':
        matched_idx = idx1 & idx2
    else:
        raise ValueError(f"Invalid loop_method: {loop_method}")

    selected_loops = loops.loc[list(matched_idx)].rename({'chrom': 'chrom1', 'chr2': 'chrom2'}, axis=1)

    sites = selected_loops[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]

    return selected_loops


def make_exclusive_combo(loops, bed_dfs, tf_method='inner', loop_method='both', vf=None):
    # search loops for either having all of beds on both sides, at least one on both sides, etc...
    loops = loops.rename({'chr1': 'chrom'}, axis=1)
    # print(loops)
    # make sites as the one/two sided intersection/union (inner/outer)
    df_final = reduce(lambda left, right: pd.merge(left, right, on=['chrom', 'start', 'end'], how=tf_method), bed_dfs)

    df_final['coord'] = df_final[['start', 'end']].mean(axis=1)
    df_final.sort_values(by='coord', inplace=True)
    print(df_final)

    loops_a1 = loops.copy()
    loops_a2 = loops.copy()
    loops_a1['loop_idx'] = loops.index
    loops_a2['loop_idx'] = loops.index
    loops_a1.sort_values(by=['anchor1'], inplace=True)
    loops_a2.sort_values(by=['anchor2'], inplace=True)

    #     print(loops_a1)
    #     print(loops_a2)
    anchor1_matches = pd.merge_asof(
        loops_a1, df_final, right_on='coord', left_on='anchor1', by='chrom', direction='nearest', tolerance=2000
    )
    anchor2_matches = pd.merge_asof(
        loops_a2, df_final, right_on='coord', left_on='anchor2', by='chrom', direction='nearest', tolerance=2000
    )
    # print(anchor1_matches)
    # Mark which anchors matched

    anchor1_matched = anchor1_matches[~anchor1_matches['start'].isna()]
    anchor2_matched = anchor2_matches[~anchor2_matches['start'].isna()]
    idx1 = set(anchor1_matched['loop_idx'])
    idx2 = set(anchor2_matched['loop_idx'])

    if loop_method == 'either':
        matched_idx = idx1 | idx2
    elif loop_method == 'both':
        matched_idx = idx1 & idx2
    else:
        raise ValueError(f"Invalid loop_method: {loop_method}")

    selected_loops = loops.loc[list(matched_idx)].rename({'chrom': 'chrom1', 'chr2': 'chrom2'}, axis=1)

    #sites = selected_loops[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]

    return selected_loops