import itertools
from pathlib import Path
import shutil
import os
import re
import pandas as pd
import ast
import numpy as np
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file, save_hdf5_file
from polychrom.contactmaps import monomerResolutionContactMapSubchains, monomerResolutionContactMap, averageContacts
import scipy as sp
from functools import partial
import time
from multiprocessing import Pool
import matplotlib.pyplot as plt


## below from james' able code (Jusuf et al., NSMB (2026))
def get_cm(input_path, block_start, block_end, block_step, **kwargs):

    chrom_size = kwargs.get('chrom_size', 70000)
    region_size = kwargs.get('region_size', 2000)
    cutoff = kwargs.get('cutoff', 3)
    nproc = kwargs.get('nproc', 8)
    balance = kwargs.get('balance', False)

    num_regions = chrom_size // region_size
    region_starts = np.arange(num_regions) * region_size

    files = list_URIs(input_path)
    print(f'+LE +EP: {len(files)} blocks')

    num_blocks = ((block_end - block_start) / block_step)
    num_regions = chrom_size / region_size

    start_time = time.time()

    if balance:
        # print("Loading...")
        # cm = monomerResolutionContactMap(filenames=files[block_start:block_end:block_step],
        #                                          cutoff=cutoff, n=nproc)
        # print("Balancing...")
        # cm_balanced = sinkhorn_knopp(cm)
        # print("Done balancing")
        # cm_repeat = averageContacts(cm_balanced, range(num_regions), region_size)
        # cm_rescaled = cm_repeat / num_blocks / num_regions
        # # get all the region_size x region_size squares and sum across them
        cm = monomerResolutionContactMapSubchains(filenames=files[block_start:block_end:block_step],
                                                  mapStarts=region_starts, mapN=region_size, cutoff=cutoff, n=nproc)

        end_time = time.time()
        print(f'done! ({int(end_time - start_time)} sec)')

        cm = cm / num_blocks / num_regions
        pad_width = region_size
        n_new = 2*pad_width + region_size
        print("Creating P(s) pad...")
        x1, y1 = get_Ps_curve_sim(cm)
        y1_ext = np.concatenate([y1, np.full(n_new - len(y1), np.median(y1[-50:]))])
        cm_padded = sp.linalg.toeplitz(y1_ext)

        #cm_padded = np.random.uniform(0, 0.0015, size=(pad_width * 2 + region_size, pad_width * 2 + region_size))
        cm_padded[pad_width:(pad_width+region_size), pad_width:(pad_width+region_size)] = cm

        plt.show()
        print("Balancing...")
        cm_balanced = sinkhorn_knopp(cm_padded)
        print("Done balancing")

        cm_rescaled = cm_balanced[pad_width:(pad_width+region_size), pad_width:(pad_width+region_size)]
        print(f"In padded matrix {np.shape(cm_balanced)}, extracting {pad_width}-{pad_width+region_size}")
    else:
        print("Not balancing")

        cm = monomerResolutionContactMapSubchains(filenames=files[block_start:block_end:block_step],
                                                  mapStarts=region_starts, mapN=region_size, cutoff=cutoff, n=nproc)
        end_time = time.time()
        print(f'done! ({int(end_time - start_time)} sec)')
        cm_rescaled = cm / num_blocks / num_regions

    return cm_rescaled



# Example

# implementation lifted from this blog post
# https://fulkast.medium.com/the-sinkhorn-knopp-algorithm-without-proof-697c9af7df7
# and then followed logic that element wise mult should work for the diagonal dots
def sinkhorn_knopp(cm, n_iterations=25):
    n = np.shape(cm)[0]
    cm_iter = cm.copy()

    D1 = np.ones(n)
    D2 = np.ones(n)

    for i in range(n_iterations):
        column_sums = cm_iter.sum(axis=0)
        row_sums = cm_iter.sum(axis=1)
        if (i % 2) == 0:
            D1 *= (1.0 / row_sums)
        else:
            D2 *= (1.0 / column_sums)

        dot_D2 = cm * D2
        cm_iter = D1[:, None] * dot_D2
    return cm_iter

def get_Ps_curve_sim(contact_map):
    '''Get the P(s) curve of a simulated monomer-resolution contact map.'''
    dist_arr = np.arange(0, contact_map.shape[0])  # array of distances to calculate P(s) for (in monomers)
    mean_contact_prob_arr = np.zeros(len(dist_arr))
    for i, dist in enumerate(dist_arr):
        contact_probs_along_diag = np.diag(contact_map, dist)
        mean_contact_prob_arr[i] = np.mean(contact_probs_along_diag)  # ignore nan's
    return dist_arr, mean_contact_prob_arr


def _process_single_block(block_file, num_regions, region_size, cutoff):
    region_starts = np.arange(num_regions) * region_size
    cm_single = monomerResolutionContactMapSubchains(filenames=[block_file],
                                                     mapStarts=region_starts, mapN=region_size, cutoff=cutoff, n=1)
    return cm_single.astype(np.float32) / num_regions


def get_cm_by_block(n_cores, n_blocks, block_files, num_regions=35, region_size=2000, cutoff=3):
    # get each CM without averaging over blocks
    # so the stack should be region_size x region_size x n_blocks

    if len(block_files) != n_blocks:
        print(f"found {len(block_files)} blocks, but requested {n_blocks}")

    worker_func = partial(
        _process_single_block,
        num_regions=num_regions,
        region_size=region_size,
        cutoff=cutoff
    )

    # split up contact map blocks over n_cores
    with Pool(processes=n_cores) as pool:
        cm_list = pool.map(worker_func, block_files)

    # Stack into a 3D array: (region_size, region_size, n_blocks)
    block_stack = np.dstack(cm_list)

    return block_stack


def jackknife_cm(block_stack, loops, w=5):
    loops.reset_index(inplace=True, drop=True)
    n_rows, n_cols, n_blocks = block_stack.shape

    # Precompute the sum of all blocks to easily calculate the Leave-One-Out (LOO) mean
    grand_sum_cm = np.sum(block_stack, axis=2)

    # Dictionary to hold the LOO metrics for each loop
    # Structure: { loop_idx: {'obs': [], 'exp': [], 'oe': []} }
    loo_data = {i: {'obs': [], 'exp': [], 'oe': []} for i in range(len(loops))}

    # For k in block_num... drop k, get the average of the remaining contact maps
    for k in range(n_blocks):

        # 1. Calculate the Leave-One-Out Mean Contact Map
        loo_cm = (grand_sum_cm - block_stack[:, :, k]) / (n_blocks - 1)

        # 2. Get expected matrix (P(s) curve) for this specific LOO map
        # Assuming get_Ps_curve_sim is available in your namespace
        x1, y1 = get_Ps_curve_sim(loo_cm)
        loo_ps_matrix = sp.linalg.toeplitz(y1)

        # 3. Faster to do this loop-by-loop: extract the 9x9 patches
        for i, loop_row in loops.iterrows():

            x_start, x_end = loop_row.left - w, loop_row.left + w
            y_start, y_end = loop_row.right - w, loop_row.right + w

            obs_dot = np.mean(loo_cm[x_start:x_end, y_start:y_end])
            exp_dot = np.mean(loo_ps_matrix[x_start:x_end, y_start:y_end])

            # Store values for this specific jackknife iteration
            loo_data[i]['obs'].append(obs_dot)
            loo_data[i]['exp'].append(exp_dot)
            loo_data[i]['oe'].append(obs_dot / exp_dot)

    # populate loop_stats_df with information
    loop_stats_list = []

    for i, loop_row in loops.iterrows():
        obs_arr = np.array(loo_data[i]['obs'])
        oe_arr = np.array(loo_data[i]['oe'])

        # Calculate Grand Means
        mean_obs = np.mean(obs_arr)
        mean_oe = np.mean(oe_arr)

        # Calculate Jackknife Standard Error: sqrt( (N-1)/N * sum( (x_i - x_mean)^2 ) )
        se_obs = np.sqrt(((n_blocks - 1) / n_blocks) * np.sum((obs_arr - mean_obs) ** 2))
        se_oe = np.sqrt(((n_blocks - 1) / n_blocks) * np.sum((oe_arr - mean_oe) ** 2))

        curr_row_dict = loop_row.to_dict()
        curr_loop_stats = {
            'loop_idx': i,
            'obs_mean': mean_obs,
            'obs_se': se_obs,
            'oe_mean': mean_oe,
            'oe_se': se_oe
        }
        append_list = {**curr_row_dict, **curr_loop_stats}
        loop_stats_list.append(append_list)

    loop_stats_df = pd.DataFrame(loop_stats_list)
    return loop_stats_df


def jackknife_cm_dummy(loops, w=5, n_blocks=2):
    loops.reset_index(inplace=True, drop=True)

    # Dictionary to hold the LOO metrics for each loop
    # Structure: { loop_idx: {'obs': [], 'exp': [], 'oe': []} }
    loo_data = {i: {'obs': [], 'exp': [], 'oe': []} for i in range(len(loops))}

    # For k in block_num... drop k, get the average of the remaining contact maps
    for k in range(n_blocks):

        for i, loop_row in loops.iterrows():

            obs_dot = 1.5
            exp_dot = 3

            # Store values for this specific jackknife iteration
            loo_data[i]['obs'].append(obs_dot)
            loo_data[i]['exp'].append(exp_dot)
            loo_data[i]['oe'].append(obs_dot / exp_dot)

    # populate loop_stats_df with information
    loop_stats_list = []

    for i, loop_row in loops.iterrows():
        obs_arr = np.array(loo_data[i]['obs'])
        oe_arr = np.array(loo_data[i]['oe'])

        # Calculate Grand Means
        mean_obs = np.mean(obs_arr)
        mean_oe = np.mean(oe_arr)

        # Calculate Jackknife Standard Error: sqrt( (N-1)/N * sum( (x_i - x_mean)^2 ) )
        se_obs = np.sqrt(((n_blocks - 1) / n_blocks) * np.sum((obs_arr - mean_obs) ** 2))
        se_oe = np.sqrt(((n_blocks - 1) / n_blocks) * np.sum((oe_arr - mean_oe) ** 2))

        curr_row_dict = loop_row.to_dict()
        curr_loop_stats = {
            'loop_idx': i,
            'obs_mean': mean_obs,
            'obs_se': se_obs,
            'oe_mean': mean_oe,
            'oe_se': se_oe
        }
        append_list = {**curr_row_dict, **curr_loop_stats}
        loop_stats_list.append(append_list)

    loop_stats_df = pd.DataFrame(loop_stats_list)
    return loop_stats_df


# add the O/E for each into sim_loops_df (expecteds should be very similar...)
def add_loops_single_cm(cm_rescaled, loops,
                        ps_matrix_base=None, w=5, plot_ind=[]):
    # compute the ps matrix from the ps curve
    x1, y1 = get_Ps_curve_sim(cm_rescaled)

    ps_matrix = sp.linalg.toeplitz(y1)

    loops.reset_index(inplace=True, drop=True)

    loop_stats_list = []
    # for each loop, make a (2w+1)kb x (2w+1)kb block centered on the dot
    for i, loop_row in loops.iterrows():
        x_start, x_end = loop_row.left - w, loop_row.left + w
        y_start, y_end = loop_row.right - w, loop_row.right + w

        obs_dot = np.mean(cm_rescaled[x_start:x_end, y_start:y_end])
        exp_dot = np.mean(ps_matrix[x_start:x_end, y_start:y_end])

        curr_row_dict = loop_row.to_dict()
        curr_loop_stats = {
            'loop_idx': i,
            'obs_mean': obs_dot,
            'obs_se': 0,
            'oe_mean': obs_dot/exp_dot,
            'oe_se': 0
        }
        append_list = {**curr_row_dict, **curr_loop_stats}
        loop_stats_list.append(append_list)

    loop_stats_df = pd.DataFrame(loop_stats_list)
    return loop_stats_df

def get_loop_df(**kwargs):
    # first, make the np array of the whole thing and save it


    CTCF_L = kwargs.get('CTCF_L')
    CTCF_R = kwargs.get('CTCF_R')
    EP = kwargs.get('EP')
    loading_loc = kwargs.get('MM')

    conv_pairs = []
    for rC in CTCF_R:
        ins_idx = np.searchsorted(CTCF_L, rC, side='right')
        pairs = [(rC, CTCF_L[i]) for i in range(ins_idx, len(CTCF_L))]

        conv_pairs.extend(pairs)
    # for rL in CTCF_L:
    #    ins_idx = np.searchsorted(CTCF_R, rL, side='left')
    #    pairs = [(CTCF_R[i-1], rL) for i in range(1, ins_idx+1)]
    #    conv_pairs.extend(pairs)

    # L within L and R within R
    all_pairs = CTCF_R + CTCF_L
    all_pairs.sort()
    nonconv_pairs = list(itertools.combinations(all_pairs, 2))

    # remove the convergent pairs
    nonconv_pairs = [x for x in nonconv_pairs if x not in conv_pairs]

    EP_pairs = list(itertools.combinations(EP, 2))

    # convergent loops directly in range of MM
    MM_inner_loops = []
    MM_outer_loops = []

    for mm in loading_loc:
        idx_l = np.searchsorted(CTCF_L, mm, side='right')
        idx_r = np.searchsorted(CTCF_R, mm, side='left')
        pair = (CTCF_R[idx_r - 1], CTCF_L[idx_l])
        MM_inner_loops.append(pair)

        outer_pair = pair = (CTCF_R[idx_r - 2], CTCF_L[idx_l + 1])
        MM_outer_loops.append(outer_pair)

    conv_df = pd.DataFrame(conv_pairs, columns=['left', 'right'])
    conv_df['loop_type'] = 'CTCF_conv'

    nonconv_df = pd.DataFrame(nonconv_pairs, columns=['left', 'right'])
    nonconv_df['loop_type'] = 'CTCF_other'

    EP_df = pd.DataFrame(EP_pairs, columns=['left', 'right'])
    EP_df['loop_type'] = 'EP'

    all_loops_df = pd.concat([conv_df, nonconv_df, EP_df], axis=0)

    all_loops_df['MM'] = None
    all_loops_df.loc[all_loops_df.set_index(['left', 'right']).index.isin(MM_inner_loops), 'MM'] = 'inner'
    all_loops_df.loc[all_loops_df.set_index(['left', 'right']).index.isin(MM_outer_loops), 'MM'] = 'outer'

    return all_loops_df


def process_folder(folder_path, arr_output, **kwargs):
    # infer block start, end, and step from folder
    block_step = kwargs.get('block_step', 1)
    block_start = kwargs.get('block_start', 0)
    block_end = kwargs.get('block_end', 0)
    overwrite = kwargs.get('overwrite', False)
    nproc = kwargs.get('nproc', 24)
    loop_width = kwargs.get('loop_width', 5)
    min_blocks = kwargs.get('min_blocks', 600)
    cutoff = kwargs.get('cutoff', 3)
    balance = kwargs.get('balance', False)
    loop_strength_method = kwargs.get('loop_strength_method', 'jackknife')

    files = list_URIs(folder_path)
    if block_end == 0: # infer
        block_end = len(files)

    num_blocks = block_end - block_start
    if num_blocks < min_blocks:
        return None

    # save into an array output folder
    if not os.path.exists(arr_output) or overwrite:
        cm_rescaled = get_cm(folder_path, block_start, block_end, block_step, balance=balance, nproc=nproc)
        np.save(arr_output, cm_rescaled, allow_pickle=False)

    # get all the loops to analyze (fast, so can always run this step)
    loop_df = get_loop_df(**kwargs)

    # jackknife and get loop strengths
    if loop_strength_method == 'jackknife':
        block_stack = get_cm_by_block(nproc, num_blocks, files, cutoff=cutoff)
        loop_df_jn = jackknife_cm(block_stack, loop_df, w=loop_width)
    elif loop_strength_method == 'single':
        loop_df_jn = add_loops_single_cm(cm_rescaled, loop_df)
    else:
        raise ValueError("Method not supported")
    return loop_df_jn

# makes sure it works before sending it to the long processing function
def process_folder_dummy(folder_path, arr_output, **kwargs):

    # infer block start, end, and step from folder
    block_step = kwargs.get('block_step', 1)
    block_start = kwargs.get('block_start', 0)
    block_end = kwargs.get('block_end', 0)
    overwrite = kwargs.get('overwrite', False)
    nproc = kwargs.get('nproc', 24)
    loop_width = kwargs.get('loop_width', 5)

    files = list_URIs(folder_path)
    if block_end == 0: # infer
        block_end = len(files)

    num_blocks = block_end - block_start
    print(f"Num blocks: {num_blocks}")

    # save into an array output folder
    if not os.path.exists(arr_output) or overwrite:
        print(f"Writes {arr_output}")


    loop_df = get_loop_df(**kwargs)
    print(loop_df.head())

    # jackknife and get loop strengths
    #block_stack = get_cm_by_block(nproc, num_blocks, files)
    loop_df_jn = jackknife_cm_dummy(loop_df,loop_width)

    return loop_df_jn
def main(root_dir, results_dir, **kwargs):
    subdirs = sorted([p for p in root_dir.iterdir() if p.is_dir()])
    overwrite = kwargs.get('overwrite', False)
    # mapping table
    mapping_rows = []

    for idx, subdir in enumerate(subdirs):

        # short unique ID
        run_id = f"run_{idx:03d}"

        simdirs = sorted([p for p in subdir.iterdir() if p.is_dir()])

        for idx2, simdir in enumerate(simdirs):

            if "sim_done.txt" not in os.listdir(os.path.join(subdir, simdir)):
                continue

            sim_id = f"sim_{idx2:02d}"
            result_file = results_dir / f"{run_id}_{sim_id}_loops.tsv"
            param_src = subdir / simdir / "run_parameters.txt"

            if not os.path.exists(result_file) or overwrite:
                result = process_folder(simdir, os.path.join(results_dir, f"{run_id}_{sim_id}_arr.npy"), **kwargs)
            else:
                result = None

            if result is not None:

                result.to_csv(result_file, sep="\t")

                if param_src.exists():
                    param_dst = results_dir / f"{run_id}_{sim_id}_run_parameters.txt"
                    shutil.copy2(param_src, param_dst)

            curr_row = {
                "run_id": run_id,
                "sim_id": sim_id,
                "original_path": str(simdir.resolve()),
                "result_file": result_file.name,
                "parameters_file": (
                    f"{run_id}_{sim_id}_parameters.txt"
                    if param_src.exists()
                    else None
                )
            }

            with open(os.path.join(subdir, simdir, "run_parameters.txt"), "r") as f:
                x = f.read().splitlines()

            param_dict = ast.literal_eval(x[0])
            del param_dict['longlived_fraction']
            del param_dict['longlived_boost_factor']
            del param_dict['dsb_boost_factor']

            if 'vel_steps' not in param_dict.keys():
                param_dict['vel_steps'] = (0.0025, 50)
            map_row = {**curr_row, **param_dict}
            mapping_rows.append(map_row)

    mapping_df = pd.DataFrame(mapping_rows)
    mapping_df.to_csv(results_dir / "mapping.tsv", index=False, sep="\t")

    print(f"Processed {len(subdirs)} folders.")

if __name__ == '__main__':
    # root_dir = Path("/mnt/md1/varshini/Blood/sim_data_velocity")
    # results_dir = Path("/mnt/md1/varshini/Blood/sim_velocity_processed_v2")
    #
    # results_dir.mkdir(exist_ok=True)
    #
    # main(root_dir, results_dir, CTCF_L=CTCF_L, CTCF_R=CTCF_R, EP=EP, MM=MM,
    #      min_blocks=2400)



    CTCF_L = [574, 694, 866, 1241, 1390, 1580, 1752, 1800]
    CTCF_R = [200, 330, 724, 1425, 1433, 1604]
    EP = [250, 372, 540, 745, 775, 833, 961, 1202, 1330, 1640, 1722]
    MM = [456, 1054, 1507]

    # root_dir = Path("/mnt/md1/varshini/Blood/sim_data_density_subset")
    #
    # results_dir = Path(f"/mnt/md1/varshini/Blood/sim_density_processed_subset_notBal")
    # results_dir.mkdir(exist_ok=True)

    #main(root_dir, results_dir, CTCF_L=CTCF_L, CTCF_R=CTCF_R, EP=EP, MM=MM, overwrite=True, balance=False)

    #results_dir = Path(f"/mnt/md1/varshini/Blood/sim_density_processed_subset_Bal")
    #results_dir.mkdir(exist_ok=True)

    #main(root_dir, results_dir, CTCF_L=CTCF_L, CTCF_R=CTCF_R, EP=EP, MM=MM, overwrite=False, balance=True,
    #      loop_strength_method='single')

    # for cutoff in [1, 2, 4, 5]:
    #     results_dir = Path(f"/mnt/md1/varshini/Blood/sim_density_processed_subset_cutoff{cutoff}")
    #     results_dir.mkdir(exist_ok=True)
    #     main(root_dir, results_dir, CTCF_L=CTCF_L, CTCF_R=CTCF_R, EP=EP, MM=MM, overwrite=False, balance=False,
    #          cutoff=cutoff)

    root_dir = Path("/mnt/md1/varshini/Blood/sim_data_compaction_fountain")

    results_dir = Path(f"/mnt/md1/varshini/Blood/sim_compaction_fountain_processed_Bal")
    results_dir.mkdir(exist_ok=True)

    main(root_dir, results_dir, CTCF_L=CTCF_L, CTCF_R=CTCF_R, EP=EP, MM=MM, overwrite=True, balance=True,
         loop_strength_method='single')



