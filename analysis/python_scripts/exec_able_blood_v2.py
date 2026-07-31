#
import numpy as np
import os
import pandas as pd
import utils
import cooltools
import multiprocessing as mp
import cooler

# import the loop quantification module
import sys
# sys.path.insert(1, LOCAL_LQ_PATH)
import looptools_local as looptools  # all the back end code is here

def calculate_and_save_Ps_curve(clr, cap_string=None, nproc=8, max_sep=3000000, output_filename=None, col = 'balanced.avg.smoothed.agg'):
    """Calculate and save the P(s) curve (average across chromosomes, weighted by chromosome size) and save to file. max_sep is the maximum genomic separation in bp to which to calculate the P(s) curve."""

    res = clr.binsize

    if cap_string is None:
        view_frame = utils.get_hg38_arms(clr, excl_ym=True)
    else:
        chr_name, st, ed = utils.get_coords(cap_string)
        # make data frame
        view_frame = pd.DataFrame([[chr_name, st, ed, chr_name]], columns=['chrom','start','end', 'name'])

    print(view_frame)
    # calculate the P(s) curve
    cvd = cooltools.expected_cis(clr=clr, view_df = view_frame, smooth=True, aggregate_smoothed=True, nproc=nproc)
    P_s_data_chr = cvd.loc[:, col]

    # save as txt file
    if output_filename is None:
        output_filename = f"P_s_{clr.filename.split('/')[-1].split('.')[0]}_{res}bp.txt"

    np.savetxt(output_filename, P_s_data_chr)
    pass


def calculate_and_save_avg_Ps_curve(clr, nproc=4, max_sep=5000000, output_filename=None):
    """Calculate and save the P(s) curve (average across chromosomes, weighted by chromosome size) and save to file. max_sep is the maximum genomic separation in bp to which to calculate the P(s) curve."""

    res = clr.binsize

    # calculate the P(s) curve
    cvd = cooltools.expected_cis(clr=clr, smooth=True, aggregate_smoothed=True, nproc=nproc)
    cvd['s_bp'] = cvd['dist'] * res

    chr_names = [chrom_name for chrom_name in clr.chromnames if len(chrom_name) > 3 and (
                chrom_name[3:].isnumeric() or chrom_name[3:] == 'X')]  # only take numbered choromosomes and chrX

    # average across chromosomes
    P_s_data_all_chrs = np.zeros((1 + max_sep // res, len(chr_names)))
    for i, chr_name in enumerate(chr_names):
        P_s_data_chr = cvd.loc[cvd['region1'] == chr_name, np.array(['s_bp', 'balanced.avg'])]
        P_s_data_chr = P_s_data_chr.loc[P_s_data_chr['s_bp'] <= max_sep]
        P_s_data_chr = P_s_data_chr.sort_values('s_bp')
        assert (np.all(P_s_data_chr['s_bp'] == np.arange(0, max_sep + 1, res)))
        P_s_data_all_chrs[:, i] = P_s_data_chr['balanced.avg']

    chrom_weights = clr.chromsizes[chr_names].values
    chrom_weights = chrom_weights / np.sum(chrom_weights)

    P_s_data_averaged_chrs = np.average(P_s_data_all_chrs, axis=1, weights=chrom_weights)

    # save as txt file
    if output_filename is None:
        output_filename = f"P_s_{clr.filename.split('/')[-1].split('.')[0]}_{res}bp.txt"
    np.savetxt(output_filename, P_s_data_averaged_chrs)

def able(PS, clr, phase):
    global _loops
    global _test_range

    if not os.path.exists(PS):
        print("Writing PS...")
        calculate_and_save_avg_Ps_curve(clr, output_filename=PS)

    P_s_values = np.loadtxt(PS)
    # print(P_s_values)

    # initialize a loop quantifier object
    lq = looptools.LoopQuantifier(clr, P_s_values,
                                  gaussian_blur_sigma_px=3, outlier_removal_radius_px=3,
                                  na_stripe_dist_to_center_px_cutoff=1)

    # output
    loops_curr = []
    toplot=False
    for i, loop in _loops.iterrows():
        if i in _test_range:
            toplot = True

        chr_name = loop.iloc[0]
        left_coord = int(np.mean([loop.iloc[1], loop.iloc[2]]))
        right_coord = int(np.mean([loop.iloc[4], loop.iloc[5]]))
        score = lq.quantify_loop(chr_name,
                                 left_coord,
                                 right_coord,
                                 local_region_size=100000,
                                 # how large of a region to use to estimate the background (bp)
                                 quant_region_size=10000,
                                 k_min=2,  # parameter for outlier detection;
                                 # lower = more sensitive, higher = less sensitive;
                                 # in a Gaussian-blurred image of the local region, consider a point an outlier if it is a local maximum and has a value at least k_min times the median value
                                 show_plot=toplot  # whether to show the plot (useful for diagnostic purposes)
                                 )
        toplot = False
        loops_curr.append(score)
        if i%1000==0:
            np.savetxt(f'/mnt/md0/varshini/Analysis/Blood/MicroC_Loops_Anchors_Merged/temp_loops_2kb{phase}.txt', np.array(loops_curr))
    return loops_curr

def able_dummy(PS, clr):
    global _loops
    global _test_range

    # print(_loops.head())
    # print(test_range)
    # print(PS)

    print(int(np.mean([_loops.iloc[0, 1], _loops.iloc[0, 2]])))

    return 0
def init_worker(loops, test_range):
    global _loops
    _loops = loops  # Assign x to the global variable in each worker

    global _test_range
    _test_range = test_range  # Assign x to the global variable in each worker

if __name__=='__main__':
    # variable definitions
    phases = ['Phase1_', 'Phase2_dedup', 'Phase3_dedup']
    resolution = 2000

    microc_cooler_path = '/mnt/coldstorage/Varshmallow/Adipose_Blood_Merged/sankaran'
    microc_cooler_names = [os.path.join(microc_cooler_path, f"Sankaran_{p}merged.250.mcool") for p in phases]
    microc_coolers = [utils.get_clr(clr_name, resolution) for clr_name in microc_cooler_names]
    ps_path = '/mnt/md0/varshini/Analysis/Blood/microc_expected'
    PS_names = [f"P_s_{f.split('/')[-1].split('.')[0]}_{resolution}bp_wavg_verify2.tsv" for f in microc_cooler_names]
    PS_names = [os.path.join(ps_path, f) for f in PS_names]

    # output_path = '/mnt/md0/varshini/Analysis/Blood/MicroC_Loops_Anchors_v2/'
    # loopcall_path = '/mnt/md0/varshini/Analysis/Blood/MicroC_Loops_Anchors_v2/reseq_calls_125kb_q0.02_2kbanch.bedpe'

    #output_path = '/mnt/md0/varshini/Analysis/Blood/MicroC_Loops_Anchors_Merged/'
    #loopcall_path = '/mnt/md0/varshini/Analysis/Blood/MicroC_Loops_Anchors_Merged/merged_calls_125kb_q0.02_5kb.bedpe'
    #loopcall_path = '/mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/merged_calls_125kb_q0.01.bedpe'

    output_path = '/mnt/md0/varshini/Analysis/Blood/complete_analyses/loops_anchors/'
    loopcall_path = '/mnt/md0/varshini/Analysis/Blood/complete_analyses/loops_anchors/merged_calls_125kb_q0.01.bedpe'
    loops = pd.read_csv(loopcall_path, sep='\s+', header=None)

    test_range = [0, 100, 500]
    loops_out = loops.copy()

    ## code
    with mp.Pool(3, initializer=init_worker, initargs=(loops,test_range,)) as p:
        out = p.starmap(able, zip(PS_names, microc_coolers, phases))

    for i, phase in enumerate(phases):
        loops_out[f'Strength_{phase}'] = out[i]
    print(loops_out)
    loops_out.to_csv(os.path.join(output_path, f'able_scores_125kb_q0.01{resolution}bp_wavg_verify.tsv'), sep='\t', header=False, index=False)
