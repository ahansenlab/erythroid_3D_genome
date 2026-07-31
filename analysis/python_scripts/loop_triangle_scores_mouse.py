# get alll the triangle scores

import acr_computation as ac
import os
import pandas as pd
import utils
import argparse
import sys

def main():
    # hard-coded inputs (can change to user input later)



    parser = argparse.ArgumentParser()
    parser.add_argument('datapath', type=str)
    parser.add_argument('subdir', type=str)
    parser.add_argument('filename', type=str)
    parser.add_argument('-p', '--perturbation', type=bool, default=False)
    parser.add_argument('-s', '--sample', type=int, default=5000)
    parser.add_argument('-d', '--dists', nargs="*", type=int, default=[50000, 100000, 200000])
    parser.add_argument('-o', '--outpath', type=str,
                        default='/mnt/md0/varshini/Analysis/Blood/MicroC_Loops_Anchors/triangle_scores_v2')


    args = parser.parse_args()

    f = args.filename
    datapath = args.datapath
    subdir = args.subdir
    is_perturbation = args.perturbation
    dists = args.dists
    sample = args.sample

    outpath=args.outpath
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    resolution = 5000
    phases = ['SMC3_asyn_auxin',
            'SMC3_asyn_no_auxin',
            'new_nipbl_auxin_merged_4h',
            'new_nipbl_untreated_merged_4h']

    microc_cooler_path = '/mnt/md1/varshini/Blood/aboreden_hic/'
    microc_cooler_names = [os.path.join(microc_cooler_path, f"{p}.mcool") for p in phases]

    microc_coolers = [utils.get_clr(clr_name, resolution) for clr_name in microc_cooler_names]

    ps_path = '/mnt/md0/varshini/Analysis/Blood/microc_expected'

    PS_names = [f"P_s_{f.split('/')[-1].split('.')[0]}_{resolution}bp_allcol.tsv" for f in microc_cooler_names]
    PS_names = [os.path.join(ps_path, f) for f in PS_names]

    PS_dfs = utils.write_get_ps(zip(PS_names, microc_coolers))
    zipped_clrs = zip(PS_dfs, microc_coolers, phases)

    cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'prob']
    cols.extend([f'strength_{cond}' for cond in phases])


    if is_perturbation:
        curr_bed = ac.read_peaks(os.path.join(datapath, subdir, f))
    else:
        bedfile = ac.read_df(os.path.join(datapath, subdir, f), sep='\t', columns=['chr', 'coord1', 'coord2'])
        bedfile['perturbation_name'] = f.split('.')[0]
        bedfile['coord'] = bedfile[['coord1', 'coord2']].mean(axis=1)
        curr_bed = bedfile.sort_values(by=['chr','coord'])

    for PS, clr, phase in zipped_clrs:
        for dist in dists:
            print(f"Processing {phase} at {dist} for {f.split('.')[0]}...", file=sys.stderr)
            print(curr_bed.head(), file=sys.stderr)
            fname = f"{phase}_{f.split('.')[0]}_{dist}_withpeaks.tsv"
            if not os.path.exists(os.path.join(outpath, fname)):
                out = ac.snip_region_withpeaks(clr, PS, curr_bed, dist, test_range=[0], verbose=False)
                out.to_csv(os.path.join(outpath, fname), sep='\t')

if __name__=='__main__':
    main()

    # for each type, make all the triangle scores by ACR
    # it should be a dataframe with three columns like phase1 triangle, phase2 triangle, phase3 triangle

    # for f in os.listdir(os.path.join(datapath, subdir)):
    #     print(f)
    # f='erythroid_tf_acrs_all_enhancernotpromoternotctcf.bed'
    # f='erythroid_tf_acrs_all_promoternotctcf.bed'

    # datapath = '/mnt/md0/varshini/Analysis/Blood/MicroC_Loops_Anchors'
    # subdir = 'classifications2'
    # f='erythroid_tf_acrs_all_enhancernotpromoternotctcf.bed'



