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
    parser.add_argument('-p', '--perturbation',  action='store_true')
    parser.add_argument('-n', '--do_sample',  action='store_true')
    parser.add_argument('-s', '--sample', type=int, default=5000)
    parser.add_argument('-d', '--dists', nargs="*", type=int, default=[50000, 100000])
    parser.add_argument('-t', '--data', type=str, default='tf_ko')
    parser.add_argument('-r', '--rand_state', type=int, default=42)
    parser.add_argument('-o', '--outpath', type=str,
                        default='/mnt/md0/varshini/Analysis/Blood/MicroC_Loops_Anchors/triangle_scores_v2')
    #parser.add_argument('-v', '--verbose', type=bool, action='store_true')
    args = parser.parse_args()

    f = args.filename
    datapath = args.datapath
    subdir = args.subdir
    is_perturbation = args.perturbation
    dists = args.dists
    sample = args.sample
    outpath = args.outpath
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    resolution = 2000

    if args.data=='diff':
        phases = ['Phase1_', 'Phase2_dedup', 'Phase3_dedup']
        microc_cooler_path = '/mnt/coldstorage/Varshmallow/Adipose_Blood_Merged/sankaran'
        microc_cooler_names = [os.path.join(microc_cooler_path, f"Sankaran_{p}merged.250.mcool") for p in phases]

    elif args.data=='tf_ko':
        phases = ['DMSO_P2', 'IKZF1_P2', 'DMSO_P3', 'IKZF1_P3', 'AAVS1_P3', 'NFE2_P3']
        microc_cooler_path = '/mnt/md1/varshini/sankaran_tfko_microc/full_merged_files/'
        microc_cooler_names = [os.path.join(microc_cooler_path, f"{p}_microc.50.mcool") for p in phases]

    elif args.data=='nfe2':
        phases = ['AAVS1_P3', 'NFE2_P3']
        microc_cooler_path = '/mnt/md1/varshini/sankaran_tfko_microc/full_merged_files/'
        microc_cooler_names = [os.path.join(microc_cooler_path, f"{p}_microc.50.mcool") for p in phases]

    elif args.data=='cohesin':
        phases = ['AAVS1_C2_G1', 'STAG1_C2_G1','STAG2_C2_G2','NIPBL_C2_G1',
                  'AAVS1_C3_G1', 'STAG1_C3_G1','STAG2_C3_G2','NIPBL_C3_G1']

        microc_cooler_path = '/mnt/md0/varshini/Analysis/sankaran_cohesin/mcools_per_rep/'
        microc_cooler_names = [os.path.join(microc_cooler_path, f"{p}_all.250.mcool") for p in phases]
        print(microc_cooler_names)
    elif args.data=='u937':
        phases = ['STAG2_KO_dmso',
                  'WT_dmso']

        microc_cooler_path = '/mnt/florence_md0/clarice/stag2_collab/v1_output/full_merged_files/'
        microc_cooler_names = [os.path.join(microc_cooler_path, f"{p}_microc.50.mcool") for p in phases]
    else:
        return 1

    microc_coolers = [utils.get_clr(clr_name, resolution) for clr_name in microc_cooler_names]

    ps_path = '/mnt/md0/varshini/Analysis/Blood/microc_expected'

    PS_names = [f"P_s_{f.split('/')[-1].split('.')[0]}_{resolution}bp_allcol.tsv" for f in microc_cooler_names]
    PS_names = [os.path.join(ps_path, f) for f in PS_names]

    PS_dfs = utils.write_get_ps(zip(PS_names, microc_coolers))
    zipped_clrs = zip(PS_dfs, microc_coolers, phases)

    cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'prob']
    cols.extend([f'strength_{cond}' for cond in phases])
    loopname = '/mnt/md0/varshini/Analysis/Blood/complete_analyses/loops_anchors/merged_calls_125kb_q0.01_classes_new.txt'
    loops = pd.read_csv(loopname, sep='\t')
    loops_formatted = ac.reformat_loops(loops, format='classes', cols_to_append=['pair'])

    if is_perturbation:
        curr_bed = ac.read_peaks(os.path.join(datapath, subdir, f))
    else:
        bedfile = ac.read_df(os.path.join(datapath, subdir, f), sep='\t', columns=['chr', 'coord1', 'coord2'])
        bedfile['perturbation_name'] = f.split('.')[0]
        bedfile['coord'] = bedfile[['coord1', 'coord2']].mean(axis=1)

        if args.do_sample:
            curr_bed = bedfile.sample(sample, random_state=args.rand_state).sort_values(by=['chr', 'coord'])
        else:
            curr_bed = bedfile.sort_values(by=['chr','coord'])

    print("loops")
    print(loops_formatted.head())

    for PS, clr, phase in zipped_clrs:
        for dist in dists:
            print(f"Processing {phase} at {dist} for {f.split('.')[0]}...", file=sys.stderr)
            print(curr_bed.head(), file=sys.stderr)
            fname = f"{phase}_{f.split('.')[0]}_loopscores_{dist}_withpeaks.tsv"
            if not os.path.exists(os.path.join(outpath, fname)):
                out = ac.loop_scores_snip(clr, PS, curr_bed, loops_formatted, dist, quantsize=5000, verbose=False)
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



