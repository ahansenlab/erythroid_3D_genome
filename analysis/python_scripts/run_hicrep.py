import argparse
import utils
import os
import itertools
from hicreppy.hicrep import h_train, genome_scc
import time
import numpy as np
from multiprocessing import Pool
import warnings
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": False})
plt.rcParams['svg.fonttype'] = 'none'
plt.rc('pdf',fonttype = 42)
plt.rcParams["ps.useafm"] = True

def process_with_progress(arguments, threads):
    results = []
    with Pool(processes=threads) as pool:
        # tqdm doesn't really work here, but that's okay
        for result in pool.imap_unordered(scccompute, arguments):
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('coolers_dir')
    parser.add_argument('outpath')
    parser.add_argument('tag')
    parser.add_argument('-d', '--max_dist',  type=int, default=10_000_000)
    parser.add_argument('-r', '--res', type=int,  default=50_000)
    parser.add_argument('-p', '--nproc', type=int, default=24)
    parser.add_argument('-m', '--h_max', type=int, default=10)
    parser.add_argument('-w', '--whitelist', nargs="*", type=str)
    parser.add_argument('-i', '--train_idxs', nargs=2, type=int, default=[0,1])
    parser.add_argument('-o', '--h_opt', type=int)
    parser.add_argument('--train', action='store_true')

    args = parser.parse_args()

    cooler_files = os.listdir(args.coolers_dir)
    cooler_names = [os.path.join(args.coolers_dir, c) for c in
                    cooler_files if 'mcool' in c]
    coolers = [utils.get_clr(clr_name, args.res) for clr_name in cooler_names]

    cooler_sums = [cooler_file.info["sum"] for cooler_file in coolers]
    downsampling_value = min(cooler_sums)

    k = len(coolers)
    pairs = list(itertools.combinations(range(len(coolers)), 2))
    n = len(pairs)

    if not args.whitelist:
        whitelist = utils.generate_all_chroms()
    else:
        whitelist = args.whitelist
    print(whitelist)
    if args.train or not args.h_opt:
        print("Training...")
        start_time = time.time()
        optimal_h = h_train(coolers[args.train_idxs[0]], coolers[args.train_idxs[1]], args.max_dist,
                            args.h_max, args.whitelist)
        end_time = time.time()
        print(f'Done training! ({int(end_time - start_time)/3600} hour)')
        h_opt = optimal_h

    else:
        h_opt = args.h_opt
        print(f"Skipping training, using provided optimal h {h_opt}...")

    arguments = list(
        zip(pairs, n * [coolers], n * [args.max_dist], n * [h_opt], n * [downsampling_value], n * [whitelist]))

    print(f"Running with optimal h {h_opt}...")
    start_time = time.time()
    results = process_with_progress(arguments, args.nproc)
    end_time = time.time()
    print(f'Done running! ({int(end_time - start_time) / 3600} hours)')

    scc_mat = np.zeros((k, k))

    for result in results:
        pair = result[1]
        scc = result[0]
        scc_mat[pair[0], pair[1]] = scc
        scc_mat[pair[1], pair[0]] = scc

    # Also fill the diagonal with 1
    for i in range(len(coolers)):
        scc_mat[i, i] = 1
    # save the scc_mat
    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)

    with open(os.path.join(args.outpath, f"{args.tag}_args.txt"), "w") as file:
        file.writelines("res\tmax_dist\th_opt\tdownsample\twhitelist")
        file.writelines(str(line) + "\t" for line in [args.res, args.max_dist, h_opt, downsampling_value, whitelist])

    np.save(os.path.join(args.outpath, f"{args.tag}_mat.npy"), scc_mat)

    scc_mat_df = pd.DataFrame(scc_mat)
    scc_mat_df.index = [cooler_names[i].split('/')[-1] for i in scc_mat_df.index]
    scc_mat_df.columns = scc_mat_df.index
    sns.clustermap(scc_mat_df, cmap=sns.color_palette("RdBu_r", as_cmap=True), annot=True)
