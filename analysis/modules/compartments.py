import cooltools
import os
import pandas as pd
import bioframe
import scipy
import numpy as np
import utils

# if cov doesn't exist at res, make it
def get_gc_cov(basedir, cov_name, res, clr=None, fasta_path=None):
    cov_fname = os.path.join(basedir, f"{cov_name}_{int(res/1000)}kb.tsv")
    if not os.path.exists(cov_fname):
        assert ((clr is not None) & (fasta_path is not None)), "these inputs are required to make the gc cov"
        print(f"generating gc cov at {int(res/1000)}kb at location {cov_fname}...")
        genome = bioframe.load_fasta(fasta_path)

        bins = clr.bins()[:]
        gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], genome)
        gc_cov.to_csv(cov_fname, index=False, sep='\t')
    else:
        gc_cov = pd.read_csv(cov_fname, sep='\t')

    return gc_cov

## TODO verify the kwargs
def get_possumm_eig(eig_df_path, eigdir=None, cov=None, res=None, **kwargs):
    if not os.path.exists(eig_df_path):
        assert ((eigdir is not None) & (cov is not None) & (res is not None)), "these inputs are required to make eig df"
        eig_df = process_possumm_eigs(cov, res, eigdir, **kwargs)
        eig_df.to_csv(eig_df_path, index=False, sep='\t')
    else:
        eig_df = pd.read_csv(eig_df_path, sep="\t")

    return eig_df

# function to prep the eigs from the output folder from possumm
# viewframe: regions on which the eigs were computed
def process_possumm_eigs(cov, res, eigdir, viewframe=None, eig_tag="oe_outfile", cov_value_col="GC"):
    if viewframe is None:
        viewframe = utils.get_hg38_arms()
        #viewframe = viewframe[~viewframe["name"].isin(["chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p"])]

    eigs_in = [np.genfromtxt(os.path.join(eigdir, f"{arm['name']}_{eig_tag}")) for _, arm in viewframe.iterrows()]
    ranges = [(reg['chrom'], reg['start'], reg['end']) for _, reg in viewframe.iterrows()]

    # gc sync them
    eigs_corr, _, _ = sync_with_gc(eigs_in, cov, ranges, value_col=cov_value_col)

    # use the view frame to create chrom start end intervals and stitch together
    start_ints = []
    end_ints = []
    chroms = []
    eigs = []

    for eig, range_tup in zip(eigs_corr, ranges):
        # arange is half-open
        # q arm will be truncated for the eig. so, the end should be the last multiple of res.
        end_rounded = range_tup[2] - (range_tup[2] % res)
        start_intervals = np.arange(range_tup[1], end_rounded, res)
        end_intervals = start_intervals + res

        assert (len(start_intervals) == len(eig) & len(end_intervals) == len(eig)), "Error with lengths"

        start_ints.extend(start_intervals)
        end_ints.extend(end_intervals)
        eigs.extend(eig)
        chroms.extend(np.repeat(range_tup[0], len(start_intervals)))

    genome_eig_df = pd.DataFrame({"chrom": chroms,
                                  "start": start_ints,
                                  "end": end_ints,
                                  "E1": eigs})

    return genome_eig_df

def get_eigen_decomp(basedir, eig_name, gc_cov, clr, res, idx=1):
    eig_out = os.path.join(basedir, f"{eig_name}_eig{idx}_{int(res/1000)}kb.tsv")
    view_df = pd.DataFrame({'chrom': clr.chromnames,
                            'start': 0,
                            'end': clr.chromsizes.values,
                            'name': clr.chromnames}
                           )

    if not os.path.exists(eig_out):
        print(f"generating eig {idx} at {int(res/1000)}kb at location {eig_out}...")
        cis_eigs = cooltools.eigs_cis(
            clr,
            gc_cov,
            view_df=view_df,
            n_eigs=idx,
        )

        cis_eigs[1].to_csv(eig_out, sep='\t', index=False)
        eigen_track = cis_eigs[1]

    else:
        eigen_track = pd.read_csv(eig_out, sep='\t')

    return eigen_track


# it's unclear to me whether gc content syncs with hetero/euchromatin at smaller binsizes
# so, use either the gc cov or the eigs themselves from 100kb to sign-adjust an eig computed at a smaller binsize
def sync_with_gc(eigs_to_sync, sync_df, ranges, value_col='GC'):
    corrs = []
    binned_eigs = []
    subs = []
    for eig_to_sync, range_tup in zip(eigs_to_sync, ranges):
        chrom, start, end = range_tup

        sync_vector_sub = sync_df[(sync_df["chrom"] == chrom)
                                  & (sync_df["start"] >= start)
                                  & (sync_df["end"] <= end)][
            value_col].values  # grab portion that corresponds to chrom start end

        eig_fin = eig_to_sync[np.isfinite(eig_to_sync)]
        int_div = int(np.floor(len(eig_fin) / len(sync_vector_sub)))
        eig_fin = eig_fin[:(len(sync_vector_sub) * int_div)]

        binned_eig = np.average(eig_fin.reshape(len(sync_vector_sub), int_div), axis=1)

        binned_eigs.append(binned_eig)

        # from here, referenced from cooltools eigdecomp.py phase_eigs_reference()
        mask = np.isfinite(binned_eig) & np.isfinite(sync_vector_sub)
        corr = scipy.stats.spearmanr(sync_vector_sub[mask], binned_eig[mask])[0]
        corrs.append(corr)
        subs.append(sync_vector_sub[mask])

    for i in range(len(eigs_to_sync)):
        eigs_to_sync[i] = np.sign(corrs[i]) * eigs_to_sync[i]

    return eigs_to_sync, binned_eigs, subs


