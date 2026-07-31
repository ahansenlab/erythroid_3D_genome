from setup_logger import logger
import os
import utils
from tabulate import tabulate
# modules
import calc_insulation
import repro
import pandas as pd

def reproducibility(cooler_names, resolution, sd=None):
    coolers = [utils.get_clr(clr_name, resolution) for clr_name in cooler_names]
    scc = repro.compute_repro(coolers, ['chr2', 'chr6', 'chr11'])
    repro.plot_repro(scc, sd)


def insulation_microc(cooler_names, resolution, regions=None, windows=None):
    coolers = [utils.get_clr(clr_name, resolution) for clr_name in cooler_names]

    if windows is None:
        windows = [3 * resolution, 5 * resolution, 10 * resolution, 25 * resolution, 20000]

    ins_tables = calc_insulation.compute_insulation(coolers, windows, regions)
    return ins_tables


def main():
    # phases = ['Phase1_', 'Phase2_dedup', 'Phase3_dedup']
    #
    # microc_cooler_path = '/mnt/coldstorage/Varshmallow/Adipose_Blood_Merged/sankaran'
    # microc_cooler_names = [os.path.join(microc_cooler_path, f"Sankaran_{p}merged.250.mcool") for p in phases]
    #
    # output_path = '/mnt/md0/varshini/Analysis/Blood/Boundaries_Merged/'
    # if not os.path.exists(output_path):
    #     os.mkdir(output_path)
    #
    # ## insulation
    # res=1000
    # ins_tables = insulation_microc(microc_cooler_names, res, windows=[5000, 10000, 20000])
    # for ins_table, p in zip(ins_tables, phases):
    #     ins_table.to_csv(os.path.join(output_path, f'{p}_microc_insulation_{res}bp.tsv'), index=False,
    #                                sep='\t')

    conds = ['DMSO_P2', 'IKZF1_P2', 'DMSO_P3', 'IKZF1_P3', 'AAVS1_P3', 'NFE2_P3', ]
    res = 1000

    microc_cooler_path = '/mnt/md1/varshini/sankaran_tfko_microc/full_merged_files'
    microc_cooler_names = [os.path.join(microc_cooler_path, f'{cond}_microc.50.mcool') for cond in conds]
    output_path = '/mnt/md0/varshini/Analysis/Blood/matchmakers/'

    ins_tables = insulation_microc(microc_cooler_names, res, windows=[5000, 10000, 20000])
    for ins_table, p in zip(ins_tables, conds):
         ins_table.to_csv(os.path.join(output_path, f'{p}_microc_insulation_{res}bp.tsv'), index=False,
                                    sep='\t')

    output_path = '/mnt/md0/varshini/Analysis/Blood/matchmakers/'
    # ins_table_merged = insulation_microc([os.path.join(microc_cooler_path,
    #                                                    'Sankaran_P123_final.merged.mcool')], res)
    #
    # ins_table_merged[0].to_csv(os.path.join(output_path, f'merged_microc_insulation_{res}.tsv'), index=False,
    #                           sep='\t')

    # id_table = pd.read_csv(os.path.join(output_path, f'merged_microc_insulation_{res}.tsv'), sep='\t')
    # match_tables=[]
    # for p in phases:
    #     match_tables.append(pd.read_csv(os.path.join(output_path, f'{p}_microc_insulation_{res}bp.tsv'), sep='\t'))
    #
    # matched_boundaries = calc_insulation.boundary_matching(id_table, match_tables, phases, window=20000, merge_method='left')
    # matched_boundaries.to_csv(os.path.join(output_path, 'merged_boundary_strengths_left.tsv'), index=False, sep='\t')

if __name__ == '__main__':
    main()