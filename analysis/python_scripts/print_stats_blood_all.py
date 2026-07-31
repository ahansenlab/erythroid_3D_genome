import repro
from tabulate import tabulate
import os
# note to self: for informative read fractions, divide by unque_mapped
stats_dir = "/mnt/md0/varshini/Analysis/Blood/complete_analyses/stats"
if not os.path.exists(stats_dir):
    os.mkdir(stats_dir)

# diff, first seq
summary, gd = repro.summary_metrics('/mnt/coldstorage/Varshmallow/Adipose_Blood_Genomics/merged_per_rep',
                                     kw_names=['Sample', 'Phase', 'Rep'], group_kws=['Phase'], filt_names='Sankaran')
print(tabulate(summary, headers = 'keys', tablefmt = 'psql'))
summary.to_csv(os.path.join(stats_dir, "stats_dec_24.tsv"), sep="\t")

# diff, reseq
summary, gd = repro.summary_metrics('/mnt/coldstorage/Varshmallow/varshinishares/Adipose_Blood_Reseq/merged_per_rep/',
                                     kw_names=['Sample', 'Phase', 'Rep'], group_kws=['Phase'])
print(tabulate(summary, headers = 'keys', tablefmt = 'psql'))
summary.to_csv(os.path.join(stats_dir, "stats_mar_25.tsv"), sep="\t")

# diff, merged (only valid stats are for dedup by run, map /unmap is not informative)
summary, gd = repro.summary_metrics('/mnt/coldstorage/Varshmallow/Adipose_Blood_Merged/sankaran/seqrun_repmerges',
                                     kw_names=['Sample', 'Phase', 'Rep'], group_kws=['Phase'])
print(tabulate(summary, headers = 'keys', tablefmt = 'psql'))
summary.to_csv(os.path.join(stats_dir, "stats_merged.tsv"), sep="\t")

# tf-ko
dir = '/mnt/md1/varshini/sankaran_tfko_microc/merged_per_rep'
summary, gd = repro.summary_metrics(dir,
                                    kw_names=['Sample', 'Condition', 'Rep'], group_kws=['Rep'])
print(tabulate(summary, headers = 'keys', tablefmt = 'psql'))
summary.to_csv(os.path.join(stats_dir, "stats_nov_25.tsv"), sep="\t")

#cohesin
dir = '/mnt/md1/varshini/sankaran_cohesin_microc/merged_per_rep'
summary, gd = repro.summary_metrics(dir,
                                    kw_names=['Sample', 'Condition', 'Guide'])
print(tabulate(summary, headers = 'keys', tablefmt = 'psql'))
summary.to_csv(os.path.join(stats_dir, "stats_jan_26.tsv"), sep="\t")