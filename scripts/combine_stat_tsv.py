import os
import numpy as np
import pandas as pd

wd = "/public/Bio/Bio_project/wgs_snv_detection/results"

# change wd
os.chdir(wd)

# read files
fastp_stat = pd.read_csv(f'{wd}/fastp_stat.tsv', sep='\t', header=0, index_col=0)
bam_stat = pd.read_csv(f'{wd}/bamstats_summary.tsv', sep='\t', header=0, index_col=0)

print(f"pre-merge: {len(fastp_stat.index)}, {len(bam_stat.index)}")

merged = pd.concat([fastp_stat, bam_stat], axis=1)
print(f"merged: {len(merged.index)}")
print(merged.head)

merged.to_csv(f'{wd}/alignment_stat_summary.csv', index=True)
