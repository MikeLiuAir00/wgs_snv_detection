import os
import glob
import json

qc_stat_dir = "/public/Bio/Bio_project/wgs_snv_detection/results/fastp/report"
outdir = "/public/Bio/Bio_project/wgs_snv_detection/results/"
# change wd
os.chdir(qc_stat_dir)

# get file list
f_list = glob.glob(r"*.json")
pre_combined_dict = {}
for j in f_list:
    samplename = j.split('.')[0]
    with open(j, 'r') as jfd:
        content = json.load(jfd)
        total_reads_before_qc = content["summary"]["before_filtering"]["total_reads"]
        total_reads_after_qc = content["summary"]["after_filtering"]["total_reads"]
        q30 = content["summary"]["after_filtering"]["q30_rate"]
        gc = content["summary"]["after_filtering"]["gc_content"]
    pre_combined_dict[samplename] = [total_reads_before_qc, total_reads_after_qc, q30, gc]

print(f"Total samples reads: {len(pre_combined_dict)}")

# read qc stat tsv file
with open(f'{outdir}fastp_stat.tsv', 'w') as f:
    f.write("sample\traw_reads\tcleaned_reads\tQ30\tGC%\n")
    for k, v in pre_combined_dict.items():
        f.write(f'{k}\t')
        f.write('\t'.join([str(s) for s in v]))
        f.write('\n')
