#!/usr/bin/bash

basepath='/public/Bio/Bio_project/wgs_snv_detection/'
bam_path='results/samtools/'
result_path='results/'

# move to target dir
cd $basepath$bam_path || exit

# iterate over bam stats files
# make header
res='sample\taverage_quality\tpercentage_of_properly_paired_reads\tBreadth_of_Coverage\n'
for f in  *.stats;
do
    sample="$(basename "$f" .stats)\t"
    # append sample name, percentage if paired reads, acerage quality
    res=$res$sample$(< "$f" grep ^SN | cut -f 2- | grep -E '(percentage of properly paired reads|average quality)' | awk -F '\t' '{print $2}' | paste -sd '\t')"\t"
    # append breadth of coverage
    res=$res$( < "$f" grep ^COV | cut -f 4 | sed '$d' | awk '{OFS="\t"}{c++; if($1 > 0) total+=1}END{print (total/c)*100}')"\n"
done

printf "$res" > $basepath$result_path/bamstats_summary.tsv


