#!/usr/bin/env bash

# extract read group id

## specification
### ID: grep -o -E ":[0-9a-zA-Z_]{9}:[0-9]{1}"
### SM: $sample
### PL: illumina
### LB: treat individual sample as LB

## two kinds of seq header exists
### ordinary seq header(illumica spec)
### uncommon format
#### grep -E '^@[VE][[:digit:]]{9}L[[:digit:]]C[[:digit:]]{3}R[[:digit:]]{3}[[:digit:]]{7}/[[:digit:]]'

shopt -s lastpipe
batch1='/public/Bio/Bio_project/wgs_snv_detection/assets/sample_batch1.txt'
batch2='/public/Bio/Bio_project/wgs_snv_detection/assets/sample_batch2.txt'
cleaned_path='/public/Bio/Bio_data/Armeniaca_Resequencing/JZP202302GM011-01/cleaned/'
outdir='/public/Bio/Bio_project/wgs_snv_detection/assets/'

cd $cleaned_path || exit

# init output string
out_str="RGID,RGLB,RGPL,RGPU,RGSM\n"
for f in $(ls *1.fastp.fastq.gz | sort -V);
do
    # extract sample name
    sample=$(echo "$f" | grep -o -E "PS[[:digit:]]{1,2}_[[:digit:]]{1,2}")
    is_ord=true

    # set flag fastq format
    if [[ -n $(zcat "$f" | head -n 1| grep -E '^@[VE][[:digit:]]{9}L[[:digit:]]C[[:digit:]]{3}R[[:digit:]]{3}[[:digit:]]{7}/[[:digit:]]') ]];
    then
        is_ord=false
    fi


    read_id=""
    # loop over headers
    zcat "${cleaned_path}$f" | grep '^@' | while read -r line;
    do
        # get RGID
        case "$is_ord" in
        true)
            read_id="$(echo "$line" | sed -n 's/.*\([0-9a-zA-Z_]\{9\}\):\([0-9]\{1\}\).*/\1\.\2/p')"
        ;;
        false)
            # uncommon fastq spec
            read_id="$(echo "$line" | sed -n 's/^@[VE][[:digit:]]\{5\}\([[:digit:]]\{4\}\)L\([[:digit:]]\).*/\1\.\2/p')"
        ;;
        *)
            exit
        esac
        # get only the metadata from first read, since multile run in sample observed
        break
    done

    ## determine batch info
    library=''
    if [[ -n $(grep $sample $batch1) ]];then
        # batch 1 sequenced data
        library='lib1'
    elif [[ -n $(grep $sample $batch2) ]];then
        # batch 2 sequence data from NC
        library='lib2'
    else
        echo "something wrong"
        exit
    fi

    # compose read group info per sample
    # RGID,RGLB,RGPL,RGPU,RGSM
    out_str="${out_str}${read_id},${library},illumina,${read_id},${sample}_T1\n"
done

printf "$out_str" > ${outdir}readgroups.csv
