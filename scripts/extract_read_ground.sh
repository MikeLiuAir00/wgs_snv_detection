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
cleaned_path='/public/Bio/Bio_data/Armeniaca_Resequencing/JZP202302GM011-01/cleaned/'
out_dir='/public/Bio/Bio_project/wgs_snv_detection/assets/read_groups/'

cd $cleaned_path || exit
for f in *.fastp.fastq.gz;
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
    ## indicate sample
    echo "$sample"
    zcat "${cleaned_path}$f" | grep '^@' | while read -r line;
    do
        # get RGID
        case "$is_ord" in
        true)
            read_id="$(echo "$line" | sed -n 's/.*:\([0-9]\{1,\}:[0-9a-zA-Z_]\{9\}\):\([0-9]\{1\}\).*/\1\.\2/p')"
        ;;
        false)
            # uncommon fastq spec
            read_id="$(echo "$line" | sed -n 's/^@[VE][[:digit:]]\{5\}\([[:digit:]]\{4\}\)L\([[:digit:]]\).*/\1\.\2/p')"
        ;;
        *)
            echo "None"
        esac
    done
    echo "$read_id"
done

