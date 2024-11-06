#!/usr/bin/bash

file=$1;  # cohort tsv file, downloaded from Terra data table

mkdir -p num_vars_cohorts
trio_denovo_vcf_col_num=$(awk -v name='trio_denovo_vcf' '{for (i=1;i<=NF;i++) if ($i==name) print i; exit}' $file)
trio_denovo_vcf_arr=( $(cat $file | cut -f$trio_denovo_vcf_col_num | tail -n +2 | tr '\n' ' ' ) )

cohort_col_num=$(awk -v name='entity:cohort_id' '{for (i=1;i<=NF;i++) if ($i==name) print i; exit}' $file)
cohort_arr=( $( awk -F"\t" -v col=$trio_denovo_vcf_col_num '{if ($col!="") print ;}' OFS="\t" $file | cut -f$cohort_col_num | tail -n +2 | tr '\n' ' ' ) )

cohort_n=0
for cohort in ${cohort_arr[@]};
do
    echo 'cohort: '$cohort
    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
    if (( $cohort_n > ${#trio_denovo_vcf_arr[@]})); then
        break
    fi
    trio_denovos=${trio_denovo_vcf_arr[$cohort_n]}
    cohort_n=$(($cohort_n+1))
    # if [ "$cohort" != "Ware_LateralityBirthDefects" ]; then
    #     continue
    # fi
    if [[ -z "$trio_denovos" ]]; then
        continue
    fi
    trio_denovos_arr=( $( echo $trio_denovos | jq -c '.[]' | xargs ) )
    parent_path=$(dirname ${trio_denovos_arr[0]})
    vcf_list=$(gsutil -m ls $parent_path | grep -v cromwell_glob_control_file | sort -u)
    num_vars_arr=()
    for vcf_file in $vcf_list;
        do
            # header_size=$(bcftools head $vcf_file | wc -l)
            # num_vars_arr+=( $(( $( gsutil cat $vcf_file | zcat | wc -l ) - $header_size )) )
            num_vars_arr+=( $( gsutil cat $vcf_file | zcat | grep -v '#' | wc -l ) )
        done
    num_vars=$( IFS=$'\n'; echo "${num_vars_arr[*]}" )
    paste <(echo "$vcf_list") <(echo "$num_vars") > num_vars_cohorts/"$cohort"_trio_denovo_num_vars.tsv
done
