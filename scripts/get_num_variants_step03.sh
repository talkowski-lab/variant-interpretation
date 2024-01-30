#!/usr/bin/bash

file=$1;  # cohort tsv file, downloaded from Terra data table

mkdir -p num_vars_cohorts

split_trio_vcf_col_num=$(awk -v name='split_trio_vcfs' '{for (i=1;i<=NF;i++) if ($i==name) print i; exit}' $file)
split_trio_vcf_arr=( $( cat $file | cut -f$split_trio_vcf_col_num | tail -n +2 | tr '\n' ' ' ) )

cohort_col_num=$(awk -v name='entity:cohort_id' '{for (i=1;i<=NF;i++) if ($i==name) print i; exit}' $file)
cohort_arr=( $( awk -F"\t" -v col=$split_trio_vcf_col_num '{if ($col!="") print ;}' OFS="\t" $file | cut -f$cohort_col_num | tail -n +2 | tr '\n' ' ' ) )

cohort_n=0
for cohort in ${cohort_arr[@]};
do
    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
    if (( $cohort_n > ${#split_trio_vcf_arr[@]})); then
        break
    fi
    split_trios=${split_trio_vcf_arr[$cohort_n]}
    cohort_n=$(($cohort_n+1))
    if [[ -z "$split_trios" ]]; then
        continue
    fi
    split_trios_arr=( $( echo $split_trios | jq -c '.[]' | xargs ) )
    num_vars_arr=()
    parent_path=$(dirname ${split_trios_arr[0]})
    vcf_list=$(gsutil -m ls $parent_path | grep -v cromwell_glob_control_file | sort -u)
    for vcf_file in $vcf_list;
        do
            num_vars_arr+=( $( bcftools view -H $vcf_file | wc -l ) )
        done
    num_vars=$( IFS=$'\n'; echo "${num_vars_arr[*]}" )
    paste <(echo "$vcf_list") <(echo "$num_vars") > num_vars_cohorts/"$cohort"_split_trio_num_vars.tsv
done
