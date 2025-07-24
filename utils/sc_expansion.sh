#!/bin/bash
declare -A PARAMS
OPTIND=1
usage() {
    echo -e "
This utility script is developed to aid in bacteria discovery in single cell RNA-sequencing
Usage:
sc_expansion.sh [-s SAMPLE] [-r REFERENCE] [-1 FQ1_PATH] [-2 FQ2_PATH] [-o OUT_DIR] [-t THREADS] [-h]
    -s SAMPLE      the name of the input sample
    -r REFERENCE   path of the hisat2 reference
    -1 FQ1_PATH    path of the R1 fastq file
    -2 FQ2_PATH    path of the R2 fastq2 file
    -o OUT_DIR     output directory
    -t THREADS     number of threads
    -h             show usage
    " 1>&2
}

while getopts ":s:r:1:2:o:t:h" opt; do
    case ${opt} in
        s) PARAMS[SAMPLE]="${OPTARG}" ;;
        r) PARAMS[REF]="${OPTARG}" ;;
        1) PARAMS[FQ1]="${OPTARG}" ;;
        2) PARAMS[FQ2]="${OPTARG}" ;;
        o) PARAMS[OUT]="${OPTARG}" ;;
        t) PARAMS[THREADS]="${OPTARG}" ;;
        h)
            usage
            exit 0
            ;;
        :) echo "Error: -${OPTARG} requires an argument" >&2; return 1 ;;
        \?) echo "Error: Invalid option -${OPTARG}" >&2; return 1 ;;
    esac
done

fq1_new="${output_dir}/${PARAMS[SAMPLE]}_removed_R1.fastq"
fq2_new="${output_dir}/${PARAMS[SAMPLE]}_removed_R2.fastq"

awk '
{
    if (NR % 4 == 1) {
        header = $0
    } else if (NR % 4 == 2) {
        sequence = $0
        prefix = substr(sequence, 1, 16)
        remaining = substr(sequence, 17)
        print header "." prefix
        print remaining
    } else if (NR % 4 == 3) {
        print $0
    } else if (NR % 4 == 0) {
        quality = $0
        remaining_quality = substr(quality, 17)
        print remaining_quality
    }
}' "${PARAMS[FQ1]}" > "${fq1_new}"

awk '
{
    if (NR % 4 == 1) {
        header = $0
    } else if (NR % 4 == 2) {
        sequence = $0
        prefix = substr(sequence, 1, 16)
        remaining = substr(sequence, 17)
        print header "." prefix
        print remaining
    } else if (NR % 4 == 3) {
        print $0
    } else if (NR % 4 == 0) {
        quality = $0
        remaining_quality = substr(quality, 17)
        print remaining_quality
    }
}' "${PARAMS[FQ2]}" > "${fq2_new}"

out_bacc="${PARAMS[OUT]}/${PARAMS[SAMPLE]}"
mkdir -p "${out_bacc}"
bacNeo --bacc -1 "${fq1_new}" -2 "${fq2_new}" -m RNA -r "${PARAMS[REF]}" -o "${out_bacc}" -l f -l g -t "${PARAMS[THREADS]}"
