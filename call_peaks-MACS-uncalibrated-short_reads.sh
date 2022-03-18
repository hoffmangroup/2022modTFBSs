#!/usr/bin/env bash

set -o pipefail -o errexit

ERR_EXIT=64

# Usage: ./call_peaks-MACS.sh <input BED directory> [output directory]
# Calls peaks, with MACS 2 using a negative (IgG) set.
#
# Does not perform spike-in normalization.
#
# Only uses fragments less than 120 bp.
#

export SHELL=$(type -p bash)

SCRIPT_DIR="$(dirname $(readlink -f ${BASH_SOURCE[0]}))"

# ----------------------------------------------------
# set and propagate expected env. variables

# slots

DEFAULT_SLOTS=1

NSLOTS=${NSLOTS:-$DEFAULT_SLOTS}

SLURM_CPUS_ON_NODE=${SLURM_CPUS_ON_NODE:-$DEFAULT_SLOTS}
SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-$SLURM_CPUS_ON_NODE}

export OMP_NUM_THREADS=$(( $SLURM_CPUS_PER_TASK > $NSLOTS ? $SLURM_CPUS_PER_TASK : $NSLOTS ))

SLOTS_PER_PAR=$(($OMP_NUM_THREADS / 4))
SLOTS_PER_SORT=$(($OMP_NUM_THREADS - ($OMP_NUM_THREADS / 4)))

# N.B. expects at least 20 GB of memory; not monitored
MAX_MEM='20GB'

# temp. dir
TMPDIR="$SLURM_TMPDIR"
TMP="$TMPDIR"
TEMP="$TMPDIR"
export SLURM_TMPDIR TMPDIR TMP TEMP
# ----------------------------------------------------
# enable env_parallel.bash
. $(which env_parallel.bash)
# ---------------------------------------------------- 

MACS_CALLER='macs2'

MACS_GENOME_SIZE='mm'

MACS_FDR=0.05

PEAK_THRESH=20

MAX_FRAG_LEN=120

NEG_SET='IgG'

INPUT_FILE_EXT='.bam'

input_dir="${1:-}"

output_dir="${2:-$input_dir}"


function rm_temp {
      rm -Rf "$TEMP_BASE_DIR"
}
trap rm_temp EXIT


function format_name {
    basename "${1:-}" $INPUT_FILE_EXT | sed -r 's/(-[12]+)[_-]S[[:digit:]]+-/\1/; s/(Mus_musculus.mm10|genome)//; s/\.sort//'
}
export -f format_name


function call_peaks {
    pos_dataset=${1:-}

    neg_dataset=${2:-}

    dataset_name="$(format_name "$pos_dataset")-$(format_name "$neg_dataset")"

    call_peaks_MACS
}
export -f call_peaks


function call_peaks_pooled {
    # N.B. exactly two of each type:
    #      pos1 pos2 neg1 neg2

    pos_datasets=("${1:-}" "${2:-}")

    neg_datasets=("${3:-}" "${4:-}")

    dataset_name="$(format_name "${1:-}" | sed -r 's/-[12]//')-$(format_name "${3:-}" | sed -r 's/-[12]//')"

    call_peaks_MACS
}
export -f call_peaks_pooled


function call_peaks_MACS {
    if [[ -n "$disable_MACS" ]]; then
        >&2 echo "Skipped MACS, as requested."
        return
    fi

    pos_datasets=("${pos_datasets[@]:-$pos_dataset}")
    neg_datasets=("${neg_datasets[@]:-$neg_dataset}")

    # * N.B. MACS instead uses the actual BAMPE output from Bowtie *
    # assumes that the BAMs have the same file naming structure
    BAM_name_portion='.bam'

    if [[ ! -f ${pos_datasets[0]/$INPUT_FILE_EXT/$BAM_name_portion} ]]; then
        BAM_name_portion='.sort.bam'

        if [[ ! -f ${pos_datasets[0]/$INPUT_FILE_EXT/$BAM_name_portion} ]]; then
               >&2 echo "Cannot find BAM files."
               exit $ERR_EXIT
        fi
    fi

    TEMP_BASE_DIR=$(mktemp -d)

    pos_datasets_TEMP="$TEMP_BASE_DIR/$(basename ${pos_datasets[0]})-filter.bam"
    neg_datasets_TEMP="$TEMP_BASE_DIR/$(basename ${neg_datasets[0]})-filter.bam"

    for TEMP_array_name in 'pos_datasets' 'neg_datasets'; do
        TEMP_array=${!TEMP_array_name}

        output_file_TEMP="$TEMP_BASE_DIR/$(basename ${TEMP_array[0]})-filter.bam"

        if [[ ${#TEMP_array[@]} -gt 1 ]]; then
            filter_input_file="${output_file_TEMP##*.}-merged.bam"
            output_file_TEMP="${output_file_TEMP##.*}-merged.bam"
            pos_datasets_TEMP="${pos_datasets_TEMP##.*}-merged.bam"
            neg_datasets_TEMP="${neg_datasets_TEMP##.*}-merged.bam"

            if [[ ! -f "$filter_input_file" ]]; then
                # first need to merge, to permit filtering
                sambamba merge -t $OMP_NUM_THREADS -l 0 "$filter_input_file" ${TEMP_array[@]}
            fi
        else
            filter_input_file="${TEMP_array[0]}"
        fi

        if [[ ! -f "$output_file_TEMP" ]]; then
            alignmentSieve -b "$filter_input_file" -o "$output_file_TEMP" --maxFragmentLength $MAX_FRAG_LEN
        fi
    done

    $MACS_CALLER callpeak --buffer-size 1000000 -t ${pos_datasets_TEMP[@]/$INPUT_FILE_EXT/$BAM_name_portion} -c ${neg_datasets_TEMP[@]/$INPUT_FILE_EXT/$BAM_name_portion} -f 'BAMPE' --gsize "$MACS_GENOME_SIZE" --outdir "$output_dir" -n "$dataset_name" --qvalue $MACS_FDR --call-summits --bdg --SPMR
}
export -f call_peaks_MACS


# Our negative controls are matched to a particular replicate; use only the matched controls.
# Could also test cross-use of each negative control with the other set of biological replicates,
# which could be used to reduce batch-effects, but not doing so for now.
#
# Do not try to call peaks on the spike-in data; skip all those files.
#
# N.B. we expect biological replicates to be denoted EXACTLY, by either "[-_]1[-_]" or "[-_]2[-_]".
#

# assumes equal number of replicates per set, but fewer in negative set than positive

# N.B. nullglob was making debugging harder; no longer needed
shopt -s extglob

pos_samples_rep_1=($input_dir/!(*$NEG_SET*)+(-|_)1+(-|_)!(*$NEG_SET*|*spike*)$INPUT_FILE_EXT)
pos_samples_rep_2=($input_dir/!(*$NEG_SET*)+(-|_)2+(-|_)!(*$NEG_SET*|*spike*)$INPUT_FILE_EXT)

neg_samples_rep_1=($input_dir/*$NEG_SET+(-|_)1+(-|_)!(*spike*)$INPUT_FILE_EXT)
neg_samples_rep_2=($input_dir/*$NEG_SET+(-|_)2+(-|_)!(*spike*)$INPUT_FILE_EXT)

if [[ ! -s ${neg_samples_rep_1[0]} ]]; then
    >&2 echo "WARNING: no first IgG replicate detected. Attempting to use single $NEG_SET to use as replicates for both."

    neg_samples_rep_1=()
    neg_samples_rep_2=()

    # "second" uses the exact same file as first
    neg_samples_rep_1=($input_dir/+(${NEG_SET})!(*spike*)$INPUT_FILE_EXT)

    neg_samples_rep_2=("${neg_samples_rep_1[@]}")

    if [[ -s "${neg_samples_rep_1[0]}" ]]; then
        >&2 echo "WARNING: using ${neg_samples_rep_1[0]} as both first and second negative sets (i.e. double-counting single or pre-pooled $NEG_SET for both replicates)."
    else
        >&2 echo "Unable to find any $NEG_SET. At least one control dataset is required for this script."
        exit $ERR_EXIT
   fi
fi

num_neg_samples_per_rep=${#neg_samples_rep_1[@]}

if [[ ${#pos_samples_rep_1[@]} -ne $num_neg_samples_per_rep || ${#pos_samples_rep_1[@]} -ne ${#pos_samples_rep_2[@]} || ${#neg_samples_rep_1[@]} -ne ${#neg_samples_rep_2[@]} ]]; then
    >&2 echo "WARNING: mis-matched datasets."
    >&2 echo -e "$num_neg_samples_per_rep\n\n${pos_samples_rep_1[@]} ${pos_samples_rep_2[@]} ${neg_samples_rep_1[@]+${neg_samples_rep_1[@]}} ${neg_samples_rep_2[@]+${neg_samples_rep_2[@]}}"
fi

if [[ $num_neg_samples_per_rep -eq 0 ]]; then
   >&2 echo "At least one control dataset is required for this script."
   exit $ERR_EXIT
fi

if [[ $num_neg_samples_per_rep -gt 1 ]]; then
   >&2 echo "More than one control per replicate set not supported."
   exit $ERR_EXIT
fi


# use env_parallel, to prevent issues with exporting the functions (in some GNU Parallel versions)

# separate biological replicates.
(env_parallel --env _ --linebuffer -k --rpl '{sample} s:.*/([^/]+[-_][12])[-_].*:\1:' --tagstring "{1sample}:" --joblog "$output_dir/call_peaks-sep-jobs.log" --resume --resume-failed --link "call_peaks {1} {2}" ::: ${pos_samples_rep_1[@]} ${pos_samples_rep_2[@]} ::: ${neg_samples_rep_1[@]} ${neg_samples_rep_2[@]}) || true &
PS1=$!

# same as above, but now pooling both replicates from each dataset (including pooling the background; MACS only)
(env_parallel --env _ --linebuffer -k --rpl '{sample} s:.*/([^/]+[-_][12])[-_].*:\1:' --tagstring "{1sample}-{3sample} (pooled):" --joblog "$output_dir/call_peaks-pooled-jobs.log" --resume --resume-failed --link "call_peaks_pooled {1} {3} {2} {4}" ::: ${pos_samples_rep_1[@]} ::: ${neg_samples_rep_1[@]} ::: ${pos_samples_rep_2[@]} ::: ${neg_samples_rep_2[@]}) || true &
PS2=$!

wait $PS1 $PS2
