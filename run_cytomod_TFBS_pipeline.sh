#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

ERR_EXIT=64

NSLOTS=${NSLOTS:-1}

CLUSTER_RESERVATION=true
cluster_res_param='n'
if [[ "$CLUSTER_RESERVATION" == true ]]; then
    cluster_res_param='y'
fi

MEME_ChIP_USER_HOLD=true
MEME_ChIP_user_hold_param=''
if [[ "$MEME_ChIP_USER_HOLD"  == true ]]; then
    MEME_ChIP_user_hold_param='-h'
fi

PE_TYPE='smp'
NUM_SORT_CORES=8
LOWER_MEM_REQ='32G'
OTHER_MEM_REQ='128G'
OTHER_EST_TIME='128:00:00'

DEFAULT_NUM_PROCS_PER_MEME_CHIP_RUN=2

MPI_ENV_ARG='-pe ompi'
num_procs_per_MEME_ChIP_run=${4:-$DEFAULT_NUM_PROCS_PER_MEME_CHIP_RUN}

par_env_arg=""
par_MEME_arg=""
if [ "$num_procs_per_MEME_ChIP_run" -gt 1 ]; then
    par_env_arg="$MPI_ENV_ARG $num_procs_per_MEME_ChIP_run"
    par_MEME_arg="-meme-p $num_procs_per_MEME_ChIP_run"
fi

IMPORTANT_FORMAT=$'\e[1m'
ERROR_FORMAT=$'[31m'
RESET_FORMAT=$'\e[0m'

# Usage: ./run_cytomod_TFBS_pipeline.sh <filename of list of BAMs> <directory containing summit BEDs to use> [output directory base=results] [num procs per MEME-ChIP run=2]
# List of BAM file format (one line per sample; tab-delimited): <WGBS>	[oxWGBS]	[TAB-seq]
# narrowPeaks directory format: <directory containing narrowPeak files (*.narrowPeak or *.narrowPeak.gz)>

# ---------------------------------- Create modified genomes ----------------------------------

DIR_PREFIX=${3:-results}

BAM_list=${1:-}

summit_BED_files_dir=${2:-}

# simulate lastpipe
readarray -t WGBS_files < <(awk '{print $1;}' "$BAM_list")
readarray -t oxWGBS_files < <(awk '{print $2;}' "$BAM_list")
readarray -t TABSeq_files < <(awk '{print $3;}' "$BAM_list")

JOB_LOG_FILE='jobs-final.log'

# removed for portability
MODULE_LOAD_CMD=''

ASSEMBLY='hg38'

ORGANISM=''
if [[ ${ASSEMBLY:0:1} == 'h' ]]; then
    ORGANISM='human'
elif [[ ${ASSEMBLY:0:1} == 'm' ]]; then
    ORGANISM='mouse'
else
    >&2 echo -e "${ERROR_FORMAT}Unknown assembly.${RESET_FORMAT}"
    exit $ERR_EXIT
fi

CYTOMOD_PATH="$HOME/cytomod/cytomod/src/cytomod.py"

BASE_FREQ_SCRIPT_NAME='makeRegionBaseFreqPlots.sh'

BASE_FREQ_SCRIPT_PATH="$PWD/$BASE_FREQ_SCRIPT_NAME"

if [[ ! -f "$BASE_FREQ_SCRIPT_PATH" && -d "$HOME/cytomod/experiments" ]]; then
    BASE_FREQ_SCRIPT_PATH="$HOME/cytomod/experiments/2014-07-10-ModifiedGenomesToDate/$BASE_FREQ_SCRIPT_NAME"
elif [[ ! -f "$BASE_FREQ_SCRIPT_PATH" ]]; then
    >&2 echo -e "${ERROR_FORMAT}Unable to find base frequency script ($BASE_FREQ_SCRIPT_NAME).${RESET_FORMAT}"
    exit $ERR_EXIT
fi

BASE_FREQ_SCRIPT_CMD_BASE="$BASE_FREQ_SCRIPT_PATH $ORGANISM $ASSEMBLY"

OUTPUT_MOD_BEDGRAPHS_PATH="$PWD/output5mC5hmCModBaseBEDsAtThreshold.py"

RUN_H_TESTS_PATH="$PWD/run_hypothesis_tests.sh"

PLOT_PATH="$PWD/processDataAndPlot.py"

BASE_QSUB_CMD='qsub -cwd -q hoffmangroup -j y -terse'

module load igenome-$ORGANISM/$ASSEMBLY

GENOME_SIZE_FILE="$REF_HOME/Sequence/WholeGenomeFasta/genome.fa.fai"

THRESHOLDS=(0.3 0.7)

EXP_DATA_DIR='processed_data'

CHIP_SEQ_DATA_DIR="$HOME/ENCODE_K562_ChIP-Seq_TFs-untreated_hg38-Mehran_processed/"

CYTOMOD_TRACK_DIR_PREFIX='mod_bedGraphs_'

MLML_SPLIT_INPUT_DIR='split_MLML_in'
ARCHIVES_BASE_DIR='archives'
RESULTS_BASE_DIR='results'
SPLIT_PREFIX='x'

BISMARK_REF="/mnt/work1/data/commondata/bismark/${ASSEMBLY}"

INITIAL_MEME_WORK_DIR="$DIR_PREFIX/MEME_workdir"

EXPANDED_NAME_STRING='expandedTo500bpRegions'

EXP_BED_DIR="$INITIAL_MEME_WORK_DIR/expanded_BEDs"

EDIT_CENTRIMO_SCRIPT="$HOME/cytomod/utils/src/editCentriMoHTML.sh" 

MEME_ChIP_CMDS_LOG='MEME-ChIP_commands_run.log'

# still include 5hmC, since heavily confounded and using x/7 for all bases
ALPHABET_FILE="$HOME/meme/src/unit-tests/data/mod_spec_no5caC_no5fC.txt"

MIN_MEME_MOTIF_WIDTH=7
MAX_MEME_MOTIF_WIDTH=12

MEME_MOTIF_DB="$HOME/meme/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme"

# Some genomes have a "full" version (e.g. mm9), in which the non-full
# version only contains a subset of chrs (e.g. missing unplaced contigs)
# Other genomes without a "full" version are generally already the
# full version (e.g. hg38), but this needs to always be checked.
if [[ -f "${BISMARK_REF}_full" ]]; then
    BISMARK_REF="${BISMARK_REF}_full"
fi

declare -a WGBS_MLML_input oxWGBS_MLML_input TABSeq_MLML_input

METHPIPE_TEMP_DIR="$DIR_PREFIX/methpipe_temp"
mkdir "$METHPIPE_TEMP_DIR"

merge_methcount_jIDs=''
merge_methcount_cmds=''
for BAM_files_array_name in 'WGBS_files' 'oxWGBS_files' 'TABSeq_files'; do
    BAM_files_array_temp_name="${BAM_files_array_name}[@]"
    BAM_files_array="${!BAM_files_array_temp_name}"

    if [[ -n "${BAM_files_array[@]// }" ]]; then
        # NB: sorting MUST be conducted by read name, since to-mr may o/w SILENTLY DISPOSE of arbitrary paired-end mates (refer to MethPipe Issue #71 on GitHub).
        # NB: need intermediate view between sort and to-mr, since forcing SAM input causes header to be read correctly
        tomr_jIDs=''
        tomr_jIDs=($(parallel -j 1 --rpl '{id} s:.*/::; s:(\.[^.]+)+$::;' --joblog "$JOB_LOG_FILE.0" "echo \"$MODULE_LOAD_CMD sambamba sort -l 0 -u -t $NUM_SORT_CORES -m $OTHER_MEM_REQ -N -o /dev/stdout {} | sambamba view -h /dev/stdin | to-mr -v -m 'bismark' -o $METHPIPE_TEMP_DIR/$BAM_files_array_name-{id}.tomr /dev/stdin\" | tee -a /dev/stderr | $(echo "$BASE_QSUB_CMD" | sed -r 's/-j y//') -N 'to-mr-$BAM_files_array_name-{id}' -R $cluster_res_param -pe $PE_TYPE $NUM_SORT_CORES -o \"$METHPIPE_TEMP_DIR/to-mr-$BAM_files_array_name.out\" -e \"$METHPIPE_TEMP_DIR/to-mr-$BAM_files_array_name.err\" -l mem_requested=$(($(echo $OTHER_MEM_REQ | sed 's/[A-Za-z]//g') / $NUM_SORT_CORES))$(echo $OTHER_MEM_REQ | sed -r 's/.*([A-Za-z])$/\1/g'),h_rt=$OTHER_EST_TIME" ::: ${BAM_files_array[@]} ))

        merge_methcount_cmds=$(parallel -j $NSLOTS --joblog "$JOB_LOG_FILE.1" --rpl '{-../} s:.*/::; s:(\.[^.]+)+$::; s:-\d+$::;' --dry-run "echo \"$MODULE_LOAD_CMD export LC_ALL=C; cat $METHPIPE_TEMP_DIR/$BAM_files_array_name-{-../}*.tomr | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 | duplicate-remover -stdin -v -S $METHPIPE_TEMP_DIR/duplicate-remover_stats.out -A -o /dev/stdout | methcounts -v -c $BISMARK_REF -o $METHPIPE_TEMP_DIR/$BAM_files_array_name-{-../}_pooled_ALL.meth /dev/stdin\" | tee -a /dev/stderr | $BASE_QSUB_CMD -hold_jid $(printf ",%s" "${tomr_jIDs[@]}") -R $cluster_res_param -N 'merge->dedup->methcount-$BAM_files_array_name-{-../}' -o \"$METHPIPE_TEMP_DIR/merge_dedup_methcount.log\" -l mem_requested=$OTHER_MEM_REQ,h_rt=$OTHER_EST_TIME" ::: ${BAM_files_array[@]})

        IFS=$'\n'
        for cmd in $merge_methcount_cmds; do
            merge_methcount_jIDs+="$(eval "$cmd"),"

            MLML_input=$(echo "$cmd" | grep -oP "$BAM_files_array_name-\w+_pooled_ALL")
            if [[ "$BAM_files_array_name" == *'oxWGBS'* ]]; then
                oxWGBS_MLML_input+=("-m $MLML_input")
            elif [[ "$BAM_files_array_name" == *'WGBS'* ]]; then
                WGBS_MLML_input+=("-u $MLML_input")
            elif [[ "$BAM_files_array_name" == *'TABSeq_files'* ]]; then
                TABSeq_MLML_input+=("-h $MLML_input")
            fi
        done
        unset IFS
    
        # remove to-mr files
        parallel -j $NSLOTS --rpl '{id} s:.*/::; s:(\.[^.]+)+$::;' "echo \"rm -f $METHPIPE_TEMP_DIR/$BAM_files_array_name-{id}.tomr\" | $BASE_QSUB_CMD -N 'clean_to-mr' -o /dev/null -hold_jid $(echo $merge_methcount_jIDs)" ::: ${BAM_files_array[@]}
    fi
done

# Job IDs to wait for prior to starting the Cytomod portion of the workflow
# This is the merge methcount job IDs, unless we are running MLML, in which case it is those job IDs
start_Cytomod_workflow_hold_jIDs="$merge_methcount_jIDs"

MLML_jIDs=''

# MLML only applies if we are estimating 5mC and 5hmC levels
if [[ ${#oxWGBS_MLML_input[@]} -ne 0 || ${#TABSeq_MLML_input[@]} -ne 0 ]]; then
    for ((i=0; i<${#WGBS_MLML_input[@]}; i++)); do
        file_ID_string_for_mlml=$(echo "$merge_methcount_cmds" | grep -oP '\w+_pooled')
        file_name="$(echo ${WGBS_MLML_input[$i]} | sed 's/-. //')"
        MLML_jIDs+=$(echo "$MODULE_LOAD_CMD mlml ${WGBS_MLML_input[$i]} ${oxWGBS_MLML_input[$i]:-} ${TABSeq_MLML_input[$i]:-} -o '$METHPIPE_TEMP_DIR/$file_name.mlml';" | tee -a /dev/stderr | $BASE_QSUB_CMD -hold_jid $(echo "$merge_methcount_jIDs") -N "MLML-$file_name" -o "$METHPIPE_TEMP_DIR/mlml.log" -l mem_requested=$OTHER_MEM_REQ,h_rt=$OTHER_EST_TIME)
    done
    start_Cytomod_workflow_hold_jIDs="$MLML_jIDs"
fi

# ---------------------------------- Cytomod and TFBS analysis ----------------------------------

for thresh in ${THRESHOLDS[@]}; do
    for ((i=0; i<${#WGBS_MLML_input[@]}; i++)); do
        file_name="$(echo ${WGBS_MLML_input[$i]} | sed 's/-. //')"

        # ---------------------------------- Cytomod ----------------------------------
        
        # if the index file for the repeat masked version is present, assume all Cytomod runs are done and skip
        if [[ ! -f $resultsDir/fullModGenome-TRFmasked.fa.fai ]]; then
            consistent_counts_file=''
            if [[ ${#oxWGBS_MLML_input[@]} -ne 0 || ${#TABSeq_MLML_input[@]} -ne 0 ]]; then
                consistent_counts_file="$METHPIPE_TEMP_DIR/$file_name.mlml"
            else
                consistent_counts_file="$METHPIPE_TEMP_DIR/$file_name.meth"
            fi

            cytomod_processing_dir="$DIR_PREFIX/$CYTOMOD_TRACK_DIR_PREFIX$(echo $file_name | sed 's/_pooled.*//')/"
            
            mkdir -p "$cytomod_processing_dir$thresh"
            
            track_creation_jID=$(echo "python $OUTPUT_MOD_BEDGRAPHS_PATH $consistent_counts_file $thresh -o \"$cytomod_processing_dir$thresh\"" | $BASE_QSUB_CMD  -hold_jid "$start_Cytomod_workflow_hold_jIDs" -N "output5mC5hmCModBaseBEDsAtThreshold_$thresh" -o "$cytomod_processing_dir/output5mC5hmCModBaseBEDsAtThreshold-$thresh.log" -l mem_requested=$LOWER_MEM_REQ)
            archiveDir="$cytomod_processing_dir/$ARCHIVES_BASE_DIR/filesForArchive-$thresh"
            mkdir -p "$archiveDir"
            resultsDir="$cytomod_processing_dir/$RESULTS_BASE_DIR/filesForArchive-$thresh"
            mkdir -p "$resultsDir"

            cytomod_NRM_log="$resultsDir/cytomod_non-RMasked.log"
            
            cytomod_NRM_jID="$($BASE_QSUB_CMD -hold_jid ${track_creation_jID:-0} -V -o $cytomod_NRM_log -l mem_requested=$LOWER_MEM_REQ -b y $CYTOMOD_PATH -v -d $REF_HOME/Sequence/Chromosomes $cytomod_processing_dir$thresh/ --archiveOutDir $archiveDir/ --BEDOutDir $resultsDir -f $resultsDir/fullModGenome.fa)"
            echo "$BASE_FREQ_SCRIPT_CMD_BASE $resultsDir/fullModGenome.fa -o $resultsDir/fullModGenome" | $(echo "$BASE_QSUB_CMD" | sed -r 's/-j y//') -hold_jid ${track_creation_jID:-0},${cytomod_NRM_jID:-0} -N 'baseFreq_mod_genome' -o "$resultsDir/fullModGenome-baseFreqsOut.txt" -e "$cytomod_NRM_log"
            echo "$MODULE_LOAD_CMD samtools faidx $resultsDir/fullModGenome.fa" | $BASE_QSUB_CMD -hold_jid ${track_creation_jID:-0},${cytomod_NRM_jID:-0} -o "$cytomod_NRM_log" -N 'faidx'

            # requires existing TRF-masked genome, linked from CWD
            cytomod_RM_log="$resultsDir/cytomod_RMasked.log"
            cytomod_RM_jID="$($BASE_QSUB_CMD -hold_jid ${track_creation_jID:-0} -V -o $cytomod_RM_log -l mem_requested=$LOWER_MEM_REQ -b y $CYTOMOD_PATH -v -b -d masked_unmod_genome $cytomod_processing_dir$thresh/ --archiveOutDir $archiveDir/ --archiveOutName TRF-masked-archive --BEDOutDir $resultsDir -f $resultsDir/fullModGenome-TRFmasked.fa)"
            echo "$MODULE_LOAD_CMD samtools faidx $resultsDir/fullModGenome-TRFmasked.fa" | $BASE_QSUB_CMD -hold_jid ${track_creation_jID:-0},${cytomod_RM_jID:-0} -o $cytomod_RM_log -N 'faidx'
        fi
        # ---------------------------------- TFBS analysis ----------------------------------

        mkdir -p "$INITIAL_MEME_WORK_DIR"
        
        mkdir -p "$EXP_BED_DIR"

        shopt -s nullglob # turn on nullglob to allow the below loop to only use actual matches
        for peakFile in "$summit_BED_files_dir"/*.bed "$summit_BED_files_dir"/*.bed.gz; do
            peakFileBName=$(echo "${peakFile##*/}" | sed -r 's/\.(bed)(\.gz)?//')
            expBEDFile="$EXP_BED_DIR/$peakFileBName-$EXPANDED_NAME_STRING.bed"
            
            if [[ -f "$expBEDFile" ]]; then  # prevent repeatedly further expanding 500 bp regions when run from within the same directory
                continue
            fi

            # NB: this assumes the input file has summit positions only as the start coordinates, with end = start + 1
            # This is not true for ENCODE narrowPeaks, where it is in column 10
            # create all 500 bp expanded BED files for summits only, using BEDTools slop
            # one less on right, since intervals are half-open and we need exactly 500 bp regions
            # only uses first three columns (no name nor score/value)
            bedtools slop -i "$peakFile" -g "$GENOME_SIZE_FILE" -l 250 -r 249 | cut -f 1-3 > "$expBEDFile"
        done
        shopt -u nullglob

        INITIAL_WORKING_DIR="$(pwd)"

        for modGenome in $resultsDir/*-TRFmasked.fa; do
            # per sample, per threshold
            base_dir_name="$(perl -pe "s|.*$CYTOMOD_TRACK_DIR_PREFIX(\w+).*/(?:filesForArchive-)?(.+)/(.+)\..+$|\1-\2_\3|" <(echo $modGenome))"

            mkdir -p "$INITIAL_MEME_WORK_DIR/$base_dir_name" && cd "$_"

            parallel -j $NSLOTS --rpl '{fasta} s:.*/::; s:\.[^/.]+$:-mod.fa:' "
                # create FASTA file from expanded BED peak intervals, corresponding to the modified genome under consideration
                fastaFromBed_jID=0
                if [[ ! -f $PWD/{fasta} ]]; then
                    fastaFromBed_jID=\$(echo \"$MODULE_LOAD_CMD fastaFromBed -fi $INITIAL_WORKING_DIR/$modGenome -bed {} -fo $PWD/{fasta}\" | $BASE_QSUB_CMD -hold_jid $cytomod_RM_jID -N make_FASTA_from_BED-{/.} -o {//}/fastaFromBed-{/.}.err)
                fi

                # run MEME-ChIP for each modified FASTA file
                MEME_ChIP_jID=\$(echo \"meme-chip -xalph $ALPHABET_FILE -meme-minw $MIN_MEME_MOTIF_WIDTH -meme-maxw $MAX_MEME_MOTIF_WIDTH $par_MEME_arg -o {/.}-mod -db $MEME_MOTIF_DB {fasta}\" | tee -a $MEME_ChIP_CMDS_LOG | $BASE_QSUB_CMD -hold_jid \$fastaFromBed_jID -N MEME-ChIP_run-{/.} -o MEME-ChIP_{/.}-mod.log -R $cluster_res_param $MEME_ChIP_user_hold_param $par_env_arg)
                $BASE_QSUB_CMD -hold_jid \$MEME_ChIP_jID -o /dev/null -b y $EDIT_CENTRIMO_SCRIPT {/.}-mod/centrimo_out/centrimo.html

                RUN_H_TESTS_jID=\$(echo \"$RUN_H_TESTS_PATH {/.}-mod {/.}-mod $ALPHABET_FILE\" | $BASE_QSUB_CMD -hold_jid \$MEME_ChIP_jID -N hypothesis_testing -o MEME-ChIP_{/.}-mod-HTests.log -R $cluster_res_param)

                " ::: $(find "$INITIAL_WORKING_DIR/$EXP_BED_DIR/" -type f -name "*$EXPANDED_NAME_STRING.bed")

            cd "$INITIAL_WORKING_DIR"
        done

        if [[ ! -f metadata-curated.tsv ]]; then
            >&2 echo -e "${ERROR_FORMAT}The current directory must contain a 'metadata-curated.tsv' file.${RESET_FORMAT}"
            exit $ERR_EXIT
        fi

        mkdir results_dir_to_plot && cd $_

        parallel -j $NSLOTS --rpl "{name} s@(?:../)+$INITIAL_MEME_WORK_DIR/(?:FASTAs-for_MEME-ChIP_input-)?(?:(?<id>\w+)-)?(?<stringency>0.\d+)(?:_fullModGenome-TRFmasked)?[-_/](?<TF>.*?)(?:_summits)?-expandedTo\d+bpRegions-mod/.*@\$+{stringency}Stringency-\$+{id}-\$+{TF}@" "mkdir {name}; ln -s ../{} {name}/" ::: $(find $INITIAL_MEME_WORK_DIR/ -type d -name 'hypothesis_testing_selected_controlledVars_top_unmod_sig_mod_from_DREME')

        cd ..

        $PLOT_PATH |& tee processDataAndPlot.log
    done
done
