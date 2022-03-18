#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# NB: all output is saved within the CWD. Change to the desired output directory first.

ERR_EXIT=64

KEEP_TEMP_FILE=true

REMOVE_ZIP_ARCHIVES=false # critical to use with current regression, of recursive ZIP inclusion

NSLOTS=${NSLOTS:-1}

THRESH_MOD_UNMOD_PREF=0

# ignore if within this amount of $THRESH_MOD_UNMOD_PREF (i.e. this close to 0)
EPSILON_THRESH_NEG=30 # for unmod pref.
EPSILON_THRESH_POS=2 # for mod pref.

NAME_OF_HYPOTHESIS_MEME_FILES='hypotheses.meme'

if [[ "$#" -lt 4 ]]; then
    >&2 echo "Usage: ./run_matrix-clustering-specific_TFs.sh <ENCODE metadata file> <name of TF> <hypothesis testing result file (CSV from Pandas dataframe, generated from plot script)> <directory containing \"hypotheses.meme\" files in subdirectories, to be read in iff they pertain to the TF being processed>"
    exit $ERR_EXIT
fi

function maybe_remove_ZIP {
    ZIP_to_rm="${1:-}"
    
    if [[ $# -gt 1 ]]; then
        >&2 echo "WARNING: multiple ZIP filenames were specified for removal. Only the first one ($ZIP_to_rm) will be removed."
    fi

    if [[ $REMOVE_ZIP_ARCHIVES == true && -f "$ZIP_to_rm" ]]; then
        rm -f "$ZIP_to_rm"
    fi  
}

ENCODE_METADATA_FILE="${1:-}"

OUTPUT_PREFIX='clustering_out'

LOWER_NCOR_THRESHOLD=0.4
LOWER_COR_THRESHOLD=0.6

NCOR_LOWER_PARAM="-lth Ncor $LOWER_NCOR_THRESHOLD"
COR_LOWER_PARAM="-lth cor $LOWER_COR_THRESHOLD"

RSAT_MATRIX_CLUSTERING_BASE_CMD="matrix-clustering -v 5 -matrix_format meme -return align_consensus,heatmap,json,newick,nb_clusters $NCOR_LOWER_PARAM $COR_LOWER_PARAM"

MATRIX_CLUSTERING_PREFIX='matrix_clustering'

MATRIX_CLUSTERING_EXT='.log'

TF_TO_CLUSTER="${2^^}"

HYPOTHESIS_TESTING_RESULTS_FILE="${3:-}"

DIR_CONTAINING_MEME_HYPOTHESES_MOTIF_SETS="${4:-}"

TMP_FILE_DIR="${TMPDIR:-/tmp}"

UNIQUE_ID="${TF_TO_CLUSTER}-$RANDOM"

# verify the ID is unique for current temporary files
# could iteratively try again with new random suffixes,
# but instead just terminate, since this is unlikely
# and the directory should be cleaned of these often
if [[ -f $TMP_FILE_DIR/*$UNIQUE_ID* ]]; then
    >&2 echo "ID generated ($UNIQUE_ID) was not unique. Clean files within $TMP_FILE_DIR and try again."
    exit $ERR_EXIT
fi

MEME_MOTIF_TEMP_FILE_FOR_MOTIF_EXTRACTION="$TMP_FILE_DIR/MEME_motifs_for_motif_extraction-$UNIQUE_ID.meme"

MEME_MOTIF_TEMP_FILE_FOR_RSAT_PATH_PREFIX="$TMP_FILE_DIR/MEME_motifs_for_RSAT-$UNIQUE_ID"

if [[ $KEEP_TEMP_FILE == false ]]; then
    trap "rm -f $MEME_MOTIF_TEMP_FILE_FOR_MOTIF_EXTRACTION ${MEME_MOTIF_TEMP_FILE_FOR_RSAT_PATH_PREFIX}*.meme" EXIT
fi

OTHER_POSSIBLE_ID_ALTERNATION_STRING="$TF_TO_CLUSTER"
ENCODE_TF_MATCHES="$(fgrep -iw $TF_TO_CLUSTER $ENCODE_METADATA_FILE | cut -f1 | perl -0777 -pe 's/\n/|/g;s/\|$//' || true)"

if [[ -n "$ENCODE_TF_MATCHES" ]]; then
    OTHER_POSSIBLE_ID_ALTERNATION_STRING="${OTHER_POSSIBLE_ID_ALTERNATION_STRING}|$ENCODE_TF_MATCHES"
fi

# NB: these Strogantsev et al. ZFP57 ChIP-seq datasets do not contain "ZFP57" within their paths, but only the mouse hybrid identifier (either BC8 or CB9)
if [[ "$TF_TO_CLUSTER" == 'ZFP57' ]]; then
    OTHER_POSSIBLE_ID_ALTERNATION_STRING="${OTHER_POSSIBLE_ID_ALTERNATION_STRING}|BC8|CB9"
fi

MEME_HYPOTHESES_MOTIF_SETS_FOR_TF="$(find -L $DIR_CONTAINING_MEME_HYPOTHESES_MOTIF_SETS/ -regextype posix-extended -iregex "$DIR_CONTAINING_MEME_HYPOTHESES_MOTIF_SETS/.*(${OTHER_POSSIBLE_ID_ALTERNATION_STRING}).*" -and -name $NAME_OF_HYPOTHESIS_MEME_FILES | tr '\n' ' ')"

# NB: only strict inequalities are used; motifs with ln(p_mod) - ln(p_unmod) == $THRESH_MOD_UNMOD_PREF + $EPSILON_THRESH are discarded

# always includes both modified and unmodified motifs

HYPOTHESIS_DATA_TO_TEST_CONSENSUS_MODIFIED_MOTIFS="$(fgrep -iw $TF_TO_CLUSTER $HYPOTHESIS_TESTING_RESULTS_FILE)"

# test all motif pairs that are not within epsilon ($EPSILON_THRESH) of 0 ($THRESH_MOD_UNMOD_PREF). Distance permitted differs for negative vs. positive values, to reflect different empirical signal strengths.
HYPOTHESIS_PAIRS_TO_TEST_CONSENSUS_MODIFIED_MOTIFS_SUBSTANTIAL_EITHER_PREF=$(echo "$HYPOTHESIS_DATA_TO_TEST_CONSENSUS_MODIFIED_MOTIFS" | awk "(\$NF < ($THRESH_MOD_UNMOD_PREF - $EPSILON_THRESH_NEG)) || (\$NF > ($THRESH_MOD_UNMOD_PREF + $EPSILON_THRESH_POS)) {print \$1}")

# combine into a single MEME motif file, using consensus IDs, and motif consensus only (non-unique) as alt. ID, to ensure meme-get-motif can look it up
meme2meme -xalph -consensus $MEME_HYPOTHESES_MOTIF_SETS_FOR_TF | sed -r 's/(MOTIF )([[:alnum:]]+)(\.[[:digit:]]+)?/\1\2\3 \2/' > "$MEME_MOTIF_TEMP_FILE_FOR_MOTIF_EXTRACTION"

RSAT_input_motif_filename="${MEME_MOTIF_TEMP_FILE_FOR_RSAT_PATH_PREFIX}.EITHER-no_ambig.meme"

parallel -j $NSLOTS -X --files --tmpdir "$TMP_FILE_DIR" "meme-get-motif -ia {} $MEME_MOTIF_TEMP_FILE_FOR_MOTIF_EXTRACTION" ::: '-id' ::: $(echo "$HYPOTHESIS_PAIRS_TO_TEST_CONSENSUS_MODIFIED_MOTIFS_SUBSTANTIAL_EITHER_PREF" | tr "-" "\n" | sort -u) | \
    parallel -X -j 1 "meme2meme -xalph" > $RSAT_input_motif_filename

matrix_clustering_log_filename="${MATRIX_CLUSTERING_PREFIX}$TF_TO_CLUSTER-EITHER-no_ambig_pref$MATRIX_CLUSTERING_EXT"

# cluster each alone
output_prefix="$OUTPUT_PREFIX-$TF_TO_CLUSTER-EITHER-no_ambig_pref"
$RSAT_MATRIX_CLUSTERING_BASE_CMD -o "$output_prefix" -matrix "$(basename ${RSAT_input_motif_filename%.*})" "$RSAT_input_motif_filename" 2> $matrix_clustering_log_filename

maybe_remove_ZIP $output_prefix*.zip

rm -Rf js/ && ln -s $RSAT/perl-scripts/lib/js 2>/dev/null # remove copy of js library files and use a symlink
