#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# Usage: ./run_hypothesis_tests.sh <MEME-ChIP input directory> [output directory] [alphabet file]

ERR_EXIT=64

MEME_INSTALL_DIR="$HOME/meme"

EDIT_CentriMo_PATH="$HOME/cytomod/utils/src/editCentriMoHTML.sh"

GEN_MOD_PERMUTATIONS_PATH="$HOME/cytomod/cytomod/src/genModPermutationsFromUnmodMotifs.py"

# assumed to be in the same directory as this script
PSEUDOCOUNT_SCRIPT="$(dirname $(readlink -f ${BASH_SOURCE[0]}))/computePseudocountForKLequality.py"

DEFAULT_PSEUDOCOUNT=0.1 # default to CentriMo's default value for its '--motif-pseudo' parameter

# was previously 4, now 5 due to the "consensus" column in CentriMo TXTs (version 4.11.4)
CENTRIMO_P_VAL_COL=5

# third argument or the alphabet with only 5mC and 5hmC modified bases
ALPHABET_FILE="${3:-$MEME_INSTALL_DIR/src/unit-tests/data/mod_spec_no5caC_no5fC.txt}"

BACKGROUND_FILE_BASENAME='background'

FIXED_CENTRAL_REGION_WIDTH=100

# trys all possible combinations of modifying Cs and Gs
# 5mC only...
POS_UNMOD_BASE='C'
NEG_UNMOD_BASE='G'
POS_STRAND_MODS_TO_TEST=('m')
NEG_STRAND_MODS_TO_TEST=('1')

MODS_TO_TEST="$(printf '%s' ${POS_STRAND_MODS_TO_TEST[@]})"

NUM_TOP_UNMOD_MOTIFS=3
THRESH_E='1e-3'  # stringent Bayesian-motivated threshold, recommended by Johnson (PNAS, 110:19313-7, 2013)
MIN_MOTIF_LEN=5

# get the hypothesis set from the top x unmodified DREME motifs, plus the unmodified versions of all significant CentriMo modified motifs
# only get those of sufficient length
CentriMo_sig_DREME_motifs=$(grep -v '^#' $1/centrimo_out/centrimo.txt | awk "{if (\$$CENTRIMO_P_VAL_COL < $THRESH_E && length(\$2) >= $MIN_MOTIF_LEN) {print \$2;}}" | grep -v '^MA' | grep -v '^[1-9]$' || true)

if [[ -z "$CentriMo_sig_DREME_motifs" ]]; then
    >&2 echo -e "No significant DREME motifs found within CentriMo results.\n"
    exit 0
fi

mod_hypotheses=$(echo "$CentriMo_sig_DREME_motifs" | grep '[a-z1-9]' || true) # may not be any
unmod_top_x=$(echo "$CentriMo_sig_DREME_motifs" | sort -g | head -$NUM_TOP_UNMOD_MOTIFS)
unmod_hypotheses=($(echo -e "$mod_hypotheses\n$unmod_top_x" | sort | uniq | tr 'a-z' 'C' | tr '1-9' 'G' | grep -P ".{$MIN_MOTIF_LEN,}"))


# getAllPureCpGModifiedIUPACMotifs <array of unmodified motifs [PASSED BY NAME]> <string of modifications to use (primary / positive strand only)>
# This uses Cytomod's genModPermutationsFromUnmodMotifs.py script to generate all motifs.
function getAllPureCpGModifiedIUPACMotifs {
    declare -a unmod_motifs=("${!1-}")
    
    local modifications="${2}"
    
    result=($($GEN_MOD_PERMUTATIONS_PATH -m "$modifications" "${unmod_motifs[@]}"))
}


# getBGFreqForBase <background file> <base whose background frequency we wish to find>
function getBGFreqForBase {
    local background_file="${1}"
    local base="${2}"
    local base_bg=$(grep -oP "^${base}\s+\K.+" "$background_file")
    result="$base_bg"
}


# computePseudocount <unmodified base frequency for a base b1> <modified base frequency for the modification of b1> [unmodified base frequency for a base b2] [modified base frequency for the modification of b2]
# NB: input frequencies are real numbers and may be in scientific notation,
#    of the "e" variety.
# Assumes that bases under consideration both have profile frequencies of 1.
# Computes the pseudo-count parameter s.t. the KL divergence, caused by
# the difference in background frequencies for the particular nucleobase
# under consideration, is equalized.
# This pseudocount value, when set for one motif, equalizes the KL divergences
# of both motifs, accounting for the pseudocount applied to the other motif.
function computePseudocount {
    num_mod_bases=$1
    unmod_bg_1=$2
    mod_bg_1=$3
    unmod_bg_2=${4:-}
    mod_bg_2=${5:-}
    
    pseudocount="$(python $PSEUDOCOUNT_SCRIPT $num_mod_bases $unmod_bg_1 $mod_bg_1 $unmod_bg_2 $mod_bg_2 $DEFAULT_PSEUDOCOUNT)"

    result="$pseudocount"
}


# the MEME-ChIP directory is the first argument
MEMEChIP_dir=$(readlink -f ${1-})  # use absolute path

background_file="$MEMEChIP_dir/$BACKGROUND_FILE_BASENAME"

# optional second argument, default to the MEME-ChIP directory (first required argument)
output_root_dir=$(readlink -f ${2-$MEMEChIP_dir}) # use absolute path

H_test_dir="$output_root_dir/hypothesis_testing_selected_controlledVars_top_unmod_sig_mod_from_DREME"

aggregated_motifs_file='hypotheses.meme'

if [[ -d "$H_test_dir" ]]; then
    >&2 echo -e "Found existing hypothesis testing directory. Removed; will be overwritten.\n"
    rm -Rf "$(readlink -f "$H_test_dir/")"
fi

mkdir "$H_test_dir" && cd $_

# create hypotheses
echo -e "The following hypotheses will be tested:\n"
# create seperately even though iupac2meme can accept all in one command,
# since we need to vary the pseudo-count parameter and
# unify the hypotheses into a single file as they are created.
for unmod_hypothesis in "${unmod_hypotheses[@]}"; do
    getAllPureCpGModifiedIUPACMotifs unmod_hypothesis "$MODS_TO_TEST"
    
    if [[ ${#result[@]} -le 0 ]]; then
        >&2 echo -e "No input motifs found.\n"
        exit $ERR_EXIT
    fi

    hypotheses=("${result[@]}")
    motif_pseudo=0
    
    for hypothesis in "${hypotheses[@]}"; do
        if [[ "$hypothesis" =~ ["${POS_STRAND_MODS_TO_TEST[@]}"] || "$hypothesis" =~ ["${NEG_STRAND_MODS_TO_TEST[@]}"] ]]; then # modified (one or more modified bases; either hemi-modified or completely)
            i=0
            for mod_base in "${POS_STRAND_MODS_TO_TEST[@]}" "${NEG_STRAND_MODS_TO_TEST[@]}"; do
                if [[ "$hypothesis" =~ $mod_base ]]; then
                    getBGFreqForBase "$background_file" "$mod_base"
                    mod_bg=${result}
                    getBGFreqForBase "$background_file" "$([[ ${POS_STRAND_MODS_TO_TEST[@]} =~ $mod_base ]] && echo $POS_UNMOD_BASE || echo $NEG_UNMOD_BASE)"
                    unmod_bg=${result}

                if [[ -z "${mod_bg_1:-}" ]]; then
                    mod_bg_1="$mod_bg"
                    unmod_bg_1="$unmod_bg"
                else
                    mod_bg_2="$mod_bg"
                    unmod_bg_2="$unmod_bg"
                    break
                fi
            fi
                ((++i))
            done
            computePseudocount "$(echo $hypothesis | tr -dc $(printf %s ${POS_STRAND_MODS_TO_TEST[@]}${NEG_STRAND_MODS_TO_TEST[@]}) | wc -c)" "$unmod_bg_1" "$mod_bg_1" "${unmod_bg_2:-}" "${mod_bg_2:-}"
            motif_pseudo=${result}
        else # unmodified, use 
            motif_pseudo=$DEFAULT_PSEUDOCOUNT
        fi
        echo -e "$hypothesis with pseudo-count $motif_pseudo\n"

        existing_motif_file=$([[ -f $aggregated_motifs_file ]] && echo "$aggregated_motifs_file" || echo '')
        meme2meme -consensus $existing_motif_file <(iupac2meme -alph "$ALPHABET_FILE" -bg "$background_file" -pseudo $motif_pseudo "$hypothesis") >> "$aggregated_motifs_file-1" && mv -f "$aggregated_motifs_file-1" "$aggregated_motifs_file"

        unset mod_bg_1 unmod_bg_1 mod_bg_2 unmod_bg_2
    done
done

# run CentriMo
# uses fixed central region widths (minreg set to one less than max, since required to be less than max)
# scans for all scores >= 0 (not default of 5)
# optimizes scores (unlike default, which uses the min)
# set motif pseudocount to 0, since already applied earlier (in iupac2meme)
cd "$MEMEChIP_dir/.."

FASTA_in=''
shopt -s nullglob extglob

# all FASTAs gzipped or not and also any file descriptors, whose names are only numbers (e.g. "63")
# match one or more and then check right after if more than one way matched
set -- "$MEMEChIP_dir"/+(*.fa?(.gz)|[1-9][0-9])

if [[ "$#" -eq 1 ]]; then
    FASTA_in="$(readlink -f $@)"
else
    >&2 echo "Unable to find the FASTA input file for CentriMo or more than one such file found."
    exit $ERR_EXIT
fi

shopt -u nullglob extglob

centrimo -o "$H_test_dir/centrimo_out" -xalph "$ALPHABET_FILE" -seqlen 500 --bfile $background_file --minreg $(($FIXED_CENTRAL_REGION_WIDTH - 1)) --maxreg $FIXED_CENTRAL_REGION_WIDTH --score 0 --optimize_score --motif-pseudo 0 <(zcat -f "$FASTA_in") "$H_test_dir/$aggregated_motifs_file"
$EDIT_CentriMo_PATH "$H_test_dir/centrimo_out/centrimo.html"

