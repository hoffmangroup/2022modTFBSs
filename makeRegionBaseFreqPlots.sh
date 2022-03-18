#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# NB: assembly-specific promoter definitions and tissue-specific enhancer definitions (merge of four histone modifications in a mouse/human dataset from doi:10.1038/nature13992, where applicable.

# Usage: ./makeRegionBaseFreqPlots.sh <organism> <assembly> <input_FASTA_file> <optional -o argument for the output prefix>

ERR_EXIT=64

BASEFREQ_BASE_PATH="$HOME/cytomod/cytomod/src/baseFrequency.py"

BASEFREQ_BASE_CMD="$BASEFREQ_BASE_PATH -Ct -r \"${region_desc:-}\""

GENCODE_FNAME='gencode.v26lift37.annotation.gtf.gz'

ORGANISM=${1:-}

ASSEMBLY=${2:-}

FASTA_IN=${3:-}

shift 3

if [[ "$@" =~ "-o " ]]; then # attempt to parse output directory
    output_dir=$(echo "${@:-}" | grep -oP -- '-o\s+\K[^\s+]+/')
fi
if [[ ! -d ${output_dir-} ]]; then
    output_dir=$(pwd)
fi


case $ORGANISM in
    'mouse')
        case $ASSEMBLY in 
        'GRCm37'|'NCBI37'|'mm9')
            GENCODE_ASSEMBLY_ID='M1'
            SEG_FILE="$HOME/cytomod/experiments/linked-2015-04-08-mESC_modGenome_Cytomod_baseFreqs_ChromHMM_Enhancers/ChromHMM_7states_4HistoneMarks_Segmentations7_Nature2014_ES-Bruce4.gz"
            ;;
        'GRCm38'|'mm10')
            GENCODE_ASSEMBLY_ID='M13'
            ;;
         ?)
            >&2 echo "Assembly $ASSEMBLY is not supported."
            exit $ERR_EXIT
        esac
        ;;
    'human')
        case $ASSEMBLY in
        'GRCh37'|'hg19')
            GENCODE_ASSEMBLY_ID='19'
            ;;
        'GRCh38'|'hg38')
            GENCODE_ASSEMBLY_ID='26'
            ;;
         ?) 
            >&2 echo "Assembly $ASSEMBLY is not supported."
            exit $ERR_EXIT
        esac
        ;;
    ?)
        >&2 echo "Only mouse or human are supported."
        exit $ERR_EXIT
        ;;
esac

GENCODE_ADDTL_PATH="release_${GENCODE_ASSEMBLY_ID}/gencode.v${GENCODE_ASSEMBLY_ID}.annotation.gtf.gz"
GENCODE_FNAME="${GENCODE_ADDTL_PATH##*/}"

# genome-wide (organism, assembly, and tissue independent)
output_prefix="$output_dir/genome-wide"
region_desc='genome-wide'
eval "$BASEFREQ_BASE_CMD -f $FASTA_IN" -o "$output_prefix" > "$output_prefix.txt"


# promoters (assembly dependent; tissue independent)
output_prefix="$output_dir/promoters"
region_desc='2 kb upstream to GENCODE known gene TSSs'
if [[ ! -e "${GENCODE_FNAME}" ]]; then
    if [[ ! -e "$HOME/${GENCODE_FNAME}" ]]; then
        wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_${ORGANISM}/${GENCODE_ADDTL_PATH}
    else
        GENCODE_FNAME="$HOME/${GENCODE_ADDTL_PATH##*/}"
    fi
fi
zcat "$GENCODE_FNAME" | gtf2bed | awk 'BEGIN{FS=OFS="\t"} {if ($8=="transcript" && /transcript_type "protein_coding"/ && /transcript_status "KNOWN";/) { split($10, info, " "); if (!seen[info[2]]) {print $1,$2-2000,$2,$4"_"gensub(/[";]|(\.[[:digit:]])/, "", "g", info[2]),$5,$6} seen[info[2]]++;}}' | bedtools getfasta -fi $FASTA_IN -bed /dev/stdin -fo /dev/stdout | eval "$BASEFREQ_BASE_CMD -f - -o $output_prefix" > "$output_prefix.txt"


# enhancers (assembly and tissue dependent)
if [[ -f ${SEG_FILE:-} ]]; then  # only if we have a defined enhancer annotation
    output_prefix="$output_dir/enhancers"
    region_desc='ChromHMM_7 Hardison Nature 2014 Enhancers'
    zcat "$SEG_FILE" | awk 'BEGIN{FS=OFS="\t";} $5==7 {print $2,$3,$4;}' | bedtools getfasta -fi $FASTA_IN -bed /dev/stdin -fo /dev/stdout | eval "$BASEFREQ_BASE_CMD -f - -o $output_prefix" > "$output_prefix.txt"
fi

