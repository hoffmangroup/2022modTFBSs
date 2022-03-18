#!/usr/bin/env python
# -*- coding: utf-8 -*-

# NB: This script only supports 5mC only or 5mC and 5hmC data.
# NB: Plotting will fail with too few TFs.

from __future__ import (with_statement, division, print_function,
                        unicode_literals)

import argparse
import codecs
import glob2
import inflect
import re

_CURRENTLY_SUPPORTED_MOD_BASES = ['m', 'h']

_DEFAULT_CSV_PATH_FILENAME_BASE = './df_for_plot-top'
_DEFAULT_CSV_PATH_FILENAME_EXT = '.csv'

_DEFAULT_ROOT_INPUT_DIR_NAME = 'result_dir_to_plot'

parser = argparse.ArgumentParser()

process_or_plot_only = parser.add_mutually_exclusive_group()

process_or_plot_only.add_argument('-p', '--no_plot', action='store_true',
                                  help="Only perform processing, to produce "
                                  "output CSV files, but do not perform "
                                  "subsequent plotting.")
process_or_plot_only.add_argument('-P', '--plot_only_from_CSVs', nargs='?',
                                  type=str,
                                  const=_DEFAULT_CSV_PATH_FILENAME_BASE,
                                  help="Use provided CSVs (from a prior run) "
                                  "for plotting, without generating any newly "
                                  "processed data. Optionally, provide a full "
                                  "path to the CSV files name, which must "
                                  "consist of the path and file name, up to, "
                                  "but exclusive of a label denoting the "
                                  "number of top hypotheses used (e.g. \"3\")"
                                  ", followed by \".csv.\", which cannot "
                                  "be altered. If this is not provided, "
                                  "it will default to: \"{}\"".
                                  format(_DEFAULT_CSV_PATH_FILENAME_BASE))

parser.add_argument('-D', '--root_input_dirname',
                    default=_DEFAULT_ROOT_INPUT_DIR_NAME,
                    help="The default root directory, containing all CentriMo "
                    "hypothesis testing, according to the prescribed "
                    "structure. This can be a symbolic link to a pre-existing "
                    "directory and defaults to \"{}\"".
                    format(_DEFAULT_ROOT_INPUT_DIR_NAME))
parser.add_argument('-K', '--keep_mixed_mod_motifs', action='store_true',
                    help="Keep modified motifs that contain more than a "
                         "single type of modified base. By default, these "
                         "motifs are removed. Examples of such motifs "
                         " include 'CAm2TG' or 'hAm1TG'."
                         "Note that this program is not presently "
                         "designed to account for such motifs when plotting. "
                         "Such motifs will be processed correctly, and "
                         "marked as mixed, via their ambiguity code being "
                         "used in their Modification column (e.g. 'x'), "
                         "but the plot will have misaligned columns and will "
                         "be poorly formatted.")
parser.add_argument('-v', '--verbose', action="count",
                    help="Increase output verbosity. Will result in the "
                    "of the actual significance decile cutoffs used, "
                    "iff this option is invoked with a value of at least 1.")

args = parser.parse_args()

if not args.no_plot:
    import matplotlib as mpl

    mpl.use('pgf')

    # switch to DejaVu fonts: sans-serif font needed for Unicode in titles,
    # whereas the others are used for consistency
    mpl.rc('font',
           **{'family': 'serif', 'serif': ['DejaVu Serif'],
              'sans-serif': ['DejaVu Sans'],
              'monospace': ['DejaVu Sans Mono']})
    # use the STIX sans-serif font for mathematics
    mpl.rc('mathtext', fontset='stixsans')

    mpl.rc('pgf', texsystem='lualatex')  # use LuaLaTeX for non-strict mem.

    from matplotlib.backends.backend_pgf import FigureCanvasPgf
    mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)

    import matplotlib.pyplot as plt

    import seaborn as sns

from scipy import stats

import numpy as np
import pandas as pd

from functools32 import lru_cache
from io import StringIO

import cytoUtils as cUtils
from cytoUtils import warn, die, assert_or_die_with_msg

from getTFFamilyInfo import getTFFamilyInfo

inf_eng = inflect.engine()

_ACTUAL_SIG_QUANTILE_CUTS_FILENAME = 'df_alldata-SigQuantiles-raw.csv'
_EXCLUDED_MIXED_MOD_MOTIFS_FILENAME_BASE = ('excluded_mixed_'
                                            'mod_motifs_df-top')

BASE_PLOT_NAME = 'plot'

_EXTENSIONS_TO_OUTPUT = ['pgf', 'pdf', 'png']

# 18 plus computed log diff and percentile columns
_NUM_ALLDATA_COLS = 20

NUM_QUANTILES = 10  # use deciles

MAX_P_VALUE = 0.05
MAX_LOG_P_VALUE = np.log(MAX_P_VALUE)

_UNMOD_COLS = slice(0, 5, None)
_MOD_COLS = slice(4, None, None)

METADATA_FILE_FULLPATH = 'metadata-curated.tsv'

H_TEST_DIR_NAME = ('hypothesis_testing_selected_controlledVars_top'
                   '_unmod_sig_mod_from_DREME')

CENTRIMO_DIR = 'centrimo_out'

CENTRIMO_FILES_NAME = 'centrimo.txt'

# using human TFClass, despite mouse data, since more comprehensive
# and all family and superfamily information applies equally to mouse
ONTOLOGY = 'TFClass/TFClass_human-curated.obo'  # from CWD

FUSED_EPITOPE_REGEX = re.compile(r'(?:[-_]|^)(FLAG|eGFP)(?:[-_]|$)',
                                 flags=re.I)

# directory or file labels, which uniquely identify a sample as ZFP57
ZFP57_ID_strings = ['MACS2_BC8', 'MACS2_CB9']

# names of Seaborn plotting methods
_SEABORN_SWARM_PLOT_TYPE = 'swarmplot'
_SEABORN_VIOLIN_PLOT_TYPE = 'violinplot'

# swarmplot scales poorly, so switch to a better-scaling plot if above this
_MAX_VALS_FOR_SWARMPLOT = 20000

PLOT_WIDTH = 36

# 1: first, 2: second, etc.
DECILE_ORD_LABELS = [inf_eng.number_to_words(inf_eng.ordinal(i + 1))
                     for i in range(10)]

PLOT_X_AXIS_LABEL = (r'$\mathrm{ln}\left(p_{unmod}\right) - '
                     r'\mathrm{ln}\left(p_{mod}\right)$')

PLOT_Y_AXIS_LABEL = 'Transcription Factor'

PRESET_SNS_CONTEXT = 'talk'

PLOT_DEFAULT_PALETTE = 'dark'

PLOT_X_TICK_SIZE = 8

PLOT_FONT_SCALE_FACTOR = 2.1

# sugg. default: 0.125 (higher => more space)
PLOT_LEFT_PARAM = 0.425

# sugg. default: 0.9 (higher => more space)
PLOT_DEFAULT_TOP_PARAM = 0.96

# sugg. default: 0.9 (lower => more space)
PLOT_RIGHT_PARAM = 0.86

PLOT_LEFT_SPINE_DROP = 40  # in points

PLOT_LEGEND_OUT_ADDTL_DIST = 0.05

PLOT_DEFAULT_LEGEND_OUT_BASE_DIST = 0.1

# amount to multiply the label offset for spans
# of additional y-axis spines
PLOT_A_MULT_FACTOR = 12

# offset of additional y-axis spines
PLOT_X_OFFSET = 800

# width of boxes containing the numerical right axis annotations
PLOT_DEFAULT_R_AN_WIDTH = 3

pd.set_option('max_colwidth', -1)  # set unlimited col width


def _iterrowsOrIdentity(df_or_object):
    """Generalization of iterrows to permit polymorphic application.

       Returns the result of calling iterrows iff the input
       is a DataFrame, the object itself, augmented with
       its indices if it was already iterable and not a Series,
       otherwise a list of two singleton lists, consisting
       respectively of '0' and the input object."""

    if isinstance(df_or_object, pd.DataFrame):
        return df_or_object.iterrows()
    elif isinstance(df_or_object, pd.Series):
        # keep entire Series together, i.e., use as a row
        return enumerate([df_or_object])
    else:
        return enumerate(df_or_object)


def _compute_plot_height(num_TFs):
    """Computes the height of the plot.
       Set based on original plot height of ~32 for ~52 TFs.
       Uses an arbitrary min below which a constant is used."""

    CONST_FACTOR = 0.62

    return CONST_FACTOR * num_TFs if num_TFs > 4 else 4


@lru_cache(maxsize=None)
def getCuratedTFClassFamilyInfo(query_TF):
    """Lookup (concise) TF family and superfamily information.
       Use the curated TFClass ontology and cache the results
       of largely redundant queries."""

    # override the query TF manually in this case
    for ZFP57_ID_str in ZFP57_ID_strings:
        if ZFP57_ID_str in query_TF:
            query_TF = 'ZFP57'
            break

    # remove fusions that cannot impact the queried TF family
    query_TF = FUSED_EPITOPE_REGEX.sub('', query_TF)

    result = getTFFamilyInfo(ONTOLOGY, query_TF)

    if result is None:  # always return a pair
        result = (None, None)
    else:
        # perform substitutions for concision
        result = tuple(re.sub(r'.+\((.+)\)(.*)', r'\1\2',
                              re.sub(r' DNA-binding domains?', r'',
                                     # need to escape, since interp. later
                                     re.sub(r'[Aa]lpha', r'$\\alpha$',
                                            re.sub(r'[Bb]eta', r'$\\beta$',
                                                   re.sub(r'zinc finger( '
                                                          r'factors)?', r'ZFs',
                                                          t)))))
                       for t in result)

    return result


def getModType(motif):
    """Returns a string describing the motif's modifications
       using the minimal number of primary base codes."""

    found_uniq_mod_bases = \
        list(set([base for base in cUtils.ALL_NON_UNMOD_BASES
                  if base in motif or
                  cUtils.complement(base) in motif]))

    assert_or_die_with_msg(found_uniq_mod_bases, 'Unable to find'
                           'cognate modified base(s).')

    # concat. of all unique bases, if nothing else possible
    result = ''.join(found_uniq_mod_bases)

    # a unique single mod, if possible
    if len(found_uniq_mod_bases) == 1:
        result = found_uniq_mod_bases[0]

    # an ambiguity code if one exists for the set
    for key, val in cUtils.AMBIG_MOD_BASES.iteritems():
        if set(found_uniq_mod_bases) == set(val):
            result = key

    return result


def _warn_if_new_missing_unmod_motif(_seen_missing_unmod_motifs_for_warn_check,
                                     mod_motif, mod_ID):
    unmod_motif = cUtils.translUnmodSeq(mod_motif)

    if unmod_motif not in _seen_missing_unmod_motifs_for_warn_check:
        warn('No matching unmodified motif (i.e. {}) was found for '
             'the (first) modified motif: {} in {}. This is likely '
             'fine, generally meaning that the unmodified motif '
             'was of insufficient significance to make it '
             'into the CentriMo output at all. All modifications '
             'of this unmodified motif are accordingly skipped, '
             'for this TF. This should be reviewed and may '
             'require curation.'.
             format(unmod_motif, mod_motif, mod_ID))

        _seen_missing_unmod_motifs_for_warn_check[unmod_motif] = None


def _get_LaTeX_subscript(text):
    """Returns LaTeX code that specifies that the provided
       text, expected to be used as a suffix, will be
       rendered in subscript font."""

    return "\\textsubscript{{{}}}".format(text)


def edit_TF_labels_for_display(str_TF_labels):
    """Returns the input list of TF labels, edited to
       be more suitable for display."""

    for label_idx, TF_label in enumerate(str_TF_labels):
        TF_label_disp = TF_label

        fused_epitope_match = FUSED_EPITOPE_REGEX.search(TF_label)

        if fused_epitope_match:
            # move the fused epitope label to the end of the string,
            # as a subscript of the main TF name
            TF_label_disp = ("{}{}{}".
                             format(TF_label[:fused_epitope_match.start()],
                                    TF_label[fused_epitope_match.end():],
                                    _get_LaTeX_subscript(fused_epitope_match.
                                    string[fused_epitope_match.start(1):
                                           fused_epitope_match.end(1)])))

        # format RNA Polymerase II subunit RPB1, phosphorylation AA pos. codes

        def _cons_phospho_repl(match):
            return '{}{}'.format(match.group(1),
                                 _get_LaTeX_subscript(match.group(2)))

        TF_label_disp = re.sub(r'(POLR2(?:[A-Z]|RPB\d+))'
                               r'(?:phospho|-)([A-Z]\d)',
                               _cons_phospho_repl, TF_label_disp, flags=re.I)

        # overwrite the previous label with the one to display
        str_TF_labels[label_idx] = TF_label_disp

    return str_TF_labels


# ----------------------------------------------------
# Adapted from: http://stackoverflow.com/a/3919443


def annotate_group(ax, name, yspan, labeloffset, colour):
    """Annotates a span of the y-axis"""
    _PLOT_SPINE_A_ADJ = 0.7

    def annotate(ax, name, left, right, x, pad, colour):
        arrow = \
            ax.annotate(name,
                        xy=(x, left), xycoords='data',
                        xytext=(x-pad, right), textcoords='data',
                        annotation_clip=False, verticalalignment='top',
                        horizontalalignment='center', linespacing=2.0,
                        arrowprops=dict(arrowstyle='-', shrinkA=0, shrinkB=0,
                                        connectionstyle='angle,angleB=180,'
                                                        'angleA=90,rad=5',
                                        color=colour),
                        fontsize='x-small')
        return arrow

    xmin = ax.get_xlim()[0] - PLOT_A_MULT_FACTOR * labeloffset
    xpad = 0.01 * np.ptp(ax.get_xlim())
    ycenter = np.mean(yspan) - 1  # minus one to align correctly

    left_arrow = annotate(ax, name, yspan[0] - _PLOT_SPINE_A_ADJ,
                          ycenter, xmin, xpad, colour)
    right_arrow = annotate(ax, name, yspan[1] - _PLOT_SPINE_A_ADJ,
                           ycenter, xmin, xpad, colour)
    return left_arrow, right_arrow


def _make_addtl_spine(ax, label, labeloffset, offset=0):
    """Makes a second bottom spine"""

    second_bottom = mpl.spines.Spine(ax, 'left', ax.spines['left']._path)
    second_bottom.set_position(('outward', offset))
    ax.spines['second_bottom'] = second_bottom

    # Make a new ylabel
    ax.annotate(label, xy=(0, 0.5), xycoords='axes fraction',
                xytext=(-labeloffset, 0), textcoords='offset points',
                verticalalignment='center', horizontalalignment='left',
                rotation=90)


# ----------------------------------------------------


def add_spine_annot(ax, name, labels):
    """Add an annotation as an additional spine of the y-axis."""

    global num_addtl_spines

    num_addtl_spines += 1

    # compute the offset for the axis label, accounting for initial y label
    approx_offset_from_y = 190
    label_adj_f = 130  # additional adjustment factor (moves label inward)
    offset_to_mult = ((0.9 * num_addtl_spines) *
                      (PLOT_X_OFFSET - approx_offset_from_y) -
                      (num_addtl_spines * label_adj_f))

    _make_addtl_spine(ax, name, approx_offset_from_y + offset_to_mult)

    pal_iter = iter(sns.color_palette("Set1", n_colors=len(labels), desat=.5))

    for name, yspan in labels:
        annotate_group(ax, name, yspan, num_addtl_spines * PLOT_X_OFFSET,
                       pal_iter.next())


def retrieve_input_data(metadata_fullpath, root_input_dirname):
    metadata = codecs.open(metadata_fullpath, encoding='UTF-8').read()

    fileset = glob2.iglob(root_input_dirname + '/**/' +
                          '/'.join((H_TEST_DIR_NAME, CENTRIMO_DIR,
                                    CENTRIMO_FILES_NAME)))

    processed_input_file = ''

    # NB: this currently does not record sex information,
    #     but just combines all per TF.
    # Replicate information is recorded as part of the ID
    # and used to correctly pair unmodified-modified hypotheses.
    for num, file_name in enumerate(fileset):
        file_name_for_label = file_name
        # manually label known ZFP57 IDs as ZFP57
        for ZFP57_id_str in ZFP57_ID_strings:
            file_name_for_label = (file_name_for_label.
                                   replace(ZFP57_id_str, "{}-ZFP57".
                                           format(ZFP57_id_str)))

        # try to figure out the ID of the TF
        try:
            # NB: there are a number of edge-cases to verify here.
            #     Small changes to the regex can compromise some cases.
            #     Any issues will propagate and result in clear failures
            #     to find TF names corresponding to IDs or malformed IDs.
            # NB: it is critical to get the *last* ID in some cases,
            #     since some ENCODE IDs (from our metadata)
            #     may be <redacted version>-<new version> and only
            #     <new version> will correctly resolve. But this
            #     is also clear if it fails, as the ID fails to resolve
            #     to a TF name.
            IDgroups = re.search(r'(?:(?:Tel_Bulge_(Nfatc1))|'  # manual
                                 r'(ZFP57)|'  # manual
                                 r'MACS2_([^\W_]+)(?!.*ZFP57)|'
                                 # below matches only *last* ENCODE accession
                                 r'(ENC[A-Z]+[0-9]+[A-Z]+)'
                                 r'(?!.*ENC[A-Z]+[0-9]+[A-Z]+)|'
                                 #
                                 r'(?:(\w+?)_)?ChIP-?[Ss]eq(?:_(\w+?)_)?.*|'
                                 r'GSM\d+_(\w+?)(?:_.*)?|'
                                 r'(ERR\d+[-_]ERR\d+)'
                                 r')(?:-mod|[-_/]|$)',
                                 file_name_for_label).groups()
        except AttributeError:  # unable to get ID: use longer str with ID
            warn("ID error for {}. Using longer string.".format(file_name))
            IDgroups = re.search(r'\W(.+?)(?:-mod|$)',
                                 file_name_for_label).groups()

        # get the last defined group, since that is the ID we want
        ID = [id for id in IDgroups if id][-1].upper()

        # get the TF name, if annot. in the metadata, else use the ID alone
        # concat. w/empty Unicode string for StringIO-compliant initial value
        TF = ''
        for line in StringIO(u'' + metadata):
            if ID in line or ID.replace('_', '-') in line:
                TF = line.split()[-1].upper()

        if (not TF):
            # manually annotate some known mappings from GEO searches
            if ID == 'RFWT' or ID == 'RMWT':
                TF = 'RXRA'
            else:
                TF = ID

        try:
            # add an additional (replicate) label to the ID, if possible
            replicate_or_addtl_ID_matches = \
                re.search(r'(?:[^\W_]+-?Stringency-?)?'
                          r'(\w+?[-_]\d+[.-]\d+)(?:.*(GSM\d+))?'
                          # above last cap. group needed, as (GSM\d+) not in ID
                          # below crucial to diff. ZFP57 reps. WRT BC8 vs. CB9
                          r'(?:[-_]?MACS2_([^\W_]+))?',
                          file_name)

            # construct additional ID by concatenation of all capture groups
            replicate_or_addtl_IDs = [id_add for id_add in
                                      replicate_or_addtl_ID_matches.groups()
                                      if id_add is not None]

            # create the final ID from the ID, plus additional hyphenated IDs
            IDs = [ID] + replicate_or_addtl_IDs
            ID = '-'.join(IDs)
        except (StopIteration, AttributeError):
            # unable to get addtl. label: warn & cont.
            warn('Unable to obtain an additional (replicate) label for . '
                 'file: {}.'
                 'This should be fine, assuming that the input dataset '
                 'contains uniquely resolved IDs (like any ENCODE ID) '
                 'and does not contain multiple replicates of the '
                 'same ChIP-seq assay across multiple biological '
                 'replicates, without changing its main ID code.'.
                 format(file_name))

        stringency = re.search(r'(\d+\.\d+|[^\W_]+)Stringency',
                               file_name).groups()[-1]

        # remove the first comment character, since it is the header
        input_file_split_comments = open(file_name).read().split('#')

        if len(input_file_split_comments) > 1:
            # take the second element, since first is not sorted warning
            assert_or_die_with_msg('WARNING' in input_file_split_comments[1],
                                   'First element does not appear to'
                                   'contain the expected unsorted warning.')

            cur_TF_data = input_file_split_comments[2].split('\n')

            if num == 0:  # handle the header
                processed_input_file += '\t'.join(('ID', 'TF', 'Stringency',
                                                   cur_TF_data[0])) + '\n'
            else:
                # add LF before each new file; terminate prev. final record
                processed_input_file += '\n'

            processed_input_file += '\n'.join([('\t'.join((ID, TF,
                                                           stringency, line)))
                                               for line in cur_TF_data[1:]
                                               if line])
    return processed_input_file


def process_data(processed_input_file, id_cols, id_cols_ordered,
                 top_H_DataFrames):
    centrimo_data = pd.read_csv(StringIO(processed_input_file),
                                delim_whitespace=True)

    # compute TF superfamily and family
    centrimo_data['TF_superfam'], centrimo_data['TF_fam'] \
        = zip(*centrimo_data['TF'].map(getCuratedTFClassFamilyInfo))

    num_data_cols = _NUM_ALLDATA_COLS + 2  # account for added cols

    # sort by the E-value
    centrimo_data.sort_values('E-value', inplace=True)
    centrimo_data.reset_index(drop=True, inplace=True)  # use the sorted index

    # masks modified motifs (True iff motif is modified)
    # used for the background for the percentile computations
    mod_mask = pd.Series([any(mod_base in motif for
                          mod_base in cUtils.MOD_MAP)
                          for motif in centrimo_data['alt']])

    # NB: the below being the number of initial input lines is sufficient
    _MAX_ARRAY_SIZE = len(centrimo_data.index)  # max processed data entries

    grouped_Hs = centrimo_data.groupby(['TF',
                                        centrimo_data['alt'].
                                        apply(cUtils.translUnmodSeq)],
                                       sort=False)

    alldata = np.zeros((_MAX_ARRAY_SIZE,
                        num_data_cols)).astype(object, copy=False)

    alldata_top_1_indices = []
    alldata_top_3_indices = []

    cur_index = 0

    for i, (name, group) in enumerate(grouped_Hs):
        # lookup table for unmodified motifs to decide to emit a warning
        _seen_missing_unmod_motifs_for_warn_check = {}

        unmod_H = cUtils.translUnmodSeq(group['alt'].head(1).tolist()[0])

        unmod_H_data = group[group['alt'] == unmod_H].reset_index(drop=True)

        # ignore hypotheses without an unmod. counterpart in CentriMo results
        if unmod_H_data.empty:
            continue

        mod_H_data = group[group['alt'] != unmod_H].reset_index(drop=True)

        if mod_H_data.empty:  # also ignore hypotheses w/o a mod. counterpart
            continue

        for top_index in xrange(0, len(mod_H_data)):
            mod_hypothesis_pairs = mod_H_data.ix[top_index, :]
            unmod_hypothesis_pairs = \
                unmod_H_data[np.logical_and(unmod_H_data['Stringency'] ==
                                            mod_hypothesis_pairs['Stringency'],
                                            unmod_H_data['ID'] ==
                                            mod_hypothesis_pairs['ID']
                                            )]

            # must have uniquely paired the modified hypothesis with its
            # unmodified counterpart
            if unmod_hypothesis_pairs.empty:
                _warn_if_new_missing_unmod_motif(
                    _seen_missing_unmod_motifs_for_warn_check,
                    mod_hypothesis_pairs.ix['alt'],
                    mod_hypothesis_pairs.ix['ID'])
                continue
            elif len(unmod_hypothesis_pairs.index) > 1:
                die('Unable to correctly (uniquely) pair an unmodified-mod. '
                    'hypothesis. {} matching unmodified motifs were found. '
                    'Manual curation is required. The unpaired data '
                    'was: \n\n\n{}\n\n{}'.format(
                        len(unmod_hypothesis_pairs.index),
                        unmod_hypothesis_pairs, mod_hypothesis_pairs))

            # convert to Series
            unmod_hypothesis_pairs = unmod_hypothesis_pairs.squeeze()

            # skip if both are not (remotely) significant
            if ((unmod_hypothesis_pairs['log_adj_p-value'] >=
                 MAX_LOG_P_VALUE) and
                (mod_hypothesis_pairs['log_adj_p-value'] >=
                 MAX_LOG_P_VALUE)):
                warn("Skipping {} due to lack of sig. "
                     "[on mod. H: {}]".format(name, top_index))
                continue

            # compute the main modification pref. metric
            log_adj_pValue_ratio_top_results = \
                (unmod_hypothesis_pairs[['log_adj_p-value']] -
                 mod_hypothesis_pairs[['log_adj_p-value']])

            # compute percentiles of motif statistical significance
            # to label each hypothesis pair with its preferred motif's
            # quantile of significance, relative to all other motif's
            # in that category (unmod. or mod.)
            percentile_of_score_top_results = []

            for (log_adj_pValue_ratio,
                 (unmod_H_top_result_idx, unmod_H_top_result_row),
                 (mod_H_top_result_idx, mod_H_top_result_row)
                 ) in zip(log_adj_pValue_ratio_top_results,
                          _iterrowsOrIdentity(unmod_hypothesis_pairs),
                          _iterrowsOrIdentity(mod_hypothesis_pairs)):

                # compute quantile labels of statistical significance of the
                # preferred 'top' motif, relative to all motifs with that pref.
                # across all transcription factors and *all* elucidated motifs
                # NB: this also includes unpaired motifs or those excluded due
                #     to neither pair being significant. This is a minority of
                #     motifs, and will consistently result in higher quantiles.
                #
                # We use the log adjusted p-value for numerical stability.
                # This is sufficient, since E-value corrections are similar
                # (i.e., number of multiple tests is similar across dataset).
                if log_adj_pValue_ratio > 0:
                    # mod. pref. => use modified motif's p-value
                    test_data_of_interest = mod_H_top_result_row

                    bg_data_of_interest = centrimo_data[mod_mask]
                else:
                    # agnostic or unmod. pref. => use unmod. motif's p-value

                    test_data_of_interest = unmod_H_top_result_row

                    bg_data_of_interest = centrimo_data[~mod_mask]

                # all values must be positive to correctly compute percentile
                cur_score = np.absolute(test_data_of_interest
                                        ['log_adj_p-value'])
                bg_stat_sigs = np.absolute(bg_data_of_interest
                                           ['log_adj_p-value'])

                # compute the percentile of each input
                # NB: this may be slow; for details and possible improvement,
                # see http://stackoverflow.com/a/28577101 (percentile comp.)
                # uses the default 'rank' average percentage ranking of score
                (percentile_of_score_top_results.
                 append(stats.percentileofscore(bg_stat_sigs, cur_score)))

            def getResultToAppend():
                """Construct the combined dataframe, including the recently computed
                   columns. Uses the unmodified columns for most values, which
                   is sufficient, since we matched unmodified-modified motif
                   pairs."""

                result_ser = pd.concat([log_adj_pValue_ratio_top_results,
                                        unmod_hypothesis_pairs.
                                        iloc[_UNMOD_COLS],
                                        mod_hypothesis_pairs.iloc[_MOD_COLS],
                                        pd.
                                        Series(percentile_of_score_top_results)
                                        ],
                                       ignore_index=True)
                return pd.DataFrame(result_ser).T

            alldata[cur_index, :] = getResultToAppend()

            if top_index == 0:
                alldata_top_1_indices.append(cur_index)

            if top_index < 3:
                alldata_top_3_indices.append(cur_index)

            cur_index += 1

    # not +1 of cur_index because it was incremented above
    alldata = alldata[:cur_index, :]

    column_names = ['diff-logAdjP-value']
    column_names.extend(unmod_H_data.columns[_UNMOD_COLS])
    column_names.extend(mod_H_data.columns[_MOD_COLS])
    column_names.extend(['SigQuantiles'])

    df_alldata = pd.DataFrame(alldata, columns=column_names)

    if args.verbose:
        pd.cut(df_alldata['SigQuantiles'],
               NUM_QUANTILES).to_csv(_ACTUAL_SIG_QUANTILE_CUTS_FILENAME,
                                     index=False)

    # discretize percentiles to q-quantiles (Seaborn hues must be categorical)
    # Like other id columns, this column is also later converted to Categorical
    df_alldata.loc[:, 'SigQuantiles'] = pd.cut(df_alldata['SigQuantiles'],
                                               NUM_QUANTILES,
                                               labels=DECILE_ORD_LABELS)

    # add a column annotating the modification type, for each motif
    df_alldata.loc[:, 'Modification'] = \
        df_alldata.iloc[:, 6].map(getModType)

    # populate result dict of DataFrames, which will be returned once fixed
    top_H_DataFrames['all'] = df_alldata

    # create view of the underlying N-hypothesis DataFrame for top 1 hypothesis
    top_H_DataFrames['1'] = df_alldata.iloc[alldata_top_1_indices, :]

    # create view of the underlying N-hypothesis DataFrame for top 3 hypotheses
    top_H_DataFrames['3'] = df_alldata.iloc[alldata_top_3_indices, :]

    # reformat data
    for num_top_Hs, cur_df in top_H_DataFrames.iteritems():
        fixed_cur_df = \
            pd.DataFrame(pd.concat([pd.to_numeric(cur_df.loc[:,
                                    ['diff-logAdjP-value']].iloc[:, 0]),
                                    (cur_df.iloc[:, 5] + '-' +
                                     cur_df.iloc[:, 6]),
                                    cur_df.loc[:, ['ID']],
                                    cur_df.loc[:, ['TF']],
                                    cur_df.loc[:, ['Stringency']],
                                    cur_df.loc[:, ['Modification']],
                                    cur_df.loc[:, ['SigQuantiles']],
                                    cur_df.loc[:, ['TF_superfam']],
                                    cur_df.loc[:, ['TF_fam']]
                                    ], axis=1))

        cols = [fixed_cur_df.columns[0]] + id_cols
        fixed_cur_df.columns = cols

        fixed_cur_df = pd.melt(fixed_cur_df, id_vars=id_cols,
                               value_vars='diff-logAdjP-value',
                               value_name='diff-logAdjP-value')
        del fixed_cur_df['variable']  # remove extra col. w/ only var. name

        # sort by TF superfamily and then by family, and finally TF name
        fixed_cur_df.sort_values(['TF_superfam', 'TF_fam', 'TF'], inplace=True)

        fixed_cur_df.reset_index(drop=True, inplace=True)  # use sorted index

        top_H_DataFrames[num_top_Hs] = fixed_cur_df

    return top_H_DataFrames


top_H_DataFrames = {'1': pd.DataFrame(), '3': pd.DataFrame(),
                    'all': pd.DataFrame()}

# format column types
id_cols = ['hypothesis', 'ID', 'TF', 'Stringency', 'Modification',
           'SigQuantiles', 'TF_superfam', 'TF_fam']

# a subset of the above id_cols
# these categoricals have a natural ordering that should be observed
id_cols_ordered = ['Stringency', 'SigQuantiles']

if args.plot_only_from_CSVs:
    # we do not create an all plot in this case
    top_H_DataFrames.pop('all')

    for num_top_Hs in top_H_DataFrames.keys():
        top_H_DataFrames[num_top_Hs] = \
            pd.DataFrame.from_csv("{}{}{}".format(
                args.plot_only_from_CSVs, num_top_Hs,
                _DEFAULT_CSV_PATH_FILENAME_EXT),
                index_col=None, sep="\t", parse_dates=False)
else:
    processed_input_file = retrieve_input_data(METADATA_FILE_FULLPATH,
                                               args.root_input_dirname)

    top_H_DataFrames = process_data(processed_input_file, id_cols,
                                    id_cols_ordered, top_H_DataFrames)

for num_top_Hs, df_for_plot in top_H_DataFrames.iteritems():
    possible_primary_mod_bases = df_for_plot.Modification.unique()

    # convert ID columns to ordered categorical types (from generic object)
    for col in id_cols:
        ordered_cat = True if col in id_cols_ordered else False

        df_for_plot[col] = df_for_plot[col].astype('category',
                                                   ordered=ordered_cat)

    # compute labels for transcription factor family annotations

    prev_sfam = ''
    prev_fam = ''

    sfams = []
    fams = []

    sfam_chng_idxs = []
    fam_chng_idxs = []

    for lab_idx, lab in enumerate(df_for_plot['TF'].unique()):
        sfam = str(df_for_plot[df_for_plot['TF'] == lab]
                   ['TF_superfam'].unique()[0])

        if sfam != prev_sfam:
            sfams.append(str(sfam))
            sfam_chng_idxs.append(lab_idx)

        fam = str(df_for_plot[df_for_plot['TF'] == lab]['TF_fam'].unique()[0])

        if fam != prev_fam:
            fams.append(str(fam))
            fam_chng_idxs.append(lab_idx)

        prev_sfam = sfam
        prev_fam = fam

    # add the final index for the last (super)family
    sfam_chng_idxs.append(lab_idx + 1)
    fam_chng_idxs.append(lab_idx + 1)

    sfam_labs = [(sfam_e, (idx_0, idx_1)) for (idx_0, idx_1), sfam_e
                 in zip(zip(sfam_chng_idxs[:-1], sfam_chng_idxs[1:]), sfams)]
    fam_labs = [(fam_e, (idx_0, idx_1)) for (idx_0, idx_1), fam_e
                in zip(zip(fam_chng_idxs[:-1], fam_chng_idxs[1:]), fams)]

    # ------------------ output data to file ------------------
    if not args.plot_only_from_CSVs:
        df_for_plot.to_csv('{}{}{}'.format(
                           _DEFAULT_CSV_PATH_FILENAME_BASE,
                           num_top_Hs,
                           _DEFAULT_CSV_PATH_FILENAME_EXT),
                           sep=b'\t', encoding='utf-8', index=False)

    if not args.no_plot:
        if not args.keep_mixed_mod_motifs:
            # remove all motifs containing mixed modifications
            mixed_mod_mask = \
                (df_for_plot['Modification'].
                 isin(cUtils.AMBIG_MOD_BASES.keys()))

            if args.verbose:
                warn("Top {}: Removed {} motifs with mixed modifications "
                     "(of {} in total)"
                     "".format(num_top_Hs, (mixed_mod_mask == 1).sum(),
                               len(df_for_plot)))

                if args.verbose > 1:
                    excluded_motifs_out_filename = \
                        '{}{}{}'.format(
                            _EXCLUDED_MIXED_MOD_MOTIFS_FILENAME_BASE,
                            num_top_Hs, _DEFAULT_CSV_PATH_FILENAME_EXT)

                    warn("The removed mixed modifications have "
                         "been saved to {}."
                         "".format(excluded_motifs_out_filename))

                    (df_for_plot.loc[mixed_mod_mask].
                     to_csv(excluded_motifs_out_filename, sep=b'\t',
                            encoding='utf-8', index=False))

            df_for_plot = df_for_plot.loc[~mixed_mod_mask]

        # ---------- compute additional right axis annotations ----------

        # compute max and min per TF, for both m and h
        df_for_plot_for_agg = df_for_plot.groupby(['TF', 'Modification'])

        max_vals = df_for_plot_for_agg.max()['diff-logAdjP-value']
        min_vals = df_for_plot_for_agg.min()['diff-logAdjP-value']

        # only use those on each side of 0 and 0 with max
        max_vals_use = max_vals.where(max_vals >= 0)
        min_vals_use = min_vals.where(min_vals < 0)

        # reformat to:
        #                         min               max
        # Modification         m        h        m       h
        val_annots = (pd.concat((min_vals_use, max_vals_use),
                                axis=1, keys=['min', 'max']).unstack().
                      sort_index(axis=1, ascending=[1, 0]).
                      reindex(df_for_plot['TF'].unique()))

        val_annots_to_plot = re.sub(r"\\quad{}\\makebox(.+?)"
                                    r"\\quad{}\\makebox",
                                    # manually re-add matched LaTeX cmds
                                    # to prevent doubling of "\"
                                    r"\quad{}\makebox\1\qquad{}\qquad{}"
                                    r"\makebox",
                                    re.sub("[\t ]+", "\quad{}",
                                           val_annots.
                                           to_string(header=False,
                                                     index=False,
                                                     float_format="@${:,.0f}$"
                                                                  "^".format).
                                           replace("@", "\makebox[{}cm][l]"
                                                   "{{".format(
                                                       PLOT_DEFAULT_R_AN_WIDTH)
                                                   ).
                                           replace("^", "}").
                                           replace("$nan$", "    "))
                                    ).split('\n')

        # ------------------ plotting ------------------

        sns.set_style("darkgrid")

        sns.set(PRESET_SNS_CONTEXT, {"xtick.major.size": PLOT_X_TICK_SIZE,
                                     "xtick.minor.size":
                                     PLOT_X_TICK_SIZE // 2},
                font_scale=PLOT_FONT_SCALE_FACTOR)

        # default to a swarm plot, to best show value distribution
        plot_type = _SEABORN_SWARM_PLOT_TYPE

        # switch to a violin plot for many points, as swarm plot scales poorly
        if len(df_for_plot) > _MAX_VALS_FOR_SWARMPLOT:
            plot_type = _SEABORN_VIOLIN_PLOT_TYPE

        plot_hues = ['SigQuantiles', 'Stringency', 'Modification']

        assert_or_die_with_msg(set(plot_hues).issubset(df_for_plot.columns))

        for hue in plot_hues:
            # global variable to track the number of spines added
            num_addtl_spines = 0

            plot_top_param = PLOT_DEFAULT_TOP_PARAM
            plot_r_an_width = PLOT_DEFAULT_R_AN_WIDTH

            plt.figure(figsize=(PLOT_WIDTH, _compute_plot_height(df_for_plot.
                                                                 TF.nunique())
                                ))

            hue_order = None
            palette = PLOT_DEFAULT_PALETTE

            if hue == 'Modification':
                hue_order = [base for base in
                             cUtils.MOD_BASE_COMPLEMENT_NUM_ORDER
                             if base in df_for_plot.Modification.unique()]

                if not hue_order:
                    warn("Provided \"Modification\" column did not contain "
                         "modified bases. Assuming replicate labels instead.")

                    df_for_plot.rename(columns={hue: 'Replicate'},
                                       inplace=True)

                    hue = 'Replicate'

                    hue_order = sorted(df_for_plot[hue].unique())
            elif hue == 'SigQuantiles':
                hue_order = [ord_label for ord_label in DECILE_ORD_LABELS
                             if ord_label in df_for_plot.SigQuantiles.unique()]

                # use a sequential cubehelix color palette
                # s.t. more sig. motifs are darker, and sequentiality
                # is preserved under B&W output or colour blindness
                palette = sns.cubehelix_palette(len(DECILE_ORD_LABELS),
                                                start=.5, rot=-.75)

            num_categories = df_for_plot[hue].nunique()

            split_hue_on_categorical_axis = False
            if num_categories == 2:
                # it only makes sense to do this when there are
                # only a couple of categories, like for two
                # stringencies or two different kinds of modifications
                split_hue_on_categorical_axis = True

            plot_ax = \
                getattr(sns, plot_type)(x="diff-logAdjP-value", y="TF",
                                        split=split_hue_on_categorical_axis,
                                        order=df_for_plot['TF'].unique(),
                                        palette=palette, hue=hue,
                                        hue_order=hue_order, data=df_for_plot)

            plot_ax.set(xlabel=PLOT_X_AXIS_LABEL, ylabel=PLOT_Y_AXIS_LABEL)

            plot_ax.get_xaxis().set_minor_locator(mpl.ticker.
                                                  AutoMinorLocator())
            plot_ax.grid(which='major', color='w', linewidth=2)
            plot_ax.grid(which='minor', color='w', linewidth=1)

            sns.despine()

            # ---------- amend display of some TF labels ----------
            TF_labs = [y_lab.get_text()for y_lab in plot_ax.get_yticklabels()]

            TF_labels_for_disp = edit_TF_labels_for_display(TF_labs)

            plot_ax.set_yticklabels(TF_labels_for_disp)

            # ---------------------------------------------
            # add TF family y-axis labels (additional spines)
            plot_ax.spines['left'].set_position(('outward',
                                                 PLOT_LEFT_SPINE_DROP))

            # add spines in order, from those closest to the y-axis
            add_spine_annot(plot_ax, 'Family', fam_labs)
            add_spine_annot(plot_ax, 'Superfamily', sfam_labs)
            # ---------------------------------------------
            annot_ax = plot_ax.twinx()

            # clone original y-axis
            annot_ax.set_ylim(plot_ax.get_ylim())
            annot_ax.set_yticks(plot_ax.get_yticks())

            # disable ticks (already have them on original axis)
            annot_ax.grid('off')

            # replace labels with the annotations
            annot_ax.set_yticklabels(val_annots_to_plot)

            add_h_header_label = True

            # base amount for legend to protrude, if placing in margin
            legend_out_base_dist = PLOT_DEFAULT_LEGEND_OUT_BASE_DIST

            # do not add the 5hmC header portion if we only have 5mC
            if (len(possible_primary_mod_bases) == 1 and
                    possible_primary_mod_bases[0] == 'm'):
                add_h_header_label = False
                plot_top_param += 0.01  # less room is needed at the top
                plot_r_an_width /= 1.7  # move min/max labels closer
            else:
                legend_out_base_dist *= 2  # move legend out more (extra cols)

            adjust_plot_right = None
            if num_categories > 4:
                # move the legend outside of the plot
                plot_ax.legend(bbox_to_anchor=(1 + PLOT_LEGEND_OUT_ADDTL_DIST +
                                               legend_out_base_dist, 1), loc=2,
                               borderaxespad=0.)
                adjust_plot_right = PLOT_RIGHT_PARAM

            # Use a tight layout and then expand for the TF family annotations
            plt.tight_layout()
            plt.subplots_adjust(left=PLOT_LEFT_PARAM, right=adjust_plot_right,
                                top=plot_top_param)

            header_label = ("\enskip{{}}\makebox[{1}cm][c]{{$\min_{{< 0}}$}}"
                            "\quad{{}}")

            if add_h_header_label:
                header_label += "\qquad{{}}"

            header_label += "\makebox[{1}cm][c]{{$\max_{{\geq 0}}$}}"

            if add_h_header_label:
                header_label += ("\n\quad{{}}\makebox[{0}cm][l]{{m}}"
                                 "\quad{{}}\makebox[{0}cm][l]{{h}}\qquad{{}}"
                                 "\qquad{{}}\enskip{{}}\makebox[{0}cm][l]{{m}}"
                                 "\quad{{}}\makebox[{0}cm][l]{{h}}")

            header_label += "\n\n"

            # add the header label for these annotations
            annot_ax.annotate(header_label.format(plot_r_an_width,
                                                  2*plot_r_an_width),
                              xy=(1, 1.001), xycoords='axes fraction')
            # ---------------------------------------------

            # remove the legend if there is only one item
            if len(plot_ax.get_legend_handles_labels()[0]) < 2:
                plot_ax.legend().set_visible(False)

            try:
                for ext in _EXTENSIONS_TO_OUTPUT:
                    plt.savefig("{}-top{}-{}.{}".format(BASE_PLOT_NAME,
                                                        num_top_Hs, hue, ext))
            except RuntimeError as ex:
                warn("Unable to save figure. The encountered LuaLaTeX "
                     "error is likely due to attempting to plot too few "
                     "TFs, which this script is not designed for. "
                     "The RuntimeError was:\n\n {}".format(ex))

                # if one fails, all formats almost certainly will
                break
