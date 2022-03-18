#!/usr/bin/env python

"""Converts MethPipe methcounts, MLML output, or other similarly
   formatted genomic interval with modification level data into
   bedGraph files, suitable for downstream processing with Cytomod.
"""

from __future__ import with_statement, division, print_function

import argparse
from collections import OrderedDict
import errno
import os
import sys
import textwrap

import numpy as np
import pandas as pd


# From: http://stackoverflow.com/a/22157136
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)


# the number of lines to read into memory at a time
_CHUNK_SIZE = 10000000

_STDIN_SPECIFIER = '-'
_MAX_CONFLICTS_VALUE = 2  # determined by the initial number of data types
_DEFAULT_OUTPUT_SUFFIX = '.bedGraph.gz'

# NB: final column should be 'i4', but integers cannot be NaN.
# This is a known limitation of Pandas (from NumPy).
FORMAT_MLML = OrderedDict([('chr', 'S12'), ('start', np.uint32),
                           ('end', np.uint32), ('5mC', np.float16),
                           ('5hmC', np.float16), ('C', np.float16),
                           ('conflicts', np.float16)])

FORMAT_METHCOUNTS_LIKE = OrderedDict([('chr', 'S12'), ('start', np.uint32),
                                      ('5mC', np.float16)])


class Range(object):
    """Adapted from: http://stackoverflow.com/a/12117089"""
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __str__(self):
        return "[{}, {}]".format(self.start, self.end)

    def __unicode__(self):
        return str(self)

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __contains__(self, other):
        return self.start <= other.start and other.end <= self.end


def _CSVAppendDF(df, file_fullpath, allow_header=True):
    """Append to a tab-separated, Gzipped, file.
       Adapted from: http://stackoverflow.com/a/30292938
    """

    header = True if allow_header else False

    if os.path.isfile(file_fullpath):
        header = False

        existing_CSV_cols = pd.read_csv(file_fullpath, nrows=1,
                                        sep='\t').columns

        if len(df.columns) != len(existing_CSV_cols):
            raise RuntimeError("Number of DataFrame ({}) and existing "
                               "CSV ({}) columns ({}) do not match"
                               ".".format(len(df.columns), file_fullpath,
                                          len(existing_CSV_cols)))
        # only compare headers if we are using a header
        elif (allow_header and not
              np.array_equal(df.columns, existing_CSV_cols)):
            raise RuntimeError("Columns/column order mismatch "
                               "(between the DataFrame ({}) and existing CSV"
                               " ({}; {}).".format(df.columns, file_fullpath,
                                                   existing_CSV_cols))

    df.to_csv(file_fullpath, sep='\t', index=False, mode='a',
              compression='gzip', header=header)


def _getFilename(nucleobase, prefix, suffix=_DEFAULT_OUTPUT_SUFFIX):
    return "{}/5{}C{}".format(prefix, nucleobase, suffix)


parser = argparse.ArgumentParser(formatter_class=SmartFormatter)

parser.add_argument('input', type=str, help=textwrap.fill(
                    "R|Tab-delimited"
                    "(optionally gzipped) input file containing the "
                    "(hydroxy)methylation data. Alternatively, a '{}' can "
                    "be used to read from STDIN. "
                    "NB: Assumes that genome coordinates are no greater "
                    "than 4294967295 (uint32) and that all float values are "
                    "no greater than 5 bits exponent, 10 bits mantissa "
                    "(float16). Assumes all reference names (i.e. \"chr..\") "
                    "fit into a five character string (S5). "
                    "The file's columns must be in one of two formats, "
                    "where '.' represents ignored columns:\n\n"
                    "chr\tstart\tend\t5mC-level\t5hmC-level\tC-level conflicts"
                    "\n\nor\n"
                    "chr\tstart\t.\t.\t5mC-level\t.".format(_STDIN_SPECIFIER),
                    55, drop_whitespace=False, replace_whitespace=False))
parser.add_argument('threshold', type=float, choices=[Range(0.0, 1.0)],
                    help="The threshold to call a modified base.")
parser.add_argument('--max-conflicts', type=int,
                    choices=range(0, _MAX_CONFLICTS_VALUE + 1), default=0,
                    help="The maximum number of MLML conflicts permitted "
                         "before setting a nucleobase as having unknown "
                         "modification state. "
                         "Ignored if only 5mC data is present "
                         "(e.g. methcounts output alone).")
parser.add_argument('-o', '--outDir', type=str, default='.',
                    help="The directory in which to place the output "
                          "bedGraph files. The output directory should "
                          "not already contain these files, since output "
                          "will be appended.")
parser.add_argument('-s', '--suffix', type=str, default=_DEFAULT_OUTPUT_SUFFIX,
                    help="A suffix to append to the output files.")
parser.add_argument('-0', '--includeBelowThreshAsZeroMethyl',
                    action='store_true',
                    help="By default, only include calls passing the "
                         "threshold within any output BED file. If this "
                         "option is set, however, include all calls."
                         "Those not passing the threshold are included "
                         "within the 5xC output file, with a score of 0.")
parser.add_argument('--expand-dinuc', action='store_true',
                    help="Assume that input intervals are the final base"
                    "of a dinucleotide (e.g. CpG), and expand it to"
                    "include the first base of the dinucleotide."
                    "This option is not strand-specific (always -1)."
                    "The use of this option, results in the start"
                    "coordinate of every output interval to be one less"
                    "than what it would have been otherwise.")

args = parser.parse_args()

input_file = sys.stdin if args.input == _STDIN_SPECIFIER else args.input

format = None
usecols = None

has_methylation_only = False

with open(input_file, 'r') as input_file_h:
    if len(input_file_h.readline().strip().split()) == 7:
        dtype = FORMAT_MLML
    else:
        dtype = FORMAT_METHCOUNTS_LIKE
        usecols = [0, 1, 4]
        has_methylation_only = True

df_chunks = pd.read_csv(input_file, usecols=usecols, engine='c',
                        dtype=dtype, chunksize=_CHUNK_SIZE, sep='\t',
                        index_col=False, names=dtype.keys())

# EFAP directory creation. Adapted from: http://stackoverflow.com/a/5032238
try:
    os.mkdir(args.outDir)
except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

for df_chunk in df_chunks:
    if has_methylation_only:
        possible_bases = ['x']

        cond_list = [df_chunk['5mC'] >= args.threshold]

        # define the end as the usual start + 1, since not in input
        df_chunk['end'] = df_chunk['start'] + 1
    else:
        possible_bases = ['z', 'x', 'hm', 'm']

        # NB: both cond_list and the given choice list must be ordered
        #     from most to least ambigous
        #     This is because the first true condition in this list
        #     will be used, even if multiple conditions are true.
        cond_list = [df_chunk['conflicts'] > args.max_conflicts,
                     np.logical_and(df_chunk['5mC'] >= args.threshold / 2,
                                    df_chunk['5hmC'] >= args.threshold / 2),
                     np.logical_and(df_chunk['5hmC'] >= args.threshold,
                                    df_chunk['5mC'] < args.threshold),
                     np.logical_and(df_chunk['5mC'] >= args.threshold,
                                    df_chunk['5hmC'] < args.threshold)]
        # NB: the known unmodified case:
        # '5mC' < args.threshold && '5hmC' < args.threshold: C
        # is ignored, since it does not require its own file.

    if args.expand_dinuc:
        # expand to include first base of dinucleotide
        # this can never be out-of-bounds, since the start
        # was the *second* base of a dinucleotide
        # NB: not strand-specific
        df_chunk['start'] -= 1

    output_file_names = [_getFilename(base, args.outDir, args.suffix)
                         for base in possible_bases]

    df_chunk['out_file'] = np.select(cond_list, output_file_names)

    df_chunk_grouped = df_chunk.groupby('out_file')

    # output to respective file for all modified loci
    for output_file_name, mod_loci in df_chunk_grouped:
        # disable false pos. SettingWithCopyWarning that would o/w result
        mod_loci.is_copy = False

        if output_file_name is not '0':
            # only file names that are not '0' rep. passing calls to output

            # Cytomod needs a final numeric column and includes all
            # intervals with non-zero values in that column.
            # These intervals should all be included, so set this to unity.
            mod_loci.loc[:, 'include'] = 1
        elif args.includeBelowThreshAsZeroMethyl:
            # intervals should not be included, so set this to nullity
            mod_loci.loc[:, 'include'] = 0

            # null regions are included only within the 5xC output
            output_file_name = _getFilename('x', args.outDir, args.suffix)
        else:
            # no output if null and not outputting null calls
            continue

        # do not permit a header, since Cytomod does not permit it
        # and it contravenes the bedGraph specification
        _CSVAppendDF(mod_loci[['chr', 'start', 'end', 'include']],
                     output_file_name, False)
