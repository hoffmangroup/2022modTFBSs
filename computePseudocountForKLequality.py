#!/usr/bin/env python

from __future__ import with_statement, division, print_function

import math
import textwrap
import sys

from scipy.optimize import fsolve

USAGE_MSG = """\
            This program computes a pseudocount s.t. its application forces
            the Kullback-Leibler (KL) divergences of two motifs, for which the
            only difference is their modification state, to be equal.
            This is done given only the provided background frequencies.

            Usage: {}
            <number of modified bases>
            <background frequency for unmodified base on one strand>
            <background frequency for modified base on one strand>
            [<background frequency for unmodified base on the other strand>]
            [<background frequency for modified base on the other strand>]

            Only the first two background frequencies are required, while
            the other two are optional, but must both be provided if one is.
            Prints (to STDOUT) the pseudocount s.t. columnar KL divergences
            will be equal, given the background.
            """.format(sys.argv[0])

PSEUDOCOUNT_GUESS = 5  # the initial value for the optimization


def lg(x):
    """Returns the base 2 log of real number x."""

    return math.log(x, 2)


def KLdivergenceOpt(pseudocount, bg1, bg2, bg3='', bg4='', one_pseud=''):
    """Assumes that bases under consideration both have profile
    frequencies of 1.

    Computes the pseudo-count parameter s.t. the KL divergence, caused by
    the difference in background frequencies for the particular nucleobase
    under consideration, is equalized.

    Assumes that modified positions all correspond to modification of the
    same nucleobase (or its complement for complete modification)
    and that the target modifications are all the same (i.e. all 5mC).
    This pseudocount value, when set for one motif (optionally specified as
    "one_pseud" or, if not provided, assumed to be part of the optimization
    (i.e. set to a free parameter like the other pseudocount, assuming both
    motifs will be set to the output pseudocount provided), equalizes the KL
    divergences of both motifs, accounting for the pseudocount applied
    to the other motif.
    """

    if not one_pseud:
        one_pseud = pseudocount

    def singleKL(bg, pseudocount=pseudocount):
        return (1 + pseudocount * bg) * lg((1 + pseudocount * bg) / bg)

    # each side has two terms to permit both complete or hemi-modification
    # wherein the complements differ in background frequency from mod bases
    lhs = singleKL(bg1) + (singleKL(bg3) if bg3 else 0)
    rhs = singleKL(bg2, one_pseud) + (singleKL(bg4, one_pseud) if bg4 else 0)
    return lhs - rhs


if len(sys.argv) < 4 or len(sys.argv) > 7:
    sys.exit(textwrap.dedent(USAGE_MSG))

# the final pseudocount is multiplied by this to ensure correct scaling for
# multiple modifications
num_mod_bases = int(sys.argv[1])

background_unmod_s1 = float(sys.argv[2])
background_mod_s1 = float(sys.argv[3])

KL_args = (background_unmod_s1, background_mod_s1)

if len(sys.argv) > 5:
    if not num_mod_bases > 1:
        ERR_MSG = """\
                  Only one set of background frequencies should be provided
                  if only a single base is being modified. Otherwise, the
                  actual number of modified bases needs to be specified."""

        sys.exit(textwrap.dedent(ERR_MSG))

    background_unmod_s2 = float(sys.argv[4])
    background_mod_s2 = float(sys.argv[5])
else:
    if len(sys.argv) == 5:
        WARN_MSG = """\
                  -----------------------------------------------------
                  Assuming that a second set of background frequencies
                  was not intended to be provided, but that the final
                  provided parameter, pertains to the pseudocount of
                  a second motif.
                  -----------------------------------------------------
                  """

        print(textwrap.dedent(WARN_MSG), file=sys.stderr)

    background_unmod_s2 = ''
    background_mod_s2 = ''

KL_args += (background_unmod_s2, background_mod_s2)

if len(sys.argv) % 2 == 0:
    KL_args += (float(sys.argv[-1]),)

sol = fsolve(KLdivergenceOpt, x0=PSEUDOCOUNT_GUESS,
             args=KL_args)

# ensure that we obtained the expected multiplicity (one root)
assert len(sol) == 1

single_pseudocount = sol[0]

# Output the final result: the computed pseudocount times
#                          the number of modified bases.
#
# This product is the pseudocount suitable for input to CentriMo
# for it to normalize the KL divergence, relative to the unmodified motif,
# whose pseudocount is expected to be 0.
print(num_mod_bases * single_pseudocount)
