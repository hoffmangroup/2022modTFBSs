#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd

MAIN_FONT = "Liberation Sans"

# ------------------------------------------------------
import matplotlib as mpl

mpl.use('pgf')

mpl.rc('font',
       **{'family': 'sans-serif',
          'sans-serif': [MAIN_FONT],
          'monospace': ['Liberation Mono']})
# use the STIX sans-serif font for mathematics
mpl.rc('mathtext', fontset='stixsans')

mpl.rc('pgf', texsystem='lualatex')  # use LuaLaTeX for non-strict mem.

from matplotlib.backends.backend_pgf import FigureCanvasPgf
mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)

import matplotlib.pyplot as plt

#from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize, SymLogNorm

import seaborn as sns

sns.set(font=MAIN_FONT)
# ------------------------------------------------------

import os
import re

_R_LIB_DIR = "/scratch/R/x86_64-redhat-linux-gnu-library/3.5"

import rpy2.robjects as ro
from rpy2.robjects.packages import importr

from rpy2.robjects import pandas2ri
pandas2ri.activate()

from rpy2.robjects.lib import grid

from rpy2.robjects.lib import grdevices

from rpy2.rinterface import NAIntegerType

base = importr(str('base'), lib_loc=_R_LIB_DIR)

scales = importr(str('scales'), lib_loc=_R_LIB_DIR)

fontLib = importr(str('fontLiberation'), lib_loc=_R_LIB_DIR)
extrafont = importr(str('extrafont'), lib_loc=_R_LIB_DIR)
extrafont.loadfonts(device='pdf', quiet=True)

# N.B. capital "grDevices", as we also need "grdevices" to call display methods
# Object is instead used to call *functions* def. in R's grdevices namespace
grDevices = importr(str('grDevices'), lib_loc=_R_LIB_DIR)

# need this before below, to find the package
ggplot = importr(str('ggplot2'), lib_loc=_R_LIB_DIR)

import rpy2.robjects.lib.ggplot2 as ggplot2

treemapify = importr(str('treemapify'), lib_loc=_R_LIB_DIR)

# N.B. applies only to the heatmap. Treemap uses a Pandas function.
# N.B. aggregating over other variables
# Aggregation functions: nanmin, nanmax, nanmean
# Using max, to emphasize modified DNA preferences.
_AGG_FUNC = np.nanmax

TM_AREA_COL_NAME = 'plot_area'

_TICK_LABEL_SIZE = 16

_COEFF_SCALE_FACTOR = 4
_AREA_SCALE_FACTOR = 20

_PLOT_SIZE = 40
_PLOT_COLS_PER_ROW = 3

_TICK_LABEL_SIZE = 20

# Name for the diff-logAdjP-value / colour bar
_PLOT_VAL_TITLE = "Score"

_PLOT_COLORBAR_CLUST_PAD = 1.2
_PLOT_COLORBAR_CLUST_SIZE = '2.5%'

_PLOT_COLORBAR_SCORE_SIZE = '1.5%'
_PLOT_COLORBAR_SCORE_PAD = 120

_PLOT_CBAR_SHRINK_FACTOR = 0.4
_PLOT_CBAR_BOT_ANCHOR = -0.16

_PLOT_FACET_SPACING = 0.4

_PLOT_EXTRA_BOT_MARGIN = 0.08

_TF_SYNONYM_MAP = {"POU5F1": "OCT4", "EGFP-": "", "KAISO": "ZBTB33"}

_TFs = ["JUND", "FEZF2", "USF1",
        "JUN", "CEBPB", "USF2",
        "OCT4", "ZNF143", "TCF12",
        "KLF4", "ZC3H11A", "ETS1",
        "ZNF384", "YY1", "CTCF",
        "ZFP57", "MYC", "ZBTB33"]

_DIRECTORIES = ["Yin_data", "ENCODE_mouse", "ENCODE_K562"]

_FILE_BASE = "main_scores_clusters"

_FILE_EXTENSION = "csv"

class OffsetSymLogNorm(SymLogNorm):
    """
    Adapted from the Matplotlib GitHub:
    https://github.com/OceanWolf/matplotlib/commit/05c1b76af7a1f3c170a3423546882ce4918f8a66
    The same as their OffsetNorm,
    except subclassing matplotlib.colors.SymLogNorm

    Normalizes data into the ``[0.0, 1.0]`` interval.
    """
    def __init__(self, linthresh, linscale=1.0,
                 vmin=None, vcenter=None, vmax=None, clip=False):
        """Normalize data with an offset midpoint

        Useful when mapping data unequally centered around a conceptual
        center, e.g., data that range from -2 to 4, with 0 as the midpoint.

        Parameters
        ----------
        linthresh: required float
            The range within which the plot is linear (to
            avoid having the plot go to infinity around zero).

        linscale : optional float
            This allows the linear range (-*linthresh* to *linthresh*)
            to be stretched relative to the logarithmic range.  Its
            value is the number of decades to use for each half of the
            linear range.  For example, when *linscale* == 1.0 (the
            default), the space used for the positive and negative
            halves of the linear range will be equal to one decade in
            the logarithmic range. Defaults to 1.
        
        vmin : optional float
            The data value that defines ``0.0`` in the normalized data.
            Defaults to the min value of the dataset.

        vcenter : optional float
            The data value that defines ``0.5`` in the normalized data.
            Defaults to halfway between *vmin* and *vmax*.

        vmax : option float
            The data value that defines ``1.0`` in the normalized data.
            Defaults to the the max value of the dataset.

        clip : optional bool (default is False)
            If *clip* is True, values beyond *vmin* and *vmax* will be set
            to ``0.0`` or ``1.0``, respectively. Otherwise, values outside
            the ``[0.0, 1.0]`` will be returned.
        """

        SymLogNorm.__init__(self, linthresh, linscale, vmin, vmax, clip)

        self.vcenter = vcenter


    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vcenter, vmax = self.vmin, self.vcenter, self.vmax
        if vmin == vmax == vcenter:
            result.fill(0)
        elif not vmin <= vcenter <= vmax:
            raise ValueError("minvalue must be less than or equal to "
                             "centervalue which must be less than or "
                             "equal to maxvalue")
        else:
            vmin = float(vmin)
            vcenter = float(vcenter)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = self._transform(result.data)

            # first scale to -1 to 1 range, then from 0 to 1.
            resdat -= self._center
            resdat[resdat > 0] /= abs(self._upper - self._center)
            resdat[resdat < 0] /= abs(self._lower - self._center)

            resdat /= 2.
            resdat += 0.5
            result = np.ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]

        return result


    def _transform_points(self):
        """
        Calculates vmin, vcenter, and vmax in the transformed system.
        """
        vmin, vcenter, vmax = self.vmin, self.vcenter, self.vmax

        arr = np.array([vmin, vcenter, vmax]).astype(np.float)

        self._upper, self._center, self._lower = self._transform(arr)


    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")

        vmin = float(self.vmin)
        vcenter = float(self.vcenter)
        vmax = float(self.vmax)

        up_diff = abs(self._upper - self._center)
        low_diff = abs(self._lower - self._center)

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val - 0.5)

            val[val > 0] *= up_diff
            val[val < 0] *= abs(self._lower - self._center)

            val += self._center

            return self._inv_transform(val)
        else:
            val = 2 * (val - 0.5)

            if val < 0:
                return self._inv_transform(val * low_diff + self._center)
            else:
                return self._inv_transform(val * up_diff + self._center)

            
    def autoscale_None(self, A):
        ' autoscale only None-valued vmin or vmax'
        if self.vmin is None and np.size(A) > 0:
            self.vmin = ma.min(A)

        if self.vmax is None and np.size(A) > 0:
            self.vmax = ma.max(A)

        if self.vcenter is None:
            self.vcenter = (self.vmax + self.vmin) * 0.5

        self._transform_points()


def format_data_for_plot(data, values, index, columns, aggfunc):
    df_for_plot = data.pivot_table(values=values, index=index,
                                   columns=columns, aggfunc=np.nanmax) 

    return df_for_plot.sort_index()


file_full_paths = (os.path.join(dirname, (_FILE_BASE +
                                os.path.extsep +
                                 _FILE_EXTENSION))
                   for dirname in _DIRECTORIES)

df = pd.concat((pd.read_csv(in_file, sep="\t")
                for in_file in file_full_paths))

# basic data check
assert pd.api.types.is_numeric_dtype(df['diff-logAdjP-value'])

try:
    assert pd.api.types.is_numeric_dtype(df['cluster'])
except Exception:
    df.cluster.unique()

# N.B. need regex enabled to allow partial match to column
df.TF.replace(_TF_SYNONYM_MAP, inplace=True, regex=True)

# filter to only analyze the specified TFs
df = df.loc[df["TF"].isin(_TFs)]

unmod_col, mod_col = zip(*df['hypothesis'].map(lambda x: x.split('-')))

df.drop("hypothesis", axis=1, inplace=True)

# N.B. treating duplicate hypotheses equally...
unmod_col = [re.sub("\.\d+", "", unmod_h) for unmod_h in unmod_col]

df.insert(0, "unmod_H", unmod_col)

df.insert(1, "mod_H", mod_col)

# N.B. needed to prevent multiple rows having the same index,
#      resulting in downstream errors or irregularities.
df.reset_index(drop=True, inplace=True)

sort_cols = ['diff-logAdjP-value']

# remove duplicates, retaining only max
key_cols = ['TF', 'mod_H'] + sort_cols
df = df.drop_duplicates(subset=key_cols, keep='first').sort_values(sort_cols)

# N.B. There exist valid tested hypotheses (i.e., within full_centrimo_df) that
#      will not match any in our filtered set for analysis (i.e., df).
#      These were simply filtered-out upstream and need not be considered here.
#      The above operations already correctly exclude these.

assert df.drop('Type', axis=1).notnull().values.all()

vmin_df = df['diff-logAdjP-value'].min()
vmax_df = df['diff-logAdjP-value'].max()

# consider scores more extreme than this to simply be at the max value
# and also ensure the scale is always set to encompass at least this value
_TRUNCATE_VAL = 4000

vmin = max(-_TRUNCATE_VAL, vmin_df)

# always set to this, since truncated to it and cannot be less than it,
# to preserve strength of positive values
vmax = _TRUNCATE_VAL
#vmax = min(_TRUNCATE_VAL, vmax_df)

_PLOT_LINSCALE_LOG_PARAM = 1

_PLOT_LIN_COLOURSCALE = 20

# ensure these positions, x, s.t. 0 < x < 100, get marked (both +/-)
_PLOT_ADDTL_TICKS_TO_ADD = np.array([10, _PLOT_LIN_COLOURSCALE])

# Light centre, with a very large portion being the midpoint colour
# (rep. no substnative preference).
# This portion is initialized to the same amount as the non-log scale portion
# of the colour bar, and then scaled (asymmetrically), by the range of values.

mid_range_ambig_colour = _PLOT_LIN_COLOURSCALE

cmap = sns.diverging_palette(220, 10, center = "light",
                             sep = mid_range_ambig_colour, as_cmap=True)

# for now, just extend the ambig. region and lighten colours near it
mid_range_ambig_colour *= 4

# create ticks >= 0, every 100, but including the additional small marks
pos_int_ticks = np.arange(0, vmax, 100)[:-1]
pos_int_ticks = np.insert(pos_int_ticks, 1, _PLOT_ADDTL_TICKS_TO_ADD)

neg_int_ticks = np.arange(np.around(vmin, -2), 0, 100)

def _del_ticks_near_extrema(ticks_extreme_last, extreme_val):
    """Prevent a bunching of labels near the extrema."""

    if (extreme_val == _TRUNCATE_VAL):
        # no need to truncate labels in this case
        return ticks_extreme_last

    extreme_val = abs(extreme_val)

    indices_to_del = [0]

    if (extreme_val < 1200 and
        (abs(ticks_extreme_last[0]) - extreme_val) < 300):
        # delete multiple ticks, if at a large magnitude (log space)

        indices_to_del += [1]

        if extreme_val < 2000:
            indices_to_del += [2, 3]

    if (abs(ticks_extreme_last[0]) - extreme_val) < 200:
        ticks_extreme_last = np.delete(ticks_extreme_last, indices_to_del)

    return ticks_extreme_last


neg_int_ticks = _del_ticks_near_extrema(neg_int_ticks, vmin)
pos_int_ticks = _del_ticks_near_extrema(pos_int_ticks[::-1], vmax)[::-1]

# including the additional small marks, for those < 0
neg_int_ticks = np.append(neg_int_ticks,
                          np.negative(_PLOT_ADDTL_TICKS_TO_ADD[::-1]))

# log scale: prevent bunching; sparsify ticks for large numbers
ticks_to_add = np.array([])

for tick in np.concatenate((neg_int_ticks, pos_int_ticks)):
    tick = int(tick)

    # successively prune more ticks as their magnitude increases
    if (abs(tick) > 100 and tick % 3 != 0 or
        abs(tick) > 2000 and tick % 8 != 0):
            continue

    ticks_to_add = np.append(ticks_to_add, tick)

# ensure the colourbar has the ticks we want:
#     0, +/-<additional intervals from 0 to 100>, every 100, and the extrema.
cbar_ticks = np.concatenate([[round(vmin, 0)], ticks_to_add, [round(vmax, 0)]])

if (cbar_ticks[1] < cbar_ticks[0]):
    cbar_ticks = np.delete(cbar_ticks, [1])

df[TM_AREA_COL_NAME] = \
    (_AREA_SCALE_FACTOR * (_COEFF_SCALE_FACTOR*(abs(df['diff-logAdjP-value']))
     - (vmin_df / _COEFF_SCALE_FACTOR)) / ((vmax_df - vmin_df) + 1))

_LABEL_ANNOT_COL = 'cluster'

tidy_df_group_cols = ['TF', 'cluster', 'root motif', 'unmod_H', 'mod_H']

# apply the same aggregation we apply in the heatmap
tidy_df = (df.groupby(tidy_df_group_cols, as_index=False)
           ['diff-logAdjP-value', TM_AREA_COL_NAME].agg(_AGG_FUNC).
           # remove duplicates, resulting from the same
           # motif/score/cluster across different datasets (e.g. replicates)
           drop_duplicates().
           sort_values('unmod_H'))  # group by the unmodified hypothesis

df_for_plot_r = pandas2ri.py2ri(tidy_df)

# order the TFs the way we defined them
df_for_plot_r[df_for_plot_r.names.index("TF")] = \
    ro.FactorVector(df_for_plot_r[df_for_plot_r.names.index("TF")],
                    levels=ro.StrVector(_TFs), ordered=True)

# replace the cluster column with the root motif
# in a second DataFrame, used for the subgroup text labels
# (i.e. label clusters with their root motif, instead of cluster number)
df_labels_roots_only = tidy_df.copy()
df_labels_roots_only.drop(labels=[_LABEL_ANNOT_COL, 'root motif'],
                          axis=1, inplace=True)
df_labels_roots_only[_LABEL_ANNOT_COL] = tidy_df['root motif']
df_labels_roots_only = pandas2ri.py2ri(df_labels_roots_only)

# Define a colourmap that we can export from matplotlib into ggplot2
# This allows the same colours to be used and, with the same transformation
# occuring on the Python end.

# Re-centre below 0, to skew toward highlighting positive results
vcenter = int(-1 * round(mid_range_ambig_colour, -1) / 8)

# Define our norm, to re-centre, and be log transformed,
# for both positive and negative values
norm = OffsetSymLogNorm(linthresh=mid_range_ambig_colour,
                        linscale=_PLOT_LINSCALE_LOG_PARAM,
                        vmin=vmin, vcenter=vcenter, vmax=vmax)

# re-centre, per our above norm
cmap_export = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)

# define the colour for every integer from the min to the max
c_vals = np.arange(vmin, vmax, 1)
colours_rgb = cmap_export.to_rgba(c_vals)

c_vals_R = ro.FloatVector(c_vals)

# convert to a usable R Vector of R (hex) colours
colours_R = ro.StrVector([grDevices.rgb(r, g, b, a)[0] for
                          r, g, b, a in (ro.FloatVector(colour_rgb)
                                         for colour_rgb in colours_rgb)])

pal_labs = [int(round(vmin)), -1000, -100, -60, -40, -30, -20, -15, -10, 0, 5,
            10, 15, 20, 30, 40, 60, 100, 1000, int(round(vmax))]

pal = cmap_export.to_rgba(pal_labs)

# Adapted from sns.palplot
num_pal_col = len(pal) - 1

f, pal_ax = plt.subplots(1, 1, figsize=(num_pal_col, 1.2))

pal_ax.imshow(np.arange(num_pal_col).reshape(1, num_pal_col),
          cmap=mpl.colors.ListedColormap(list(pal)),
          interpolation="nearest", aspect="auto")

pal_ax.tick_params(labelsize=_TICK_LABEL_SIZE) 

pal_ax.set_xticks(np.arange(len(pal)) - .5)
pal_ax.set_yticks([-.5, .5])

pal_labs_disp = [str(pal_lab).replace("-", u"−") for pal_lab in pal_labs]
pal_labs_disp[0] = u"≤" + pal_labs_disp[0]
pal_labs_disp[-1] = u"≥" + pal_labs_disp[-1]

# easiest way to simply move first/last labels, to prevent overlap
pal_labs_disp[0] = pal_labs_disp[0] + "        "
pal_labs_disp[-1] =  "  " + pal_labs_disp[-1]

pal_ax.set_xticklabels(pal_labs_disp)
pal_ax.set_yticklabels([])

plt.tight_layout()

for ext in ("png", "pdf"):
    plt.savefig("colour_scale.{}".format(ext))

cbar_ticks_zero_idx = np.where(cbar_ticks == 0)[0][0]
cbar_ticks_to_rm = [cbar_ticks_zero_idx - i for i in range(-3, 4) if i != 0]

# remove the six values nearest 0 (no space for them on this colourbar)
cbar_ticks_adj = np.rint(np.delete(cbar_ticks, cbar_ticks_to_rm)).astype(int)

cbar_ticks_breaks = cbar_ticks_adj.copy()

cbar_ticks_R = ro.IntVector(cbar_ticks_adj)

cbar_ticks_breaks_R = ro.IntVector(cbar_ticks_breaks)

cbar_ticks_R_labs = ro.StrVector(cbar_ticks_breaks)

if cbar_ticks[0] == vmin:
    cbar_ticks_R_labs[0] = u"≤ " + str(vmin)

if cbar_ticks[-1] == vmax:
    cbar_ticks_R_labs[-1] = u"≥ " + str(vmax)

_ENABLE_SUBGROUP_LABELS = False


gp = ggplot2.ggplot(df_for_plot_r)

pp = gp + \
     ggplot2.aes_string(area=TM_AREA_COL_NAME, fill='diff.logAdjP.value',
                        label='mod_H', subgroup='cluster') + \
     treemapify.geom_treemap() + \
     treemapify.geom_treemap_subgroup_border(colour="grey35", size=ro.r.rel(4)) + \
     treemapify.geom_treemap_subgroup_text(data=df_labels_roots_only, place='centre', grow=True,
                                           alpha=0.33 if _ENABLE_SUBGROUP_LABELS else 0,
                                           colour='black', fontface='italic') + \
     treemapify.geom_treemap_text(colour='white', place='topleft', reflow=True, grow=True) + \
     ggplot2.scale_fill_gradientn(colours=colours_R, values=scales.rescale(c_vals_R),
                                  name=_PLOT_VAL_TITLE, breaks=cbar_ticks_breaks_R, labels=cbar_ticks_R_labs,
                                  # N.B. below, limit must be exactly min/max for correct colour mapping
                                  # all values above or below min/max are set to the min/max ("squish")
                                  limits=ro.IntVector([vmin - 1, vmax]), oob=scales.squish) + \
     ggplot2.facet_wrap('TF', ncol=3) + \
     ggplot2.theme(legend_position='bottom',
                   legend_key_width=ro.lib.grid.unit(32, 'lines'),
                   strip_text_x=ggplot2.element_text(size=32, colour="blue",
                                                     family=MAIN_FONT),
                   legend_text=ggplot2.element_text(size=ro.r.rel(2),
                                                    family=MAIN_FONT),
                   legend_title=ggplot2.element_text(size=ro.r.rel(2.5),
                                                     family=MAIN_FONT))

_LR_OUT_NAME = 'treemap.png'
_HR_OUT_NAME = 'treemap.pdf'

_LR_DPI = 300  # raster
_HR_DPI = 600  # unused; vector output

LR_DEV = 'png'
HR_DEV = 'pdf'

_OUT_W = 130
_OUT_H = 90

for out_name, out_dpi, out_dev in \
    ((_LR_OUT_NAME, _LR_DPI, LR_DEV),
     (_HR_OUT_NAME, _HR_DPI, HR_DEV)):
    pp.save("{}".format(out_name), device=out_dev,
            width=_OUT_W, height=_OUT_H, units="cm",
            dpi=out_dpi, limitsize=False)
