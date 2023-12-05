"""Plot the results from the IES case study."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu, David Shu
import pickle

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{eurosym}')

# define colors
colors = {'blue': (2/255, 61/255, 107/255),
          'lightblue': (173/255, 189/255, 227/255),
          'green': (185/255, 210/255, 95/255),
          'yellow': (250/255, 235/255, 90/255),
          'violet': (175/255, 130/255, 185/255),
          'red': (235/255, 95/255, 115/255),
          'orange': (250/255, 180/255, 90/255)}
# allocate colors to component types
colormap = {'AC': colors['lightblue'],
            'B': colors['red'],
            'CC': colors['lightblue'],
            'CHP': colors['orange'],
            'HP': colors['violet'],
            'PV': colors['yellow'],
            'BAT': colors['yellow'],
            r'TES$_\mathrm{c}$': colors['lightblue'],
            r'TES$_\mathrm{h}$': colors['red']}


def get_data(dvs, ovs):
    """Get capacity data for each iteration."""
    all_caps = dvs['value'].unstack(level=1)
    cap_names = [name for name in all_caps.columns if name.endswith('nom')]
    caps = dvs['value'].unstack(level=1)[cap_names].rename(columns=lambda n:
                                                           n[:5].split('_')[0])
    # Cumulate CHPs and PV
    cum_caps = caps.groupby(caps.columns, axis=1).sum()

    # Split conversion and storage units
    conv = ['AC', 'B', 'CC', 'CHP', 'HP', 'PV']
    sto = ['BAT', 'STC', 'STH']
    mean_conv_caps = cum_caps[conv].mean().sort_values(ascending=False)
    conv_order = mean_conv_caps[mean_conv_caps > 0].index
    mean_sto_caps = cum_caps[sto].mean().sort_values(ascending=False)
    sto_order = mean_sto_caps[mean_sto_caps > 0].index

    return (cum_caps[conv_order],
            cum_caps[sto_order].rename(columns={'STC': r'TES$_\mathrm{c}$',
                                                'STH': r'TES$_\mathrm{h}$'}))


def stacked_plot(pos, heights, width, offset=0, ax=None, legend=True):
    """Create a stacked bar plot at the given positions."""
    if ax is None:
        ax = plt.subplot()

    bottom = pd.Series(0, heights.index)
    for col in heights.columns:
        height = heights[col].values
        ax.bar(pos + offset, height, width, bottom=bottom,
               color=colormap[col], label=col if legend else None,
               edgecolor='k',
               alpha=.99, linewidth=0.25)
        # last two options: https://stackoverflow.com/a/59389823/7416115
        bottom += height


# Lin data
from comando.utility import get_latest
filename = get_latest('*_results.pickle')
with open(filename, 'rb') as f:
    obj_vals, dvs, ovs = pickle.load(f)

conv_data, sto_data = get_data(dvs, ovs)

gwi_data = obj_vals.gwi/1000
tac_data = obj_vals.tac/1000

# NL Data
filename_nl = get_latest('*_results_nl.pickle')
with open(filename_nl, 'rb') as f:
    nl_obj_vals, nl_dvs, nl_ovs = pickle.load(f)

nl_conv_data, nl_sto_data = get_data(nl_dvs, nl_ovs)

nl_gwi_data = nl_obj_vals.gwi/1000
nl_tac_data = nl_obj_vals.tac/1000

# PLOT Nonlinear results
NL = True

# plot
fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(3.5, 4),
                        gridspec_kw={'hspace': 0, 'height_ratios': [2, 1.2]})

if NL:
    # width = (gwi_data.max() - gwi_data.min()) / len(obj_vals) / 4
    width = (gwi_data.max() - gwi_data.min()) / len(obj_vals) / 6
else:
    width = (gwi_data.max() - gwi_data.min()) / len(obj_vals) / 2.5

# Conversion data
stacked_plot(gwi_data, conv_data, width, -width * 1/2, axs[0])
# NOTE: Pandas seems incapable to specify bar position (just order!)
# conv_data['gwi'] = gwi_data
# sto_data['gwi'] = gwi_data
# conv_data.plot(x='gwi', position=1.1, width=0.3, color=colormap, kind='bar',
#                stacked=True, ax=axs[0])
if NL:
    stacked_plot(gwi_data, nl_conv_data, width, -width * 7/4, axs[0],
                 legend=False)

# Storage data
ax2 = axs[0].twinx()
stacked_plot(gwi_data, sto_data, width, width * 1/2, ax2)

if NL:
    stacked_plot(gwi_data, nl_sto_data, width, width * 7/4, ax2, legend=False)

# Hatching storages since they have the same colors as conversion components
bars = ax2.patches
hatches = sum(([h] * len(sto_data) for h in ['\\\\\\\\', '////', '.....']), [])
for bar, hatch in zip(bars, 2 * hatches):
    bar.set_hatch(hatch)

# Merging the two legends
leg1 = axs[0].legend()
leg2 = ax2.legend()
l_diff = len(conv_data.columns) - len(sto_data.columns)
patches = [*leg1.get_patches(), *leg2.get_patches()]
labels = [text.get_text() for text in leg1.get_texts() + leg2.get_texts()]
plt.legend(patches, labels, loc='upper right', borderpad=.5, labelspacing=0.25)
#          fancybox=True, framealpha=1, shadow=True,
leg1.remove()
axs[0].set_ylabel(r'conversion capacity [MW]')
ax2.set_ylabel('storage capacity [MWh]')
y_max = 3.9  # axs[0].get_ylim()[1] * 1.05
axs[0].set_ylim(0, y_max)

# Now we want to make visible how pareto points correspond to capacities
# OPTION 1: Vertical lines from pareto points to corresponding capacities
axs[1].set_ylim(0.25, 3.9)
axs[1].vlines(gwi_data, tac_data, axs[1].get_ylim()[1], color='k', alpha=0.3,
              linestyles='dotted', linewidth=1)
axs[1].scatter(gwi_data, tac_data, marker='x', color='k',
               label='MILP solution')
if NL:
    axs[1].scatter(nl_gwi_data, nl_tac_data, marker='+', color='k',
                   label='NLP solution')
    for gwi, tac, dgwi, dtac in zip(gwi_data, tac_data, nl_gwi_data - gwi_data,
                                    nl_tac_data - tac_data):
        axs[1].arrow(gwi, tac, dgwi, dtac, lw=0.25)

rx = 0.04
ry = 0.45
for i, pos in enumerate(zip(gwi_data, tac_data)):
    a = np.pi/180 * 45  # (180 + 90/7 * i)
    textpos = [rx * np.cos(a) + pos[0], ry * np.sin(a) + pos[1]]
    axs[1].annotate(8-i, pos, xytext=textpos, ha='center', va='center',
                    bbox=dict(boxstyle="Round, pad=0.15", fc="w", ec="k",
                              lw=.5))

from matplotlib.ticker import FormatStrFormatter
axs[0].yaxis.set_major_formatter(FormatStrFormatter('%d'))

axs[1].set_xlabel(r'GWI [kt$_{\mathrm{CO}_2}$/a]')
axs[1].set_ylabel(r'TAC [Mio. {\euro}]')
axs[1].legend(loc='upper right', framealpha=1, handlelength=0.5, )

plt.tight_layout()
plt.show()
if NL:
    fig.savefig(f"CS1_NL_{filename[:16]}_{filename_nl[:16]}.pdf")
    print(obj_vals.join(nl_obj_vals, rsuffix='_nl'))
else:
    fig.savefig(f"CS1_NL_{filename[:16]}_{filename_nl[:16]}.pdf")
