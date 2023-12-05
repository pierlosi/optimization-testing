"""Result plots for ORC case study.

Based on:

@InCollection{huster2019impact,
  author    = {Wolfgang R. Huster and Artur M. Schweidtmann and Alexander
               Mitsos},
  booktitle = {Computer Aided Chemical Engineering},
  title     = {Impact of accurate working fluid properties on the globally
               optimal design of an organic rankine cycle},
  doi       = {10.1016/b978-0-12-818597-1.50068-0},
  pages     = {427--432},
  publisher = {Elsevier},
  volume    = {47},
  year      = {2019},
}
"""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import pickle

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import numpy as np
import pandas as pd

import comando
comando.set_backend('symengine')

from ANN import make_ann

plt.rcParams.update({
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    'hatch.linewidth': 0.5   # reduce hatch linewidth
    })


# Dummy symbols as placeholders
p = comando.Symbol('p')
h = comando.Symbol('h')
s = comando.Symbol('s')
pressures = np.linspace(1.5e5, 36.29e5, 50)

s_liq_func = comando.lambdify([p, h], [make_ann('s_liq')(p, h)])
s_sat_liq_func = comando.lambdify([p], [make_ann('s__sat_liq')(p)])
s_sat_vap_func = comando.lambdify([p], [make_ann('s__sat_vap')(p)])
s_vap_func = comando.lambdify([p, h], [make_ann('s_vap')(p, h)])

T_liq_func = comando.lambdify([p, h], [make_ann('T_liq')(p, h)])
T_sat_func = comando.lambdify([p], [make_ann('T__sat')(p)])
T_vap_func = comando.lambdify([p, h], [make_ann('T_vap')(p, h)])
T_sat = T_sat_func(pressures)
# T3 = T_sat(p2)

h_liq_func = comando.lambdify([p, s], [make_ann('h_liq')(p, s)])
h_sat_liq_func = comando.lambdify([p], [make_ann('h__sat_liq')(p)])
h_sat_vap_func = comando.lambdify([p], [make_ann('h__sat_vap')(p)])

# Two-phase region
h_sat_liq = h_sat_liq_func(pressures)
h_sat_vap = h_sat_vap_func(pressures)


def get_data(res):
    """Get working fluid state data based on the given results."""
    p1 = res['p1']
    p2 = res['p2']
    h2 = res['h2']
    h2r = res['h2r']
    h3 = res['h3']
    h4 = res['h4']
    h6 = res['h6']
    h6r = res['h6r']
    h_pinch = h_sat_vap_func(p1)
    data = {'1': {'P': p1, 'T': res['T1'], 'H': res['h1'], 'S': res['s1']},
            '2': {'P': p2, 'T': res['T2'], 'H': h2, 'S': s_liq_func(p2, h2)},
            '2r': {'P': p2, 'T': res['T2r'], 'H': h2r,
                   'S': s_liq_func(p2, h2r)},
            # '2s': {'P': p2, 'T': T2s},
            '3': {'P': p2, 'T': res['T3'], 'H': h3, 'S': s_sat_liq_func(p2)},
            '4': {'P': p2, 'T': res['T4'], 'H': h4, 'S': s_sat_vap_func(p2)},
            '5': {'P': p2, 'T': res['T5'], 'H': res['h5'], 'S': res['s5']},
            '6': {'P': p1, 'T': res['T6'], 'H': h6, 'S': s_vap_func(p1, h6)},
            # '6s': {'P': p1, 'T': T6s},
            '6r': {'P': p1, 'T': res['T6r'], 'H': h6r,
                   'S': s_vap_func(p1, h6r)},
            'pinch': {'P': p1, 'T': T_sat_func(p1), 'H': h_pinch,
                      'S': s_vap_func(p1, h_pinch)},
            }
    return pd.DataFrame(data)


def plot_ph(ax, res=None, **kwargs):
    """Plot a pressure enthalpy diagram."""
    if res is None:

        ax.plot(np.hstack([h_sat_liq, np.flip(h_sat_vap)])/1e3,
                np.hstack([pressures, np.flip(pressures)])/1e3,
                c='k', **kwargs)
        return

    data = get_data(res)
    cycle = data.copy()
    cycle['1_'] = data['1']
    ax.plot(cycle.loc['H']/1e3, cycle.loc['P']/1e3, c='k',  **kwargs)
    return data.loc[['H', 'P']]/1e3


def plot_Th(ax, res=None, **kwargs):
    """Plot a temperature enthalpy diagram."""
    if res is None:
        p_bound_enths = np.linspace(h_min, h_max)

        def get_temps(p):
            h_lim1 = h_sat_liq_func(p)
            h_lim2 = h_sat_vap_func(p)
            T_sat = T_sat_func(p)
            return np.array([T_liq_func(p, h) if h <= h_lim1 else T_sat
                             if h <= h_lim2 else T_vap_func(p, h)
                             for h in p_bound_enths])
        p1_bound_temps = [get_temps(5e5), get_temps(2e5)]
        p2_bound_temps = [get_temps(35e5), get_temps(6e5)]

        ax.fill_between(p_bound_enths/1000, *p1_bound_temps, hatch='\\' * 5,
                        facecolor="none", edgecolor='k', alpha=0.4, lw=0.5,
                        label='$p_1$ region')
        ax.fill_between(p_bound_enths/1000, *p2_bound_temps, hatch='/' * 5,
                        facecolor="none", edgecolor='k', alpha=0.4, lw=0.5,
                        label='$p_2$ region')

        saturation_enths = np.array([*h_sat_liq, *np.flip(h_sat_vap)])/1e3
        saturation_temps = np.array([*T_sat, *np.flip(T_sat)])
        ax.plot(saturation_enths, saturation_temps, c='k', lw=1, **kwargs)
        return
    data = get_data(res)
    p1 = res['p1']
    p2 = res['p2']

    def isobar(state1, state2, func, n):
        enths = np.linspace(*data.loc['H'][[state1, state2]], n)
        temps = np.array([func(h) for h in enths])
        return enths, temps

    # enths23, temps23 = isobar('2', '3', lambda x: T_liq_func(p2, x), 10)
    enths22r, temps22r = isobar('2', '2r', lambda x: T_liq_func(p2, x), 4)
    enths2r3, temps2r3 = isobar('2r', '3', lambda x: T_liq_func(p2, x), 4)
    enths45, temps45 = isobar('4', '5', lambda x: T_vap_func(p2, x), 3)
    enths66r, temps66r = isobar('6', '6r', lambda x: T_vap_func(p1, x), 5)

    process_enths = np.array([data.loc['H']['1'], *enths22r, *enths2r3,
                              *enths45, *enths66r,
                              *data.loc['H'][['pinch', '1']]])/1e3
    process_temps = np.array([data.loc['T']['1'], *temps22r, *temps2r3,
                              *temps45, *temps66r,
                              *data.loc['T'][['pinch', '1']]])

    # positions of process states 1, 2, 2r, 3, 4, 5, 6, 6r and pinch
    pro = [0, 1, 5, 8, 9, 11, 12, 16, 17]
    kwargs['zorder'] = 1
    ax.plot(process_enths, process_temps, c='k', markevery=pro, **kwargs)
    kwargs['label'] = None
    kwargs['ls'] = ''
    kwargs['zorder'] = 2
    ax.plot(process_enths, process_temps, c='k', markevery=pro, **kwargs)

    # ax.scatter(process_enths[pro], process_temps[pro], **kwargs)
    return data.loc[['H', 'T']]/np.array([[1e3, 1]]).T


# plot boundaries
h_min = 5e3
h_max = 550e3
p_min = 2e5
p_max = 38e5
T_min = 270
T_max = 420


with open('results.pickle', 'rb') as f:
    results_baron, results_pyomo, results_maingo = pickle.load(f)

# TH plot
INCH_PER_POINT = 1 / 72.27
w = 275 * INCH_PER_POINT  # Note Latex reports 252 points but image is smaller
h = 1.5 * 0.618 * w
fig, ax = plt.subplots(figsize=(w, h))
plot_Th(ax, ls='dotted', label='saturation curve')
states_baron = plot_Th(ax, results_baron, label='process (BARON)', marker='^',
                       markerfacecolor='#eb5f73', markeredgecolor='#eb5f73',
                       ls='dashed', lw=1)
states_maingo = plot_Th(ax, results_maingo, label='process (MAiNGO)',
                        markerfacecolor='#023d6b', markeredgecolor='#023d6b',
                        marker='v', lw=1)
rx = 40  # lable distance in x
ry = 20  # lable difference in y
label_angles = [
    270,  # 1
    170,  # 2
    115,  # 2r
    120,  # 3
    145,  # 4
    -45,  # 5
    0,   # 6
    -45,  # 6r
    145  # pinch
]
states_mean = 0.5 * (states_baron + states_maingo)
# plt.scatter(res['h6s']/1e3, res['T6s'], label='6s')
for i, (index, data_point) in enumerate(states_mean.items()):
    a = np.pi/180 * (label_angles[i])
    xy_mean = data_point.values
    xytext = xy_mean + [rx * np.cos(a), ry * np.sin(a)]
    # plt.scatter(*xytext, marker='x', c='r')
    ax.plot(*np.vstack([states_baron.values[:, i], xytext]).T, c='k', lw=.5)
    ax.plot(*np.vstack([states_maingo.values[:, i], xytext]).T, c='k', lw=.5)
    ax.annotate(index, xy_mean, xytext=xytext, ha='center',
                bbox=dict(boxstyle="Round, pad=0.2", fc="w", ec="k", lw=0.5))


ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.xlabel('specific enthalpy [kJ/kg]')
plt.ylabel('temperature [K]')
plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
           mode="expand", borderaxespad=0, ncol=2, handlelength=2.88)
ax.set_ylim(T_min, T_max)
ax.set_xlim(h_min/1e3, h_max/1e3)
plt.tight_layout()
plt.savefig("TH.pgf")
