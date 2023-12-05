"""Case study for the optimal operating point of an organic Rankine cycle.

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
from contextlib import nullcontext
import sys

from time import time
import comando
comando.set_backend('symengine')


def run_ORC_case_study(tol=1e-3, verbose=True, reuse=False):
    """Run the organic Rankine cycle case study."""
    from .ORC_components import HeatExchanger, Pump, Turbine, CoolingSystem
    from .ANN import make_ann
    t1 = time()
    ES = comando.System('ORC')

    # Variables  (All in SI units)
    mdot = ES.make_operational_variable("mdot", bounds=(100, 2000),
                                        init_val=408.3575967)  # kg/s
    h2r = ES.make_operational_variable("h_2r", bounds=(30e3, 400e3),
                                       init_val=59785.3)  # J/kg
    h6s = ES.make_operational_variable("h_6s", bounds=(350e3, 700e3),
                                       init_val=350e3)  # J/kg
    p1 = ES.make_operational_variable("p_1", bounds=(2e5, 5e5),
                                      init_val=2.439e5)  # Pa = kg/[m*s^2]
    p2 = ES.make_operational_variable("p_2", bounds=(6e5, 35e5),
                                      init_val=15.950538e5)  # Pa = kg/[m*s^2]

    T_cw_in = ES.make_parameter('T_cw_in', 288)  # K
    T_br_in = ES.make_parameter('T_br_in', 408)  # K
    T_br_out = ES.make_parameter('T_br_out', 357)  # K
    mdot_cp_br = ES.make_parameter('mdot_cp_br', 3627e3)  # W/K
    dTmin_cond = ES.make_parameter('dTmin_cond', 10)  # K
    dTmin_rec = ES.make_parameter('dTmin_rec', 15)  # K
    dTmin_eva = ES.make_parameter('dTmin_eva', 15)  # K
    eta_p = ES.make_parameter('eta_p', 0.8)
    etat = ES.make_parameter('etat', 0.9)
    eps = ES.make_parameter('eps', 1e-5)

    Qdot_in = ES.add_expression('Qdot_in', mdot_cp_br * (T_br_in - T_br_out))

    s_sat_liq = make_ann('s__sat_liq')
    s1 = s_sat_liq(p1)

    T_sat = make_ann('T__sat')  # depth : 13
    T1 = T_sat(p1)
    T3 = T_sat(p2)

    h_sat_liq = make_ann('h__sat_liq')  # depth : 13
    h1 = h_sat_liq(p1)
    h3 = h_sat_liq(p2)

    h_sat_vap = make_ann('h__sat_vap')  # depth : 13
    h7 = h_sat_vap(p1)
    h4 = h_sat_vap(p2)

    h_liq = make_ann('h_liq')  # depth : 25
    h2s = h_liq(p2, s1)

    T_liq = make_ann('T_liq')  # depth : 13
    T2r = T_liq(p2, h2r)

    h5 = h2r + Qdot_in / mdot
    T_vap = make_ann('T_vap')  # depth : 13
    T5 = T_vap(p2, h5)

    s_vap = make_ann('s_vap')  # depth : 13
    s5 = s_vap(p2, h5)

    T6s = T_vap(p1, h6s)

    s6s = s_vap(p1, h6s)

    pump = Pump('pump', mdot, p1, p2, h1, h2s, eta_p)
    h2 = pump['h_out']

    T2 = T_liq(p2, h2)

    economizer = HeatExchanger('eco', dTmin_eva,
                               T_h_out=T_br_out, T_c_in=T2r, T_c_out=T3,
                               h_c_in=h2r, h_c_out=h3, mdot_c=mdot,
                               mdot_cp_h=mdot_cp_br)
    T_br_3 = economizer['T_h_in']

    T4 = T3  # isothermal evaporation

    superheater = HeatExchanger('sup', dTmin_eva, T_h_in=T_br_in,
                                T_c_in=T4, T_c_out=T5, h_c_in=h4, h_c_out=h5,
                                mdot_c=mdot, mdot_cp_h=mdot_cp_br)
    T_br_4 = superheater['T_h_out']

    evaporator = HeatExchanger('eva', dTmin_eva, T_h_in=T_br_4,
                               T_h_out=T_br_3, T_c_in=T3, T_c_out=T4,
                               h_c_in=h3, h_c_out=h4, mdot_c=mdot,
                               mdot_cp_h=mdot_cp_br)
    # NOTE: By fully specifying both sides of the evaporator, we introduce the
    #       additional constraint that the heat flows computed based on the
    #       corresponding quantities at both sides must be equal!
    T_br_3 = evaporator['T_h_out']

    turbine = Turbine('tur', mdot, p2, p1, h5, h6s, etat)
    h6 = turbine['h_out']
    T6 = T_vap(p1, h6)
    T6r = T_vap(p1, h6 - (h2r - h2))

    recuperator = HeatExchanger('rec', dTmin_rec, T_h_in=T6,
                                T_h_out=T6r, T_c_in=T2, T_c_out=T2r, h_h_in=h6,
                                h_c_in=h2, h_c_out=h2r, mdot_h=mdot,
                                mdot_c=mdot)
    h6r = recuperator['h_h_out']  # == h6 - (h2r - h2)

    # pinch inside of the condenser, right after saturation of the vapor phase
    T_cw_7 = T1 - dTmin_cond
    dT_cw_7_calc = T_cw_7 - T_cw_in
    dT_cw_7 = comando.Max(eps, dT_cw_7_calc)

    Q_17 = mdot * (h7 - h1)
    mdot_cp_cw = Q_17/(dT_cw_7)

    condenser = HeatExchanger('cond', dTmin_cond, T_h_in=T6r,
                              T_h_out=T1, T_c_in=T_cw_in, h_h_in=h6r,
                              h_h_out=h1, mdot_h=mdot, mdot_cp_c=mdot_cp_cw)
    condenser.add_ge_constraint(dT_cw_7_calc, 0, 'dT_cw_7 >= 0')
    T_cw_out = condenser['T_c_out']

    cooling_sys = CoolingSystem('cooling_sys', mdot_cp_cw, T_cw_out, T_cw_in)

    # NOTE: The constraints generated from the connections are all recognized
    #       to be redundant and thus eliminated during problem generation,
    #       because of the sequential way in which we defined the components!
    connections = {'1': (condenser.Hdot_h_out, pump.Hdot_in),
                   '2': (pump.Hdot_out, recuperator.Hdot_c_in),
                   '2r': (recuperator.Hdot_c_out, economizer.Hdot_c_in),
                   '3': (economizer.Hdot_c_out, evaporator.Hdot_c_in),
                   '4': (evaporator.Hdot_c_out, superheater.Hdot_c_in),
                   '5': (superheater.Hdot_c_out, turbine.Hdot_in),
                   '6': (turbine.Hdot_out, recuperator.Hdot_h_in),
                   '6r': (recuperator.Hdot_h_out, condenser.Hdot_h_in),
                   'br_4': (superheater.Hdot_h_out, evaporator.Hdot_h_in),
                   'br_3': (evaporator.Hdot_h_out, economizer.Hdot_h_in),
                   'cw_in': (cooling_sys.Hdot_out, condenser.Hdot_c_in),
                   'cw_out': (condenser.Hdot_c_out, cooling_sys.Hdot_in)
                   }
    components = [pump, recuperator, economizer, evaporator, superheater,
                  turbine, condenser, cooling_sys]
    for comp in components:
        ES.add(comp)
    for bus_id, bus in connections.items():
        ES.connect(bus_id, bus)

    P_net = turbine['P'] - pump['P'] - cooling_sys['P']
    ES.add_ge_constraint(P_net * 1e-7, 0, 'P_net >= 0')
    ES.add_expression('efficiency', P_net / Qdot_in)

    ES.add_eq_constraint((s5 - s6s) * 1e-3, 0, 's5 == s6s')
    ES.add_le_constraint((T1 - T6s) * 1e-2, 0, 'T1 <= T6s')

    # Maximize net power production, considering the nominal operating point
    P = ES.create_problem(0, -P_net * 1e-7, timesteps=(['nominal'], 1))
    print(f'Problem has {P.num_vars} variables and {P.num_cons} constraints.')
    print(f'\nCreating ORC model took: {time() - t1} s')

    def get_results():
        from comando.utility import evaluate
        return {'mdot': evaluate(mdot, idx='nominal'),
                'p1': evaluate(p1, idx='nominal'),
                'p2': evaluate(p2, idx='nominal'),
                'T1': evaluate(T1, idx='nominal'),
                'T2': evaluate(T2, idx='nominal'),
                'T2r': evaluate(T2r, idx='nominal'),
                'T3': evaluate(T3, idx='nominal'),
                'T4': evaluate(T4, idx='nominal'),
                'T5': evaluate(T5, idx='nominal'),
                'T6s': evaluate(T6s, idx='nominal'),
                'T6': evaluate(T6, idx='nominal'),
                'T6r': evaluate(T6r, idx='nominal'),
                'h1': evaluate(h1, idx='nominal'),
                'h2s': evaluate(h2s, idx='nominal'),
                'h2': evaluate(h2, idx='nominal'),
                'h2r': evaluate(h2r, idx='nominal'),
                'h3': evaluate(h3, idx='nominal'),
                'h4': evaluate(h4, idx='nominal'),
                'h5': evaluate(h5, idx='nominal'),
                'h6s': evaluate(h6s, idx='nominal'),
                'h6': evaluate(h6, idx='nominal'),
                'h6r': evaluate(h6r, idx='nominal'),
                's1': evaluate(s1, idx='nominal'),
                's5': evaluate(s5, idx='nominal'),
                's6s': evaluate(s6s, idx='nominal'),
                'P_P': evaluate(pump['P'], idx='nominal'),
                'P_T': evaluate(turbine['P'], idx='nominal'),
                'P_F': evaluate(cooling_sys['P'], idx='nominal'),
                'P_net': evaluate(P_net, idx='nominal')}

    ###########################################################################
    # Solve with BARON...
    ###########################################################################

    end = '\n' if verbose else ''

    baron_opts = dict(epsr=tol, epsa=tol)
    # ...via our own interface
    print('Writing custom .bar file and solving...', end=end, flush=True)
    custom_t0 = time()

    from comando.interfaces.baron import solve, baron_str_map

    # Teach baron how to handle the Max function (the tanh function is handled
    # as lambda arg: f'(1 - 2/(exp(2 * ({arg})) + 1))'} by our interface)
    baron_str_map['Max'] = lambda a, b: \
        f'(0.5 * ({a} + {b} + (({a} - {b} + 1e-04)^2)^0.5))'
    solve(P, 'custom.bar', silent=~verbose, reuse=reuse, **baron_opts)
    print(f' took {time() - custom_t0} s')
    results_baron = get_results()

    # ...indirectly via Pyomo (to make sure our interface isn't flawed)
    print('Translating to pyomo and solving...', end=end, flush=True)
    pyomo_t0 = time()
    import pyomo
    from comando.interfaces.pyomo import to_pyomo, pyomo_op_map
    # Teach pyomo how to handle the Max function
    pyomo_op_map['Max'] = comando.utility.smooth_max
    # Replace tanh(arg) with 1 - 2/(exp(2 * arg) + 1) so baron can handle it
    pyomo_op_map['tanh'] = lambda arg: 1 - 2/(pyomo.core.exp(2 * arg) + 1)

    m = to_pyomo(P)

    from comando.utility import silence, get_latest
    log = 'pyomo.baron.log'

    context = nullcontext if verbose else silence
    if reuse:
        from pyomo.solvers.plugins.solvers import BARON
        shell = BARON.BARONSHELL()
        with context():
            res = shell.solve('pyomo.bar', type='baron', tee=True)
    else:
        import tempfile
        with context():
            res = m.solve('baron', tee=True, keepfiles=True,
                          options=baron_opts, logfile=log)
        # Retrieve the last .bar file from the system's TMP folder
        newest_bar_file = get_latest(tempfile.gettempdir() + '/*.bar')
        from shutil import move
        import pathlib
        move(newest_bar_file, pathlib.Path('.').resolve() / 'pyomo.bar')
    print(f' took {time() - pyomo_t0} s')

    results_pyomo = get_results()

    ###########################################################################
    # Solve with MAiNGO
    ###########################################################################
    print('Solving via MAiNGO...', end=end, flush=True)
    maingo_t0 = time()
    from comando.interfaces.maingo_api import MaingoProblem

    mp = MaingoProblem(P)
    sol, res = mp.solve(epsilonR=tol, epsilonA=tol)
    print(f' took {time() - maingo_t0} s')
    results_maingo = get_results()

    results = results_baron, results_pyomo, results_maingo
    import pickle
    with open('results.pickle', 'wb') as f:
        pickle.dump(results, f)

    return results


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser()
    # Add the arguments to the parser
    ap.add_argument("-tol", "-t", required=False, default=1e-3,
                    type=float, help="absolute/relative optimality tolerance")
    ap.add_argument("-reuse", "-r", required=False, default=False,
                    action='store_true', help="reuse custom BARON input file")
    ap.add_argument("-verbose", "-v", required=False, default=False,
                    action='store_true', help="give verbose output")

    if __package__ is None:
        from pathlib import Path
        file = Path(__file__).resolve()
        sys.path.append(str(file.parents[1]))
    import ORC
    __package__ = ORC.__name__
    results = run_ORC_case_study(**vars(ap.parse_args()))
    interfaces = ['BARON', 'Pyomo-BARON', 'MAiNGO-API']
    from pandas import DataFrame
    print(DataFrame(results, interfaces).T)
