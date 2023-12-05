"""Industrial energy system (IES) case study - design benchmark."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu, David Shu


def IES_model():
    """Create a model for an industrial energy system."""
    from comando import System
    from components.example_components import Grid, Demand
    from .IES_components import CombinedHeatAndPower, Boiler, HeatPump, \
        AbsorptionChiller, CompressionChiller, PhotovoltaicModule, \
        Battery, CoolWaterStorage, HotWaterStorage
    ES = System('ES')
    T_amb = ES.make_parameter('T_amb')

    # Power components:
    POWER = Grid('Electricity', consumption_limit=100, feedin_limit=100,
                 single_connector=True, constrain_flow=True)
    ED = Demand('Electric')
    CHP = CombinedHeatAndPower('CHP')

    # Heating components:
    # constant gas price in [kEUR/MWh] and constant emissions factor in [t/MWh]
    GAS = Grid('Gas', price=0.06, co2_factor=0.244)
    HD = Demand('Heat')
    B = Boiler('B')
    HP = HeatPump('HP', T_amb)

    # Cooling components:
    CD = Demand('Cooling')
    AC = AbsorptionChiller('AC')
    CC = CompressionChiller('CC')

    # Solar components
    PV_OB = PhotovoltaicModule('PV_OB')
    PV_EF = PhotovoltaicModule('PV_EF')

    # Storage units
    STH = HotWaterStorage('STH')
    STC = CoolWaterStorage('STC')
    BAT = Battery('BAT')

    # Set cyclic initial conditions
    # TODO: We should provide a cleaner method to do this (or sth equivalent)!
    from comando import cyclic
    for comp in STH, STC, BAT:
        comp._states_dict[comp['soc']] = \
            cyclic, *comp._states_dict[comp['soc']][1:]

    components = [POWER, GAS,
                  HD, CD, ED,
                  CHP, B, HP,
                  AC, CC,
                  PV_OB, PV_EF,
                  STH, STC, BAT]
    connections = {
        'Power_Bus': [POWER.CONSUMPTION, ED.IN, HP.IN, CC.IN, BAT.IN,
                      CHP.POWER_OUT, PV_OB.OUT, PV_EF.OUT, BAT.OUT],
        'Gas_Bus': [CHP.IN, B.IN, GAS.CONSUMPTION],
        'Heat_Bus': [HD.IN, AC.IN, STH.IN,
                     CHP.HEAT_OUT, B.OUT, HP.OUT, STH.OUT],
        'Cooling_Bus': [CD.IN, STC.IN, AC.OUT, CC.OUT, STC.OUT]}

    for c in components:
        ES.add(c)
    for bus_id, connectors in connections.items():
        ES.connect(bus_id, connectors)

    # add expressions to energy system
    for i in ['investment_costs', 'fixed_costs', 'variable_costs',
              'emissions']:
        ES.add_expression(i, ES.aggregate_component_expressions(i))

    # limit PV nominal size based on available rooftop surface
    A_OB = 5261  # [m^2]
    A_EF = 205  # [m^2]
    pel_area = PV_OB.pel_area  # [MW/m^2] nominal output per m^2
    ES.add_le_constraint(PV_OB['Pel_out_nom'], A_OB * pel_area, 'PV_OB_cap')
    ES.add_le_constraint(PV_EF['Pel_out_nom'], A_EF * pel_area, 'PV_EF_cap')

    return ES


def run_IES_case_study(num_sol=8, tolerance=0.01, seed=123, verbose=False,
                       nl_correction_only=False):
    """Run the integrated energy system case study."""
    from datetime import datetime
    import pickle
    import pathlib

    from pandas import read_csv

    from comando.utility import make_tac_objective
    from . import user_defined_algorithm
    from .user_defined_algorithm import make_multiobjective, \
        solve_multiobjective

    ES = IES_model()  # Build model

    # Define Objective
    # economic and environmental parameters and expressions
    # TAC with 4 years amortization and 8% interest (design and variable cost)
    dc, vc = make_tac_objective(ES, n=4, i=0.08)
    co2_tot = ES.get_expression('emissions')

    # Load and format data
    data = read_csv(pathlib.Path(__file__).parent / 'data.csv', sep=',')
    data.set_index(['s', 't'], inplace=True)
    scenarios = data['period_weight'].groupby('s').first().rename('pi')
    timesteps = data['dt']

    # Define problem
    P_mo = make_multiobjective(ES, timesteps, scenarios, data,
                               tac=(dc, vc), gwi=(0, co2_tot))

    print(f'\nUsing relative solver tolerance {tolerance} and seed {seed}\n')
    user_defined_algorithm.TOL = tolerance
    user_defined_algorithm.SEED = seed

    if not nl_correction_only:
        # First we solve a linearization (MILP)
        obj_vals, dvs, ovs = solve_multiobjective(P_mo, num_sol,
                                                  (3, 'convex_combination'),
                                                  silent=~verbose)
        now = datetime.now().strftime("%Y_%m_%d_%H_%M")
        with open(f'{now}_results.pickle', 'wb') as f:
            pickle.dump((obj_vals, dvs, ovs), f)
    else:
        from comando.utility import get_latest
        res_file = get_latest('*_results.pickle')  # Retrieve the last results
        print(f'Using data from latest results file: {res_file}')
        with open(res_file, 'rb') as f:
            obj_vals, dvs, ovs = pickle.load(f)

    # Reset bounds
    for obj_var in P_mo.obj_vars:
        obj_var.bounds = None, None
    gwi_max = obj_vals.gwi.max()
    step = (obj_vals.gwi.max() - obj_vals.gwi.min()) / (num_sol - 1)
    ub_MILP = [gwi_max - step * it for it in range(num_sol)]
    ub_MILP[-1] = None  # In the last iteration gwi is minimized!

    # Now we define a callback that initializes all variables with values from
    # a given iteration and then fixes the binary ones, giving us an NLP...
    def init_and_fix(P, it, lb, ub):
        """Initialize variables and fix binaries with values from it."""
        print('\t\tinitializing...')
        P.design = dvs.value.loc[it]
        P.operation = ovs.loc[it]
        for dv in P.design_variables:
            if dv.is_binary:
                dv.fix()
                print(f'\t\t\t{dv}', bool(dv.value), dv.bounds)
        for ov in P.operational_variables:
            if ov.is_binary:
                for ovi in ov:
                    ovi.fix()

        print('\t\tinitial objectives:')
        for name, obj_var, expr in zip(P.obj_names, P.obj_vars, P.obj_exprs):
            obj_var.value = expr.value
            print(f'\t\t\t{name} = {obj_var.value}')

        print('\t\tmaximum constraint violation:')
        vio = P.get_constraint_violations()
        max_vio = vio.max()
        max_vio_con = vio[vio == max_vio]
        if any(max_vio_con):
            print('\t\t\t' + ': '.join(str(e) for e in
                                       next(iter(max_vio_con.items()))))
        else:
            print('\t\t\tNone!')

        ub_index = 0 if it == num_sol - 1 else 1
        ub[ub_index] = ub_MILP[it]  # set bounds to MILP compatible values

    # .. which we can again solve using the same algorithm
    nl_obj_vals, nl_dvs, nl_ovs = solve_multiobjective(P_mo, num_sol,
                                                       silent=~verbose,
                                                       callback=init_and_fix)

    now = datetime.now().strftime("%Y_%m_%d_%H_%M")
    with open(f'{now}_results_nl.pickle', 'wb') as f:
        pickle.dump((nl_obj_vals, nl_dvs, nl_ovs), f)

    print(obj_vals.join(nl_obj_vals, rsuffix='_nl'))


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser()
    # Add the arguments to the parser
    ap.add_argument("-num_sol", "-num", required=False, default=8,
                    type=int,
                    help="number of solutions in multiobjective optimization")
    ap.add_argument("-tolerance", "-tol", required=False,
                    default=0.01, type=float,
                    help="default solver tolerance")
    ap.add_argument("-seed", required=False, default=123, type=int,
                    help="random seed for MILP solver")
    ap.add_argument("-verbose", "-v", required=False, default=False,
                    action='store_true', help="give verbose output")
    ap.add_argument("-nl_correction_only", "-nl", required=False,
                    default=False, action='store_true',
                    help="whether to only do the nonlinear correction")

    if __package__ is None:
        import sys
        from pathlib import Path
        file = Path(__file__).resolve()
        sys.path.append(str(file.parents[1]))
    import IES
    __package__ = IES.__name__
    run_IES_case_study(**vars(ap.parse_args()))
