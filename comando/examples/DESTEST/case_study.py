"""Case study for low temperature district heating network design."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu, Dominik Hering
import pathlib

import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import ConvexHull
from scipy.spatial.qhull import QhullError
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

from comando.core import System
from comando.interfaces.gurobi import to_gurobi
from comando.utility import make_tac_objective


data_path = pathlib.Path(__file__).parent / 'data'


SEED = 123


def C2K(temperature):
    """Convert temperature from Celsius to Kelvin."""
    return temperature + 273.15


def cluster(data, n_clusters=11, rs=42):
    """Cluster the chosen base data and plot total_data against clusters."""
    cm = plt.get_cmap('gist_earth')
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(
        color=[cm(1. * i / n_clusters) for i in range(n_clusters)])

    tot_name = 'SimpleDistrict'
    dem_name = [column for column in data.columns
                if column.startswith('SimpleDistrict_')]
    consumption_data = data.where(data[tot_name] >= 0.1).dropna()
    tot_data = consumption_data[['T_amb', tot_name]]

    def get_extreme_point(n_clusters, base_data):
        """Return the extreme point matching the base data."""
        if base_data is consumption_data:
            extreme_point = consumption_data.max().to_frame().T
            extreme_point['T_amb'] = -12
            return extreme_point
        return pd.DataFrame([-12, consumption_data[dem_name].max().sum()],
                            ['T_amb', tot_name]).T

    scaler = StandardScaler()
    scaler.fit(tot_data)
    scaled_data = scaler.transform(tot_data)

    res = KMeans(n_clusters - 1, random_state=rs).fit(scaled_data)
    membership = res.predict(scaled_data)

    # Total demands (needed for plot)
    tot_map = {0: get_extreme_point(n_clusters, tot_data),
               **{c + 1: tot_data[membership == c]
                  for c in range(n_clusters - 1)}}
    # Consumer demands (needed for calculation)
    cluster_map = {0: get_extreme_point(n_clusters, consumption_data),
                   **{c + 1: consumption_data[membership == c]
                      for c in range(n_clusters - 1)}}

    # Plot resulting clusters
    for members in tot_map.values():
        plt.plot(members['T_amb'], members[tot_name] / 1000, '+', ms=5.5)
        center = members.mean().values
        try:
            hull = ConvexHull(members)
            for simplex in hull.simplices:
                plt.plot(members.values[simplex, 0],
                         members.values[simplex, 1] / 1000, 'k-')
        except QhullError:
            pass
        plt.plot(center[0], center[-3:].sum() / 1000, 'o', ms=5.5, color='k')
    ldata = len(data)
    cluster_centers = pd.DataFrame(
        {c: members.mean().append(
            pd.Series([len(members) / ldata], ['weight']))
         for c, members in cluster_map.items()}
    ).T.sort_values('T_amb', ignore_index=True)

    cluster_centers.loc[0, 'weight'] = 0  # set weight of extreme demand to 0

    cluster_centers.to_csv(data_path / 'cluster_centers.csv', index=False)
    return cluster_centers


def calc_heat_curve(t_amb, t_min, t_max, t_amb_design=C2K(-12),
                    t_amb_threshold=C2K(20)):
    """Calculate linear heat demand curve.

    Arguments
    ---------
    t_amb: pd.Series
        Series of ambient air temperature data, [K]
    t_min: float
        Flow temperature at threshold ambient air temperature, [K]
    t_max: float
        Flow temperature at design ambient air temperature, [K]
    t_amb_design: float
        Design ambient air temperature, [K], default = 261.15 K
    t_amb_threshold: float
        Threshold ambient air temperature, [K], default = 293.15 K
    """
    m = (t_max - t_min) / (t_amb_design - t_amb_threshold)
    y = t_max - m * t_amb_design
    t_hc = t_amb * m + y
    return t_hc


def create_energysystem(consumers):
    """Create the an energy system model based on demand data.

    Arguments
    ---------
    consumers: Iterable
        Contains the names of consumers to be considered
    """
    # import components
    from .DESTEST_components import DummySource, Network, WasteHeatMflowT, \
        LinkingComponent, Consumer
    # Instantiate central components
    source_el = DummySource('Source_el')
    source_gas = DummySource('Source_gas')
    source_wasteheat = WasteHeatMflowT('Source_WH')
    # Load pipe data
    # Pipe network data is used from https://github.com/ibpsa/project1
    pipe_data = pd.read_csv(data_path / 'pipe_data.csv')
    net = Network('Network', pipe_data=pipe_data)
    lc_central = LinkingComponent('central',
                                  t_in_a=source_wasteheat['T_flow'],
                                  t_out_a=source_wasteheat['T_return'],
                                  t_in_b=net['T_out_a'],
                                  t_out_b=net['T_in_a'])
    lc_central.extend_connection('IN_mflow_a')
    lc_central.extend_connection('OUT_mflow_b')

    # Add central components to comp
    comp = [source_el,
            source_gas,
            source_wasteheat,
            lc_central,
            net,
            ]

    # Set connections of central components
    conn = {
        # # Source - Central LC
        'mflow_centLC_source': [source_wasteheat.OUT_mflow,
                                lc_central.IN_mflow_a],
        # # Central LC - net
        'mflow_centLC_net': [net.IN_mflow, lc_central.OUT_mflow_b],
        # mflow net
        'mflow_net_HP_bus': [net.OUT_mflow],

        # HP - Power
        'el_Bus': [source_el.OUT,
                   lc_central.IN_P,
                   ],
        'gas_Bus': [source_gas.OUT]
    }

    # Add consumers and set connections of decentral components
    for consumer in consumers:
        cons = Consumer(f'C_{consumer}',
                        b_connect=net[f'b_{consumer}'],
                        t_in_a=net['T_out_b'],
                        t_out_a=net['T_in_b'])
        comp.append(cons)
        # Net - con
        conn['mflow_net_HP_bus'].append(cons.IN_mflow_a)
        cons.extend_connection('IN_P_el')
        conn['el_Bus'].append(cons.IN_P_el)
        conn['gas_Bus'].append(cons.IN_P_gas)

    # Instatiate energy system
    ES = System(label="ES", components=comp, connections=conn)

    return ES


def run_destest_case_study(validate_result=False):
    """Run the DESTEST case study.

    Arguments
    ---------
    validate_result: Boolean
        If True, an existing result file (design.csv) is used to create a
        reduced energysystem. The result of this reduced energysystem can
        be compared to the optimistion result, which uses clustered demand
        data.
        If False, a design optimisation with clustered demand data is
        started.
    """
    # Demand data is used from https://github.com/ibpsa/project1
    data = pd.read_csv(data_path / 'data.csv', index_col='time')
    consumer_groups = [col.strip('_Power_W') for col in data.columns
                       if col.startswith('SimpleDistrict_')]
    ES = create_energysystem(consumers=consumer_groups)
    if validate_result:
        print("\n\n\n=============================================")
        print("Solving operational problem for fixed design.")
        print("=============================================\n\n\n")
        # create reference sim based on result_dict using original data
        data['weight'] = 1 / len(data)
        used_data = data.where(data['SimpleDistrict'] >= 0.1).dropna()
        scenarios = used_data['weight']
    else:
        print("\n\n\n=======================")
        print("Solving design problem.")
        print("=======================\n\n\n")
        used_data = cluster(data, n_clusters=11)
        # create Scenarios
        scenarios = dict(used_data[['T_amb', 'weight']].values)

    # Add expressions to ES
    for expr_id in ['investment_costs', 'fixed_costs', 'variable_costs']:
        ES.add_expression(expr_id, ES.aggregate_component_expressions(expr_id))

    # Collect parameters of all components
    params = dict()

    params['Source_el_price'] = 0.1598  # [€/kWh]
    params['Source_gas_price'] = 0.0634  # [€/kWh]
    params['Network_T_amb'] = C2K(used_data['T_amb'])  # °C
    # Waste Heat net
    params['Source_WH_T_flow'] = C2K(30)  # [K]
    params['Source_WH_T_return'] = C2K(20)  # [K]
    params['Source_WH_Q_dot'] = 3000  # [kW]

    # Consumer settings
    t_dict = dict(
        SimpleDistrict_2_3_5_6=[C2K(40), C2K(35)],
        SimpleDistrict_9_12_13_14=[C2K(50), C2K(40)],
        SimpleDistrict_1_4_7_8=[C2K(70), C2K(50)],
        SimpleDistrict_10_11_15_16=[C2K(85), C2K(60)],
     )

    for consumer_group in consumer_groups:
        # Heat Pump Settings
        params[f'HP_LC_C_{consumer_group}_n'] = 4
        params[f'HX_LC_C_{consumer_group}_n'] = 4
        # BES settings
        params[f'BES_C_{consumer_group}_Qdot'] = \
            used_data[f'{consumer_group}_Power_W'] / 1000  # [kW]
        hc = calc_heat_curve(
            t_amb=C2K(used_data['T_amb']),
            t_max=t_dict[consumer_group][0], t_min=t_dict[consumer_group][1])
        params[f'BES_C_{consumer_group}_T_flow'] = hc
        params[f'BES_C_{consumer_group}_dT'] = 15
        params[f'HS_el_C_{consumer_group}_dT'] = 15
        params[f'HS_el_C_{consumer_group}_n'] = 4
        params[f'HS_el_C_{consumer_group}_p_spec'] = 10
        params[f'HS_el_C_{consumer_group}_p_fix'] = 100
        params[f'HS_gas_C_{consumer_group}_dT'] = 15
        params[f'HS_gas_C_{consumer_group}_n'] = 4
        params[f'HS_gas_C_{consumer_group}_p_spec'] = 111
        params[f'HS_gas_C_{consumer_group}_p_fix'] = 4300

    # create Problem
    P = ES.create_problem(
        *make_tac_objective(ES, n=30, i=0.04),
        timesteps=([1], 8760),
        scenarios=scenarios,
        data=params,
        name='DESTEST'
    )

    print()
    print(f'Problem has {P.num_cons} constraints and {P.num_vars} variables.')
    print()

    if validate_result:
        # set design variables according to design.csv
        result_df = pd.read_csv('design.csv', index_col='name')
        for dv in P.design_variables:
            dv.value = result_df.loc[dv.name]['value']
            dv.fix()
            print(f'set {dv.name} to {dv.value}')
        log_name = 'DESTEST_validation.log'
        design_name = 'design_validation.csv'
        operation_name = 'operation_validation.csv'
    else:
        log_name = 'DESTEST.log'
        design_name = 'design.csv'
        operation_name = 'operation.csv'

    m = to_gurobi(P)
    print('Solving...')
    options = dict(  # Options assuming Gurobi 9.1.1
        Seed=SEED,
        NonConvex=2,
        MIPGap=0,
        MIPFocus=1,
        LogFile=log_name,
        OutputFlag=1,
    )
    m.solve(**options)

    print(P.design)
    P.design.to_csv(design_name)
    P.operation.to_csv(operation_name)

    tac = m.objVal
    print(f'\nExpected TAC {tac} €')
    demands = sum(p for p in P.parameters if p.name.endswith('Qdot'))
    tot_demand = P.weighted_sum(demands, symbolic=False) / 1e3
    print(f'\nAnnual heat demand: {tot_demand} MWh')
    print(f'\ncorresponds to: {tac / tot_demand} €/MWh')

    return ES, P, m


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser()

    ap.add_argument('-validate_result', '-vr', action='store_true',
                    default=False)

    if __package__ is None:
        import sys
        from pathlib import Path

        file = Path(__file__).resolve()
        sys.path.append(str(file.parents[1]))
    import DESTEST

    __package__ = DESTEST.__name__
    run_destest_case_study(**vars(ap.parse_args()))
