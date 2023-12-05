"""Case study for building demand response."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Florian Joseph Baader, Marco Langiu
import pickle

from comando.core import System
from comando.utility import make_mayer_objective
from comando.interfaces.pyomo_dae import to_pyomo_dae


def run_BDR_case_study():
    """Run the building demand response case study."""
    from .BDR_components import AirHeatPump, Mass, HeatTransport, Grid

    ES = System('ES')
    T_amb = ES.make_parameter('T_amb')

    # components of the building model:
    #                         [kg/m³]   [J/kg K]  [m³]        [K]
    AirMass = Mass('AirMass', rho=1.19, C_p=1007, V=2.5*13*6, T_0=294.15,
                   T_bounds=True)
    WallMass = Mass('WallMass', rho=218, C_p=1000, V=0.435*2.5*(13+2*6),
                    T_0=284.15)
    CoreMass = Mass('CoreMass', rho=218, C_p=1000, V=0.125*(13*6), T_0=294.15)
    #                                           [W]    number of outflows
    AHP = AirHeatPump('AHP', T_amb, Q_max=3000, n_outflows=2)

    # alpha [W /(m² K)]; Area [m²]; T_A [K]; T_B [K]
    HTCoreAir = HeatTransport('HTCoreAir', alpha=2.7, Area=13*6,
                              T_A=CoreMass['T'], T_B=AirMass['T'])
    HTAirWall = HeatTransport('HTAirWall', alpha=2.7, Area=2.5*(13+2*6),
                              T_A=AirMass['T'], T_B=WallMass['T'])
    HTWallAmbient = HeatTransport('HTWallAmbient', alpha=5, Area=2.5*(13+2*6),
                                  T_A=WallMass['T'], T_B=T_amb)
    HTCoreWall = HeatTransport('HTCoreWall', alpha=2.7, Area=0.125*(13+2*6),
                               T_A=CoreMass['T'], T_B=WallMass['T'])

    PG = Grid('PG')

    # make a list with all components
    components = [AirMass, WallMass, CoreMass, AHP, HTAirWall, HTWallAmbient,
                  HTCoreAir, HTCoreWall, PG]
    # define connections
    connections = {'HeatToCoreMass': [CoreMass.Q_IN, HTCoreAir.A,
                                      AHP.Qdot_out_1, HTCoreWall.A],
                   'HeatToAirMass': [AirMass.Q_IN, HTCoreAir.B, HTAirWall.A,
                                     AHP.Qdot_out_2],
                   'HeatToWallMass': [WallMass.Q_IN, HTAirWall.B,
                                      HTWallAmbient.A, HTCoreWall.B],
                   'Power': [AHP.P_IN, PG.OUT]}
    for c in components:
        ES.add(c)
    for bus_id, connectors in connections.items():
        ES.connect(bus_id, connectors)
    # max = 24 h; min = 1/StepsPerHour;
    # should be n*(1/StepsPerHour) with an element N
    hours = 24
    # number of intervals per hour for which data are given
    StepsPerHour = 4
    ncp = 4
    # number of finite elements for which controls are set picewise constant
    control_factor = 1
    nfe = int(hours*StepsPerHour*control_factor)

    P = ES.create_problem(operational_objective=make_mayer_objective(ES),
                          timesteps=([3600*i/StepsPerHour for i in
                                      range(1, hours*StepsPerHour+1)],
                                     hours*3600))
    # Quarter-hourly data
    ambient_T_data = [278.32, 278.38, 278.4, 278.51, 278.43, 278.51, 278.55,
                      278.61, 278.63, 278.72, 278.42, 278.45, 278.52, 278.38,
                      278.38, 278.41, 278.27, 277.89, 277.8, 277.81, 277.97,
                      278.06, 277.94, 277.99, 278.04, 278.16, 277.96, 277.58,
                      277.79, 277.91, 277.86, 277.89, 277.25, 276.87, 277.4,
                      277.65, 277.7, 277.48, 277.35, 277.18, 277.04, 276.98,
                      277.09, 277.29, 277.32, 277.77, 278.44, 278.45, 278.41,
                      278.5, 278.2, 278.29, 278.54, 278.63, 278.82, 279.11,
                      279.24, 279.59, 279.33, 279.2, 279.47, 279.48, 279.64,
                      279.51, 279.33, 279.64, 280.11, 280.22, 279.79, 280.15,
                      280, 279.16, 277.28, 276.69, 277.31, 278.04, 277.58,
                      277.92, 278.14, 279.42, 279.7, 279.34, 279.6, 281.37,
                      280.72, 280.52, 280.58, 280.34, 279.39, 281.01, 281.69,
                      281.75, 281.82, 281.73, 282.81, 282.47]

    # Hourly data
    electricity_Price_data = [0.002135, 0.002098, 0.002119, 0.002064, 0.002113,
                              0.002178, 0.002973, 0.00365, 0.004044, 0.003839,
                              0.003457, 0.003055, 0.002579, 0.002309, 0.002313,
                              0.0023, 0.002245, 0.002428, 0.003191, 0.00339,
                              0.003585, 0.00313, 0.002433, 0.002603]

    # time series of ambient temperature [K]
    P['ES_T_amb'] = ambient_T_data
    # time series of electricity prices [ct./Wh]
    P['PG_ElectricityPrice'] = [p for p in electricity_Price_data
                                for repetition in range(StepsPerHour)]
    P['AirMass_T_min'] = 293.15
    P['AirMass_T_max'] = 295.15

    P.controls = {AHP['b_on']}
    pyomo_dae_options = dict(method='dae.collocation', ncp=ncp, nfe=nfe)
    dynamic_options = dict(controls=P.controls,
                           control_factor=control_factor,
                           pyomo_dae_options=pyomo_dae_options)

    m = to_pyomo_dae(P, **dynamic_options)

    from pyomo.common.errors import ApplicationError
    try:
        m.solve(solver='gurobi', keepfiles=False, tee=True)
    except (ApplicationError, ValueError) as e:
        # If we got a demo version pyomo will state solver status as 'aborted'!
        if isinstance(e, ValueError) and 'aborted' not in str(e):
            raise  # This is a different ValueError, raise it!
        m.solve('gurobi', remote=True, options={'mipgap': 0})

    with open('results.pickle', 'wb') as f:
        pickle.dump((P.data.getter(), P.operation.T), f)

    return m.objective()


def plot(P):
    """Generate a plot of the results."""
    from matplotlib import pyplot as plt

    plt.rcParams.update({
        "font.family": "serif",  # use serif/main font for text elements
        "text.usetex": True,     # use inline math for ticks
        "pgf.rcfonts": False,    # don't setup fonts from rc parameters
        'hatch.linewidth': 0.5   # reduce hatch linewidth
        })

    with open('results.pickle', 'rb') as f:
        data, operation = pickle.load(f)

    temp = data[['ES_T_amb', 'AirMass_T_min', 'AirMass_T_max']]
    temp = temp.join(operation[['AirMass_T', 'CoreMass_T', 'WallMass_T']])
    temp = temp - 273.15  # Kelvin to Celsius
    temp.columns = [
        'Environment',
        r'$T_{\mathrm{air}}^{\mathrm{min}}$',
        r'$T_{\mathrm{air}}^{\mathrm{max}}$',
        'Air',
        'Core',
        'Wall'
    ]
    nrg = (operation['AHP_Q_HP'] * 1e-3)
    cost = (data['PG_ElectricityPrice'] * 1e3)

    fig, axs = plt.subplots(3, 1, sharex=True, sharey=False,
                            figsize=(8/2.54, 14/2.54),
                            gridspec_kw={'hspace': 0,
                                         'height_ratios': [6, 4, 4]})
    temp.plot(ax=axs[0])
    axs[0].set_xlim(0, 24 * 3600)  # time limits in s
    xticks = range(0, 25, 4)
    axs[0].set_xticks([tic * 3600 for tic in xticks])
    axs[0].set_xticklabels(xticks)  # time labels in hours
    axs[0].set_xlabel('$t$ [h]')
    axs[0].set_ylabel('$T$ [°C]')
    axs[0].legend(ncol=3, labelspacing=.1, columnspacing=.5, handlelength=.8,
                  loc=(0.03, 0.35))
    nrg.plot(ax=axs[1])
    axs[1].set_ylabel(r'$\dot Q_{\mathrm{HP}}$ [kW]')  # energy use in kW
    cost.plot(ax=axs[2])
    axs[2].set_ylabel(r'$C^{\mathrm{elec}}$ [ct/kWh]')  # cost in € ct / kWh
    plt.show()
    plt.savefig('results.png')


if __name__ == '__main__':
    if __package__ is None:
        import sys
        from pathlib import Path
        file = Path(__file__).resolve()
        sys.path.append(str(file.parents[1]))
    import BDR
    __package__ = BDR.__name__
    run_BDR_case_study()
