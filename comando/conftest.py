"""Setup for automatic testing."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHOR: Marco Langiu
from contextlib import contextmanager
import os

import pytest


@pytest.fixture
def reset_backend():
    """Reset the backend to the one before the test after it is done."""
    import comando
    old_backend = comando.get_backend()
    yield
    comando.set_backend(old_backend.__name__)


@pytest.fixture
def clear_components():
    """Clean al defined components after test."""
    yield
    from comando import System
    print('Deleting components:\n' + '\n'.join(System.existing_components))
    System.existing_components.clear()


@pytest.fixture()
def voll_system():
    """Create an example energy system for functional tests.

    Component models taken from:
    @PhdThesis{voll2014automated,
      author            = {Voll, Philip},
      title             = {Automated optimization based synthesis of
                           distributed energy supply systems; 1. Aufl.},
      institution       = {Lehrstuhl für Technische Thermodynamik und Institut
                           für Thermodynamik},
      year              = {2014},
      isbn              = {978-3-86130-474-6},
      url               = {http://publications.rwth-aachen.de/record/228954},
      series            = {Aachener Beiträge zur Technischen Thermodynamik},
      volume            = {1},
    }
    """
    from comando import System
    from components.example_components import Source, Sink, Demand
    from examples.IES.IES_components import CombinedHeatAndPower, Boiler, \
        AbsorptionChiller, CompressionChiller

    # Power:
    GRID = Source('Electricity', price=0.16)
    FEED = Sink('Feedin', compensation=0.1)
    CHP1 = CombinedHeatAndPower('CHP')

    # Heating:
    GAS = Source('Gas', price=0.06)  # €/[kWh gas drawn from grid]
    B1 = Boiler('B1')
    B2 = Boiler('B2')
    HD = Demand('Heat')

    # Cooling: AC, CC
    AC1 = AbsorptionChiller('AC')
    CC1 = CompressionChiller('CC1')
    CC2 = CompressionChiller('CC2')
    CD = Demand('Cooling')

    comp = [GRID, FEED, CHP1, GAS, B1, B2, HD, AC1, CC1, CC2, CD]
    conn = {
        'Power_Bus': [GRID.OUT, FEED.IN, CHP1.POWER_OUT,
                      CC1.IN, CC2.IN],
        'Gas_Bus': [GAS.OUT, B1.IN, B2.IN, CHP1.IN],
        'Heat_Bus': [B1.OUT, B2.OUT, CHP1.HEAT_OUT, HD.IN,
                     AC1.IN],
        'Cooling_Bus': [AC1.OUT, CC1.OUT, CC2.OUT, CD.IN]}

    ES = System('Energy_System', comp, conn)

    ES.parameters
    ES.design_variables
    ES.operational_variables
    ES.expressions
    ES.states

    for id in ['investment_costs', 'fixed_costs', 'variable_costs']:
        ES.add_expression(id, ES.aggregate_component_expressions(id))

    # Breaking symmetries
    ES.add_le_constraint(ES['B1_Qdot_out_nom'], ES['B2_Qdot_out_nom'],
                         'break_boiler_symmetry')
    ES.add_le_constraint(ES['CC1_Qdot_out_nom'], ES['CC2_Qdot_out_nom'],
                         'break_compression_chiller_symmetry')
    yield ES
    System.existing_components.clear()


@pytest.fixture()
def det_data():
    """Return data for a deterministic `VollSystem` formulation.

    Data taken from:
    @PhdThesis{voll2014automated,
      author            = {Voll, Philip},
      title             = {Automated optimization based synthesis of
                           distributed energy supply systems; 1. Aufl.},
      institution       = {Lehrstuhl für Technische Thermodynamik und Institut
                           für Thermodynamik},
      year              = {2014},
      isbn              = {978-3-86130-474-6},
      url               = {http://publications.rwth-aachen.de/record/228954},
      series            = {Aachener Beiträge zur Technischen Thermodynamik},
      volume            = {1},
    }
    """
    data = {  # kW
        'Cooling_demand': [1200, 1300, 2600, 1900],
        'Heat_demand': [2400, 1500, 700, 1500]}
    return data


@pytest.fixture()
def det_voll_problem(voll_system, det_data):
    """Combine the Voll energy system with deterministic data."""
    # voll_system = voll_system()
    ts = ['spring', 'summer', 'autumn', 'winter']
    n = 10
    i = 0.08
    af = ((1 + i) ** n * i) / ((1 + i) ** n - 1)  # annuity factor
    ic = voll_system.get_expression('investment_costs')
    fc = voll_system.get_expression('fixed_costs')
    vc = voll_system.get_expression('variable_costs')
    P = voll_system.create_problem(af*ic + fc, vc, (ts, 8760), name='Test')
    for parameter_name, parameter_value in det_data.items():
        P[parameter_name] = parameter_value
    yield P
    voll_system.existing_components.clear()


@pytest.fixture()
def cleanup_pyomo_files():
    """Remove pyomo- and solver-generated files."""
    import os
    import re
    for f in os.listdir('.'):
        if re.search("tmp.*(pyomo.nl|neos.sol|neos.log)", f):
            os.remove(os.path.join('.', f))


@pytest.fixture()
@contextmanager
def run_in_tmpdir(tmpdir):
    """Run in temporary directory then return where you've been."""
    oldpwd = os.getcwd()  # we start in oldpwd
    os.chdir(tmpdir)  # then we go to tmpdir
    print('RUNNING IN:', os.getcwd())
    try:
        yield
    finally:
        os.chdir(oldpwd)  # and return to oldpwd


@pytest.fixture()
def test_problem(scenarios, timesteps):
    """Create a simple instance of the Rosenbrok problem for testing."""
    import comando

    x = comando.Variable('x', bounds=(-1.5, 1.5))
    y = comando.VariableVector('y', bounds=(-0.5, 2.5))
    c = comando.Parameter('c', 2)  # scalar parameter - constant
    pv = comando.Parameter('p')  # parameter vector - one entry per (s, t)
    do = (1 - x) ** 2
    oo = pv * (y - x ** 2) ** 2
    constraints = {'c1': (x - 1) ** 3 - y + 1 <= 0,
                   'c2': x + y - c <= 0}
    P = comando.Problem(do, oo, constraints, timesteps=timesteps,
                        scenarios=scenarios, name='Rosenbrok'
                        + ('_stoch' if scenarios else ''))
    pv.value = y.value * 100

    P.c = comando.cyclic  # Just so we test if this is can be serialized!

    return P
