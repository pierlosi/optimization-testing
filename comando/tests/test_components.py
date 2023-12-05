"""Tests for components."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import pytest

from comando.core import System
from comando.utility import make_tac_objective
from comando.interfaces.pyomo import to_pyomo
from components.example_components import Source, Demand
from examples.IES.IES_components import Boiler, Battery


def test_Boiler_use(clear_components):
    G = Source('Gas')
    B = Boiler('Boiler', optional=False)
    H = Demand('Heat')
    connections = {'gas supply': [G.OUT, B.IN],
                   'demand satisfaction': [B.OUT, H.IN]}
    ES = System('Sys', components=[G, B, H], connections=connections)
    P = ES.create_problem(*make_tac_objective(ES), (range(1, 4), 8760))

    # TODO: Create problem
    demand = (1, 1.5, 2)  # in kW
    data = {'Gas_price': 0.3, 'Heat_demand': demand}
    for param, value in data.items():
        P[param] = value

    m = to_pyomo(P)

    # DEBUG: Necessary when pyomo looses track of its logger for some reason...
    global logger
    logger = lambda: None  # Dummy logger, just need a namespace...
    logger.warning = print
    import pyomo
    pyomo.repn.util.logger = logger
    # END DEBUG

    from pyomo.common.errors import ApplicationError
    try:
        m.solve('baron', keepfiles=False)
    except (ApplicationError, ValueError) as e:
        # If we got a demo version of baron we'll get a licensing error and
        # pyomo will complain that solver status is 'aborted'!
        if isinstance(e, ValueError) and 'aborted' not in str(e):
            raise  # This is a different ValueError, raise it!
        from tests import skip_long
        if skip_long:
            pytest.skip()
        from pyomo.opt.parallel.manager import ActionManagerError
        try:
            m.solve('knitro', remote=True, options={'maxtime_real': 2})
        except ActionManagerError as e:
            if '[]' in str(e):
                pytest.xfail('NEOS is unavailable')

    # Test constraint violations
    for i in m.constraints:
        con = m.constraints[i]
        lb_violation = 0 if con.lower is None else max(0, con.lower() - con())
        ub_violation = 0 if con.upper is None else max(0, con() - con.upper())
        assert lb_violation == pytest.approx(0, abs=1e-5)
        assert ub_violation == pytest.approx(0, abs=1e-5)

    # We need to meet demand...
    assert tuple(ES['Boiler_output'].value) == pytest.approx(demand)
    # ...and consequently need to install power that can meet maximum demand!
    assert ES['Boiler_Qdot_out_nom'].value == pytest.approx(max(demand),
                                                            abs=1e-4)

    # Finally consumed gas resources should match the boiler input
    diff = G['use'].value - B.get_expression('input').value
    for val in diff.values:
        assert val == pytest.approx(0.0, rel=1e-6)


def test_Storage_use(clear_components):
    # Here we create a (nonsensical) system in which we only compensate the
    # self discharge of a battery.
    G = Source('Grid')
    Battery.min_cap = 0.15  # 150 kWh
    Battery.c_loss = 4.2e-5  # small self discharge
    BAT = Battery('Battery')
    min_soc = 0.015
    BAT['soc'].lb = min_soc
    BAT.parameters_dict['Battery_soc_init'].value = min_soc
    D = Demand('Power')
    connections = {'supply': [G.OUT, BAT.IN],
                   'demand': [BAT.OUT, D.IN]}
    ES = System('Sys', components=[G, BAT, D], connections=connections)
    P = ES.create_problem(*make_tac_objective(ES), (range(2), 8760))
    data = {'Grid_price': 0.3, 'Power_demand': 0}
    for param, value in data.items():
        P[param] = value

    # FIXME: P.design_variables.pop() pops from the actual design variables!

    m = to_pyomo(P)

    from pyomo.common.errors import ApplicationError
    try:
        m.solve('baron', keepfiles=False)
    except (ApplicationError, ValueError) as e:
        # If we got a demo version of baron we'll get a licensing error and
        # pyomo will complain that solver status is 'aborted'!
        if isinstance(e, ValueError) and 'aborted' not in str(e):
            raise  # This is a different ValueError, raise it!
        from tests import skip_long
        if skip_long:
            pytest.skip()
        from pyomo.opt.parallel.manager import ActionManagerError
        try:
            m.solve('knitro', remote=True, options={'maxtime_real': 2})
        except ActionManagerError as e:
            if '[]' in str(e):
                pytest.xfail('NEOS is unavailable')
    exp_cap = Battery.min_cap  # 150 kWh
    exp_SOC = min_soc  # 15 kWh
    exp_discharge = Battery.c_loss * exp_SOC  # 0.000623 kWh / h = kW
    exp_in = exp_discharge / BAT.charge_eff  # 0.0006833 kW

    assert P['Battery_capacity_nom'].value == pytest.approx(exp_cap)
    assert all(v == pytest.approx(1) for v in P['Battery_b_in'].value)
    assert all(v == pytest.approx(exp_in) for v in P['Battery_input'].value)
    assert all(v == pytest.approx(exp_in) for v in P['Grid_use'].value)
    assert all(v == pytest.approx(exp_SOC) for v in P['Battery_soc'].value)
    assert all(v == pytest.approx(0) for v in P['Battery_output'].value)
    assert all(v == pytest.approx(0) for v in P['Battery_socdot'].value)
    assert all(v == pytest.approx(0) for v in P['Battery_b_out'].value)

# TODO: Similar tests for other components
