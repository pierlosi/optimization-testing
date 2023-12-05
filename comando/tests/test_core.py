"""Tests for the COMANDO core."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import pytest

import comando
from comando import Domain, Symbol, Parameter, Variable, VariableVector, \
    Component, System, Problem, INF


def test_slots():
    if comando.sympy.__version__ == '1.7':
        pytest.xfail('See https://github.com/sympy/sympy/issues/20567')
    for i, ty in enumerate([Symbol, Parameter, Variable, VariableVector]):
        s = ty(f's{i}')
        with pytest.raises(AttributeError):
            s.__dict__  # should not have __dict__ due to slots


def test_variable():
    var = Variable('v', bounds=(3.141, 6.282), init_val=5)
    assert var.name == 'v'
    assert var.domain == Domain.REAL
    assert var.is_integer is False
    assert var.bounds == (3.141, 6.282)
    assert var.indexed is False
    assert var.value == 5
    data = {'case0': 3.1, 'case1': 6.2}
    with pytest.raises(TypeError):
        var.value = data
    assert var.indexed is False
    var.fix()
    assert all(val == 5 for val in (var.value, *var.bounds))
    var.fix(4)
    assert all(val == 4 for val in (var.value, *var.bounds))
    var.fix(4.5)
    assert all(val == 4.5 for val in (var.value, *var.bounds))
    var.unfix()
    assert var.bounds == (3.141, 6.282)
    assert var.value == 4.5
    with pytest.raises(ValueError):
        var.fix(1)  # not within original bounds!

    bin_var = Variable('bin', domain=Domain.BINARY)
    assert bin_var.domain == Domain.BINARY
    assert bin_var.bounds == (0, 1)
    assert bin_var.is_integer is True
    assert bin_var.indexed is False
    # assert bin_var.indexed is True

    # INFO: old indexed Variable features
    # TODO: This should fail (data is not integral)!
    # with pytest.xfail():
    #     with pytest.raises(ValueError):
    #         bin_var.value = data
    #         for key, value in bin_var.value.items():
    #             assert data[key] == value
    # for i in range(2):
    #     data[f'case{i}'] = i
    # bin_var.value = data
    # for key, value in bin_var.value.items():
    #     assert data[key] == value

    bin_var.fix(0)
    assert all(val == 0 for val in (bin_var.value, *bin_var.bounds))
    bin_var.fix(1)
    assert all(val == 1 for val in (bin_var.value, *bin_var.bounds))
    with pytest.raises(ValueError):
        bin_var.fix(2)  # not within original bounds!
    with pytest.warns(UserWarning):
        bin_var.fix(0.5)  # will warn but round to nearest even number (here 0)
    with pytest.raises(ValueError):
        bin_var.fix(1.00001)  # outside of bounds: will raise
    bin_var.unfix()
    assert bin_var.bounds == (0, 1)
    assert bin_var.value == 0


def test_variable_vector():
    y = comando.VariableVector('y', init_val=42)
    assert y.indexed is True
    assert y.is_expanded is False
    with pytest.raises(RuntimeError):
        y.value = 1
    with pytest.raises(IndexError):
        y[1]

    # Bounds in the unexpanded state
    assert y._bounds == y.bounds == (y.lb, y.ub) == (-INF, INF)
    y.bounds = (0, 100)
    assert y._bounds == y.bounds == (y.lb, y.ub) == (0, 100)
    y.ub = 200
    y.lb = 0
    assert y._bounds == y.bounds == (y.lb, y.ub) == (0, 200)
    y.fix(4)
    assert y.bounds == (4, 4)
    y.fix(5)
    assert y.bounds == (5, 5)
    y.unfix()
    assert y.bounds == (0, 200)
    assert y.value is None
    with pytest.raises(ValueError):
        y.fix(300)  # not within original bounds!

    # TODO: Test with other index types, e.g., float, string, datetime, ...
    index = [1, 2, 3]
    y.instantiate(index)
    assert y.is_expanded is True
    assert y[1] == y.expansion[1]  # getitem works...
    with pytest.raises(IndexError):
        y[0]  # ... but only for elements of the index!

    # y is properly expanded and init_val is taken as initial value
    assert [*y.value] == [42, 42, 42]
    y.value = 1  # Setting value with scalar
    assert [*y.value] == [1, 1, 1]
    y.value = {i: i for i in index}  # Setting value with Mapping
    assert [*y.value] == index
    y.value = y.value + 1  # Setting value with Series
    assert [*y.value] == [2, 3, 4]

    # Bounds in the expanded state
    y.ub = 50
    assert y._bounds == (0, 50)
    assert all(ub == 50 for ub in y.ub)
    y.bounds = None, {1: 150}
    assert y._bounds == (-INF, 150)
    y.lb = {3: 5}
    y.ub = {3: None}
    assert y._bounds == (-INF, INF)
    assert all(val == exp for val, exp in zip(y.lb, [-INF, -INF, 5]))
    assert all(val == exp for val, exp in zip(y.ub, [150, 50, INF]))

    with pytest.raises(ValueError):
        y.fix()  # value is still 2, 3, 4 but lb[3] is 5 by now!
    y.value = 42
    y.fix()
    assert all(lbi == vi == ubi for lbi, vi, ubi in zip(y.lb, y.value, y.ub))
    y.fix(10)
    assert all(lbi == vi == ubi for lbi, vi, ubi in zip(y.lb, y.value, y.ub))
    y.unfix()
    assert all(val == 10 for val in y.value)
    assert all(val == exp for val, exp in zip(y.lb, [-INF, -INF, 5]))
    assert all(val == exp for val, exp in zip(y.ub, [150, 50, INF]))

    y[2].value = 9  # __getitem__ + setting value
    y[2] = 5  # __setiem__: equivalent to y[2].value = value
    assert y.value[2] == 5

    assert [*y.expansion] == [*y]  # __iter__: Allows iterating directly over y

    # Sanity reminder: Look no further
    Note: lambda: lambda: lambda: lambda: lambda: lambda: lambda: All is well

    # def mymatmul(self, other):
    #     try:
    #         return self.expansion @ other.value
    #     except AttributeError:
    #         return (self.expansion * other).sum()
    #
    # comando.VariableVector.__matmul__ = mymatmul
    #
    # # Implementing multiple timestep refrences via matrix multiplication
    # def_value = {i: 0 for i in index}
    # time = [0]  # Initial Value
    # e = dict()
    # for i in index:
    #     time.append(comando.Parameter(f't{i}'))
    #     time[-1].value = def_value
    #     time[-1][i] = 1
    #     e[f'switching_@_{i}'] = y @ time[i] - y @ time[i-1]
    # e
    #
    # assert not any(l - r for l, r in \
    #                zip(e.values(), [y[1], y[2] - y[1], y[3] - y[2]]))

    # Test setting bounds with iterables of length 2.
    # (This has caused incorrect behavior before)
    y.instantiate([1, 2])
    y.lb = 1
    assert all(val == 1 for val in y.lb)
    assert all(val == INF for val in y.ub)


def test_component(clear_components):
    C = Component('Component')

    with pytest.raises(RuntimeError) as e:
        Component('Component')  # Labels can only be used once
    assert 'already been defined' in str(e)

    eta = C.make_parameter('efficiency')
    demand = C.make_parameter('demand')
    max_output = C.make_design_variable('max_output', domain=Domain.REAL,
                                        bounds=(40, 100), init_val=40)
    output = C.make_operational_variable('output', domain=Domain.REAL,
                                         bounds=(40, None), init_val=40)
    C.add_eq_constraint(output, demand)
    active = C.make_operational_variable('active', domain=Domain.INTEGER,
                                         bounds=(0, 1), init_val=1)
    C.add_le_constraint(output, max_output * active)
    C.add_expression('investment_costs', max_output * 1000)
    resource_price = C.make_parameter('resource_price')
    C.add_expression('variable_costs', output / eta * resource_price)

    assert C.parameters == {eta, demand, resource_price}
    expected_parameter_names = ['Component_demand',
                                'Component_efficiency',
                                'Component_resource_price']
    assert sorted(par.name for par in C.parameters) == expected_parameter_names

    assert C.design_variables == {max_output}
    assert C.operational_variables == {output, active}

    constraints = [con for con in C.constraints_dict.values()]
    assert str(constraints[0].__class__).endswith("Equality'>")
    # Symengine order seems to be lexical!
    assert set(constraints[0].args) == {output, demand}
    assert str(constraints[1].__class__).endswith("LessThan'>")
    assert constraints[1].lts == output
    assert constraints[1].gts == max_output * active

    Tin = C.make_operational_variable('Tin', bounds=(0, None))
    Tout = C.make_operational_variable('Tout', bounds=(0, None))
    C.add_le_constraint(Tin, Tout)  # Constraint without a name
    # Internally one is created automatically
    assert 'Component_Tin ≤ Component_Tout' in C._constraints_dict
    # Externally it is prefixed with the Component label (as all others)
    assert 'Component_Component_Tin ≤ Component_Tout' in C.constraints_dict

    state, der = C.make_state('Total_output', rate_of_change=output)
    assert state.name == 'Component_Total_output'
    state_init, der_, roc = C.states_dict[state]
    assert isinstance(state_init, comando.Parameter)
    assert state_init.value is None
    assert der_ == der
    assert roc == output

    C.BLA = 'Some stuff you stick to your Component for whatever reason'
    inp = output/eta
    with pytest.raises(RuntimeError) as e:
        C.add_connectors('BLA', inp)  # First possible signature (1 Connector)
    assert 'already has an attribute called BLA' in str(e)
    with pytest.raises(RuntimeError) as e:
        C.add_connectors(IN=inp, BLA=inp*10)  # Second signature (n Connectors)
    assert 'already has an attribute called BLA' in str(e)
    # Even though A was valid, it was not added, since BLA was not!
    assert C.connectors == {}
    C.add_input('IN', inp)
    # Connectors may also be used to represent 'flow of information', note
    # however that the user has to take care to combine such connectors in
    # pairs only and either taking care of the sign manually (as done here),
    # or using input and output connectors.
    C.add_connectors(TIN=Tin, TOUT=-Tout)
    assert len(C.connectors) == 3

    # TODO: C.__getitem__

    # IDEA: Make the above sample component a fixture. We could then use it in
    #       integration tests with the EnergySystem and scenario generation.
    # # Independent events...
    # efficiency_scenarios = {'high_efficiency': 0.25,
    #                         'medium_efficiency': 0.75}
    # demand_scenarios = {'high_demand': 0.3,
    #                     'medium_demand': 0.2,
    #                     'low_demand': 0.5}
    # # ...are combined to obtain scenarios with associated probabilities
    # scenarios = {es + ' & ' + ds: p_es * p_ds
    #              for es, p_es in efficiency_scenarios.items()
    #              for ds, p_ds in demand_scenarios.items()}


def test_system(clear_components):
    # TODO: Replace the example components with fixtures such as the sample
    #       Component defined in test_component above!
    from components.example_components import Source, Demand
    from examples.IES.IES_components import Boiler
    G = Source('Gas')
    B = Boiler('Boiler')

    # specifying that the gas source is hooked up to the boiler
    connections = {'gas supply': [G.OUT, B.IN]}

    # creating the energy system
    ES = System('Sys', components=[G, B], connections=connections)
    # iterating over System gives components, but not necessarily in order!
    assert {*ES} == {G, B}

    # The boiler has an unconnected output!
    assert len(ES.get_open_connectors()) == 1

    # Create, add and connect a heat demand
    H = Demand('Heat')
    ES.add(H)
    assert H in ES.components

    ES.connect('demand satisfaction', [B.OUT, H.IN])
    assert len(ES.connections) == 2
    assert len(ES.get_open_connectors()) == 0

    # Detach B_OUT from the demand_satisfection bus
    ES.detach('demand satisfaction', [B.OUT])
    with pytest.raises(KeyError) as e:
        ES.detach('demand satisfaction', [B.OUT])  # Can't attach again
    assert 'not part of bus' in str(e)
    assert 'demand satisfaction' in ES.connections
    assert ES.constraints_dict['Sys_demand satisfaction'] == \
        comando.Eq(H.IN.expr, 0)
    ES.detach('demand satisfaction', [H.IN])
    assert 'demand satisfaction' not in ES.connections

    ES.connect('demand satisfaction', [B.OUT, H.IN])
    assert 'demand satisfaction' in ES.connections
    assert ES.constraints_dict['Sys_demand satisfaction'] == \
        comando.Eq(B.OUT.expr + H.IN.expr, 0)
    ES.detach('demand satisfaction')
    assert 'demand satisfaction' not in ES.connections
    assert 'demand satisfaction' not in ES.constraints_dict


def test_nested_system(clear_components):

    class C(Component):
        def __init__(self, label):
            super().__init__(label)
            C_in = self.make_operational_variable('in')
            C_out = self.make_operational_variable('out')
            self.add_eq_constraint(C_in, 0.5 * C_out, 'transformation')
            self.add_input('IN', C_in)
            self.add_output('OUT', C_out)

    class Subsystem(System):
        def __init__(self, label):
            super().__init__(label)
            comps = {l: C(f'{label}_{l}') for l in 'ABC'}
            for comp in comps.values():
                self.add(comp)
            self.connect('internal', [comps['A'].OUT,
                                      comps['B'].OUT,
                                      comps['C'].IN])
            self.expose_connector(comps['A'].IN, 'A_IN')
            self.expose_connector(comps['B'].IN, 'B_IN')
            self.expose_connector(comps['C'].OUT)

    subsystems = [Subsystem(f'S{i}') for i in range(2)]
    S = System('S', subsystems)

    S.connect('01A', [subsystems[0].OUT, subsystems[1].A_IN])
    S.connect('10B', [subsystems[1].OUT, subsystems[0].B_IN])
    S.close_connector(subsystems[1].B_IN)
    subsystems[1].extend_connection('internal')
    S.connect('1int0A', [subsystems[1].internal, subsystems[0].A_IN])

    S.connections
    S.constraints_dict
    S0_A_in = S.operational_variables_dict['S0_A_in']
    S0_A_out = S.operational_variables_dict['S0_A_out']
    S0_B_in = S.operational_variables_dict['S0_B_in']
    S0_B_out = S.operational_variables_dict['S0_B_out']
    S0_C_in = S.operational_variables_dict['S0_C_in']
    S0_C_out = S.operational_variables_dict['S0_C_out']
    S1_A_in = S.operational_variables_dict['S1_A_in']
    S1_A_out = S.operational_variables_dict['S1_A_out']
    S1_B_in = S.operational_variables_dict['S1_B_in']
    S1_B_out = S.operational_variables_dict['S1_B_out']
    S1_C_in = S.operational_variables_dict['S1_C_in']
    S1_C_out = S.operational_variables_dict['S1_C_out']
    Eq = comando.Eq
    expected = {'S_01A': Eq(-S0_C_out + S1_A_in, 0),
                'S_10B': Eq(S0_B_in - S1_C_out, 0),
                'S_IN': Eq(S1_B_in, 0),
                'S0_internal': Eq(-S0_A_out - S0_B_out + S0_C_in, 0),
                'S_1int0A': Eq(S0_A_in - S1_A_out - S1_B_out + S1_C_in, 0),
                'S0_A_transformation': Eq(S0_A_in, 0.5*S0_A_out),
                'S0_B_transformation': Eq(S0_B_in, 0.5*S0_B_out),
                'S0_C_transformation': Eq(S0_C_in, 0.5*S0_C_out),
                'S1_C_transformation': Eq(S1_C_in, 0.5*S1_C_out),
                'S1_A_transformation': Eq(S1_A_in, 0.5*S1_A_out),
                'S1_B_transformation': Eq(S1_B_in, 0.5*S1_B_out)}
    assert S.constraints_dict == expected


# TODO: The syntax in the following two tests seems more natural
@pytest.mark.skip(reason="dataframe getitem is by column, not by row!")
def test_problem_design_operation_getitem():
    dv = comando.Variable('dv')
    ov = comando.VariableVector('ov')
    P = comando.Problem(dv, ov, timesteps=(['test'], 1))
    # NOTE: This currently doesn't work
    assert dv.value == P.design['dv']
    assert ov.value == P.operation['ov']

    # NOTE: This currently works
    with pytest.raises(KeyError):
        P.design.T['dv']
    with pytest.raises(KeyError):
        P.operation.T['ov']


@pytest.mark.skip(reason="dataframe setitem is by column, not by row!")
def test_problem_design_operation_setitem():
    import pandas as pd
    dv = comando.Variable('dv')
    ov = comando.VariableVector('ov')
    P = comando.Problem(dv, ov, timesteps=(['test'], 1))
    P.design['dv'] = 1
    # NOTE: This currently doesn't work
    assert dv.value == 1
    for data in {'test': 1}, [2, ], 3:
        P.operation['ov'] = data
        assert all(ov.value == pd.Series(data, ['test']))


def test_problem(clear_components):
    """Test the creation of a simple problem."""
    from collections.abc import Mapping
    from pandas import Series

    timesteps_data = [(range(1, 5), 5), (['t1', 't2', 't3', 't4'], 0),
                      {'t1': 1, 't2': 2, 't3': 3, 't4': 4}]
    scenarios_data = [None, ['s1', 's2', 's3'], range(1, 4),
                      {'s1': 1/3, 's2': 1/3, 's3': 1/3}]
    for scenarios in scenarios_data:
        for timesteps in timesteps_data:
            P = Problem(timesteps=timesteps, scenarios=scenarios)
            t_len = len(timesteps) if isinstance(timesteps,
                                                 (Mapping, Series)) \
                else len(timesteps[0])
            s_len = 1 if scenarios is None else len(scenarios)
            assert len(P.index) == t_len * s_len

    from components.example_components import Source, Demand
    from examples.IES.IES_components import Boiler
    from comando.utility import make_tac_objective, indexed, split
    G = Source('Gas', price=0.06)
    B = Boiler('Boiler')
    H = Demand('Heat')
    # specifying that the gas source is hooked up to the boiler
    connections = {'gas supply': [G.OUT, B.IN],
                   'demand satisfaction': [B.OUT, H.IN]}

    # creating the energy system
    ES = System('Sys', components=[G, B, H], connections=connections)

    scenarios = ['s1', 's2', 's3', 's4']
    P = ES.create_problem(*make_tac_objective(ES), timesteps=(range(4), 4),
                          scenarios=scenarios)
    P['Heat_demand'] = 1, 2, 3, 4
    assert all(all(P['Heat_demand'].value[s] == (1, 2, 3, 4))
               for s in scenarios)
    with pytest.raises(ValueError) as e:
        P['Heat_demand'] = 1, 2, 3
    assert str(e.exconly()) == 'ValueError: Value must be a scalar, a ' \
        'Mapping or an Iterable with appropriate length (4 for ' \
        'time-dependent data, or 16 for data depending on both scenario ' \
        'and time)!'
    scenarios = ['s1', 's2', 's3']
    P = ES.create_problem(*make_tac_objective(ES), timesteps=(range(4), 4),
                          scenarios=scenarios)
    P['Heat_demand'] = 1, 2, 3, 4  # time-dependent data
    assert all(all(P['Heat_demand'].value[s] == (1, 2, 3, 4))
               for s in scenarios)
    P['Heat_demand'] = 1, 2, 3  # scenario dependent data
    assert all(all(P['Heat_demand'].value[s] == val)
               for s, val in zip(scenarios, [1, 2, 3]))
    with pytest.raises(ValueError) as e:
        P['Heat_demand'] = 1, 2
    assert str(e.exconly()) == 'ValueError: Value must be a scalar, a ' \
        'Mapping or an Iterable with appropriate length (4 for ' \
        'time-dependent data, 3 for scenario-dependent data or 12 for data ' \
        'depending on both scenario and time)!'

    # Post-initialization change of timesteps and scenarios
    timesteps = P.timesteps['s1']
    timesteps[4] = 1
    P.timesteps = timesteps
    P.data
    assert all(all(P.timesteps[s] == timesteps) for s in scenarios)
    assert all(P['Heat_demand'].value.index == P.index)
    P.scenarios = scenarios[:2]
    assert all(a == b for a, b in zip(P.scenarios, scenarios[:2]))
    assert all(P['Heat_demand'].value.index == P.index)

    P = ES.create_problem(*make_tac_objective(ES), timesteps=(range(4), 4))
    for i, (t, timesteps) in enumerate(P.timesteps.items()):
        assert i == t
        assert timesteps == 1

    P['Heat_demand'] = data = 2.2, 4.2, 3, 2.5  # set all via iterable
    P['Gas_price'] = {1: 0.05, 3: 0.07}  # change only selected values via dict

    problem_data = P.data
    assert all(problem_data['Heat_demand'].values == data)
    assert all(P['Heat_demand'].value == problem_data['Heat_demand'])
    assert all(problem_data['Gas_price'] == [0.06, 0.05, 0.06, 0.07])
    assert all(P['Gas_price'].value == problem_data['Gas_price'])

    new_problem_data = problem_data.copy()
    new_problem_data['Gas_price'] *= 2
    new_problem_data['Heat_demand'] += 0.2
    P.data = new_problem_data
    assert all(P['Heat_demand'].value == new_problem_data['Heat_demand'])
    assert all(P['Gas_price'].value == new_problem_data['Gas_price'])

    # Some exemplary initialization using P's design and operation properties
    # TODO: we should improve the syntax for individual assignments!
    assert P.design.value['Boiler_Qdot_out_nom'] == 14
    # Making Boiler large enough to supply head demand
    max_dem = new_problem_data['Heat_demand'].max()
    P.design = {'Boiler_Qdot_out_nom': max_dem}
    assert P['Boiler_Qdot_out_nom'].value == max_dem
    # This affects the Variable in the boiler!
    assert B['Qdot_out_nom'].value == max_dem
    Qnom = B['Qdot_out_nom']

    # TODO: Decide on identical orientation of data (index/symbol or
    #       symbol/index) for both parameter and operational data!
    Qrel = B['Qdot_out_rel']
    # Corresponding part load can be computed
    Qrel.value = comando.utility.evaluate(P['Heat_demand']/Qnom.value)

    # TODO: figure out why solve causes error on CI
    # solve input-output relation with assumed Q_nom for input variable
    # expr = comando.solve(B.constraints_dict['Boiler_input_output_relation'],
    #                   B['Qdot_in'])[0]
    # guess = comando.utility.evaluate(expr)
    # # then set values for each timestep
    # op_data = P.operation
    # for t in P.timesteps.keys():
    #     op_data[t]['Boiler_Qdot_out_rel'] = P['Heat_demand'][t]/Qnom.value
    #     op_data[t]['Gas_use'] = guess[t]
    #     op_data[t]['Boiler_Qdot_in'] = guess[t]
    # old_vio = P.get_constraint_violations()
    # P.operation = op_data
    # new_vio = P.get_constraint_violations()
    # # Our initialization reduced constraint violation at the initial point
    # assert old_vio.max() > new_vio.max()

    sys_cons = [c_id for c_id in P.constraints if c_id.startswith('Sys')]
    assert len(sys_cons) == 2
    # Checking number of constraints and variables
    cons = ([], [])
    vars = ([], [])
    for comp in G, B, H:
        for i, cd in enumerate(split(comp.constraints, indexed)):
            cons[i].extend(cd)
        vars[0].extend(comp.design_variables)
        vars[1].extend(comp.operational_variables)
    cons[1].extend(sys_cons)
    assert P.num_cons == len(cons[0]) + len(cons[1]) * 4
    assert P.num_vars == len(vars[0]) + len(vars[1]) * 4
    i = 0.08
    n = 10
    af = ((1 + i) ** n * i) / ((1 + i) ** n - 1)  # annuity factor
    assert P.design_objective == B['investment_costs'] * af + B['fixed_costs']
    assert P.operational_objective == G['variable_costs']

    # Post-initialization change of timesteps and scenarios
    timesteps = P.timesteps
    timesteps[4] = 1
    P.timesteps = timesteps
    assert all(P.timesteps == timesteps)
    assert all(P['Gas_price'].value.index == P.index)
    P.scenarios = scenarios[:2]
    assert all(a == b for a, b in zip(P.scenarios, scenarios[:2]))
    assert all(P['Heat_demand'].value.index == P.index)


def test_time_dependent_scenarios():
    """Concept for advanced indexing."""
    import pandas as pd

    p1 = p2 = .5
    # Scenario-dependent timesteps:
    # Here a steady state and a time-series are considered
    scenarios = {'SS1': (p1, 1),
                 'D1': (p2, {'t1': 1, 't2': 3, 't3': 2, 't4': 1})}
    # this input may be disallowed completely or transformed to the following:
    scenarios = {'SS1': (p1, {'t1': 5}),
                 'D1': (p2, {'t1': 1, 't2': 3, 't3': 2, 't4': 1})}

    def get_indices(scenarios):
        for s, (p, timesteps) in scenarios.items():
            try:
                for t, dt in timesteps.items():
                    yield (s, t)
            except AttributeError:
                if not isinstance(timesteps, (int, float)):
                    raise ValueError("The timesteps need to be specified as a "
                                     "single int\\float or a dictionary!")
                yield (s, 't')

    index = pd.MultiIndex.from_tuples(get_indices(scenarios))
    ser = pd.Series(index=index, dtype='O')

    data_variants = [
            {'SS1': 0,
             ('D1', 't1'): 1,
             ('D1', 't2'): 2,
             ('D1', 't3'): 3,
             ('D1', 't4'): 4},
            {('SS1', 't'): 0,  # Note: incorrect time index (should be 't1')
             ('D1', 't1'): 1,
             ('D1', 't2'): 2,
             ('D1', 't3'): 3,
             ('D1', 't4'): 4},
    ]
    for data in data_variants:
        ser[:] = 99
        for k, v in data.items():
            if k not in ser.index:
                print('Oh oh!')
                continue
            ser[k] = v
        print(ser)

    ser[:] = 99
    for k, v in zip(ser.index, data.values()):
        ser[k] = v
    print(ser)

    # time dependent
    # for s in ['SS1', 'D1']:
    #     ser[s] = 1, 2, 3, 4
    # print(ser)

    ser[:] = 99
    # scenario dependent
    value = 0, 1
    for i, (s, (p, timesteps)) in enumerate(scenarios.items()):
        try:
            for t in timesteps:
                ser[s, t] = value[i]
        except TypeError:
            ser[s, t] = value[i]
    print(ser)


def test_remove(clear_components):
    """Test correct behavior of a component's removal from the system."""
    comps = [Component(f'A{i}') for i in range(2)]
    s = System('B')
    for comp in comps:
        d = comp.make_design_variable('d')
        o = comp.make_operational_variable('o')
        comp.add_expression('de', 2 * d)
        comp.add_expression('oe', 2 * o)
        p = comp.make_parameter('p')
        comp.add_le_constraint(d * o, p, 'c')
        comp.add_connectors('C', o)
        comp.declare_state(d, 0.05 * d)
        s.add(comp)
    s.connect('Conn', [getattr(comp, 'C') for comp in comps])

    assert set(s.constraints_dict['B_Conn'].args) == \
        {0, sum(comp['o'] for comp in comps)}
    s.remove(comps[0])
    assert set(s.constraints_dict['B_Conn'].args) == {0, comp['o']}
