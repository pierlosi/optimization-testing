"""A simple interface to Pyomo.DAE."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Florian Joseph Baader, Marco Langiu
import pyomo.core
from pyomo.core import ConcreteModel, Var, Objective, ConstraintList
from pyomo.core import Reals as pyomo_Reals
from pyomo.core import Integers as pyomo_Integers
from pyomo.core import Binary as pyomo_Binary
from pyomo.environ import SolverManagerFactory, SolverFactory, \
    TransformationFactory
from pyomo.dae import ContinuousSet, DerivativeVar

import comando
import sympy
from comando.utility import parse

import copy
import numpy as np
from pandas import Series


# Explicitly registering numeric types that may occur in our expressions
pyomo.core.expr.numvalue.RegisterNumericType(comando.sympy.Float)
pyomo.core.expr.numvalue.RegisterIntegerType(comando.sympy.Integer)


pyomo_op_map = comando.base_op_map.copy()
for name in ['exp', 'log', 'sqrt', 'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh',
             'asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh']:
    pyomo_op_map[name] = getattr(pyomo.core, name)


def _translate(comando_var, pyomo_var, index=None):
    if comando_var.is_binary:
        pyomo_var.domain = pyomo_Binary
    elif comando_var.is_integer:
        pyomo_var.domain = pyomo_Integers
    else:
        pyomo_var.domain = pyomo_Reals
    pyomo_var.setlb(comando_var[index].lb)
    pyomo_var.setub(comando_var[index].ub)
    value = comando_var.value if index is None else comando_var.value[index]
    try:
        pyomo_var.set_value(value)
    except ValueError as err:
        if 'is not in domain' in str(err):
            pyomo_var.set_value(round(value))


def to_pyomo_dae(problem, **dynamic_options):
    """Parse the optimization problem in energy_system to Pyomo.DAE.

    Pyomo.DAE allows to solve dynamic problem formulations with advanced time
    discretization approaches.

    Note
    ----
    The objective value resulting form optimization of the returned pyomo model
    - which can be computed with `m.objective()` - will generally NOT be the
    same as the objective value that can be computed with the optimal solution
    using the original COMANDO Problem!

    Arguments
    ---------
    problem : comando.Problem
    dynamic_options : dict
        dict with 3 entries:

        1. 'controls': iterable of comando.Variable objects acting as controls
        2. 'control_factor': number of finite elements for which controls are
           set picewise constant
        3. 'pyomo_dae_options': options for pyomo.dae

    Returns
    -------
    m : pyomo.ConcreteModel
    """
    P = problem
    if P.timesteps is None:
        raise RuntimeError('Problem without timesteps cannot be handled via '
                           'the Pyomo.DAE interface!')
    if P.scenarios is not None:
        raise NotImplementedError('Problem with scenarios can currently not '
                                  'be handled by the Pyomo.DAE interface!')
    if 'pyomo_dae_options' not in dynamic_options:
        raise ValueError("dynamic_options must contain 'pyomo_dae_options'!")
    pyomo_dae_options = dynamic_options['pyomo_dae_options']

    if P.scenarios:
        raise NotImplementedError('scenarios in combination with pyomo.dae '
                                  'have not been implemented yet!')
    # so far only option for dynamic opt #!!!!!!!!!!!!Todo
    # assert P.objective_type == 'Mayer term objective'
    # we have to change this, if we want to do integrated design and operation

    pars = P.parameters
    dv = P.design_variables
    ov = copy.copy(P.operational_variables)

    # define set of state variables sv
    # derivatives of states are stored in a dictionary so they can later be
    #  mapped to m.der['statename']
    sv = set()
    derv = {}
    for state in P.states:
        sv.add(state)
        derivative = P.states[state][1]
        derv.update({derivative: state})
        ov.remove(state)
        ov.remove(derivative)

    m = ConcreteModel()
    t_final = P.T
    m.t = ContinuousSet(bounds=(0, t_final))
    t_0 = 0

    if dynamic_options['controls']:
        cv = dynamic_options['controls']
        for c in cv:
            ov.remove(c)
        if dynamic_options['control_factor']:
            control_factor = dynamic_options['control_factor']
        else:
            control_factor = 1
        nfe_controls = pyomo_dae_options['nfe'] / control_factor

        t_control = np.linspace(t_0 + (t_final - t_0) / nfe_controls, t_final,
                                int(nfe_controls))
        P.t_control = t_control
        m.c = Var((v.name for v in cv), m.t)
    else:
        cv = set()

    # variable costs have to be determined via integration if discretization
    # is flexible
    IntegratedVariableObjective = \
        comando.VariableVector('IntegratedVariableObjective')
    sv.add(IntegratedVariableObjective)

    m.x = Var(v.name for v in dv)
    m.y = Var((v.name for v in ov), m.t)
    m.s = Var((v.name for v in sv), m.t)
    m.der = DerivativeVar(m.s, withrespectto=m.t)
    assert (pyomo_dae_options['method'] == 'dae.finite_difference'
            or pyomo_dae_options['method'] == 'dae.collocation')
    if pyomo_dae_options['method'] == 'dae.finite_difference':
        assert (pyomo_dae_options['scheme'] == 'BACKWARD'
                or pyomo_dae_options['scheme'] == 'CENTRAL'
                or pyomo_dae_options['scheme'] == 'FORWARD')
        assert isinstance(pyomo_dae_options['nfe'], int)
        discretizer = TransformationFactory(pyomo_dae_options['method'])
        discretizer.apply_to(m, wrt=m.t, scheme=pyomo_dae_options['scheme'],
                             nfe=pyomo_dae_options['nfe'])
    elif pyomo_dae_options['method'] == 'dae.collocation':
        assert isinstance(pyomo_dae_options['nfe'], int)
        assert isinstance(pyomo_dae_options['ncp'], int)
        # TODO: 'LAGRANGE-RADAU' 'LAGRANGE-LEGENDRE'
        discretizer = TransformationFactory(pyomo_dae_options['method'])
        discretizer.apply_to(m, wrt=m.t, nfe=pyomo_dae_options['nfe'],
                             ncp=pyomo_dae_options['ncp'])

    m.comando_time_map = {}
    idx = 0
    for pyomo_t in m.t:
        if pyomo_t > P.index[idx]:
            idx += 1
        m.comando_time_map[pyomo_t] = P.index[idx]

    m.control_time_map = {}
    idx = 0
    for pyomo_t in m.t:
        if pyomo_t > P.t_control[idx]:
            idx += 1
        m.control_time_map[pyomo_t] = P.t_control[idx]

    m.sym_map = {}
    for p_comando_t in pars:
        if p_comando_t.indexed:
            p_pyomo_t = Series(index=m.t, name=p_comando_t.name, dtype=float)
            for t_pyomo in m.t:
                p_pyomo_t[t_pyomo] = p_comando_t.value.get(
                    m.comando_time_map[t_pyomo]
                )
            m.sym_map[p_comando_t] = p_pyomo_t
        else:
            m.sym_map[p_comando_t] = p_comando_t.value
    for v in dv:
        pvar = m.x[v.name]
        _translate(v, pvar)
        m.sym_map[v] = pvar
    for v in ov:
        pvars = {}
        for t_pyomo in m.t:
            pvar = m.y[v.name, t_pyomo]
            _translate(v, pvar, index=m.comando_time_map[t_pyomo])
            pvars[t_pyomo] = pvar
        m.sym_map[v] = pvars
    for v in sv:
        if v.name == 'IntegratedVariableObjective':
            pvars = {}
            for t in m.t:
                pvar = m.s[v.name, t]
                pvars[t] = pvar
            m.sym_map[v] = pvars
        else:
            pvars = {}
            for t_pyomo in m.t:
                pvar = m.s[v.name, t_pyomo]
                _translate(v, pvar, index=m.comando_time_map[t_pyomo])
                pvars[t_pyomo] = pvar
            m.sym_map[v] = pvars
    for v in derv:
        pvars = {}
        for t_pyomo in m.t:
            pvar = m.der[derv[v].name, t_pyomo]
            _translate(v, pvar, index=m.comando_time_map[t_pyomo])
            pvars[t_pyomo] = pvar
        m.sym_map[v] = pvars

    for v in cv:
        pvars = {}
        for t_pyomo in m.t:
            pvar = m.c[v.name, m.control_time_map[t_pyomo]]
            _translate(v, pvar, index=m.comando_time_map[t_pyomo])
            pvars[t_pyomo] = pvar
        m.sym_map[v] = pvars

    # Raise an error if any of the parameters has no value
    miss = '\n\t'.join(p.name for p in P.parameters if m.sym_map[p] is None)
    if miss:
        raise ValueError(f"Lacking data for parameter(s):\n\t{miss}")

    def parse2pyomo(expr, idx=None):
        """Parse expressions from comando to pyomo expanding indexed Symbols."""
        return parse(expr, m.sym_map, pyomo_op_map, idx=idx)

    # for later use
    m.parse2pyomo = parse2pyomo

    m.constraints = ConstraintList()

    for conKey in P.constraints:
        if any(v in P.constraints[conKey].free_symbols
               for v in ov | sv | cv | derv.keys()):  # indexed constraints
            for t in m.t:
                m.constraints.add(parse2pyomo(P.constraints[conKey], t))
        else:  # unindexed constraints
            m.constraints.add(parse2pyomo(P.constraints[conKey]))

    # define mayer term objective
    phi_dot = P.operational_objective
    d = parse2pyomo(P.design_objective)
    Phi = parse2pyomo(IntegratedVariableObjective, m.t.last())
    m.objective = Objective(expr=d+Phi)

    m.constraints.add(m.s[IntegratedVariableObjective.name, t_0] == 0)
    for t in m.t:
        m.constraints.add(m.der[IntegratedVariableObjective.name, t]
                          == parse2pyomo(phi_dot, t))

    for state in P.states:
        initial_state = P.states[state][0]
        if initial_state is not None:
            m.constraints.add(parse2pyomo(sympy.Eq(state, initial_state), t_0))

    def solve(solver='cplex', remote=False, writeback=True, **kwargs):
        """Solve pyomo model `m`.

        The solution is calculated using the specified solver. The `remote`
        flag determines whether this is done locally or remotely on the 'NEOS'
        server.
        If no removal of the generated solver input, logs and output files is
        desired, the `keepfiles` flag can be set to `True`.
        Any additional pyomo options such as `tee=True` or
        `symbolic_solver_labels=True` can be specified via the other kwargs.
        """
        # import os
        # keepfiles = False
        # if not kwargs.get('keepfiles', False):
        #     keepfiles = True
        #     dir_contents = {*os.listdir('.')}
        if remote:
            with SolverManagerFactory('neos') as manager:
                opt = SolverFactory(solver)
                for option, val in kwargs.pop('options', {}).items():
                    opt.options[option] = val
                results = manager.solve(m, opt=opt, **kwargs)
        else:
            solver = SolverFactory(solver, options=kwargs.pop('options', {}),
                                   executable=kwargs.pop('executable', None))
            results = solver.solve(m, **kwargs)
        results.write()

        # write back
        time = [t for t in m.t]
        new_timesteps = {}
        new_timesteps[time[0]] = 0
        for i in range(1, len(time)):
            new_timesteps[time[i]] = time[i] - time[i - 1]
        from pandas import Series
        new_timesteps = Series(new_timesteps)
        P.timesteps = new_timesteps

        for v in dv:
            v.value = m.sym_map[v]()
        for vv in ov | sv | cv | derv.keys():
            if not vv.name == 'IntegratedVariableObjective':
                vv.value = {idx: pvar() for idx, pvar in m.sym_map[vv].items()}
        for p in pars:
            if p.indexed:
                for t in m.t:
                    p[t].value = m.sym_map[p][t]

    m.solve = solve
    return m
