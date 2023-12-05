"""A simple interface to Pyomo."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import pyomo.core
from pyomo.core import ConcreteModel, Set, Var, Param, Objective, \
    Constraint, Any
from pyomo.core import Reals as pyomo_Reals
from pyomo.core import Integers as pyomo_Integers
from pyomo.core import Binary as pyomo_Binary
from pyomo.environ import SolverManagerFactory, SolverFactory

import comando
from comando.utility import indexed, parse, split


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
    pyomo_var.setlb(comando_var.lb)
    pyomo_var.setub(comando_var.ub)
    value = comando_var.value if index is None else comando_var.value[index]
    try:
        pyomo_var.set_value(value)
    except ValueError as err:
        if 'is not in domain' in str(err):
            pyomo_var.set_value(round(value))


def to_pyomo(problem):
    """Parse the optimization problem in energy_system to pyomo.

    Arguments
    ---------
    problem : comando.Problem

    Returns
    -------
    m : pyomo.ConcreteModel
    """
    P = problem

    m = ConcreteModel()
    m.idx = Set(initialize=P.index)

    pars = P.parameters
    dp, op = split(pars, indexed)
    dv = P.design_variables
    ov = P.operational_variables
    # TODO: it seems like better care when creating P should make this obsolete
    fake_dv, dv = split(dv, lambda dvi: dvi.parent is None)
    ov.update({fdv.parent for fdv in fake_dv})

    m.dp = Param((p.name for p in dp), mutable=True)
    m.op = Param((p.name for p in op), m.idx, mutable=True)
    m.x = Var(v.name for v in dv)
    m.y = Var((v.name for v in ov), m.idx)
    m.sym_map = {}
    for p in dp:
        ppar = m.dp[p.name]
        ppar.value = p.value
        m.sym_map[p] = ppar
    for pv in op:
        ppars = {}
        for i, p in pv.expansion.items():
            ppar = m.op[pv.name, i]
            ppar.value = p.value
            m.sym_map[p] = ppar
            ppars[i] = ppar
        m.sym_map[pv] = ppars
    for v in dv:
        pvar = m.x[v.name]
        _translate(v, pvar)
        m.sym_map[v] = pvar
    for vv in ov:
        pvars = {}
        for i, v in vv.expansion.items():
            pvar = m.y[vv.name, i]
            _translate(v, pvar)  # ov is VariableVector
            m.sym_map[v] = pvar
            pvars[i] = pvar
        m.sym_map[vv] = pvars

    # Raise an error if any of the parameters has no value
    miss = '\n\t'.join(p.name for p in P.parameters if m.sym_map[p] is None)
    if miss:
        raise ValueError(f"Lacking data for parameter(s):\n\t{miss}")

    def parse2pyomo(expr, idx=None):
        """Parse expressions from comando to pyomo."""
        return parse(expr, m.sym_map, pyomo_op_map, idx=idx)
    m.parse2pyomo = parse2pyomo

    # Adding objective
    # TODO: The goal should be:
    #       obj = parse2pyomo(P.objective)
    #       but currently parsing replaces all parameters with their values!
    #       Consequently updated values for parameter in the objective function
    #       would not be reflected!
    do = parse2pyomo(P.design_objective)
    oo = 0
    ts = P.timesteps
    if ts is None:
        oo = sum(p * parse2pyomo(P.operational_objective, s)
                 for s, p in P.scenario_weights.items())
    elif P.scenarios is None:
        oo = sum(dt * parse2pyomo(P.operational_objective, t)
                 for t, dt in ts.items())
    else:
        oo = sum(p * dt * parse2pyomo(P.operational_objective, (s, t))
                 for s, p in P.scenario_weights.items()
                 for t, dt in ts[s].items())

    m.objective = Objective(expr=do + oo)

    # Add constraints defined by components and their connections
    constraints = P.constraints
    m.constraints = Constraint(Any)

    def add_constraint(name, con, idx=None):
        """Add constraints to `m` if satisfaction cannot be pre-determined."""
        # NOTE: We cannot use implicit conversion to bool because pyomo
        #       constraint objects will use the current variable values.
        if con == True:  # constraint can never be violated -> ignore it
            return
        if con == False:
            raise comando.ImpossibleConstraintException
        if idx is None:
            m.constraints[name] = parse2pyomo(con)
        else:
            m.constraints[f'{name}_{idx}'] = parse2pyomo(con, idx)

    for con_id, con in constraints.items():
        if indexed(con):
            for idx in m.idx:
                try:
                    add_constraint(con_id, con, idx)
                except comando.ImpossibleConstraintException:
                    # FIXME: Should be a proper warning
                    print(f'WARNING: Dropping constraint {con_id}[{idx}], as '
                          'it will always be violated!')
        else:
            try:
                add_constraint(con_id, con)
            except comando.ImpossibleConstraintException:
                # FIXME: Should be a proper warning
                print(f'WARNING: Dropping constraint {con_id}[{idx}], as '
                      'it will always be violated!')

    if P.states:
        print('WARNING: There are differential constraints that will be '
              'discretized using implicit Euler discretization, if you want '
              'an advanced time discretization you must reformulate your '
              'problem prior to passing it to the Pyomo interface!')
        skipped_indices = set()
        Delta_t = P.Delta_t
        if Delta_t not in op:
            # P.parameters.add(Delta_t)
            m.Delta_t = Param(m.idx, mutable=True)
            ppars = {}
            for i, p in Delta_t.expansion.items():
                ppar = m.Delta_t[i]
                ppar.value = p.value
                m.sym_map[p] = ppar
                ppars[i] = ppar
            m.sym_map[Delta_t] = ppars

        from comando.utility import get_previous

    for state, state_data in P.states.items():
        initial_state, derivative, state_expression = state_data

        prev = get_previous(state, initial_state)
        prev_map = {}
        for i, v in prev.expansion.items():
            try:
                prev_map[i] = m.sym_map[v]
            except KeyError:
                prev_map[i] = m.sym_map[initial_state]
        m.sym_map[prev] = prev_map
        con = comando.Eq((state - prev) / Delta_t, derivative)
        for i in m.idx:
            if i in skipped_indices:
                continue
            if ts[i] == 0:
                print(f'WARNING: index {i} corresponds to a length of 0... '
                      'skipping!')
                skipped_indices.add(i)
                continue
            add_constraint(state.name, con, i)

    def writeback():
        for v in dv:
            v.value = m.sym_map[v]()
        for vv in ov:
            vv.value = {idx: pvar() for idx, pvar in m.sym_map[vv].items()}
    m.writeback = writeback

    def update():
        """Update the variables and parameters withdata from COMANDO.

        When doing multiple solves or manually manipulating the COMANDO problem
        this method can be used to update the corresponding Pyomo objects.

        Note
        ----
        The structure of the COMANDO problem is assumed to remain the same.
        If the structure changes a new Pyomo model should be created!
        """
        for p in P.parameters:
            try:  # operational parameter
                for par, ppar in zip(p, m.sym_map[p].values()):
                    ppar.value = par.value
            except AttributeError:  # design parameter
                m.sym_map[p].value = p.value
        for dv in P.design_variables:
            _translate(dv, m.sym_map[dv])
        for ov in P.operational_variables:
            for var, pvar in zip(ov, m.sym_map[ov].values()):
                _translate(var, pvar)
    m.update = update

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
        import os
        keepfiles = kwargs.get('keepfiles', False)
        if not keepfiles:
            dir_contents = {*os.listdir('.')}
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

        # Load results back into COMODO Variables
        if writeback:
            m.writeback()
        if not keepfiles:
            for f in {*os.listdir('.')} - dir_contents:
                os.remove(os.path.join('.', f))
        return results

    m.solve = solve
    return m
