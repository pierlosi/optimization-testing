"""Simple user defined algorithms for the IES case study."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu, David Shu
from pandas import DataFrame, Index, MultiIndex

import comando
from comando.interfaces.pyomo import to_pyomo
from comando.utility import silence


SEED = 123  # Random seed for repeatability, unfortunately only for GUROBI
TOL = 0.01  # Default solver tolearance


def make_multiobjective(ES, timesteps, scenarios, data, **objectives):
    """Create a problem with a weighting between two objectives.

    Arguments
    ---------
    ES : System
    timesteps : 2-tuple / Mapping / Series
    scenarios : sequence / Mapping / Series
    data : Mapping / DataFrame
        Parameter data corresponding to the chosen scenarios and time steps
    objectives : names and 2-tuple
        2-tuples containing the design objective and operational objective
    """
    if len(objectives) != 2:
        raise NotImplementedError('Currently only bi-objective problems can '
                                  'be solved!')
    obj_names = [*objectives]
    obj_vars = [comando.Variable(obj_name) for obj_name in obj_names]
    obj_switch = comando.Parameter('obj_switch')  # 0 or 1
    obj = (1 - obj_switch) * obj_vars[0] + obj_switch * obj_vars[1]

    # NOTE: The order of the next four lines is important:
    #       The problem needs to be created, symbols that are in the objectives
    #       but not in any other part of the problem need to be added
    #       explicitly, available data needs to be set so time-variable
    #       parameters are expanded appropriately and finally the objective
    #       expressions can be formed.
    #       If the order is changed, some parameters may still lack values and
    #       thus be scalar when forming the expressions. Since expressions do
    #       not update automatically, expanding some of these parameters later
    #       on would wrongly make the expressions operational, as they would
    #       now contain parameter vectors.
    P = ES.create_problem(design_objective=obj,
                          operational_objective=0,
                          timesteps=timesteps,
                          scenarios=scenarios,
                          name=f"min {'|'.join(objectives)}")

    # NOTE: We need sympify to ensure things like 0 are converted to exprs!
    P.add_symbols(set().union(*(comando.sympify(obj_part).free_symbols
                                for obj_parts in objectives.values()
                                for obj_part in obj_parts)))
    P.data = data
    P.obj_names = obj_names
    P.obj_vars = obj_vars
    P.obj_switch = obj_switch
    P.obj_exprs = [dv + P.weighted_sum(ov) for dv, ov in objectives.values()]

    # TODO: Create a wrapper for adding constraints that takes care of this
    def add_constraint(P, con, con_id):
        """Manually add a constraint to Problem P."""
        P.add_symbols(con.free_symbols)
        P.constraints[con_id] = con

    for obj_name, obj_var, obj_expr in zip(obj_names, obj_vars, P.obj_exprs):
        add_constraint(P, comando.Eq(obj_var, obj_expr), obj_name)

    print(f'Problem {P.name} has {P.num_vars} variables and {P.num_cons} '
          'constraints!')

    return P


def _setup_solve(P, linearization, solver, options, silent):
    """Return the Problem to be solved, a Pyomo model and a solve method."""
    if options is None:
        options = {}
    if linearization is not None:
        from comando.linearization import linearize
        _P = P_lin = linearize(P, *linearization)
        print(f'Linearization of problem {P.name} has {P_lin.num_vars} '
              f'variables and {P_lin.num_cons} constraints!')
        m = to_pyomo(P_lin)
        if solver is None:
            solver = 'gurobi'
            options = {**dict(mipgap=TOL, timelimit=86400, threads=1,
                              seed=SEED),
                       **options}
        if silent:
            def solve(logname):
                parts = logname.split('.')
                parts[-2] += ('_lin')
                logname = '.'.join(parts)
                parts
                with silence():
                    res = m.solve(solver, options=options,
                                  warmstart=True, tee=True,
                                  logfile=logname, keepfiles=True)
                return res
        else:
            def solve(logname):
                parts = logname.split('.')
                parts[-2] += ('_lin')
                logname = '.'.join(parts)
                parts
                return m.solve(solver, options=options, warmstart=True,
                               tee=True, logfile=logname,
                               keepfiles=True)
    else:
        _P = P
        m = to_pyomo(P)
        if solver is None:
            solver = 'baron'
            options = {**dict(EpsR=TOL, MaxTime=3600), **options}
        if silent:
            def solve(logname):
                with silence():
                    res = m.solve(solver, options=options,
                                  logfile=logname, keepfiles=True)
                return res
        else:
            def solve(logname):
                return m.solve(solver, options=options, tee=True,
                               logfile=logname, keepfiles=True)
    return _P, m, solve


def solve_multiobjective(P, num_sol=8, linearization=None, solver=None,
                         options=None, silent=True, callback=None):
    """Solve a multi-objective problem based on the given energy system.

    Solve the given multi-objective problem via the epsilon-constraint method.
    In each iteration one of the objectives selected as the primary objective.
    Then minimization of this primary objective is performed with limitations
    imposed on the secondary one and afterwards a correction is performed by
    minimizing the secondary objective while limiting the primary one.
    After establishing the best and worst values achievable for both
    objectives, multiple iterations are performed with the first objective as
    the primary one, iteratively reducing the upper bound of the second one.

    Arguments
    ---------
    P : Problem
        a multi-objective problem
    num_sol : int
        number of points on the pareto curve to be generated
    linearization : None / 2-tuple
        Options for linearization (n_bp, method)
    solver : str
        name of solver executable (must be on path and compatible with Pyomo)
    options : dict
        keys and values of valid options for the chosen solver
    callback : Callable
        callback that is executed before each iteration and is handed P and it

    Returns
    -------
    obj_vals : DataFrame
        Values of the objectives for each iteration
    dvs : DataFrame
        Design variable values of the original or linearized problem for each
        iteration
    ovs : DataFrame
        Operational variable values of the original or linearized problem for
        each iteration
    """
    _P, m, solve = _setup_solve(P, linearization, solver, options, silent)

    lb = [None, None]
    ub = [None, None]

    index = Index(range(num_sol), name='iteration')
    obj_vals = DataFrame(index=index, columns=P.obj_names, dtype='d')
    dvs = DataFrame(index=MultiIndex.from_product([index, _P.design.index]),
                    columns=['value'])
    ovs = DataFrame(index=MultiIndex.from_product([index, _P.operation.index]),
                    columns=_P.index)

    def run(it, obj_num, callback=None):
        """Do an iteration minimizing the objective with the given number."""
        print('\n\n' + '#' * 79)
        print(f'Iteration {it}:')
        print(f'\tminimizing {P.obj_names[obj_num]}...')

        if callback:
            callback(P, it, lb, ub)

        P.obj_switch.value = obj_num

        # Remove constraints from other objective
        for obj_var, ubi in zip(P.obj_vars, ub):
            obj_var.ub = ubi
        m.update()
        for n in P.obj_names:
            print(f'\t\t{m.x[n].lb} ≤ {n} ≤ {m.x[n].ub}')

        res = solve(f'{it}_min_{P.obj_names[obj_num]}_{SEED}_{TOL}.log')
        print(f"\t...took {res['Solver'][0]['Time']}s of CPU time")
        if P is not _P:
            for n in P.obj_names:
                print(f'\t\tlinearized {n} = {m.x[n].value}')
        for n, v in zip(P.obj_names, P.obj_exprs):
            print(f'\t\t{n} = {v.value}')

        lb[obj_num] = P.obj_vars[obj_num].lb = \
            max(P.obj_vars[obj_num].lb, res['Problem'][0]['Lower bound'])

        other = 1 - obj_num
        print(f'\n\tcorrecting {P.obj_names[other]}...')

        P.obj_switch.value = other

        # Use current value as upper bound
        for obj_var, expr in zip(P.obj_vars, P.obj_exprs):
            val = expr.value
            try:
                obj_var.ub = val
            except ValueError:
                print(f'\t\tMinimization of {obj_var.name} was exact (LB=UB), '
                      'but the lower bound appears to be incorrect!\n\t\t'
                      f'Adjusting from {obj_var.lb} to {val}')
                obj_var.lb = val
                obj_var.ub = val

        m.update()
        for n in P.obj_names:
            print(f'\t\t{m.x[n].lb} ≤ {n} ≤ {m.x[n].ub}')

        res = solve(f'{it}_correct_{P.obj_names[other]}_{SEED}_{TOL}.log')

        print(f"\t...took {res['Solver'][0]['Time']}s of CPU time")
        if P is not _P:
            for n in P.obj_names:
                print(f'\t\tlinearized {n} = {m.x[n].value}')
        for n, v in zip(P.obj_names, P.obj_exprs):
            print(f'\t\t{n} = {v.value}')

        ub[other] = P.obj_exprs[other].value

        P.obj_switch.value = obj_num
        m.update()
        *obj_vals.loc[it], dvs.loc[it], ovs.loc[it] = \
            *(e.value for e in P.obj_exprs), _P.design.values, \
            _P.operation.values

    # use first objective (best/worst solution for first/second objective)
    run(0, 0, callback)

    # use second objective (worst/best solution for first/second objective)
    run(num_sol-1, 1, callback)

    # Establish steps for second objective
    ub_max = ub[1]
    step = (ub_max - P.obj_exprs[1].value) / (num_sol - 1)
    for iter in range(1, num_sol - 1):  # do num_sol-2 iters (already did 2)!
        ub[1] = ub_max - step * iter
        run(iter, 0, callback)

    return obj_vals, dvs, ovs


def get_vio(P):
    """Get constraint id and amount of largest constraint violation."""
    vio = P.get_constraint_violations()
    max_vio = vio.max()
    max_vio_con = vio[vio == max_vio]
    if any(max_vio_con):
        return next(iter(max_vio_con.items()))
    return None, max_vio
