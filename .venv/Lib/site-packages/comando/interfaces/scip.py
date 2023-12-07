"""COMANDO-SCIP interface."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Alexander Holtwerth, Marco Langiu
from pyscipopt import Model, quicksum, quickprod, exp, log, sqrt

from comando import base_op_map
from comando.utility import parse, indexed


def scip_sum(*args):
    """Use SCIP's quicksum for addition."""
    return quicksum(args)


def scip_prod(*args):
    """Use SCIP's quickprod for multiplication."""
    return quickprod(args)


scip_op_map = base_op_map.copy()
scip_op_map['Add'] = scip_sum
scip_op_map['Mul'] = scip_prod
scip_op_map['exp'] = exp
scip_op_map['log'] = log
scip_op_map['sqrt'] = sqrt


def _get_vtype(var):
    if var.domain.value == 1:  # REAL -> 'CONTINUOUS' (OR SIMPLY 'C')
        return 'C'
    return var.domain.name  # 'INTEGER' AND 'BINARY' are accepted by SCIP


class ScipProblem(Model):
    """An extension of the SCIP model for COMANDO."""

    def __init__(self, P):
        super().__init__(P.name)
        self.P = P
        index = P.index
        dvs = P.design_variables
        pars = P.parameters
        ovs = P.operational_variables

        self.sym_map = sym_map = {}
        for p in pars:
            sym_map[p] = p.value
        for v in dvs:
            sym_map[v] = self.addVar(v.name, _get_vtype(v), *v.bounds, v.value)
        for vv in ovs:
            svars = {}
            for idx, v in vv.items:
                svar = self.addVar(v.name, _get_vtype(v), *v.bounds, v.value)
                sym_map[v] = svars[idx] = svar
            sym_map[vv] = svars

        # Raise an error if any of the parameters has no value
        miss = '\n\t'.join(p.name for p in P.parameters if sym_map[p] is None)
        if miss:
            raise ValueError(f"Lacking data for parameter(s):\n\t{miss}")

        def parse2scip(expr, idx=None):
            """Parse expressions from COMANDO to SCIP."""
            return parse(expr, sym_map, scip_op_map, idx=idx)

        # Adding objective
        do = parse2scip(P.design_objective)
        ts = P.timesteps
        if ts is None:
            oo = quicksum(p * parse2scip(P.operational_objective, s)
                          for s, p in P.scenario_weights.items())
        elif P.scenarios is None:
            oo = quicksum(dt * parse2scip(P.operational_objective, t)
                          for t, dt in ts.items())
        else:
            oo = quicksum(
                p * quicksum(dt * parse2scip(P.operational_objective, (s, t))
                             for t, dt in ts[s].items())
                for s, p in P.scenario_weights.items()
            )

        # NOTE: SCIP doesn't support using nonlinear objectives directly!
        self.obj_var = self.addVar('obj_var', lb=None, ub=None)
        self.setObjective(self.obj_var)
        self.addCons(do + oo <= self.obj_var)

        # Add constraints defined by components and their connections
        constraints = P.constraints

        for con_id, con in constraints.items():
            if indexed(con):
                for idx in index:
                    self.addCons(parse2scip(con, idx), con_id)
            else:
                self.addConstr(parse2scip(con), con_id)

    def solve(self, **options):
        """Solve the SCIP model.

        Arguments
        ---------
        options : dict
            SCIP options (also called 'Parameters') with '/' replaced by '_',
            i.e., a value for the option `limits/time` is given with the key
            `'limits/time'=...`. For a full list of available options see
            https://www.scipopt.org/doc-7.0.2/html/PARAMETERS.php
            or the equivalent for the SCIP version you have installed.
        """
        self.redirectOutput()  # To make output work in interactive sessions
        for option, value in options.items():
            self.setParam(option.replace('_', '/'), value)
        self.optimize()
        self.write_back()

    def write_back(self, P=None):
        """Write back results from gurobi model to a COMANDO problem."""
        if P is None:
            P = self.P
        sol = self.getBestSol()
        for v in P.design_variables:
            v.value = sol[self.sym_map[v]]
        for vv in P.operational_variables:
            for i, v in vv.items:
                v.value = sol[self.sym_map[vv][i]]