"""Code to generate an ALE input file for MAiNGO.

ALE is a library for Algebraic Logical Expressions that can be used to generate
MAiNGO problems from human-readable input files.
This module allows to generate such input files based on a COMANDO problem.
"""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import re

import comando
from comando.utility import split, indexed, StrParser

inf = float('inf')
INF = 1e10


maingo_str_map = {
   '()': lambda expr: f'({expr})',
   'Add': lambda *args: " + ".join(args),
   'Mul': lambda *args: " * ".join(f'{arg}' for arg in args),
   'Pow': lambda base, exponent: f"pow({base}, {exponent})",
   'LessThan': lambda lhs, rhs: f'{lhs} <= {rhs}',
   'GreaterThan': lambda lhs, rhs: f'{lhs} >= {rhs}',
   'Equality': lambda lhs, rhs: f'{lhs} = {rhs}'
}


def maingo_pow_callback(parser, expr, idx):
    """Handle special pow calls in MAiNGO."""
    base, exponent = expr.args
    if base == comando.E:
        return parser.str_map['exp'](*parser.parse_args((exponent, ), idx))
    # if exponent.is_Number:
    #     if exponent < 0:
    #         return parser.str_map['Inv'](base ** -exponent, idx)
    #     from math import log2
    #     log2_exp = log2(float(exponent))
    #     if log2_exp % 1 == 0:  # exponent is a multiple of 2 or 0.5
    #         # NOTE: the case log2_exp == 0 is simplified to 1 by the backend!
    #         if log2_exp > 0:
    #             func = parser.str_map['sqr']
    #             res = func(*parser.parse_args((base, ), idx))
    #             while True:
    #                 log2_exp -= 1
    #                 if log2_exp == 0:
    #                     return res
    #                 res = func(res)
    #         if log2_exp < 0:
    #             func = parser.str_map['sqrt']
    #             res = func(*parser.parse_args((base, ), idx))
    #             while True:
    #                 log2_exp += 1
    #                 if log2_exp == 0:
    #                     return res
    #                 res = func(res)
    return None  # No special case, handle normally


class AleParser(StrParser):  # pylint: disable=too-few-public-methods
    """A class for parsing comando expressions to baron Syntax."""

    def __init__(self, sym_map):
        super().__init__(sym_map, maingo_str_map,
                         pow_callback=maingo_pow_callback)


# TODO: Actively check that no two labels coincide, also could be generalized!
def _normalize(label, replacement='__'):
    return re.sub('[^0-9a-zA-Z]+', replacement, label)


def _ale_var_rep(var, name):
    lb = -INF if var.lb == -inf else var.lb
    ub = INF if var.ub == inf else var.ub
    if var.is_integer:
        if (lb, ub) == (0, 1):
            return f'binary {name};\n'
        return f'integer {name} in [{lb}, {ub}];\n'
    return f'real {name} in [{lb}, {ub}];\n'


def _write_cons(file, cons, parse, suffixes, prefix=''):
    dcons, ocons = split(cons, indexed)
    file.write(f'\n{prefix}constraints:\n')
    for c_id, con in dcons.items():
        file.write(f'{parse(con)} "{c_id}";\n')
    for c_id, con in ocons.items():
        for sfx in suffixes:
            file.write(f'{parse(con, sfx)} "{c_id}_{sfx}";\n')


def _write_ale_file(P, file_name, relaxation_only_constraints,
                    squashing_constraints, cse, outputs,
                    add_intermediates_as_output,
                    sym_map, suffixes, dvs, ovs, parse):
    # Shorthands
    do = P.design_objective
    oo = P.operational_objective
    cons = P.constraints
    ro_cons = relaxation_only_constraints
    sq_cons = squashing_constraints

    with open(file_name, 'w') as file:
        # Variables section
        file.write('definitions:\n\n# Variables\n')
        for dv in dvs:
            file.write(_ale_var_rep(dv, sym_map[dv]))
        for ov in ovs:
            for idx in suffixes.keys():
                file.write(_ale_var_rep(ov[idx], sym_map[ov][idx]))

        if cse:  # Using Common Subexpression Elimination to reduce file size
            from comando.utility import evaluate

            file.write('\n# Intermediate variables:\n')
            orig_exprs = [do, oo, *cons.values()]
            if ro_cons:
                orig_exprs.extend(ro_cons.values())
            if sq_cons:
                orig_exprs.extend(sq_cons.values())
            if outputs:
                orig_exprs.extend(outputs.values())
            reps, exprs = comando.cse(orig_exprs)
            defs = {}
            expr_defs = {}
            for sym, rep in reps:
                e = rep.subs(defs)
                if indexed(e):
                    x = comando.VariableVector(sym.name)
                    x.instantiate(P.index)
                    x_map = sym_map[x] = {}
                    for idx, sfx in suffixes.items():
                        n = x_map[idx] = _normalize(x[idx].name)
                        file.write(f'real {n} := {parse(e, idx)};\n')
                else:
                    x = comando.Variable(sym.name)
                    sym_map[x] = _normalize(x.name)
                    file.write(f'real {x} := {parse(e)};\n')
                x.value = evaluate(e)
                defs[sym] = x
                expr_defs[x] = e

            do = exprs[0].subs(defs)
            oo = exprs[1].subs(defs)
            con_end = 2 + len(cons)
            cons = {con_id: expr.subs(defs) for con_id, expr in
                    zip(P.constraints, exprs[2:con_end])}

            if ro_cons:
                ro_con_end = con_end + len(ro_cons)
                ro_cons = {con_id: expr.subs(defs) for con_id, expr in
                           zip(ro_cons, exprs[con_end:ro_con_end])}
            else:
                ro_con_end = con_end

            if sq_cons:
                sq_con_end = ro_con_end + len(sq_cons)
                sq_cons = {con_id: expr.subs(defs) for con_id, expr in
                           zip(sq_cons, exprs[ro_con_end:sq_con_end])}

            P.cse = [do, oo, cons]
            P.expr_defs = expr_defs
            if outputs:
                outputs = {k: v.subs(defs) for k, v in
                           zip(outputs, exprs[-len(outputs):])}
            else:
                outputs = {}
            if add_intermediates_as_output:
                outputs.update({aux.name: aux for aux in defs.values()})

        file.write('\n# Initial point\n')
        for dv in dvs:
            file.write(f'{sym_map[dv]}.init <- {dv.value};\n')
        for ov in ovs:
            for idx, val in zip(suffixes, ov.value):
                file.write(f'{sym_map[ov][idx]}.init <- {float(val)};\n')

        # Constraints sections
        if cons:
            _write_cons(file, cons, parse, suffixes)
        if ro_cons:
            _write_cons(file, ro_cons, parse, suffixes, 'relaxation only ')
        if sq_cons:
            _write_cons(file, sq_cons, parse, suffixes, 'squashing ')

        # Adding objective
        file.write('\nobjective:\n')
        # file.write(parse(P.objective))
        file.write(parse(do))
        if oo != 0:
            file.write(' + ')
            ts = P.timesteps
            if ts is None:  # Assume scenarios is not None
                file.write(' + '.join(f'{p} * ({parse(oo, s)})'
                                      for s, p in P.scenario_weights.items()))
            elif P.scenarios is None:
                file.write(' + '.join(f'{dt} * ({parse(oo, t)})'
                                      for t, dt in ts.items()))
            else:
                file.write(
                    ' + '.join(
                        f"""{p} * ({' + '.join(f'{dt} * ({parse(oo, (s, t))})'
                                               for t, dt in ts[s].items())})"""
                        for s, p in P.scenario_weights.items()
                    )
                )
        file.write(';\n')
        if outputs:
            from comando.utility import get_index
            file.write('\noutputs:\n')
            for k, v in outputs.items():
                index = get_index(v)
                if index is None:
                    file.write(f'  {parse(v)} "{k}";\n')
                else:
                    for i in index:
                        file.write(f'  {parse(v, i)} "{k}[{i}]";\n')


def write_ale_file(P, file_name, relaxation_only_constraints=None,
                   squashing_constraints=None, cse=True, outputs=None,
                   add_intermediates_as_output=False, reuse=False):
    """Write the problem in ALE syntax to a file or stdout."""
    suffixes = {i: str(i) for i in P.index} \
        if P.scenarios is None or P.timesteps is None \
        else {ii: '_'.join(str(i) for i in ii) for ii in P.index}
    dvs = sorted(P.design_variables, key=lambda x: x.name)
    ovs = sorted(P.operational_variables, key=lambda x: x.name)
    consts, params = split(P.parameters, indexed)

    sym_map = {c: str(c.value) for c in consts}
    for dv in dvs:
        sym_map[dv] = _normalize(dv.name)
    for ov in ovs:
        sym_map[ov] = {idx: _normalize(ov[idx].name)
                       for idx, sfx in suffixes.items()}
    for p in params:
        sym_map[p] = {idx: str(val) for idx, val in zip(suffixes, p.value)}

    parser = AleParser(sym_map)
    parse = parser.cached_parse

    if not reuse:
        _write_ale_file(P, file_name, relaxation_only_constraints,
                        squashing_constraints, cse, outputs,
                        add_intermediates_as_output, sym_map, suffixes, dvs,
                        ovs, parse)
    return sym_map, suffixes


# TODO: Can probably be generalized and moved to utility.py
def write_settings_file(options, settings_name='MAiNGOSettings.txt'):
    """Generate a settings file with the given options."""
    with open(settings_name, 'w') as f:
        f.writelines(f'{option} {value}\n'
                     for option, value in options.items())


def call_maingo(file_name, settings_name=None, silent=False):
    """Call the maingo executable with a problem and possibly settings file."""
    from comando.utility import syscall
    if settings_name:
        return syscall('MAiNGO', file_name, settings_name, silent=silent)
    # Will try to use MAiNGOSettings.txt
    return syscall('MAiNGO', file_name, silent=silent)


def get_results(results_file_name='MAiNGOresult.txt'):
    """Code for parsing MAiNGO results files."""
    import itertools
    import re
    p = re.compile(r" +(\S+) +\S+ +(\S+)")

    val_map = {}
    with open(results_file_name, 'r') as f:
        # Advance two lines
        for line in itertools.islice(f, 2, None):
            try:
                var_name, val = p.search(line).groups()
            except AttributeError:
                break
            val_map[var_name] = val
    return val_map


def solve(P, file_name=None, relaxation_only_constraints=None,
          squashing_constraints=None, silent=False, cse=True, outputs=None,
          add_intermediates_as_output=False, reuse=None, **options):
    """Solve poblem P using MAiNGO."""
    from comando.utility import canonical_file_name, check_reuse_or_overwrite

    base_name, file_name = canonical_file_name(P.name, '.ale', file_name)
    check_reuse_or_overwrite(file_name, reuse)

    sym_map, suffixes = write_ale_file(P, file_name,
                                       relaxation_only_constraints,
                                       squashing_constraints, cse, outputs,
                                       add_intermediates_as_output, reuse)
    if options:
        settings_name = f'{base_name}_Settings.txt'
        write_settings_file(options, settings_name)
    else:
        settings_name = None
    ret = call_maingo(file_name, settings_name, silent)
    if ret != 0:
        raise RuntimeError('ERROR: Solver returned nonzero exit status!')

    with open('MAiNGO.log', 'r') as f:
        content = f.read()
    INF_MSG = "Problem is infeasible!"
    if INF_MSG in content:
        print('Log says:', INF_MSG)
    else:
        vals = get_results()

        for dv in P.design_variables:
            dv.value = vals[sym_map[dv]]
        for ov in P.operational_variables:
            for idx, name in sym_map[ov].items():
                ov[idx] = vals[name]
    return ret
