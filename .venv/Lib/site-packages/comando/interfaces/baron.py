"""Input file generation for BARON."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
from functools import partial
from math import isinf
import sys

import comando
from comando import Eq, Le

from comando.utility import split, indexed, get_index, get_vars, StrParser

# In some newer macOS versions 'System Integrity Protection' appears to prevent
# interactive Python sessions to access certain environment variables like
# DYLD_LIBRARY_PATH which is where shared libraries such as for CPLEX are
# stored by default.
# As Baron relies on DYLD_LIBRARY_PATH to be set correclty, we do this here
if sys.platform == 'darwin':
    # NOTE: Make sure the shared library for CPLEX 12.10 can be found!
    import os
    import subprocess

    def search_cplex_library(major='12', minor='10', patch='0', *ignored):
        """Locate the shared library file for the given CPLEX version."""
        cplex_library = f"libcplex{major}{minor}{patch}.dylib"
        print(f'searching for {cplex_library}...')
        if ignored:
            print('ignoring:', ignored)
        res = subprocess.check_output(f'locate {cplex_library}', shell=True)
        try:
            loc = res.decode().split()[0]
            print(f'...found at {loc}!')
            print('Temporarily adding to DYLD_LIBRARY_PATH!')
            os.environ['DYLD_LIBRARY_PATH'] = os.path.dirname(loc)
            return True
        except IndexError:
            print('... nothing found!')
            return False

    if not search_cplex_library():
        prompt = ('Please enter the CPLEX version you have installed\n'
                  'in the form major.minor.patch, e.g., 12.10.0 or hit\n'
                  "enter if you do not have CPLEX or don't want to use it\n")
        while True:
            version = input(prompt).split('.')
            if version[0] == '' or search_cplex_library(*version):
                break

baron_str_map = {'Neg': lambda arg: f'(-{arg})',  # unary negation requires ()!
                 'LessThan': lambda lts_minus_gts, ZERO:
                     f'{lts_minus_gts} <= 0',
                 'GreaterThan': lambda lts_minus_gts, ZERO:
                     f'{lts_minus_gts} <= 0',
                 # NOTE: symengine stores both a == 0 and 0 == a as the latter!
                 'Equality': lambda ZERO, lhs_minus_rhs:
                     f'{lhs_minus_rhs} == 0',  # BARON needs a == 0
                 'Pow': lambda base, exponent: f'{base} ^ {exponent}',
                 'tanh': lambda arg: f'(1 - 2/(exp(2 * ({arg})) + 1))'}


def baron_pow_callback(parser, expr, idx):
    """Handle special pow calls in BARON."""
    base, exponent = expr.args
    if base == comando.E:
        return parser.str_map['exp'](*parser.parse_args((exponent, ), idx))
    # NOTE: BARON does not allow x^y, where x and y are both variables.
    #       It is permissible to have either x or y as a variable in this
    #       case but not both. The following reformulation can be used
    #       around this: x^y = exp(y * log(x))
    if all(get_vars(arg) for arg in expr.args):
        arg = base * comando.log(exponent)
        return parser.str_map['exp'](*parser.parse_args((arg, ), idx))
    return None  # No special case, handle normally


class BaronParser(StrParser):  # pylint: disable=too-few-public-methods
    """A class for parsing comando expressions to baron Syntax."""

    def __init__(self, sym_map):
        super().__init__(sym_map, baron_str_map,
                         pow_callback=baron_pow_callback)


def options_section(options):
    """Write the OPTIONS section."""
    opts = '\n'.join(f'{opt}: {val};' for opt, val in options.items())
    return f'OPTIONS {{\n{opts}\n}}\n\n' if options else ''


def _name(prefix, n=1, i=None):
    """Generate n unused names of the form {prefix}{i} starting from i."""
    if i is None:
        try:
            i = _name.i[prefix]
        except (AttributeError, KeyError):
            i = 0
    for _ in range(n):
        _name.i[prefix] = i + 1
        yield f'{prefix}{i}'
        i += 1


_name.i = {}

var_name = partial(_name, 'x')
con_name = partial(_name, 'c')

# TESTS
# y = comando.Variable('y', Integers, (0, 1), 0, True)
# y.value = {1: 0, 2: 0}
# x1 = comando.Variable('x1', bounds=(0, None))
# x3 = comando.Variable('x3', Integers)
# x4 = comando.Variable('x4', bounds=(0, 100))
# x5 = comando.Variable('x5', bounds=(-300, None))
# x6 = comando.Variable('x6', bounds=(0, None))
# x7 = comando.Variable('x7', Integers, bounds=(10, None))
# k = comando.Variable('k', bounds=(2, 50), indexed=True, init_val=42)
# k.value = {1: 42, 2: 42}
# variables = {y, x1, x3, x4, x5, x6, x7, k}
# prios = {x3: 10, x5: 0}


def variables_section(var_map, prios=None):
    """Write the (BINARY/INTEGER/POSITIVE) VARIABLES sections."""
    all_vars = set(var_map)
    # continuous vs. integer
    c_vars, i_vars = split(all_vars, lambda v: v.is_integer)
    # general integer vs. binary
    i_vars, b_vars = split(i_vars, lambda v: v.is_binary)
    # general vs. 'positive' (lb = 0 for BARON)
    _, p_vars = split(all_vars - b_vars,
                      lambda v: all(v.lb == 0) if v.indexed else v.lb == 0)
    gc_vars = c_vars - p_vars  # 'general continuous' variables
    variable_groups = {'BINARY_VARIABLES': b_vars,
                       'INTEGER_VARIABLES': i_vars,
                       'POSITIVE_VARIABLES': p_vars,
                       'VARIABLES': gc_vars}
    res = ''
    for group, vars in variable_groups.items():  # sections for each group
        if vars:
            var_reps = (', '.join(var_map[v].values()) if v.indexed
                        else var_map[v] for v in vars)
            res += f"""{group} {', '.join(var_reps)};\n"""

    for var, repr in var_map.items():
        if var.indexed:
            for i, var_i_repr in repr.items():
                res += f'// {var_i_repr}: {var}[{i}]\n'
        else:
            res += f'// {repr}: {var}\n'

    lbs = i_vars.union(c_vars) - p_vars
    # lbnds_gen = ('\n'.join(f'{r}: {v[i].lb};' for i, r in var_map[v].items()
    #                        if not isinf(v[i].lb)) if v.indexed
    #              else ('' if isinf(v.lb) else f'{var_map[v]}: {v.lb}; ')
    #              for v in lbs)
    # lbnds = '\n'.join(lbnds_gen)

    def lbnds_gen():
        for v in lbs:
            if v.indexed:
                yield from (f'{r}: {v[i].lb};' for i, r in var_map[v].items()
                            if not isinf(v[i].lb))
            elif not isinf(v.lb):
                yield f'{var_map[v]}: {v.lb}; '

    lbnds = '\n'.join(lbnds_gen())

    ubs = lbs.union(p_vars)
    # ubnds_gen = ('\n'.join(f'{r}: {v[i].ub};' for i, r in var_map[v].items()
    #                        if not isinf(v[i].ub)) if v.indexed
    #              else ('' if isinf(v.ub) else f'{var_map[v]}: {v.ub}; ')
    #              for v in ubs)
    # ubnds = '\n'.join(ubnds_gen)

    def ubnds_gen():
        for v in ubs:
            if v.indexed:
                yield from (f'{r}: {v[i].ub};' for i, r in var_map[v].items()
                            if not isinf(v[i].ub))
            elif not isinf(v.ub):
                yield f'{var_map[v]}: {v.ub}; '

    ubnds = '\n'.join(ubnds_gen())

    for b_type, bnds in ('LOWER', lbnds), ('UPPER', ubnds):
        if bnds:
            res += f"\n{b_type}_BOUNDS {{\n{bnds}\n}}\n"

    if prios:
        prio_gen = ('\n'.join(f'{r}: {p};' for r in var_map[v].values())
                    if v.indexed
                    else f'{var_map[v]}: {p};' for v, p in prios.items())
        prio_section = '\n'.join(prio_gen)
        if prio_section:
            res += f"\nBRANCHING_PRIORITIES {{\n{prio_section}\n}}\n"
    return res


def constraints_section(con_map, rel_only_cons, convex_cons, parse):
    """Write the (RELAXATION_ONLY/CONVEX) EQUATIONS sections."""
    # Assuming constraints is a dict con->name and rel_only_cons and
    # convex_cons are sets of cons
    ro_cons = {con: con_map[con] for con in rel_only_cons}
    conv_cons = {con: con_map[con] for con in convex_cons}
    constraint_groups = {'EQUATIONS': con_map,
                         'RELAXATION_ONLY_EQUATIONS': ro_cons,
                         'CONVEX_EQUATIONS': conv_cons}

    res = ''
    for group, cons in constraint_groups.items():
        if cons:
            res += f"""{group} {', '.join(', '.join(name.values())
                                          if indexed(con) else name
                                          for con, name in cons.items())};\n"""

    con_defs = ('\n'.join(f'{n[i]}: {parse(c, i)};' for i in get_index(c))
                if indexed(c) else f'{n}: {parse(c)};'
                for c, n in con_map.items())

    res += '\n\n' + '\n'.join(con_defs)
    return res


def objective_section(P, parse):
    do = parse(P.design_objective)
    ooe = P.operational_objective
    if ooe == 0:
        return f"\nOBJ: minimize {do};"
    ts = P.timesteps
    if ts is None:
        oo = ' + '.join(f'{p} * ({parse(ooe, s)})'
                        for s, p in P.scenario_weights.items())
    elif P.scenarios is None:
        oo = ' + '.join(f'{dt} * ({parse(ooe, t)})' for t, dt in ts.items())
    else:
        oo = " + ".join(f"""{p} * ({' + '.join(f'{dt} * ({parse(ooe, (s, t))})'
                                              for t, dt in ts[s].items())})"""
                        for s, p in P.scenario_weights.items())
    return f"\nOBJ: minimize {do} + {oo};"


def start_section(var_map):
    s_gen = ('\n'.join(f'{r}: {val};'
                       for r, val in zip(var_map[v].values(), v.value))
             if v.indexed else f'{var_map[v]}: {v.value};' for v in var_map)
    _start_section = '\n'.join(s_gen)
    return f"\nSTARTING_POINT {{\n{_start_section}\n}}\n"
#
# par = comando.Parameter('p')
# par.value = {0: 40, 1: 50}
# constraints = {
#     5*x3 + y - 3*x5 ** 3 >= 1: ['e11', 'e12'],
#     Eq(y + 2*x4 - 2*x7, 25.7): ['e21', 'e22'],
#     -20 <= x4 + 2*y*x3 + x6: ['e3lb1', 'e3lb2'],
#     x4 + 2*y*x3 + x6 <= par: ['e3ub1', 'e3ub2']
# }
# rel_only_cons = {Eq(y + 2*x4 - 2*x7, 25.7)}
# convex_cons = {5*x3 + y - 3*x5 ** 3 >= 1}
#
# sym_map = {v: {i: f'{v.name}{i}' for i in index} if v.indexed else v.name
#            for v in set().union(*(get_vars(con) for con in constraints))}

def _spread_bounds(lb, ub, frac=0.1):
    if lb < 0:
        lb *= (1 + frac)
    else:
        lb *= (1 - frac)
    if ub < 0:
        ub *= (1 - frac)
    else:
        ub *= (1 + frac)
    return lb, ub


# DEBUG:
# For very long expressions BARON seems to have issues, terminating without
# having reached the specified termination criteria and complaining that:
#     User did not provide appropriate variable bounds.
#     Some model expressions are unbounded.
# or incorrectly reporting problems as infeasible when intermediate expressions
# are introduced as variables with appropriate bounds.
# One issue might be expressions whose bounds are identical, however, even
# spreading the bounds of such expressions results in the reported
# infeasibility.
# Using extreme bounds that are certain to be valid again results in the
# complaint about unbounded expressions...
# b = -1e10, 1e10
def apply_cse(P, var_map):
    """Replace reoccurring expressions with variables."""
    do = P.design_objective
    oo = P.operational_objective
    cons = P.constraints
    reps, exprs = comando.cse((do, oo, *cons.values()))
    defs = {}
    cons = {}
    n = len(P.index)
    for sym, rep in reps:
        e = rep.subs(defs)
        # print(defs)
        index = get_index(e)
        if index is None:
            b = comando.utility.bounds(e)
            # try:
            # except ValueError:
            #     b = -1e10, 1e10
            if b[1] - b[0] < 1e-15:
                b = _spread_bounds(b[0], b[1])
            #     # print(e, b)
            #     # for s in e.free_symbols:
            #     #     print(type(s), s, end=' ')
            #     #     if s.is_Parameter:
            #     #         print(s.value)
            #     #     else:
            #     #         print(s.bounds)
            #     # const_val = (b[0] + b[1]) * 0.5
            #     # defs[sym] = const_val
            #     # print(sym, const_val)
            #     # continue
            x = comando.Variable(f'intermediate{sym.name[1:]}', bounds=b)
            var_map[x] = next(var_name())
        else:
            b = comando.utility.bounds(e)
            b = (min(b[0]), max(b[1]))
            # try:
            # except ValueError:
            #     b = -1e10, 1e10
            if b[1] - b[0] < 1e-15:
                b = _spread_bounds(b[0], b[1])
            #     # print(e, b)
            #     # for s in e.free_symbols:
            #     #     print(type(s), s, end=' ')
            #     #     if s.is_Parameter:
            #     #         print(s.value)
            #     #     else:
            #     #         print(s.bounds)
            #     # const_val = (b[0] + b[1]) * 0.5
            #     # defs[sym] = const_val
            #     # print(sym, const_val)
            #     # continue
            x = comando.VariableVector(f'intermediate{sym.name[1:]}', bounds=b)
            x.instantiate(get_index(e))
            var_map[x] = dict(zip(index, var_name(len(index))))
        # print(sym, b)
        cons[f'{x.name}_def'] = Eq(x, e)
        defs[sym] = x

    P2 = comando.Problem(exprs[0].subs(defs), exprs[1].subs(defs),
                         {**{con_id: expr.subs(defs) for con_id, expr in
                             zip(P.constraints, exprs[2:])}, **cons}, None,
                         P.timesteps, P.scenarios,
                         name=f'CSE_reformulation_of_{P.name}')
    return P2


def handle_tanh(expr, tanh_definitions, var_map):
    """Search the expression for occurrences of tanh and substitute them."""
    tanh_occurrences = {}
    for atom in expr.atoms(comando.Function):
        if (comando.utility.get_type_name(atom)) == 'tanh':
            arg = atom.args[0]
            index = get_index(arg)
            if index is None:
                tanh_var = comando.Variable(f'tanh{handle_tanh.n}',
                                            bounds=(-1, 1))
                var_map[tanh_var] = next(var_name())
            else:
                tanh_var = comando.VariableVector(f'tanh{handle_tanh.n}',
                                                  bounds=(-1, 1))
                tanh_var.instantiate(index)
                var_map[tanh_var] = dict(zip(index, var_name(len(index))))
            handle_tanh.n += 1
            tanh_occurrences[atom] = tanh_var
            tanh_definitions[tanh_var] = 1 - 2/(comando.exp(2 * arg) + 1)
    return expr.xreplace(tanh_occurrences)


handle_tanh.n = 0


def normalize(con):
    """Bring constraints to a normal form baron can handle."""
    try:
        return Le(con.lts - con.gts, 0)
    except AttributeError:
        return Eq(0, con.lhs - con.rhs)


# TODO: This may be general enough to move it to utility
def discretize(P, sym_map):
    """Discretize differential equations in P using implicit Euler."""
    skipped_indices = set()
    timesteps = P.timesteps
    Delta_t = P.Delta_t

    from comando.utility import get_previous

    sc = {}  # state constraints
    for state, state_data in P.states.items():
        initial_state, derivative, state_expression = state_data

        prev = get_previous(state, initial_state)
        prev_map = {}
        for i, v in prev.expansion.items():
            try:
                prev_map[i] = sym_map[v]
            except KeyError:
                prev_map[i] = sym_map[initial_state]
        sym_map[prev] = prev_map
        con = comando.Eq((state - prev) / Delta_t, derivative)
        for i in get_index(con):
            if i in skipped_indices:
                continue
            if timesteps[i] == 0:
                print(f'WARNING: index {i} corresponds to a length of 0... '
                      'skipping!')
                skipped_indices.add(i)
                continue
            sc[f'{derivative.name}_{i}'] = con
    return sc


def write_bar_file(P, file_name, options=None, cse=False, reuse=False):
    """Write a baron input file for problem P."""
    if options is None:
        options = {}
    prios = options.pop('branching_priorities', None)
    n = len(P.index)
    _name.i = {}
    # var_map = {v: dict(zip(P.index, var_name(n))) if v.indexed
    #            else next(var_name())
    #            for v in P.design_variables.union(P.operational_variables)}
    var_map = {}
    for v in P.design_variables.union(P.operational_variables):
        if v.indexed:
            index = v.expansion.index
            var_dict = dict(zip(index, var_name(len(index))))
            for v_i, name_i in zip(v, var_dict.values()):
                var_map[v_i] = name_i
            var_map[v] = var_dict
        else:
            var_map[v] = next(var_name())

    if cse:
        P = apply_cse(P, var_map)

        # tanh handling
        # TODO: Do not overwrite but use copies of constraints and objective!
        tanh_definitions = {}
        P.design_objective = handle_tanh(P.design_objective,
                                         tanh_definitions, var_map)
        P.operational_objective = handle_tanh(P.operational_objective,
                                              tanh_definitions, var_map)
        for c_id, con in P.constraints.items():
            P.constraints[c_id] = handle_tanh(con, tanh_definitions, var_map)
        for tanh_var, tanh_def in tanh_definitions.items():
            P.constraints[f'{tanh_var.name}_def'] = Eq(tanh_var, tanh_def)

    if P.states:
        print('WARNING: There are differential constraints that will be '
              'discretized using implicit Euler discretization, if you want '
              'an advanced time discretization you must reformulate your '
              'problem prior to passing it to the BARON interface!')
        # NOTE: state discretization will create constraints that directly
        #       contain state VeriableVector elements so we add them here!
        for state in P.states:
            var_map.update(zip(state.expansion, var_map[state].values()))
        cons = {**P.constraints, **discretize(P, var_map)}
    else:
        cons = P.constraints

    parser = BaronParser(var_map)
    parse = parser.cached_parse

    def check(c_id, con):
        """Check whether a constraint has variables, if not evaluate."""
        if get_vars(con):
            return True
        # If constraint doesn't have any variables...
        res = comando.utility.parse(con,  # ...test whether it is satisfied
                                    {p: p.value for p in con.free_symbols})
        if not all(res) if indexed(con) else not res:  # if there's a violation
            from comando import ImpossibleConstraintException
            raise ImpossibleConstraintException(
                f'Constraint "{c_id}" contains only parameters and '
                'evaluates to False!')
        else:  # else continue (returns None which is equivalent to False!)
            print(f'INFO: Constraint "{c_id}" contains only parameters and '
                  'evaluates to True! Skipping...')

    con_map = {normalize(c): dict(zip(P.index, con_name(n))) if indexed(c)
               else next(con_name())
               for c_id, c in cons.items() if check(c_id, c)}

    if not reuse:
        with open(file_name, 'w') as f:
            f.write(options_section(options))
            f.write(variables_section(var_map, prios))
            f.write(constraints_section(con_map, {}, {}, parse))
            f.write(objective_section(P, parse))
            f.write(start_section(var_map))
    return var_map, con_map

# TESTS f = open(input_file_name, 'w')
# str_parse(-b, sym_map)
# str_parse(a + -b, sym_map)
# str_parse(1 + a + -b, sym_map)
# str_parse(a - (b + 1), sym_map)
# str_parse(a + 1 - (b + 1), sym_map)
# str_parse(2 ** a + -b, sym_map)
# str_parse(a + -(b + 1), sym_map)
# str_parse(a + -b - c, sym_map, idx=1)
# str_parse(a + -log(b), sym_map)
# str_parse(a + -log(b + 1), sym_map)
# str_parse(a + -exp(b + 1), sym_map)
# str_parse(a ** (log(b) + 1), sym_map)
# str_parse(a ** -(log(b) + 1), sym_map)
# str_parse(-(b + 1), sym_map)
# str_parse(a ** -1.2, sym_map)
#
# str_parse(a - b, sym_map)
# str_parse(a + -(1 + b), sym_map)
# str_parse(a + c**-(1 + b), sym_map, idx=1)
# str_parse(a + c**log(1 + b), sym_map, idx=1)
# str_parse(log(b + 1), sym_map, idx=1)


def get_results(results_file_name='res.lst'):
    """Code for parsing baron results files."""
    import re
    p = re.compile(r"  (\S+)\s+\S+\s+(\S+)")

    val_map = {}
    with open(results_file_name, 'r') as f:
        # Advance until the results section
        # while not f.readline().startswith('  variable'):
        #     continue
        ready = False
        for line in f.readlines():
            if not ready:
                ready = line.startswith('  variable')
                continue
            try:
                var_name, val = p.search(line).groups()
            except AttributeError:
                pass

            val_map[var_name] = float(val)
    return val_map


def solve(P, file_name=None, silent=False, cse=False,
          reuse=None, **options):
    """Solve the problem specified in the input_file with baron."""
    from comando.utility import syscall, canonical_file_name, \
        check_reuse_or_overwrite

    base_name, file_name = canonical_file_name(P.name, '.bar', file_name)
    check_reuse_or_overwrite(file_name, reuse)

    for out_name in 'ResName', 'SumName', 'TimName':
        if out_name not in options:
            options[out_name] = f'"{base_name}.{out_name[:3].lower()}.lst"'
    for bool_option in ['results', 'summary', 'times']:
        if bool_option in options:
            options[bool_option] = int(options[bool_option])  # bool to int!

    var_map, con_map = write_bar_file(P, file_name, options, cse, reuse)
    result_file_name = options['ResName'][1:-1]

    log_name = f'{base_name}.baron.log'
    ret = syscall('baron', file_name, log_name=log_name, silent=silent)
    if ret and not silent:
        print(f'BARON terminated with return code {ret}!')
        return ret
    val_map = get_results(result_file_name)
    if not val_map:
        return -1  # no results (e.g. problem infeasible, no license)
    for comando_var, baron_var in var_map.items():
        if comando_var.indexed:
            comando_var.value = {i: val_map[vi] for i, vi in baron_var.items()}
        else:
            comando_var.value = val_map[baron_var]

    return ret


def get_times_and_bounds(baron_log_file):
    """Code for parsing baron logs for time, and bounds."""
    import re
    p = re.compile(r"[ \*] +\S+ +\S+ +(\S+) +(\S+) +(\S+)")

    times = []
    lb_data = []
    ub_data = []
    lines = iter(baron_log_file.split('\n'))
    baron_log_file.split('\n')[0].startswith("  Iteration")
    for line in lines:
        if not line.startswith("  Iteration"):
            continue
        break
    else:
        raise RuntimeError
    lines = [*lines]
    for line in lines:
        try:
            t, lb, ub = p.search(line).groups()
        except AttributeError:
            break
        times.append(float(t))
        lb_data.append(float(lb))
        ub_data.append(float(ub))
    return times, lb_data, ub_data
