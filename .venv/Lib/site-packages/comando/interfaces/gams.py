"""Routines for translation to GAMS syntax."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import comando.utility
from comando.utility import indexed, split, _str_map, StrParser

gams_str_map = {'LessThan': lambda lhs, rhs: f'{lhs} =L= {rhs}',
                'GreaterThan': lambda lhs, rhs: f'{lhs} =G= {rhs}',
                'Equality': lambda lhs, rhs: f'{lhs} =E= {rhs}',
                'power': lambda base, exponent: f'power({base}, {exponent})',
                'exp': lambda arg: f'EXP({arg})',
                'log': lambda arg: f'LOG({arg})'}
for key in '()', 'Add', 'Neg', 'Sub', 'Mul', 'Div', 'Inv', 'Pow':
    gams_str_map[key] = _str_map[key]


def gams_pow_callback(parser, expr, idx):
    """Handle special pow calls in GAMS."""
    base, exponent = expr.args
    if base == comando.E:
        return parser.str_map['exp'](*parser.parse_args((exponent, ), idx))
    # NOTE: GAMS provides a special function for integer powers.
    #       the notation x ** y, is handled as exp(y * log(x))
    if exponent.is_Integer and exponent.is_positive:
        return parser.str_map['power'](*parser.parse_args((base, exponent, ),
                                                          idx))
    return None  # No special case, handle normally


class GamsParser(StrParser):  # pylint: disable=too-few-public-methods
    """A class for parsing comando expressions to GAMS Syntax."""

    def __init__(self, sym_map):
        super().__init__(sym_map, gams_str_map,
                         pow_callback=gams_pow_callback)


def _entry(name, value, i_name=None, field='', doc=''):
    """Create a string representation of a gams section entry."""
    entry = f'\t{name}'
    if i_name:
        entry += f"({i_name}) {doc}"
        if value is not None:
            delim = pad = '\n\t'
            data = delim.join(f'{idx}{field}\t{v}'
                              for idx, v in zip(_entry.defs[i_name], value))
            entry += f'{pad}/{pad}{data}{pad}/'
    else:
        entry += f' {doc}'
        if value is not None:
            delim = ', '
            pad = ' '
            try:  # iterable of values
                data = delim.join(str(v) for v in value)
            except TypeError:  # scalar
                data = str(value)
            entry += f'{pad}/{pad}{data}{pad}/'
        if not hasattr(_entry, 'defs'):
            _entry.defs = {}
        _entry.defs[name] = value
    return entry + '\n'


def literal(i):
    """Get a gams representation for the literal index i in parentheses."""
    if isinstance(i, tuple):  # scenario-timestep pairs
        idx = f"'{'_'.join(str(e) for e in i)}'"
    else:  # only timesteps
        idx = f"'{i}'"
    return f'({idx})'


GAMS_HEADER = '* GAMS cannot handle too many digits in numbers.\n' \
    '* This tells it to ignore those extra digits.\n' '$OFFDIGIT\n\n'


def write_sets_section(P):
    """Generate a string representation of a SETS section."""
    ts = P.timesteps
    if ts is None:
        sets = f"SETS\n\n{_entry('s', P.scenarios, doc='scenarios')}"
    elif P.scenarios is None:
        sets = f"SETS\n\n{_entry('t', ts.keys(), doc='timesteps')}"
    else:
        sets = f"SETS\n\n{_entry('s', P.scenarios, doc='scenarios')}"
        timesteps = ['_'.join(str(i) for i in st) for st in P.index]
        sets += f"\n{_entry('t', timesteps, doc='timesteps')}"
    return sets + ';'


def write_parameters_section(P, iname):
    """Generate a string representation of a PARAMETERS section."""
    entries = []
    if P.timesteps is not None:
        entries.append(_entry('Delta_t', P.timesteps, 't'))
    if P.scenarios is not None:
        entries.append(_entry('prob', P.scenario_weights, 's'))
    param_entries = (*entries,
                     *(_entry(p.name, p.value, iname) if p.indexed
                       else _entry(p.name, p.value) for p in P.parameters))
    return f"PARAMETERS\n\t\n{''.join(param_entries)};"


def write_variables_section(P, iname):
    """Write the different 'VARIABLES' sections and variable bounds."""
    var_types = {'BINARY ': [], 'INTEGER ': [], '': []}  # last == continuous
    all_vars = P.design_variables.union(P.operational_variables)
    for v in all_vars:  # Sort the variabes according to their domain
        if v.is_binary:
            var_types['BINARY '].append(v)
        elif v.is_integer:
            var_types['INTEGER '].append(v)
        else:
            var_types[''].append(v)

    def _var_values(var):  # A helper method for variable entries
        """Convert a possibly indexed variable into GAMS representation."""
        if var.indexed:
            return (v.value, iname)
        else:
            return (None, None) if var.value is None \
                else ([f'L {var.value}'], None)
    variables = ""
    for prefix, vars in var_types.items():
        if vars:  # Only print this section if there are variables of this type
            # var_entries = (_entry(v.name, *_var_values(v)) for v in vars)
            # NOTE: For some odd reason neither of the following works:
            #       [_entry(v.name, *_var_values(v), '.L') for v in vars]
            #
            #       def var_entries():
            #           for v in vars:
            #               yield _entry(v.name, *_var_values(v), '.L')
            var_entries = []
            for v in vars:
                var_entries.append(_entry(v.name, *_var_values(v), '.L'))
            variables += f"{prefix}VARIABLES\n\n{''.join(var_entries)};\n\n"

    # Variable bounds
    def _bounds(v):  # A helper method for variable bound entries
        """Return the string to set the bounds of variable v."""
        b_types = ['LO', 'UP']
        if v.indexed:
            # print the default bounds and the indexed ones that are different:
            # e.g. for a positive vector foo -> foo(t, s).LO = 0;
            res = ' '.join(f'{v.name}.{b_type}({iname}) = {val};'
                           for b_type, val in zip(b_types, v._bounds)
                           if val is not None) + '\n'
            # if in scenario 's1' foo is bounded by 1 from above this gives:
            # foo('s1_t1').UB = 1;
            # foo('s1_t2').UB = 1
            # ...
            for i, iv in v.expansion.items():
                data = zip(b_types, iv.bounds, v._bounds)
                iv_entry = ' '.join(f'{v.name}.{b_type}{literal(i)} = {val};'
                                    for b_type, val, default in data
                                    if val not in [default, None])
                if iv_entry:
                    res += f'{iv_entry}\n'
            return res
        return ' '.join(f'{v.name}.{b_type} = {val};' for b_type, val in
                        zip(b_types, v.bounds) if val is not None) + '\n'

    bounds = ''.join(_bounds(v) for v in var_types['INTEGER '] + var_types[''])
    if bounds:
        variables += f'* BOUNDS\n{bounds}'
    return variables


def write_equations_section(P, iname, parse):
    """Generate a string representation of a EQUATIONS section."""
    def eq(id):
        """Normalize the equation id if it doesn't fit GAMS syntax."""
        import re
        if re.match(r'^\w+$', id):
            eq_id = id
        else:
            eq_id = f'constraint_{eq.n}'
            eq.n += 1
        eq.map[eq_id] = id
        return eq_id
    eq.n = 0
    eq.map = {}

    cons = {}
    for c_id, c in P.constraints.items():
        if comando.utility.indexed(c):
            cons[f'{eq(c_id)}({iname})'] = c
        else:
            cons[eq(c_id)] = c

    declarations = "\n".join(cons)
    if declarations:
        equations = f'EQUATIONS\n\n{declarations};\n\n* INSTANTIATIONS\n' + \
            "\n".join(f'{c_id}.. {parse(c)};'
                      for c_id, c in cons.items())
    return equations


def write_objective(P, parse):
    """Create the string for the objective function of P."""
    return parse(P.objective)  # simple solution
    # explicit flat parsing
    # ts = [*P.timesteps.items()]
    # if P.scenarios is None:
    #     oo = ' + '.join(f'{dt} * ({parse(ooe)})' for t, dt in ts)
    # else:
    #     oo = " + ".join(f"""{p} * ({' + '.join(f'{dt} * ({parse(ooe, (s, t))})'
    #                                            for t, dt in ts)})"""
    #                     for s, p in P.scenarios.items())
    # return f"{do} + {oo}"

# DEBUG:
# S.existing_components.clear()
# S = comando.System('TEST')
# p = S.make_parameter('p')
# q = S.make_parameter('q')
# x = S.make_design_variable('x')
# y = S.make_operational_variable('y')
# S.add_eq_constraint(x ** 2, p, 'dc')
# S.add_eq_constraint(y ** 0.5, -x * q, 'oc')
# P = S.create_problem(x, y, timesteps=[1,2,3])
# P.set_values({'TEST_p': 0.5, 'TEST_q': {1: -42, 2: -15, 3: 0}})
# P.get_data()
# DEBUG
def populate_sym_map(iname, dvars=(), ovars=(), pars=(), sym_map=None):
    """Create a GAMS symbol map or populate an existing one with new symbols.

    Arguments
    ---------
    dvars : iterable
        design variables
    ovars : iterable
        operational variables
    pars : iterable
        parameters

    Returns
    -------
    sym_map : dict
        Symbol map with GAMS representation for each passed symbol.
    """
    from itertools import chain
    import re
    p = re.compile(r"(\w+)\[\(?'?(\w+)'?(?:, '?(\w+)'?)*\)?\]")

    if sym_map is None:
        sym_map = {}
    dpars, opars = split(pars, indexed)
    for sym in chain(dpars, dvars):
        sym_map[sym] = sym.name
    for sym in chain(opars, ovars):
        sym_map[sym] = f'{sym.name}({iname})'
        for i, symi in sym.expansion.items():
            name, *indices = p.match(symi.name).groups()
            sym_map[symi] = f"{name}('{'_'.join(i for i in indices if i)}')"  # f'{sym.name}{literal(i)}'
    return sym_map


def write_gms_file(P, file_name, model_type='MINLP'):  # , flat=False):
    """Write a GAMS file based on the COMANDO Problem."""
    name = P.name
    iname = 't' if P.index.name is None else P.index.name

    # Parsing setup
    sym_map = populate_sym_map(iname, P.design_variables,
                               P.operational_variables,
                               P.parameters)
    if P.states:
        from comando.utility import implicitEuler
        from comando import is_trivial
        state_constraints, prev_states = implicitEuler(P)
        # P.constraints.update(state_constraints)
        # populate_sym_map(ovars=prev_states, pars=[P.Delta_t], sym_map=sym_map)
        populate_sym_map(pars=[P.Delta_t], sym_map=sym_map)
        for c_id, con in state_constraints.items():
            for i, idx in enumerate(P.index):
                # BUG: Gives False for some constraints!
                # con_i = comando.utility.parse(con, idx=idx)
                con_i_lhs = comando.utility.parse(con.lhs, idx=idx)
                con_i_rhs = comando.utility.parse(con.rhs, idx=idx)
                con_i = comando.Eq(con_i_lhs, con_i_rhs)
                c_name = f'{c_id}_{i}'
                if is_trivial(con_i):
                    print(f'WARNING: Constraint {c_name} is trivially '
                          'satisfied, skipping...')
                else:
                    P.constraints[c_name] = con_i

    parse = GamsParser(sym_map)
    # OBJECTIVE
    # model_type = 'MINLP'  # TODO: Automatic selection
    # + (f"OPTION {model_type} = {solver};\n" if solver else "")
    objective = "VARIABLE OBJ;\nEQUATION OBJ_EXPRESSION;\n" \
        f"OBJ_EXPRESSION.. OBJ =E= {write_objective(P, parse)};\n\n" \
        f"MODEL {name} /all/;\n" \
        f"SOLVE {name} USING {model_type} minimising OBJ;\n"
    with open(file_name, 'w') as f:
        f.write(GAMS_HEADER)
        f.write('\n\n')
        f.write(write_sets_section(P))
        f.write('\n\n')
        f.write(write_parameters_section(P, iname))
        f.write('\n\n')
        f.write(write_variables_section(P, iname))
        f.write('\n\n')
        f.write(write_equations_section(P, iname, parse))
        f.write('\n\n')
        f.write(objective)
    # return sym_map


def get_results(results_file_name):
    """Code for parsing gams results files."""
    import re
    p1 = re.compile(r"---- VAR (\S+)\s+\S+\s+(\S+)")
    p2 = re.compile(r"---- VAR (\S+)\s+\n")
    p3 = re.compile(r"(\S+)\s+\S+\s+(\S+)")

    val_map = {}

    with open(results_file_name, 'r') as f:
        lines = iter(f.readlines())

        def skip_header(line):
            if line.startswith('\x0cGAMS '):
                next(lines)
                next(lines)
                return next(lines)
            return line

        def advance(n=1):
            for i in range(n):
                line = skip_header(next(lines))
            return line

        for line in lines:
            try:  # Match scalar variable
                var_name, val = p1.search(line).groups()
                val_map[var_name] = 0 if val == '.' else float(val)
                advance()
                continue
            except AttributeError:
                pass
            try:  # Match indexed variable
                var_name = p2.search(line).group(1)
                line = advance(3)
                values = {}
                while True:
                    try:
                        line = skip_header(next(lines))
                        idx, value = p3.search(line).groups()
                        values[tuple(idx.split('.'))] = 0 if val == '.' \
                            else float(value)
                    except AttributeError:
                        val_map[var_name] = values
                        advance()
                        break

            except AttributeError:
                pass
            if line.startswith('**** REPORT SUMMARY'):
                break
    return val_map


def solve(P, input_file_name=None, silent=False, **options):
    """Solve the problem specified in the input_file using GAMS.

    Apart from the usual GAMS options, the model type (e.g. LP/MINLP) may be
    specified explixitly using the model_type option, the default is MINLP.
    """
    from comando.utility import syscall
    import os

    if input_file_name is None:
        base_name = P.name
        input_file_name = f'{P.name}.gms'
    elif input_file_name.endswith('.gms'):
        base_name = input_file_name[:-4]
    else:
        base_name = input_file_name
        input_file_name = base_name + '.gms'

    while os.path.isfile(input_file_name):
        yn = input(f"A file '{input_file_name}' already exists, overwrite "
                   "(y/n)?")
        if yn.lower() == 'y':
            break
        if yn.lower() == 'n':
            print('Aborting...')
            return 1
    model_type = options.pop('model_type', 'MINLP')
    write_gms_file(P, input_file_name, model_type)

    ret = syscall('gams', input_file_name,
                  *(f'{option}={value}' for option, value in options.items()),
                  log_name=f'{base_name}.gams.log', silent=silent)
    if ret != 0:
        print('Something went wrong during solving!')
    else:
        val_map = get_results(f'{base_name}.lst')
        for v in P.design_variables.union(P.operational_variables):
            v.value = val_map[v.name]

    return ret


######
# def is_num(x):
#     """Test whether x is numeric."""
#     return isinstance(x, (int, float))


# def allnums(*args):
#     """Test if all arguments passed are numerical."""
#     return all(is_num(arg) for arg in args)

#
# def is_neg(x):
#     """Test if x is a negated summand and store the result and negation."""
#     if x in is_neg.cache:
#         return is_neg.cache[x][0]
#     try:
#         if x[0] == '-':
#             is_neg.cache[x] = [True, x[1:]]
#             return True
#     except TypeError:
#         if x < 0:
#             is_neg.cache[x] = [True, str(-x)]
#             return True
#     return False
#
#
# is_neg.cache = {}
#
#
# def gams_sum(*args):
#     """Create a string representation of `sum(args)` in GAMS syntax."""
#     summands, negs = comando.utility.split(args, is_neg)
#     res = ' + '.join(str(s) for s in summands)
#     if negs:  # negs is nonempty!
#         if len(negs) > 1:
#             subtrahends = f"[{' + '.join(is_neg.cache[n][1] for n in negs)}]"
#         else:
#             subtrahends = is_neg.cache[negs[0]][1]
#         return f'{res} - {subtrahends}'
#     return res
#
#
# def is_inv(x):
#     """Test if x is an inverted factor and store the result and inversion."""
#     if x in is_inv.cache:
#         return is_inv.cache[x][0]
#     try:
#         if x[:2] == '1/':
#             is_inv.cache[x] = [True, x[2:]]
#             return True
#     except TypeError:
#         pass
#     return False
#
#
# is_inv.cache = {}
#
#
# def gams_prod(*args):
#     """Create a string representation of `product(args)` in GAMS syntax."""
#     negate = False
#     if args[0] == -1:
#         negate = True
#         args = args[1:]
#     factors, invs = comando.utility.split(args, is_inv)
#     res = ' * '.join(f'({f})' if '+' in f else str(f) for f in factors)
#     if invs:  # invs is nonempty!
#         if len(invs) > 1:
#             dividends = f"[{' + '.join(is_inv.cache[i][1] for i in invs)}]"
#         else:
#             dividends = is_inv.cache[invs[0]][1]
#         res = f'{res} / {dividends}'
#     if negate:
#         return f"-[{res}]" if ' ' in res else f'-{res}'
#     return res
#
#
# def gams_pow(*args):
#     """Create a string representation of `pow(*args)` in GAMS syntax."""
#     if args[1] == -1:
#         return f'1/[{args[0]}]'
#     if args[1] == 2:
#         return f'SQR({args[0]})'
#     try:  # test if the exponent is integer
#         if int(args[1]) == args[1]:
#             return f"POWER({args[0]}, {args[1]})"
#     except ValueError:
#         pass
#     # real power (equivalent to EXP[n*LOG(x)])
#     return f"RPOWER({args[0]}, {args[1]})"
#
#
# gams_map = {
#     'Add': gams_sum,
#     'Mul': gams_prod,
#     'Pow': gams_pow,
#     'LessThan': lambda x, y: f"{x} =L= {y}",
#     'GreaterThan': lambda x, y: f"{x} =G= {y}",
#     'Equality': lambda x, y: f"{x} =E= {y}"}
#
#
# def experimental():
#     from gams import GamsWorkspace
#     import os
#     import sys
#     cwd = os.getcwd()
#     filename = 'test_if_else.gms'
#     ws = GamsWorkspace()  # system_directory='/Applications/GAMS25.0/sysdir/')
#     filename = "Voll.gms"
#     t1 = ws.add_job_from_file(os.path.join(cwd, filename))
#     t1.run(output=sys.stdout)
#
#     for o in t1.out_db:
#         print(o.name)
#         for record in o:
#             print('\t', record)
#         # try:
#         #
#         #         # if record.keys:
#         #         #     keys = ", ".join(str(key) for key in record.keys)
#         #         #     print(f'\n\t({keys}): level = {record.level}, '
#         #         #           f'marginal = {record.marginal}')
#         #         # else:
#         #         #     print(f'\n\tlevel = {record.level}, '
#         #         #           f'marginal = {record.marginal}')
#         # except:
#         #     for record in o:
#         #         print('\t', record)
