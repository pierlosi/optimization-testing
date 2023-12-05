"""Utility functions used in various parts of COMANDO."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
from collections import defaultdict
import contextlib
import math
import os
import operator
from functools import reduce, partial

import sympy
import numpy
from numpy import minimum, maximum

import comando


def canonical_file_name(base_name, suffix, file_name=None) -> (str, str):
    """Get a canonical file name and the corresponding base name.

    Arguments
    ---------
    base_name : str
        the fallback if filename is None
    suffix : str
        desired suffix
    file_name : str or None
        a desired name with or without the desired suffix or None

    Returns
    -------
    base_name : str
        the canonical base_name
    file_name : str
        the canonical file_name
    """
    if file_name is None:
        file_name = f'{base_name}{suffix}'
    elif file_name.endswith(suffix):
        base_name = file_name[:-4]
    else:
        base_name = file_name
        file_name = base_name + suffix
    return base_name, file_name


def check_reuse_or_overwrite(file_name, reuse) -> None:
    """Check whether a file exists for reuse or overwriting.

    Arguments
    ---------
    file_name : str
        the name of the file to reuse or overwrite
    reuse : bool or None
        whether to reuse the given file. There are three cases:

        1. reuse is False -> do nothing
        2. reuse is True -> check if file_name is a valid file
        3. reuse is None -> if file_name is a file, ask before overwriting
    """
    if reuse is False:
        return
    if reuse is True:
        if not os.path.isfile(file_name):
            raise FileNotFoundError(file_name)
    else:  # reuse is None, explicitly ask
        while os.path.isfile(file_name):
            yn = input(f"A file '{file_name}' already exists, "
                       "overwrite (y/n)?")
            if yn.lower() == 'y':
                break
            if yn.lower() == 'n':
                raise FileExistsError(f'File {file_name} exists! '
                                      'Aborting...')


def syscall(executable, *args, log_name=os.devnull, silent=False):
    """Issue a system call.

    Arguments
    ---------
    executable : str
        the name of the executable to be called. Needs to be an executable in
        one of the directories listed in the PATH environment variable.
    args : list of str
        the arguments passed to the executable
    log_name : str
        the location of a file to which the output of the system call is
        mirrored in parralel to execution (default: os.devnull)
    silent : bool
        whether the output should be hidden from stdout (default: False)

    Returns
    -------
    ret : int
        the return code of the system call
    """
    from subprocess import Popen, PIPE, STDOUT
    import sys

    if silent:
        with open(log_name, 'w') as f:
            process = Popen([executable, *args], stderr=STDOUT, stdout=f)
    else:
        process = Popen([executable, *args], stderr=STDOUT, stdout=PIPE)
        with open(log_name, 'w') as f:
            while True:
                output = process.stdout.readline()
                if process.poll() is not None and not output:
                    break
                line = output.decode('UTF-8')
                sys.stdout.write(line)
                f.write(line)
    process.wait()  # wait until executable has finished
    ret = process.returncode
    return ret


def implicitEuler(P):
    """Discretize differential equations in P using implicit Euler."""
    from comando.utility import get_previous

    P.states
    Delta_t = P.Delta_t

    sc = {}  # state constraints
    previous_states = []
    skipped_indices = set()
    for state, state_data in P.states.items():
        initial_state, derivative, state_expression = state_data

        prev = get_previous(state, initial_state)
        previous_states.append(prev)
        # prev_map = {}
        # for i, v in prev.expansion.items():
        #     try:
        #         prev_map[i] = sym_map[v]
        #     except KeyError:
        #         prev_map[i] = sym_map[initial_state]
        # new_syms[prev] = prev_map
        con = comando.Eq((state - prev) / Delta_t, derivative)
        for i in get_index(con):
            if i in skipped_indices:
                continue
            if Delta_t.value[i] == 0:
                print(f'WARNING: Index {i} corresponds to a time-step of zero '
                      'length, so it is skipped!')
                skipped_indices.add(i)
                continue
            sc[f'{derivative.name}_{i}'] = con
        return sc, previous_states


@contextlib.contextmanager
def silence():
    """Execute code without output to stdout."""
    import sys
    with open(os.devnull, 'w') as f:
        save_stdout = sys.stdout
        sys.stdout = f
        yield
        sys.stdout = save_stdout


@contextlib.contextmanager
def cleanup():
    """Take a snapshot of all files on current path and remove new ones."""
    import os
    dir = os.getcwd()
    dir_contents = {*os.listdir('.')}
    yield
    os.chdir(dir)  # ensure we didn't change the directory!
    for f in {*os.listdir('.')} - dir_contents:
        os.remove(os.path.join('.', f))


def get_latest(path_glob):
    """Get the path to the latest file matching the path_glob."""
    import glob
    import os
    matches = glob.glob(path_glob)
    if matches:
        return __builtins__['max'](matches, key=os.path.getctime)
    else:
        raise FileNotFoundError(f'No file in {os.getcwd()} matches the glob '
                                f"'{path_glob}'!")


def depth(expr):
    """Compute the depth of the given expression."""
    if isinstance(expr, (int, float, comando.Number, comando.Symbol)):
        return 0
    else:
        return 1 + __builtins__['max']([depth(arg) for arg in expr.args])


def fuzzy_not(v):
    """Return None if `v` is None else `not v`."""
    if v is None:
        return v
    else:
        return not v


def smooth_abs(x, delta=1e-4):
    """Smooth approximation of abs(x)."""
    with sympy.evaluate(False):
        res = ((x + delta) ** 2) ** 0.5
    return res


def smooth_max(a, b, delta=1e-4):
    """Smooth approximation of max(a, b)."""
    return 0.5 * (a + b + smooth_abs(a-b, delta))


def smooth_min(a, b, delta=1e-4):
    """Smooth approximation of min(a, b)."""
    return 0.5 * (a + b - smooth_abs(a-b, delta))


# def min(iterable):
#     """Generalize min function to array-like structures."""
#     return reduce(minimum, iterable)
#
#
# def max(iterable):
#     """Generalize max function to array-like structures."""
#     return reduce(maximum, iterable)


# TODO: 'algebraic expression' has a different meaning in general (e.g. it
#       excludes transcendental expressions or numbers). Also see the arcticle
#       https://en.wikipedia.org/wiki/Algebraic_expression.
#       A better name might be _is_numeric
def _assert_algebraic(expr):
    """Assert that `expr` is an 'algebraic expression'.

    An algebraic expression is considered to be an inctance of type
    `comando.Expr` that is not of type `comando.Boolean`.
    An exception to this rule are symbols, which are always considered to be
    algebraic expressions.
    """
    # NOTE: Due to the imports in sympy's __init__.py, the type
    #       hierarchy does not correspond to the attribute hierarchy!
    assert isinstance(expr, (int, float)) or expr.is_Symbol or \
        isinstance(expr, comando.Expr) and not \
        isinstance(expr, comando.Boolean), \
        "expression should be 'algebraic sympy expressions'" \
        " (i.e., instances of type `comando.Expr` that are" \
        " not of type `comando.Boolean`)!"


# TODO: Rename to is_indexed for consistency
def indexed(expr):
    """Test whether the given expression is indexed."""
    return any(sym.indexed for sym in expr.free_symbols)


def is_negated(expr):
    """Check if expr is a negated expression or a negative number."""
    return (expr.is_Mul and (expr.args[0].is_Number
                             and bool(expr.args[0] < 0))) or \
        (expr.is_Number and bool(expr < 0)) or \
        expr.is_Add and all(is_negated(arg) for arg in expr.args)


def get_vars(expr):
    """Get a set of all variables in `expr`."""
    return {sym for sym in expr.free_symbols
            if isinstance(sym, (comando.Variable, comando.VariableVector))}


def get_pars(expr):
    """Get a set of all parameters in `expr`."""
    return {sym for sym in expr.free_symbols
            if isinstance(sym, (comando.Parameter))}


def get_previous(state, initial_state):
    """Create a VariableVector for the state at previous timesteps."""
    if initial_state == comando.cyclic:
        try:
            prev_data = state.expansion.groupby('s').shift(1)
            prev = comando.VariableVector(state.name + '__prev')
            prev.expansion = prev_data
            for s in prev.expansion.index.levels[0]:
                prev.expansion[s].iloc[0] = state.expansion[s].iloc[-1]
            return prev
        except KeyError:  # no scenarios
            init = state.expansion.tail(1)
    elif isinstance(initial_state, comando.Variable):
        init = initial_state
    elif isinstance(initial_state, comando.Parameter):
        init = comando.Variable(f'{state}__init')
        init.fix(initial_state.value)
    elif initial_state.is_Number:
        init = comando.Variable(f'{state}__init')
        init.fix(initial_state)
    else:
        raise NotImplementedError(f'The initial state of {state} must be '
                                  'a Variable, Parameter or Number, but '
                                  f'is a {type(initial_state)}!')
    prev = comando.VariableVector(state.name + '__prev')
    try:
        prev_data = state.expansion.groupby('s').shift(1, fill_value=init)
    except KeyError:  # no scenarios
        prev_data = state.expansion.shift(1, fill_value=init)
    prev.expansion = prev_data
    return prev


def split(iterable, condition, c_type=None):
    """Sort `iterable`'s elements into two containers based on `condition`.

    Arguments
    ---------
    iterable : iterable
        iterable to be sorted
    condition : callable
        function used for sorting, return value must convert to bool
    c_type : type
        type of the resulting containers; if unset will attempt to use
        `type(iterable)` to create the containers or default to `list` if
        that fails. If after this both `isinstance` and `c_type()` are
        mappings, the new containers will be submappings of `iterable`
        and `condition` is only applied to the values of `iterable`.

    Returns
    -------
    cond_true, cond_false : tuple
        tuple containing:
        - the container of elements that do **not** satisfy `condition`
        - the container of elements that do satisfy `condition`
    """
    if c_type is None:
        c_type = type(iterable)
    try:  # Try to deduce the container type from iterable
        res = [c_type(), c_type()]
    except TypeError:
        c_type = list
        res = [c_type(), c_type()]

    if c_type is list:
        for x in iterable:
            res[condition(x)].append(x)
    elif c_type is set:
        for x in iterable:
            res[condition(x)].add(x)
    elif c_type is tuple:
        for x in iterable:
            res[condition(x)] += (x, )
    elif c_type is dict:
        if isinstance(iterable, dict):
            for k, v in iterable.items():
                res[condition(v)][k] = v
        else:
            for x in iterable:
                res[condition(x)][x[0]] = x[1]
    else:
        raise NotImplementedError(f'No insertion known for {c_type}!')
    return res[False], res[True]


# NOTE: Precedence is already handled by the tree datastructures in the object
#       oriented setting, but it becomes important when writing to strings!
#       Thats where we will need to be able to handle the '()' entry in the
#       op_map. For the object oriented case, '()' simply maps to identity.
def identity(expr):
    """Return expr."""
    return expr


def _sum(*args):
    """Return the sum of the elements in args."""
    return sum(args)
    return numpy.sum(args)


def prod(*args):
    """Return the product of the elements in args."""
    return reduce(operator.mul, args)
    return numpy.prod(args)


# def pow(*args):
#     """Return the result of args[0] to power of args[1]."""
#     return 1/args[0] if args[1] == -1 else args[0] ** args[1]

# comando_op_map = {
#     '()': identity,
#     'Add': _sum,
#     'Neg': lambda arg: -arg,
#     'Mul': prod,
#     'Div': operator.truediv,
#     'Pow': operator.pow,
#     'Inv': partial(operator.truediv, 1),
#     'LessThan': operator.le,
#     'GreaterThan': operator.ge,
#     'Equality': operator.eq}

numpy_op_map = {
    '()': identity,
    'Add': _sum,
    'Neg': lambda arg: -arg,
    'Mul': prod,
    'Div': operator.truediv,
    'Pow': operator.pow,
    'Inv': partial(operator.truediv, 1),
    'LessThan': operator.le,
    'GreaterThan': operator.ge,
    'Equality': operator.eq,
    'Min': lambda *args: reduce(minimum, args),
    'Max': lambda *args: reduce(maximum, args),
    'Abs': numpy.abs}

# exponential_functions = {'exp', 'log'}
# nonsmooth_function_names = {'sign', 'Abs', 'Min', 'Max', 'ceiling', 'floor'}
# trigonometric_function_names = {'sin', 'cos', 'tan', 'cot', 'sec', 'csc'}
# trigonometric_inverse_function_names = {'asin', 'acos', 'atan', 'acot', 'asec',
#                                         'acsc'}
# hyperbolic_function_names = {'sinh', 'cosh', 'tanh', 'coth', 'sech', 'csch'}
# hyperbolic_inverse_function_names = {'asinh', 'acosh', 'atanh', 'acoth',
#                                      'asech', 'acsch'}
# comando_functions = set().union(
#     exponential_functions, nonsmooth_function_names,
#     trigonometric_function_names, trigonometric_inverse_function_names,
#     hyperbolic_function_names, hyperbolic_inverse_function_names
# )
# NOTE: All of these functions depend on the backend and therefore need to be
#       looked up dynamically to allow for backend switches!
# def _inject_funcdefs():
#     for f_name in comando_functions:
#         comando_op_map[f_name] = \
#             lambda arg: getattr(comando.get_backend(), f_name)(arg)
# _inject_funcdefs()

for name in 'exp', 'log', 'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh':
    numpy_op_map[name] = getattr(numpy, name)
for name in {'asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh'}:
    numpy_op_map[name] = getattr(numpy, f'arc{name[1:]}')


def make_str_func(name):
    """Create a default string representation of functions.

    Arguments
    =========
    name : str
        name of the function to be represented

    Returns
    =======
    func : callable
        function that returns a string representation for the function with
        given string arguments
    """
    return lambda *args: f"{name.lower()}({', '.join(args)})"


class DefaultStringMap(defaultdict):
    """A default dict implementation whose factory is called with the key."""

    def __missing__(self, key):
        """Call the factory with the missing key."""
        if self.default_factory is None:
            raise KeyError((key, ))
        self[key] = value = self.default_factory(key)
        return value


_str_map = DefaultStringMap(make_str_func, {
    '()': lambda expr: f'({expr})',
    'Add': lambda *args: " + ".join(args),
    'Sub': lambda minu, subtr: f'{minu} - {subtr}',
    'Neg': lambda arg: f'-{arg}',
    'Mul': lambda *args: " * ".join(args),  #
    'Div': lambda numer, denom: f'{numer} / {denom}',
    'Inv': lambda arg: f'1 / {arg}',
    'LessThan': lambda lhs, rhs: f'{lhs} <= {rhs}',
    'GreaterThan': lambda lhs, rhs: f'{lhs} >= {rhs}',
    'Equality': lambda lhs, rhs: f'{lhs} == {rhs}',
    'Pow': lambda base, exponent: f'{base} ** {exponent}',
    #'exp': lambda arg: f'exp({arg})',
    #'log': lambda arg: f'log({arg})',
    #'tanh': lambda arg: f'(1 - 2/(exp(2 * ({arg})) + 1))'
}
)
# for name in 'exp', 'log', 'sqrt', 'tan', 'sinh', 'tanh':
#     numpy_op_map[name] = getattr(numpy, name)
# for name in 'asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh':
#     numpy_op_map[name] = getattr(numpy, f'arc{name[1:]}')


def get_index(expr):
    """Get the index of the current expression or None if it is not indexed."""
    dummy = sum((sym.value[:] * 0 for sym in expr.free_symbols if sym.indexed),
                0)
    return None if isinstance(dummy, int) else dummy.index


def parse(expr, sym_map=None, op_map=None, idx=None, const_fun=float):
    """Parse the sympy expression `expr`, replacing symbols and operators.

    The `sym_map` provides the possibility to replace symbols (Parameters and
    Variables) with new objects, if the symbol map is empty or there is no
    mapping for a given symbol, it is replaced with its `value` attrubute; if
    `value` is `None`, the symbol is not replaced.
    The `op_map` provides the possibility to replace sympy's operator objects
    with other operators. If no op map is provided, the default python
    operators are used.
    """
    if sym_map is None:
        if op_map is None and idx is None and const_fun is float:
            return expr  # Quick exit, nothing changes
        sym_map = {}

    for sym in expr.free_symbols:  # Complete the sym map with missing symbols
        if sym not in sym_map:
            sym_map[sym] = sym

    # Based on default or user-provided op_map:
    if op_map is None:
        op_map = comando.op_map
    else:
        op_map = {**comando.base_op_map, **op_map}  # override default with user ops

    if idx is None:
        return _parse(expr, sym_map, op_map, const_fun)
    return _idx_parse(expr, sym_map, op_map, idx, const_fun)


def _parse(expr, sym_map, op_map, const_fun):
    """Parse expr, all symbols/operators must be in sym_map/op_map."""
    if expr.args:
        try:
            f_name = get_type_name(expr)
            op = op_map[f_name]
        except KeyError:
            raise NotImplementedError(f"The function '{f_name}' is not defined"
                                      " in the operational map of the "
                                      "interface you are using. Either provide"
                                      " your own definition or reformulate.")
        args = tuple(_parse(arg, sym_map, op_map, const_fun)
                     for arg in expr.args)
        return op(*args)

    if expr.is_Number:
        return const_fun(expr)
    return sym_map[expr]


def _idx_parse(expr, sym_map, op_map, idx, const_fun):
    """Parse expr at idx, all symbols/operators must be in sym_map/op_map."""
    if expr.args:
        try:
            f_name = get_type_name(expr)
            op = op_map[f_name]
        except KeyError:
            raise NotImplementedError(f"The function '{f_name}' is not defined"
                                      " in the operational map of the "
                                      "interface you are using. Either provide"
                                      " your own definition or reformulate.")
        args = tuple(_idx_parse(arg, sym_map, op_map, idx, const_fun)
                     for arg in expr.args)
        return op(*args)

    if expr.is_Number:
        return const_fun(expr)
    rep = sym_map[expr]

    if indexed(expr):
        return rep[idx]
    else:
        return rep


def evaluate(expr, idx=None):
    """Evaluate the given expression at the symbols' current values."""
    index = get_index(expr)

    if index is None:  # scalar expression
        return _parse(expr, {sym: sym.value for sym in expr.free_symbols},
                      numpy_op_map, float)
    # indexed expression
    if idx is None:  # No particular index desired
        import pandas
        res = _parse(expr, {sym: sym.value.values if sym.indexed else sym.value
                            for sym in expr.free_symbols}, numpy_op_map, float)
        return pandas.Series(res, index, float)
    # only values for index=idx desired
    return float(expr.subs({sym: sym[idx].value if sym.indexed else sym.value
                            for sym in expr.free_symbols}))
    # TODO: Since we assume unique index elements, this should only come into
    #       play when we introduce slicing, e.g. idx is slice(None)
    # if len(res) > 1:
    #    index = index[index == idx]
    # ...
    # res = _idx_parse(expr, {sym: sym[idx].value if sym.indexed else sym.value
    #                         for sym in expr.free_symbols}, numpy_op_map, idx,
    #                  float)
    # return Series(res, index, float)


def precedence(expr):
    """Get the operator precedence of the expression.

    If the argument of an expression has a higher precedence than the
    expression itself, precedence needs to be represented explicitly in string
    representations. An example is (a + b) * c, where dropping the parentheses
    changes the result. We define precedences as:

        ±x, f(x): 0
        x ** y: 1
        x * y, x/y: 2
        x ± y: 3
    """
    return {comando.Add: 3, comando.Mul: 2, comando.Pow: 1}.get(type(expr), 0)


def as_numer_denom(expr):
    """Normalize as_numer_denom for sympy and symengine.

    Also see:
    https://github.com/symengine/symengine/issues/681#issuecomment-724039779
    """
    numer, denom = expr.as_numer_denom()
    return numer.subs(1., 1), denom.subs(1., 1)


def get_type_name(expr):
    """Get the type name of an expression."""
    f_name = type(expr).__name__
    if f_name == 'FunctionSymbol':  # User-defined functions in symengine
        f_name = expr.get_name()
    return f_name


# TODO: We can save a lot of time for indexed expressions if we parse the non
#       changing part (which for large expressions will be the majority of the
#       expression) only once. This would result in a string template instead
#       of a string. The template arguments can then be replaced with the
#       string corresponding to the index. This principle could also be nested!
class StrParser:
    """A class for creating interface-specific string representations.

    Callbacks must have a signature consisting of three arguments: The parser
    the expression to be parsed and an optional index
    Arguments:
    =========

    sym_map
    str_map
    add_callback
    mul_callback
    pow_callback
    """

    def __init__(self, sym_map, str_map=None, override_default_str_map=True,
                 add_callback=None, mul_callback=None, pow_callback=None):
        self.sym_map = sym_map

        # Based on default or user-provided str_map:
        if str_map is None:
            self.str_map = _str_map.copy()
        elif override_default_str_map:
            self.str_map = _str_map.copy()
            self.str_map.update(str_map)
        else:
            self.str_map = str_map

        if add_callback is None:
            self.parse_add = self._parse_add
        else:
            def parse_add(expr, idx):
                res = add_callback(self, expr, idx)
                if res is None:
                    return self._parse_add(expr, idx)
                return res
            self.parse_add = parse_add

        if mul_callback is None:
            self.parse_mul = self._parse_mul
        else:
            def parse_mul(expr, idx):
                res = mul_callback(self, expr, idx)
                if res is None:
                    return self._parse_mul(expr, idx)
                return res
            self.parse_mul = parse_mul

        if pow_callback is None:
            self.parse_pow = self._parse_pow
        else:
            def parse_pow(expr, idx):
                res = pow_callback(self, expr, idx)
                if res is None:
                    return self._parse_pow(expr, idx)
                return res
            self.parse_pow = parse_pow
        self.cache = {}

    def _parse_add(self, expr, idx):
        p_args, n_args = split(expr.args, is_negated)
        if p_args:
            if n_args:  # (a -b + c - d) -> a + c - (b + d)
                m_expr = expr.func(*p_args)  # minuend
                s_expr = expr.func(*(-a for a in n_args))  # subtrahend
                # While the minuend can be handled normally, the subtrahend
                # must be parenthesized
                args = (*self.parse_args((m_expr, ), idx, 3),
                        *self.parse_args((s_expr, ), idx, 2.5))
                return self.str_map['Sub'](*args)
            return self.str_map['Add'](*self.parse_args(expr.args, idx, 3))
        # This handles -a -b -> -(a + b)
        # NOTE: We cannot just use -expr as this can result in an infinite
        #       recursion with the 'Neg' cases in expr.is_Mul!
        neg_expr = expr.func(*(-arg for arg in n_args))
        return self.str_map['Neg'](*self.parse_args((neg_expr, ), idx, 2))

    def _parse_mul(self, expr, idx):
        a0 = expr.args[0]
        if a0.is_Number and a0 < 0:
            return self.str_map['Neg'](*self.parse_args((-expr, ), idx, 2))
        numer, denom = as_numer_denom(expr)
        if denom == 1 or denom == 1.0:
            return self.str_map['Mul'](*self.parse_args(expr.args, idx, 2))
        # In contrast to sympy, we don't want to distinguish between cases
        # like -(a + b) / c, -( (a + b) / c) and (a + b)/-c
        if is_negated(numer):
            return self.str_map['Neg'](*self.parse_args((-numer/denom, ),
                                       idx, 2))
        if is_negated(denom):
            return self.str_map['Neg'](*self.parse_args((numer/-denom, ),
                                       idx, 2))
        # While the numerator can be handled normally, the denominator
        # must be parenthesized if it's equivalent to a multiplication
        # to correctly represent forms like a / (b * c)
        args = (*self.parse_args((numer, ), idx, 2),
                *self.parse_args((denom, ), idx, 1.5))
        return self.str_map['Div'](*args)

    def _parse_pow(self, expr, idx):
        base, exponent = expr.args
        if is_negated(exponent):
            return self.str_map['Inv'](*self.parse_args((base ** -exponent, ),
                                                        idx, 1.5))
        return self.str_map['Pow'](*self.parse_args(expr.args, idx, 1))

    def __call__(self, expr, idx=None):
        """Execute the StrParser."""
        if expr.is_Add:
            return self.parse_add(expr, idx)
        if expr.is_Mul:
            return self.parse_mul(expr, idx)
        if expr.is_Pow:
            return self.parse_pow(expr, idx)
        # prec == 0
        if expr.args:  # A univariate function
            return self.str_map[get_type_name(expr)](
                *self.parse_args(expr.args, idx)
            )
        if expr.is_Symbol:  # Variables or Parameters
            replacement = self.sym_map.get(expr, None)
            if replacement is None:
                replacement = expr.value
            if idx is not None and expr.indexed:
                replacement = replacement[idx]
            return str(replacement)  # replacement specified
        # Number
        if expr.is_rational and not expr.is_integer:
            # we parenthesize rationals to avoid problems with, e.g., a * b
            # where b is a rational
            return self.str_map['()'](str(expr))
        return str(expr)

    def parse_args(self, args, idx, parent_prec=0):
        """Parse arguments of an expression with parentheses as appropriate."""
        if parent_prec == 0:
            return tuple(self(arg, idx) for arg in args)
        else:
            return tuple(self.str_map['()'](self(arg, idx))
                         if precedence(arg) > parent_prec else
                         self(arg, idx) for arg in args)

    def cached_parse(self, expr, idx=None):
        """Look up previously parsed expressions in self.cache."""
        # try:
        if (expr, idx) in self.cache:
            return self.cache[expr, idx]
        else:
            res = self(expr, idx)
            self.cache[expr, idx] = res
            return res
        # except RecursionError as e:
        #     if not hasattr(e, 'expr'):
        #         e.expr = expr
        #     raise


def str_parse(expr, sym_map=None, str_map=None, idx=None):
    """Parse the sympy expression `expr`, replacing symbols and operators.

    The `sym_map` provides the possibility to replace symbols (Parameters and
    Variables) with their string representations, if the symbol map is empty
    or there is no mapping for a given symbol, it is replaced with its `value`
    attribute.
    The `str_map` provides the possibility to replace sympy's operator objects
    with string representations.
    """
    # DEBUG
    # try:
    # DEBUG
    if sym_map is None:
        sym_map = {}

    # Based on default or user-provided str_map:
    if str_map is None:
        str_map = _str_map.copy()
    else:  # override default with user map
        __str_map = _str_map.copy()
        __str_map.update(str_map)
        str_map = __str_map

    def rec_str_parse(expr):
        """Recursively call to str_parse."""
        return str_parse(expr, sym_map, str_map, idx)

    def parse_args(args, parent_prec=0):
        """Parse arguments of an expression with parentheses as appropriate."""
        if parent_prec == 0:
            return tuple(str_parse(arg, sym_map, str_map, idx) for arg in args)
        else:
            return tuple(str_map['()'](str_parse(arg, sym_map, str_map, idx))
                         if precedence(arg) > parent_prec else
                         str_parse(arg, sym_map, str_map, idx) for arg in args)

    if expr.is_Add:
        p_args, n_args = split(expr.args, is_negated)
        if p_args:
            if n_args:  # (a -b + c - d) -> a + c - (b + d)
                m_expr = expr.func(*p_args)  # minuend
                s_expr = expr.func(*(-a for a in n_args))  # subtrahend
                # While the minuend can be handled normally, the subtrahend
                # must be parenthesized
                args = (*parse_args((m_expr, ), 3),
                        *parse_args((s_expr, ), 2.5))
                return str_map['Sub'](*args)
            return str_map['Add'](*parse_args(expr.args, 3))
        # This handles -a -b -> -(a + b)
        # NOTE: We cannot just use -expr as this can result in an infinite
        #       recursion with the 'Neg' cases in expr.is_Mul!
        neg_expr = expr.func(*(-arg for arg in n_args))
        return str_map['Neg'](*parse_args((neg_expr, ), 2))
    if expr.is_Mul:
        a0 = expr.args[0]
        if a0.is_Number and a0 < 0:
            return str_map['Neg'](*parse_args((-expr, ), 2))
        numer, denom = as_numer_denom(expr)
        if denom == 1 or denom == 1.0:
            return str_map['Mul'](*parse_args(expr.args, 2))
        # In contrast to sympy, we don't want to distinguish between cases
        # like -(a + b) / c, -( (a + b) / c) and (a + b)/-c
        if is_negated(numer):
            return str_map['Neg'](*parse_args((-numer/denom, ), 2))
        if is_negated(denom):
            return str_map['Neg'](*parse_args((numer/-denom, ), 2))
        # While the numerator can be handled normally, the denominator
        # must be parenthesized if it's equivalent to a multiplication
        # to correctly represent forms like a / (b * c)
        args = (*parse_args((numer, ), 2), *parse_args((denom, ), 1.5))
        return str_map['Div'](*args)
    if expr.is_Pow:
        base, exponent = expr.args
        if is_negated(exponent):
            return str_map['Inv'](*parse_args((base ** -exponent, ), 1.5))
        return str_map['Pow'](*parse_args(expr.args, 1))
    # prec == 0
    if expr.args:  # A univariate function
        return str_map[get_type_name(expr)](*parse_args(expr.args))
    if expr.is_Symbol:  # Variables or Parameters
        replacement = sym_map.get(expr, None)
        if replacement is None:
            replacement = expr.value
        if idx is not None and expr.indexed:
            replacement = replacement[idx]
        return str(replacement)  # replacement specified
    # Number
    if expr.is_rational and not expr.is_integer:
        # we parenthesize rationals to avoid problems with, e.g., a * b where
        # b is a rational
        return str_map['()'](str(expr))
    return str(expr)

    # DEBUG
    # except RecursionError as e:
    #     if not hasattr(e, 'expr'):
    #         e.expr = expr
    #     raise
    # DEBUG


# Method for tighter bounds of univariate continuous functions

class RootFindingError(Exception): pass


def cont_univ_bounds(expr, var):
    """Compute tight bounds for bounded continuous univariate expressions.

    A bounded continuous expression of a single variable can be bounded exactly
    by looking at its function values at 'critical points', i.e., the union of
    its bounds and the points at which its derivative is zero that are in the
    interval implied by the bounds.
    """
    from sympy.solvers import solveset
    from sympy import Interval
    # lb = var.lb if var.lb is not None else float('-inf')
    # ub = var.ub if var.ub is not None else float('inf')
    lb, ub = bounds(var)
    crit_points = {lb, ub}
    try:
        crit_points.update(solveset(expr.diff(var), var, Interval(lb, ub)))
    except Exception:
        raise RootFindingError
    crit_values = [expr.subs(var, point) for point in crit_points]
    return min(crit_values), max(crit_values)


# Functions for interval arithmetic
# NOTE: The below definitions assume that for any bound interval (lb, ub) with
#       -∞ < lb == ub < ∞. In any case we should probably warn or throw an
#       error if there are any infinite bounds.
def sum_bounds(*bounds):
    """Return the sum of a sequence of intervals."""
    _lb, _ub = 0, 0
    for lb, ub in bounds:
        _lb += lb
        _ub += ub
    return (_lb, _ub)


def mul_bounds(a, b):
    """Return the result of a * b by interval arithmetic."""
    # NOTE: If a and b come from the same expression, such as in the case
    #       x * x, the bound obtained by this function is larger than
    #       necessary as the multiplication is equivalent to x ** 2 and in
    #       general [lb, ub] ** 2 ⊆ [lb, ub] * [lb, ub] holds! Luckily
    #       sympy does a fairly good job of simplifying expressions, so
    #       this problem should not occur too often.
    # NOTE: Edge cases of 0 * ±inf are treated as 0 because infinite bounds
    #       should never be binding, otherwise the solution would be unbounded
    inf = float('inf')
    S = [i * j if {i, j} not in [{0, -inf}, {0, inf}] else 0
         for i in a for j in b]  # all combinations of products from a and b
    return (min(S), max(S))


def prod_bounds(*bounds):
    """Return the product of a sequence of intervals."""
    return reduce(mul_bounds, bounds)


def pow_bounds(a, b):
    """Return the result of a ** b if a[0] == a[1] or b[0] == b[1].

    The power function is monotonous on any domain that excludes zero!
    We therefore only need to consider special cases explicitly and can
    return the sorted tuple of powers of the base boundaries otherwise.
    """
    # NOTE: Infinities give rise to another special case, but currently we're
    #       not interested in this, so we save ourselves the headache...
    if b[0] in [None, -float('inf'), float('inf')] or b[0] != b[1]:
        if a[0] == a[1]:
            a = a[0]
            if a < -1:   # increasing and only valid for integral b
                # vals = integer_values in [b[0], b[1]]
                # # will be the last two elements of vals
                # return (min(vals), max(vals))
                raise NotImplementedError
            elif a < 0:  # decreasing and only valid for integral b
                # vals = integer_values in [b[0], b[1]]
                # # will be the first two elements of vals
                # return (min(vals), max(vals))
                raise NotImplementedError
            elif a > 1:   # decreasing
                return (a ** b[1], a ** b[0])
            elif a > 0:  # 0 < a <= 1
                return (a ** b[0], a ** b[1])
            else:  # a == 0
                return (0, 0)
        else:
            raise ValueError("Interval power requires real-valued singleton "
                             f"exponents (with -∞ < lb == ub < ∞), got {b}")
    b = float(b[0])
    if b == 0:
        return (1, 1)  # x ** 0 == 1 for any float x, including 0, ±∞ and nan!
    lb, ub = a
    if (lb <= 0) and (ub >= 0):  # 0 is contained in base set a
        if b < 0:  # negative exponent with 0 base -> Division by zero!
            raise ValueError("Base of 0 in power with negative exponent")
        if not b % 2:  # b is an even power -> a ** b is symmetric around 0!
            return (0, max([lb ** b, ub ** b]))
    if lb < 0 and (b % 1) != 0:  # b is not integral -> a ** b is complex!
        raise ValueError("Negative base with a non-integer exponent")
    return tuple(sorted([lb ** b, ub ** b]))


def bLe(x, y):
    """Return the truth value of x <= y, if it can be determined, else None."""
    if x[1] <= y[0]:
        return True
    if x[0] > y[1]:
        return False


def bGe(x, y):
    """Return the truth value of x >= y, if it can be determined, else None."""
    if x[0] >= y[1]:
        return True
    if x[1] < y[0]:
        return False


def bEq(x, y):
    """Return the truth value of x == y, if it can be determined, else None."""
    if x[0] == y[0] == x[1] == y[1]:
        return True
    if x[1] <= y[0] or y[1] <= x[0]:
        return False


def bounds_cost_turton(x, c1, c2, c3):
    """Compute bounds for the the Guthrie cost function.

    cost_turton(x, c1, c2, c3) \
        = pow(10, c1 + c2 * log10(x) + c3 * log10(x) ** 2)
    """
    from numpy import log10
    # bounds for log10(x)
    log_bounds = log10(x[0]), log10(x[1])
    # bounds for c2 * log10(x)
    c2_bounds = mul_bounds(c2, log_bounds)
    # bounds for log10(x) ** 2
    log_squared_bounds = tuple(sorted([log_bounds[0] ** 2,
                                       log_bounds[1] ** 2]))
    # bounds for c3 * log10(x) ** 2
    c3_bounds = mul_bounds(c3, log_squared_bounds)
    # bounds for c1 + c2 * log10(x) + c3 * log10(x) ** 2
    exponent_bounds = sum_bounds(c1, c2_bounds, c3_bounds)
    return (10 ** exponent_bounds[1], 10 ** exponent_bounds[0])


def lmtd_bounds(dT1, dT2):
    """Calculate the logarithmic mean temperature difference.

    lmtd(dT1, dT2) = (dT1 - dT2) / log(dT1 / dT2)
    """
    from numpy import log
    # assuming dT2 > 0!
    inv_dT2_bounds = (1/dT2[1], 1/dT2[0])
    # bounds for dT1 / dT2
    arg_bounds = mul_bounds(dT1, inv_dT2_bounds)
    # bounds for log(dT1 / dT2) and its inverse
    log_bounds = log(arg_bounds[0]), log(arg_bounds[1])
    inv_log_bounds = (1/log_bounds[1], 1/log_bounds[0])
    # bounds for dT1 - dT2
    diff_bounds = sum_bounds(dT1, (-dT2[1], -dT2[0]))
    return mul_bounds(inv_log_bounds, diff_bounds)


def rlmtd_bounds(dT1, dT2):
    """Calculate bounds for the inverse of the LMTD.

    rlmtd(dT1, dT2) = log(dT1 / dT2) / (dT1 - dT2)
    """
    from numpy import log
    # assuming dT2 > 0!
    inv_dT2_bounds = (1/dT2[1], 1/dT2[0])
    # bounds for dT1 / dT2
    arg_bounds = mul_bounds(dT1, inv_dT2_bounds)
    # bounds for log(dT1 / dT2)
    log_bounds = log(arg_bounds[0]), log(arg_bounds[1])
    # bounds for dT1 - dT2 and its inverse
    diff_bounds = sum_bounds(dT1, (-dT2[1], -dT2[0]))
    inv_diff_bounds = (1/diff_bounds[1], 1/diff_bounds[0])
    return mul_bounds(log_bounds, inv_diff_bounds)


def floor_substitute_bounds(x, LB, UB):
    from numpy import floor
    return floor(LB[0]), floor(UB[1])


bounds_op_map = {
    'Add': sum_bounds,
    'Mul': prod_bounds,
    'Pow': pow_bounds,
    'LessThan':  bLe,
    'GreaterThan': bGe,
    'Equality': bEq,
    'cost_turton': bounds_cost_turton,
    'lmtd': lmtd_bounds,
    'rlmtd': rlmtd_bounds,
    'floor_substitute': floor_substitute_bounds
    }


# EXPERIMENTAL!!!
def _cont_univ_bounds(func, lb, ub):
    import comando
    from sympy.solvers import solveset
    from sympy import Interval
    DUMMY = comando.Symbol('DUMMY')
    f = func(DUMMY)
    crit_points = {lb, ub}
    try:
        crit_points.update(solveset(f.diff(DUMMY), DUMMY, Interval(lb, ub)))
    except Exception:
        raise RootFindingError
    crit_values = [f.subs(DUMMY, point) for point in crit_points]
    return min(crit_values), max(crit_values)


def _mon_inc(name, b):
    func = getattr(math, name)
    return func(b[0]), func(b[1])


# monotonously increasing functions
for name in ['exp', 'log', 'sqrt', 'tan', 'sinh', 'tanh',
             'asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh']:
    bounds_op_map[name] = partial(_mon_inc, name)

# symmetric functions
# bounds_op_map['cosh'] = ...


# for name in ['sin', 'cos']:
#     func = getattr(comando, name)
#     bounds_op_map[name] = lambda bounds: _cont_univ_bounds(func, *bounds)

bounds_op_map['Min'] = lambda *b: (min(bi[0] for bi in b),
                                   min(bi[1] for bi in b))
bounds_op_map['Max'] = lambda *b: (max(bi[0] for bi in b),
                                   max(bi[1] for bi in b))


def _normalize_bounds(bounds):
    lb = float('-inf') if math.isnan(bounds[0]) else bounds[0]
    ub = float('inf') if math.isnan(bounds[1]) else bounds[1]
    return(lb, ub)


def bounds(expr):
    """Propagate variable bounds through the given expression."""
    from pandas import Series
    index = get_index(expr)
    if index is None:
        return _bounds(expr)
    lb = Series(-comando.INF, index)
    ub = Series(comando.INF, index)
    for i in index:
        lb[i], ub[i] = _bounds(parse(expr, idx=i))
    return lb, ub


def _bounds(expr):
    """Propagate variable bounds through the given expression."""
    if expr.args:
        op = bounds_op_map[get_type_name(expr)]
        args = tuple(_bounds(arg) for arg in expr.args)
        try:
            return op(*args)
        except ValueError as e:
            raise ValueError(f'{e}: {expr}')

    # if expr.is_Number:  false for symengine.E
    try:
        return float(expr), float(expr)
    except (TypeError, RuntimeError):
        pass

    try:  # treating expr as a Variable by getting its bounds
        return expr.bounds
    except AttributeError:
        # expr is a Parameter it has a value attribute...
        val = expr.value
        return _normalize_bounds((val, val))


def make_tac_objective(system, ic_labels=None, fc_labels=None, vc_labels=None,
                       n=10, i=0.08):
    """Create the objective components corresponding to total annualized cost.

    Arguments
    ---------
    ic_labels : iterable of `str`
        labels for the investment cost expressions
    fc_labels : iterable of `str`
        labels for the fixed cost expressions
    vc_labels : iterable of `str`
        labels for the variable cost expressions
    n : `int`
        number of repetitions of the time-period specified by `timesteps`.
    i : `float`
        interest_rate

    Returns
    -------
    Expressions for design_objective and operational_objective forming part
    of the total annualized costs.
    """
    if ic_labels is None:
        ic_labels = ['investment_costs']
    if fc_labels is None:
        fc_labels = ['fixed_costs']
    if vc_labels is None:
        vc_labels = ['variable_costs']
    ic = sum(system.aggregate_component_expressions(ic_label)
             for ic_label in ic_labels)
    fc = sum(system.aggregate_component_expressions(fc_label)
             for fc_label in fc_labels)
    vc = sum(system.aggregate_component_expressions(vc_label)
             for vc_label in vc_labels)
    af = ((1 + i) ** n * i) / ((1 + i) ** n - 1)  # annuity factor
    return af*ic + fc, vc


def make_mayer_objective(system):
    """Create a Mayer term objective.

    An option for dynamic optimization
    phi(t_f) is minimized and phi_dot = f(operation)
    """
    phi_dot = \
        system.aggregate_component_expressions('variable Mayer objective')

    system.objective_type = 'Mayer term objective'
    # no design objective, just operational part
    return phi_dot


def lambdify(expr, vars=None, eval_params=True):
    """Create a function for evaluation of expr with numpy.

    A variant of sympy's lambdify which creates a function for fast numerical
    evaluation of an expression.

    Arguments
    ---------

    expr : Expression
        the expression to be turned into a function

    eval_params : bool
        Whether to replace parameters by their values (if expr contains indexed
        parameters, evaluating the function will thus return a Series)
        default: True

    Returns
    -------
    f : callable
        a function returning the value of the expression for given variable
        values
    """
    if vars is None:
        vars = get_vars(expr)
    from collections import OrderedDict
    v_reps = OrderedDict((v, sympy.Symbol(f'VALUE_OF_{v.name}')) for v in vars)
    p_reps = {}
    if eval_params:
        p_reps = {p: p.value for p in get_pars(expr)}
        if p_reps and get_index(sum(p_reps)) is not None:  # multiple expressions
            from pandas import Series
            exprs = parse(expr, {**p_reps, **v_reps})
            funcs = Series((sympy.lambdify(v_reps.values(), ei)
                            for ei in exprs), exprs.index, 'O')

            def func(*args):
                return funcs.apply(lambda f: f(*args))

            return func
    return sympy.lambdify(v_reps.values(), parse(expr, {**p_reps, **v_reps}))


def define_function(name, implementation):
    """Define a new symbolic function with a name and an implementation.

    The given implementation will be used when the function is called with
    nonsymbolic arguments only.

    Arguments
    ---------
    name : str
        the function's name
    implementation : callable
        a callable that serves as the numerical implementation and is used when
        no symbols are within the arguments. Must have a fixed number of
        positional arguments.

    Returns
    -------
    new_function : sympy.FunctionClass
        the newly defined function
    """

    from inspect import signature
    sig = signature(implementation)
    if any(arg.kind.value > 2 for arg in sig.parameters.values()):
        raise ValueError('Arguments of implementation must be either '
                         'POSITIONAL_ONLY, POSITIONAL_OR_KEYWORD or '
                         'VAR_POSITIONAL!')

    ###########################################################################
    # Needed for sympy evaluation
    @classmethod
    def eval(cls, *args):
        for arg in args:
            try:
                if any(arg.free_symbols):
                    break
            except AttributeError:  # e.g. pandas.Series
                pass
        else:
            return implementation(*args)

    new_func = comando.Function(name)
    new_func.eval = eval
    new_func.n_args = len(sig.parameters)

    # allow evaluation
    comando.utility.numpy_op_map[name] = implementation
    ###########################################################################
    # Needed for symengine evaluation
    if not hasattr(comando.Function, 'implementations'):
        comando.Function.implementations = {}
    comando.Function.implementations[name] = implementation

    def _wrapped(*args):
        # print('wrapped call')
        new_func_instance = new_func(*args)

        def __float__(self):
            # print(f'wrapped float of {self.get_name()}')
            implementation = \
                comando.Function.implementations[self.get_name()]
            # return float(implementation(*self.args))
            return float(implementation(*map(float, self.args)))
        type(new_func_instance).__float__ = __float__

        return new_func_instance

    comando.op_map[name] = _wrapped
    return _wrapped


def make_function(expr):
    """Create a callable function from the given expression.

    The user can either specify values for all Variables in the order specified
    by the returned function's '__doc__' (in alphabetical order), or pass
    keyword arguments in the form of identifier-value pairs. The values are
    taken as alternatives for the default values used in the evaluation for the
    Variables or Parameters with matching identifiers.
    """
    syms = expr.free_symbols
    pars = sorted((sym for sym in syms if isinstance(sym, comando.Parameter)),
                  key=lambda x: x.name)
    vars = sorted((sym for sym in syms if isinstance(sym, comando.Variable)),
                  key=lambda x: x.name)

    def function(*user_var_values, **user_values):
        nargs = len(vars)
        if user_var_values:
            if len(user_var_values) != nargs:
                raise ValueError(f"0 or {nargs} arguments required!")
            sym_map = {var: user_var_values[i] for i, var in enumerate(vars)}
        else:
            sym_map = {var: var.value for var in vars}
        par_values = {par: par.value for par in pars}
        sym_map.update(par_values)

        for id, val in user_values.items():
            known_symbols = iter(sym_map)
            while True:
                try:
                    sym = next(known_symbols)
                except StopIteration:
                    names = [sym.name for sym in sym_map]
                    raise ValueError(f'Expression only has symbols {names}!')
                if sym.name == id:
                    break
            sym_map[sym] = val

        par_values.update(user_values)

        return parse(expr, sym_map)

    # Manually passing docstring to specify argument order
    doc = f"Call {expr} providing {vars} and optionally {pars} as arguments."
    function.__doc__ = doc

    return function

# NOTE: We can now do things like the following
# my_pars = [Parameter(f"p{i}", value=i**2) for i in range(3)]
# my_vars = [Variable(f"v{i}", value=10-i) for i in range(2)]
# expr = my_vars[0] ** my_vars[1] + sum(my_pars)
# f = make_function(expr)
# f.__doc__
# f(1, 2)
#
# from pandas import Series
# s = Series([1,3,5])
# f(s, 2)
