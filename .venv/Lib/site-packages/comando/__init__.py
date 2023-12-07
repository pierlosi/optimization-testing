"""Configuration of the backends for COMANDO."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import os
from enum import Enum
from collections.abc import Iterable, Mapping
import operator
from functools import reduce, partial

from pandas import Series, DataFrame

os.environ['SYMPY_CACHE_SIZE'] = 'None'
import sympy
from sympy.logic.boolalg import Boolean


sympy.Boolean = Boolean
# sympy.compatibility.NotIterable = NotIterable
backends = {'sympy': sympy}
# BACKEND = symengine
BACKEND = sympy


class Domain(Enum):
    """Simple Enum for variable domains, specify other types via bounds."""

    REAL = 1
    INTEGER = 2
    BINARY = 3


REAL, INTEGER, BINARY = Domain.REAL, Domain.INTEGER, Domain.BINARY
NAN = float('nan')
INF = float('inf')
UNBOUNDED = {-INF, INF}
EPS = 1e-9


def identity(expr):
    """Return expr."""
    return expr


def _sum(*args):
    """Return the sum of the elements in args."""
    return sum(args)


def prod(*args):
    """Return the product of the elements in args."""
    return reduce(operator.mul, args)


base_op_map = {
    '()': identity,
    'Add': _sum,
    'Neg': lambda arg: -arg,
    'Mul': prod,
    'Div': operator.truediv,
    'Pow': operator.pow,
    'Inv': partial(operator.truediv, 1),
    'LessThan': operator.le,
    'GreaterThan': operator.ge,
    'Equality': operator.eq}


exponential_function_names = {'exp', 'log'}
nonsmooth_function_names = {'sign', 'Abs', 'Min', 'Max', 'ceiling', 'floor'}
trigonometric_function_names = {'sin', 'cos', 'tan', 'cot', 'sec', 'csc'}
trigonometric_inverse_function_names = {'asin', 'acos', 'atan', 'acot', 'asec',
                                        'acsc'}
hyperbolic_function_names = {'sinh', 'cosh', 'tanh', 'coth', 'sech', 'csch'}
hyperbolic_inverse_function_names = {'asinh', 'acosh', 'atanh', 'acoth',
                                     'asech', 'acsch'}
# NOTE: All of these functions depend on the backend and therefore need to be
#       looked up dynamically to allow for backend switches!
comando_functions = set().union(
    exponential_function_names, nonsmooth_function_names,
    trigonometric_function_names, trigonometric_inverse_function_names,
    hyperbolic_function_names, hyperbolic_inverse_function_names
)


class SlotSerializationMixin(object):
    """A Mixin to make classes with slots serializable."""

    __slots__ = ()

    def __getstate__(self):
        """Get the state."""
        # Collect all data stored in state and slots of mro hierarchy
        state = {}
        for ty in reversed(type(self).__mro__):  # reversed to overwrite
            try:
                state.update(super(ty).__getstate__())
            except AttributeError:
                pass
            if hasattr(ty, '__slots__'):
                for slot in ty.__slots__:
                    if hasattr(self, slot):
                        state[slot] = getattr(self, slot)
        return state
        # try:
        #     super_state = super().__getstate__()
        # except AttributeError:
        #     super_state = {}
        # return {**super_state, **{slot: getattr(self, slot)
        #         for slot in self.__slots__ if hasattr(self, slot)}}

    def __setstate__(self, state):
        """Set the state."""
        for slot, value in state.items():
            setattr(self, slot, value)


class set_backend:
    """Set the symbolic backend used by COMANDO."""

    def __init__(self, backend_name):
        import comando
        self.previous_backend = comando.get_backend()
        self.backend_name = backend_name
        self._set(self.backend_name)

    def _set(self, backend_name):
        import comando
        if backend_name in comando.backends:
            comando.BACKEND = comando.backends[backend_name]
            self.update_symbols()
            comando.op_map = base_op_map.copy()
            for f_name in comando_functions:
                comando.op_map[f_name] = getattr(comando.get_backend(), f_name)

            def __getattr__(name):
                """Get attributes that can't be found from the backend."""
                from comando import core
                try:
                    return getattr(core, name)
                except AttributeError:
                    return getattr(get_backend(), name)
            comando.__getattr__ = __getattr__
        else:
            raise NotImplementedError(f'{backend_name} is not a known '
                                      'backend!\n Available backends are '
                                      f'{[*comando.backends]}!')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        import comando
        comando.BACKEND = self.previous_backend
        self._set(self.previous_backend.__name__)

    def update_symbols(self):
        """Update the Symbol definitions for the new backend."""
        import comando

        BE = comando.get_backend()

        comando.Zero = BE.S(0)  # TODO: Obsolete

        class Symbol(SlotSerializationMixin, BE.Symbol):
            """A placeholder for a value which can occur within expressions."""

            __slots__ = ('_value', '_indexed', '__qualname__', '_newargs')

            __qualname__ = 'Symbol'

            def __getnewargs__(self):
                return self._newargs

            def __getnewargs_ex__(self):
                return self._newargs, {}

            def __hash__(self):
                """Return the symbol's id, since it is unique."""
                return id(self)

            @property
            def indexed(self):
                """Check if the Symbol is indexed."""
                return self._indexed

            @property
            def value(self):
                """Get the Symbol's value."""
                return self._value

            @property
            def is_Parameter(self):
                return False

            @property
            def is_Variable(self):
                return False
        comando.Symbol = Symbol

        class Parameter(Symbol):
            """A `Symbol` representing a parameter whose value is known."""

            __slots__ = ('expansion', '_parent')

            __qualname__ = 'Parameter'

            def __new__(cls, name, value=NAN, parent=None):
                self = BE.Symbol.__new__(cls, name)
                # self = super().__new__(cls, name)

                self.value = value  # TODO: Add safeguards
                self._parent = parent

                self._newargs = name, value, parent

                return self

            @property
            def indexed(self):
                """Check whether the Parameter is indexed or not."""
                return self.expansion is not None

            @property
            def value(self):
                """Return the value or values of the Parameter."""
                return self._value if self.expansion is None \
                    else self.expansion.apply(lambda e: getattr(e, 'value'))

            @value.setter
            def value(self, data):
                """Set the value of the Parameter.

                In contrast to the Variable, a Parameter can be made indexed by
                simply providing some Mapping or a pandas.Series that imply
                both an index and values.
                """
                if isinstance(data, (Mapping, Series)):
                    self.expand(data)
                    self._value = None
                else:
                    self._value = None if data is None else float(data)
                    self.expansion = None

            @property
            def indices(self): return self.expansion.keys()

            @property
            def elements(self): return self.expansion.values

            @property
            def items(self): return self.expansion.items()

            @property
            def parent(self):
                """Return the parent of this parameter."""
                return self._parent

            def expand(self, data):
                """Expand the `Parameter` with indexed data."""
                if self._parent is not None:
                    raise RuntimeError(f'Attempted to expand parameter {self} '
                                       f'which is a member of {self._parent}!')
                self.expansion = Series((Parameter(f'{self.name}[{i}]', v,
                                                   self)
                                         for i, v in data.items()),
                                        data.keys(), dtype='O')

            def __getitem__(self, index):
                try:  # TODO: Custom errors still includes long stack trace...
                    return self.expansion[index]
                except KeyError:
                    if self.expansion is None:
                        raise TypeError(f'Parameter {self.name} is not '
                                        'indexed!')
                    raise IndexError("Parameter index out of range")

            def __setitem__(self, index, value):
                """Set the value of the element corresponding to the index."""
                self.expansion[index].value = value

            def __iter__(self):
                """Iterate over the elements of this Parameter."""
                if self.expansion is None:
                    raise TypeError(f'{self} is scalar.')
                for elem in self.elements:
                    yield elem

            @property
            def is_Parameter(self):
                return True
        comando.Parameter = Parameter

        class Variable(Symbol):
            """A `Symbol` representing a variable whose value is unknown."""

            __slots__ = ('_domain', '_bounds', '__bounds', '_init_val',
                         '_parent')

            __qualname__ = 'Variable'

            def __new__(cls, name, domain=REAL, bounds=(None, None),
                        init_val=None, indexed=False, parent=None):
                if domain is REAL:
                    assumptions = {'real': True}
                elif domain in {INTEGER, BINARY}:
                    assumptions = {'integer': True}
                else:
                    raise ValueError(f'Domain must be either {REAL}, '
                                     f'{INTEGER} or {BINARY}!')

                self = BE.Symbol.__new__(cls, name, **assumptions)
                # self = super().__new__(cls, name, **assumptions)
                # WITH DEFAULT KWARG: domain='REAL'
                # if domain == 'Binaries':
                #     self._bounds = comando.Interval(0, 1)
                #     self._domain = comando.INTEGER
                # domain = comando.S(domain)
                # if domain not in (comando.INTEGER, comando.REAL):
                #     raise ValueError("Kwarg 'domain' must be one of 'REAL',"
                #                      "' INTEGER' or 'Binaries'!")
                self._domain = domain
                self.bounds = bounds
                self.__bounds = None  # storing bounds when Variable is fixed
                # TODO: This is a second option to create time-indexed
                #       variables, but at some point we should decide wheter
                #       we simply distinguish by value or via a different class
                self._indexed = indexed
                self.init_val = init_val
                # NOTE: indexed Variables need an index that is specified later
                if not indexed:
                    self.value = self.init_val
                self._parent = parent

                self._newargs = name, domain, bounds, init_val, indexed, parent

                return self

            @property
            def parent(self):
                """Return the parent of this variable."""
                return self._parent

            @property
            def domain(self): return self._domain

            @property
            def bounds(self): return self._bounds

            @bounds.setter
            def bounds(self, bounds):
                if self.domain is BINARY:
                    lb = 0 if bounds[0] is None else bounds[0]
                    ub = 1 if bounds[1] is None else bounds[1]
                    if lb not in {0, 1} or ub not in {0, 1}:
                        raise ValueError('Bounds for binary variables may '
                                         'only be [0, 0], [0, 1] or [1, 1], '
                                         f'but are {[lb, ub]}!')
                else:
                    lb = -INF if bounds[0] is None else bounds[0]
                    ub = INF if bounds[1] is None else bounds[1]

                if lb <= ub:
                    self._bounds = lb, ub  # comando.Interval(lb, ub)
                else:
                    raise ValueError('Lower bound of', self.name,
                                     'must be less than or '
                                     'equal to upper bound!')

            @property
            def lb(self): return self.bounds[0]  # self.bounds.inf

            @lb.setter
            def lb(self, lb):
                lb = -INF if lb is None else lb
                if lb <= self.ub:
                    self.bounds = lb, self.ub  # comando.Interval(lb, self.ub)
                else:
                    raise ValueError('Lower bound of', self.name,
                                     'must be less than or '
                                     'equal to upper bound!')

            @property
            def ub(self): return self.bounds[1]  # self.bounds.sup

            @ub.setter
            def ub(self, ub):
                ub = INF if ub is None else ub
                if self.lb <= ub:
                    self.bounds = self.lb, ub  # comando.Interval(self.lb, ub)
                else:
                    raise ValueError('Upper bound of', self.name,
                                     'must be greater than or '
                                     'equal to lower bound!')

            @Symbol.value.setter
            def value(self, data):
                """Set the value of the Variable.

                A variable is declared to be either indexed or not at creation;
                according to this specification the value is a scalar or
                pandas.Series.
                If the Variable is indexed, the first time the value is set, it
                must be specified via a Mapping or pandas.Series, which imply
                an index.
                After this the values can be changed by simply providing
                iterables of appropriate length and the index and length can be
                adapted by again specifying the value via a Mapping or Series.
                """
                if self._indexed:
                    if isinstance(data, (Mapping, Series)):
                        self._value = Series(data)
                    if isinstance(data, (Iterable)):
                        try:
                            self._value = Series(data, self._value.index)
                        except AttributeError:
                            raise AttributeError("The variable's index has "
                                                 "not been specified yet, set "
                                                 "the value using a Mapping "
                                                 "or a pandas.Series!")
                        except ValueError:
                            raise ValueError("The variable's index does not "
                                             "match the length of the "
                                             "provided value!")
                            # data & index don't match -> infer from data
                            self._value = Series(data)
                else:
                    self._value = NAN if data is None else float(data)

            def fix(self, value=None):
                """Fix the variable by setting both bounds to `value`."""
                if value is None:
                    try:
                        value = self._value
                    except AttributeError:  # uninitialized indexed Variable
                        value = self._init_val
                if self.is_integer and not float(value).is_integer():
                    lb, ub = self.__bounds if self.__bounds else self.bounds
                    if lb <= value <= ub:
                        from warnings import warn
                        round_val = round(value)
                        warn(f'Fixing value of {self.name} to {round_val} '
                             f'instead of {value}!')
                        value = round_val
                    else:
                        kind = self.domain.name.lower()
                        raise ValueError(f'Cannot fix {kind} variable '
                                         f'"{self.name}" with bounds '
                                         f'{lb, ub} to non-integer value '
                                         f'{value}!')
                # Get original bounds (either the previously fixed or the
                # current ones)
                lb, ub = self.__bounds if self.__bounds else self.bounds
                if value <= lb - EPS or ub + EPS <= value:
                    raise ValueError(f'Value {value} is not within original '
                                     f'bounds [{lb}, {ub}]')
                self.__bounds = lb, ub  # NoOp or setting current bounds
                self.value = value
                self.bounds = (value, value)

            def unfix(self):
                """Recover the original bounds."""
                if self.__bounds:
                    self.bounds = self.__bounds
                    self.__bounds = None

            @property
            def init_val(self):
                return self._init_val

            @init_val.setter
            def init_val(self, val):
                from math import isnan, isinf, ceil, floor
                lb, ub = self._bounds
                if val is None:
                    if self.is_binary:
                        self._init_val = 0
                        return
                    tmp = (ub + lb) * 0.5
                    if isnan(tmp) or isinf(tmp):
                        if not isnan(ub) and not isinf(ub):
                            self._init_val = \
                                floor(ub) if self.is_integer else ub
                        elif not isnan(lb) and not isinf(lb):
                            self._init_val = \
                                ceil(lb) if self.is_integer else lb
                        else:
                            self._init_val = 0
                        return
                    val = tmp
                if self.is_integer:  # make adjustments for the integer case
                    # ... round val to nearest int
                    if not float(val).is_integer():
                        l_val = floor(val)
                        u_val = ceil(val)
                        val = l_val if val - l_val <= u_val - val else u_val
                    # TODO: could be moved to the bounds property?
                    lb = ceil(lb)
                    ub = floor(ub)
                if lb <= val <= ub:
                    self._init_val = val
                    return
                raise ValueError(f'init_val {val} is incompatible with bounds '
                                 f'[{lb}, {ub}]!')

            @property
            def is_integer(self):
                return self.domain in {BINARY, INTEGER}

            @property
            def is_binary(self):
                return self.domain is BINARY

            # TODO: Might need EPS for comparisons as well
            @property
            def is_positive(self):
                """Check if all possible values of the variable are positive.

                We can assert positivity if the lower bound is positive,
                otherwise we can assert nonpositivity if the upper bound is
                nonpositive. If we cannot assert either of these facts, the
                variable may contain both positive and negative values. To
                reflect this we return None.
                """
                return True if self.lb > 0 else False if self.ub <= 0 else None

            @property
            def is_negative(self):
                """Check if all possible values of the variable are negative.

                We can assert negativity if the upper bound is negative,
                otherwise we can assert nonnegativity if the lower bound is
                nonnegative. If we cannot assert either of these facts, the
                variable may contain both positive and negative values. To
                reflect this we return None.
                """
                return True if self.ub < 0 else False if self.lb >= 0 else None

            @property
            def is_nonnegative(self):
                """Check if all possible values of the variable are negative.

                This is the fuzzy not of self.is_negative
                """
                # True if self.lb >= 0 else False if self.ub < 0 else None
                return comando.utility.fuzzy_not(self.is_negative)

            @property
            def is_nonpositive(self):
                """Check if all possible values of the variable are negative.

                This is the fuzzy not of self.is_negative
                """
                # True if self.ub <= 0 else False if self.lb > 0 else None
                return comando.utility.fuzzy_not(self.is_positive)

            @property
            def is_Variable(self):
                return True
        comando.Variable = Variable

        class VariableVector(Symbol):
            """A `Symbol` representing a vector of `Variables`."""

            __slots__ = ('_domain', '_bounds', '__bounds', '_init_val',
                         'expansion')

            __qualname__ = 'VariableVector'

            def __new__(cls, name, domain=REAL, bounds=(None, None),
                        init_val=None):
                if domain is REAL:
                    assumptions = {'real': True}
                elif domain in {INTEGER, BINARY}:
                    assumptions = {'integer': True}
                else:
                    raise ValueError(f'Domain must be either {REAL}, '
                                     f'{INTEGER} or {BINARY}!')

                self = BE.Symbol.__new__(cls, name, **assumptions)
                # self = sympy.Symbol.__xnew__(cls, name, **assumptions)
                self._domain = domain
                self.expansion = Series(dtype='O')
                # self.bounds = bounds
                self.bounds = bounds
                # TODO
                self.__bounds = None  # for storing bounds
                # self._bounds = bounds
                self.init_val = init_val

                self._newargs = name, domain, bounds, init_val

                return self

            # TODO: Rename (here and elsewehere) to is_indexed
            @property
            def indexed(self):
                return True

            @property
            def is_expanded(self):
                return bool(len(self.expansion))

            @property
            def domain(self): return self._domain

            @property
            def indices(self): return self.expansion.keys()

            @property
            def elements(self): return self.expansion.values

            @property
            def items(self): return self.expansion.items()

            def _get_property(self, property):
                """Get a `Series` with property values from `self.elements`."""
                return (self.expansion.apply(lambda e: getattr(e, property))
                        .rename(property))

            def _set_property(self, property, scalar_or_mapping):
                """Set the property values of `self.elements`."""
                # Attempting to treat `scalar_or_mapping` as a mapping...
                if not self.is_expanded:
                    raise RuntimeError(f'VariableVector {self.name} has not '
                                       'been instantiated yet.')
                if isinstance(scalar_or_mapping, (Mapping, Series, DataFrame)):
                    for i, elem in self.expansion.items():
                        try:
                            setattr(elem, property, scalar_or_mapping[i])
                        except KeyError:
                            continue  # Leave existing value
                    return
                try:  # treating as iterable
                    for (i, elem), val in zip(self.expansion.items(),
                                              scalar_or_mapping):
                        setattr(elem, property, val)
                    return
                except TypeError:  # scalar_or_mapping not iterable -> scalar
                    for elem in self.elements:
                        setattr(elem, property, scalar_or_mapping)
                    return
                raise RuntimeError(f'Could not set {property} property of '
                                   f'{self.name} with {scalar_or_mapping}!')

            @property
            def value(self):
                values = self._get_property('value')
                return values if len(values) else None

            @value.setter
            def value(self, values):
                self._set_property('value', values)

            @property
            # def bounds(self): return self._bounds
            def bounds(self):
                return (self._get_property('lb'), self._get_property('ub')) \
                    if self.is_expanded else self._bounds

            @bounds.setter
            def bounds(self, bounds):
                if self.is_expanded:
                    # NOTE: We expect bounds to be a 2-tuple of scalars,
                    #       Iterables or Mappings
                    _lb, _ub = (0, 1) if self.domain is BINARY else (-INF, INF)
                    lb, ub = self._get_property('lb'), self._get_property('ub')
                    if isinstance(bounds[0], (Iterable, Mapping, Series)):
                        lb.update(Series(bounds[0]).fillna(_lb))
                    else:
                        lb[:] = _lb if bounds[0] is None else bounds[0]
                    if isinstance(bounds[1], (Iterable, Mapping, Series)):
                        ub.update(Series(bounds[1]).fillna(_ub))
                    else:
                        ub[:] = _ub if bounds[1] is None else bounds[1]
                    if all(lb <= ub):
                        self._bounds = lb.min(), ub.max()
                        self._set_property('bounds', lb.to_frame().join(ub).T)
                    else:
                        raise ValueError('Lower bound of ' + self.name
                                         + ' must be less than or '
                                         'equal to upper bound!')
                else:
                    # NOTE: We expect bounds to be a 2-tuple of scalars
                    if self.domain is BINARY:
                        lb = 0 if bounds[0] is None else bounds[0]
                        ub = 1 if bounds[1] is None else bounds[1]
                        if lb not in {0, 1} or ub not in {0, 1}:
                            raise ValueError('Bounds for binary variables may '
                                             'only be [0, 0], [0, 1] or '
                                             f'[1, 1], but are {[lb, ub]}!')
                    else:
                        lb = -INF if bounds[0] is None else bounds[0]
                        ub = INF if bounds[1] is None else bounds[1]

                    if lb <= ub:
                        self._bounds = lb, ub  # comando.Interval(lb, ub)
                    else:
                        raise ValueError('Lower bound of ' + self.name
                                         + ' must be less than or '
                                         'equal to upper bound!')

            @property
            # def lb(self): return self._bounds[0]
            def lb(self):
                return self._get_property('lb') if self.is_expanded \
                    else self._bounds[0]

            @lb.setter
            def lb(self, lb):
                try:  # Scalar lb
                    lb = -INF if lb is None else float(lb)
                except (ValueError, TypeError):  # assume intended for elements
                    # setting lb by iterable
                    old_lb = self.lb
                    try:
                        self._set_property('lb', lb)
                    except ValueError as e:
                        self._set_property('lb', old_lb)
                        raise e
                    self._bounds = min(self.lb), self._bounds[1]
                    return
                ub = self.ub
                if self.is_expanded:
                    if all(lb <= ub):
                        self._set_property('lb', lb)
                        self._bounds = (self.lb.min(), self._bounds[1])
                    else:
                        raise ValueError('Lower bound of ' + self.name
                                         + ' must be less than or '
                                         'equal to upper bound!')
                else:
                    if lb <= ub:
                        self._bounds = lb, ub
                    else:
                        raise ValueError('Lower bound of ' + self.name
                                         + ' must be less than or '
                                         'equal to upper bound!')

            @property
            # def ub(self): return self._bounds[1]
            def ub(self):
                return self._get_property('ub') if self.is_expanded \
                    else self._bounds[1]

            @ub.setter
            def ub(self, ub):
                try:  # Scalar ub
                    ub = INF if ub is None else float(ub)
                except (ValueError, TypeError):  # assume intended for elements
                    # setting ub by iterable
                    old_ub = self.ub
                    try:
                        self._set_property('ub', ub)
                    except ValueError as e:
                        self._set_property('ub', old_ub)
                        raise e
                    self._bounds = self._bounds[0], max(self.ub)
                    return
                lb = self.lb
                if self.is_expanded:
                    if all(lb <= ub):
                        self._set_property('ub', ub)
                        self._bounds = self._bounds[0], self.ub.max()
                    else:
                        raise ValueError('Upper bound of ' + self.name
                                         + ' must be greater than or '
                                         'equal to lower bound!')
                else:
                    if lb <= ub:
                        self._bounds = lb, ub
                    else:
                        raise ValueError('Upper bound of ' + self.name
                                         + ' must be greater than or '
                                         'equal to lower bound!')
                #
                # self._set_property('ub', ub)
                # ub = INF if ub is None else ub
                # if self.lb <= ub:
                #     self.bounds = self.lb, ub

            @property
            def init_val(self):
                return self._init_val

            @init_val.setter
            def init_val(self, val):
                from math import isnan, isinf, ceil, floor
                lb, ub = self._bounds
                if val is None:
                    tmp = (ub + lb) * 0.5
                    if isnan(tmp) or isinf(tmp):
                        if not isnan(ub) and not isinf(ub):
                            self._init_val = \
                                floor(ub) if self.is_integer else ub
                        elif not isnan(lb) and not isinf(lb):
                            self._init_val = \
                                ceil(lb) if self.is_integer else lb
                        else:
                            self._init_val = 0
                        return
                    val = tmp
                if self.is_integer:  # make adjustments for the integer case
                    if not float(val).is_integer():  # ... round to nearest int
                        l_val = floor(val)
                        u_val = ceil(val)
                        val = l_val if val - l_val <= u_val - val else u_val
                    # TODO: could be moved to the bounds property?
                    lb = ceil(lb)
                    ub = floor(ub)
                if lb <= val <= ub:
                    self._init_val = val
                    return
                raise ValueError(f'init_val {val} is incompatible with bounds '
                                 f'[{lb}, {ub}]!')

            def fix(self, value=None):
                """Fix `self.elements` by setting both bounds to `value`."""
                if value is None:
                    try:
                        value = self.value
                    except AttributeError:  # uninitialized
                        value = self._init_val

                lb, ub = self.__bounds if self.__bounds else self.bounds
                if self.is_expanded:
                    if any(value <= lb - EPS) or any(ub + EPS <= value):
                        raise ValueError(f'Value {value} is not within '
                                         f'original bounds [{lb}, {ub}]')
                    self.value = value
                elif value <= lb - EPS or ub + EPS <= value:
                    raise ValueError(f'Value {value} is not within '
                                     f'original bounds [{lb}, {ub}]')
                self.__bounds = lb, ub  # NoOp or setting current bounds
                self.bounds = (value, value)

            def unfix(self):
                """Recover the original bounds."""
                self.bounds = self.__bounds
                self.__bounds = None

            @property
            def is_integer(self):
                return self.domain in {BINARY, INTEGER}

            @property
            def is_binary(self):
                return self.domain is BINARY

            @property
            def is_positive(self):
                """Check if all possible values of the variable are positive.

                We can assert positivity if the lower bound is positive,
                otherwise we can assert nonpositivity if the upper bound is
                nonpositive. If we cannot assert either of these facts, the
                variable may contain both positive and negative values. To
                reflect this we return None.
                """
                try:
                    return True if self.lb > 0 else \
                        False if self.ub <= 0 else None
                except ValueError:
                    return True if all(self.lb > 0) else \
                        False if all(self.ub <= 0) else None

            @property
            def is_negative(self):
                """Check if all possible values of the variable are negative.

                We can assert negativity if the upper bound is negative,
                otherwise we can assert nonnegativity if the lower bound is
                nonnegative. If we cannot assert either of these facts, the
                variable may contain both positive and negative values. To
                reflect this we return None.
                """
                try:
                    return True if self.ub < 0 else \
                        False if self.lb >= 0 else None
                except ValueError:
                    return True if all(self.ub < 0) else \
                        False if all(self.lb >= 0) else None

            @property
            def is_nonnegative(self):
                """Check if all possible values of the variable are negative.

                This is the fuzzy not of self.is_negative
                """
                # True if self.lb >= 0 else False if self.ub < 0 else None
                return comando.utility.fuzzy_not(self.is_negative)

            @property
            def is_nonpositive(self):
                """Check if all possible values of the variable are negative.

                This is the fuzzy not of self.is_negative
                """
                # True if self.ub <= 0 else False if self.lb > 0 else None
                return comando.utility.fuzzy_not(self.is_positive)

            def instantiate(self, index):
                """Create a `Variable` instance for every element in `index`."""
                self.expansion = Series((Variable(f'{self.name}[{i}]',
                                                  self.domain, self._bounds,
                                                  self.init_val, parent=self)
                                         for i in index), index, dtype='O')

            def __getitem__(self, index):
                try:  # TODO: Custom errors still includes long stack trace...
                    return self.expansion[index]
                except KeyError:
                    raise IndexError("Vector index out of range")

            def __setitem__(self, index, value):
                """Set the value of the element corresponding to the index."""
                self.expansion[index].value = value

            def __iter__(self):
                """Iterate over the elements of this VariableVector."""
                for elem in self.elements:
                    yield elem
        comando.VariableVector = VariableVector

        try:
            Singleton = BE.core.singleton.Singleton
        except AttributeError:
            # NOTE: raises TypeError: metaclass conflict: the metaclass of a
            #                         derived class must be a (non-strict)
            #                         subclass of the metaclasses of all its
            #                         bases
            # Singleton = BE.lib.symengine_wrapper.Singleton
            Singleton = object

        class Cyclic(BE.Basic, metaclass=Singleton):
            """A Placeholder class for cyclic initial conditions."""

            __qualname__ = 'Cyclic'
        comando.Cyclic = Cyclic
        comando.cyclic = Cyclic()


try:
    import symengine
    backends['symengine'] = symengine
    from symengine.lib import symengine_wrapper
    symengine.Boolean = symengine_wrapper.Boolean

    def _0(self):
        return self.args[0]

    def _1(self):
        return self.args[1]

    symengine_wrapper.Relational.lhs = property(_0)
    symengine_wrapper.Relational.rhs = property(_1)
    symengine_wrapper.Le.lts = property(_0)
    symengine_wrapper.Le.gts = property(_1)
except ModuleNotFoundError:
    pass
finally:
    def get_backend():
        """Get the symbolic backend used by COMANDO."""
        return BACKEND

    set_backend('sympy')
    import comando
    comando.get_backend = get_backend
    comando.backend = set_backend


from .core import *
from comando import utility

# Teaching all sympy expressions about their value (Won't work with symengine!)
sympy.Expr.value = property(lambda self: comando.utility.evaluate(self))

# NOTE: must be last line of file!
__version__ = '1.0.2'
