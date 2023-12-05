"""Tests for various utility functions."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import pytest
import comando
from comando import Variable
from comando.utility import sum_bounds, mul_bounds, prod_bounds, pow_bounds, \
    _normalize_bounds, bounds, str_parse, StrParser
inf = float('inf')


class TestBounding():

    @pytest.mark.parametrize('args, result',
                             [([(1, 2), (3, 4), (0, 0)], (4, 6)),
                              ([(1, 2), (-3, 4), (0, 0)], (-2, 6)),
                              ([(1, 2), (3, -4), (0, 0)], (4, -2)),
                              ([(1, 2), (-3, -4), (0, 0)], (-2, -2)),
                              ([(-1, 2), (3, 4), (0, 0)], (2, 6)),
                              ([(-1, 2), (-3, 4), (0, 0)], (-4, 6)),
                              ([(-1, 2), (3, -4), (0, 0)], (2, -2)),
                              ([(-1, 2), (-3, -4), (0, 0)], (-4, -2)),
                              ([(1, -2), (3, 4), (0, 0)], (4, 2)),
                              ([(1, -2), (-3, 4), (0, 0)], (-2, 2)),
                              ([(1, -2), (3, -4), (0, 0)], (4, -6)),
                              ([(1, -2), (-3, -4), (0, 0)], (-2, -6)),
                              ([(-1, -2), (3, 4), (0, 0)], (2, 2)),
                              ([(-1, -2), (-3, 4), (0, 0)], (-4, 2)),
                              ([(-1, -2), (3, -4), (0, 0)], (2, -6)),
                              ([(-1, -2), (-3, -4), (0, 0)], (-4, -6)),
                              ([(-inf, 1), (1, 2)], (-inf, 3)),
                              ([(-2, 0), (-inf, 2)], (-inf, 2)),
                              ([(-2, 1), (3, inf)], (1, inf)),
                              ([(3, 1), (-inf, inf)], (-inf,  inf))])
    def test_sum_bounds(self, args, result):
        assert sum_bounds(*args) == result

    @pytest.mark.parametrize('a, b, result',
                             [((1, 2), (3, 4), (3, 8)),
                              ((1, 2), (-3, 4), (-6, 8)),
                              ((1, 2), (3, -4), (-8, 6)),
                              ((1, 2), (-3, -4), (-8, -3)),
                              ((-1, 2), (3, 4), (-4, 8)),
                              ((-1, 2), (-3, 4), (-6, 8)),
                              ((-1, 2), (3, -4), (-8, 6)),
                              ((-1, 2), (-3, -4), (-8, 4)),
                              ((1, -2), (3, 4), (-8, 4)),
                              ((1, -2), (-3, 4), (-8, 6)),
                              ((1, -2), (3, -4), (-6, 8)),
                              ((1, -2), (-3, -4), (-4, 8)),
                              ((-1, -2), (3, 4), (-8, -3)),
                              ((-1, -2), (-3, 4), (-8, 6)),
                              ((-1, -2), (3, -4), (-6, 8)),
                              ((-1, -2), (-3, -4), (3, 8))])
    def test_mul_bounds(self, a, b, result):
        assert mul_bounds(a, b) == result

    @pytest.mark.parametrize('args, result',
                             [([(1, 2), (3, 4), (1, 1)], (3, 8)),
                              ([(1, 2), (-3, 4), (1, 1)], (-6, 8)),
                              ([(1, 2), (3, -4), (1, 1)], (-8, 6)),
                              ([(1, 2), (-3, -4), (1, 1)], (-8, -3)),
                              ([(-1, 2), (3, 4), (1, 1)], (-4, 8)),
                              ([(-1, 2), (-3, 4), (1, 1)], (-6, 8)),
                              ([(-1, 2), (3, -4), (1, 1)], (-8, 6)),
                              ([(-1, 2), (-3, -4), (1, 1)], (-8, 4)),
                              ([(1, -2), (3, 4), (1, 1)], (-8, 4)),
                              ([(1, -2), (-3, 4), (1, 1)], (-8, 6)),
                              ([(1, -2), (3, -4), (1, 1)], (-6, 8)),
                              ([(1, -2), (-3, -4), (1, 1)], (-4, 8)),
                              ([(-1, -2), (3, 4), (1, 1)], (-8, -3)),
                              ([(-1, -2), (-3, 4), (1, 1)], (-8, 6)),
                              ([(-1, -2), (3, -4), (1, 1)], (-6, 8)),
                              ([(-1, -2), (-3, -4), (1, 1)], (3, 8))])
    def test_prod_bounds(self, args, result):
        assert prod_bounds(*args) == result

    def test_pow_bounds_singleton_exponent(self):
        with pytest.raises(ValueError) as e:
            pow_bounds((-1, 1), (1, 2))
            exception_msg = e.value.args[0]
            assert exception_msg == "Interval power requires real-valued " \
                "singleton exponents (with -∞ < lb == ub < ∞)!"

    @pytest.mark.parametrize('base',
                             [(-inf, -1),
                              (-inf, 0),
                              (-inf, 1),
                              (-inf, inf),
                              (-1, inf),
                              (0, inf),
                              (1, inf),
                              (-1, 1),
                              (0, 0)])
    def test_pow_bounds_with_0_exponent(self, base):
        assert pow_bounds(base, (0, 0)) == (1, 1)

    @pytest.mark.parametrize('base, expo',
                             [((-inf, 0), -1),
                              ((-inf, 0), -1.1),
                              ((-inf, 1), -1),
                              ((-inf, 1), -1.1),
                              ((-inf, inf), -1),
                              ((-inf, inf), -1.1),
                              ((-1, inf), -1),
                              ((-1, inf), -1.1),
                              ((0, inf), -1),
                              ((0, inf), -1.1),
                              ((-1, 1), -1),
                              ((-1, 1), -1.1),
                              ((0, 0), -1),
                              ((0, 0), -1.1)])
    def test_pow_bounds_with_0_in_base_and_negative_exponent(self, base, expo):
        with pytest.raises(ValueError) as e:
            pow_bounds(base, (expo, expo))
            exception_msg = e.value.args[0]
            assert exception_msg == "Base of 0 in power with negative exponent"

    @pytest.mark.parametrize('base, expo, res',
                             [((-1, 2), 4, (0, 2 ** 4)),
                              ((-1, 2), 3, ((-1) ** 3, 2 ** 3)),
                              ((1, 2), 4, (1 ** 4, 2 ** 4)),
                              ((1, 2), -4, (2 ** -4, 1 ** -4)),
                              ((1, 2), 5, (1 ** 5, 2 ** 5)),
                              ((1, 2), -5, (2 ** -5, 1 ** -5)),
                              ((-2, -1), 4, ((-1) ** 4, (-2) ** 4)),
                              ((-2, -1), -4, ((-2) ** -4, (-1) ** -4)),
                              ((-2, -1), 5, ((-2) ** 5, (-1) ** 5)),
                              ((-2, -1), -5, ((-1) ** -5, (-2) ** -5))])
    def test_pow_bounds_with_integer_exponent(self, base, expo, res):
        assert pow_bounds(base, (expo, expo)) == res

    @pytest.mark.parametrize('base, expo',
                             [((-1, 1), 1.1),
                              ((-1, -0.1), -1.1)])
    def test_pow_bounds_with_negative_base_and_real_exponent(self, base, expo):
        with pytest.raises(ValueError) as e:
            pow_bounds(base, (expo, expo))
            exception_msg = e.value.args[0]
            assert exception_msg == "Negative base with a non-integer exponent"

    @pytest.mark.parametrize('base, expo, res',
                             [((0, 1), 1.1, (0 ** 1.1, 1 ** 1.1)),
                              ((0.5, 2), 1.1, (0.5 ** 1.1, 2 ** 1.1)),
                              ((0.5, 2), -1.1, (2 ** -1.1, 0.5 ** -1.1))])
    def test_pow_bounds_with_real_exponent(self, base, expo, res):
        assert pow_bounds(base, (expo, expo)) == res


def test_bounds():
    # NOTE: For some reason just using 'from comando import Parameter' causes a
    #       problem in testing. While manually executing it seems to work fine.
    #       Therefore we use this approach instead...
    Parameter = comando.Parameter
    assert Parameter == comando.Parameter

    # Parameters
    pr = Parameter('pr', value=3.2)  # positive real number
    nr = Parameter('nr', value=-2.3)  # negative real number
    poi = Parameter('poi', value=3)  # positive odd integer
    noi = Parameter('noi', value=-5)  # negative odd integer
    pei = Parameter('pei', value=2)  # positive even integer
    nei = Parameter('nei', value=-4)  # negative even integer
    zero = Parameter('zero', value=0)  # zero
    undef = Parameter('undef')  # None value

    positive_parameters = {pr, poi, pei}
    negative_parameters = {nr, noi, nei}
    parameters = positive_parameters.union(negative_parameters)
    special_parameters = {zero, undef}
    for p in parameters.union(special_parameters):
        # import comando
        # print(type(p), comando.Parameter)
        assert bounds(p) == _normalize_bounds((p.value, p.value))

    for p in parameters:
        for i in [-inf, -2, 1.5, -1, -0.5, -0.0, 0, 0.0, 0.5, 1, 1.5, 2, inf]:
            assert bounds(p + i) == (p.value + i, p.value + i)
            assert bounds(p * i) == (p.value * i, p.value * i)

    for a in [-inf, -2, 1.5, -1, -0.5, -0.0, 0, 0.0, 0.5, 1, 1.5, 2, inf]:
        for b in [-2, 1.5, -1, -0.5, -0.0, 0, 0.0, 0.5, 1, 1.5, 2]:
            if (a < 0 and not float(b).is_integer()) or (a == 0 and b < 0):
                continue  # Error cases handled by 'TestBounding' unit tests
            res = a ** b
            dummy = Parameter(str(a), value=a)
            assert bounds(dummy ** b) == (res, res)

    assert bounds(undef + 1) == (-inf, inf)  # undefined value -> no bounds!

    # TODO: Variables
    u = Variable('u', bounds=(None, None))  # unbounded...
    unn = Variable('unn', bounds=(0, None))  # & nonnegative
    unp = Variable('unp', bounds=(None, 0))  # & nonpositive
    un = Variable('un', bounds=(1, None))  # & negative
    up = Variable('up', bounds=(None, -1))  # & positive

    b = Variable('b', bounds=(-10, 10))  # bounded...
    bnn = Variable('bnn', bounds=(0, 10))  # & nonnegative
    bnp = Variable('bnp', bounds=(-10, 0))  # & nonpositive
    bn = Variable('bn', bounds=(1, 10))  # & positive
    bp = Variable('bp', bounds=(-10, -1))  # & negative
    unbounded_vars = {u, unn, unp, un, up}
    bounded_vars = {b, bnn, bnp, bn, bp}
    vars = unbounded_vars.union(bounded_vars)
    for v in vars:
        assert bounds(v) == _normalize_bounds(v.bounds)
    for v in bounded_vars:
        lb, ub = v.bounds
        for i in [-inf, -2, 1.5, -1, -0.5, -0.0, 0, 0.0, 0.5, 1, 1.5, 2, inf]:
            assert bounds(v + i) == (lb + i, ub + i)
            v.value


        for i in [-2, 1.5, -1, -0.5, -0.0, 0, 0.0, 0.5, 1, 1.5, 2]:
            assert bounds(v * i) == tuple(sorted([lb * i, ub * i]))
        # TODO: Test case for 0 * ±inf

    # TODO: Multiplication for unbounded vars

    # TODO: Pow for vars

    # TODO: Some more general expressions


@pytest.mark.xfail
def test_parse():
    assert False  # TODO


def test_str_parse():
    a = Variable('a')
    b = Variable('b')
    c = Variable('c')
    sym_map = {v: str(v) for v in [a, b, c]}

    for sp in StrParser(sym_map), lambda expr: str_parse(expr, sym_map):
        assert sp(a + 1) == '1 + a'  # constants come first
        assert sp(-(a + 1)) == '-(1 + a)'
        assert sp(a + b) == 'a + b'
        assert sp(1 - a + b - c) == '1 + b - (a + c)'
        assert sp(-(a + b)) == '-(a + b)'
        assert sp(-a - b) == '-(a + b)'
        assert sp(a * b) == 'a * b'
        assert sp(-a * b) == '-a * b'
        assert sp(a * -b) == '-a * b'
        assert sp(a * 2) == '2 * a'
        assert sp(a * b + c) == 'a * b + c'
        # obsolete parentheses are skipped
        assert sp((a * b) + c) == 'a * b + c'
        assert sp(a * (b + c)) == '(b + c) * a'
        assert sp(a + -(b + c)) == 'a - (b + c)'
        assert sp(a + -(b * c)) == 'a - b * c'
        assert sp(a / (b + c)) == 'a / (b + c)'
        assert sp(-a / (b + c)) == '-a / (b + c)'
        assert sp(a / -(b + c)) == '-a / (b + c)'
        assert sp(-(a / (b + c))) == '-a / (b + c)'
        assert sp(a / b / c) == 'a / (b * c)'
        assert sp(a / (b * c)) == 'a / (b * c)'
        assert sp(-a / (b * c)) == '-a / (b * c)'
        assert sp(a / -(b * c)) == '-a / (b * c)'
        assert sp(-(a / (b * c))) == '-a / (b * c)'
        assert sp(1/a * 1/b * 1/c) == '1 / (a * b * c)'
        assert sp((a + b) / c) == '(a + b) / c'
        assert sp(a * b / c) == 'a * b / c'
        assert sp(-(a + b) / c) == '-(a + b) / c'
        assert sp(-(a + b) / (b + c)) == '-(a + b) / (b + c)'
        assert sp((a + b) / -(b + c)) == '-(a + b) / (b + c)'
        assert sp(-((a + b) / (b + c))) == '-(a + b) / (b + c)'
        assert sp(-(a * b) / (2 * c)) == '-a * b / (2 * c)'
        assert sp((a * b) / -(2 * c)) == '-a * b / (2 * c)'
        assert sp(-((a * b) / (2 * c))) == '-a * b / (2 * c)'
        assert sp(1/a) == '1 / a'
        assert sp(1/(1 + a)) == '1 / (1 + a)'
        assert sp(1/(a * b)) == '1 / (a * b)'
        assert sp(1/a**2) == '1 / a ** 2'
        assert sp(a**1/2) == 'a / 2'
        # Exponents that are python floats produce long string representations
        assert sp(a**0.5) == 'a ** 0.500000000000000'
        assert sp(a**(1/2)) == 'a ** 0.500000000000000'
        # exact sympy rationals are wrapped in parentheses
        assert sp(a**(comando.S(1)/2)) == 'a ** (1/2)'
        assert sp(comando.sqrt(a)) == 'a ** (1/2)'
        assert sp(-a ** b) == '-a ** b'
        assert sp(-(a ** b)) == '-a ** b'
        assert sp((-a) ** b) == '(-a) ** b'
        assert sp(a ** -b) == '1 / a ** b'
        assert sp(a ** (-b)) == '1 / a ** b'
        assert sp(1 / a ** b) == '1 / a ** b'

    try:
        import symengine
        backends = 'sympy', 'symengine'
    except ModuleNotFoundError:
        backends = 'sympy',

    for backend in backends:
        comando.set_backend(backend)
        x = comando.Variable('x')
        assert str_parse(-(1 + x), {x: x.name}) == '-(1 + x)'

    # DEBUG:
    # from operator import add, sub, mul, truediv, pow
    # args =  {comando.S(42), -comando.S(42), a, -a, a + b, a - b, a * b, a / b, a ** b}
    # for arg1 in args:
    #     for arg2 in (args - {arg1}):
    #         arg2 = arg2.subs(a, c)
    #         for op in add, sub, mul, truediv, pow:
    #             # print(arg1, op.__name__, arg2, '==')
    #             expr = op(arg1, arg2)
    #             # print(expr)
    #             # print()
    #             # print()
    #             z = f"({expr}) - ({str_parse(expr, sym_map).replace('^', '**')})"
    #             if comando.S(z) != 0:
    #                 print(expr)
    #
    # expr = -c*(a - b)
    # reint = comando.S(str_parse(expr, sym_map).replace('^', '**'))
    # expr
    # reint
    # DEBUG
