"""Components for the ORC case study.

Based on:

@InCollection{huster2019impact,
  author    = {Wolfgang R. Huster and Artur M. Schweidtmann and Alexander
               Mitsos},
  booktitle = {Computer Aided Chemical Engineering},
  title     = {Impact of accurate working fluid properties on the globally
               optimal design of an organic rankine cycle},
  doi       = {10.1016/b978-0-12-818597-1.50068-0},
  pages     = {427--432},
  publisher = {Elsevier},
  volume    = {47},
  year      = {2019},
}
"""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import numpy

import comando
from comando import Min, Max, tanh, log

from comando.utility import define_function


def add_scaled_var(component, name, lb, ub, init=None, design=False):
    """Add a new variable scaled to the range [0, 1] to the model.

    Returns
    -------
    unscaled_var : Expression
        the expression corresponding to the unscaled variable
    scaled_var : Variable
        the scaled variable with bounds [0, 1]
    """
    if init is not None:
        init = (init - lb) / (ub - lb)
    Var = component.make_design_variable if design \
        else component.make_operational_variable
    scaled_var = Var(name + '_scaled', bounds=(0, 1), init_val=init)
    unscaled_var = lb + scaled_var * (ub - lb)
    return component.add_expression(name, unscaled_var), scaled_var


def make_linar_func(ref, **dependencies):
    """Define a linear function.

    For each dependency a reference value and sensitivity can be given.
    The resulting function then has the form:
        ref + sum(dep_sens * (dep_val - dep_ref))

    Arguments
    ---------
    ref : float
        the reference value
    dependencies : dict[str] -> (float, float)
        mapping from name of the dependency to a 2-tuple containing the
        associated reference value and the sensitivity
    """

    def func(**kwargs):
        return ref + sum((dep_val - dependencies[dep][0])
                         * dependencies[dep][1]
                         for dep, dep_val in kwargs.items())

    return func


def A_ST(d_s):
    """Calculate the heat transfer area for a shell & tube exchanger.

    The assumption is that tube parameters are fixed and length, baffle
    spacing, etc., are only dependent on the shell diameter.
    """
    return 25396.6478670648 \
        + 1174.5329443002 * tanh(0.419223942469735
                                 - 0.911727986258998 * d_s) \
        - 26510.3978290862 * tanh(2.19067185114821
                                  - 0.645904206987058 * d_s)


CEPCI_HIST = {
    2020: 596.2,
    2019: 607.5,
    2018: 603.1,
    2017: 567.5,
    2016: 541.7,
    2015: 556.8,
    2014: 576.1,
    2013: 567.3,
    2012: 584.6,
    2011: 585.7,
    2010: 550.8,
    2009: 521.9,
    2008: 575.4,
    2007: 525.4,
    2006: 499.6,
    2005: 468.2,
    2004: 444.2,
    2003: 402.0,
    2002: 395.6,
    2001: 394.3,
    2000: 394.1,
}


def CEPCI_extrapolation(year):
    """Return an extrapolation of the CEPCI for the given year.

    The result is based on a linar regression from the 15 most recent years in
    CEPCI_HIST.
    """
    import numpy as np

    x = np.array(list(CEPCI_HIST.keys()))
    y = np.array(list(CEPCI_HIST.values()))

    N = min(len(x), 15)
    # indices of elements between the Nth largest and the largest element
    i_N = np.argpartition(x, -N)[-N:]
    # sort the indices according to the x value
    i_N = i_N[np.argsort(x[i_N])]
    x_N = x[i_N]
    y_N = y[i_N]

    # coefficients of the linear regression
    # NOTE: There is no need to solve an LP, we can just use KKT conditions of
    #       the minimization problem for mean squared error as second
    #       derivatives are nonnegative constants
    m = ((N * np.sum(x_N * y_N) - np.sum(x_N) * np.sum(y_N))
         / (N * np.sum(x_N * x_N) - np.sum(x_N) ** 2))
    b = np.sum(y_N - m * x_N) / N

    return m * year + b


def CEPCI_correction(base_CEPCI, CEPCI=None, year=None):
    """Correct purchase costs from a base year to a given year."""
    if CEPCI is None:
        if year is None:
            raise ValueError('Need either CEPCI value or year for '
                             'extrapolation!')
        CEPCI = CEPCI_extrapolation(year)
    return CEPCI / base_CEPCI


USD_EUR_2014 = 0.7541  # amount of US $ we got per € in 2014
COST_CONVERSION_ASTOLFI = USD_EUR_2014 \
    * CEPCI_correction(CEPCI_HIST[2014], year=2021)


def guthrie_cost_function(x, c1, c2, c3):
    """Compute costs or cost coefficients using the Guthrie cost function."""
    from numpy import log10
    return pow(10, c1 + c2 * log10(x) + c3 * log10(x) ** 2)


cost_function = define_function('cost_turton', guthrie_cost_function)


def lmtd_implementation(dT1, dT2):
    """Calculate the logarithmic mean temperature difference."""
    from numpy import log
    try:
        return (dT1 - dT2) / log(dT1 / dT2)
    except RuntimeWarning:
        return dT1


lmtd = define_function('lmtd', lmtd_implementation)


def rlmtd_implementation(dT1, dT2):
    """Calculate the inverse of the logarithmic mean temperature difference."""
    from numpy import log
    try:
        return log(dT1 / dT2) / (dT1 - dT2)
    except RuntimeWarning:
        return 1/dT1


rlmtd = define_function('rlmtd', rlmtd_implementation)


# NOTE: A substitute for the floor function that can be implemented using step
#       functions
def floor_substitute_implementation(x, LB, UB):
    """Calculate the floor of x for lb <= x <= ub.

    Returns lb for values below lb and x-1 for values above ub.
    """
    from numpy import floor, heaviside
    lb = int(floor(LB))
    ub = int(floor(UB))

    s = sum(heaviside(x - i, 1) for i in range(lb + 1, ub + 1))
    return max(x - 1, lb + s)

floor_substitute = define_function('floor_substitute',
                                   floor_substitute_implementation)

# FIXME: user-defined functions within other expressions cannot be evaluated!
# float(rlmtd(1, 2))  # works because it uses __float__
# float(1 + rlmtd(1, 2))  # "RuntimeError: Not Implemented" because it uses
# FIX:
# def eval(expr):
#     try:
#         return float(expr)
#     except RuntimeError:
#         return float(type(expr)(*map(eval, expr.args)))
# eval(rlmtd(1, 2))
# eval(1 + rlmtd(1, 2))


def F_even(R, S):
    """The F factor for an even number of tube passes."""
    eps = 1e-5
    S = Max(eps, Min(1-eps, S))
    R = Max(eps, Min(1/S, R))
    a = (R ** 2 + 1) ** 0.5
    numer = a * log((1 - S)/Max(eps, 1 - R * S))
    # numer = a * log(Max(eps, (1 - S)/(1 - R * S)))
    # numer = a * log(Max(eps, (1/Max(eps, S) - 1) / Max(eps, 1 - R)))
    # denom = (R - 1) * log((2 - S * (R + 1 - a)) / (2 - S * (R + 1 + a)))
    # denom = (R - 1) * log((2 - S * (R + 1 - a)) / Max(eps, 2 - S * (R + 1 + a)))
    denom = (R - 1) * log(Max(eps, (2 - S * (R + 1 - a))
                              / Max(eps, 2 - S * (R + 1 + a))))
    return numer / Max(eps, denom)


def F_sym(RS, S):
    """Symmetric form of F_even"""
    # from numpy import log
    # if abs(S) < 1e-5 or abs(RS) < 1e-5:
    #     return 1
    # if abs(S - 1) < 1e-5 or abs(RS - 1) < 1e-5 or S >= S_max(RS/S):
    #     return float('-inf')
    # if abs(RS - S) < 1e-5:
    #     sqrt2 = 2 ** 0.5
    #     return S * sqrt2 / (1 - S) / log((2 - S * (2 - sqrt2)) / (2 - S * (2 + sqrt2)))
    R = RS/S
    a = (R ** 2 + 1) ** 0.5
    log_denom = (1 - RS)
    numer = a * comando.log((1 - S)/log_denom)
    denom = (R - 1) * comando.log((2 - S * (R + 1 - a)) / (2 - S * (R + 1 + a)))
    return numer / denom


def F90_ANN_sym(RS, S):
    """Compute F for values of RS and S satisfying the 90% limit for F.

    The temperature correction factor F is used for adjusting the logarithmic
    mean temperature differnce to the mean temperature difference suitable for
    a given heat exchanger geometry.
    The 90% limit for F is given by the equation:

        RS <= (10 * S - 9)/(5 * (S - 2))

    or equivalently.

        S <= (10 * RS - 9)/(5 * (RS - 2))

    The present correlation was obtained from training an ANN with one hidden
    layer containing three neurons on data satisfying the 90% limit.
    For data at this limit, the maximum error is around 1%, while for data
    outside it, errors can get larger.
    """
    return (1.98485878182041
            - 1.04254 * tanh(3.24816 + 2.041391 * RS - 2.089098 * S)
            - 1.04254 * tanh(3.24816 - 2.089098 * RS + 2.041391 * S)
            + 1.09428753535601 * tanh(3.63440625035811
                                      - 2.48316976 * RS
                                      - 2.48316976 * S))


def F_inv_90_ANN_sym(RS, S):
    """ANN fit using 3 Neurons of the """
    return (1.80795318691684
            + 0.298978006443871 * tanh(3.07419713556056
                                       - 2.79020407672056 * RS
                                       + 1.50401518707256 * S)
            + 0.298978008519322 * tanh(3.07419713810265
                                       + 1.50401518532696 * RS
                                       - 2.79020407480942 * S)
            - 1.40088081913597 * tanh(3.98604272647425
                                      - 2.82401213955892 * RS
                                      - 2.82401213961128 * S))


def F_nrtp(p, q, n=2):
    """Temperature correction factor for chross flow with n tube rows & passes.

    Arguments
    ---------
    p = (T_{t,in} - T_{t,out}) / (T_{t,in} - T_{air,in})
    q = (T_{air,out} - T_{air,in}) / (T_{t,in} - T_{air,in})

    Publication:
        Roetzel and Nicole 1975: "Mean Temperature Difference for Heat
        Exchanger Design - A General Approximate Explicit Equation"
    """
    from numpy import log as ln, sin, arctan
    R = p / q
    r = (p - q) / ln((1 - q) / (1 - p))

    # k-1, i-1
    try:
        alpha = {
            1: [
                [-4.62e-1, -3.13e-2, -1.74e-1, -4.20e-2],
                [+5.08e+0, +5.29e-1, +1.32e+0, +3.47e-1],
                [-1.57e+1, -2.37e+0, -2.93e+0, -9.53e-1],
                [+1.72e+1, +3.18e+0, +1.99e+0, +6.49e-1]
            ],
            2: [
                [-2.35e-1, -7.73e-2, -5.98e-2, +5.25e-3],
                [+2.28e+0, +6.32e-1, +3.64e-1, -1.27e-2],
                [-6.44e+0, -1.63e+0, -6.13e-1, -1.14e-2],
                [+6.24e+0, +1.35e+0, +2.76e-1, +2.72e-2]
            ],
            3: [
                [-8.43e-1, +3.02e-2, +4.80e-1, +8.12e-2],
                [+5.85e+0, -9.64e-3, -3.28e+0, -8.34e-1],
                [-1.28e+1, -2.28e-1, +7.11e+0, +2.19e+0],
                [+9.24e+0, +2.66e-1, -4.90e+0, -1.69e+0]
            ],
            4: [
                [-3.39e-1, +2.77e-2, +1.79e-1, -1.99e-2],
                [+2.38e+0, -9.99e-2, -1.21e+0, +4.00e-2],
                [-5.26e+0, +9.04e-2, +2.62e+0, +4.94e-2],
                [+3.90e+0, -8.45e-4, -1.81e+0, -9.81e-2]
            ]
        }[n]
    except KeyError:
        raise ValueError(f"Parameter n must be 1, 2, 3, or 4, was {n}!")
    return 1 - sum(sum(alpha[k][i] * (1 - r) ** (k+1)
                       * sin(2 * (i + 1) * arctan(R)) for k in range(4))
                   for i in range(4))


# NOTE: p == RS, q == S
def F_Schedwill(p, q):
    """Scheddwill correlation for a singele finned tube in cross flow.

    Taken from
    Roetzel and Neubert 1978: "Calculation of Mean Temperature Difference in
    Air-CooIed Cross-Flow Heat Exchangers"

    Arguments
    ---------
    p = (T_{t,in} - T_{t,out}) / (T_{t,in} - T_{air,in})
    q = (T_{air,out} - T_{air,in}) / (T_{t,in} - T_{air,in})
    """
    from numpy import log as ln
    r = - q / ln(1 + q / p * ln(1 - p))
    rln = (p - q) / ln((1 - q) / (1 - p))
    return r / rln


def F_Schedwill_lim(p, fac=1):
    """Schedwill limit in terms of the tube-side effectiveness.

    Describes the boundary at which F_Schedwill approaches negative infinity.

    Arguments
    ---------
    p = (T_{t,in} - T_{t,out}) / (T_{t,in} - T_{air,in})
    fac a reduction factor in (0, 1] to shift the limit towards the origin
    """
    from numpy import log as ln
    return -p / ln(1 - p / fac)


def F_Schedwill_ANN_inv(p, q):
    """Fit of 1/F_Schedwill for single finned tube in cross-flow.

    Valid for data in the 90% region of the Schedwill limit with an accuracy of
    approximately 1%.
    """
    return (0.601270552853066
            - 1.98431035170102 * tanh(-2.31031107792008 * p
                                      - 1.69787794706641 * q
                                      + 3.15952494927165)
            + 1.98431035170102 * tanh(2.01956880315572 * p
                                      - 2.38366995223538 * q
                                      + 3.80961342263366)
            - 0.224468430494661 * tanh(-6.51019861125896 * p
                                       - 9.34554507220272 * q
                                       + 9.60235524498431)
            + 0.619159197037758 * tanh(-2.76228171305818 * p
                                       - 1.21562742194533 * q
                                       + 2.85878849610653))


class HeatExchanger(comando.Component):
    """A model for a counterflow heat exchanger.

    We describe the heat flows experienced by the hot and cold fluid either via
    mass flow and specifiv enthalpies at in and outlet or heat capacity flow
    and temperatures at the in and outlet.
    As we assume an adiabatic heat exchanger, the two heat flows must be
    identical. We therefore allow for one of the six unknowns describing the
    two heat flows to be expressed in terms of the other five. If all six are
    given, we introduce an equality constraint forcing heat flow identity.

    .. note:
    With an enthalpy-based description, temperatures are also assumed to be
    given in order to impose constraints on temperature differences!

    Arguments
    ---------
    label : `str`
        A unique sting that serves as an identifier

    dT_min : expression
        The minimum temperature difference that needs to be upheld

    T_h_in : expression
        The temperature at the input of the hot side

    T_h_out : expression
        The temperature at the output of the hot side

    T_c_in : expression
        The temperature at the input of the cold side

    T_c_out : expression
        The temperature at the output of the cold side

    h_h_in : expression
        The specific enthalpy at the input of the hot side

    h_h_out : expression
        The specific enthalpy at the output of the hot side

    h_c_in : expression
        The specific enthalpy at the input of the cold side

    h_c_out : expression
        The specific enthalpy at the output of the cold side

    mdot_h : expression
        The mass flow passing through the hot side

    mdot_c : expression
        The mass flow passing through the cold side

    mdot_cp_h : expression
        The mean heat capacity flow passing through the hot side

    mdot_cp_c : expression
        The mean heat capacity flow passing through the cold side

    F_p : expression
        Coefficient for investment cost calculations
    """

    def _set_delta(self, q, q1, q2, qmin=0, create_con=True):
        dq_calc = q1 - q2
        if create_con:
            try:
                qminstr = float(qmin)
            except RuntimeError:
                qminstr = f'd{q}_min'
            if q1 == q2 or dq_calc == qmin:
                print(f'INFO: Skipping trivial delta constraint {self.label}_d{q}!')
            else:
                self.add_ge_constraint(dq_calc, qmin, f'd{q} >= {qminstr}')
        return self.add_expression(f'd{q}', Max(qmin, dq_calc))

    def __init__(self, label, dT_min, **kwargs):
        super().__init__(label)

        exprs = self._expressions_dict
        _set_delta = self._set_delta
        self.dT_min = dT_min
        def alldefined(*args):
            """Check if all arguments are specified, if yes store them."""
            if all(arg in kwargs for arg in args):
                for arg in args:
                    exprs[arg] = kwargs[arg]
                return True

        def T_based_Hdot(x):
            """Set hot/cold enthalpy flows based on T and return dT."""
            if x == 'h':
                return _set_delta('T_h', exprs['T_h_in'], exprs['T_h_out'])
            else:
                return _set_delta('T_c', exprs['T_c_out'], exprs['T_c_in'])

        def h_based_Hdot(x):
            """Set hot/cold enthalpy flows based on h and return dh."""
            if not alldefined(f'T_{x}_in', f'T_{x}_out'):
                raise RuntimeError(f'Temperatures T_{x}_in and T_{x}_out must '
                                   'be given!')
            if x == 'h':
                _set_delta('T_h', exprs['T_h_in'], exprs['T_h_out'])
                return _set_delta('h_h', exprs['h_h_in'], exprs['h_h_out'])
            else:
                _set_delta('T_c', exprs['T_c_out'], exprs['T_c_in'])
                return _set_delta('h_c', exprs['h_c_out'], exprs['h_c_in'])

        def hotQ(undefined):
            if alldefined('T_h_in', 'T_h_out', 'mdot_cp_h'):
                Qdot_h = exprs['Qdot'] = exprs['mdot_cp_h'] * T_based_Hdot('h')
            elif alldefined('h_h_in', 'h_h_out', 'mdot_h', 'T_h_in',
                            'T_h_out'):
                Qdot_h = exprs['Qdot'] = exprs['mdot_h'] * h_based_Hdot('h')
            else:
                raise RuntimeError(f'If {undefined} is not given, you must '
                                   'specify T_h_in and T_h_out, and either '
                                   'mdot_cp_h or h_h_in, h_h_out and mdot_h!')
            return Qdot_h

        def coolQ(undefined):
            if alldefined('T_c_in', 'T_c_out', 'mdot_cp_c'):
                Qdot_c = exprs['Qdot'] = exprs['mdot_cp_c'] * T_based_Hdot('c')
            elif alldefined('h_c_in', 'h_c_out', 'mdot_c', 'T_c_in',
                            'T_c_out'):
                Qdot_c = exprs['Qdot'] = exprs['mdot_c'] * h_based_Hdot('c')
            else:
                raise RuntimeError(f'If {undefined} is undefined, you must '
                                   'specify T_c_in and T_c_out, and either '
                                   'mdot_cp_c or h_c_in, h_c_out and mdot_c!')
            return Qdot_c

        # We start looking for undefined quantities on the hot side
        if 'h_h_in' not in kwargs and alldefined('h_h_out', 'mdot_h'):
            Qdot_c = coolQ()
            exprs['h_h_in'] = exprs['h_h_out'] + Qdot_c / exprs['mdot_h']
            h_based_Hdot('h')
        elif 'h_h_out' not in kwargs and alldefined('h_h_in', 'mdot_h'):
            Qdot_c = coolQ('h_h_out')
            exprs['h_h_out'] = exprs['h_h_in'] - Qdot_c / exprs['mdot_h']
            h_based_Hdot('h')
        elif 'mdot_h' not in kwargs and alldefined('h_h_in', 'h_h_out'):
            Qdot_c = coolQ('mdot_h')
            dh_h = h_based_Hdot('h')
            exprs['mdot_h'] = Qdot_c / dh_h
        elif 'T_h_in' not in kwargs and alldefined('T_h_out', 'mdot_cp_h'):
            Qdot_c = coolQ('T_h_in')
            try:  # if cold side is defined via mdot
                exprs['T_h_in'] = (exprs['T_h_out']
                                   + exprs['mdot_c'] / exprs['mdot_cp_h']
                                   * exprs['dh_c'])
            except KeyError:  # cold side is defined via mdot_cp
                exprs['T_h_in'] = (exprs['T_h_out']
                                   + exprs['mdot_cp_c'] / exprs['mdot_cp_h']
                                   * exprs['dT_c'])
            T_based_Hdot('h')
        elif 'T_h_out' not in kwargs and alldefined('T_h_in', 'mdot_cp_h'):
            Qdot_c = coolQ('T_h_out')
            try:  # if cold side is defined via mdot
                exprs['T_h_out'] = (exprs['T_h_in']
                                    - exprs['mdot_c'] / exprs['mdot_cp_h']
                                    * exprs['dh_c'])
            except KeyError:  # cold side is defined via mdot_cp
                exprs['T_h_out'] = (exprs['T_h_in']
                                    - exprs['mdot_cp_c'] / exprs['mdot_cp_h']
                                    * exprs['dT_c'])
            T_based_Hdot('h')
        elif all(k not in kwargs for k in ('mdot_cp_h', 'h_h_in', 'h_h_out')) \
                and alldefined('T_h_in', 'T_h_out'):
            Qdot_c = coolQ('mdot_cp_h')
            dT_h = T_based_Hdot('h')
            exprs['mdot_cp_h'] = Qdot_c / dT_h
        # Now the cold side...
        elif 'h_c_in' not in kwargs and alldefined('h_c_out', 'mdot_c'):
            Qdot_h = hotQ('h_c_in')
            exprs['h_c_in'] = exprs['h_c_out'] - Qdot_h / exprs['mdot_c']
            h_based_Hdot('c')
        elif 'h_c_out' not in kwargs and alldefined('h_c_in', 'mdot_c'):
            Qdot_h = hotQ('h_c_out')
            exprs['h_c_out'] = exprs['h_c_in'] + Qdot_h / exprs['mdot_c']
            h_based_Hdot('c')
        elif 'mdot_c' not in kwargs and alldefined('h_c_in', 'h_c_out'):
            Qdot_h = hotQ('mdot_c')
            dh_c = h_based_Hdot('c')
            exprs['mdot_c'] = Qdot_h / dh_c
        elif 'T_c_in' not in kwargs and alldefined('T_c_out', 'mdot_cp_c'):
            Qdot_h = hotQ('T_c_in')
            try:  # if hot side is defined via mdot
                exprs['T_c_in'] = (exprs['T_c_out']
                                   - exprs['mdot_h']/exprs['mdot_cp_c']
                                   * exprs['dh_h'])
            except KeyError:  # hot side is defined via mdot_cp
                exprs['T_c_in'] = (exprs['T_c_out']
                                    - exprs['mdot_cp_h']/exprs['mdot_cp_c']
                                    * exprs['dT_h'])
            T_based_Hdot('c')
        elif 'T_c_out' not in kwargs and alldefined('T_c_in', 'mdot_cp_c'):
            Qdot_h = hotQ('T_c_out')
            try:  # if hot side is defined via mdot
                exprs['T_c_out'] = (exprs['T_c_in']
                                    + exprs['mdot_h']/exprs['mdot_cp_c']
                                    * exprs['dh_h'])
            except KeyError:  # hot side is defined via mdot_cp
                exprs['T_c_out'] = (exprs['T_c_in']
                                    + exprs['mdot_cp_h']/exprs['mdot_cp_c']
                                    * exprs['dT_h'])
            T_based_Hdot('c')
        elif all(k not in kwargs for k in ('mdot_cp_c', 'h_c_in', 'h_c_out')) \
                and alldefined('T_c_in', 'T_c_out'):
            Qdot_h = hotQ('mdot_cp_c')
            dT_c = T_based_Hdot('c')
            exprs['mdot_cp_c'] = Qdot_h / dT_c
        else:  # Assume all quantities are defined
            Qdot_c = coolQ('h_h_out')
            Qdot_h = hotQ('mdot_cp_c')
            self.add_eq_constraint((Qdot_c - Qdot_h), 0, 'heat_flow')

    def set_shell_and_tube_geometry(self, p_t, d_o, d_i, B_per_d_s=None, L_per_d_s=None, n_tp=1, arr='triangular', d_s_min=0.5, d_s_max=2.5):
        d_s, d_s_scaled = add_scaled_var(self, 'd_s', d_s_min, d_s_max,
                                         design=True)
        self.kind = 'S&T'
        self.d_s = self.add_expression('d_s', d_s)
        self.d_o = self.add_expression('d_o', d_o)
        self.d_i = self.add_expression('d_i', d_i)
        self.p_t = self.add_expression('p_t', p_t)  # tube pitch [m]
        # B/d_s should lie between 0.2 to 1 (TEMA, HEDH)
        if B_per_d_s is None:
            B_per_d_s, _ = add_scaled_var(self, 'B_per_d_s', 0.2, 1,
                                          design=True)
        self.B = self.add_expression('B', B_per_d_s * d_s)

        # tube bundle diameter [m]
        bundle_gap = 0.005 + 0.012 / d_s
        self.d_tb = d_tb = self.add_expression('d_tb', d_s - bundle_gap)

        if arr == 'triangular':
            c = 3 ** 0.5 / 2
        elif arr == 'square':
            c = 1
        # equivalent wetted shell diameter
        self.d_e = 4 * c * p_t ** 2 / numpy.pi / d_o - d_o
        d_ctl = d_tb - d_o
        self.n_t = n_t = self.add_expression('n_t',
                                             0.78 * (d_ctl / p_t) ** 2 / c)

        self.arr = arr

        self.n_tp = self.add_expression('n_tp', n_tp)
        # HEDH: assume L/d_s = 8 initially
        if L_per_d_s is None:
            L_per_d_s, _ = add_scaled_var(self, 'L_per_d_s', 4, 12,
                                          design=True)
        self.L_t = L_t = self.add_expression('L_t', L_per_d_s * d_s)
        self.A = self.add_expression('A', numpy.pi * d_o * L_t * n_t)  # A_HX
        self.A_tcs = self.add_expression(
            'A_tcs', numpy.pi * d_i ** 2 / 4 * (n_t / n_tp))

        def limit_tube_velocity(Vdot_max, vel_max):
            """Add a constraint to limit tube velocity."""
            self.add_expression('vel_t_max', Vdot_max / self.A_tcs)
            self.add_le_constraint(Vdot_max, self.A_tcs * vel_max,
                                   f'vel_t <= {vel_max} m/s')

        self.limit_tube_velocity = limit_tube_velocity

        # 6 / 254 in serth2007process
        self.A_scs = self.add_expression(
            'A_scs', self.B * (bundle_gap + (d_tb - d_o) * (1 - d_o / p_t)))
            # TODO: Simplification
            # 'A_scs', self.B * ((d_s - d_tb) + (d_tb - d_o) * (1 - d_o / p_t)))
            # 'A_scs', self.B * ((d_s - d_tb) + (d_tb - d_o) - (d_tb - d_o) * (d_o / p_t)))
            # 'A_scs', self.B * ((d_s - d_o) - (d_tb - d_o) * (d_o / p_t)))
            # 'A_scs', self.B * (d_s - d_o * (1 + (d_tb - d_o) / p_t)))

        def limit_shell_velocity(Vdot_max, vel_max):
            """Add a constraint to limit shell velocity."""
            self.add_expression('vel_s_max', Vdot_max / self.A_scs)
            self.add_le_constraint(Vdot_max, self.A_scs * vel_max,
                                   f'vel_s <= {vel_max} m/s')

        self.limit_shell_velocity = limit_shell_velocity

        self.set_F()

        return self.A

    def set_ACC_geometry(self, A_min, A_max, T_c_pinch, T_h_pinch, Qdot_con, eta_F=0.7):
        """Specify the heat exchanged to be an ACC."""
        self.kind = 'ACC'
        self.eta_F = eta_F
        A, A_scaled = add_scaled_var(self, 'A', A_min, A_max, design=True)

        self.n_tp = n_tp = self.add_expression('n_tp', 4)
        self.arr = arr = 'triangular'

        self.n_r = n_tp  # number of tube rows = number of tube passes

        # ghasemi2013modeling
        n_fans_per_bay = 1

        L_t_val = 10.97  # == 36 feet
        # L_t, L_t_scaled = add_scaled_var(self, 'L_t', 4, 20, design=True)
        self.L_t = L_t = self.add_expression('L_t', L_t_val)
        self.p_t = p_t = self.add_expression('p_t', 0.06985)  # 2.75 in
        self.d_f = d_f = self.add_expression('d_f', 0.0635)  # 2.5 in
        self.d_o = d_o = self.add_expression('d_o', 0.03175)  # 1.25 in
        self.d_i = d_i = self.add_expression('d_i', 0.02967)  # BWG 14
        self.t_f = t_f = self.add_expression('t_f', 0.00041)
        self.s_f = s_f = self.add_expression('s_f', 0.0019)  # fin spacing
        # number of fins per length
        self.n_f = n_f = self.add_expression('n_f', 1/(s_f + t_f))

        # Serth Table 12.1 (1 in tube 2 in fin)
        #     d_f = 0.0508  # 2 in
        #     d_o = 0.0254  # 1 in
        #     d_i = 0.021184  # BWG 14
        #     p_t = 0.05715  # 2.25 in
        #     t_f = 0.0003302  # 0.013 in
        #     n_f = 9 / 0.0254  # 9 per inch
        #     arr='triangular'

        self.n_t = n_t = self.add_expression(
            'n_t',
            A / (L_t * numpy.pi * (d_o + n_f * ((d_f ** 2 - d_o ** 2) / 2
                                                + t_f * (d_f - d_o))))
        )

        self.A_o = self.add_expression('A_o', numpy.pi * d_o * L_t * n_t)

        d_f_o = d_f / d_o
        self.A_per_A_o = A_per_A_o = self.add_expression(
            'A_per_A_o', 1 + n_f * ((d_f_o ** 2 - 1) * d_o / 2
                                    + (d_f_o - 1) * t_f)
        )

        d_f_i = d_f / d_i
        d_o_i = d_o / d_i
        self.A_per_A_i = self.add_expression(
            'A_per_A_i', d_o_i + n_f * ((d_f_i ** 2 - d_o_i ** 2) * d_i / 2
                                        + (d_f_i - d_o_i) * t_f)
        )

        self.A_tcs = self.add_expression(
            'A_tcs',
            numpy.pi * d_i ** 2 / 4 / n_tp * self.n_t
        )

        def limit_tube_velocity(Vdot_max, vel_t_max, suffix=''):
            """Add a constraint to limit tube velocity."""
            self.add_expression(f'vel_t{suffix}_max', Vdot_max / self.A_tcs)
            self.add_le_constraint(Vdot_max, self.A_tcs * vel_t_max,
                                   f'vel_t{suffix} <= {vel_t_max} m/s')

        self.limit_tube_velocity = limit_tube_velocity

        if arr == 'triangular':
            self.W = self.add_expression('W', 3 ** 0.5 / 2 * n_t / n_tp * p_t)
        elif arr == 'square':
            self.W = self.add_expression('W', n_t / n_tp * p_t)
        else:
            raise ValueError

        self.A_face = self.add_expression('A_face', self.W * L_t)


        from numpy import pi, floor, ceil

        # Bay size, assuming d_F <= W_B <= L_B
        # L_B = L_t
        # L_B * W_B / n_fans_per_bay * 0.4 <= pi * d_F ** 2 / 4  # smallest bay size
        # L_B * W_B / n_fans_per_bay * 0.4 <= pi * W_B ** 2 / 4  # min W_B @ d_F == W_B <= L_B
        # L_B / n_fans_per_bay * 1.6 / pi <= W_B
        W_bay_min = L_t / n_fans_per_bay * 1.6 / pi
        # L_B * W_B / n_fans_per_bay * 0.4 <= pi * L_B ** 2 / 4  # max W_B @ d_F == L_B <= W_B
        # W_B * 0.4 <= 10 * pi * L_B * n_fans_per_bay
        W_bay_max = 10 * pi * L_t * n_fans_per_bay

        self.n_banks = n_banks = 2

        W_bank = self.add_expression('W_bank', self.W / self.n_banks)

        W_bank_min = A_min / A_per_A_o / L_t / n_banks
        W_bank_max = A_max / A_per_A_o / L_t / n_banks
        n_bays_per_bank_min = max(4, ceil(W_bank_min / W_bay_max))
        n_bays_per_bank_max = floor(W_bank_max / W_bay_min)

        print(W_bay_min, W_bay_max)
        print(n_bays_per_bank_min, n_bays_per_bank_max)
        # exit()

        W_bay_approx = L_t / n_fans_per_bay
        n_bay_per_bank = floor_substitute(W_bank / W_bay_approx,
                                          n_bays_per_bank_min,
                                          n_bays_per_bank_max)

        self.n_bays = n_bays = self.add_expression('n_bays',
                                                   n_bay_per_bank
                                                   * self.n_banks)
        self.n_fans = self.add_expression('n_fans', n_bays * n_fans_per_bay)

        W_bay = self.add_expression('W_bay', W_bank / n_bay_per_bank)
        self.add_le_constraint(W_bay_min, W_bay, f'{W_bay_min} m <= W_bay')
        self.add_le_constraint(W_bay, W_bay_max, f'W_bay <= {W_bay_max} m')

        A_bay = self.add_expression('A_bay', W_bay * L_t)
        A_F = self.add_expression('A_F', 0.4 * A_bay / n_fans_per_bay)
        d_F = self.add_expression('d_F', (4 / pi * A_F) ** 0.5)

        self.add_le_constraint(n_fans_per_bay * d_F, L_t,
                               f'n_fans_per_bay * d_F <= {L_t} m')

        self.set_F(T_c_in=T_c_pinch, T_h_out=T_h_pinch, suffix='_DES')
        self.set_F(T_c_out=T_c_pinch, T_h_in=T_h_pinch, suffix='_CON')
        Qdot = self.get_expression('Qdot')
        self.add_expression('Qdot_CON', Qdot_con)
        self.add_expression('Qdot_DES', Max(0, Qdot - Qdot_con))

    def set_F(self, T_c_in=None, T_c_out=None, T_h_in=None, T_h_out=None, suffix=''):
        # Deltas for temperatures
        if T_c_in is None:
            T_c_in = self.get_expression('T_c_in')
        if T_c_out is None:
            T_c_out = self.get_expression('T_c_out')
        if T_h_in is None:
            T_h_in = self.get_expression('T_h_in')
        if T_h_out is None:
            T_h_out = self.get_expression('T_h_out')

        try:
            dT_min_val = f'{float(self.dT_min)} K'
        except RuntimeError:
            dT_min_val = 'dT_min'

        dT_c = self._set_delta(f'T_c{suffix}', T_c_out, T_c_in)

        dT_h = self._set_delta(f'T_h{suffix}', T_h_in, T_h_out)

        dT_h_in = self._set_delta(f'T_h_in{suffix}', T_h_in, T_c_out,
                                  self.dT_min)

        dT_h_out = self._set_delta(f'T_h_out{suffix}', T_h_out, T_c_in,
                                   self.dT_min)

        dT_lim = self._set_delta(f'T_lim{suffix}', T_h_in, T_c_in,
                                 self.dT_min, False)

        # Heat capacity ratio
        R = self.add_expression(f'R{suffix}', dT_h / dT_c)  # can be 1 / 0 !!!
        # Hot effectiveness
        RS = self.add_expression(f'RS{suffix}', Max(0, Min(1, dT_h / dT_lim)))
        # Cold effectiveness
        S = self.add_expression(f'S{suffix}', Max(0, Min(1, dT_c / dT_lim)))

        # NOTE: The analytical correlations for F approach negative infinity
        #       for certain combinations of S and RS. To avoid steep slopes and
        #       the associated uncertainty in the value of F, we restrict
        #       ourselves to 90% of this limiting value.
        fac = 0.90
        isothermal = T_c_in == T_c_out or T_h_in == T_h_out
        if isothermal:
            F_inv_calc = 1
        elif self.kind == 'ACC':
            # NOTE: Here we assume that the desuperheating mostly takes place
            #       within the first tube row and consequently use Schedwill's
            #       correlation for a single finned tube in crossflow.
            F_inv_calc = F_Schedwill_ANN_inv(RS, S)

            # Update RS
            self.add_le_constraint(dT_h, fac * dT_lim,
                                   f'dT_h <= {fac} * dT_lim')
            self.add_le_constraint(S,
                                   -RS / comando.log(Max(1e-9, Min(1 - 1e-9, 1 - RS / fac))),
                                   "S <= 90_percent_Schedwill_limit")
        elif self.kind == 'S&T':  # We always assume a single shell pass!
            if self.n_tp == 1:  # Assuming cross-flow configuration
                F_inv_calc = 1
            elif self.n_tp % 2 == 1:  # Correlation for even # of tube passes
                F_inv_calc = F_inv_90_ANN_sym(RS, S)

                # NOTE: As suggested in smith2005chemical instead of limiting
                #       F, we can derive limits on S and RS.
                self.add_le_constraint(2 * RS + 2 * S - RS * S, 2 * fac,
                                       f'F{suffix} <= {fac} F_lim')
                # Additional individual limits
                # 2 * RS + 2 * S - RS * S ≤ 2 * fac
                # 2 * RS + S * (2  - RS) ≤ 2 * fac
                # S * (2  - RS) ≤ 2 * (fac - RS)
                # S ≤ 2 * (fac - RS) / (2 - RS)
                self.add_le_constraint(S, 2 * (fac - RS) / (2 - RS),
                                       'S <= Smax')
                self.add_le_constraint(RS, 2 * (fac - S) / (2 - S),
                                       'RS <= RSmax')
        else:
            raise NotImplementedError("No F correlation implemented for heat "
                                      f"exchanger kind '{self.kind}'!")

        F_inv = self.add_expression(f'F_inv{suffix}',
                                    Max(1, Min(10/7, F_inv_calc)))
        self.add_expression(f'F{suffix}', 1/F_inv)

        self.add_expression(f'lmtd{suffix}',
                            lmtd(Max(self.dT_min, dT_h_in),
                                 Max(self.dT_min, dT_h_out)))
        self.add_expression(f'rlmtd{suffix}',
                            rlmtd(Max(self.dT_min, dT_h_in),
                                  Max(self.dT_min, dT_h_out)))

        # Average temperature of cold fluid, hot fluid, and wall and average
        # temperature difference
        T_c = self.add_expression(f'T_c{suffix}', 0.5 * (T_c_in + T_c_out))
        T_h = self.add_expression(f'T_h{suffix}', 0.5 * (T_h_in + T_h_out))
        T_w = self.add_expression(f'T_w{suffix}', 0.5 * (T_c + T_h))

        self._set_delta(f'T{suffix}', T_h, T_c, self.dT_min)

    def calc_U_inv(self, k_t=16, k_f=200, R=1.3e-4, R_con=0, evaporating=False, condensing=False, **properties):
        """Calculate the inverse of U.

        It is assumed that the hot fluid is always in the tube-side.

        Arguments
        ---------
        k_t : Expression
            heat conductivity of tube material [W/m/K]
        R : Expression
            shell-side fouling factor [m²K/W]
        evaporating : bool
            whether the cold fluid is being evaporated
        properties : dict
            fluid properties, must include:

            - For hot fluid:
                mdot_h, v_h, mu_h, Pr_h, k_h
                if condensing additionally: p_h, p_crit_h
            - For cold fluid:
                if condensing: mdot_c, v_in_c, rho_c, mu_c, Pr_c, k_c, vel_face_max_c
                if evaporating: p_c, p_crit_c, T_sat_c
                else: mdot_c, mu_c, mu_w_c, Pr_c, k_c

        Returns
        -------
        U_inv : Expression
            inverse of total heat transfer coefficient
        """
        pi = numpy.pi
        # Shorthands
        g = self  # geometry
        from types import SimpleNamespace
        c = SimpleNamespace()
        h = SimpleNamespace()
        for prop, val in properties.items():
            if prop.endswith('_c'):
                setattr(c, prop[:-2], val)
            elif prop.endswith('_h'):
                setattr(h, prop[:-2], val)
            else:
                raise ValueError(f'Unexpected property {prop}')

        # Hot fluid (tube side)
        mdot_t = h.mdot * g.n_tp / g.n_t
        Vdot_t = mdot_t * h.v
        vel_t = self.add_expression(
            'vel_t_calc', Vdot_t / (pi * g.d_i ** 2 / 4))

        if condensing:
            # NOTE: Re_t and h.Pr should be evaluated for saturated liqud!
            # Re = d_i * vel_t * rho / mu
            #    = d_i * 4 * Vdot_t / (pi * d_i ** 2) * rho / mu
            #    = 4 * Vdot_t / (pi * d_i) * rho / mu
            #    = 4 * mdot_t / (pi * d_i) / mu
            Re_t_CON_calc = self.add_expression(
                'Re_t_CON',
                4 * mdot_t / (pi * g.d_i) * 1 / Max(1e-9, h.mu_CON)
            )
            self.add_ge_constraint(Re_t_CON_calc, 2100, 'Re_t_ >= 2100')
            Re_t_CON = Max(2100, Re_t_CON_calc)

            self.add_ge_constraint(h.Pr_CON, 0.6, 'Pr_h_CON >= 0.6')
            Pr_h_CON = Max(0.6, h.Pr_CON)

            # coefficient from [Shah 1979]
            # C_int = (0.55 + 2.09 / p_r ** 0.38)
            # Evaluation of simplification suggested in [Shah 1979]
            # C_int = (0.555555555555556 + 2.10401844532279 / p_r ** 0.38)
            # Analytical integration of C assuming x decreases linearly
            C_int = (2.0433943842548 / (h.p / h.p_crit) ** 0.38 + 5/9)
            Nu_t_CON = self.add_expression(
                'Nu_t_CON',
                0.023 * Re_t_CON ** 0.8 * Pr_h_CON ** 0.4 * C_int
            )
            h_t_CON_calc = self.add_expression('alpha_t_CON',
                                               Nu_t_CON * h.k_CON / g.d_i)
            # self.add_ge_constraint(h_t_CON_calc, 100, 'h_t_CON >= 100')
            h_t_CON = Max(100, h_t_CON_calc)

            Re_t_DES_calc = self.add_expression(
                'Re_t_DES',
                4 * mdot_t / (pi * g.d_i) * 1 / Max(1e-9, h.mu_DES)
            )
            self.add_ge_constraint(Re_t_DES_calc, 2100, 'Re_t_DES >= 2100')
            Re_t_DES = Max(2100, Re_t_DES_calc)

            self.add_ge_constraint(h.Pr_DES, 0.6, 'Pr_h_DES >= 0.6')
            Pr_h_DES = Max(0.6, h.Pr_DES)
            # NOTE: Gnielinski correlation is only valid for fully turbulent
            #       flow (2100 ≤ Re_t ≤ 1e6) and 0.6 ≤ Pr ≤ 2000
            f_over_8 = Max(1.087, 0.782 * comando.log(Re_t_DES) - 1.51) ** -2 / 8
            Nu_t_DES = self.add_expression(
                'Nu_t_DES',
                f_over_8 * (Re_t_DES - 1000) * Pr_h_DES
                / (1 + 12.7 * f_over_8 ** 0.5 * (Pr_h_DES ** (2/3) - 1))
            )
            h_t_DES_calc = self.add_expression('alpha_t_DES',
                                           Nu_t_DES * h.k_DES / g.d_i)
            # self.add_ge_constraint(h_t_DES_calc, 100, 'h_t_DES >= 100')
            h_t_DES = Max(100, h_t_DES_calc)
        else:  # Single phase (no condensation)
            # Re = d_i * vel_t * rho / mu
            #    = d_i * 4 * Vdot_t / (pi * d_i ** 2) * rho / mu
            #    = 4 * Vdot_t / (pi * d_i) * rho / mu
            #    = 4 * mdot_t / (pi * d_i) / mu
            Re_t_calc = self.add_expression(
                'Re_t',
                4 * mdot_t / (pi * g.d_i) * 1 / Max(1e-9, h.mu)
            )
            self.add_ge_constraint(Re_t_calc, 2100, 'Re_t >= 2100')
            Re_t = Max(2100, Re_t_calc)

            self.add_ge_constraint(h.Pr, 0.6, 'Pr_h >= 0.6')
            Pr_h = Max(0.6, h.Pr)

            f_over_8 = Max(1.087, 0.782 * comando.log(Re_t) - 1.51) ** -2 / 8

            # NOTE: Gnielinski correlation is only valid for fully turbulent
            #       flow (2100 ≤ Re_t ≤ 1e6) and 0.6 ≤ Pr ≤ 2000
            Nu_t = self.add_expression(
                'Nu_t_Gnielinski',
                f_over_8 * (Re_t - 1000) * Pr_h * (1 + g.d_i / g.L_t) ** (2/3)
                / (1 + 12.7 * f_over_8 ** 0.5 * (Pr_h ** (2/3) - 1))
            )

            # NOTE: Seider, Tate, & Hausen correlation (Re >= 1e4)
            # ratio for reasonable mu values won't be too large
            # mu_h_frac = Max(0.1, Max(1e-9, h.mu) / Max(1e-9, h.mu_w))
            # self.add_expression('mu_h', h.mu)
            # self.add_expression('mu_h_w', h.mu_w)
            # self.add_expression('mu_h_frac', h.mu / h.mu_w)
            #
            # Nu_t = self.add_expression(
            #     'Nu_t_SeiderTateHausen',
            #     0.023 * Re_t ** 0.8 * Pr_h ** (1 / 3) * mu_h_frac ** 0.14
            # )
            h_t_calc = self.add_expression('alpha_t', Nu_t * h.k / g.d_i)
            # self.add_ge_constraint(h_t_calc, 100, 'h_t >= 100')
            h_t = Max(100, h_t_calc)

        # Cold fluid
        if condensing:  # (air over finned tubes)
            # mdot_cp_air = c.mdot * c.cp
            Vdot_air_in = self.add_expression('Vdot_air_in',
                                              c.mdot * c.v_in)
            # Face inlet velocity
            vel_face = self.add_expression('vel_face',
                                           Vdot_air_in / g.A_face)
            # self.add_le_constraint(vel_face, c.vel_face_max,
            #                        f'vel_face <= {c.vel_face_max}')

            vel_air_max = self.add_expression(
                'vel_air_max',
                vel_face * g.p_t / (g.p_t - g.d_o - g.n_f * g.d_f * g.t_f)
            )
            Re_air_calc = self.add_expression(
                'Re_air',
                c.mdot * c.v_in * c.rho / (g.A_face * Max(1e-6, c.mu))
                * g.d_o * g.p_t / (g.p_t - g.d_o - g.n_f * g.d_f * g.t_f)
            )
            self.add_ge_constraint(Re_air_calc, 1000, 'Re_air >= 1000')
            Re_air = Max(1000, Re_air_calc)

            # Intermediate calculation for air side pressure drop
            a = (g.p_t - g.d_f) / g.d_o
            b = (g.d_f - g.d_o) / 2  # fin height (from tube to fin tip)
            Re_eff_inv = 1 / Re_air * b / g.s_f
            # Re ≤ 18e3
            # -> Re_eff = Re * s_f / b ≤ 50e3 * s_f / b
            # -> 1/(50e3 * s_f / b) = b / (50e3 * s_f) ≤ Re_eff_inv
            self.add_le_constraint(b / (18e3 * g.s_f), Re_eff_inv,
                                   'Re_air <= 18e3')
            Re_eff_inv = Max(Re_eff_inv, b / (18e3 * g.s_f))
            # Friction factor of Ganguli et al.
            f = (1 + 2 * numpy.exp(-a / 4) / (1 + a)) \
                * (0.021 + 27.2 * Re_eff_inv + 0.29 * Re_eff_inv ** 0.2)
            # Robinson & Briggs
            # f = 18.93 * Re_air_inv ** (0.316) * (p_t / d_o) ** (-0.927)
            self.add_expression(
                'Delta_p_air', 2 * f * g.n_r * vel_air_max ** 2 * c.rho
            )  # should be below about 500 Pa
            self.add_expression('Vdot_air', c.mdot * c.v)
            # +20% of the bundle pressure drop as a conservative estimate
            self.add_expression(  # P_req = 1.2 * Vdot_air * Delta_p_air / eta_F
                'P_req',
                1.2 * 2 * f * g.n_r * vel_air_max ** 2 * c.mdot / g.eta_F
            )

            self.add_ge_constraint(c.Pr, 0.6, 'Pr_c >= 0.6')
            Pr_c = Max(0.6, self.add_expression('Pr_air', c.Pr))

            Nu_air = self.add_expression(
                'Nu_air',
                0.38 * Re_air ** 0.6 * Pr_c ** (1/3) * g.A_per_A_o ** -0.15
            )
            h_air_calc = self.add_expression('alpha_air', Nu_air * c.k / g.d_o)
            h_air = Max(1, h_air_calc)

            psi = (g.d_f + g.t_f - g.d_o) \
                / 2 * (1 + 0.35 * numpy.log((g.d_f + g.t_f) / g.d_o))
            m = (2 * h_air / k_f / g.t_f) ** 0.5
            arg = m * psi
            eta_f = comando.tanh(arg) / arg

            # weighted fin efficiency (eliminating pi * L_t)
            # A_fin = 2 * pi / 4 * (d_f ** 2 - d_o ** 2) + pi * d_f * t_f
            # A_fins = A_fin * n_f * L_t
            #        = A_fins_piL * pi * L_t
            A_fins_piL = g.n_f * ((g.d_f ** 2 - g.d_o ** 2) / 2
                                  + g.d_f * g.t_f)

            # A_prime = pi * d_o * L_t * (1 - n_f * t_f)
            #         = A_prime_piL * pi * L_t
            A_prime_piL = g.d_o * (1 - g.n_f * g.t_f)

            # A_tot = A_prime + A_fins
            #       = A_tot_piL * pi * L_t
            A_tot_piL = A_prime_piL + A_fins_piL

            # eta_w = (A_prime + eta_f * A_fins) / A_tot
            eta_w = (A_prime_piL + eta_f * A_fins_piL) / A_tot_piL

            # Air cooled condenser with finned tubing
            U_CON_inv = self.add_expression(
                'U_CON_inv',
                (1 / h_t_CON) * self.A_per_A_i
                + A_tot_piL * numpy.log(g.d_o / g.d_i) / (2 * k_t)
                + R_con * self.A_per_A_o
                + 1 / h_air / eta_w
            )
            self.add_expression('U_CON', 1 / U_CON_inv)

            U_DES_inv = self.add_expression(
                'U_DES_inv',
                (1 / h_t_DES) * self.A_per_A_i
                + A_tot_piL * numpy.log(g.d_o / g.d_i) / (2 * k_t)
                + R_con * self.A_per_A_o
                + 1 / h_air / eta_w
            )
            self.add_expression('U_DES', 1 / U_DES_inv)

            return U_CON_inv, U_DES_inv

        if evaporating:  # (twophase mixture in shell side)
            # Nucleate boiling coefficinet according to Mostinski
            p_r = c.p / c.p_crit
            F_p = 1.8 * p_r ** 0.17 + 4 * p_r ** 1.2 + 10 * p_r ** 10
            DelatT_e = Max(0.1, self.get_expression('T_w') - c.T_sat)
            h_nb = self.add_expression(
                'alpha_nb',
                (1.4692e-15 * c.p_crit ** 2.3
                 * DelatT_e ** (7/3) * F_p ** (10/3))
            )
            # Correction for convective and bundle effects
            C = 3 ** 0.5 / 2 if g.arr == 'triangular' else 1  # else square
            F_b = 1 + 0.1 * (0.785 * g.d_tb * g.d_o
                             / (C * g.p_t ** 2) - 1) ** 0.75
            h_s_calc = self.add_expression('alpha_s', h_nb * F_b + 250)
        else:  # (single phase fluid in shell)
            Re_s_calc = self.add_expression(
                'Re_s',
                g.d_e * c.mdot / g.A_scs * 1 / Max(1e-9, c.mu)
            )
            # self.add_ge_constraint(Re_s_calc, 2100, 'Re_s >= 2100')
            Re_s = Max(2100, Re_s_calc)

            # Colburn factor
            j_H = 0.5 * (1 + g.B / g.d_s) * (
                    0.08 * Re_s ** 0.6821 + 0.7 * Re_s ** 0.1772)

            # ratio for reasonable mu values won't be too large
            mu_frac = Max(0.1, Max(1e-9, c.mu) / Max(1e-9, c.mu_w))
            self.add_expression('mu_c', c.mu)
            self.add_expression('mu_c_w', c.mu_w)
            self.add_expression('mu_c_frac', c.mu/c.mu_w)

            # self.add_ge_constraint(c.Pr, 0.6, 'Pr_c >= 0.6')
            Pr_c = Max(0.6, self.add_expression('Pr_c', c.Pr))

            # Nusselt number according to the simplified Delaware method
            Nu_s = self.add_expression('Nu_s',
                                       j_H * Pr_c ** (1/3) * mu_frac ** 0.14)
            h_s_calc = self.add_expression('alpha_s', Nu_s * c.k / g.d_e)
        h_s = Max(100, h_s_calc)
        self.add_ge_constraint(h_s_calc, 100, 'h_s >= 100')

        # shell & tube exchanger
        U_inv = self.add_expression(
            'U_inv',
            g.d_o / g.d_i / h_t
            + g.d_o * numpy.log(g.d_o / g.d_i) / 2 / k_t
            + 1 / h_s + R
        )
        self.add_expression('U', 1 / U_inv)

        return U_inv

    def set_U(self, U=None, **dependencies):
        """Set the value for U and relate A to thermal quantities."""
        if U is None:
            U_inv = self.calc_U_inv(**dependencies)
        else:
            try:
                U = float(U)
                assert not dependencies, ('When a fixed U value is given, no '
                                          'dependencies should be specified')
            except ValueError:  # interpret U as a path to a file
                import json
                with open(U, 'r') as f:
                    ref, model = json.load(f)

                for dependency_name, dependency in dependencies.items():
                    self.add_expression(dependency_name, dependency)
                # parse string representation to expression
                U_expr = comando.S(model)
                # determine symbols acting as placeholders
                deps = U_expr.free_symbols
                # Replace placeholder symbols with corresponding expressions
                U = Max(10, U_expr.subs({dep: self.get_expression(dep.name)
                                         for dep in deps}))
                self.add_ge_constraint(U, 1, 'U >= 1 W/m2/K')
            self.add_expression('U', U)
            U_inv = 1 / U

        A = self.get_expression('A')
        if self.kind == 'ACC':
            U_CON_inv, U_DES_inv = U_inv
            # F = 1 for isothermal condensation
            A_CON = (self.get_expression('Qdot_CON')
                     * self.get_expression('rlmtd_CON') * U_CON_inv)
            A_DES = (self.get_expression('Qdot_DES')
                     * self.get_expression('rlmtd_DES')
                     * self.get_expression('F_inv_DES') * U_DES_inv)
            self.add_expression('L_DES', self.L_t * self.n_tp * A_DES / A)
            self.add_expression('L_CON', self.L_t * self.n_tp * A_CON / A)
            A_req = self.add_expression('A_req', A_CON + A_DES)
        else:
            A_req = self.add_expression('A_req',
                                        self.get_expression('Qdot') * U_inv
                                        * self.get_expression('rlmtd')
                                        * self.get_expression('F_inv'))

        self.add_eq_constraint(A, A_req, 'A == (Qdot/LMTD/F/U)')

    def set_investment(self, A, p=None, year=2021, kind='fixed_tubesheet'):
        """Calculate the investment costs.

        The corresponding cost correlations are taken from:

        ```bibtex
        @Book{turton2018analysis,
          author    = {Turton, Richard and
                       Shaeiwitz, Joseph A. and
                       Bhattacharyya, Debangsu and
                       Whiting, Wallace B.},
          date      = {2018-08-01},
          title     = {Analysis, Synthesis and Design of Chemical Processes},
          isbn      = {0134177401},
          publisher = {Prentice Hall},
          ean       = {9780134177403}
        }
        ```

        giving costs for September 2001 (CEPCI = 397) and are corrected to the
        given year via an extrapolation of the CEPCI.

        Arguments
        ---------
        A : the relevant heat exchange area [m2]
        p : maximum pressure [Pa]
        year : year for cost correction
        kind : str
            either 'fixed_tubesheet', 'kettle_reboiler', or 'ACC'
        """
        # Material factor Table A3 p. 1301 and
        #      carbon steel  stainless steel
        # FmA       1              2.75        FmA = 2.75
        FmA = 2.75

        # Base burchase cost
        # turton2018analysis Table A1 p. 1286
        #                  k1A    k2A     k3A    A_min A_max
        # Fixed tube:      4.3247 -0.303  0.1634 10    1000
        # kettle reboiler: 4.4646 -0.5277 0.3955 10    10
        if kind == 'fixed_tubesheet':
            Cp0 = cost_function(A, 4.3247, -0.303, 0.1634)
        elif kind == 'kettle_reboiler':
            Cp0 = cost_function(A, 4.4646, -0.5277, 0.3955)
        elif kind == 'U_tube':
            Cp0 = cost_function(A, 4.1884, -0.2503, 0.1974)
        elif kind == 'ACC':
            Cp0 = cost_function(A, 4.0336, 0.2341, 0.0497)
            FmA = 1

        # Pressure correction
        Fp = 1
        if p is not None:
            # turton2018analysis Table A2 p. 1298
            # Both fixed tube and kettle reboiler:
            # c1A(0.03881), c2A(-0.11272), c3A(0.08183)
            # NOTE: pressure correction is based on bar gauge (barg)
            #       (1e5 Pa = 1 bar = 0 barg)
            Fp = cost_function(p/1e5 - 1, 0.03881, -0.11272, 0.08183)

        # Bare module cost factors Table A4 p. 1303
        # B1A(1.63), B2A(1.66)
        FbmA = Cp0 * (1.63 + 1.66 * FmA * Fp)

        # Bare module cost corrected to current year
        C_BM = CEPCI_correction(CEPCI_HIST[2001], year=2021) * FbmA

        if kind == 'ACC':
            inv_con = self.add_expression('investment_condenser', C_BM)

            P_F_nom_min = 37e3
            P_F_nom_max = 200e3
            P_F_min_rel = 0.05
            P_F_min = P_F_nom_min * P_F_min_rel

            P_F_nom, P_F_nom_scaled = add_scaled_var(self, 'P_F_nom',
                                                     P_F_nom_min, P_F_nom_max,
                                                     design=True)

            P_F, P_F_scaled = add_scaled_var(self, 'P_F', 0, P_F_nom_max)
            # self.add_ge_constraint(P_F, P_F_min_rel * P_F_nom,
            #                        f'P_F >= {P_F_min_rel} * P_F_nom')
            self.add_ge_constraint(P_F_nom, P_F, 'P_F_nom >= P_F')

            P = self.add_expression('P', P_F * self.n_fans)
            P_req = self.get_expression('P_req')
            self.add_eq_constraint(P, P_req, 'P == P_req')

            Vdot_air_F = self.add_expression(
                'Vdot_air_F', self.get_expression('Vdot_air') / self.n_fans
            )
            self.add_expression(
                'vel_F', Vdot_air_F / self.get_expression('A_F')
            )
            # FbmF = cost_function(Vdot_air_nom, 3.1761, -0.1373, 0.3414)
            # inv_F_Turton = self.add_expression(
            #     'investment_fan_Turton', CEPCI_correction(CEPCI_HIST[2001], year=2021) * FbmF * self.n_bays
            # )
            # self.add_le_constraint(Vdot_air_F, 100, 'Vdot_air_F <= 100 m3/s')

            CEPCI_01_2000 = 391.1
            inv_F_Smith = self.add_expression(
                'investment_fan_Smith', CEPCI_extrapolation(2021) / CEPCI_01_2000
                * 1.23e4 * (P_F_nom / 50e3) ** 0.76 * self.n_fans
            )
            inv_F = self.add_expression('investment_fan', inv_F_Smith)
            self.add_expression('investment', inv_con + inv_F)
        else:
            self.add_expression('investment', C_BM)


class Pump(comando.Component):
    """A simple model for an electrical pump.

    Arguments
    ---------
    label : `str`
        A unique sting that serves as an identifier
    mdot: expression
        The mass flow passing through the pump
    p_in : expression
        The pressure at the input of the pump
    p_out : expression
        The pressure at the output of the pump
    h_in : expression
        The specific enthalpy at the input of the pump
    h_out_is : expression
        The specific enthalpy at the output of the pump for isentropic
        compression
    eta_is : expression
        The isentropic efficiency of the pump
    """

    def __init__(self, label, mdot, p_in, p_out, h_in, h_out_is, eta_is,
                 eta_m=0.95):
        super().__init__(label)
        exprs = self.expressions_dict
        exprs.update({'mdot': mdot, 'p_in': p_in, 'p_out': p_out,
                      'h_in': h_in, 'h_out_is': h_out_is})
        exprs['w'] = w = (exprs['h_out_is'] - exprs['h_in']) / eta_is
        # h_out = h_in + w
        #       = h_in + (h_out_is - h_in) / eta_is
        #       = h_out_is / eta_is - (1/eta_is - 1) * h_in
        exprs['h_out'] = exprs['h_out_is'] / eta_is \
            - (1 / eta_is - 1) * exprs['h_in']

        # exprs['P'] = P = exprs['mdot'] * w / eta_m
        # self.add_ge_constraint(P, 0, 'P >= 0')

        P, P_scaled = add_scaled_var(self, 'P', 1e5, 2e6)
        self.add_eq_constraint(self['P'], exprs['mdot'] * w / eta_m,
                               'P == mdot * dh / eta_m')

        # Design variable
        P_nom, P_nom_scaled = add_scaled_var(self, 'P_nom', 1e5, 2e6,
                                             design=True)
        self.add_ge_constraint(P_nom, P, 'P_nom >= P')
        # silveira2003thermoeconomic
        # exprs['investment'] = 3540 * pow(P_nom * 1e-3, 0.71)
        # astolfi2014binaryb
        inv_EUR_2014 = 14000 * pow(P_nom / 200e3, 0.67)
        exprs['investment'] = COST_CONVERSION_ASTOLFI * inv_EUR_2014

        self.add_input('Hdot_in', exprs['mdot'] * exprs['h_in'])
        self.add_output('Hdot_out', exprs['mdot'] * exprs['h_out'])

        for id, expr in exprs.items():
            self.add_expression(id, expr)

    # def set_investment(self):
    #     P_nom = self.get_expression('P_nom')
    #     P = self.get_expression('P')
    #     index = comando.utility.get_index(P)
    #     P_max = self.add_expression('P_max',
    #                                 Max(*(comando.utility.parse(P, idx=i)
    #                                               for i in index)))
    #     self.add_eq_constraint(P_nom, P_max, 'P_nom = P_max')


class Turbine(comando.Component):
    """A simple model for a vapor turbine.

    Arguments
    ---------
    label : A unique sting that serves as an identifier
    mdot: The mass flow passing through the turbine
    p_in : The pressure at the input of the turbine
    p_out : The pressure at the output of the turbine
    h_in : The specific enthalpy at the input of the turbine
    h_out_is : The specific enthalpy at the output of the turbine for
        isentropic expansion
    eta_is : The isentropic efficiency of the turbine
    """

    def __init__(self, label, mdot, p_in, p_out, h_in, h_out_is, v_vap, eta_is,
                 eta_g=0.95, units=1, stages=1):
        super().__init__(label)
        exprs = self.expressions_dict
        self.eta_g = eta_g
        self.v_vap = v_vap
        self.units = self.add_expression('units', units)
        self.stages = self.add_expression('stages', stages)
        mdot_t = mdot / units
        exprs.update({'mdot': mdot_t, 'p_in': p_in, 'p_out': p_out,
                      'h_in': h_in, 'h_out_is': h_out_is})

        self.eta_is = eta_is
        dh_is = exprs['dh_is'] = Max(exprs['h_in'] - h_out_is, 10e3)
        self.add_ge_constraint(dh_is, 10e3, 'dh_is >= 10e3')
        exprs['dh'] = dh = dh_is * eta_is
        h_out = exprs['h_out'] = exprs['h_in'] - dh

        # exprs['P'] = P = exprs['mdot'] * dh * eta_g
        # self.add_ge_constraint(P, 1, 'P >= 1 MW')
        # self.add_le_constraint(P, 15e6, 'P <= 15 MW')

        P, P_scaled = add_scaled_var(self, 'P', 1e6, 15e6)
        self.add_eq_constraint(self['P'], exprs['mdot'] * dh * eta_g,
                               'P == mdot * dh * eta_g')
        # Design variable
        P_nom, P_nom_scaled = add_scaled_var(self, 'P_nom', 1e6, 15e6,
                                             design=True)
        self.add_ge_constraint(P_nom, P, 'P_nom >= P')

        v_in = v_vap(p_in, h_in)
        Vdot_in = exprs['Vdot_in'] = Max(mdot_t * v_in, 1)
        v_out = v_vap(p_out, h_out)
        Vdot_out = exprs['Vdot_out'] = Max(mdot_t * v_out, 1)
        v_out_is = v_vap(p_out, h_out_is)
        Vdot_out_is = exprs['Vdot_out_is'] = Max(mdot_t * v_out_is, 1)
        # self.add_ge_constraint(Vdot_out, 1, 'Vdot_out <= 1')

        # Current VR := Vdot_out_is / Vdot_in
        self.add_expression('current_VR', Vdot_out_is / Vdot_in)

        # Ratio between output and input volume flow
        self.add_expression('actual_volume_ratio', Vdot_out / Vdot_in)

        # Pressure ratio
        PR = self.add_expression('PR', Min(0.999, p_out / p_in))
        self.add_le_constraint(p_out, p_in, 'p_out <= p_in')

        # reduced mass flow
        phi = self.add_expression('phi', mdot_t / (p_in / v_in) ** 0.5)
        kappa = 1.08
        a = (2 / (kappa + 1)) ** (self.stages * kappa / (kappa - 1))

        K_Stodola, _ = add_scaled_var(self, 'K_Stodola', 0.01, 0.05, design=True)
        self.add_eq_constraint(K_Stodola,
                               mdot_t / (p_in / v_in * (1 - (Max(0, PR - a) / (1 - a)) ** 2)) ** 0.5,
                               "Stodola's ellipse law")

        # K_Stodola, _ = add_scaled_var(self, 'K_Stodola', 0.001, 5000, design=True)
        # g = 9.81
        # def f(PR):
        #     return PR ** (2 / kappa) - PR ** ((kappa + 1) / kappa)
        #
        # # With this correlation K_Stodola corresponds to nozzle area!
        # self.add_eq_constraint(K_Stodola,
        #                        mdot_t / (2 * g * kappa / (kappa - 1)
        #                                  * p_in / v_in \
        #                                  * comando.Max(f(PR), 0.001)) ** 0.5,
        #                        "Single isentropic nozzle")

        self.set_investment(P_nom, dh_is, Vdot_out_is)

        self.add_input('Hdot_in', exprs['mdot'] * exprs['h_in'])
        self.add_output('Hdot_out', exprs['mdot'] * exprs['h_out'])

        for id, expr in exprs.items():
            self.add_expression(id, expr)

    def set_investment(self, P_nom, dh_is_des=None, Vdot_out_is_des=None):
        """Set the investment costs.

        Note that the real flow rate can be used for a conservative estimate of
        costs, as it is larger than the flow rate under isentropic expansion.

        Arguments
        ---------
        P_nom : nominal electical power [W]
        dh_is_des : isentropic enthalpy drop at the design point [J/kg]
        Vdot_out_is_des : outlet flow rate under isentropic expansion [m3/s]
        n : number of stages
        """
        # Astolfi
        if dh_is_des is None:
            dh_is_des, _ = add_scaled_var(self, 'dh_is_des',
                                          10e3 * self.stages,
                                          65e3 * self.stages,
                                          design=True)

        if Vdot_out_is_des is None:
            Vdot_out_is_des, _ = add_scaled_var(self, 'Vdot_out_is_des',
                                                10, 100, design=True)

        n = self.stages
        dh_is_des_stage = dh_is_des / n
        self.add_expression('NS', Vdot_out_is_des ** 0.5 / dh_is_des_stage ** 0.75)
        SP = self.add_expression('SP',
                                 Vdot_out_is_des ** 0.5 / dh_is_des_stage ** 0.25)
        C0 = 1.23e6  # €
        n0 = 2  # -
        SP0 = 0.18  # m
        C_T = self.add_expression('investment_turbine',
                                  COST_CONVERSION_ASTOLFI * self.units
                                  * C0 * (n / n0) ** 0.5 * (SP / SP0) ** 1.1)

        C_G0 = 200e3  # €
        P_nom0 = 5e6  # W
        C_G = self.add_expression('investment_generator',
                                  COST_CONVERSION_ASTOLFI * self.units
                                  * C_G0 * (P_nom / P_nom0) ** 0.67)
        self.add_expression('investment', (C_T + C_G))

    def consider_off_design_Ghasemi(self, v_vap):
        mdot = self['mdot']
        p_out = self['p_out']
        dh_rel = self.make_operational_variable('dh_rel', bounds=(0.2, 1.2))
        Vdot_out_rel = self.make_operational_variable('Vdot_out_rel',
                                                      bounds=(0.2, 1.2))

        # Original correlation
        # NOTE: This has a bad relaxation due to the polynomial form, for a
        #       particular instance, the two fits below need about 600
        #       iterations, while no termination is reached after 15000
        #       iterations for the original form!
        # r = self.add_expression('r', r_Ghasemi(dh_rel, Vdot_out_rel))

        # single ANN for r
        # r = self.add_expression(
        #     'r',
        #     0.788244315940139
        #     - 0.125538123499497 * tanh(-5.05288676247691 * Vdot_out_rel
        #                                + 0.342038486516288 * dh_rel
        #                                + 2.56368525799966)
        #     + 0.0826184905352226 * tanh(0.0967477215040203 * Vdot_out_rel
        #                                 + 5.99878176443436 * dh_rel
        #                                 - 3.7011654313655)
        # )

        # # separate ANNs for rh and rVdot (for arguments from 0.2 to 1.2)
        rh = (0.213946609079737
              + 0.770558293837478 * tanh(0.0641545397993912
                                         + 1.71396121263698 * dh_rel)
              + 0.100292835001918 * tanh(2.82762125246941
                                         - 2.16281940159924 * dh_rel))
        rVdot = (0.704715907112635
                 - 0.275815866961846 * tanh(0.396299444020688
                                            - 3.79822439333504 * Vdot_out_rel)
                 - 0.0200169522504285 * tanh(3.17864205120355
                                             - 5.09840252639156 * Vdot_out_rel))

        self.add_ge_constraint(rh, 0.61, 'rh >= 0.61')
        self.add_ge_constraint(rVdot, 0.78, 'rVdot >= 0.78')
        r = self.add_expression('r', rh * rVdot)
        self.add_ge_constraint(r, 0.4758, 'r >= 0.4758')

        dh = self.add_expression('dh', self['dh_is'] * self.eta_is * r)

        P, P_scaled = add_scaled_var(self, 'P', 1e6, 15e6)
        self.add_eq_constraint(self['P'], mdot * dh * self.eta_g,
                               'P == mdot * dh * eta_g')

        h_out = self.add_expression('h_out', self['h_in'] - dh)

        n = self.stages
        dh_des, dh_des_scaled = add_scaled_var(self, 'dh_des',
                                               10e3 * n, 65e3 * n,
                                               design=True)

        self.add_eq_constraint(dh_des, dh / dh_rel, 'dh_des == dh / dh_rel')

        # Update Volume flow
        Vdot_out = self.add_expression('Vdot_out', mdot * v_vap(p_out, h_out))
        Vdot_out_des, _ = add_scaled_var(self, 'Vdot_out_des', 1, 50,
                                         design=True)
        self.add_eq_constraint(Vdot_out_des, Vdot_out / Vdot_out_rel,
                               'Vdot_out_des == Vdot_out / Vdot_out_rel')

        # Update existing constraint with new expression for P
        P_nom = self['P_nom']
        self.add_ge_constraint(P_nom, P, 'P_nom >= P')

        dh_is_des = self.add_expression('dh_is_des', dh_des / self.eta_is)

        # Update investment calculations
        # (Use Vdot_out_des as conservative estimate of Vdot_out_is_des)
        self.set_investment(P_nom, dh_is_des, Vdot_out_des)
        Vdot_in = self.get_expression('Vdot_in')
        Vdot_out_is = self.add_expression('Vdot_out_is',
                                          mdot * v_vap(p_out,
                                                       self['h_out_is']))

        # Current VR := Vdot_out_is / Vdot_in
        self.add_expression('current_VR', Vdot_out_is / Vdot_in)

        # Ratio between output and input volume flow
        self.add_expression('actual_volume_ratio', Vdot_out / Vdot_in)


def r_Ghasemi(dh_rel, Vdot_rel):
    """Calculate the off-design factor for turbine efficiciency.

    Source: http://dx.doi.org/10.1016/j.energy.2012.10.039

    Note that the correlation given in the above source has a typo, the correct
    correlation can be found in the associated preprint on Researchgate.

    Arguments
    ---------
    dh_rel : enthalpy drop divided by maximum possible enthalpy drop
    Vdot_rel : outlet flow rate divided by maximum outlet flow rate

    Returns
    -------
    r : off-design correction factor for isentropic turbine efficiency
    """
    rT = (dh_rel) ** 0.5
    rh = (((1.398 * rT - 5.425) * rT + 6.274) * rT - 1.866) * rT + 0.619
    rVT = Vdot_rel
    # Note that in the official version, published in Energy the last
    #      coefficient below is incorrectly given as 0.38, i.e., it is missing
    #      one additional zero!
    rv = (((-0.21 * rVT + 1.117) * rVT - 2.533) * rVT + 2.588) * rVT + 0.038
    return rh * rv


def vel_max_gas(p, M=58.12, F_mat=1):
    """Estimate the maximum allowable gas velocity in steel tubes.

    Arguments
    ---------
    p : pressure [Pa]
    M : molecular mass [g/mol]
    F_mat : material factor (1.5 for stainless steel, 5 for titanium)

    Returns
    -------
    vel_max : maximum allowable velocity [m/s]
    """
    p_psi = p / 6894.76  # convert to psi
    return F_mat * 548.64 / (M * p_psi) ** 0.5

