"""Case study for the optimal operating point if an organic rankine cycle.

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
# Copyright © 2020 Institute of Energy and Climate Research
# Energy Systems Engineering (IEK-10)
# Forschungszentrum Jülich GmbH
# Tel.: +49 2461 61-96307
# http://www.fz-juelich.de/iek/iek-10/EN
# 52425 Jülich, Germany
#
# AUTHORS: Marco Langiu
import pickle
import sys
from time import time

from numpy import array, atleast_1d, minimum, maximum
from pandas import Series

import CoolProp.CoolProp as CP
CP.set_reference_state('isobutane', 'NBP')

import comando
comando.set_backend('symengine')


def p_sat(T):
    """Calculate the saturation pressure of isobutane given temperature."""
    from comando import tanh
    return (17799254.1384071
            + 331236.710292479 * tanh(2.67857499525001
                                      - 0.0106734099609901 * T)
            - 18147149.5826235 * tanh(4.19335080960359
                                      - 0.00781580113529208 * T))


p_crit_br = 22064000.0  # Pa
p_crit_ib = 3629000.0  # Pa


def v_br(T_br):
    """Calculate the specific volume of brine, given temperature.

    This fit was created assuming properties of water at p = 897 kPa.
    """
    return (0.00104034947240847
            + 8.42239252864916e-7 * (-369.639795918367 + T_br))


def mu_br(T_br):
    """Calculate the kinematic viscosity of brine, given temperature.

    This fit was created assuming properties of water at p = 897 kPa.
    """
    return (0.000306228501649275
            - 2.60156065738765e-06 * (-365.476530612245 + T_br))


def Pr_br(T_br):
    """Calculate the Prandtl number of brine, given temperature.

    This fit was created assuming properties of water at p = 897 kPa.
    """
    return 1.90974767090314 - 0.0163723946967896 * (-365.476530612245 + T_br)


def k_br(T_br):
    """Calculate the thermal conductivity of brine, given temperature.

    This fit was created assuming properties of water at p = 897 kPa.
    """
    return 0.000261524442034341 * T_br + 0.578798774654223


def rho_air(T):
    """Calculate air density for a given temperature."""
    return 1.2045 - 0.004299 * (T - 293.15)


def v_air(T):
    """Calculate air specific volume for a given temperature."""
    return 0.8302 + 0.002841 * (T - 293.15)


def mu_air(T):
    """Calculate air visosity for a given temperature."""
    return 1.821e-05 + 4.903e-08 * (T - 293.15)


def Pr_air(T):
    """Calculate air Prandtl number for a given temperature."""
    return 0.708 + -0.0001391 * (T - 293.15)


def k_air(T):
    """Calculate air thermal conductivity for a given temperature."""
    return 0.02587 + 7.527e-05 * (T - 293.15)


dTmin_eva = 1

T_air_max = 363.15  # 90ºC
dT_air_in_pinch_max = 75
cp_air = 1006
P_AIR = 101325  # Pa (1 atm)

vel_fac = 1
vel_max_liq = 3.0 * vel_fac # maximum allowed speed for liquids
vel_max_gas = 20.0 * vel_fac  # maximum allowed speed for gases
vel_face_max = 3.5 * vel_fac  # maximum allowed velocity of air at ACC inlet

# MDOT_BR = 567.73
MDOT_BR = 660

P_BR = 897e3
cp_br = 4100  # J/kg/K
T_BR_IN = 408.15

# Variable bounds
T_br_out_min, T_br_out_max = 333.15, 403.15

dTmin_con_min = 1
dTmin_con_max = 55

A_min_con = 1e4
A_max_con = 1e6

mdot_lb, mdot_ub = 50, 1500

P1_LB = 1.1e5
P1_UB = 20e5
P2_LB = 3e5
P2_UB = 25e5

# operational cost as % of annualized initial investment (huster2019uimpact)
op_cost_fraction = 0.06
# Balance of plant cost as % of equipment cost (macchi2017organic)
BOP_fraction = 0.3

n = 20  # project lifetime
i = 0.06  # discount relative_change
annuity_factor = i / (1 - 1 / (1 + i) ** n)
var_cost = 4  # US$ / MWh
tot_operation = 8000  # h/a
C_el = 80  # US$ / MWh

h_liq_hp_min = CP.PropsSI('H', 'P', P2_LB, 'Q', 0, 'isobutane')
h_liq_hp_max = CP.PropsSI('H', 'P', P2_UB, 'Q', 0, 'isobutane')
dh_eva_p1_min = CP.PropsSI('H', 'Q', 1, 'P', P1_UB, 'isobutane') \
    - CP.PropsSI('H', 'Q', 0, 'P', P1_UB, 'isobutane')


def run_ORC_case_study(silent=True, tol=1e-3, maxTime=None, T_c=288,
                       use_superheater=False, use_recuperator=False,
                       obj='TAC',
                       subdir='.',
                       etat=0.88,  # da Lio (VR ≈ 4, SP ≈ 0.2)
                       off_design=None, branching_priorities=0, plot=False,
                       start_from=None, fix_design=False):
    """Run the organic Rankine cycle case study."""
    from .ORC_components import HeatExchanger, Pump, Turbine
    from .ANN import ann_from_json
    from pathlib import Path
    ann_dir = Path(__file__).parent.resolve() / Path('old_anns')

    if use_superheater:
        model_directory = 'use_superheater'
    if use_recuperator:
        model_directory = 'use_recuperator'
    if not use_superheater and not use_recuperator:
        model_directory = 'heat_exchanger'

    t1 = time()
    ES = comando.System('ORC')

    # Normalize T_c, to either stay a float or convert to a Series if multiple
    # values are given or it was explicitly requested by setting
    # off_design=True
    if isinstance(T_c, (int, float)):
        scenarios = Series(1, [T_c])
    else:
        off_design = True
        if isinstance(T_c, list):
            scenarios = Series(1/len(T_c), T_c)
        elif isinstance(T_c, dict):
            scenarios = Series(T_c)
            T_c = list(T_c)

    print('Considered scenarios and their weights:\n\n')
    print(scenarios)
    print('\noptions:')
    print(f'{off_design = }, {use_superheater = }, {use_recuperator = }, '
          f'{tol = }')
    # Variables (scaled to bounds in range 0, 1)

    def add_scaled_var(name, lb, ub, init=None, design=False):
        """Add a new variable scaled to the range [0, 1] to the model.

        Returns
        -------
        unscaled_var : Expression
            the expression corresponding to the unscaled (original) variable
        scaled_var: Variable or VariableVector
            the scaled variable with bounds 0, 1
        """
        if init is not None:
            init = (init - lb) / (ub - lb)
        Var = ES.make_design_variable if design \
            else ES.make_operational_variable
        scaled_var = Var(name + '_scaled', bounds=(0, 1), init_val=init)
        unscaled_var = lb + scaled_var * (ub - lb)
        return ES.add_expression(name, unscaled_var), scaled_var

    mdot_br = ES.add_expression('mdot_br', MDOT_BR)
    # mdot_br, mdot_br_scaled = add_scaled_var('mdot_br', MDOT_BR/10, MDOT_BR)

    p_air = ES.add_expression('p_air', P_AIR)

    # NOTE: Limiting lower pressure bound for individual scenarios to
    #       T_c + DeltaT_min
    T_sat_lb = array(T_c) + dTmin_con_min
    p1_lbs = atleast_1d(CP.PropsSI('P', 'Q', 1, 'T', T_sat_lb,
                                   'isobutane'))
    T_sat_ub = minimum(T_air_max,
                       array(T_c) + dTmin_con_max + dT_air_in_pinch_max)
    p1_ubs = minimum(P1_UB, atleast_1d(CP.PropsSI('P', 'Q', 1, 'T', T_sat_ub,
                                                  'isobutane')))
    p1_lbs = maximum(P1_LB, p1_lbs)

    p1_range = P1_UB - P1_LB
    p1_lbs_scaled = (p1_lbs - P1_LB) / p1_range
    # p1_ubs_scaled = (p1_ubs - p1_lb) / p1_range
    p1_ubs_scaled = 1

    p2_lbs = maximum(P2_LB, p1_lbs + 2e5)
    p2_range = P2_UB - P2_LB
    p2_lbs_scaled = (p2_lbs - P2_LB) / p2_range

    p1, p1_scaled = add_scaled_var('p_1', P1_LB, P1_UB)
    ES.add_expression('p1_lb', P1_LB)
    ES.add_expression('p1_ub', P1_UB)

    p2, p2_scaled = add_scaled_var('p_2', P2_LB, P2_UB)
    p2_max, p2_max_scaled = add_scaled_var('p_2_max', P2_LB, P2_UB,
                                           design=True)
    ES.add_ge_constraint(p2_max_scaled, p2_scaled, 'p2_max >= p2')
    ES.add_expression('p2_lb', P2_LB)
    ES.add_expression('p2_ub', P2_UB)

    # mdot_br, _ = add_scaled_var('mdot_br', MDOT_BR * 0.2, MDOT_BR)
    mdot_cp_br = mdot_br * cp_br
    ES.add_expression('p_br', P_BR)
    T_br_in = ES.make_parameter('T_br_in', T_BR_IN)  # K
    T_br_out, T_br_out_scaled = add_scaled_var('T_br_out',
                                               T_br_out_min, T_br_out_max)  # K

    dTmin_con, dTmin_con_scaled = add_scaled_var('dTmin_con',
                                                 dTmin_con_min, dTmin_con_max)

    eta_p = ES.make_parameter('eta_p', 0.8)
    etat = ES.make_parameter('etat', etat)

    Qdot_in = ES.add_expression('Qdot_in', mdot_cp_br * (T_br_in - T_br_out))

    s_sat_liq = ann_from_json(ann_dir / 's_sat_liq_p_6_6.json')
    s_vap = ann_from_json(ann_dir / 's_vap_p_h_6_6.json')

    T_liq = ann_from_json(ann_dir / 'T_liq_p_h_6_6.json')
    T_sat = ann_from_json(ann_dir / 'T_sat_p_6_6.json')
    T_vap = ann_from_json(ann_dir / 'T_vap_p_h_6_6.json')

    h_liq = ann_from_json(ann_dir / 'h_liq_p_s_6_6.json')
    h_sat_liq = ann_from_json(ann_dir / 'h_sat_liq_p_6_6.json')
    h_sat_vap = ann_from_json(ann_dir / 'h_sat_vap_p_6_6.json')
    h_vap = ann_from_json(ann_dir / 'h_vap_p_s_6_6.json')

    Pr_liq = ann_from_json(ann_dir / 'Pr_liq_p_h_6_6.json')
    Pr_vap = ann_from_json(ann_dir / 'Pr_vap_p_h_6_6.json')
    k_liq = ann_from_json(ann_dir / 'k_liq_p_h_6_6.json')
    k_vap = ann_from_json(ann_dir / 'k_vap_p_h_6_6.json')

    # Specific volume (1/density)
    v_liq = ann_from_json(ann_dir / 'v_liq_p_h_6_6.json')

    # NOTE: We're introducing maxima here to assure the solver that no division
    #       by zero can occur
    def v_vap(p, h):
        return comando.Max(0.0023,
                           ann_from_json(ann_dir / 'v_vap_p_h_6_6.json')(p, h))

    # NOTE: only ANNs for the inverse of the viscosity were trained
    muinv_liq = ann_from_json(ann_dir / 'muinv_liq_p_h_6_6.json')

    def mu_liq(p, h):
        return 1/comando.Max(3407.6, muinv_liq(p, h))

    def mu_vap(p, h):
       return 1/comando.Max(14737.0,
                            ann_from_json(ann_dir
                                          / 'muinv_vap_p_h_6_6.json')(p, h))

    s1 = s_sat_liq(p1)
    T1 = T_sat(p1)
    h1 = h_sat_liq(p1)
    h7 = h_sat_vap(p1)
    h4 = h_sat_vap(p2)
    h2s = h_liq(p2, s1)
    h_hot_pinch = h_sat_liq(p2)
    T4 = T_sat(p2)

    mdot, mdot_scaled = add_scaled_var('mdot', mdot_lb, mdot_ub)  # kg/s

    pump = Pump('pump', mdot, p1, p2, h1, h2s, eta_p)

    if use_recuperator:
        h2 = pump['h_out']  # INFO: With recuperator
        h2r, _ = add_scaled_var('h2r', h_liq_hp_min, h_liq_hp_max)
        ES.add_le_constraint(h2, h2r, 'h2 <= h2r')
    else:
        h2 = h2r = pump['h_out']  # INFO: Without recuperator

    T2 = T_liq(p2, h2)
    T2r = T_liq(p2, h2r)

    if use_superheater:
        h5 = Qdot_in / mdot + h2r
        T5 = T_vap(p2, h5)
    else:
        T5 = T4
        h5 = h4
        ES.add_eq_constraint(h2r, h5 - Qdot_in / mdot,
                             'h2r == h5 - Qdot_in / mdot')
        ES.add_expression('Qdot_in_ib', mdot * (h5 - h2r))
        ES.add_expression('Qdot_in_br', Qdot_in)

    s5 = s_vap(p2, h5)

    # In off-design the output of the economizer is generally not saturated
    if off_design:
        h3, _ = add_scaled_var('h_ib_eco_out', h_liq_hp_min, h_liq_hp_max)
        ES.add_ge_constraint(h3, h_hot_pinch * 0.95,
                             'h3 >= h_hot_pinch * 0.95')
        ES.add_le_constraint(h3, h_hot_pinch, 'h3 <= h_hot_pinch')
        T3 = T_liq(p2, h3)
        ES.add_le_constraint(T3, T4, 'T3 <= T4')
    else:
        h3 = h_hot_pinch
        T3 = T_sat(p2)

    if use_recuperator:
        ES.add_le_constraint(h2r, h3, 'h2r <= h3')

    ES.add_expression('dTmin_eva', dTmin_eva)
    economizer = HeatExchanger('eco', dTmin_eva,
                               T_h_out=T_br_out, T_c_in=T2r, T_c_out=T3,
                               h_c_in=h2r, h_c_out=h3, mdot_c=mdot,
                               mdot_cp_h=mdot_cp_br)
    T_br_3 = economizer['T_h_in']

    T_br_pinch = T_br_out + mdot * (h_hot_pinch - h2r) / mdot_cp_br
    A_eco = economizer.set_shell_and_tube_geometry(p_t=0.02064, d_o=0.01588,
                                                   d_i=0.01423,
                                                   )

    T_br_eco = economizer.get_expression('T_h')
    T_w_eco = economizer.get_expression('T_w')
    # NOTE: fluid properties are normally evaluated at averages for p and T
    #       however, we only have ANNs ready for 1/mu(p,h), so we instead use
    #       'mean enthalpies' corresponding to the mean temperatures.
    #       For the wall viscosity, we use the maximum enthalpy in lack of one
    #       corresponding to T_w
    # T_ib_eco = economizer.get_expression('T_c')
    h_ib_eco = 0.5 * (h2r + h3)
    h_w_eco = h3
    economizer.set_U(
        mdot_h=mdot_br, v_h=v_br(T_br_eco), mu_h=mu_br(T_br_eco),
        Pr_h=Pr_br(T_br_eco), k_h=k_br(T_br_eco),
        mdot_c=mdot, mu_c=mu_liq(p2, h_ib_eco), mu_w_c=mu_liq(p2, h_w_eco),
        Pr_c=Pr_liq(p2, h_ib_eco), k_c=k_liq(p2, h_ib_eco))

    economizer.set_investment(A_eco, p2_max)

    # Limit economizer area such that brine velocity <= vel_max_liq m/s
    economizer.limit_tube_velocity(mdot_br * v_br(T_br_3), vel_max_liq)
    Vdot_liq_2 = mdot * v_liq(p2, h2)
    economizer.limit_shell_velocity(Vdot_liq_2, vel_max_liq)

    if use_superheater:
        superheater = HeatExchanger('sup', dTmin_eva, T_h_in=T_br_in,
                                    T_c_in=T4, T_c_out=T5,
                                    h_c_in=h4, h_c_out=h5,
                                    mdot_c=mdot, mdot_cp_h=mdot_cp_br)
        T_br_4 = superheater['T_h_out']

        A_sup = superheater.set_shell_and_tube_geometry(p_t=0.02064,
                                                        d_o=0.01588,
                                                        d_i=0.01423,
                                                        )

        T_br_sup = superheater.get_expression('T_h')
        T_w_sup = superheater.get_expression('T_w')
        h_ib_sup = 0.5 * (h4 + h5)
        h_w_sup = h5
        superheater.set_U(
            mdot_h=mdot_br, v_h=v_br(T_br_sup), mu_h=mu_br(T_br_sup),
            Pr_h=Pr_br(T_br_sup), k_h=k_br(T_br_sup),
            mdot_c=mdot, mu_c=mu_vap(p2, h_ib_sup), mu_w_c=mu_vap(p2, h_w_sup),
            Pr_c=Pr_vap(p2, h_ib_sup), k_c=k_vap(p2, h_ib_sup))

        superheater.set_investment(A_sup, p2_max)

        # Limit superheater area such that brine velocity <= vel_max_liq m/s
        superheater.limit_tube_velocity(mdot_br * v_br(T_br_in), vel_max_liq)
        Vdot_vap_5 = mdot * v_vap(p2, h5)
        superheater.limit_shell_velocity(Vdot_vap_5, vel_max_gas)
    else:
        T_br_4 = T_br_in

    # NOTE: We set T_c_in = T4 instead of T3 here, because LMTD calculations
    #       for evaporators should be based on saturation temperature only!
    evaporator = HeatExchanger('eva', dTmin_eva, T_h_in=T_br_4,
                               T_h_out=T_br_3, T_c_in=T4, T_c_out=T4,
                               h_c_in=h3, h_c_out=h4, mdot_c=mdot,
                               mdot_cp_h=mdot_cp_br)
    # NOTE: By fully specifying both sides of the evaporator, we introduce the
    #       additional constraint that the heat flows computed based on the
    #       corresponding quantities at both sides must be equal!
    #       This constraint is redundant, as it is an equivalent form of
    #       h5 = h2r + Qdot_in / mdot.
    #       Nevertheless, removing it leads to slowdown, (17.976 s, 3235 it)
    #       even when adding it as relaxation-only (16.941 s 3237 it)!
    eva_heat_flow_con = evaporator._constraints_dict.pop('heat_flow')

    A_eva = evaporator.set_shell_and_tube_geometry(p_t=0.02064, d_o=0.01588,
                                                   d_i=0.01423, n_tp=2)

    # NOTE: As a conservative estimate, we calculate lmtd and U based on T_sat
    T_br_eva = 0.5 * (T_br_pinch + T_br_4)
    T_w_eva = 0.5 * (T_br_eva + T4)
    evaporator.set_U(
        mdot_h=mdot_br, v_h=v_br(T_br_eva), mu_h=mu_br(T_br_eva),
        Pr_h=Pr_br(T_br_eva), k_h=k_br(T_br_eva),
        p_c=p2, p_crit_c=p_crit_ib, T_sat_c=T4, evaporating=True)

    evaporator.set_investment(A_eva, p2_max, kind='U_tube')

    # Limit evaporator area such that brine velocity <= vel_max_liq m/s
    evaporator.limit_tube_velocity(mdot_br * v_br(T_br_4), vel_max_liq)

    # NOTE: We can either express h6s as h_vap(p1, s5), or we introduce a tear
    #       variable to reset relaxations
    h6s = h_vap(p1, s5)
    # h6s, h6s_scaled = add_scaled_var('h_6s', 400e3, 550e3, 475e3)

    T6s = T_vap(p1, h6s)
    # s6s = s_vap(p1, h6s)
    # Ensure tear enthalpies are consistent
    # ES.add_eq_constraint((h6s - h_vap(p1, s5)) * 1e-5, 0, 'h6s == h_vap(p1, s5)')
    # Ensure entropies based on tear variable is consistent with calculated
    # ES.add_eq_constraint((s5 - s6s), 0, 's5 == s6s')

    ###########################################################################
    # NOTE: macchi2017axial states that maximum electric power output per
    #       turbine is typically limited to 15 MWel. Also the real plant
    #       investigated by ghasemi2013modeling had two turbines.
    ###########################################################################
    turbine = Turbine('tur', mdot, p2, p1, h5, h6s, v_vap, etat, units=2)
    turbine.consider_off_design_Ghasemi(v_vap)
    h6 = turbine['h_out']
    T6 = T_vap(p1, h6)

    if use_recuperator:
        dTmin_rec = ES.make_parameter('dTmin_rec', 10)  # K
        T6r = T_vap(p1, h6 - (h2r - h2))
        mdot_rec = ES.add_expression('mdot_rec', mdot/turbine.units)
        recuperator = HeatExchanger('rec', dTmin_rec, T_h_in=T6,
                                    T_h_out=T6r, T_c_in=T2, T_c_out=T2r,
                                    h_h_in=h6, h_c_in=h2, h_c_out=h2r,
                                    mdot_h=mdot_rec, mdot_c=mdot_rec)
        h6r = recuperator['h_h_out']  # == h6 - (h2r - h2)

        A_rec = recuperator.set_shell_and_tube_geometry(p_t=0.03969,
                                                        d_o=0.03175,
                                                        d_i=0.02964,
                                                        )

        h_ib_h_rec = 0.5 * (h6 + h6r)
        h_ib_c_rec = 0.5 * (h2 + h2r)
        h_w_rec = 0.5 * (h_ib_h_rec + h_ib_c_rec)
        recuperator.set_U(
            mdot_h=mdot_rec, v_h=v_vap(p1, h_ib_h_rec),
            mu_h=mu_vap(p1, h_ib_h_rec),
            Pr_h=Pr_vap(p1, h_ib_h_rec), k_h=k_vap(p1, h_ib_h_rec),
            mdot_c=mdot_rec, mu_c=mu_liq(p2, h_ib_c_rec),
            mu_w_c=mu_liq(p2, h_w_rec), Pr_c=Pr_liq(p2, h_ib_c_rec),
            k_c=k_liq(p2, h_ib_c_rec)
        )

        recuperator.set_investment(A_rec, p2_max)

        # Multiply investment by number of recuperators
        recuperator.add_expression('investment',
                                   recuperator.get_expression('investment')
                                   * turbine.units)

        # Limit recuperator area such that isobutane vel <= vel_max_gas m/s
        recuperator.limit_tube_velocity(turbine.get_expression('Vdot_out'),
                                        vel_max_gas)

        Vdot_liq_2r = mdot * v_liq(p2, h2r)
        recuperator.limit_shell_velocity(Vdot_liq_2r, vel_max_liq)
    else:
        T6r, h6r = T6, h6

    ###########################################################################
    # Direct air cooling
    T_air_in = ES.add_expression('T_air_in', ES.make_parameter('T_air_in'))
    T7 = T1  # condensation

    # air pinch
    T_air_pinch = T1 - dTmin_con
    dT_air_p_calc = T_air_pinch - T_air_in

    dh_eva_p1 = h7 - h1
    dh_eva_p1 = comando.Max(dh_eva_p1_min, dh_eva_p1)
    ES.add_ge_constraint(dh_eva_p1, dh_eva_p1_min,
                         'dh_eva_p1 >= dh_eva_p1_min')
    Q_17 = mdot * dh_eva_p1
    mdot_cp_air = Q_17 / comando.Max(dTmin_con_min, dT_air_p_calc)
    ACC = HeatExchanger('ACC', dTmin_con, T_h_in=T6r,
                        T_h_out=T1, T_c_in=T_air_in, h_h_in=h6r,
                        h_h_out=h1, mdot_h=mdot, mdot_cp_c=mdot_cp_air)
    T_air_out = ES.add_expression('T_air_out', ACC['T_c_out'])

    ES.add_le_constraint(T_air_out, T_air_max, f'T_air_out <= {T_air_max} K')

    mdot_air = ES.add_expression('mdot_air', mdot_cp_air / cp_air)

    ACC.set_ACC_geometry(A_min_con, A_max_con, T_air_pinch, T7, Qdot_con=Q_17)

    # NOTE: As a conservative estimate we again calculate LMTD and U based
    #       on the conditions of the saturated liquid
    T_air = 0.5 * (T_air_in + T_air_out)
    h_DES = 0.5 * (h7 + h6r)
    ACC.set_U(
        mdot_h=mdot,
        v_h=v_liq(p1, h1),
        mu_CON_h=mu_liq(p1, h1),
        mu_DES_h=mu_vap(p1, h_DES),
        Pr_CON_h=Pr_liq(p1, h1),
        Pr_DES_h=Pr_vap(p1, h_DES),
        k_CON_h=k_liq(p1, h1),
        k_DES_h=k_vap(p1, h_DES),
        p_h=p1, p_crit_h=p_crit_ib,
        mdot_c=mdot_air, v_in_c=v_air(T_air_in), v_c=v_air(T_air), rho_c=rho_air(T_air),
        mu_c=mu_air(T_air), Pr_c=Pr_air(T_air), k_c=k_air(T_air),
        vel_face_max_c=vel_face_max, condensing=True
    )

    Vdot_gas_max = mdot * v_vap(p1, h6r)
    ACC.limit_tube_velocity(Vdot_gas_max, vel_max_gas, '_gas')
    Vdot_liq_max = mdot * v_liq(p1, h1)
    ACC.limit_tube_velocity(Vdot_liq_max, vel_max_liq, '_liq')

    ACC.set_investment(ACC.A_o, kind='ACC')

    components = [
        pump,
        economizer,
        evaporator,
        turbine,
        ACC
    ]
    if use_superheater:
        components.append(superheater)
    if use_recuperator:
        components.append(recuperator)
    for comp in components:
        ES.add(comp)

    P_net = ES.add_expression('P_net', turbine.units * turbine['P']
                              - pump['P'] - ACC['P'])

    T_cold = T_air_in + dTmin_con_min
    T_hot = T_br_in - dTmin_eva
    eta_carnot = ES.add_expression('eta_carnot', 1 - T_cold/T_hot)
    ES.add_le_constraint(P_net, eta_carnot * Qdot_in, 'P_net <= P_carnot')

    # NOTE: Scaling this constraint appears to be detrimental!
    ES.add_ge_constraint(P_net, 0, 'P_net >= 0')
    ES.add_le_constraint((T1 - T6s), 0, 'T1 <= T6s')

    ES.add_expression('p1', p1)
    ES.add_expression('p2', p2)

    ES.add_expression('T1', T1)
    ES.add_expression('T2', T2)
    ES.add_expression('T2r', T2r)
    ES.add_expression('T3', T3)
    ES.add_expression('T4', T4)
    ES.add_expression('T5', T5)
    ES.add_expression('T6', T6)
    ES.add_expression('T6r', T6r)
    ES.add_expression('T7', T7)

    ES.add_expression('s1', s1)
    ES.add_expression('s5', s5)

    ES.add_expression('h1', h1)
    ES.add_expression('h2', h2)
    ES.add_expression('h2r', h2r)
    ES.add_expression('h3', h3)
    ES.add_expression('h_hp', h_hot_pinch)
    ES.add_expression('h4', h4)
    ES.add_expression('h5', h5)
    ES.add_expression('h6', h6)
    ES.add_expression('h6r', h6r)
    ES.add_expression('h7', h7)

    ES.add_expression('efficiency', P_net / Qdot_in)
    try:
        ES.add_expression('T_cw_in', T_cw_in)
        ES.add_expression('T_cw_p', T_cw_p)
        ES.add_expression('T_cw_out', T_cw_out)
    except NameError:
        pass
    ES.add_expression('T_air_in', T_air_in)
    ES.add_expression('T_air_cp', T_air_pinch)
    ES.add_expression('T_air_out', T_air_out)

    ES.add_expression('T_br_in', T_br_in)
    ES.add_expression('T_br_4', T_br_4)
    ES.add_expression('T_br_hp', T_br_pinch)
    ES.add_expression('T_br_3', T_br_3)
    ES.add_expression('T_br_out', T_br_out)

    # Fixed costs for exploration, drilling, etc. in $
    exploration_and_development_costs = ES.add_expression('E&D_cost', 15e6)
    equipment_cost = ES.add_expression(
        'equipment_cost',
        ES.aggregate_component_expressions('investment')
    )
    ES.add_expression('BOP_cost', BOP_fraction * equipment_cost)
    inv = ES.add_expression('investment', exploration_and_development_costs
                            + (1 + BOP_fraction) * equipment_cost)

    # Maximize net power production, considering the nominal operating point
    T_c_rep = '_'.join(f'{T:g}' for T in scenarios.index)
    if obj == 'P_net':
        P = ES.create_problem(0, -P_net * 1e-6, scenarios=scenarios,
                              name=f'ORC_off-design_{T_c_rep}')

        ES.add_expression(
            'P_net_mean', P.weighted_sum(P_net)
        )
    elif obj == 'LCOE':
        P_net_mean, P_net_mean_scaled = add_scaled_var('P_net_mean',
                                                       0.1e6, 50e6,
                                                       design=True)

        P_net_mean_MW = comando.Max(1, P_net_mean * 1e-6)
        LCOE = ES.add_expression('LCOE',
                                 (1 + op_cost_fraction) * annuity_factor * inv
                                 / (P_net_mean_MW * tot_operation) + var_cost)

        P = ES.create_problem(LCOE, 0, scenarios=scenarios,
                              name=f'ORC_off-design_{T_c_rep}')

        P.constraints['P_net_mean = mean(P_net)'] = \
            comando.Eq(P_net_mean, P.weighted_sum(P_net))
        P.add_symbols({P_net_mean_scaled})
    elif obj == 'TAC':
        CAPEX = ES.add_expression(
            'CAPEX', (1 + op_cost_fraction) * annuity_factor * inv
        )
        OPEX = ES.add_expression(
            'OPEX', -C_el * (P_net * 1e-6 * tot_operation)
        )
        P = ES.create_problem(CAPEX, OPEX, scenarios=scenarios,
                              name=f'ORC_off-design_{T_c_rep}')
        ES.add_expression(
            'TAC', P.objective
        )

        ES.add_expression(
            'P_net_mean', P.weighted_sum(P_net)
        )

    else:
        raise ValueError("obj needs to be one of 'P_net' or 'LCOE'!")

    P['ORC_T_air_in'] = T_c

    # Adding auxiliary constraints forcing P_nom to be the maximum value of P
    # for comp in pump, turbine, ACC:
    #     index = comando.utility.get_index(comp['P'])
    #     expr = comp['P']
    #     sym_maps = [{sym: sym[idx] for sym in expr.free_symbols if sym.indexed}
    #                 for idx in index]
    #     P.constraints[f'{comp.label}_P_nom == {comp.label}_P_max'] \
    #         = comando.Eq(comp['P_nom'],
    #                      comando.Max(*(expr.subs(reps) for reps in sym_maps)))

    # Update bounds
    p1_scaled.lb = p1_lbs_scaled
    p1_scaled.ub = p1_ubs_scaled
    p2_scaled.lb = p2_lbs_scaled

    # Manual bound adjustments for optimizations with reduced ranges
    # dv_bnds =  {
    ###########################################################################
    # normal pump
    ###########################################################################
    #     'eco_d_s_scaled': (0.2614032374806455, 0.48107557959943426),
    #     'eva_d_s_scaled': (0.34446300152438314, 0.40310754186843545),
    #     'sup_d_s_scaled': (0.3606909464087589, 0.6924217162238628),
    #     'rec_d_s_scaled': (0.33821083249926187, 1),
    #     'eco_L_per_d_s_scaled': (0.9999999999804658, 0.9999999999804658),
    #     'eva_L_per_d_s_scaled': (0.7761308185610004, 1),
    #     'sup_L_per_d_s_scaled': (0, 0.2060942118735339),
    #     'rec_L_per_d_s_scaled': (0, 0.04187234438507546),
    #     'eco_B_per_d_s_scaled': (0.3503875485709095, 0.9085720890879553),
    #     'sup_B_per_d_s_scaled': (0.9999999999931103, 0.9999999999931103),
    #     'rec_B_per_d_s_scaled': (0, 0.672810810277909),
    #     'ACC_A_scaled': (0.9054310909404975, 1),
    #     'ORC_p_2_max_scaled': (0.4231340951167861, 0.7549456246734089),
    #     'ACC_P_F_nom_scaled': (1.5182702362391052e-12, 1.5182702362391052e-12),
    #     'pump_P_nom_scaled': (0.9999999999998577, 0.9999999999998577),
    #     'tur_P_nom_scaled': (0.34880572005153226, 1),
    #     'tur_Vdot_out_des_scaled': (0.14693962779671452, 0.8005967438984161),
    #     'tur_dh_des_scaled': (0.5347588183267318, 1),
    #     'tur_K_Stodola_scaled': (0.1739912565497725, 0.6706411708995028)
    ###########################################################################
    # large pump
    ###########################################################################
    # 'eco_d_s_scaled': (0.43487561321569135, 0.5242417787647259),
    # 'eva_d_s_scaled': (0.3463639027337039, 0.41587634715989563),
    # 'sup_d_s_scaled': (0.3101767055317472, 0.6708124363522635),
    # 'rec_d_s_scaled': (0.3288586971664608, 1),
    # 'eco_L_per_d_s_scaled': (0.9999999999913637, 0.9999999999913637),
    # 'eva_L_per_d_s_scaled': (0.9756610457941579, 1),
    # 'sup_L_per_d_s_scaled': (4.642471019458101e-13, 4.642471019458101e-13),
    # 'rec_L_per_d_s_scaled': (1.7630206193201963e-12, 1.7630206193201963e-12),
    # 'eco_B_per_d_s_scaled': (0.2582357497165164, 0.3744696225126911),
    # 'sup_B_per_d_s_scaled': (0.9999999999990369, 0.9999999999990369),
    # 'rec_B_per_d_s_scaled': (0, 0.6863372993278605),
    # 'ACC_A_scaled': (0.904483627495233, 1),
    # 'ORC_p_2_max_scaled': (0.4394675736145944, 0.8389670190162987),
    # 'ACC_P_F_nom_scaled': (1.0321943260465258e-12, 1.0321943260465258e-12),
    # 'pump_P_nom_scaled': (0.4764643278129342, 0.6300048439026376),
    # 'tur_P_nom_scaled': (0.36789262957806884, 1),
    # 'tur_Vdot_out_des_scaled': (0.14116745537453995, 0.8005967438982184),
    # 'tur_dh_des_scaled': (0.5706260367141801, 1),
    # 'tur_K_Stodola_scaled': (0.12080290555455185, 0.632137395012611)
    # }
    # for dv in P.design_variables:
    #     try:
    #         lb, ub = dv_bnds[dv.name]
    #     except KeyError:
    #         print('No bound update for', dv)
    #     dv.lb = max(dv.lb, lb)
    #     dv.ub = min(dv.ub, ub)
    # ov_bnds = {
    ###########################################################################
    # normal pump
    ###########################################################################
    #     'ORC_mdot_scaled': (0.17640703331282157, 0.258040331550651),
    #     'ORC_p_1_scaled': (0.04078034326702827, 0.3348365187061884),
    #     'ORC_p_2_scaled': (0.4231340951018073, 0.7549456246441316),
    #     'ORC_h2r_scaled': (0.24518879250817438, 0.5914866427738358),
    #     'ORC_h_ib_eco_out_scaled': (0.5559109392254025, 0.8299895428253424),
    #     'ACC_P_F_scaled': (0, 0.19499999999884482),
    #     'pump_P_scaled': (0.7202400209747096, 1),
    #     'tur_P_scaled': (0.25151691180335684, 1),
    #     'tur_Vdot_out_rel': (0.5551259440096945, 1.2),
    #     'tur_dh_rel': (0.45699661201382147, 1.2),
    #     'ORC_T_br_out_scaled': (0.07529781078109858, 0.4270159649991394),
    #     'ORC_dTmin_con_scaled': (0.043325194190139624, 0.0728825210523157)
    ###########################################################################
    # large pump
    ###########################################################################
    # 'ORC_mdot_scaled': (0.1752227863456328, 0.25472986529325065),
    # 'ORC_p_1_scaled': (0.04107747669865423, 0.331929466399127),
    # 'ORC_p_2_scaled': (0.43946757359926747, 0.8389670189894092),
    # 'ORC_h2r_scaled': (0.20572441769508173, 0.5228053910441782),
    # 'ORC_h_ib_eco_out_scaled': (0.5707661297093068, 0.8909875680496007),
    # 'ACC_P_F_scaled': (0, 0.19499999997409367),
    # 'pump_P_scaled': (0.47346201416228745, 0.6300048436120712),
    # 'tur_P_scaled': (0.3059981836178343, 1),
    # 'tur_Vdot_out_rel': (0.49975250762139845, 1.2),
    # 'tur_dh_rel': (0.5008867260287917, 1.2),
    # 'ORC_T_br_out_scaled': (0.08476858099656023, 0.44559918372944035),
    # 'ORC_dTmin_con_scaled': (0.03676397944661167, 0.07232451711336958)
    # }
    # for ov in P.operational_variables:
    #     try:
    #         lb, ub = ov_bnds[ov.name]
    #     except KeyError:
    #         print('No bound update for', ov)
    #     for ovi in ov:
    #         ov.lb = max(ovi.lb, lb)
    #         ov.ub = min(ovi.ub, ub)

    print()
    print('Bounds for p1_scaled:\n')
    print(p1_scaled.lb)
    print(p1_scaled.ub)
    print()

    print('Lower bounds for p2_scaled:\n')
    print(p2_scaled.lb)
    print()

    print(f'\nCreating ORC model took: {time() - t1} s')

    print(f'Problem has {P.num_vars} variables and {P.num_cons} constraints.')

    ###########################################################################
    # initialization
    ###########################################################################
    if start_from is not None:
        with open(start_from, 'rb') as f:
            results = pickle.load(f)
        for v in P.design_variables.union(P.operational_variables):
            try:
                v.value = results[v.name]
                print(f"found initial value of {v.value} for '{v}'")
            except KeyError:
                print(f"found no initial value for '{v}'")
                pass
    exceptions = {
        'ORC_P_net_mean_scaled',
    }
    if fix_design:
        for v in list(P.design_variables):
            if v.name in exceptions:
                print('not fixing', v.name)
            else:
                v.fix()
                print('fixed', v.name, 'to', v.value)

    ###########################################################################
    # Priority strategies
    ###########################################################################
    from comando.utility import get_vars, indexed
    from collections import Counter

    # NOTE: Helper methods to descend into a given level of the expression tree
    #       and to obtain the depth of a given symbol within an expression
    def descend(expr, depth):
        """Return a generator of subexpressions at the given depth.

        Note: Even works with negative depths!
        """
        if depth == 0:
            yield expr
        elif depth == 1:
            try:
                yield from expr.args
            except AttributeError:
                yield from []
        else:
            for arg in expr.args:
                yield from descend(arg, depth - 1)

    def symdepth(expr, sym, min_depth=0):
        """Obtain the maximum depth of `sym' within `expr'."""
        # print(expr)
        if min_depth > 0:
            for arg in expr.args:
                min_depth = max((min_depth,
                                 *(min_depth + symdepth(subexpr, sym)
                                   for subexpr
                                   in descend(arg, min_depth - 1))))
            return min_depth
        if expr.args:
            if sym in expr.free_symbols:
                for arg in expr.args:
                    # min_depth = 1 + symdepth2(arg, sym, min_depth - 1)
                    min_depth = max(min_depth,
                                    1 + symdepth(arg, sym, min_depth - 1))
                return min_depth
            else:
                return -1
        return (expr == sym) - 1

    # TESTS
    # symdepth(x, y) == 0
    # symdepth(x, x) == 1
    # symdepth(1 + x, x) == 2
    # symdepth(2 * (1 + x), x) == 3
    # symdepth((2 + 2 * (x)) * (1 + 3 * x ** 2), x)
    # symdepth(expr, x)
    # symdepth(expr, y)

    # [(ORC_p_1_scaled, 56), (ORC_mdot_scaled, 31), (ORC_h_6s_scaled, 30), (ORC_p_2_scaled, 44)]

    if branching_priorities == 0:
        print('Using presence-based branching priorities!')
        pp = Counter()
        for expr in (P.design_objective + P.operational_objective,
                     *P.constraints.values()):
            pp.update({var: len(scenarios) if indexed(expr) else 1
                       for var in get_vars(expr)})
        # default bp = 0
        # prios = {var: int(1 + pp_val) for var, pp_val in pp.items()}

        branching_priorities = '0-squared'
        prios = {var: int(1 + pp_val ** 2) for var, pp_val in pp.items()}

    elif branching_priorities == 1:
        print('Using default branching priorities!')
        prios = {var: 1
                 for var in P.design_variables.union(P.operational_variables)}
    elif branching_priorities == 2:
        print('Using branching priorities based on maximum depth!')
        # NOTE: Currently better than default but worse than presence-based
        #       Also current implementation is extremely slow, should analyze
        #       subexpressions from CSE instead of original expressions.
        from collections import defaultdict
        maxdepth = defaultdict(int)
        i = 1
        for expr in (P.design_objective + P.operational_objective,
                     *P.constraints.values()):
            print(i)
            i += 1
            for var in get_vars(expr):
                var_depth = maxdepth[var]
                maxdepth[var] = symdepth(expr, var, var_depth)
        prios = {var: depth for var, depth in maxdepth.items()}
    elif branching_priorities % 100 == 0:
        des_prio = branching_priorities // 100
        print(f'Using elevated branching priorities of {des_prio} for design variables!')
        prios = {var: des_prio for var in P.design_variables}
    else:
        raise NotImplementedError(f'Strategy {branching_priorities} for '
                                  'branching priority selection has not been '
                                  'implemented!')

    # NOTE: Here you can add variables to be eliminated completely from
    #       branching. In general this will put convergence at risk.
    auxiliaries = [
        p2_max_scaled,  # will be largest value of p2
        turbine['P_nom_scaled'],  # will be largest value of tur_P
        pump['P_nom_scaled'],  # will be largest value of pump_P
        ACC['P_F_nom_scaled'],  # will be largest value of ACC_P_F
        turbine['P_scaled'],  # implied by tur_P == mdot * dh * eta_T * eta_g
        pump['P_scaled'],  # implied by pump_P == mdot * dh / eta_P / eta_m
        ACC['P_F_scaled'],  # will be equal to P_req
        turbine['K_Stodola_scaled'],  # implied by tur_Stodola's ellipse law
        # P_net_mean_scaled
    ]
    if off_design:
        auxiliaries.extend([
            # economizer['d_s_scaled'],
            # economizer['L_per_d_s_scaled'],
            # evaporator['d_s_scaled'],
            ACC['A_scaled'],  # implied by ACC_A == (Qdot/LMTD/F/U)
            # turbine['Vdot_out_des_scaled'],  # either des or rel value implied
            turbine['Vdot_out_rel'],  # either des or rel value implied
            # turbine['dh_des_scaled'],  # either des or rel value implied
            turbine['dh_rel'],  # either des or rel value implied
        ])
        # for name in ['P_scaled', 'dh_is_rel', 'Vdot_out_rel']:
        #     try:
        #         expr = turbine[name]
        #         if isinstance(expr, comando.Symbol):
        #             auxiliaries.append(expr)
        #     except KeyError:
        #         pass

        # if use_superheater:
        #     auxiliaries.append(superheater['d_s_scaled'])
        #     # auxiliaries.append(superheater['L_per_d_s_scaled'])
        #     # auxiliaries.append(superheater['B_per_d_s_scaled'])
        # if use_recuperator:
        #     auxiliaries.append(recuperator['d_s_scaled'])

    var_prios = dict(sorted({**prios,
                             **{var: 1 for var in auxiliaries}}.items()))

    # Ensure design variables have higher priority than operational ones
    max_op_prio = max(var_prios.get(var, 1) for var in P.operational_variables)
    for var in P.design_variables:
        if var not in auxiliaries:
            var_prios[var] += max_op_prio

    print('Variables and branching priorities:')
    print('\n  Design:')
    from operator import attrgetter
    namegetter = attrgetter('name')

    for var in sorted(P.design_variables, key=namegetter):
        print(f'    {var}: {var_prios.get(var, 1)}')
    print('\n  Operation:')
    for var in sorted(P.operational_variables, key=namegetter):
        print(f'    {var}: {var_prios.get(var, 1)}')

    import json
    with open(subdir + '/bp.json', 'w') as f:
        json.dump({v.name: p for v, p in var_prios.items()}, f)

    print('Constraints\n\t' + '\n\t'.join(P.constraints.keys()))

    ro_cons = {
        # 'T2 - T1 >= 0': (T2 - T1) >= 0,
        # 'T5 - T6 >= 0': (T5 - T6) >= 0,
        # 'h2 - h1 >= 0': (h2 - h1) >= 0,
        # 'h5 - h6 >= 0': (T5 - T6) >= 0,
        # 'eva_heat_flow': eva_heat_flow_con,
        # 'ACC_P >= 0': P.constraints.pop('ACC_P >= 0'),
    }

    # Using CSE to obtain the intermediate expressions
    obj_proxy = P.design_objective + P.operational_objective
    intermediates, new_exprs = comando.cse((*ES.expressions_dict.values(),
                                            *P.constraints.values(),
                                            obj_proxy))

    defs = {}
    for sym, expr in intermediates:
        e = expr.subs(defs)
        index = comando.utility.get_index(e)
        try:
            b = comando.utility.bounds(e)
        except ValueError as e:
            print(e)
            if index is None:
                b = float('-inf'), float('inf')
            else:
                b = [float('-inf')], [float('inf')]

        if index is None:
            x = comando.Variable(sym.name, bounds=b)
        else:
            b = (min(b[0]), max(b[1]))
            x = comando.VariableVector(sym.name, bounds=b)
            x.instantiate(index)
        defs[sym] = x

    expr_map = {}
    for expr_id, maybe_sym in zip(ES.expressions_dict, new_exprs):
        if maybe_sym.is_Symbol:
            expr_map[maybe_sym] = expr_id

    with open(subdir + '/expressions.txt', 'w') as f, open(subdir + '/expressions.csv', 'w') as csv:
        for sym, expr in intermediates:
            lb, ub = defs[sym]._bounds
            gap = ub - lb
            data = [expr_map.get(sym, ''), str(sym),
                    *map(str, (lb, ub, gap)), str(expr)]
            csv.write('; '.join(data) + '\n')
        f.write('\n# INTERMEDIATE EXPRESSIONS\n')
        f.writelines(f'{sym} in {list(defs[sym]._bounds)} := {expr}\n'
                     for sym, expr in intermediates)
        exprs = iter(new_exprs)
        f.write('\n# EXPRESSIONS\n')
        f.writelines(f'"{e_id}" {e}\n'
                     for e_id, e in zip(ES.expressions_dict, exprs))
        f.write('\n# CONSTRAINTS\n')
        f.writelines(f'"{e_id}" {e}\n'
                     for e_id, e in zip(P.constraints, exprs))
        f.write('\n# OBJECTIVE\n')
        f.writelines(f'{next(exprs)}\n')

    # create a graphical representation of the problem
    if plot is True:
        print('Preparing plot...')

        from comando.linearization import _is_linear

        def is_linear(expr, var):
            try:
                return _is_linear(expr, [var])
            except RuntimeError:
                return False

        # Option 1: reuse above cse
        # con_exprs = new_exprs[len(ES.expressions_dict):-1]

        # Option 2: only consider constraints and objective
        intermediates, new_exprs = comando.cse((*P.constraints.values(),
                                                obj_proxy))
        con_exprs = new_exprs[:-1]

        import networkx as nx

        level_0 = [
            *P.design_variables,
            *P.operational_variables,
            # *P.parameters  # comment if interested only in variables
        ]
        ignore = {*P.parameters}
        levels = [level_0]

        seen = set()

        i = 0
        g = nx.DiGraph()
        g.add_nodes_from(level_0, level=i)

        intermediates = {k: v for k, v in intermediates}
        interms = {*intermediates}
        fzjcolors = {'blue': (2/255, 61/255, 107/255),
                     'lightblue': (173/255, 189/255, 227/255),
                     'green': (185/255, 210/255, 95/255),
                     'yellow': (250/255, 235/255, 90/255),
                     'violet': (175/255, 130/255, 185/255),
                     'red': (235/255, 95/255, 115/255),
                     'orange': (250/255, 180/255, 90/255)}
        dv_color = fzjcolors['orange']
        ov_color = fzjcolors['red']
        p_color = fzjcolors['green']
        f_color = fzjcolors['violet']
        con_color = fzjcolors['lightblue']
        obj_color = fzjcolors['blue']
        colors = {}
        node_colors = []
        labels = {}
        replacements = {}
        styles = []
        x_i = 0
        p_i = 0
        for sym in level_0:
            if sym.is_Parameter:
                labels[sym] = f'p{p_i}'
                colors[sym] = p_color
                node_colors.append(p_color)
                p_i += 1
            else:
                labels[sym] = f'x{x_i}'
                x_i += 1
                if sym.indexed:
                    colors[sym] = ov_color
                    node_colors.append(ov_color)
                else:
                    colors[sym] = dv_color
                    node_colors.append(dv_color)

        f_i = 1

        # First determine all ignored intermediate expressions (static)
        size = len(interms) + 1
        while len(interms) < size:  # As long as size changes...
            size = len(interms)
            # ... find nodes that only have dependencies on ignored quantities
            for xi in interms:
                syms = intermediates[xi].free_symbols
                if all(sym in ignore for sym in syms):
                    ignore.add(xi)
            interms.difference_update(ignore)

        while interms:
            seen.update(levels[i])
            i += 1
            next_level = {}
            # find nodes that only have dependencies in previous layers
            for xi in interms:
                syms = intermediates[xi].free_symbols
                if all(sym in seen or sym in ignore for sym in syms):
                    next_level[xi] = syms - ignore
                    name = f'f{f_i}'
                    labels[xi] = name
                    replacements[xi] = comando.Symbol(name)
                    colors[xi] = f_color
                    node_colors.append(f_color)
                    f_i += 1
            interms.difference_update(next_level)

            levels.append(next_level)
            g.add_nodes_from(next_level.keys(), level=i)
            for node, dependencies in next_level.items():
                g.add_edges_from(([(node, dependency)
                                   for dependency in dependencies]))
                subexpression = intermediates[node]
                styles.extend(
                    [':' if is_linear(subexpression, dependency) else '-'
                     for dependency in dependencies]
                )

        c_i = 0

        i += 1
        next_level = {}

        # add node corresponding to objective
        next_level[new_exprs[-1]] = new_exprs[-1].free_symbols - ignore
        labels[new_exprs[-1]] = 'obj'
        colors[new_exprs[-1]] = obj_color
        node_colors.append(obj_color)

        # add nodes corresponding to constraints
        for con_expr in con_exprs:
            next_level[con_expr] = con_expr.free_symbols - ignore
            labels[con_expr] = f'c{c_i}'
            colors[con_expr] = con_color
            node_colors.append(con_color)
            c_i += 1

        levels.append(next_level)
        g.add_nodes_from(next_level.keys(), level=i)
        for (node, dependencies), subexpression \
                in zip(next_level.items(), new_exprs):
            g.add_edges_from(([(node, dependency)
                               for dependency in dependencies]))
            styles.extend(
                [':' if is_linear(subexpression, dependency) else '-'
                 for dependency in dependencies]
            )

        # color = [subset_color[data["layer"]] for v, data in G.nodes(data=True)]
        from matplotlib import pyplot as plt

        width = max(len(level) for level in levels)

        Lx = i + 1
        Ly = width
        plt.figure(figsize=(Lx, Ly))

        pos = nx.multipartite_layout(g, subset_key="level")

        dy = 2
        dx = Lx / Ly * dy * (width - 1) / i
        pos = {n: [dx * (l - i / 2),
                   dy * ((len(level) - 1) / 2 - k)]
               for l, level in enumerate(levels)
               for k, n in enumerate(level)}

        import matplotlib.pyplot as plt

        nx.draw(g, pos, node_color=node_colors)
        ax = plt.gca()
        for patch, style in zip(ax.patches, styles):
            patch.set_linestyle(style)
        nx.draw_networkx_labels(g, pos, labels, font_size=8)
        plt.axis("equal")
        # plt.show()
        plt.savefig('subexpressions.pdf')

        # extract relevant subexpressions with intermediate symbols renamed as
        # in figures
        f = {}
        for xi, subexpression in intermediates.items():
            if xi in labels:
                f[labels[xi]] = subexpression.subs(replacements)
        # print in order
        for fi in sorted(f, key=lambda fi: int(fi[1:])):
            print(fi, f[fi])

        exit()

    def get_results():
        """Efficient evaluation of expressions using common subexpressions."""
        results = dict()
        syms = {*P.design_variables, *P.operational_variables, *P.parameters}
        for sym in syms:
            results[sym.name] = sym.value

        # Get a CSE representation based on all defined expressions
        intermediates, new_exprs = comando.cse(ES.expressions_dict.values())

        # evaluate all intermediate expressions that occur more than once
        values = {}
        for sym, expr in intermediates:
            par = comando.Parameter(sym.name)
            par.value = evaluate(expr.subs(values))
            values[sym] = par

        # use representation of original expressions in terms of intermediates
        for e_id, e in zip(ES.expressions_dict, new_exprs):
            results[e_id] = evaluate(e.subs(values))
        results['scenarios'] = scenarios
        return results

    ###########################################################################
    # Constraint scaling
    ###########################################################################
    from comando.utility import get_index, _parse, _idx_parse, numpy_op_map

    def neval(expr):
        """Auxiliary method to help evaluation of user-defined functions."""
        try:
            return float(expr)
        except RuntimeError:
            return float(type(expr)(*map(neval, expr.args)))

    def evaluate(expr, idx=None):
        """Evaluate the given expression at the symbols' current values."""
        index = get_index(expr)

        if index is None:  # scalar expression
            return neval(_parse(expr, {sym: sym.value
                                       for sym in expr.free_symbols},
                                numpy_op_map, float))
        # indexed expression
        if idx is None:  # No particular index desired
            sym_map = {sym: sym.value.values if sym.indexed else sym.value
                       for sym in expr.free_symbols}
            res = [neval(expr.subs({sym: val[i] if sym.indexed else val
                                    for sym, val in sym_map.items()}))
                   for i in range(len(index))]
            return Series(res, index, float)
        # only values for index=idx desired
        return neval(expr.subs({sym:
                                sym[idx].value if sym.indexed else sym.value
                                for sym in expr.free_symbols}))

    # def scale_cons_by_value(cons):
    #     for c_id, c in cons.items():
    #         try:
    #             expr = c.lts - c.gts
    #         except AttributeError:
    #             expr = c.lhs - c.rhs
    #         current_val = abs(evaluate(expr, 288))
    #         if 1e-4 < current_val:
    #             print(f'scaling constraint\n\t"{c_id}"\n by {current_val}\n')
    #             scaled_con = expr/current_val
    #             cons[c_id] = type(c)(scaled_con, 0)
    #
    # scale_cons_by_value(P.constraints)
    # scale_cons_by_value(ro_cons)

    # Example of manual scaling by name
    # cons = P.constraints
    # for c_id, c in cons.items():
    #     try:
    #         expr = c.lts - c.gts
    #     except AttributeError:
    #         expr = c.lhs - c.rhs
    #     if 'P' in c_id:
    #         expr *= 1e-6
    #     elif 'dh' in c_id:
    #         expr *= 1e-4
    #     elif 'A' in c_id:
    #         expr *= 1e-3
    #     elif 'dT' in c_id:
    #         expr *= 0.1
    #     elif 'T' in c_id:
    #         expr *= 0.01
    #     cons[c_id] = type(c)(expr, 0)

    for dv in P.design_variables:
        print('\t'.join(map(str, (dv, dv.lb, dv.value, dv.ub))))

    # P.constraints = {}  # deleting all constraints
    ###########################################################################
    # Solve with MAiNGO
    ###########################################################################
    print('Solving via MAiNGO...', end='' if silent else '\n', flush=True)
    maingo_t0 = time()
    maingo_options = dict(
        branching_priorities={**prios,
                              **{var: 0 for var in auxiliaries}},
        epsilonR=tol,
        epsilonA=tol,
        BAB_nodeSelection=0,  # use node with lowest lower bound in the tree
        BAB_branchVariable=1,   # use dimension with largest relative diameter
        LBP_solver=2,  # CPLEX
        PRE_maxLocalSearches=10,
    )
    if maxTime:
        maingo_options['maxTime'] = maxTime

    from comando.interfaces.maingo_api import MaingoProblem, FEASIBLE_POINT, GLOBALLY_OPTIMAL
    mp = MaingoProblem(P, ro_cons)

    # from comando.interfaces.maingo_api import MAiNGO, LANG_ALE
    # NOTE: The MAiNGO API just writes the problem to disc without using CSE,
    #       as a consequence the file size is enourmous and writing it takes
    #       forever! Instead, we can use our own implementation, that extracts
    #       all reoccurring subexpressions as redintermediate expressions.
    # solver = MAiNGO(mp)
    # solver.write_model_to_file_in_other_language(LANG_ALE,
    #                                              f'{P.name}_{T_cw_rep}.ale')

    # from comando.interfaces.maingo_ale import write_ale_file, write_settings_file
    # write_settings_file(maingo_options, f'{P.name}_settings.txt')
    # write_ale_file(P, f'{P.name}.ale', ro_cons)
    # exit(0)

    def fn(suffix):
        return subdir + '/' + P.name + suffix
    sol, res = mp.solve(**maingo_options, writeCsv=True,
                        iterations_csv_file_name=
                        fn(f'_iterations_{branching_priorities}.csv'),
                        log_file_name=fn('.log'),
                        result_file_name=fn('_results.txt'),
                        solution_and_statistics_csv_file_name=fn('_stat.csv'),
                        )

    try:
        i = [FEASIBLE_POINT, GLOBALLY_OPTIMAL].index(res)
    except ValueError:
        print('No feasible point found!')
        exit(0)

    print('Terminated with ' + ' '.join(res.name.split('_')).lower()
          + ('!' if i == 0 else ' solution!'))

    print(f' took {time() - maingo_t0} s')
    results = get_results()

    ###########################################################################
    # Solve with SCIP
    ###########################################################################
    # NOTE: SCIP appears to be inefficient for the present formulation
    # scip_t0 = time()
    # from comando.interfaces.scip import ScipProblem, scip_op_map, exp
    # scip_op_map['Max'] = comando.utility.smooth_max
    # scip_op_map['tanh'] = lambda arg: 1 - 2 / (exp(2 * arg) + 1)
    # sp = ScipProblem(P)
    # sp.solve(limits_gap=tol, limits_absgap=tol, numerics_epsilon=1e-10,
    #          numerics_feastol=1e-10, numerics_dualfeastol=1e-10)
    # print(f' took {time() - scip_t0} s')
    # results = get_results()

    ###########################################################################
    # Solve with BARON
    ###########################################################################
    # NOTE: BARON appears to have problems with the representation of the ANNs
    #       and reports that no feasible point can be found
    # baron_t0 = time()
    # from comando.interfaces.baron import solve, baron_str_map
    # Teach baron how to handle various functions we use
    # Note that the tanh function is handled as
    #     lambda arg: f'(1 - 2/(exp(2 * ({arg})) + 1))'}
    # by our interface
    # from functools import reduce
    #
    # smoother = 1e-08
    #
    # def baron_binary_smooth_min(a, b):
    #     """Obtain a string representation for a smooth min in BARON syntax."""
    #     return f'(0.5 * ({a} + {b} - (({a} - {b} + {smoother})^2)^0.5))'
    #
    # def baron_binary_smooth_max(a, b):
    #     """Obtain a string representation for a smooth max in BARON syntax."""
    #     return f'(0.5 * ({a} + {b} + (({a} - {b} + {smoother})^2)^0.5))'
    #
    # def baron_cost_func(x, _, c1, c2, c3):
    #     """Obtain a string representation for cost_function in BARON syntax."""
    #     log10x = f'log({x})/log(10)'
    #     return f'(10 ^ ({c1} + {c2} * {log10x} + {c3} * {log10x} ^ 2))'
    #
    # def baron_lmtd(dT1, dT2):
    #     """Obtain a string representation for lmtd in BARON syntax."""
    #     return f'(({dT1} - {dT2}) / log({dT1} / {dT2}))'
    #
    # def baron_rlmtd(dT1, dT2):
    #     """Obtain a string representation for rlmtd in BARON syntax."""
    #     return f'(log({dT1} / {dT2}) / ({dT1} - {dT2}))'
    #
    # baron_str_map['Min'] = lambda *args: reduce(baron_binary_smooth_min, args)
    # baron_str_map['Max'] = lambda *args: reduce(baron_binary_smooth_max, args)
    # baron_str_map['cost_function'] = baron_cost_func
    # baron_str_map['lmtd'] = baron_lmtd
    # baron_str_map['rlmtd'] = baron_rlmtd
    #
    # tol = 1e-3
    # baron_opts = dict(
    #     epsr=tol,
    #     epsa=tol,
    #     branching_priorities=prios,
    #     dolocal=0,
    #     numloc=0
    # )
    #
    # # DEBUG: eliminate all constraints
    # backup = {**P.constraints}
    # P.constraints = {}
    # solve(P, 'ORC_off-design_cse.bar', cse=True, reuse=False,
    #       **baron_opts)
    #
    # print(f' took {time() - baron_t0} s')
    # results = get_results()

    ###########################################################################
    # Print constraints value
    ###########################################################################
    def get_active():
        from comando.utility import get_index
        active = set()
        for c_id, c in P.constraints.items():
            expr = c.rhs - c.lhs
            index = get_index(expr)
            if index is None:
                val = evaluate(expr)
                print(c_id, val)
                if abs(val) <= 1e-3:
                    active.add(c_id)
            else:
                for i in index:
                    val = evaluate(expr, i)
                    print(c_id, i, val)
                    if abs(val) <= 1e-3:
                        active.add(c_id)
        return active

    # active = get_active()

    with open(f'{subdir}/results_{T_c_rep}.pickle', 'wb') as f:
        pickle.dump(results, f)

    return results


if __name__ == '__main__':
    import argparse

    def parse_dict(args):
        res = {}
        for arg in args:
            pos = arg.find(':')
            if pos == -1:
                raise ValueError
            res[float(arg[:pos])] = float(arg[pos+1:])
        return res

        return dict(map(lambda arg: (float(arg), float(arg)), args))

    def parse_list(args):
        return list(map(float, args))

    class T_Parser(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            try:
                result = parse_dict(values)
            except ValueError:
                try:
                    result = parse_list(values)
                    if len(result) == 1:
                        result = result[0]
                except ValueError:
                    raise ValueError(f'Could not understand argument {values} '
                                     'for option -T')
            setattr(namespace, self.dest, result)

    ap = argparse.ArgumentParser()
    ap.add_argument("-silent", "-s", required=False, default=False,
                    action='store_true', help="silencing output")
    ap.add_argument("-tol", "-t", required=False, default=1e-2,
                    type=float, help="absolute/relative optimality tolerance")
    ap.add_argument("-maxTime", "-time", required=False, default=None,
                    type=float, help="maximum time")
    ap.add_argument('-T_c', '-T', nargs='+', required=False, default=288,
                    action=T_Parser, help="Temperature of cooling water")
    ap.add_argument("-off-design", "-od", required=False, default=None,
                    action='store_true', help="whether to consider off-design"
                    "for single-temperature optimizations")
    ap.add_argument("-branching_priorities", "-bp", required=False, default=0,
                    type=int, help="select strategy for branching priorities")

    ap.add_argument("-use_superheater", "-sup", required=False, default=False,
                    action='store_true', help="whether to use a superheater")
    ap.add_argument("-use_recuperator", "-rec", required=False, default=False,
                    action='store_true', help="whether to use a recuperator")

    ap.add_argument("-plot", "-p", required=False, default=False,
                    action='store_true', help="whether to create a plot of "
                    "the subexpression graph of the resulting problem")

    ap.add_argument("-subdir", required=False, default='.',
                    type=str, help="subdirectory in which to store results")
    ap.add_argument("-obj", required=False, default='TAC',
                    type=str, help="choice for objective function")

    ap.add_argument("-start_from", "-start", required=False, default=None,
                    type=str, help="pickle file with start values")
    ap.add_argument("-fix_design", "-fix", required=False, default=False,
                    action='store_true', help="whether to fix design")
    if __package__ is None:
        from pathlib import Path
        file = Path(__file__).resolve()
        sys.path.append(str(file.parents[1]))

    import ORC_off_design
    __package__ = ORC_off_design.__name__
    results = run_ORC_case_study(**vars(ap.parse_args()))
    interfaces = ['MAiNGO-API']