"""Example Components for the COMANDO framework."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu, David Shu, Florian Joseph Baader
import comando
from comando.core import Component, BINARY


class Link(Component):
    """A directed link between two buses."""

    def __init__(self, label, lb=None, ub=None):
        """Initialize the Link with a label and transfer bounds."""
        super().__init__(label)

        bounds = [lb, ub]
        expr_bounds = [None, None]
        for i, bound in enumerate(bounds):
            if isinstance(bound, comando.Expression):
                expr_bounds[i], bounds[i] = bound, None
        transfer = self.make_operational_variable('transfer', bounds=bounds)

        expr_lb, expr_ub = expr_bounds
        if expr_lb is not None:
            self.add_le_constraint(expr_lb, transfer, 'transfer_lb')
        if expr_ub is not None:
            self.add_le_constraint(transfer, expr_ub, 'transfer_ub')

        self.add_input('IN', transfer)
        self.add_output('OUT', transfer)


class Grid(Component):
    """A component representing a generic grid connection."""

    def __init__(self, label, consumption_limit=None, feedin_limit=0, price=0,
                 compensation=0, single_connector=False, constrain_flow=False,
                 co2_factor=0):
        """Initialize the Grid.

        Arguments
        ---------
        - label : str
            Unique string that serves as an identifier of this Grid.
        - consumption_limit : None, numeric data or expression (default: None)
            Limit the amount of commodity that can be supplied by the grid.
        - feedin_limit : None, numeric data or expression (default: 0)
            Limit the amount of commodity that can be fed into the grid.
        - price : numeric data
            Price per unit of commodity that is consumed.
        - compensation : numeric data
            Price per unit of commodity that is fed into the grid.
        - constrain_flow: bool (default: False)
            Enforce a constraint on the input/output flows
        - single_connector: bool (default: False)
            Specify whether a single, bidirectional connector should be used in
            favor of two directed connectors (consumption/feed-in)
        """
        super().__init__(label)

        if consumption_limit == 0:  # We don't need to look at consumption
            consumption = price = spendings = 0
        else:
            if consumption_limit is None \
                    or isinstance(consumption_limit, (int, float)):
                consumption = self.make_operational_variable(
                    'consumption', bounds=(0, consumption_limit))
            else:  # consumption_limit is an expression
                consumption = self.make_operational_variable('consumption',
                                                             bounds=(0, None))
                self.add_le_constraint(consumption, consumption_limit,
                                       'consumption_limit')
            price = self.make_parameter('price', price)
            spendings = self.add_expression('spendings', price * consumption)

        if feedin_limit == 0:  # We don't need to look at feedin
            feedin = compensation = earnings = 0
        else:
            if feedin_limit is None \
                    or isinstance(feedin_limit, (int, float)):
                feedin = self.make_operational_variable('feedin',
                                                        bounds=(0,
                                                                feedin_limit))
            else:  # feedin_limit is an expression
                feedin = self.make_operational_variable('feedin',
                                                        bounds=(0, None))
                self.add_le_constraint(feedin, feedin_limit, 'feedin_limit')
            compensation = self.make_parameter('compensation', compensation)
            earnings = self.add_expression('earnings', compensation * feedin)

        self.add_expression('variable_costs', spendings - earnings)

        co2_factor = self.make_parameter('co2_factor', value=co2_factor)
        emissions = (consumption - feedin) * co2_factor
        self.add_expression('emissions', emissions)

        if constrain_flow and feedin != 0 and consumption != 0:
            if feedin_limit and consumption_limit:  # Neither is 0 or None
                # Big-M formulation
                consuming = self.make_operational_variable('consuming',
                                                           domain=BINARY,
                                                           init_val=1)
                self.add_le_constraint(consumption,
                                       consuming * consumption_limit,
                                       'consumption_limit')
                self.add_le_constraint(feedin,
                                       (1 - consuming) * feedin_limit,
                                       'feeding_in')
            else:
                # Complementarity formulation
                self.add_eq_constraint(feedin * consumption, 0,
                                       'flow_complementarity')
            if single_connector:
                self.add_connector('CONSUMPTION', feedin - consumption)
                return
        if feedin != 0:
            self.add_input('FEEDIN', feedin)
        if consumption != 0:
            self.add_output('CONSUMPTION', consumption)


class Consumer(Component):
    """A consumer of one or more commodities."""

    def __init__(self, label, *variable_commodities, **fixed_commodities):
        """Specify consumer demands for commodities.

        Arguments
        ---------
        - label : str
            Unique string that serves as an identifier of this Consumer.
        - commodities : dict mapping str to numeric data
            Demands for the different commodities.
        """
        super().__init__(label)
        from itertools import repeat
        commodities = {**dict(zip(variable_commodities, repeat(None))),
                       **fixed_commodities}
        for commodity, demand in commodities.items():
            if not isinstance(demand, comando.Expr):
                demand = self.make_parameter(f'{commodity}_demand', demand)
            self.add_input(commodity.upper(), demand)


class Source(Component):
    """A source for an arbitrary commodity."""

    def __init__(self, label, price=0, co2_factor=0):
        """Initialize the Source.

        Arguments
        ---------
        - label : str
            Unique string that serves as an identifier of this Source.
        - price : numeric data
            Price per unit of provided commodity, can be a scalar or a
            pandas Series.
        - co2_factor : numeric data
            CO2 emission per unit of provided commodity, can be a scalar or a
            pandas Series.
        """
        super().__init__(label)
        price = self.make_parameter('price', price)
        use = self.make_operational_variable('use', bounds=(0, None))
        self.add_expression('variable_costs', price * use)
        co2_factor = self.make_parameter('co2_factor', value=co2_factor)
        emissions = self['use'] * co2_factor
        self.add_expression('emissions', emissions)
        self.add_output('OUT', use)


class Sink(Component):
    """A sink for an arbitrary commodity."""

    def __init__(self, label, compensation=0, co2_factor=0):
        """Initialize the Sink.

        Arguments
        ---------
        - label : str
            Unique string that serves as an identifier of this Sink.
        - price : numeric data
            Compensation per unit of consumed commodity, can be a scalar
            or a pandas Series.
        """
        super().__init__(label)
        compensation = self.make_parameter('compensation', compensation)
        sink = self.make_operational_variable('sink', bounds=(0, None))
        self.add_expression('variable_costs', -compensation * sink)
        co2_factor = self.make_parameter('co2_factor', value=co2_factor)
        emissions = self['sink'] * co2_factor
        self.add_expression('emissions', -emissions)
        self.add_input('IN', sink)


class Demand(Component):
    """A demand requiring a known amount of an arbitrary commodity."""

    def __init__(self, label, data=0):
        """Initialize the Demand.

        Arguments
        ---------
        - label : str
            Unique string that serves as an identifier of this Demand.
        - data : numeric data
            Amount of consumed commodity, can be a scalar or a pandas
            Series.
        """
        super().__init__(label)
        demand = self.make_parameter('demand', data)
        self.add_input('IN', demand)


class Resource(Component):
    """A resource providing a known amount of an arbitrary commodity."""

    def __init__(self, label, data=0):
        """Initialize the Resource.

        Arguments
        ---------
        - label : str
            Unique string that serves as an identifier of this Resource.
        - data : numeric data
            Amount of provided commodity, can be a scalar or a pandas
            Series.
        """
        super().__init__(label)
        resource = self.make_parameter('resource', data)
        self.add_output('OUT', resource)


class LinearConverter(Component):
    """Notation from Chapter 3 of source: voll2013automated.

    Arguments
    =========
    label : str
        name of the component
    investment_costs : Mapping
        mapping from values of nominal size to investment costs
    part_load : Mapping
        mapping from values of relative output to relative input
    c_m : Expression
        percentage of investment costs incurred annually as maintenance costs
    eta_N : Expression
        efficiency when providing nominal output

    Examples
    ========

    investment_costs coefficients
    -----------------------------

    voll2013automated (Table A.2, p. 149)

    technology V_N_1 [kW] I_1 [k€] V_N_2 [kW] I_2 [k€] V_N_3 [kW] I_3 [k€]
    ----------------------------------------------------------------------
    Boiler            100   34.343      14000  379.580         -         -
    CHP engine        500  230.022        712  278.644      3200   850.563
    AC                 50   68.493        750  154.012      6500   522.651
    TC                400   89.006      10000 1572.302         -         -

    sass2020model (Table C.2, p. 17)

    Comp. V_N_1      I_1 V_N_2      I_2  V_N_3       I_3 V_N_4       I_4
           [kW]     [k€]  [kW]     [k€]   [kW]      [k€]  [kW]      [k€]
    BOI     100  27.6807  2000   86.955      –         –     –        –
    CHP 1   100 138.090   1400  436.869      –         –     –        –
    CHP 2  1400 436.869   2300  643.716      –         –     –        –
    CHP 3  2300 643.716   3200  850.563      –         –     –        –
    AC      100  71.404    711  155.9175  2000  243.252      –        –
    CC      400  89.006  10000 1572.302      –         –     –        –
    PV        5  32.8583    55  295.6515   550 2187.492      –        –
    HP        5   5.1137    27   14.9994    83   31.2601   200   55.486
    BAT       0   0         40   47.323    120  118.368   2000 1238.006
    TES_c     0   0         20   14.492   1000   35.264  25000  543.968
    TES_h     0   0         20    1.929  23175  502.079 115000 2086.885

    part_load coefficients
    ----------------------

    sass2020model (note that u and v are swapped!)

    Component u1  v1      u2      v2      u3 v3
    -------------------------------------------
    BOI       0.2 0.22608 1       1       –  –
    CHP (th)  0.5 0.46035 1       1       –  –
    CHP (el)  0.1 0.20251 1       1       –  –
    AC        0.2 0.25006 0.60778 0.48792 1  1
    CC        0.2 0.31204 0.70497 0.59543 1  1
    HP        0.2 0.21584 1       1       –  –

    voll2013automated (Table A.3, p. 149)

    technology v1  u1     v2  u2     v3  u3
    -------------------------------------------
    Boiler     0.2 0.2184 1.0 1.0004 -   -
    CHP engine 0.5 0.4790 1.0 0.9815 -   -
    AC         0.2 0.2722 0.6 0.4833 1.0 0.9833
    TC         0.2 0.3185 0.7 0.5936 1.0 0.9828

    c_m
    ---

    voll2013automated (Table 3.1, p. 32)

    technology thermal power investment cost maintenance cost      ηN
               range [MW]    range [k€]      [% investment cost]   []
    -----------------------------------------------------------------
    Boiler        0.1 - 14.0        34 - 380                1.5  0.9
    CHP engine     0.5 - 3.2       230 - 850               10.0  0.87
    AC             0.1 - 6.5        75 - 520                4.0  0.67
    TC            0.4 - 10.0       89 - 1570                1.0  5.54
    """

    def __init__(self, label, investment_costs, part_load, c_m, eta_N):
        super().__init__(label)
        from comando.linearization import piecewise_linear
        # TODO: handle special case of \en(elements) == 1

        # A mapping from nominal size to investment costs
        V_N_bp, I0, dIdVN = piecewise_linear(investment_costs)
        # V_N_bp = sorted(investment_costs.keys())
        # I = [investment_costs[V_N_bp_i] for i in V_N_bp]
        # dIdVn = [(I_ub - I_lb) / (V_N_ub_i - V_N_lb_i)
        #          for V_N_lb_i, V_N_ub_i, I_lb, I_ub
        #          in zip(V_N_bp, V_N_bp[1:], I, I[1:])]
        # I0 = []
        i_elements = range(1, len(V_N_bp))

        # A mapping from relative output (V/V_N) to relative input (U/U_N)
        v_bp, u0, dudv = piecewise_linear(part_load)
        # v = sorted(part_load.keys())
        # u = [part_load[v_i] for v_i in v]
        # dudv = [(u_ub - u_lb) / (v_ub - v_lb)
        #         for v_lb, v_ub, u_lb, u_ub
        #         in zip(v, v[1:], u, u[1:])]
        # u0 = [u_i - dudv_i * v_i for u_i, dudv_i, v_i in zip(u, dudv, v)]
        elements = range(1, len(v_bp))

        V_N = self.make_design_variable('out_nom')
        V_N_ = [self.make_operational_variable(f'out_nom_{i}')
                for i in elements]  # NOTE: ξ_t in source!
        gamma = [self.make_design_variable(f'gamma_{i}', BINARY)
                 for i in elements]
        # y = self.make_design_variable('existence', BINARY)
        y = self.add_expression('exists', sum(gamma))

        delta = self.make_operational_variable('on', BINARY)
        V = [self.make_operational_variable(f'out_{i}')
             for i in elements]
        V_t = self.make_operational_variable('out')

        # Enforce disjunction
        self.add_le_constraint(sum(delta), 1, 'one_output')
        self.add_eq_constraint(V_t, sum(V), 'output_matches')
        # NOTE: Alternatively V_t = sum(V)

        inv = sum(gamma_i * I0_i + V_N * dIdVN_i
                  for gamma_i, I0_i, dIdVN_i in zip(gamma, I0, dIdVN))  # 3.6
        # set (investment) costs
        self.add_expression('investment_costs', inv)
        # add maintenance costs to fixed costs
        self.add_expression('fixed_costs', c_m * inv)

        self.add_le_constraint(y, 1, 'one_investment_size')  # 3.7
        for i, gamma_i, V_N_lb_i, V_N_ub_i \
                in zip(i_elements, gamma, V_N_bp, V_N_bp[1:]):  # 3.8 & 3.9
            self.add_le_constraint(gamma_i * V_N_lb_i, V_N, f'V_N_lb_{i}')
            self.add_le_constraint(V_N, gamma_i * V_N_ub_i, f'V_N_ub_{i}')

        # U = u * V_N / eta_nom
        # U0 = [u0_i * V_N_i for u0_i, V_N_i in zip(u0, V_N_)]

        # 3.16 after reformulation of δ_t * V_N -> ξ_t (here V_N_)
        U_t = sum(u0_i * V_N_i + dudv_i * V_i
                  for u0_i, V_N_i, dudv_i, V_i in zip(u0, V_N_, dudv, V)) \
            / eta_N

        v_lb = v_bp[0]
        v_ub = v_bp[-1]
        V_N_lb = V_N_bp[0]
        V_N_ub = V_N_bp[-1]
        for i, V_N_i in zip(elements, V_N_):
            # 3.17
            self.add_le_constraint(V_N_i * v_lb, V_t, f'min_out_{i}')
            self.add_le_constraint(V_t, V_N_i * v_ub, f'max_out_{i}')
            # 3.19 (Glover 1 & 2)
            self.add_le_constraint(delta * V_N_lb, V_N_i, f'V_N_{i}_lb')
            self.add_le_constraint(V_N_i, delta * V_N_ub, f'V_N_{i}_ub')
            # 3.20 (Glover 3 & 4)
            self.add_le_constraint((1 - delta) * V_N_lb, V_N - V_N_i,
                                   f'V_N_{i}_V_N_lb')
            self.add_le_constraint(V_N - V_N_i, (1 - delta) * V_N_ub,
                                   f'V_N_{i}_V_N_ub')
        # set connectors
        self.add_input('IN', U_t)
        self.add_output('OUT', V_t)
