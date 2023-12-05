"""Component models for the low temperature district heating network."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Dominik Hering, Marco Langiu
from comando.core import Component, BINARY, System


class BESMFlowTFixDT(Component):
    """A model for a building energy system.

    using mass flow and temperature connectors. Temperature difference is fixed

    Parameters
    ----------
    Qdot:
        The building's heat load, [kW]
    cp:
        Heating capacity, default=4.12 [kJ/kgK]
    dT:
        Temperature difference between t_flow and t_return, [K]
    t_flow:
        Flow temperature, [K]

    Connectors
    ----------
    IN_mflow:      Mass flow rate
    """

    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        qdot = self.make_parameter('Qdot')  # [W]
        cp = self.make_parameter('cp', 4.12)
        d_t = self.make_parameter('dT')
        self.make_parameter('T_flow')

        m_flow = qdot / (d_t * cp)

        self.add_connector('IN_mflow', m_flow)


class DummySource(Component):
    """A dummy resource that serves as a source for an arbitrary commodity.

    Parameters
    ----------
    price:
        float, defines the price per energy unit of the dummy source.

    Connectors
    ----------
    OUT:
        use of resource

    Expressions
    -----------
    variable_costs:
        Cost of for use
    """

    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        price = self.make_parameter('price')

        #######################################################################
        # Operational Variables
        #######################################################################
        use = self.make_operational_variable('use', bounds=(0, None),
                                             init_val=0)
        self.add_output('OUT', use)

        self.add_expression('variable_costs', price * use)


class HeatPump(Component):
    """Quadratic Heat pump model.

    Based on the following temperatures:
    - t_in_evap:
        Incoming temperature at evaporator [K]
    - t_out_evap:
        Outgoing temperature at evaporator [K]
    - t_out_cond:
        Outgoing temperature at condenser [K]
    - t_in_cond:
        Incoming temperature at condenser [K]

    Parameters
    ----------
    cop_nom : Float, COP efficiency.
        Determines the performance of the heat pump relative to Carnot in
        percent, default=0.6
    n:
        Integer, Number of HP instances, default = 1
    cp:
        Heat capacity of fluid, default = 4.12 kJ/(kgK)

    Design Variables
    ----------------
    b_build:
        Boolean build decision of HP
    Qdot_design:
        Maximum heating capacity [kW]

    Operational Variables
    ---------------------
    b_op:
        Boolean for operation (1 for on, 0 for off)
    mflow_evap:
        Mass flow rate at evaporator [kg/s]
    mflow_cond:
        Mass flow rate at condenser [kg/s]
    Qdot_cond:
        Thermal power at condenser [kW]
    P_HP:
        Electrical power demand of HP [kW]

    Connectors
    ----------
    IN_mflow_evap:
        mass flow rate at evaporator [kg/s]
    IN_P:
        Electrical Power consumption [kW]
    OUT_mflow_cond:
        mass flow rate at Condenser [kg/s]

    Expressions
    -----------
        investment_costs:    Investment cost of HP
    """

    Q_max = 400  # [kW]
    mdot_max = 10  # [kg/s]

    def __init__(self, label, t_in_evap, t_out_evap, t_in_cond,
                 t_out_cond, b_connect=None):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        cop_nom = self.make_parameter('cop_nom', 0.6)
        cp = self.make_parameter('cp', 4.12)
        n = self.make_parameter('n', 1)

        #######################################################################
        # Design Variables
        #######################################################################
        qdot_design = self.make_design_variable('Qdot_design',
                                                bounds=(0, self.Q_max))
        b_build = self.make_design_variable('b_build', domain=BINARY,
                                            bounds=(0, 1))
        self.add_le_constraint(qdot_design, self.Q_max * b_build, 'HP_Q_max')

        #######################################################################
        # Operational Variable
        #######################################################################
        b_op = self.make_operational_variable('b_op', domain=BINARY)
        self.add_le_constraint(b_op, b_build, 'HP_operation')
        # Evaporator
        mflow_evap = self.make_operational_variable('mflow_evap',
                                                    bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_evap, b_op * self.mdot_max,
                               'HP_mlfow_evap')
        self.add_le_constraint(mflow_evap, b_build * self.mdot_max,
                               'HP_mlfow_evap')
        self.add_input(identifier='IN_mflow_evap', expr=mflow_evap * n)

        qdot_evap = mflow_evap * (t_in_evap - t_out_evap) * cp

        # Condenser
        mflow_cond = self.make_operational_variable('mflow_cond',
                                                    bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_cond, b_op * self.mdot_max,
                               'HP_mflow_cond_op')
        self.add_le_constraint(mflow_cond, b_build * self.mdot_max,
                               'HP_mflow_cond_op')
        self.add_output(identifier='OUT_mflow_cond', expr=mflow_cond * n)

        # HP
        qdot_cond = self.make_operational_variable('Qdot_cond',
                                                   bounds=(0, self.Q_max))
        self.add_le_constraint(qdot_cond, self.Q_max * b_op, 'Qdot_HP_op')

        p_hp = self.make_operational_variable('P_HP', bounds=(0, self.Q_max))
        self.add_le_constraint(p_hp, self.Q_max * b_op, 'P_HP_op')
        self.add_eq_constraint(p_hp * t_out_cond * cop_nom,
                               (t_out_cond - t_out_evap) * qdot_cond,
                               'input_output_relation')
        # n_HP Heat pumps are connected to the grid in parallel
        self.add_input(identifier='IN_P', expr=p_hp * n)
        self.add_le_constraint(qdot_cond, qdot_design, 'output_limit')

        qdot_out = mflow_cond * cp * (t_out_cond - t_in_cond)
        self.add_eq_constraint(qdot_cond, qdot_out, 'heat_flow')
        self.add_eq_constraint(qdot_evap + p_hp, qdot_cond, 'energy_balance')

        #######################################################################
        # Misc
        #######################################################################
        # Investment costs
        if b_connect is None:
            inv_costs = 500 * qdot_design * n
        else:
            inv_costs = 620 * qdot_design * n
            self.add_le_constraint(b_build, b_connect, 'connection')
        self.add_expression('investment_costs', inv_costs * n)


class HeatExchanger(Component):
    """Heat exchanger model.

    Based on the following temperatures:
    - t_in_a:
        Incoming temperature at side a [K]
    - t_out_a:
        Outgoing temperature at side a [K]
    - t_out_b:
        Outgoing temperature at side b [K]
    - t_in_b:
        Incoming temperature at side b [K]

    Parameters
    ----------
    n:
        Integer, Number of HP instances, default = 1
    cp:
        Heat capacity of fluid, default = 4.12 kJ/(kgK)

    Design Variables
    ----------------
    b_build:
        Boolean for build decision of HX
    Qdot_design:
        Maximum heating capacity [kW]

    Operational Variables
    ---------------------
    b_op:
        Boolean operational variable.
        b_op=0 for no mass flow and no temperature constraints
        b_op=1 for heat exchanger operation
    mflow_a:
        Incoming mass flow rate at side a [kg/s]
    mflow_b:
        Outgoing mass flow rate at side b [kg/s]

    Connectors
    ----------
    IN_mflow_a:
        mass flow rate at side a [kg/s]
    OUT_mflow_b:
        mass flow rate at side b [kg/s]

    Expressions
    -----------
        investment_costs:    Investment cost of HX
    """

    Q_max = 400  # [kW]
    mdot_max = 10  # [kg/s]

    def __init__(self, label, b_connect=None, b_build_hp=None,
                 t_in_a=None, t_out_a=None, t_in_b=None, t_out_b=None):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        cp = self.make_parameter('cp', 4.12)
        n = self.make_parameter('n', 1)

        #######################################################################
        # Design Variables
        #######################################################################
        qdot_design = self.make_design_variable('Qdot_design',
                                                bounds=(0, self.Q_max))
        b_build = self.make_design_variable('b_build', domain=BINARY)
        if b_build_hp is not None:
            # Opposite of HP build decision
            self.add_le_constraint(b_build, 1 - b_build_hp)
        self.add_le_constraint(qdot_design, 400 * b_build, 'Qdot_max')
        if b_connect is not None:
            self.add_le_constraint(b_build, b_connect, 'connection')

        #######################################################################
        # Operational Variable
        #######################################################################
        b_op = self.make_operational_variable('b_op', domain=BINARY)
        self.add_le_constraint(b_op, b_build, 'operation')
        # Side a
        mflow_a = self.make_operational_variable('mflow_a',
                                                 bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_a, b_op * self.mdot_max, 'b_hx_mflow_a')
        self.add_le_constraint(mflow_a, b_build * self.mdot_max,
                               'bbuild_hx_mflow_a')
        self.add_input(identifier='IN_mflow_a', expr=mflow_a * n)

        qdot_in = mflow_a * (t_in_a - t_out_a) * cp

        # Side b
        mflow_b = self.make_operational_variable('mflow_b',
                                                 bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_b, b_op * self.mdot_max, 'b_hx_mflow_b')
        self.add_le_constraint(mflow_b, b_build * self.mdot_max,
                               'bbuild_hx_mflow_b')
        self.add_output(identifier='OUT_mflow_b', expr=mflow_b * n)

        qdot_out = mflow_b * cp * (t_out_b - t_in_b)

        self.add_le_constraint(qdot_out, qdot_design, 'design_Limit')
        self.add_le_constraint(t_out_b * b_op, t_in_a * b_op,
                               'max_temp_increase')

        self.add_eq_constraint(qdot_in, qdot_out, 'heat_flow')

        #######################################################################
        # Misc
        #######################################################################
        # Investment costs
        inv_costs = 90 * qdot_design * n
        self.add_expression('investment_costs', inv_costs)


class HeatSourceDecentral(Component):
    """A model for a generic decentral heat source.

    Using mass flow and temperature connectors. Temperature difference is
    fixed.

    Parameters
    ----------
    cp:
        Heating capacity, default=4.12 [kJ/kgK]
    dT:
        Temperature difference between t_flow and t_return, [K]
    efficiency:
        Efficiency
    n:
        number of parallel units
    p_spec:
        Specific price of component  [€/kW]
    p_fix:
        Fix price of component [€]

    Design Variables
    ----------------
    b_build:
        Boolean for build decision of HS

    Connectors
    ----------
    OUT_mflow:
        Mass flow rate
    IN_P:
        Input for power consumption

    Expressions
    -----------
    investment_costs:    Investment cost of Heat source
    """

    Q_max = 20  # [kW]
    mdot_max = 10  # [kg/s]

    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        cp = self.make_parameter('cp', 4.12)
        spec_price = self.make_parameter('p_spec', 100)  # €/kW
        fix_price = self.make_parameter('p_fix', 100)  # €
        d_t = self.make_parameter('dT')
        eff = self.make_parameter('efficiency', value=1)
        n = self.make_parameter('n', value=1)

        qdot_design = self.make_design_variable('Qdot_max',
                                                bounds=(0, self.Q_max))
        b_build = self.make_design_variable('b_build', domain=BINARY)
        self.add_le_constraint(qdot_design, b_build * self.Q_max, 'Qdot_max')
        #######################################################################
        # Operational Variables
        #######################################################################
        m_flow = self.make_operational_variable('mflow',
                                                bounds=(0, self.mdot_max))
        self.add_le_constraint(m_flow, b_build * self.mdot_max, 'b_mflow')
        self.add_output('OUT_mflow', m_flow * n)
        p_in = m_flow * (d_t * cp)
        self.add_le_constraint(p_in, qdot_design, 'p_max')
        self.add_input('IN_P', p_in * n / eff)

        # Investment costs
        inv_costs = (spec_price * qdot_design + b_build * fix_price) * n
        self.add_expression('investment_costs', inv_costs)


class Network(Component):
    """A network component.

    Collects and distributes energy flows. Temperature losses are calculated.

    Parameters
    ----------
    U_avg:
        Average U value of network [W/K]
    cp:
        Heating capacity of medium, default=4.12 [kJ/kgK]
    Design_dT:
        Design temperature difference of network, default=15 [K]
    T_amb:
        Ambient air temperature. Is used to generate T_flow curve, [K]
    T_amb_design:
        Design ambient air temprature of lowest point in heating curve,
        default= -16 °C
    T_amb_tresh:
        Treshold temperature when no heating is required any more,
        default= 20 °C

    Design Variables
    ----------------
    b_{nodes}:
        Boolean decision variables for each possible pipe segment
    Return_Temp_max:
        Design Return temperature at T_amb_design [K]
    Return_Temp_min:
        Design Return temperature at T_amb_tresh [K]

    Operational Variables
    ---------------------
    mflow:
        Mass flow rate in network [kg/s]
    b_op:
        Boolean for operation (1 for On, 0 for off)
    T_in_a:
        Incoming Temperature at side a, [K]
    T_out_a:
        Outgoing Temperature at side a, [K]
    T_in_b:
        Incoming Temperature at side b, [K]
    T_out_b:
        Outgoing Temperature at side b, [K]

    Connectors
    ----------
    IN_mflow:
        Mass flow rate through network [kg/s]
    OUT_mflow:
        Mass flow rate through network [kg/s]

    Expressions
    -----------
    investment_costs:    Investment cost of pipes
    """

    def __init__(self, label, pipe_data):

        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        u_avg = self.make_parameter('U_avg', 0.035)  # [W/mK]
        cp = self.make_parameter('cp', 4.12)
        # Temperature diff at source
        design_dt = self.make_parameter('Design_dT', 15)
        t_amb = self.make_parameter('T_amb')  # [K]
        t_amb_design = self.make_parameter('T_amb_design', 273.15-12)  # [K]
        t_amb_tresh = self.make_parameter('T_amb_tresh', 273.15+20)  # [K]

        #######################################################################
        # Design Variables
        #######################################################################
        # Flow temperature curve
        t_return_max = self.make_design_variable('Return_temp_max',
                                                 bounds=(273.15, 373.15))
        t_return_min = self.make_design_variable('Return_temp_min',
                                                 bounds=(273.15, 373.15))
        self.add_le_constraint(t_return_min, t_return_max,
                               'Design_temperature')

        for b_node in pipe_data['Beginning Node']:
            self.make_design_variable(f'b_{b_node}', domain=BINARY)
            if 'SimpleDistrict' in b_node:
                self.add_input(identifier=f'IN_b_{b_node}',
                               expr=self[f'b_{b_node}'])
        # set dependencies
        pipe_ref_data = pipe_data[['Beginning Node', 'Depend', 'cost',
                                   'Length [m]']].set_index('Beginning Node')
        for b_node in pipe_data['Beginning Node']:
            if type(pipe_ref_data['Depend'][b_node]) is str:
                dependencies_list = list(pipe_ref_data['Depend'][b_node])
                for dependency in dependencies_list:
                    self.add_le_constraint(self[f'b_{b_node}'],
                                           self[f'b_{dependency}'],
                                           'b_connect_dependencies')

        #######################################################################
        # Operational Variables
        #######################################################################
        mflow = self.make_operational_variable("mflow", bounds=(0, 10))
        b_op = self.make_operational_variable("b_op", domain=BINARY)
        self.add_le_constraint(mflow, 10 * b_op, 'bop_mflow')
        self.add_input(identifier='IN_mflow', expr=mflow)
        self.add_output(identifier='OUT_mflow', expr=mflow)

        # temperature curve
        # Slope and y-axis of linear temperature correlation
        m = (t_return_max - t_return_min) / (t_amb_design - t_amb_tresh)
        y = t_return_max - m * t_amb_design

        t_in_a = self.make_operational_variable('T_in_a',
                                                bounds=(273.15, 373.15))
        t_out_a = self.make_operational_variable('T_out_a',
                                                 bounds=(273.15, 373.15))
        self.add_eq_constraint(t_out_a,
                               (m * t_amb + y) * b_op + 281.15 * (1 - b_op),
                               'Treturn_heating_curve')
        self.add_eq_constraint(t_in_a,
                               (t_out_a + design_dt) * b_op
                               + 281.15 * (1 - b_op),
                               'Tflow_heating_curve')

        # Temperature with losses
        t_out_b = self.make_operational_variable('T_out_b',
                                                 bounds=(273.15, 373.15))
        t_in_b = self.make_operational_variable('T_in_b',
                                                bounds=(273.15, 373.15))
        self.add_le_constraint(t_out_b, t_in_a, 'dT_flow')
        self.add_le_constraint(t_in_a - t_out_b, 5 * b_op, 'dT_flow_max')
        self.add_ge_constraint(t_in_a - t_out_b, 0, 'dT_flow_min')
        self.add_le_constraint(t_out_a, t_in_b, 'dT_return')
        self.add_le_constraint(t_in_b - t_out_a, 5 * b_op, 'dT_return_max')
        self.add_ge_constraint(t_in_b - t_out_a, 0, 'dT_return_min')
        length = sum(
            pipe_ref_data['Length [m]'][b_node] * self[f'b_{b_node}']
            for b_node in pipe_data['Beginning Node'])
        self.add_eq_constraint(mflow * cp * 1000 * (t_in_a - t_out_b),
                               (t_in_a - 281.15) * length * u_avg,
                               'flow_losses')
        self.add_eq_constraint(mflow * cp * 1000 * (t_in_b - t_out_a),
                               (t_in_b - 281.15) * length * u_avg,
                               'return_losses')

        #######################################################################
        # Investment Costs
        #######################################################################
        inv_costs = sum(
            pipe_ref_data['cost'][b_node] * self[f'b_{b_node}']
            for b_node in pipe_data['Beginning Node'])
        self.add_expression('investment_costs', inv_costs)


class WasteHeatMflowT(Component):
    """A resource that serves as a source for waste heat.

    Parameters
    ----------
    T_flow:
        Flow Temperature [K]
    T_return:
        return temperature of waste heat [K]

    Operational Variables
    ---------------------
    m_flow
        Mass flow rate, [kg/s]

    Connectors
    ----------
    OUT_mflow:
        Outgoing mass flow
    """

    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        # Temperatures
        self.make_parameter('T_flow')
        self.make_parameter('T_return')
        # Mass flow rate
        mflow = self.make_operational_variable('mflow', bounds=(0, None))
        self.add_output(identifier='OUT_mflow', expr=mflow)


class LinkingComponent(System):
    """Create a linking component, modeled as a COMANDO system.

    This component contains a Heat Pump and a Heat Exchanger model.
    """

    def __init__(self, label, t_in_a, t_out_a,
                 t_in_b, t_out_b, b_connect=None):
        super().__init__(label)
        hp = HeatPump(f'HP_{label}', b_connect=b_connect,
                      t_in_evap=t_in_a, t_out_evap=t_out_a,
                      t_in_cond=t_in_b, t_out_cond=t_out_b)
        hx = HeatExchanger(f'HX_{label}', b_build_hp=hp['b_build'],
                           t_in_a=t_in_a, t_out_a=t_out_a,
                           t_in_b=t_in_b, t_out_b=t_out_b,
                           b_connect=b_connect)
        for comp in [hp, hx]:
            self.add(comp)
        # Side a
        self.connect('IN_mflow_a', [hp.IN_mflow_evap,
                                    hx.IN_mflow_a])
        # Side b
        self.connect('OUT_mflow_b', [hp.OUT_mflow_cond,
                                     hx.OUT_mflow_b])

        self.expose_connector(hp.IN_P, 'IN_P')


class Consumer(System):
    """Create a consumer group, modeled as a COMANDO system.

    This component contains a linking component, two heat sources and one
    demand.
    """

    def __init__(self, label, t_in_a, t_out_a, b_connect=None):
        super().__init__(label)
        hs_el = HeatSourceDecentral(f'HS_el_{label}')
        hs_gas = HeatSourceDecentral(f'HS_gas_{label}')
        self.add_le_constraint(hs_el['b_build'] + hs_gas['b_build'], 1)
        building = BESMFlowTFixDT(f'BES_{label}')
        lc = LinkingComponent(f'LC_{label}', b_connect=b_connect,
                              t_in_a=t_in_a, t_out_a=t_out_a,
                              t_in_b=building['T_flow'] - building['dT'],
                              t_out_b=building['T_flow'])
        lc.extend_connection('OUT_mflow_b')
        lc.extend_connection('IN_mflow_a')

        for comp in [lc, hs_el, hs_gas, building]:
            self.add(comp)
        self.connect('mflow_to_bes', [lc.OUT_mflow_b,
                                      hs_el.OUT_mflow,
                                      hs_gas.OUT_mflow,
                                      building.IN_mflow])
        self.connect('IN_P_el', [lc.IN_P,
                                 hs_el.IN_P])
        self.expose_connector(hs_gas.IN_P, 'IN_P_gas')

        # Net connection
        self.expose_connector(lc.IN_mflow_a, 'IN_mflow_a')
