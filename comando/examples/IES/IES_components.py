"""Components for the IES case study."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu, David Shu, Florian Joseph Baader
from comando.core import System, Component, BINARY


class IESComponent(Component):
    """Component of an industrial energy system.

    Arguments
    ---------
    - label : str
        Unique sting that serves as an identifier of this Demand.
    - nom_ref: The nominal reference size of the component
    - c_ref: coefficient of investment costs
    - c_m: factor of investment costs corresponding to maintenance costs
          TODO: This implies calculations with a year as the time series!
    - M: (1) exponent of investment costs
    - nom_min: (0) lower bound for nominal size
    - nom_max: (None) upper bound for nominal size
    - min_part_load: (0) minimum allowed value for relative output
    - base_eff: the base efficiency of the component
    - fit_params: dict of power, coefficient items for polynomial
        approximation of efficiency
    - in_name: name of input commodity
    - out_name: name of output commodity
    - in_connector_name: name of input connector
    - out_connector_name: name of output connector
    - exists: bool
        Specifies whether the Component is currently installed; if set to
        True, its nominal size will be a user-specified parameter. If set
        to False, the nominal size of the component is a design variable.
    - optional: bool
        Specifies whether the Component may be added or removed.
        If exists is True, setting optional to True allows the Component to
        be sold for a user-specified price.
        If exists is False, setting optional to False requires the
        component to be installed with a minimal nominal size of at least
        `nom_min`, while setting it to True adds the posibility of not
        installing it.
    """

    def __init__(self, label, nom_ref, c_ref, c_m, M=1, nom_min=0,
                 nom_max=None, min_part_load=0, base_eff=1,
                 fit_params_nom=None, fit_params_den=None,
                 in_name='input', out_name='output', in_connector_name='IN',
                 out_connector_name='OUT', exists=False, optional=True):
        # If no fit_params are given we assume no dependency on relative output
        if not fit_params_nom:
            fit_params_nom = {0: 1}
        if not fit_params_den:
            fit_params_den = {0: 1}

        super().__init__(label)

        # nominal (installed) output power
        nom_name = out_name + '_nom'
        available = 1
        if exists:
            if optional:  # Component exits but may be sold
                installed = self.make_parameter(nom_name)
                self.add_le_constraint(nom_min, installed,
                                       name=nom_name + '_min')
                self.add_le_constraint(installed, nom_max,
                                       name=nom_name + '_max')
                keep = self.make_design_variable('keep', BINARY, init_val=1)
                available = keep
                out_nom = keep * installed

                # IDEA: Add a discount factor, that would make more sense!
                # TODO: Is that a bad formulation?
                revenue = c_ref * ((installed * (1 - keep)) / nom_ref) ** M
                # revenue from selling counts as negative investment costs
                self.add_expression('investment_costs', -revenue)
            else:  # Component exists and may NOT be sold
                # TODO: This shows that we might want to follow GAMS, Pyomo,
                #       etc. and introduce bounds on parameters!
                out_nom = self.make_parameter(nom_name)
                self.add_le_constraint(nom_min, out_nom,
                                       name=nom_name + '_min')
                self.add_le_constraint(out_nom, nom_max,
                                       name=nom_name + '_max')

            # investment cost of the component
            inv = c_ref * (out_nom / nom_ref) ** M
            # add maintenance costs to fixed costs
            self.add_expression('fixed_costs', c_m * inv)
        else:
            # TODO: If the component does NOT exist we might add a method to
            #       the component to specify a finite number of possible
            #       nominal sizes.
            out_nom = self.make_design_variable(nom_name,
                                                bounds=(0, nom_max),
                                                init_val=nom_max)
            if optional:  # Component doesn't exist but may be installed or not
                exists = self.make_design_variable('exists', BINARY,
                                                   init_val=1)
                available = exists
                self.add_le_constraint(exists * nom_min, out_nom,
                                       name=nom_name + '_min')
                self.add_le_constraint(out_nom, exists * nom_max,
                                       name=nom_name + '_max')
            else:
                out_nom.lb = nom_min
            # NOTE: if optional is false, the component must be installed

            # investment cost of the component
            inv = c_ref * (out_nom / nom_ref) ** M
            # fixed cost of the component
            fc = c_m * inv

        self.add_expression('investment_costs', inv)
        self.add_expression('fixed_costs', fc)

        # relative output
        out_rel = self.make_operational_variable(out_name + '_rel',
                                                 bounds=(0, 1), init_val=1)

        # If minimal part load is considered, we need to decide whether it's
        # operational or not to determine the lower bound for the part load.
        self.add_expression('min_part_load', min_part_load)
        if min_part_load == 0:
            operating = 1
        else:
            # binary indicating whether the component is operational
            operating = self.make_operational_variable('operating',
                                                       domain=BINARY,
                                                       init_val=1)
            if optional:  # available is variable
                self.add_le_constraint(operating, available,
                                       name='availablility')
            self.add_ge_constraint(out_rel, operating * min_part_load,
                                   name='min_part_load')
            self.add_le_constraint(out_rel, operating,
                                   name='max_part_load')

        # efficiency in terms of relative output
        numerator = sum(coeff * out_rel ** power
                        for power, coeff in fit_params_nom.items())
        denominator = sum(coeff * out_rel ** power
                          for power, coeff in fit_params_den.items())
        eff = base_eff * numerator/denominator
        self.add_expression('eff', eff)

        out = self.make_operational_variable(out_name,
                                             bounds=(0, nom_max),
                                             init_val=nom_max)
        self.add_eq_constraint(out, out_rel * out_nom,
                               name='component_output')
        inp_max = nom_max/eff.subs(out_rel, 1)
        try:
            init_ub = float(inp_max)
        except (ValueError, TypeError):  # base_eff contains symbols!
            init_ub = None
        inp = self.make_operational_variable(in_name, bounds=(0, init_ub),
                                             init_val=init_ub)
        # Redundant constraints (add as relaxation only?)
        self.add_le_constraint(inp, out_nom/eff.subs(out_rel, 1),
                               name='input_limit')
        self.add_le_constraint(inp, operating * inp_max,
                               name='input_upper_limit')
        self.add_eq_constraint(inp * eff, out, name='input_output_relation')

        self.add_expression('input', inp)
        self.add_expression('output', out)

        # set connectors
        self.add_input(in_connector_name, inp)
        self.add_output(out_connector_name, out)


class Boiler(IESComponent):
    """Boiler parameterized as in Sass2020.

    design variable: nominal output power to be installed
    operational variable: heat output during a given time step
    The investment costs 'c_inv' is a nonlinear function of the nominal output
    'Qdot_nom':
    c_inv = c_inv,ref * (Qdot_nom/Qdot_ref)^M
    The input 'Qdot_in' can be determined via the efficiency relation:
    Qdot_in =  Qdot_out / eff
    The efficiency 'eff' is determined via the product of the nominal
    efficiency and a polynomial fitting function of the load fraction 'q':
    q = Qdot_out/Qdot_nom
    """

    Qdot_ref = .001  # [MW], reference nominal power
    c_ref = 2701.6/1000  # [k€], reference cost
    M = 0.4502  # [-], cost exponent
    c_m = 0.015  # [-], maintenance coefficient, (fraction of investment cost)
    Qdot_min = 100/1000  # [kW], minimal nominal power allowed for the model
    Qdot_max = 14000/1000  # [kW], maximal nominal power allowed for the model
    qdot_min = 0.2  # [-] minimum modeled thermal output part load
    eff_nom = 0.9  # [-], nominal efficiency
    fit_params_nom = {0: -0.07557,
                      1: 1.39731,
                      2: -7.00130,
                      3: 21.75378}
    fit_params_den = {0: 0.03487,
                      1: 0.67774,
                      2: -5.34196,
                      3: 20.66646}

    def __init__(self, label, exists=False, optional=True):
        super().__init__(label, self.Qdot_ref, self.c_ref, self.c_m, M=self.M,
                         nom_min=self.Qdot_min, nom_max=self.Qdot_max,
                         min_part_load=self.qdot_min, base_eff=self.eff_nom,
                         fit_params_nom=self.fit_params_nom,
                         fit_params_den=self.fit_params_den, in_name='Qdot_in',
                         out_name='Qdot_out', exists=exists, optional=optional)


class AbsorptionChiller(IESComponent):
    """Absorption chiller (AC) parameterized as in Sass2020.

    design variable: nominal cooling power to be installed
    operational variable: operational output cooling for each time step
    The investment costs 'c_inv' is a nonlinear function of the nominal output
    Qdot_nom:
    c_inv = c_inv,ref * (Qdot_nom/Qdot_ref)^M
    The input 'Qdot_in' can be determined via the efficiency relation:
    Qdot_in =  Qdot_out / eff
    The efficiency 'eff' is determined via the product of the nominal
    efficiency and a polynomial fitting function of the load fraction 'q':
    q = Qdot_out/Qdot_nom
    """

    Qdot_ref = .001  # [MW], reference nominal power
    c_ref = 8847.5/1000  # [k€], reference cost
    M = 0.4345  # [-], cost exponent
    c_m = 0.01  # [-], maintenance coefficient, (fraction of investment cost)
    Qdot_min = 100/1000  # [MW], minimal nominal power allowed for the model
    Qdot_max = 2000/1000  # [MW], maximal nominal power allowed for the model
    qdot_min = 0.2  # [-] minimum output part load
    cop_nom = 0.67  # [-], nominal efficiency
    fit_params_nom = {1: 1}
    fit_params_den = {0: 0.24999,
                      1: -0.08330,
                      2: 0.83330}

    def __init__(self, label, exists=False, optional=True):
        super().__init__(label, self.Qdot_ref, self.c_ref, self.c_m, M=self.M,
                         nom_min=self.Qdot_min, nom_max=self.Qdot_max,
                         min_part_load=self.qdot_min, base_eff=self.cop_nom,
                         fit_params_nom=self.fit_params_nom,
                         fit_params_den=self.fit_params_den, in_name='Qdot_in',
                         out_name='Qdot_out', exists=exists, optional=optional)


class CompressionChiller(IESComponent):
    """Compression-chiller (CC) parameterized as in Sass2020.

    design variable: nominal cooling power to be installed
    operational variable: operational output cooling for each time step
    The investment costs 'c_inv' is a nonlinear function of the nominal output
    'Qdot_nom':
    c_inv = c_inv,ref * (Qdot_nom/Qdot_ref)^M
    The input 'P_IN' can be determined via the efficiency relation:
    P_IN =  Qdot_out / eff
    The efficiency 'eff' is determined via the product of the nominal
    efficiency and a polynomial fitting function of the load fraction 'q'
    q = Qdot_out/Qdot_nom
    """

    Qdot_ref = .001  # [MW], reference nominal power
    c_ref = 444.3/1000  # [k€], reference cost
    M = 0.8632  # [-], cost exponent
    c_m = 0.04  # [-], maintenance coefficient, (fraction of investment cost)
    Qdot_min = 400/1000  # [MW], minimal nominal power allowed for the model
    Qdot_max = 10000/1000  # [MW], maximal nominal power allowed for the model
    qdot_min = 0.2  # [-] minimum output part load
    cop_nom = 5.54  # [-], nominal efficiency
    # coefficients of rational function describing the part load behavior
    fit_params_nom = {0:  0.0126,
                      1:  3.679,
                      2: -3.5494,
                      3: 0.8615}

    # maintenance coefficient, fraction of the investment cost
    def __init__(self, label, exists=False, optional=True):
        super().__init__(label, self.Qdot_ref, self.c_ref, self.c_m, M=self.M,
                         nom_min=self.Qdot_min, nom_max=self.Qdot_max,
                         min_part_load=self.qdot_min, base_eff=self.cop_nom,
                         fit_params_nom=self.fit_params_nom, in_name='P_in',
                         out_name='Qdot_out', exists=exists, optional=optional)


class HeatPump(IESComponent):
    """Heat pump (HP) parameterized as in Sass2019.

    design variable: nominal heating power to be installed
    operational variable: operational output heating for each time step
    The investment costs 'c_inv' is a nonlinear function of the nominal output
    'Qdot_nom':
    c_inv = c_inv,ref * (Qdot_nom/Qdot_ref)^M
    The input 'P_IN' can be determined via the efficiency relation:
    P_IN =  Qdot_out / eff
    The efficiency 'eff' is determined via the product of the nominal
    efficiency and a polynomial fitting function of the load fraction 'q'
    q = Qdot_out/Qdot_nom
    """

    Qdot_ref = .001  # [MW], reference nominal power
    c_ref = 1654.7/1000  # [k€], reference cost
    M = 0.6611  # [-], cost exponent
    c_m = 0.01  # [-], maintenance coefficient, (fraction of investment cost)
    Qdot_min = 5/1000  # [MW], minimal nominal power allowed for the model
    Qdot_max = 200/1000  # [MW], maximal nominal power allowed for the model
    qdot_min = 0.2  # [-] minimum modeled thermal output part load

    def __init__(self, label, T_eva, T_con=273.15 + 60, eff_exer=0.36,
                 exists=False, optional=True):
        eff_carnot = 1 - T_eva/T_con
        super().__init__(label, self.Qdot_ref, self.c_ref, self.c_m, M=self.M,
                         nom_min=self.Qdot_min, nom_max=self.Qdot_max,
                         min_part_load=self.qdot_min,
                         base_eff=eff_exer/eff_carnot, in_name='P_in',
                         out_name='Qdot_out', exists=exists, optional=optional)


class CHP_sub(IESComponent):
    """Combined-heat-and-power (CHP) engine parameterized as in Sass2020.

    design variable: nominal heating power to be installed
    operational variable: operational output (heat & elec.) for each time step
    The investment costs 'c_inv' is a nonlinear function of the nominal thermal
    output 'Qdot_nom':
        c_inv = c_inv,ref * (Qdot_nom/Qdot_ref)^M
    The input 'Qdot_in' can be determined via the efficiency relation:
        Qdot_in =  Qdot_out / eff_th
    The thermal efficiency 'eff_th' is determined via the product of nominal
    thermal efficiency and a polynomial fitting function f(q), with the load
    fraction q = Qdot_out/Qdot_nom, giving the relationship:
        eff_th = f(q) * eff_th_nom
    The electric efficiency 'eff_el' is determined via the product of nominal
    electrical efficiency and f(q), giving the relationship:
        eff_el = f(q) * eff_el_nom
    Thus the output power can be computed as:
        P_out = Qdot_in * eff_el = Qdot_out / eff_th * eff_el
    """

    # investment model parameters
    Qdot_ref = .001  # [MW], reference nominal power
    c_ref = 9332.6/1000  # [k€], reference cost
    M = 0.539  # [-], cost exponent
    c_m = 0.1  # [-], maintenance coefficient, (fraction of investment cost)
    Qdot_min = 100/1000  # [MW], minimal nominal power allowed for the model
    Qdot_max = 3200/1000  # [MW], maximal nominal power allowed for the model
    qdot_min = 0.5  # [-] minimum output part load
    # coefficients for polynomial fits of thermal
    fit_params_th = {0: 1.0960, 1: -0.0199, 2: -0.0768}
    fit_params_el = {0: 0.5868, 1: 0.6743, 2: -0.2611}

    def __init__(self, label, exists=False, optional=True,
                 Qdot_min=.1, Qdot_max=3.2):
        # NOTE: Since the CHP module hase two outputs (heat and power) we use
        #       one of them (heat) as the base commodity for which the part
        #       load restrictions are imposed and then add the other one, along
        #       with the related constraints and connector.
        #       Also, since the base efficiency is assumed to depend on the
        #       nominal size, we remove the standard input output relation and
        #       create a custom one.
        super().__init__(label, self.Qdot_ref, self.c_ref, self.c_m, M=self.M,
                         nom_min=Qdot_min, nom_max=Qdot_max,
                         min_part_load=self.qdot_min, in_name='Qdot_in',
                         out_name='Qdot_out', out_connector_name='HEAT_OUT',
                         exists=exists, optional=optional)

        # size-dependent base efficiencies for custom input output relation
        # approximate the base efficiency by using the mean value
        # [-], nominal thermal eff.
        eff_th_nom = 0.498 - 3.55e-2 * (Qdot_min + Qdot_max)/2
        self.add_expression('eff_th_nom', eff_th_nom)
        # [-], nominal electrical eff.
        eff_el_nom = 0.372 + 3.55e-2 * (Qdot_min+Qdot_max)/2
        self.add_expression('eff_el_nom', eff_el_nom)

        # operation-dependent efficiencies (depends on q and Qdot_nom)
        out_rel = self['Qdot_out_rel']
        eff_th = eff_th_nom * sum(coeff * out_rel ** power for power, coeff
                                  in self.fit_params_th.items())
        self.add_expression('eff_th', eff_th)
        eff_el = eff_el_nom * sum(coeff * out_rel ** power for power, coeff
                                  in self.fit_params_el.items())
        self.add_expression('eff_el', eff_el)

        inp = self.get_expression('input')
        inp.ub = Qdot_max/eff_th_nom
        out_nom = self['Qdot_out_nom']
        self.add_le_constraint(inp, out_nom / eff_th_nom, name='input_limit')
        # update input_output_relation constraint with new efficiency
        out = self.get_expression('output')
        self.add_eq_constraint(inp * eff_th, out, 'input_output_relation')
        # add second output
        out_pow = self.make_operational_variable('power_output',
                                                 bounds=(0, (Qdot_max
                                                             / eff_th_nom
                                                             * eff_el_nom)))
        self.add_eq_constraint(inp * eff_el, out_pow, 'input_output_power')

        self.add_expression('power_output', out_pow)
        self.add_output('POWER_OUT', out_pow)


class CombinedHeatAndPower(System):
    """Combined-heat-and-power (CHP) engine parameterized as in Sass2020.

    design variable: nominal heating power to be installed
    operational variable: operational output (heat and elec.) for each timestep
    The investment costs 'c_inv' is a nonlinear function of the nominal thermal
    output 'Qdot_nom':
    c_inv = c_inv,ref * (Qdot_nom/Qdot_ref)^M
    The input 'Qdot_in' can be determined via the efficiency relation:
    Qdot_in =  Qdot_out / eff_th
    The thermal efficiency 'eff_th' is determined via the product of nominal
    thermal efficiency and a polynomial fitting function f(q), with the load
    fraction q = Qdot_out/Qdot_nom, giving the relationship:
    eff_th = f(q) * eff_th_nom
    The electric efficiency 'eff_el' is determined via the product of nominal
    electrical efficiency and f(q), giving the relationship:
    eff_el = f(q) * eff_el_nom
    Thus the output power can be computed as:
    P_out = Qdot_in * eff_el = Qdot_out / eff_th * eff_el
    """

    def __init__(self, label, choose=True):
        CHP_small = CHP_sub(label+'__small', Qdot_min=.1, Qdot_max=1.4)
        CHP_medium = CHP_sub(label+'__medium', Qdot_min=1.4, Qdot_max=2.3)
        CHP_large = CHP_sub(label+'__large', Qdot_min=2.3, Qdot_max=3.2)

        comps = [CHP_small, CHP_medium, CHP_large]
        conns = {conn: [getattr(comp, conn) for comp in comps]
                 for conn in ['IN', 'HEAT_OUT', 'POWER_OUT']}
        super().__init__(label, comps, conns)

        self.add_expression('Qdot_out_nom', sum(c['Qdot_out_nom']
                            for c in comps))
        # Join connectors of the individual components
        for conn in conns:
            self.extend_connection(conn)

        # Store aggregate expressions for the overall CHP
        for expr_name in ['input', 'output']:
            expr = self.aggregate_component_expressions(expr_name)
            self.add_expression(expr_name, expr)
        # NOTE: For expressions that are aggregated again at a later stage
        #       storing them under the same name as in the subcomponents would
        #       result in doubling!
        for expr_name in ['investment_costs', 'fixed_costs']:
            expr = self.aggregate_component_expressions(expr_name)
            self.add_expression(f'aggregate_{expr_name}', expr)

        if choose:
            # ensure only one of the submodels exists
            self.add_le_constraint(CHP_small['exists'] + CHP_medium['exists']
                                   + CHP_large['exists'], 1, name='selection')


class PhotovoltaicModule(IESComponent):
    """Photovoltaic modules (PV) parameterized as in Sass2019.

    design variable: nominal power to be installed
    operational variable: operational output power for each time step
    The investment costs 'c_inv' is a nonlinear function of the nominal output
    'Pel_nom':
    c_inv = c_inv,ref * (Pel_nom/Pel_ref)^M
    The output 'P_el' is limited by the solar irradiance ´i_solar´ and the PV
    output capacity ´Pel_nom´:
    P_el <= eff * Pel_nom/pel_area * i_solar
    """

    Pel_ref = .001  # [MW], reference nominal power
    c_ref = 4254.3/1000  # [k€], reference cost
    c_m = 0.01  # [-], maintenance coefficient, (fraction of investment cost)
    M = 0.9502  # [-], cost exponent
    Pel_min = 5/1000  # [MW], minimal nominal power allowed for the model
    Pel_max = 550/1000  # [MW], maximal nominal power allowed for the model
    pel_min = 0.0  # [-] minimum output part load
    pel_area = 0.171 / 1000  # [MW/m^2] nominal capacity per square meter

    # maintenance coefficient, fraction of the investment cost
    def __init__(self, label, exists=False, optional=True):
        super().__init__(label, self.Pel_ref, self.c_ref, self.c_m, M=self.M,
                         nom_min=self.Pel_min, nom_max=self.Pel_max,
                         min_part_load=self.pel_min,
                         out_name='Pel_out', exists=exists, optional=optional)

        # remove input connection
        del(self.connectors['IN'])
        del self._expressions_dict['input']
        del self._operational_variables_dict['input']

        # overwrite the input-output equality with an inequality
        eff = 0.19  # [-], efficiency
        i_solar = self.make_parameter('i_solar')  # [MW/m^2] irradiance
        self.add_ge_constraint(1/self.pel_area*eff*self['Pel_out_nom']*i_solar,
                               self['output'],  name='input_output_relation')


class Storage(Component):
    """General storage unit.

    design variable: nominal capacity to be installed
    operational variable: charge/discharge rate in operation
    The investment costs 'c_inv' is a nonlinear function of the capacity
    'C_nom':
    c_inv = c_inv,ref * C_nom^M
    The storage state of charge is modeled using state variables, for which
    time coupling constraints are created after the parametrization of the
    problem.

    Arguments
    ---------
    - label : str
        Unique sting that serves as an identifier of the storage unit.
    - nom_ref: The nominal reference size of the storage unit
    - c_ref: coefficient of investment costs
    - c_m: factor of investment costs corresponding to maintenance costs
          TODO: This implies calculations with a year as the time series!
    - M: (1) exponent of investment costs
    - min_cap: (0) lower bound for nominal size
    - max_cap: (None) upper bound for nominal size
    - charge_eff: (1) the charging efficiency of the storage unit
    - discharge_eff: (1) the discharging efficiency of the storage unit
    - c_loss: (0) the relative self-discharge of the unit (% SOC/h)
    - in_min: (0) minimum charging rate relative to nominal capacity
    - in_max: (1) maximum charging rate relative to nominal capacity
    - out_min: (0) minimum discharging rate relative to nominal capacity
    - out_max: (1) mxaimum discharging rate relative to nominal capacity

    - in_name: name of input commodity
    - out_name: name of output commodity
    - in_connector_name: name of input connector
    - out_connector_name: name of output connector
    - exists: bool
        Specifies whether the Component is currently installed; if set to
        True, its nominal size will be a user-specified parameter. If set
        to False, the nominal size of the component is a design variable.

    """

    def __init__(self, label, nom_ref, c_ref, c_m, M=1, min_cap=0,
                 max_cap=None, charge_eff=1, discharge_eff=1, c_loss=0,
                 in_min=0, in_max=1, out_min=0, out_max=1, in_name='input',
                 out_name='output', in_connector_name='IN',
                 out_connector_name='OUT', exists=False):
        super().__init__(label)

        # nominal (installed) output power
        nom_name = 'capacity_nom'
        if exists:
            cap = self.make_parameter(nom_name)
            self.add_le_constraint(min_cap, cap, name=nom_name + '_min')
            self.add_le_constraint(cap, max_cap, name=nom_name + '_max')
            # investment cost of the component
            inv = c_ref * (cap / nom_ref) ** M
        else:
            cap = self.make_design_variable(nom_name, bounds=(0, max_cap),
                                            init_val=max_cap)
            exists = self.make_design_variable('exists', BINARY, init_val=1)
            # bound component size
            self.add_le_constraint(exists * min_cap, cap,
                                   name=nom_name + '_min')
            self.add_le_constraint(cap, exists * max_cap,
                                   name=nom_name + '_max')
            # investment cost of the component
            inv = c_ref * (cap / nom_ref) ** M
            # set (investment) costs
            self.add_expression('investment_costs', inv)
        # add maintenance costs to fixed costs
        self.add_expression('fixed_costs', c_m * inv)

        # define operational variables
        soc = self.make_operational_variable('soc', bounds=(0, max_cap))
        self.add_le_constraint(soc, exists * cap, name=nom_name + '_max')

        # input at timestep t
        inp = self.make_operational_variable(in_name, bounds=(0, max_cap))
        # binary indicating whether the storage is charging
        b_in = self.make_operational_variable('b_in', domain=BINARY,
                                              init_val=0)
        # the upper/lower bound for input depends on the nominal output, but
        # since that's a variable  we have to write this bound as a constraint:
        self.add_le_constraint(inp, b_in * cap * in_max, name='input_max')
        self.add_ge_constraint(inp, b_in * cap * in_min, name='input_min')
        # output at timestep t
        out = self.make_operational_variable(out_name, bounds=(0, max_cap))
        # binary indicating whether the storage is discharging
        b_out = self.make_operational_variable('b_out', domain=BINARY,
                                               init_val=1)
        # the upper/lower bound for output depends on the nominal output, but
        # since that's a variable  we have to write this bound as a constraint:
        self.add_le_constraint(out, b_out * cap * out_max, name='output_max')
        self.add_ge_constraint(out, b_out * cap * out_min, name='output_min')

        dissipation = self.add_expression('dissipation',  soc * c_loss)
        # define storage constraints and declare state variables
        state_change = inp * charge_eff - out / discharge_eff \
            - dissipation
        self.declare_state(soc, state_change, 0)

        # make sure, storage units cannot charge and discharge at the same time
        self.add_le_constraint(b_in + b_out, 1, name='charge_discharge_switch')

        self.add_expression('input', inp)
        self.add_expression('output', out)

        # set connectors
        self.add_input(in_connector_name, inp)
        self.add_output(out_connector_name, out)


class Battery(Storage):
    """Battery unit (BAT) parametrized as in Sass2020.

    design variable: nominal capacity to be installed
    operational variable: charge/discharge rate in operation
    The investment costs 'c_inv' is a nonlinear function of the capacity
    'C_nom':
    c_inv = c_inv,ref * C_nom^M
    The storage volume is modeled using state variables, for which time
    coupling constraints are created post-parametrization.
    """

    # cost parameters
    nom_ref = .001  # [MWh], reference nominal capacity
    c_ref = 2116.1/1000  # [k€], reference cost
    M = 0.8382  # [-], cost exponent
    c_m = 0.025  # [-], maintenance coefficient, (fraction of investment cost)
    # bounds of component size
    min_cap = 0  # [MW], minimal capacity allowed for the model
    max_cap = 2000/1000  # [MW], maximum capacity allowed for the model
    # efficiencies
    charge_eff = 0.920  # [-], charging efficiency
    discharge_eff = 0.926  # [-], discharging efficiency
    c_loss = 0  # [1/h], self-discharge loss
    # operational parameters
    in_min = 0  # [1/h], minimum charging rate relative to nominal power
    in_max = 0.36  # [1/h], maximum charging rate relative to nominal power
    out_min = 0  # [1/h], minimum discharging rate relative to nominal power
    out_max = 0.36  # [1/h], maximum discharging rate relative to nominal power

    def __init__(self, label, exists=False):
        super().__init__(label, self.nom_ref, self.c_ref, self.c_m,
                         self.M, self.min_cap, self.max_cap, self.charge_eff,
                         self.discharge_eff, self.c_loss, self.in_min,
                         self.in_max, self.out_min, self.out_max,
                         in_name='input', out_name='output',
                         in_connector_name='IN', out_connector_name='OUT',
                         exists=exists)


class CoolWaterStorage(Storage):
    """Cooling storage unit (STC) parametrized as in Sass2020.

    design variable: nominal capacity to be installed
    operational variable: charge/discharge rate in operation
    The investment costs 'c_inv' is a nonlinear function of the capacity
    'C_nom':
    c_inv = c_inv,ref * C_nom^M
    The storage volume is modeled using state variables, for which time
    coupling constraints are created post-parametrization.
    """

    # cost parameters
    nom_ref = .001  # [MWh], reference nominal capacity
    c_ref = 57.5/1000  # [k€], reference cost
    M = 0.9037  # [-], cost exponent
    c_m = 0.01  # [-], maintenance coefficient, (fraction of investment cost)
    # bounds of component size
    min_cap = 0  # [MW], minimal capacity allowed for the model
    max_cap = 25000/1000  # [MW], maximum capacity allowed for the model
    # efficiencies
    charge_eff = 0.95  # [-], charging efficiency
    discharge_eff = 0.95  # [-], discharging efficiency
    c_loss = 0.005  # [1/h], self-discharge loss
    # operational parameters
    in_min = 0  # [1/h], minimum charging rate relative to nominal power
    in_max = 1  # [1/h], maximum charging rate relative to nominal power
    out_min = 0  # [1/h], minimum discharging rate relative to nominal power
    out_max = 1  # [1/h], maximum discharging rate relative to nominal power

    def __init__(self, label, exists=False):
        super().__init__(label, self.nom_ref, self.c_ref, self.c_m,
                         self.M, self.min_cap, self.max_cap, self.charge_eff,
                         self.discharge_eff, self.c_loss, self.in_min,
                         self.in_max, self.out_min, self.out_max,
                         in_name='input', out_name='output',
                         in_connector_name='IN', out_connector_name='OUT',
                         exists=exists)


class HotWaterStorage(Storage):
    """Heating storage unit (STH) parametrized as in Sass2020.

    design variable: nominal capacity to be installed
    operational variable: charge/discharge rate in operation
    The investment costs 'c_inv' is a nonlinear function of the capacity
    'C_nom':
    c_inv = c_inv,ref * C_nom^M
    The storage volume is modeled using state variables, for which time
    coupling constraints are created post-parametrization.
    """

    # cost parameters
    nom_ref = .001  # [MWh], reference nominal capacity
    c_ref = 83.8/1000  # [k€], reference cost
    M = 0.8663  # [-], cost exponent
    c_m = 0.01  # [-], maintenance coefficient, (fraction of investment cost)
    # bounds of component size
    min_cap = 0  # [MW], minimal capacity allowed for the model
    max_cap = 115000/1000  # [MW], maximum capacity allowed for the model
    # efficiencies
    charge_eff = 0.95  # [-], charging efficiency
    discharge_eff = 0.95  # [-], discharging efficiency
    c_loss = 0.005  # [1/h], self-discharge loss
    # operational parameters
    in_min = 0  # [1/h], minimum charging rate relative to nominal power
    in_max = 1  # [1/h], maximum charging rate relative to nominal power
    out_min = 0  # [1/h], minimum discharging rate relative to nominal power
    out_max = 1  # [1/h], maximum discharging rate relative to nominal power

    def __init__(self, label, exists=False):
        super().__init__(label, self.nom_ref, self.c_ref, self.c_m,
                         self.M, self.min_cap, self.max_cap, self.charge_eff,
                         self.discharge_eff, self.c_loss, self.in_min,
                         self.in_max, self.out_min, self.out_max,
                         in_name='input', out_name='output',
                         in_connector_name='IN', out_connector_name='OUT',
                         exists=exists)
