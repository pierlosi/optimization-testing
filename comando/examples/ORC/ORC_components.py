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
import comando


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

    def __init__(self, label, dT_min, **kwargs):
        super().__init__(label)

        exprs = self.expressions_dict

        def alldefined(*args):
            """Check if all arguments are specified, if yes store them."""
            if all(arg in kwargs for arg in args):
                for arg in args:
                    exprs[arg] = kwargs[arg]
                return True

        def T_based_Hdot(x):
            """Set hot/cold enthalpy flows based on T and return dT."""
            exprs[f'Hdot_{x}_in'] = exprs[f'mdot_cp_{x}'] * exprs[f'T_{x}_in']
            exprs[f'Hdot_{x}_out'] = exprs[f'mdot_cp_{x}'] \
                * exprs[f'T_{x}_out']
            if x == 'h':
                exprs['dT_h'] = dT = exprs['T_h_in'] - exprs['T_h_out']
            else:
                exprs['dT_c'] = dT = exprs['T_c_out'] - exprs['T_c_in']
            return dT

        def h_based_Hdot(x):
            """Set hot/cold enthalpy flows based on h and return dh."""
            if not alldefined(f'T_{x}_in', f'T_{x}_out'):
                raise RuntimeError('When defining the hot side via enthalpies')
            exprs[f'Hdot_{x}_in'] = exprs[f'mdot_{x}'] * exprs[f'h_{x}_in']
            exprs[f'Hdot_{x}_out'] = exprs[f'mdot_{x}'] * exprs[f'h_{x}_out']
            if x == 'h':
                exprs['dh_h'] = dh = exprs['h_h_in'] - exprs['h_h_out']
                exprs['dT_h'] = exprs['T_h_in'] - exprs['T_h_out']
            else:
                exprs['dh_c'] = dh = exprs['h_c_out'] - exprs['h_c_in']
                exprs['dT_c'] = exprs['T_c_out'] - exprs['T_c_in']
            self.add_ge_constraint(dh/1e5, 0, f'dh_{x} >= 0')
            return dh

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
            Qdot_c = coolQ('h_h_in')
            exprs['h_h_in'] = exprs['h_h_out'] + Qdot_c/exprs['mdot_h']
            h_based_Hdot('h')
        elif 'h_h_out' not in kwargs and alldefined('h_h_in', 'mdot_h'):
            Qdot_c = coolQ('h_h_out')
            exprs['h_h_out'] = exprs['h_h_in'] - Qdot_c/exprs['mdot_h']
            h_based_Hdot('h')
        elif 'mdot_h' not in kwargs and alldefined('h_h_in', 'h_h_out'):
            Qdot_c = coolQ('mdot_h')
            exprs['mdot_h'] = Qdot_c / (exprs['h_h_in'] - exprs['h_h_out'])
            h_based_Hdot('h')
        elif 'T_h_in' not in kwargs and alldefined('T_h_out', 'mdot_cp_h'):
            Qdot_c = coolQ('T_h_in')
            exprs['T_h_in'] = exprs['T_h_out'] + Qdot_c/exprs['mdot_cp_h']
            T_based_Hdot('h')
        elif 'T_h_out' not in kwargs and alldefined('T_h_in', 'mdot_cp_h'):
            Qdot_c = coolQ('T_h_out')
            exprs['T_h_out'] = exprs['T_h_in'] - Qdot_c/exprs['mdot_cp_h']
            T_based_Hdot('h')
        elif all(k not in kwargs for k in ('mdot_cp_h', 'h_h_in', 'h_h_out')) \
                and alldefined('T_h_in', 'T_h_out'):
            Qdot_c = coolQ('mdot_cp_h')
            exprs['mdot_cp_h'] = Qdot_c / (exprs['T_h_in'] - exprs['T_h_out'])
            T_based_Hdot('h')
        # Now the cold side...
        elif 'h_c_in' not in kwargs and alldefined('h_c_out', 'mdot_c'):
            Qdot_h = hotQ('h_c_in')
            exprs['h_c_in'] = exprs['h_c_out'] - Qdot_h/exprs['mdot_c']
            h_based_Hdot('c')
        elif 'h_c_out' not in kwargs and alldefined('h_c_in', 'mdot_c'):
            Qdot_h = hotQ('h_c_out')
            exprs['h_c_out'] = exprs['h_c_in'] + Qdot_h/exprs['mdot_c']
            h_based_Hdot('c')
        elif 'mdot_c' not in kwargs and alldefined('h_c_in', 'h_c_out'):
            Qdot_h = hotQ('mdot_c')
            exprs['mdot_c'] = Qdot_h / (exprs['h_c_out'] - exprs['h_c_in'])
            h_based_Hdot('c')
        elif 'T_c_in' not in kwargs and alldefined('T_c_out', 'mdot_cp_c'):
            Qdot_h = hotQ('T_c_in')
            exprs['T_c_in'] = exprs['T_c_out'] - Qdot_h/exprs['mdot_cp_c']
            T_based_Hdot('c')
        elif 'T_c_out' not in kwargs and alldefined('T_c_in', 'mdot_cp_c'):
            Qdot_h = hotQ('T_c_out')
            exprs['T_c_out'] = exprs['T_c_in'] + Qdot_h/exprs['mdot_cp_c']
            T_based_Hdot('c')
        elif all(k not in kwargs for k in ('mdot_cp_c', 'h_c_in', 'h_c_out')) \
                and alldefined('T_c_in', 'T_c_out'):
            Qdot_h = hotQ('mdot_cp_c')
            exprs['mdot_cp_c'] = Qdot_h / (exprs['T_c_out'] - exprs['T_c_in'])
            T_based_Hdot('c')
        else:  # Assume all quantities are defined
            Qdot_c = coolQ('h_h_out')
            Qdot_h = hotQ('mdot_cp_c')
            self.add_eq_constraint((Qdot_c - Qdot_h) * 1e-6, 0, 'heat_flow')

        # Deltas for temperatures
        self.add_ge_constraint(exprs['dT_h']/100, 0, 'dT_h >= 0')
        self.add_ge_constraint(exprs['dT_c']/100, 0, 'dT_c >= 0')
        exprs['dT_h_in'] = dT_h_in = exprs['T_h_in'] - exprs['T_c_out']
        self.add_ge_constraint(dT_h_in/100, dT_min/100, f'dT_h_in >= {dT_min}')
        exprs['dT_h_out'] = dT_h_out = exprs['T_h_out'] - exprs['T_c_in']
        self.add_ge_constraint(dT_h_out/100, dT_min/100,
                               f'dT_h_out >= {dT_min}')

        # Connectors
        self.add_input('Hdot_h_in', exprs['Hdot_h_in'])
        self.add_output('Hdot_h_out', exprs['Hdot_h_out'])
        self.add_input('Hdot_c_in', exprs['Hdot_c_in'])
        self.add_output('Hdot_c_out', exprs['Hdot_c_out'])

        for id, expr in exprs.items():
            self.add_expression(id, expr)


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

    h_out : expression
        The specific enthalpy at the output of the pump or the specific,
        isentropic enthalpy if eta_s is also given

    eta_s : expression
        The isentropic efficiency of the pump

    The mass flow, pressures and specific enthalpies are optional; if they are
    not given, variables are created for the corresponding quantities.
    """

    def __init__(self, label, mdot=None, p_in=None, p_out=None, h_in=None,
                 h_out=None, eta_s=None):
        super().__init__(label)
        exprs = self.expressions_dict
        exprs.update({'mdot': mdot, 'p_in': p_in, 'p_out': p_out,
                      'h_in': h_in, 'h_out': h_out})
        # Creation of variables for missing mass flow, pressures and enthalpies
        for e_name, e in exprs.items():
            # TODO: is lower bound of 0 for enthalpy problematic?
            exprs[e_name] = e if e is not None else \
                self.make_operational_variable(e_name, bounds=(0, None))

        if eta_s is None:  # h_out is the actual specific output enthalpy
            exprs['w'] = w = exprs['h_out'] - exprs['h_in']
        else:  # h_out is the isentropic specific output enthalpy
            exprs['h_out_is'] = exprs['h_out']
            exprs['w'] = w = (exprs['h_out_is'] - exprs['h_in']) / eta_s
            exprs['h_out'] = exprs['h_in'] + w

        exprs['P'] = P = exprs['mdot'] * w
        self.add_ge_constraint(P * 1e-5, 0, 'P >= 0')

        self.add_input('Hdot_in', exprs['mdot'] * exprs['h_in'])
        self.add_output('Hdot_out', exprs['mdot'] * exprs['h_out'])

        for id, expr in exprs.items():
            self.add_expression(id, expr)


class Turbine(comando.Component):
    """A simple model for a vapor turbine.

    Arguments
    ---------
    label : `str`
        A unique sting that serves as an identifier

    mdot: expression
        The mass flow passing through the turbine

    p_in : expression
        The pressure at the input of the turbine

    p_out : expression
        The pressure at the output of the turbine

    h_in : expression
        The specific enthalpy at the input of the turbine

    h_out : expression
        The specific enthalpy at the output of the turbine or the specific,
        isentropic enthalpy if eta_s is also given

    eta_s : expression
        The isentropic efficiency of the turbine

    The mass flow, pressures and specific enthalpies are optional; if they are
    not given, variables are created for the corresponding quantities.
    """

    def __init__(self, label, mdot=None, p_in=None, p_out=None, h_in=None,
                 h_out=None, eta_s=None):
        super().__init__(label)

        exprs = self.expressions_dict
        exprs.update({'mdot': mdot, 'p_in': p_in, 'p_out': p_out,
                      'h_in': h_in, 'h_out': h_out})
        # Creation of variables for missing mass flow, pressures and enthalpies
        for e_name, e in exprs.items():
            # TODO: is lower bound of 0 for enthalpy problematic?
            exprs[e_name] = e if e is not None else \
                self.make_operational_variable(e_name, bounds=(0, None))

        if eta_s is None:  # h_out is the actual specific output enthalpy
            exprs['w'] = w = exprs['h_in'] - exprs['h_out']
        else:  # h_out is the isentropic specific output enthalpy
            exprs['w'] = w = (exprs['h_in'] - exprs['h_out']) * eta_s
            exprs['h_out'] = exprs['h_in'] - w

        exprs['P'] = P = exprs['mdot'] * w
        self.add_ge_constraint(P * 1e-7, 0, 'P >= 0')

        self.add_input('Hdot_in', exprs['mdot'] * exprs['h_in'])
        self.add_output('Hdot_out', exprs['mdot'] * exprs['h_out'])

        for id, expr in exprs.items():
            self.add_expression(id, expr)


class CoolingSystem(comando.Component):
    """A simple model for a fan-based cooling system.

    Arguments
    ---------
    label : `str`
        A unique sting that serves as an identifier

    mdot_cp_fluid: expression
        The mean heat capacity flow of the fluid being cooled

    T_fluid_in : expression
        The temperature at the fluid inlet

    T_fluid_out : expression
        The temperature at the fluid outlet

    Delta_p_f : expression
        The pressure increase over the fans

    rho_air : expression
        The density of air passing through the fans

    cp_air : expression
        The heat capacity of air passing through the fans

    The heat capacity flow and temperatures of the fluid are optional; if they
    are not given, variables are created for the corresponding quantities.
    """

    def __init__(self, label, mdot_cp_fluid=None, T_fluid_in=None,
                 T_fluid_out=None, Delta_p_f=170, rho_air=1.2,
                 cp_air=1e3, eta_f=0.65):
        super().__init__(label)

        exprs = self.expressions_dict
        exprs.update({'mdot_cp_fluid': mdot_cp_fluid, 'T_fluid_in': T_fluid_in,
                      'T_fluid_out': T_fluid_out})

        # Creation of variables for missing mass flow, pressures and enthalpies
        for e_name, e in exprs.items():
            # TODO: is lower bound of 0 for enthalpy problematic?
            exprs[e_name] = e if e is not None else \
                self.make_operational_variable(e_name, bounds=(0, None))

        mdot_cp_air = mdot_cp_fluid = exprs['mdot_cp_fluid']  # ASSUMPTION!
        exprs['mdot_air'] = mdot_air = mdot_cp_air / cp_air
        exprs['Vdot_air'] = Vdot_air = mdot_air / rho_air
        exprs['P'] = P = Vdot_air * Delta_p_f / eta_f
        # NOTE: Results in P = Q_17 \
        #           / (dT_cw_7) / cp_air / rho_air * Delta_p_f / eta_f
        # NOTE: The following constraint can be proven to always be satisfied
        #       using interval arithmetic
        self.add_ge_constraint(P * 1e-6, 0, 'P >= 0')
        exprs['Qdot_diss'] = mdot_cp_fluid * (exprs['T_fluid_in']
                                              - exprs['T_fluid_out'])

        self.add_input('Hdot_in', mdot_cp_fluid * exprs['T_fluid_in'])
        self.add_output('Hdot_out', mdot_cp_fluid * exprs['T_fluid_out'])

        for id, expr in exprs.items():
            self.add_expression(id, expr)
