"""Components for the building demand response case study."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Florian Joseph Baader, Marco Langiu
from comando.core import Component, INTEGER


class Mass(Component):
    """A mass with a temperature that varies due to heat fluxes."""

    def __init__(self, label, rho, C_p, V, T_0, T_bounds=False):
        super().__init__(label)

        # create heat flows
        # could be in and out
        Qdot_in = self.make_operational_variable('Qdot_in')
        # connect heat flow
        self.add_connectors(Q_IN=-Qdot_in)

        T, der_T = self.make_state('T', Qdot_in/rho/C_p/V, init_state=T_0)

        if T_bounds is True:
            T_min = self.make_parameter('T_min')
            T_max = self.make_parameter('T_max')

            self.add_ge_constraint(T, T_min)
            self.add_le_constraint(T, T_max)


class HeatTransport(Component):
    """A heat transport between component A and B."""

    # temperature of A and B is given in that oder
    # Heat transport is symmetric, just make sure that tempearture A and
    # connector A are connected to the same component
    def __init__(self, label, alpha, Area, T_A, T_B):
        super().__init__(label)

        # the heat flows to A and to B which can both be negative or positive

        # calculate heat flows
        Qdot_A = alpha * Area * (T_B - T_A)
        Qdot_B = - Qdot_A

        # 2 connectors to be connected to component A and B
        self.add_connectors(A=Qdot_A)
        self.add_connectors(B=Qdot_B)


class Grid(Component):
    """A simple electricity grid."""

    def __init__(self, label):
        super().__init__(label)
        # define objective, which is to minimize Pel*costs
        electricityPrice = self.make_parameter('ElectricityPrice')
        Pel_out = self.make_operational_variable('Pel', bounds=(0, None))
        self.add_expression('variable Mayer objective',
                            Pel_out * electricityPrice / 3600)
        self.add_output('OUT', Pel_out)


class AirHeatPump(Component):
    """An air heat pump with COP dependent on outside air temperature."""

    # give variable of outside air temperature
    def __init__(self, label, T_ambient, Q_max, n_outflows):
        super().__init__(label)
        self.Q_HP_max = Q_max  # 30 * 13 * 6  # maximum power of heat pump [W]
        self.zeta_HP = 0.5  # exergetic efficiency [-]
        # assumption that hot supply-temperature (deutsch: "Vorlauftempratur")
        # is constant
        self.T_HP_warm = 313.15  # warm temperature [K]
        # cold temperature is assumed to be 5 K below the ambient
        self.T_HP_cold = T_ambient - 5  # cold temperature [K]
        # define Q dot
        Q_HP = self.make_operational_variable('Q_HP',
                                              bounds=(0, self.Q_HP_max))
        # roundabout Qdot_HP_max / COP_max with COP_max = 4.5
        self.Pel_HP_max = self.Q_HP_max / 4.5
        # minimal partial load of heat pump
        self.Pel_HP_min = 0.2 * self.Pel_HP_max
        Pel_HP = self.make_operational_variable('Pel_HP',
                                                bounds=(0, self.Pel_HP_max))

        # indicate if heat pump is on or off
        b_on = self.make_operational_variable('b_on', domain=INTEGER,
                                              bounds=(0, 1))
        self.add_le_constraint(Pel_HP, b_on * self.Pel_HP_max)
        self.add_ge_constraint(Pel_HP, b_on * self.Pel_HP_min)

        # define relation between Q_HP and Pel dependent on COP
        self.add_eq_constraint(Q_HP, Pel_HP * self.zeta_HP
                               * (self.T_HP_warm
                                  / (self.T_HP_warm - self.T_HP_cold)))
        self.add_connectors('P_IN', Pel_HP)

        sumQdot_out = 0
        for i in range(1, n_outflows + 1):
            name = f'Qdot_out_{i}'
            Qdot_i = self.make_operational_variable(name, bounds=(0, None))
            sumQdot_out += Qdot_i
            self.add_connectors(id=name, expr=Qdot_i)

        self.add_eq_constraint(Q_HP, sumQdot_out)
