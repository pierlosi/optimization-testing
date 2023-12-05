# %%

"""
This is the implementation of the pseudo code. The dispatcher has been implemented as a function.
To be able to run the dispatcher I used the electricity price profile of zone SE3 in Sweden in year 2021 and a fictious fixed hydrogen price.
The project design technical variables are the same as the one used in the file provided in the Case Challenge.
"""

import pyomo.environ as pyo
from pyomo.opt import SolverFactory

from math import sqrt
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# %% Here I define the dispatcher function


def dispatcher(design):
    time_horizon = len(design["P_PV_data"])

    # Transoforming dataframes into dictionaries so that pyomo can read them
    P_PV_dict = dict()
    for i in range(time_horizon):
        P_PV_dict[i] = (
            design["P_PV_data"]["P_PV_kW"].values[i] * design["solar_pv_size"]
        )

    electricity_price_dict = dict()
    for i in range(time_horizon):
        electricity_price_dict[i] = (
            design["Electricity_price_data"]["Grid_Price"].values[i] / 100
        )  # [€/kWh]

    # INITIALIZE the optimization framework
    m = pyo.ConcreteModel()

    # INITIALIZATION of SETS
    # m.iIDX is the time index, which set the optimization horizon of the problem in an hourly fashion
    m.iIDX = pyo.Set(initialize=range(time_horizon))

    # INITIALIZATION of PARAMETERS

    """importing data in the pyomo framewrok"""
    m.P_PV = pyo.Param(m.iIDX, initialize=P_PV_dict)  # [kW]
    m.el_price = pyo.Param(m.iIDX, initialize=electricity_price_dict)  # [€/kWh]
    m.h2_price = pyo.Param(initialize=design["Hydrogen_price_data"])  # [€/kg]

    m.m_h2_demand = pyo.Param(initialize=design["Hydrogen_demand"])  # [kg/h]

    m.SEC_elec = pyo.Param(initialize=design["SEC_nom"])
    m.SEC_compr = pyo.Param(initialize=design["comp_SEC_nom"])

    # here I compute the upper and lower bound for the electrolyzer power
    P_elec_min = design["Pmin"] / 100 * design["elec_size"] * 1e3  # [kW]
    P_elec_max = design["Pmax"] / 100 * design["elec_size"] * 1e3  # [kW]
    m.P_elec_min = pyo.Param(initialize=P_elec_min)
    m.P_elec_max = pyo.Param(initialize=P_elec_max)

    # INITIALIZATION of VARIABLES

    m.P_grid = pyo.Var(
        m.iIDX, domain=pyo.NonNegativeReals
    )  # [kW] power from the PV to the grid
    m.P_PV_to_elec = pyo.Var(
        m.iIDX, domain=pyo.NonNegativeReals
    )  # [kW] power from the PV to the electrolyzer
    m.P_curt = pyo.Var(
        m.iIDX, domain=pyo.NonNegativeReals
    )  # [kW] power from the PV to the curtailment
    m.P_elec = pyo.Var(
        m.iIDX, domain=pyo.NonNegativeReals
    )  # [kW] power from the PV to the electrolyzer
    m.bin_elec = pyo.Var(
        m.iIDX, domain=pyo.Binary
    )  # state of the electrolyzer -> binary variable -> electrolyzer either on or off

    m.m_elec = pyo.Var(
        m.iIDX, domain=pyo.NonNegativeReals
    )  # [kg/h] hydrogen mass flow out of the electrolyzer

    # DEFINITION of the optimization problem

    # definition of the objective function
    def obj_funct(m):
        return sum(
            (-m.el_price[i] * m.P_grid[i] + m.h2_price * m.m_elec[i]) for i in m.iIDX
        )

    m.obj = pyo.Objective(rule=obj_funct, sense=pyo.maximize)

    # power equilibrium constraint
    def f_power_equilibrium(m, i):
        return m.P_PV_to_elec[i] + m.P_grid[i] == m.P_elec[i]

    m.cstr_power_equilibrium = pyo.Constraint(m.iIDX, rule=f_power_equilibrium)

    # PV power constraint
    def f_PV_power(m, i):
        return m.P_PV_to_elec[i] + m.P_curt[i] == m.P_PV[i]

    m.cstr_PV_power = pyo.Constraint(m.iIDX, rule=f_PV_power)

    # here I implement the constraints so that when the electolyzer is on it can only work between P_min and P_max
    # the binary commitment variable is necessary to leave the possibility to the electrolyzer to be switched off
    # (bin_elec = 1 when electrolyzer is on and bin_elec = 0 when it is off)

    def f_elec_power_commit_sup(m, i):
        return m.P_elec[i] <= m.bin_elec[i] * m.P_elec_max

    m.cstr_elec_power_commit_sup = pyo.Constraint(m.iIDX, rule=f_elec_power_commit_sup)

    def f_elec_power_commit_inf(m, i):
        return m.P_elec[i] >= m.bin_elec[i] * m.P_elec_min

    m.cstr_elec_power_commit_inf = pyo.Constraint(m.iIDX, rule=f_elec_power_commit_inf)

    # electrolyzer hydrogen production constraint

    def f_hydrogen_production(m, i):
        return m.m_elec[i] == m.P_elec[i] / m.SEC_elec

    m.cstr_hydrogen_production = pyo.Constraint(m.iIDX, rule=f_hydrogen_production)

    def f_demand_satisfaction(m, i):
        return m.m_elec[i] == m.m_h2_demand

    m.cstr_demand_satisfaction = pyo.Constraint(m.iIDX, rule=f_demand_satisfaction)

    # selection of the optimization solver (minding that it is suited for the kind of problem)
    opt = pyo.SolverFactory("cbc")

    # resolution of the problem
    opt.solve(m)

    # Storing the data in arrays that can later be exported (maybe in a dataframe) and/or be displayed
    P_elec = np.array([pyo.value(m.P_elec[i]) for i in m.iIDX])
    P_grid = np.array([pyo.value(m.P_grid[i]) for i in m.iIDX])
    P_PV = np.array([pyo.value(m.P_PV[i]) for i in m.iIDX])
    el_price = np.array([pyo.value(m.el_price[i]) for i in m.iIDX])
    m_elec = np.array([pyo.value(m.m_elec[i]) for i in m.iIDX])

    data_time = pd.DataFrame(
        {
            "Time": design["P_PV_data"]["Time"],
            "P_elec": P_elec,
            "P_grid": P_grid,
            "P_PV": P_PV,
            "el_price": el_price,
            "m_elec": m_elec,
        }
    )

    return data_time


# %% READING input data

P_PV_data = pd.read_csv(
    "PVdata_1MWp_1year_unit_kW.csv", names=["P_PV_kW"], nrows=8760
)  # [kW] ([kWh/h])

# Creating a datetime type column
start_date = "2022-01-01 00:00:00"
start_datetime = datetime.datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
datetime_array = [
    start_datetime + datetime.timedelta(hours=i) for i in range(len(P_PV_data))
]

P_PV_data["Time"] = datetime_array

# The electricity data refers to Sweden zone SE3 for year 2021
electricity_price_data = pd.read_csv(
    "PriceCurve_SE3_2021.csv", nrows=8760, header=0, sep=";"
)  # [cents/kWh]

# Random price chosen for hydrogen sale
hydrogen_price_data = 10  # [€/kg]

""" 
It is preferrable to sell electricity when its price is higher than 
the price of hydrogen divided by the specific energy consumption of the electrolyzer

Hence, supposing that the price of hydrogen is 20 €/kg and the SEC_elec is 40 kWh/kg,
electricity will be sold to the grid if the price of electricity will be higher than 20/40 = 0.5 €/kWh,
otherwise it will be used to produce hydrogen and sell it

"""


# Here I define what is the design of the project. I also add the dataframes with PV power production, electricity prices and hydrogen prices
design = {
    "elec_size": 1,  # MW
    "solar_pv_size": 1.8,  # MW
    "Pmax": 100,  # %
    "Pmin": 10,  # %
    "horizon": 1,  # year
    "deg_power": 0.5,  # %/year
    "SEC_nom": 45,  # kWh/kgH2       #Specific Energy Consumption of the electrolyzer system as a whole
    "comp_size": 1,  # MW
    "comp_to_pressure": 200,  # bar
    "comp_from_pressure": 20,  # bar
    "comp_SEC_nom": 0.1,  # kWh/kgH2
    "P_PV_data": P_PV_data,  # [kW]
    "Electricity_price_data": electricity_price_data,  # [cents/kWh]
    "Hydrogen_price_data": hydrogen_price_data,  # [€/kg]
    "Hydrogen_demand": 15,
}

# %% I run the simulation

data_time = dispatcher(design)

# %%
day_start_display = 0  # which day of the year do you want to start? cardinal number
days_display = 5  # how many days do you want to display?

fig, ax_pow = plt.subplots()
fig.suptitle("Power")

ax_pow.plot(
    data_time["P_elec"][
        24 * day_start_display : 24 * (day_start_display + days_display)
    ],
    color="green",
    label="P_elec",
)
# ax_pow.plot(P_compr[24*day_start_display:24*(day_start_display+days_display)], color = 'blue', label = 'P_compr')
ax_pow.plot(
    data_time["P_grid"][
        24 * day_start_display : 24 * (day_start_display + days_display)
    ],
    color="red",
    label="P_grid",
)
# ax_pow.plot(P_PV[24*day_start_display:24*(day_start_display+days_display)], color = 'black', label = 'P_PV')
ax_pow.plot(
    data_time["P_PV"][24 * day_start_display : 24 * (day_start_display + days_display)],
    color="purple",
    label="P_PV",
)  # [kW]
ax_pow.set_ylabel("Power [kW]")
ax_pow.set_xlabel("Time [h]")

ax_pow.legend(loc="best", bbox_to_anchor=(1.3, 1))

ax_price = ax_pow.twinx()

ax_price.plot(
    data_time["el_price"][
        24 * day_start_display : 24 * (day_start_display + days_display)
    ],
    color="orange",
    label="El_price",
)
ax_price.set_ylabel("Price [€/kWh]")
ax_price.legend(loc="best", bbox_to_anchor=(1.3, 0.8))

plt.grid(True)
plt.show()
plt.close()
# %%
