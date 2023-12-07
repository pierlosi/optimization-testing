# %%
import comando
import datetime
import numpy as np
import pandas as pd

from comando import Component


# %%


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
        demand = self.make_parameter("demand", data)
        self.add_input("IN", demand)


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
        resource = self.make_parameter("resource", data)
        self.add_output("OUT", resource)


class Sink(Component):
    """A sink for an arbitrary commodity."""

    def __init__(self, label):
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
        sink = self.make_operational_variable("sink", bounds=(0, None))
        self.add_input("IN", sink)


class Source(Component):
    """A source for an arbitrary commodity."""

    def __init__(self, label):
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
        source = self.make_operational_variable("source", bounds=(0, None))
        self.add_output("OUT", source)


class Electrolyzer(Component):
    """An electrolyzer for hydrogen production."""

    def __init__(self, label, sec, size):
        """Initialize the Electrolyzer.

        Arguments
        ---------
        - label : str
            Unique string that serves as an identifier of this Electrolyzer.
        - data : numeric data
            Amount of provided commodity, can be a scalar or a pandas
            Series.
        """
        super().__init__(label)

        sec = self.make_parameter("resource", sec)

        # Create operational variables
        P_in = self.make_operational_variable("P_in")

        # Expressions
        Q_out = self.add_expression("Q_out", P_in / sec)

        # Add constraints
        self.add_le_constraint(P_in, size, "P_max")

        # Add connectors
        self.add_input("IN", P_in)
        self.add_output("OUT", Q_out)


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

# %%
Component.existing_components = dict()
H2_CONS = Demand("H2_CONS", design["Hydrogen_demand"])
# PV = Resource("PV", design["P_PV_data"]["P_PV_kW"])
CURT = Sink("CURT")
GRID = Source("GRID")
PEM = Electrolyzer("PEM", design["SEC_nom"], design["elec_size"])


# %%
