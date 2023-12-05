# Simultaneous Optimization of Design and Operation of an Air-Cooled Geothermal ORC under Consideration of Multiple Operating Points

This directory contains code, data and results for the publication
```bibtex
@Article{langiu2021simultaneous,
  title   = {Simultaneous Optimization of Design and Operation of an Air-CooledGeothermal ORC under Consideration of Multiple Operating Points},
  journal = {Computers \& Chemical Engineering},
  volume  = {???},
  pages   = {??????},
  year    = {2021},
  issn    = {????-????},
  doi     = {https://doi.org/10.1016/j.compchemeng.2022.??????},
  author  = {Marco Langiu and Manuel Dahmen and Alexander Mitsos}
}
```
Available at [TODO: PREPRINT]()
[TODO: Journal](https://doi.org/10.1016/j.compchemeng.2021.??????).

Here we consider the design of an a air-cooled geothermal organic Rankine cycle (ORC).
We take into account the annual ambient air temperature distribution at the considered construction site, and thus the resulting operating points.
In this way, we are able to obtain a design that is robust in the sense that it allows for operation in all considered ambient conditions and optimal in the sense that its expected annual return is maximized.

## Code
The code implementing the considered system is split into `ORC_components.py`, where component models are implemented and `case_study.py`, where the system model is formulated based on the component models and considered temperatures, an optimization problem is formulated and solved, and the results are written.
The latter file can be executed with various command-line options, using any python version >= 3.8, for an overview see
```
python case_study.py --help
```

## ANN models
Code used for handling ANNs can be found in the `ANN` subdirectory.
The parameters for the ANNs used for calculating fluid properties of isobutane from
```bibtex
@Article{huster2019impact,
  author    = {Huster, Wolfgang R. and Schweidtmann, Artur M. and Mitsos, Alexander},
  pages     = {427--432},
  title     = {{Impact of Accurate Working Fluid Properties on the Globally Optimal Design of an Organic Rankine Cycle}},
  volume    = {47},
  doi       = {10.1016/b978-0-12-818597-1.50068-0},
  journal   = {Comput. Aided Chem. Eng.},
  publisher = {Elsevier},
  year      = {2019},
}
```
are encoded as json files in the `old_anns` subdirectory.
The new ANNs for the off-design efficiency factors for ORC-turbines according to
```bibtex
@Article{ghasemi2013modeling,
  author    = {Ghasemi, Hadi and Paci, Marco and Tizzanini, Alessio and Mitsos, Alexander},
  pages     = {412--428},
  title     = {Modeling and optimization of a binary geothermal power plant},
  volume    = {50},
  doi       = {10.1016/j.energy.2012.10.039},
  journal   = {Energy},
  publisher = {Elsevier},
  year      = {2013},
}
```
and the inverse of the temperature correction factor are implemented in `ORC_components.py`

## Results
The results of the performed calculations from Section 3 of the paper can be found in the following subdirectories:
- 3.1
  - Design optimization for average ambient temperature: `AD`

    `python case_study.py -sup -rec -od -T=289 -subdir=results/AD -time=86400`

  - Operational optimizations with fixed design from AD: `AO`

    `python run_individual.py -subdir=results/AO -start=results/AD/results_289.pickle -fix -n=11`

  - Operational optimizations with fixed design from AD and a **relaxed fluid velocities and variable brine mass flow**: `AO_relaxed`

    `python run_individual.py -subdir=results/AO_relaxed -start=results/AD/results_289.pickle -fix -n=11`

  - Similarly, `results/HD` and `results/HO` contain results for a maximum-temperature design and off-design operation, respectively.

- 3.2
  - Multiple-temperature Design optimization for operational scenarios corresponding to 11 ambient temperatures: `MD`

    `python case_study.py -T 263.15:0.00029 268.15:0.01852 273.15:0.07207 278.15:0.13581 283.15:0.19003 288.15:0.159 293.15:0.13603 298.15:0.12674 303.15:0.09895 308.15:0.05504 313.15:0.00295 -sup -rec -subdir=results/MD -time=86400`

  - Operational optimizations with fixed design from MD: `MO`

    `python run_individual.py -subdir=results/MO -start=results/MD/results_263.15_268.15_273.15_278.15_283.15_288.15_293.15_298.15_303.15_308.15_313.15.pickle -fix -T=-10_41_1 -n=11 -time=86400`

  - Design optimizations for individual temperatures: `SD`

    `python run_individual.py -subdir=results/SD -n=11 -time=86400`

- 3.3
  - Design optimization combining operational scenarios for 11 ambient temperatures and reduced variable ranges: `results/MDRR`
  - Design optimizations for individual temperatures with reduced variable ranges: `results/SDRR`
