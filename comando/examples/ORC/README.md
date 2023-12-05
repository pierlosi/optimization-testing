<!-- This file is part of the COMANDO project which is released under the MIT
license. See file LICENSE for full license details. -->
# ORC case Study

This case study considers the maximization of net power production for an organic Rankine cycle (ORC) fueled by geothermal brine and cooled by a cooling water cycle whose heat is dissipated through fans.

This case study is based on:

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

The directory `Network Databank` contains parameters for Artificial Neural Networks (ANNs) describing thermodynamic properties of isobutane and was kindly provided by Wolfgang R. Huster and Artur M. Schweidtmann.
The ANNs were trained with data from [CoolProp](http://www.coolprop.org/) (last accessed February 1st 2021) using the reference state `NBP`.

## Requirements

This case study requires the following additional features to be installed:
- pyomo
- cpp-backend

This can be done from the base directory of this repository using
```shell
pip install .[pyomo,cpp-backend]
```

Further requirements are the following solvers
- BARON 20.10.16
- CPLEX 12.10
- MAiNGO 0.3