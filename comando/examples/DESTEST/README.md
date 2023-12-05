<!-- This file is part of the COMANDO project which is released under the MIT
license. See file LICENSE for full license details. -->
# DESTEST case study

This case study is based on

    @InProceedings{saelens2020towards,
        author    = {Saelens, Dirk and de Jaeger, Ina and B{\"u}nning, Felix and Mans, Michael and Maccarini, Alessandro and Garreau, Enora and R{\o}nneseth, {\O}ystein and Sartori, Igor and Vandermeulen, Annelies and {van der Heijde}, Bram and Helsen, Lieve},
        booktitle = {{Proc. 16th Int. Conf. Build. Simul.}},
        date      = {2019-09-02},
        title     = {{Towards a DESTEST: a District Energy Simulation Test Developed in IBPSA Project 1}},
        doi       = {10.26868/25222708.2019.210806},
        editor    = {Corrado, Vincenzo and Fabrizio, Enrico and Gasparella, Andrea and Patuzzi, Francesco},
        pages     = {3569--3577},
        publisher = {IBPSA},
        series    = {Building Simulation Conference proceedings},
        address   = {Rome, Italy},
        year      = {2020},
    }

The file `data/pipe_data.csv` contains network data, computed based on the file `Pipe_data.csv`, available at [https://github.com/ibpsa/project1-destest/tree/master/Networks/CE_1/input_data](https://github.com/ibpsa/project1-destest/tree/master/Networks/CE_1/input_data) (last accessed February 1st 2021) and the file `data/data.csv` contains demand data, computed based of the file `simple_district.csv`, available at [https://github.com/ibpsa/project1/tree/master/wp_3_1_destest/Buildings/SimpleDistrict/Results/SimpleDistrict_IDEAS](https://github.com/ibpsa/project1/tree/master/wp_3_1_destest/Buildings/SimpleDistrict/Results/SimpleDistrict_IDEAS) (last accessed February 1st 2021).
Both `Pipe_data.csv` and `simple_district.csv` are licensed under an [adapted BSD 3-clause license](https://github.com/ibpsa/project1-destest/blob/master/license.md).
