# Interfaces to solvers and algebraic modeling languages (AMLs)
<!-- This file is part of the COMANDO project which is released under the MIT
license. See file LICENSE for full license details. -->

COMANDO itself allows to model energy systems and formulate optimization
problems based on the resulting models.
To solve these optimization problems, they need to be passed to a suitable
solver.
This can either be done by creating an appropriate input file and passing it to
a solver or interfacing a solver's application programming interface (API) if
present.
Alternatively we can communicate the problem formulation to an algebraic
modeling language (AML) - again via some input file or via API - which takes
care of the interaction with the solver.
The benefit of using an AML is that multiple solvers can be tested easily, a
minor downside is that an additional representation (the one for the AML) needs
to be created.

Both solver and AML interfaces are either text-based (creating input files) or
use an API.
Using any interface generally results in the following steps to be performed:

- translate a given problem to an appropriate representation
- instructs the solver to solve the problem with the given options
- reads back the solution from result files or API-objects
- updates the COMANDO variables with the optimal values

## Available Solver interfaces

We provide the following interfaces:
<!--
### CPLEX (TODO)

Register and download CPLEX solver from [here](https://www-03.ibm.com/isc/esd/dswdown/searchPartNumber.wss?partNumber=CJ6BPML)

From a commandline navigate to
- Linux: `/opt/ibm/ILOG/CPLEX?Studio1210/python`
- Mac: `/Applications/CPLEX_Studio1210/python`
- Windows `...`
and run

```shell
python setup install
```

To check if everything worked properly, run python and do:
```python
import docplex
``` -->

### GUROBI (API)

Problem classes:

- (MI)LP
- (MI)QP
- (MI)QCQP

#### Installation

Download and install the GUROBI solver (for academics you can create an account and obtain a free license for noncommercial use [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/))

From a commandline navigate to

- Linux: `/opt/gurobi902/linux64`
- Mac: `/Library/gurobi902/mac64`
- Windows `C:\gurobi\win64`

and run

```shell
python setup install
```

On Linux you may have to manually set additional paths (e.g. in
`~/.bash_profile`)

```shell
export GUROBI_HOME="/opt/gurobi902/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
```

To check if everything worked properly, run python and do:

```python
import gurobipy
```

### BARON (Text-based)

Problem classes:

- (MI)LP
- (MI)QP
- (MI)QCQP
- (MI)NLP

### SCIP (API)

Problem classes:

- (MI)LP
- (MI)QP
- (MI)QCQP
- (MI)NLP

#### Installation
Either download and install the prebuild binaries or download and compile the source from [the SCIPOPT website](https://www.scipopt.org/index.php#download).

Then install the pyscipopt python interface following [these](https://github.com/scipopt/PySCIPOpt/blob/master/INSTALL.md) instructions.

### MAiNGO (Text-based & API)

Problem classes:

- (MI)LP
- (MI)QP
- (MI)QCQP
- (MI)NLP


#### Installation

Using the MAiNGO interface requires the MAiNGO solver, which must currently
build from source code.
For this you will need CMake, compilers for C++11 and Fortran77.

<!-- We tested: -->
<!-- - Linux (Ubuntu): gcc / gfortran 9.3.0 (with 8 even specifying -lstdc++fs doesn't work!) -->
<!-- - Mac: clang, gfortran -->
<!-- - Windows: VisualStudio, intel fortran compiler -->

An optional dependency is CPLEX as a MILP subsolver, alternatively the open
source solver CBC is used.

<!--
CPLEX_BIN_DIR                    /opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux                                                                        
CPLEX_CONCERT_INCLUDE_DIR        /opt/ibm/ILOG/CPLEX_Studio1210/concert/include                                                                               
CPLEX_CONCERT_LIBRARY            /opt/ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/libconcert.a                                              
CPLEX_ILOCPLEX_LIBRARY           /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libilocplex.a                                               
CPLEX_INCLUDE_DIR                /opt/ibm/ILOG/CPLEX_Studio1210/cplex/include                                                                                 
CPLEX_LIBRARY                    /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libcplex.a                                                  
CPLEX_ROOT_DIR                   /opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux
-->

To obtain the source code, we recommend using git:

```shell
git clone https://git.rwth-aachen.de/avt.svt/public/maingo.git
cd maingo
git submodule init
git submodule update -j 1
git submodule foreach git checkout master
git submodule foreach git pull
```

When using the commandline you can configure the cmake run via

```shell
mkdir build && cd build
ccmake ..
```

Now you may need to point CMake to the location of CPLEX if available.
At this step you must activate the `MAiNGO_build_python_interface` switch if
you want to use the API interface.
Hit configure and generate to create the required files for building MAiNGO.

On Linux/Mac you can build the software with

```shell
make
```

On Windows you may use Visual Studio

The build process should produce the MAiNGO executable and pymaingo shared library
`pymaingo.<your_python_version_number>-<your_platform>.<so OR pyd>`
which needs to be added to the `PYTHONPATH` environment variable to be
discoverable by python.
For this add the following line to your `~/.bash_profile` or similar:

```shell
export PYTHONPATH="<PATH_TO...>/maingo/build:$PYTHONPATH"
```

To check if everything worked properly, run python and do:

```python
import pymaingo
```

#### Usage

```python
from comando.interfaces.maingo_api import MaingoProblem

mp = MaingoProblem(P)
solver, ret = mp.solve(epsilonA=1e-12, outstreamVerbosity=1)
# `solver` is a maingo::MAiNGO solver object that can be queried for
# solve-related information, also see `help(pymaingo.MAiNGO)`.
# ret is a maingo::RETCODE, also see `help(pymaingo.RETCODE)`.
```

## Available AML interfaces

### GAMS (Text-based)

The GAMS interface relies on the GAMS AML, which you can get (along with a trial license) [here](https://www.gams.com/download/).

To use the COMANDO interface to GAMS, simply insyall the software and ensure it
is found when running a commandline with

```shell
gams
```

#### Usage

```python
# With a comando.Problem P...

from comando.interfaces.pyomo import to_pyomo

# NOTE: In contrast to us, pyomo refers to objects representing optimization
#       problems as 'models', hence 'm' is typically used!
m = to_pyomo(P)
# Now different locally installed solvers can be used, e.g., BARON (Requires
# standalone solver and license, not via GAMS!)
# The tee=True option results in the solver output to be displayed.
# Other, solver-specific options can be added, such as BARON's `MaxTime` here.
res = m.solve('baron', tee=True, MaxTime=300)
# res is a Pyomo results object, see the Pyomo documentation for details
```

### Pyomo, Pyomo.DAE (API)

The Pyomo and Pyomo.DAE interfaces rely on the python package `pyomo`, which
can be installed by running

```shell
pip install pyomo
```

#### Usage

```python
# With a comando.Problem P...

from comando.interfaces.gams import solve

# Create a .gms file and pass it to GAMS.
# For a problem class other than MINLP, the `model_type` must be specified
# explicitly.
# Any additional GAMS options can be specified.
ret = solve(P, model_type='NLP', NLP='baron', optCR=1e-3)
# ret is the GAMS return code, also see:
# https://www.gams.com/latest/docs/UG_GAMSReturnCodes.html
if ret == 0:
    print('Εὕρηκα!')
elif ret == 7:
    raise RuntimeError("GAMS couldn't find a valid license!")
```
