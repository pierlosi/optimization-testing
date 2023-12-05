.. This file is part of the COMANDO project which is released under the MIT
.. license. See file LICENSE for full license details.
..
.. AUTHOR: Marco Langiu
.. _installation:

############
Installation
############

Currently COMANDO is tested and should thus work with python 3.7 and 3.8.
COMANDO itself only has :mod:`sympy` and :mod:`pandas` as its dependencies,
however, while it can be used on its own for structural analysis,
reformulation, or evauation of models and problems, actually solving
optimization problems requires the use of an interface.

The available interfaces and their installation instructions can be found in
:ref:`interfaces`.


For the solution of optimization problems COMANDO provides interfaces to
different solvers or relies on the interfaces provided by algebraic modeling
languages (AMLs).
Interfaces are either **text-based** (i.e., they allow for the generation of an
input file for a solver or an AML) or **api-based** (i.e., they make use of a
Python interface provided by the solver or AML).

Currently we provide the following interfaces.

- text-based:

  - BARON (solver)
  - GAMS (AML)
  - MAiNGO (solver)

- API-based:

  - Pyomo / Pyomo.DAE (AML)
  - Gurobi (solver)
  - MAiNGO (solver)

Until COMANDO is listed on pypi you have the following options to install
COMANDO:

- Download the GitLab repository as a zip file

- clone the GitLab repository

Once you did that, you can install via

  .. code-block:: shell

    # from the parent directory of this repository...
    pip install .

Uninstallation
--------------

.. code-block:: shell

  # from anywhere
  pip uninstall comando