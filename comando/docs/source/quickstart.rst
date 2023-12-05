.. This file is part of the COMANDO project which is released under the MIT
.. license. See file LICENSE for full license details.
..
.. AUTHOR: Marco Langiu
.. _quickstart:

################
Quickstart Guide
################

This section contains code examples for the use of COMANDO

.. code-block:: python

  # With some energy system model ES = ...

  design_objective = ...  # $F_1$ expression
  operational_objective = ...  # $\dot{F}_2$ expression

  P = ES.create_problem(design_objective,
                        operational_objective,
                        timesteps=timesteps,
                        name=f"min {'|'.join(objectives)}")

  from comando.linearization import linearize

  # Automatically linearize:
  # - find all nonlinear expressions
  # - evaluate them on an equidistant grid with 3 breakpoints per variable
  # - create a triangulation of the expression values
  # - encode the triangulation as linear constraints using the
  #   convex-combination method
  P_lin = linearize(P, n_bp=3, method='convex_combination')

  # Use the Pyomo interface to interact with solvers
  from comando.interfaces.pyomo import to_pyomo

  # convert P to Pyomo 'model' objects
  m = to_pyomo(P)
  m_lin = to_pyomo(P_lin)

  # Solve with open source solvers (need to be installed beforehand)...
  res_lin = m_lin.solve('cbc')  # ...the MILP approximation globally
  res = m.solve('ipopt')  # ...locally using a nonlinear solver
  res_glob = m.solve('couenne')  # ...globally
