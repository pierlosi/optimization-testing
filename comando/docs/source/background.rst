.. This file is part of the COMANDO project which is released under the MIT
.. license. See file LICENSE for full license details.
..
.. AUTHOR: Marco Langiu
.. _background:

##########
Background
##########

This section gives a brief overview of the concept behind COMANDO, the problem
formulations that can be addressed with it, as well as its usage.

*******
Concept
*******

The aim of COMANDO is to enable object-oriented modeling of energy systems and their constituting components for the purpose of creating optimization problems related to the system's design and operation.
Unlike other energy system modeling frameworks, COMANDO does not impose restrictions to maintain the resulting optimization problems in a particular class such as linear programming (LP) or mixed-integer linear programming (MILP).
Instead, component and system models may contain a wide range of nonlinear expressions, including nonconvex and nonsmooth functions.
Furthermore ordinary differential equations describing system dynamics can be included directly, without manual discretization.
In this way COMANDO allows users to create component and system models in a natural fashion, without having to obscure them with optimization-specific implementation details.

By specifying the connections between components, different component models can be combined into a system model.
Based on a finished system model, different optimization problems can be created by specifying objective functions, the operational scenarios and their time-structure to be considered, as well as associated parameter data.

To bring the resulting problem formulations into a form that solvers understand, certain reformulations may be necessary.
COMANDO provides routines for substitution and reformulation of expressions, as well as for automatic approximation via linearization and time-discretization.
Once an acceptable problem formulation is created, it can be translated to different algebraic modeling languages (AMLs) or be passed directly to an appropriate solver for solution.

The details of the different steps of this process are explained in the following.

*****
Usage
*****

The usage of COMANDO can be split into three phases

1. Modeling phase

   - Component model creation
   - System model creation

2. Problem formulation phase

   - problem generation

     - objective selection
     - time-structure selection
     - scenario-structure selection
     - providing data
   - problem reformulation

     - time discretization
     - linearization
     - ...

3. Problem solution phase

   - via a solver interface (e.g., GUROBI, BARON)
   - via an algebraic modeling language (AML) interface (e.g., GAMS, PYOMO)
   - via custom algorithm

     - initialization
     - decomposition
     - ...


Modeling phase
==============

In the modeling phase the system behavior is specified in terms of the
component behavior and system structure.


Components
----------

In COMANDO a component is the basic building block of an energy system.
It represents a model of a generic real-world component, specified via a
collection of mathematical expressions that are given in a symbolic form.

Take for example the following simple description of some generic component `C`:

.. math::
   :nowrap:

   \begin{align}
     C_{\mathrm{out}} &= C_{\mathrm{in}} C_{\eta} \\
     C_{\eta} &= (C_0 + C_1 \, C_{\mathrm{pl}} + C_2 \, C_{\mathrm{pl}}^2) \\
     C_{\mathrm{pl}} &= C_{\mathrm{out}} / C_{\mathrm{out,max}}
   \end{align}

It contains different *symbols* that represent different design and operational quantities of `C` and equations that correlate these symbols with each other.
Here, we may assume that $C_{\mathrm{out,max}}$ is a design quantity, while the symbols $C_{\mathrm{out}}$, $C_{\mathrm{in}}$, $C_{\eta}$, and $C_{\mathrm{pl}}$ are operational quantities and $C_0$, $C_1$, and $C_2$ are some numerical coefficients.

We can create a COMANDO model of this component by deciding which role the different symbols play for some design or operation process.
In COMANDO a distinction is made between quantities that are given as input data and those that are to be determined by the use of the model, i.e., by formulating some optimization problem incorporating the model and solving that problem.

Symbols that act as placeholder for input data can be modeled as a :class:`comando.Parameter`, while those that are to be determined can be modeled as a :class:`comando.Variable` or a  :class:`comando.VariableVector`, depending on whether they represent design- or an operational-related quantity.

.. note::

  Objects of type :class:`.Variable` / :class:`.VariableVector` always represent scalar / vector values while those of type :class:`.Parameter` may represent both, depending on the data that is assigned to them.
  In the context of design and operation, design quantities are always modeled by symbols or composed expressions that represent scalar values, while operational quantities are always modeled by symbolss or composed expressions that represent vector values, with each entry representing the value of the corresponding quantity for a certain operational scenario and timestep.

We may for example model the coefficients $C_0$, $C_1$, and $C_2$ as objects `C_0`, `C_1`, and `C_2` of type :class:`.Parameter`, $C_{\mathrm{out,max}}$ as a :class:`.Variable` `C_out_max` and $C_{\mathrm{out}}$, $C_{\mathrm{in}}$ as objects `C_out`, `C_in` of type :class:`.VariableVector`.

For the definition of component COMANDO provides the :class:`.Component` class.
This class defines methods to create and add the corresponding symbols to :class:`.Component` instances:

.. code:: python

  # Define a new type of component
  class C(comando.Component):
      """A generic component."""

      def __init__(self, label):
          """Initialize the component."""
          super().__init__(label)  # Initialize the parent class (Component)

          # Create parameters
          C_coeffs = [self.make_parameter(i) for i in range(3)]

          # Create design and operational variables
          C_out_max = self.make_design_variable('out_max')
          C_out = self.make_operational_variable('out')
          C_in = self.make_operational_variable('in')

          # Expressions formed with the defined symbols
          C_pl = self.add_expression(C_out / C_out_max, 'part_load')
          C_eta = sum(C_i * C_pl ** i for i, C_i in enumerate(C_coeffs))

          # Add constraints
          self.add_eq_constraint(C_out, C_in * C_eta, 'conversion')
          self.add_le_constraint(C_out, C_out_max, 'output_limit')

          # Add connectors
          self.add_input('INPUT', C_in)
          self.add_output('OUTPUT', C_out)

  # Create an instance of our new class
  c = comando.Component('C')
  c_pl = c.get_expression('part_load')  # Access previously defined expression

Note how instead of introducing additional :class:`.VariableVector` instances for $C_{\eta}$ and $C_{\mathrm{pl}}$, we directly used their definitions and create expressions for these quantities using overloaded Python functions!
The expression for `C_pl` is considered interesting and is therefore stored in `C` under the name 'part_load'.
At a later point, this expression can be accessed via the :func:`.Component.get_expression` method.

We may have similarly defined `C_out` as an expression, but instead we decided to make it an operational variable and introduce the correlation between input and output as an equality constraint.
Similarly, we introduced an inequality constraint specifying that the component's output must be smaller than it's maximum value.

It is also possible to declare operational variables to be 'states', i.e.,
quantities whose time-derivative is given by some algebraic expression.
This allows the consideration of dynamic effects.

.. todo::

  Examples for :func:`.Component.declare_state` and :func:`.Component.make_state`

Another thing to notice are the 'INPUT' and 'OUTPUT' connectors that can be used to connect the component to others within a :class:`comando.System` model.
Components may also assign individual algebraic expressions to connectors, allowing them to be interfaced with other components.
Connectors may be specified as inputs, outputs or bidirectional connectors.
The former two restrict the value of the corresponding expression to be nonnegative and nonpositive, respectively, while the latter does not impose additional restrictions.

Systems
-------

Systems are modeled as collections of interconnected components, i.e., a system
model consists of a set of components and the specification of how their
connectors are connected.
For this purpose COMANDO provides the :class:`comando.System` class.
As the :class:`.System` class inherits from the :class:`.Component` class, systems can define additional expressions, constraints and connectors and be incorporated as subsystems within other systems.

.. todo::

  Example for system creation

Problem formulation phase
=========================

Given a system model, COMANDO can currently be used to create a :class:`comando.Problem`
object, representing a mathematical optimization problem (OP) of the form:

.. \require{boldsymbol}
.. math::
   :nowrap:

   \[
   % DEFINITIONS
   \renewcommand{\bm}[1]{\boldsymbol{\mathbf{#1}}}
   \renewcommand{\min}[1][]{
	   {\underset{#1}{\text{min}\,}}
   }
   \renewcommand{\st}{\mathrm{s.\,t.}}
   \renewcommand{\dv}{\bm{x}}
   \renewcommand{\ov}{\bm{y}_{s}(\cdot)}  % operational variables
   \renewcommand{\ovfunc}{\bm{y}_{s}(t)}  % operational variables
   \renewcommand{\derov}{\dot{\bm{y}}^\mathrm{d}_s(t)}  % derivative of ov
   \renewcommand{\roc}{\bm{f}}                          % specified rate of change for ov
   \renewcommand{\iv}{\bm{y}^\mathrm{d}_{s}(t=0)}       % initial value
   \renewcommand{\design}{I}
   \renewcommand{\operation}{{II}}
   \renewcommand{\of}{F}                       % objective function
   \renewcommand{\dobj}{\of_\design}           % design objective
   \renewcommand{\oobj}{\of_{\operation,s}}    % operational objective
   \renewcommand{\moobj}{\dot{\of}_\operation} % momentary operational objective
   \renewcommand{\myargs}{\big(\dv, \ovfunc, \bm{p}_s(t)\big)}
   %
   \begin{alignedat}{4}
     & \min[\dv] \;\;
       && \dobj(\dv) + \sum_{s \in \mathcal{S}} w_s \, \oobj^*(\dv) \\[-1mm]
     & \underset{\hphantom{\dv, \bm{y}}}{\st}
       && \bm{g}_\design(\dv) \leq \bm{0} \\[-2mm]
     & && \bm{h}_\design(\dv) = \bm{0} \\
     & && \left.\kern-0.75ex
            \begin{array}{lll}
              \oobj^*(\dv) =
              & \min[\bm{y}_s(\cdot)]
                & \oobj(\dv, \ov) = \!\!
                  \int_{\mathcal{T}_s} \!\!
                  \moobj\myargs \, \text{d}t \\
              & \st
                & \bm{y}^\text{d}_{s}(\cdot) \in \bm{y}_{s}(\cdot) \\
              & & \iv = \bm{y}^\text{d}_{s, 0} \\
              & & \left.\kern-0.3ex
                  \begin{array}{ll}
                    \derov = \roc\myargs \!\!\!&\\
                    \bm{g}_\operation\myargs \leq \bm{0} &\\
                    \bm{h}_\operation\myargs = \bm{0} &\\
                    \bm{y}_{s}(t) = [\bm{y}^\text{d}_{s}(t), ...] \\
                    \bm{y}_{s}(t) \in \mathcal{Y}_{s}(t) \subset \mathbb{R}^{n_y} \!\times \mathbb{Z}^{m_y} \!\!\!\!\!\!\!\!\! \\
                  \end{array}
                  \right\} \; \forall t \in \mathcal{T}_s \\
              & & \mathcal{T}_s = \left[0, T_s\right] \\[-2mm]
            \end{array}
          \right\}
            \; \forall s \in \mathcal{S} \\[1mm] %= \{s_1, s_2, \cdots, s_{\card{\mathcal{S}}}\}\\
      & && \dv \in \mathcal{X} \subset \mathbb{R}^{n_x} \!\times \mathbb{Z}^{m_x} \\[-2mm]
      & && \mathcal{S} = \{s_1, s_2, \cdots, s_{|\mathcal{S}|}\}
   \end{alignedat}
   \]

Where $\dv$ and $\ov$ are the vectors of design- and operation-variables, respectively.

Problem generation
------------------

$F_I$ and $F_{II}$ are user-specified scalar and indexed expressions corresponding to
one-time and momentary costs, $\mathcal{T}_s$ are time horizons for scenarios $s$ in from the set $\mathcal{S}$, with corresponding probabilities $w_s$.

The constraints for the OP are automatically generated from the system model,
i.e., all scalar relational expressions are taken as contraints and all indexed
relational expressions are taken as constraints, parametrized by t, and s.

The dependence of the objective and constraint functions on time and scenario can be expressed in terms of the parameter values, which are user input.
This data can currently be given only in discrete form, i.e., for a discrete time and scenario.
The time-steps at which data is available

**If any states were defined in the model the corresponding differential
equations are currently discretized by default.**
An exception is the use of the :mod:`.pyomo_dae` interface; here the time-continuous representation is passed and discretization via collocation can be performed.
Data for the parameter values and initial guesses for the variable values can be provided based on the user-chosen sets T and S.

.. note::

  The time- and scenario structure is only specified during :class:`.Problem` creation and is therefor not available during the *modeling phase*!
  While this may be unintuitive at first, it allows for a clean separation between component and system behavior on the one hand, and problem-specific imlementation details on the other.
  An important consequence of this is, that certain aspects of component behavior such as, e.g., ramping constraints must either be defined in a more general way (e.g., via limiting the derivative of the ramped quantity), or must be deferred and carried out after :class:`.Problem` creation (e.g. via appropriate callbacks).

.. todo::

  Example for :class:`.Problem` creation

.. todo::

  Examples for creating problem-specific constraints such as:

  - Constraints that restrict cumulative quantities (time-integrals of quantities)
  - Constraints that only consider a certain time range
  - Constraints explicitly referencing quantities corresponding to multiple time-points and/or scenarios

Problem reformulation
---------------------

It is possible to use manual or automated reformulations of the original
problem formulation.

.. todo::

  Example reformulations, :mod:`.linearization` module and default discretization scheme in the interfaces.

Problem solution
================

Given a Problem, the user can chose to pass it directly to a solver capable of handling the corresponding problem type, transform the COMANDO Problem to a representation in an AML (currently Pyomo or GAMS), or work directly with the COMANDO representation in a custom algorithm to preprocess or solve the
Problem.
Currently available interfaces as well as their installarion requirements are described in :ref:`interfaces`.

If a solution is obtained, it is loaded back into the COMANDO Problem after the solver terminates.
