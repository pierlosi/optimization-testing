"""Methods to generate (possibly mixed-integer) linear problem formulations."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: David Shu, Marco Langiu
import itertools

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, Delaunay
from scipy.spatial.qhull import QhullError
from sympy import Rel, Le, Eq

import comando


class Namer:
    """Generator for names of the form {prefix}{i}."""

    def __init__(self, prefix, i_0=0):
        """Create a namer with a given prefix."""
        self.prefix = prefix
        self.i = i_0

    def get(self):
        """Get the next generated name."""
        res = f'{self.prefix}{self.i}'
        self.i += 1
        return res


def piecewise_linear(points):
    """Approximate a univariate function via a discrete mapping.

    compute the axis intersects y0_i and slopes dydx_i for the linear segments
    connecting the given points in sorted order.

    Arguments
    ---------
    points : Mapping
        mapping from input to output values

    Returns
    -------
    x : list
        sorted list of x values
    y0 : list
        len(points) - 1 axis intersections for line segments
    dydx : list
        len(points) - 1 slopes for line segments
    """
    x = sorted(points.keys())
    y = [points[x_i] for x_i in x]
    dydx = [(y_ub - y_lb) / (x_ub - x_lb)
            for x_lb, x_ub, y_lb, y_ub
            in zip(x, x[1:], y, y[1:])]
    y0 = [y_i - dydx_i * x_i for y_i, dydx_i, x_i in zip(x, dydx, x)]
    return x, y0, dydx


def _is_linear(expr, vars=None):
    """Test whether the expression 'expr' is linear in variables 'vars'.

    Approach: The derivative of the expression with respect to all variables in
    'vars' may no longer contain any element of 'vars'.

    Arguments
    ---------
    expr: comando.Expression
        expression to be tested
    vars: set
        variables with respect to which to test for linearity

    Returns
    -------
    is_linear : boolean
        Boolean reflecting whether the expression is linear in all variables
    """
    if vars is None:
        vars = comando.utility.get_vars(expr)
    for var in vars:
        if expr.diff(var).free_symbols.isdisjoint(vars):
            continue
        else:
            return False
    return True


def cached(reformulator):
    """Equip the `reformulator` function with a cache."""
    def wrapped(*args):  # args are `expr, n_bp, method`
        if args in wrapped.cache:
            # NOTE: We only want to return the modifications the first time we
            #       encounter `expr`, on subsequent passes we only need its
            #       reformulation as the modifications were already registered!
            return wrapped.cache[args][0], [], [], []
        reformulation, constraints, dv, ov = reformulator(*args)
        if reformulation is not args[0]:  # reformulation differs, cache it
            wrapped.cache[args] = reformulation, constraints, dv, ov
        return reformulation, constraints, dv, ov

    wrapped.cache = {}
    return wrapped


def linearize(P, n_bp, method):
    """Return a (possibly mixed-integer) linear approximation of problem P.

    In the linearization process, all expressions that will form part of an
    optimization problem will be analyzed and reformulated into linear forms if
    necessary.
    In this process of reformulation new variables and constraints may need to
    be introduced.
    To refrain from modifying the passed energy system instance the new set of
    variables, a dictionaty of constraints and an objective expression are
    generated and returned.

    Arguments
    ---------
    P : comando.Problem
        The problem to be linearized
    n_bp : number
        the number of breakpoints used per variable
    method : str
        method used to encode the triangulation over the generated simplices
        currently either 'convex_combination' or 'multiple_choice'.

    Returns
    -------
    P_lin: comando.Problem
        A linear reformulation of problem P possibly with additional variables
        and constraints.
    """
    # Prepare a new cache for nonlinear expressions the caller can access later
    linearize_expr.cache = linearize.cache = {}

    # linearize constraints
    # UNUSED: variables = P.design_variables.union(P.operational_variables)
    constraints = {}
    new_constraints = set()

    def add_constraints(lin_cons):
        """Add the linear constraints if they are new."""
        for lin_con in lin_cons:
            # TODO: We might be able to drop this test due to caching!
            if lin_con in new_constraints:
                continue
            if not comando.is_trivial(lin_con, "Linearization constraint"):
                new_constraints.add(lin_con)
                constraints[f'lin_constr_{add_constraints.n}'] = lin_con
                add_constraints.n += 1

    add_constraints.n = 0  # Counter for auxiliary linear constraints

    # Note: MILP solvers can be affected by the order of constraints, even when
    #       their random seed is fixed! To allow for consistent results we
    #       process constraints in lexicographical order.
    for con_id, con in sorted(P.constraints.items(), key=lambda item: item[0]):
        con_sides = [con.lhs, con.rhs]
        for i, xhs in enumerate(con_sides):
            if not _is_linear(xhs):
                # TODO: We do not need to return the new variables explicitly
                #       as any relevant variable will be part of an expression
                #       and thus can be determined at the end of this method!
                #       The only problem with this is the objective, which may
                #       contain individual elements of operational variables,
                #       that the suggested approach would wrongly count as
                #       design variables.
                #       A solution for this is to only check for operational
                #       variables and instantiate them!
                coeff, expr = xhs.as_coeff_Mul()
                expr, lin_cons, lin_dvars, lin_ovars = linearize_expr(
                    expr, n_bp, method)
                con_sides[i] = coeff * expr
                add_constraints(lin_cons)
        if con_sides is [con.lhs, con.rhs]:
            constraints[con_id] = con
        else:
            lin_con = Rel(*con_sides, con.rel_op)
            con_repr = f"Linearization of constraint {con_id}"
            if not comando.is_trivial(lin_con, con_repr):
                constraints[con_id] = lin_con

    # linearize objective
    objective = []
    for expr in P.design_objective, P.operational_objective:
        if _is_linear(expr):
            objective.append(expr)
        else:
            lin_expr, lin_cons, lin_dvars, lin_ovars = \
                linearize_expr(expr, n_bp, method)
            objective.append(lin_expr)

            # register auxiliary variables
            # UNUSED: variables.update(lin_dvars, lin_ovars)

            # register auxiliary linearization constraints
            add_constraints(lin_cons)

    # TODO: STATES
    if P.states:
        from warnings import warn
        warn("Linearization has not yet been implemented for states!")

    # NOTE: Already happens in P_lin!
    # for lin_ovar in lin_ovars:
    #     lin_ovar.instantiate(P.index)

    return comando.Problem(*objective, constraints, P.states, P.timesteps,
                           P.scenario_weights, name=f'{P.name}_linearized')


@cached
def linearize_expr(expr, n_bp, method):
    """Linearize the expression 'expr' with `n_bp` breakpoints using `method`.

    Arguments
    ---------
    expr: comando.Expression
        expression to be linearized
    n_bp : number
        the number of breakpoints used per variable
    method : str
        method used to encode the triangulation over the generated simplices
        currently either 'convex_combination' or 'multiple_choice'.

    Returns
    -------
    lin_expr : comando.Expression
        linearized expression to replace `expr`
    lin_constr : list
        newly introduced constraints
    lin_dvars : list
        newly introduced design variables
    lin_ovars : list
        newly introduced operational variables
    """
    # Share cache with callees
    glover_reformulation.cache = linearize_expr.cache
    piecewise_linear_reformulation.cache = linearize_expr.cache

    # determine all variables in expr
    if expr.is_Add:  # linearize every summand separately
        res = [0, [], [], []]  # [lin_expr, lin_cons, lin_dvars, lin_ovars]
        for arg in expr.args:
            ret = linearize_expr(arg, n_bp, method)
            for i, element in enumerate(ret):
                res[i] += element
        return res
    else:
        if _is_linear(expr):  # expr is constant or linear
            return expr, [], [], []
        else:  # expr is nonlinear
            # `expr` may be a product coeff * _expr, where coeff is a constant
            # and _expr is a nonlinear function.
            # In the linearization functions `create_piecewise_linear_f` and
            # `glover_reformulation`, we linearize _expr only (without the
            # coefficient coeff), to improve the use of caching. The
            # coefficient coeff is multiplied to the result after linearization
            coeff, _expr = expr.as_coeff_Mul()

            if _expr.is_Mul:
                try:
                    lin_expr, lin_cons, lin_dvars, lin_ovars = \
                        glover_reformulation(_expr, n_bp, method)
                    return lin_expr * coeff, lin_cons, lin_dvars, lin_ovars
                except NoGloverReformulationPossible:
                    pass
            lin_expr, lin_cons, lin_dvars, lin_ovars = \
                piecewise_linear_reformulation(_expr, n_bp, method)
            lin_expr = lin_expr * coeff
    return lin_expr, lin_cons, lin_dvars, lin_ovars


class NoGloverReformulationPossible(Exception): pass


glover_namer = Namer('glover_var_')

# @lru_cache(maxsize=cache_size)
@cached
def glover_reformulation(expr, n_bp, method):
    """Create a glover reformulation of the expression `b_var * g`.

    The Glover reformulation is a linear reformulation for the product of a
    binary variable `b_var` and an algebraic expression `g`.
    The reformulation substitutes the expression `b_var * g` with `lin_expr`.
    If `b_var` is 0, additional constraints ensure that `lin_expr` is bounded
    above and below by 0. If `b_var` is 1, `lin_expr` is either equivalent to
    `g` or to its linearization if `g` is nonlinear.

    **NOTE:** The expression `g` must be independent of `b_var`!

    Arguments
    ---------
    expr: comando.Expression
        expression to be linearized
    n_bp : number
        the number of breakpoints used per variable
    method : str
        method used to encode the triangulation over the generated simplices
        currently either 'convex_combination' or 'multiple_choice'.

    Returns
    -------
    lin_expr : comando.Expression
        linearized expression to replace `expr`
    lin_constr : list
        newly introduced constraints
    lin_dvars : list
        newly introduced design variables
    lin_ovars : list
        newly introduced operational variables

    Glover, Fred (1975): Improved Linear Integer Programming Formulations of
    Nonlinear Integer Problems. In: Management Science 22 (4), S. 455–460.
    DOI: 10.1287/mnsc.22.4.455.
    """
    # Share cache with callee
    linearize_expr.cache = glover_reformulation.cache

    if not expr.is_Mul:  # ensuring that we are dealing with a multiplication
        raise NoGloverReformulationPossible

    # testing whether glover reformulation is actually possible
    def selector(arg):
        """Decide if an expression is a binary variable."""
        try:
            return bool(arg.is_binary)
        except AttributeError:
            return False

    others, b_vars = comando.utility.split(expr.args, selector)
    for b_var in b_vars:
        # NOTE: We rely on the fact that sympy aggregates expressions, i.e., we
        #       assume here that `b_vars` never contains two instanced of the
        #       same `b_var`. Only if this holds we can compose `g` as below!
        if all(b_var not in arg.free_symbols for arg in others):
            other_b_vars = {*b_vars} - {b_var}
            g = expr.func(*other_b_vars, *others)
            break  # exiting here makes `b_var` and `g` the glover form
    else:  # We did not find `g`
        raise NoGloverReformulationPossible

    if _is_linear(g):  # if g is linear we don't need to reformulate
        g_lin, lin_constr, lin_dvars, lin_ovars = g, [], [], []
    else:
        g_lin, lin_constr, lin_dvars, lin_ovars = \
            linearize_expr(g, n_bp, method)

    # linearization constraints and variables
    g_min, g_max = comando.utility.bounds(g_lin)
    # var_name = f'xi_{{{id(b_var)}}}_{{{id(g)}}}'  # TODO: Naming!
    var_name = glover_namer.get()

    # Ensure that the lin_expr may take on a value of 0 if b_var = 0
    var_lb = min(0, g_min)
    var_ub = max(0, g_max)
    if comando.utility.indexed(expr):
        lin_expr = comando.VariableVector(var_name, bounds=(var_lb, var_ub))
        lin_ovars.append(lin_expr)
    else:
        lin_expr = comando.Variable(var_name, bounds=(var_lb, var_ub))
        lin_dvars.append(lin_expr)
    # (1 - b_var) * g_min <= g_lin - lin_expr <= (1 - b_var) * g_max)
    lin_constr.append(Le((1 - b_var) * g_min, g_lin - lin_expr))
    lin_constr.append(Le(g_lin - lin_expr, (1 - b_var) * g_max))
    # b_var * g_min <= lin_expr <= b_var * g_max
    if g_min != 0:
        lin_constr.append(Le(b_var * g_min, lin_expr))
    if g_max != 0:
        lin_constr.append(Le(lin_expr, b_var * g_max))
    # NOTE: Since b_var is binary we only need to consider two cases:
    #       1) b_var=0: From the second set of constraints follows lin_expr = 0
    #                   This and the first set of constraints imply that
    #                   g_min <= g_lin <= g_max which is True since g_min and
    #                   g_max are bounds computed by interval arithmetic.
    #       2) b_var=1: From the first set of constraints follows that
    #                   g_lin = lin_expr, from the second set of constraints
    #                   follows that g_min <= lin_expr <= g_max, which is True
    #                   by the definition of lin_expr's bounds.
    return lin_expr, lin_constr, lin_dvars, lin_ovars


simplex_namer = Namer('Simplex_var_')
lambda_namer = Namer('Lambda_var_')
aux_namer = Namer('Aux_var_')

# @lru_cache(maxsize=cache_size)
@cached
def piecewise_linear_reformulation(expr, n_bp, method):
    """Create a piecewise linear approximation of the expression `expr`.

    Arguments
    ---------
    expr: comando.Expression
        expression to be linearized
    n_bp : number
        the number of breakpoints used per variable
    method : str
        method used to encode the triangulation over the generated simplices
        currently either 'convex_combination' or 'multiple_choice'.

    Returns
    -------
    lin_expr : comando.Expression
        linearized expression to replace `expr`
    lin_constr : list
        newly introduced constraints
    lin_dvars : list
        newly introduced design variables
    lin_ovars : list
        newly introduced operational variables
    """
    assert method in ['convex_combination', 'multiple_choice'], \
        'Linearization method, {}, not implemented.'.format(method)

    if _is_linear(expr):  # quick exit
        return expr, [], [], []
    else:
        #######################################################################
        # Setup
        #######################################################################
        vars = comando.utility.get_vars(expr)
        d = len(vars)
        # unique identifier for f, used in names of auxiliary variables
        _id = id(expr)
        # Determine which type of variable is to be used: design variables are
        # currently implemented as vectors.
        # NOTE: The assumption is that all auxiliary variables will be of the
        #       same type, based on whether f is indexed or not!
        indexed = comando.utility.indexed(expr)
        # Var = comando.Variable  # ov is Variable
        if indexed:  # ov is VariableVector
            Var = comando.VariableVector  # ov is VariableVector
        else:  # ov is VariableVector
            Var = comando.Variable  # ov is VariableVector

        # linearization constraints
        lin_constr = []

        # New variables (lin_vars[0]: design, lin_vars[1]: operational)
        lin_vars = ([], [])

        #######################################################################
        # Vertices
        #######################################################################
        # create unique identifiers for the vertices:
        # vertex_id = list(range(n_bp ** d))
        # calculate the breakpoints of each variable
        #   e.g. with n_bp = 3 if `f.free_symbols == {x, y}`, x is in [0, 10]
        #   and y is in [5, 10], then bp = [[0.0, 5.0, 10.0], [5.0, 7.5, 10.0]]
        bp = list()
        for v in vars:
            # FIXME: We probably want a more rigorous way to handle bounds:
            #        The big constants used as fallbacks for unbounded
            #        variables may be worse than throwing an error and asking
            #        the user to provide meaningful bounds.
            #        Also operational variables may have different bounds for
            #        each index for which we currently have no implementation!
            bounds = v._bounds if comando.utility.indexed(v) else v.bounds
            lb = -1e10 if bounds[0] in (None, -comando.INF) else bounds[0]
            ub = 1e10 if bounds[1] in (None, comando.INF) else bounds[1]
            bp.append([lb + i * (ub - lb) / (n_bp - 1) for i in range(n_bp)])

        # calculate the x values of the vertices for the linearization
        vertex_x = np.array([*itertools.product(*bp)], 'd')
        # calculate the function values of the vertices
        f = comando.utility.lambdify(expr, eval_params=False)

        # vertex_f = np.fromiter((f(*vert) for vert in vertex_x), 'd')
        vertex_f = [f(*vert) for vert in vertex_x]

        # create dataframe for vertices
        data = {var: [v[i] for v in vertex_x] for i, var in enumerate(vars)}
        # data.update({'vertex_id': vertex_id, 'vertex_f': vertex_f})
        data['vertex_f'] = vertex_f
        vertices = pd.DataFrame(data)

        #######################################################################
        # Simplices
        #######################################################################
        # get vertices for each simplex
        if d == 1:
            vertices_per_simplex = [[i, i + 1] for i in range(n_bp - 1)]
        else:
            try:
                tri = Delaunay(vertex_x, qhull_options='QJ QbB')
            except QhullError as e:
                if "Use option 'Qz'" in str(e):
                    tri = Delaunay(vertex_x, qhull_options='QJ QbB Qz')
                else:
                    raise
            simps = tri.simplices
            keep = np.empty(len(simps), dtype=bool)
            for i, t in enumerate(simps):
                # A = np.hstack((vertex_x[t], np.ones([1, d + 1]).T))
                # keep[i] = abs(np.linalg.det(A)) >= 1E-10
                # Note: Provideds the tolerance is chosen adequately, the above
                #       test is more general than the one we use below!
                #       However in high dimension the above test can be very
                #       expensive and in the uniform axis-aligned grids that we
                #       currently consider, we can test for coplanarity simply
                #       by checking for vertices that share the same coordinate
                #       value in any dimension!
                keep[i] = not any((row == row[0]).all()
                                  for row in vertex_x[t].T)
            vertices_per_simplex = simps[keep]

            # NOTE: In 2-D the triangulation can be visualized as follows:
            #       import matplotlib.pyplot as plt
            #       plt.triplot(vertex_x[:, 0], vertex_x[:,1], tri.simplices)
            #       plt.plot(vertex_x[:, 0], vertex_x[:,1], 'o')
            #       plt.show()

        # TODO: Should be implemented more like this...
        # simplices = pd.DataFrame(vertices_per_simplex)
        # create dataframe for simplices
        data = {'simplex_id': range(len(vertices_per_simplex)),
                'vertices': list(vertices_per_simplex)}
        simplices = pd.DataFrame(data)

        # TODO: The 'vertices' column in `simplices` is a Series of (list)
        #       **objects**, which defies much of the Series's purpose. A better
        #       alternative would be to stores the numbers in vertices:
        # simplices = pd.DataFrame(vertices_per_simplex)

        # create a binary `b_simplex` for every simplex (independent of method)

        # ov is VariableVector
        b_vars = [Var(simplex_namer.get(), domain=comando.BINARY)
                  for i in simplices.index]

        lin_vars[indexed].extend(b_vars)
        simplices['b_simplex_vars'] = b_vars
        # only one of the simplices can be active, i.e., 1
        lin_constr.append(Le(simplices.b_simplex_vars.sum(), 1))

        #######################################################################
        # Linearization methods
        #######################################################################
        if method == 'convex_combination':
            # create a convex combination variable 'lambda' for each vertex

            # ov is Variable
            # c_vars = [Var(f'lambda_{_id}_{i}', bounds=(0, 1), indexed=indexed)
            #           for i in vertices.index]

            # ov is VariableVector
            c_vars = [Var(lambda_namer.get(), bounds=(0, 1))
                      for i in vertices.index]

            vertices['lambda_vars'] = c_vars

            # convex combination of all variables, x = sum_v(lambda_v * x_v)
            lin_constr.extend(Eq(
                    vertices.lambda_vars.mul(vertices[v]).sum(), v)
                for v in vars)

            # convex combination of the function values
            # f(x) = sum_v(lambda_v*f(x_v))
            lin_expr = vertices.lambda_vars.mul(vertices.vertex_f).sum()

            # coupling of linearization variables:
            # lambdas can only be greater than one if one of adjacent simplices
            # is active (b_simplex = 1)
            for n_v in vertices.index:
                # sum over all contributions from simplices at vertex n_v
                s = simplices.loc[[n_v in vtcs for vtcs in simplices.vertices],
                                  'b_simplex_vars'].sum()
                # HOTFIX for older numpy versions where [].sum() = False
                if s is False:
                    s = 0
                lin_constr.append(Le(vertices.at[n_v, 'lambda_vars'], s))

            # the sum of all lambda has to be equal to the sum of b_simplex
            lin_constr.append(Eq(vertices.lambda_vars.sum(),
                                 simplices.b_simplex_vars.sum()))

            # classify linearization variables as design or operational type
            lin_vars[indexed].extend(c_vars)
        elif method == 'multiple_choice':
            # TODO: c_vars for convex_combination above is added to vertices
            #       the c_vars here are added to simplices. Also the created
            #       variables are named and created differently. This probably
            #       means we should use different names for the two!

            # create copies of each variable for each simplex.
            # NOTE: The bounds of the auxiliary variables are not important as
            # their values are always limited by the new constraints.
            # We can provide meaningful bounds, by using the bounds of the
            # original variable, however, we always have to include 0, as all
            # auxiliary variables can be 0, if the simplex is inactive.
            c_vars = {}
            for v in vars:

                bounds = min(v._bounds[0], 0), max(v._bounds[1], 0)
                # ov is Variable
                # simplices[v] = c_vars[v] = [Var(f'aux_{v}_{_id}_{i}',
                #                                 bounds=bounds, indexed=indexed)
                #                             for i in simplices.index]

                # ov is VariableVector
                simplices[v] = c_vars[v] = [Var(aux_namer.get(),
                                                bounds=bounds)
                                            for i in simplices.index]

                # the variable value is the sum of all the auxiliary variables
                lin_constr.append(Eq(v, simplices[v].sum()))

            # use linear equation to approximate the function value f within a
            # simplex s:
            # f_s(x_s) = m_s * x_s + c_s (m_s and x_s are constants, x_s is the
            # vector of all variables)
            # For multiple simplices, we want f_s(x_s) = 0 and x_s = 0, if x is
            # not in the simplex s. Hence, change f_s(x) to  m_s * x_s +
            # b_s * c_s, where b_s is a shorthand for b_simplex_vars.
            # The parameters in m_s and c_s are determined by solving a linear
            # system of equations using the values for x and f(x) at the
            # vertices. f(x) = m_s * x + c_s is brought in the form
            # A * [m_s, c_s] = b <=> A * c = b, which is solved for c.
            # Finally sum over all simplices: f(x) = sum_s(f_s(x))
            lin_expr = 0  # initialize

            for n_s in simplices.index:
                vtcs = simplices.at[n_s, 'vertices']  # vertices in a simplex
                vertex_vals = vertices.loc[vtcs, vars]  # values for all vars
                b_simplex_var = simplices.at[n_s, 'b_simplex_vars']
                s_vars = simplices.loc[n_s, vars]  # auxiliary simplex vars
                A = np.hstack((vertex_x[vtcs], np.ones([1, d + 1]).T))
                b = vertices.loc[vtcs, 'vertex_f'].values
                c = np.linalg.solve(A, b)
                lin_expr += np.dot(c, [*s_vars, b_simplex_var])

                # describe the simplices using sets of linear equations in the
                # form: A * x_s <= b. This can be reformulated as:
                # lin_eq * [x_s, 1] <= 0
                # The equation parameters 'lin_eq' are calculated using
                # scipy.spatial.ConvexHull if x_s contains more than 1 variable
                if d == 1:
                    b_min = float(vertex_vals.min())
                    b_max = float(vertex_vals.max())
                    con1 = Le(s_vars[0], b_max * b_simplex_var)
                    con2 = Le(b_min * b_simplex_var, s_vars[0])
                    for con in con1, con2:
                        if con.is_Boolean:
                            if not con:
                                raise comando.ImpossibleConstraintException(
                                    'An impossible constraint ocurred during '
                                    f'linearization of expression\n {expr}!')
                            continue  # This constraint is always satisfied!
                        lin_constr.append(con)
                else:
                    convex_h = ConvexHull(vertex_vals)
                    for n_eq, eq in enumerate(convex_h.equations):
                        con = Le(np.dot(eq[:-1], s_vars), -eq[-1] * b_simplex_var)
                        if con.is_Boolean:
                            if not con:
                                raise comando.ImpossibleConstraintException(
                                    'An impossible constraint ocurred during '
                                    f'linearization of expression\n {expr}!')
                            continue  # This constraint is always satisfied!
                        lin_constr.append(con)
            lin_vars[indexed].extend(e for l in c_vars.values() for e in l)
        return lin_expr, lin_constr, lin_vars[0], lin_vars[1]
