"""Collection of routines for visualizing data, expressions and results."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
from matplotlib import pyplot as plt
import numpy as np

import comando
from comando.utility import get_vars, get_pars, indexed


def convergence_plot(times, lb_data, ub_data, color='b', label=''):
    """Make a simple performance plot."""
    plt.plot(times, lb_data, linewidth=4, color=color, label=label)
    plt.plot(times, ub_data, linewidth=4, color=color)


def plot_expr(expr, label=None, axis=0, fig=None, cutoff=1E10, x=None, kind='surface', **kwargs):
    """Plot an expression."""
    if expr.is_Relational:
        lhs_vars = get_vars(expr.lhs)
        rhs_vars = get_vars(expr.rhs)
        if len(lhs_vars) != 1 or len(rhs_vars) != 1:
            raise NotImplementedError('Can only plot univariate relational '
                                      'expressions!')
        lhs_var = lhs_vars.pop()
        rhs_var = rhs_vars.pop()
        if lhs_var is rhs_var:  # same var, use x axis
            plot_expr(expr.rhs)
            plot_expr(expr.lhs)
        else:  # different vars, use y1 and y2 axes
            plot_expr(expr.rhs, label, 1, fig, cutoff)
            plot_expr(expr.lhs, label, 2, fig, cutoff)
    else:
        vars = get_vars(expr)
        n_vars = len(vars)
        if n_vars == 2:
            # Make data.
            if x is not None:  # Axis order is specified
                vars = [*vars]
                if vars[0].name == x:
                    x, y = vars
                elif vars[1].name == x:
                    y, x = vars
                else:
                    raise ValueError(f'No variable with name {x}!')
            else:  # Axis order is determined by free_symbols
                x = vars.pop()
                y = vars.pop()
            lb = float(max(-cutoff, x._bounds[0]))
            ub = float(min(cutoff, x._bounds[1]))
            X = np.array(range(int(lb), int(ub) + 1)) if x.is_integer \
                else np.linspace(lb, ub, 20)
            lb = float(max(-cutoff, y._bounds[0]))
            ub = float(min(cutoff, y._bounds[1]))
            Y = np.array(range(int(lb), int(ub) + 1)) if y.is_integer \
                else np.linspace(lb, ub, 20)

            # TODO: Taking the mean as fallback here might be confusing...
            data = {p: p.value.mean() if indexed(p) else p.value
                    for p in get_pars(expr)}
            func = comando.utility.lambdify(expr.subs(data), [x, y])
            XX, YY = np.meshgrid(X, Y)
            Z = func(XX, YY)

            if fig is None:
                fig = plt.figure()
            if kind == 'surface':
                from mpl_toolkits.mplot3d import Axes3D
                ax = fig.gca(projection='3d')
                plt.surf = ax.plot_surface(XX, YY, Z, label=label, **kwargs)
            elif kind == 'contour':
                ax = fig.gca()
                CS = ax.contour(X, Y, Z, **kwargs)
                ax.clabel(CS, inline=1, fontsize=10)
            else:
                raise ValueError(f'"{kind}" is not a recognized plot kind')

            if not label:
                label = str(expr)
            plt.title(label)
            plt.xlabel(str(x))
            plt.ylabel(str(y))
            plt.show()
        elif n_vars == 1:
            x = vars.pop()
            lb = float(max(-cutoff, x._bounds[0]))
            ub = float(min(cutoff, x._bounds[1]))
            X = np.array(range(int(lb), int(ub) + 1)) if x.is_integer \
                else np.linspace(lb, ub)
            data = {p: p.value.mean() if indexed(p) else p.value
                    for p in get_pars(expr)}

            # NOTE: More intuitive variant fails for Min and Max functions
            # _expr = expr.subs(data)
            # func = comando.utility.lambdify(_expr)
            # F = func(X)

            F = np.empty(X.shape)
            for i, val in enumerate(X):
                data[x] = val
                try:
                    F[i] = float(expr.subs(data))
                except TypeError:
                    F[i] = float('nan')

            if not label:
                label = str(expr)
            plt.title(label)
            if axis == 0:
                plt.plot(X, F, label=label, **kwargs)
                plt.xlabel(str(x))
            elif axis == 1:
                plt.plot(F, X, label=label, **kwargs)
                plt.ylabel(str(x))
            else:
                plt.plot(F, X, label=label, **kwargs)
                plt.gca().twinx().set_ylabel(str(x))
        elif n_vars == 0:
            if not label:
                label = str(expr)
            plt.title(label)
            plt.axhline(float(expr), label=label)
        else:
            raise NotImplementedError(f'Cannot plot expression with {n_vars} '
                                      'vars!')


def plot_incidence(exprs, vars=None):
    """Visualize which variables are contained in which expression."""
    if vars is None:  # If no vars are specified use all, sorted by name
        vars = sorted(set.union(*(get_vars(expr) for expr in exprs)), key=str)
    mat = np.empty((len(exprs), len(vars)))
    for i, expr in enumerate(exprs):
        mat[i, :] = [var in expr.free_symbols for var in vars]
    plt.imshow(mat, cmap='Greys')
    plt.xticks(range(len(vars)), [var.name for var in vars])
    plt.xticks(rotation=90)
    plt.yticks(range(len(exprs)), [str(expr) for expr in exprs])
