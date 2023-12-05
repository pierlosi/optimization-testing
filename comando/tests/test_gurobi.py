"""Tests for the GUROBI interface."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import pytest
import comando

from contextlib import contextmanager
import os
from shutil import which


@contextmanager
def cwd(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


@pytest.mark.skipif(which("gurobi") is None, reason="GUROBI is not installed")
@pytest.mark.parametrize('scenarios, timesteps',
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_gurobi_solve(scenarios, timesteps):
    try:
        from comando.interfaces.gurobi import to_gurobi
    except ModuleNotFoundError as e:
        if 'gurobipy' in str(e):
            print('Module gurobipy cannot be found, you may have forgotten to '
                  'install it!')
            pytest.xfail('Module gurobipy cannot be found')
        raise
    from comando.utility import silence

    x = comando.Variable('x', bounds=(-1.5, 1.5))
    x_squared = comando.Variable('x_squared', bounds=(-2.25, 2.25))
    y = comando.VariableVector('y', bounds=(-0.5, 2.5))
    # y_squared = comando.VariableVector('y_squared', bounds=(-0.25, 6.25))
    p = comando.Parameter('p', 2)
    do = 1 - 2 * x + x_squared
    oo = 100 * (y ** 2 - 2 * y * x_squared + x_squared ** 2)
    constraints = {'c1': x * x_squared - 3 * x_squared + 3 * x - y <= 0,
                   'c2': x + y - p <= 0,
                   'x_squared': comando.Eq(x_squared, x ** 2)}
    P = comando.Problem(do, oo, constraints, timesteps=timesteps,
                        scenarios=scenarios, name='Rosenbrok')
    gm = to_gurobi(P)
    with silence():
        gm.solve(NonConvex=2)

    assert x.value == pytest.approx(1)
    for yi in y:
        assert yi.value == pytest.approx(1)
    assert comando.utility.evaluate(P.objective) == pytest.approx(0, abs=1e-6)
