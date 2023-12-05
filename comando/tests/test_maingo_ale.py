"""Tests for the MAiNGO ALE interface."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import os
from shutil import which

import pytest

import comando
from comando.interfaces.maingo_ale import solve
from comando import exp, tanh
from comando.utility import define_function


@pytest.mark.skipif(which("MAiNGO") is None, reason="MAiNGO is not installed")
@pytest.mark.parametrize('scenarios, timesteps',  # for test_problem
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_maingo_solve(test_problem, run_in_tmpdir):
    tp = test_problem

    # testing the inclusion of user defined functions into Maingo by adding a
    # trivial constraint which includes MAiNGO's lb_func
    lb_func = define_function('lb_func', lambda x, lb: comando.Max(x, lb))
    x = tp['x']  # getting variable x ∈ [-1.5, 1.5]
    # adding a trivial constraint (always satisfied)
    expr = (2 + x) ** 4
    tp.constraints['c3'] = expr / lb_func(expr, comando.EPS) <= 1
    # another trivial constraint to test tanh implementation
    tp.constraints['c4'] = comando.Eq(1 - 2/(exp(2 * x) + 1), tanh(x))

    name = tp.name
    with run_in_tmpdir:
        ret = solve(tp, reuse=False, epsilonA=1e-5)
        assert ret == 0
        file_name = f'{name}.ale'
        assert os.path.isfile(file_name)
        with open(file_name, 'r') as f:
            content = f.read()
        assert 'tanh' in content
        assert 'lb_func' in content
        settings_name = f'{name}_Settings.txt'
        assert os.path.isfile(settings_name)
        with open(settings_name, 'r') as f:
            assert f.readline().split() == ['epsilonA', '1e-05']
    assert x.value == pytest.approx(1, abs=1e-3)
    for yi in tp['y']:
        assert yi.value == pytest.approx(1, abs=1e-3)
    assert comando.utility.evaluate(tp.objective) \
        == pytest.approx(0, abs=1e-5)


