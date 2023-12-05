"""Tests for the GAMS interface."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
from shutil import which

import pytest
import comando


@pytest.mark.skipif(which("gams") is None, reason="GAMS is not installed")
@pytest.mark.parametrize('scenarios, timesteps',  # for test_problem
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_gams_solve(test_problem, run_in_tmpdir):
    from comando.interfaces.gams import solve
    from comando.utility import silence

    solver = 'couenne' if which("baron") is None else 'baron'
    with silence():
        with run_in_tmpdir:
            res = solve(test_problem, model_type='NLP', NLP=solver, silent=True)
            if res == 7:
                pytest.xfail("GAMS couldn't find a valid license!")

    assert test_problem['x'].value == pytest.approx(1)
    for yi in test_problem['y']:
        assert yi.value == pytest.approx(1)
    assert comando.utility.evaluate(test_problem.objective) \
        == pytest.approx(0)
