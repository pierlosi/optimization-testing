"""Tests for the BARON interface."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import os
from shutil import which

import pytest

import comando


@pytest.mark.skipif(which("baron") is None, reason="BARON is not installed")
@pytest.mark.parametrize('scenarios, timesteps',  # for test_problem
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_baron_solve(test_problem, run_in_tmpdir):
    from comando.interfaces.baron import solve

    name = test_problem.name
    with run_in_tmpdir:
        solve(test_problem, reuse=False, MaxTime=1, times=False, epsa=1e-9)
        assert os.path.isfile(f'{name}.res.lst')
        assert os.path.isfile(f'{name}.sum.lst')
        assert not os.path.isfile(f'{name}.tim.lst')  # explicitly turned off!
    assert test_problem['x'].value == pytest.approx(1)
    for yi in test_problem['y']:
        assert yi.value == pytest.approx(1)
    assert comando.utility.evaluate(test_problem.objective) \
        == pytest.approx(0)
