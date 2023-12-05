"""Tests for the Pyomo interface."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
from shutil import which

import pytest

import comando
from comando.interfaces.pyomo import to_pyomo

# TODO:
# class TestPyomo(object):
#     """Test the pyomo interface."""
#
#     def test_model_creation(self, test_problem):
#         self.m = to_pyomo(test_problem)
#
#     def test_solution_with_couenne(self):
#         m.solve('couenne')
#
#     def test_solution_with_baron(self):
#         m.solve('baron')


@pytest.mark.parametrize('scenarios, timesteps',  # for test_problem
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_pyomo_problem(test_problem):
    """Test the pyomo interface."""
    m = to_pyomo(test_problem)

    if which("baron") is None:
        pytest.skip("BARON is not installed")
    m.solve('baron', options={'epsa': 1e-9})
    assert test_problem['x'].value == pytest.approx(1)
    for i in test_problem.index:
        assert test_problem['y'][i].value == pytest.approx(1)
    assert comando.utility.evaluate(test_problem.objective) \
        == pytest.approx(0)


def test_pyomo_voll_problem(det_voll_problem):
    """Test the pyomo interface."""
    m = to_pyomo(det_voll_problem)
    assert len(m.x) + len(m.y) == det_voll_problem.num_vars
    assert len(m.constraints) == det_voll_problem.num_cons
    # TODO: Actually test whether the resulting object works as expected!
