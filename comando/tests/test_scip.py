"""Tests for the SCIP interface."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import pytest
import comando

from shutil import which


@pytest.mark.skipif(which("scip") is None, reason="SCIP is not installed")
@pytest.mark.parametrize('scenarios, timesteps',  # for test_problem
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_scip_solve(test_problem, run_in_tmpdir, scenarios):
    try:
        from comando.interfaces.scip import ScipProblem
    except ModuleNotFoundError as e:
        if 'pyscipopt' in str(e):
            print('Module pyscipopt cannot be found, you may have forgotten to'
                  ' install it!')
            pytest.xfail('Module pyscipopt cannot be found')
        raise
    sp = ScipProblem(test_problem)
    with run_in_tmpdir:
        # TODO: Why do we need to set such extreme tolerances to get reasonable
        #       results with SCIP?
        sp.solve(numerics_epsilon=1e-10,
                 numerics_feastol=1e-10,
                 numerics_dualfeastol=1e-10)

    assert test_problem['x'].value == pytest.approx(1, abs=1e-5)
    for yi in test_problem['y']:
        assert yi.value == pytest.approx(1, abs=1e-5)
    assert comando.utility.evaluate(test_problem.objective) \
        == pytest.approx(0, abs=1e-6)
