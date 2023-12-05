"""Tests for the MAiNGO API interface."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import pytest
import comando

from tests import missing_module


@pytest.mark.skipif(missing_module('maingopy') and missing_module('pymaingo'),
                    reason="maingopy/pymaingo is not installed")
@pytest.mark.parametrize('scenarios, timesteps',  # for test_problem
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_maingo_solve(test_problem):
    from comando.interfaces.maingo_api import MaingoProblem
    from comando.utility import define_function

    tp = test_problem

    # testing the inclusion of user defined functions into Maingo by adding a
    # trivial constraint which includes MAiNGO's lb_func
    lb_func = define_function('lb_func', lambda x, lb: comando.Max(x, lb))
    x = tp['x']  # getting variable x ∈ [-1.5, 1.5]
    # adding a trivial constraint (always satisfied)
    tp.constraints['c3'] = (2 + x) / lb_func(2 + x, comando.EPS) <= 1

    mp = MaingoProblem(tp)
    solver, ret = mp.solve(epsilonA=1e-5)  # ), outstreamVerbosity=1)
    assert ret == comando.interfaces.maingo_api.GLOBALLY_OPTIMAL
    assert solver.get_cpu_solution_time() \
        <= solver.get_wallclock_solution_time()
    assert solver.get_final_abs_gap() < 1e-5
    assert solver.get_final_LBD() == 0
    assert solver.get_final_rel_gap() == 1  # as LBD is 0 this is  UBD / UBD
    assert solver.get_iterations() >= solver.get_LBP_count()
    assert solver.get_iterations() >= solver.get_UBP_count()
    assert solver.get_iterations() >= solver.get_max_nodes_in_memory()
    sp = solver.get_solution_point()
    assert sp[0] == pytest.approx(1, abs=1e-3) == tp['x'].value
    for yi, spi in zip(tp['y'], sp[1:]):
        assert spi == pytest.approx(1, abs=1e-3) == yi.value
    vals = solver.evaluate_model_at_solution_point()
    assert vals[0] == pytest.approx(0, abs=1e-3) == tp.objective.value
    exp_vals = iter(vals[1:])
    for con_val, exp_val in zip(tp.constraints['c1'].lhs.value, exp_vals):
        assert con_val == pytest.approx(exp_val)
    for con_val, exp_val in zip(tp.constraints['c2'].lhs.value, exp_vals):
        assert con_val == pytest.approx(exp_val)
    assert len(solver.evaluate_additional_outputs_at_solution_point()) == 0
    assert solver.get_objective_value() == pytest.approx(0, abs=1e-5) \
        == tp.objective.value
