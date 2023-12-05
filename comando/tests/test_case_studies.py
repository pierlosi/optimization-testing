"""Testing the case studies."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
from shutil import which

import pytest

from tests import skip_long, missing_module


# @pytest.mark.long
@pytest.mark.skipif(skip_long, reason='Skipping long running tests...')
@pytest.mark.skipif(which("gurobi") is None, reason="GUROBI is not installed")
@pytest.mark.skipif(which("baron") is None, reason="BARON is not installed")
def test_IES_case_study(run_in_tmpdir, clear_components):
    """Execute the industrial energy system case study.

    Takes about 2200 s on an i7-8700 CPU @ 3.20GHz, 3192 Mhz, 6 Core(s)
    """
    from examples.IES.case_study import run_IES_case_study
    with run_in_tmpdir:
        run_IES_case_study()
        import pickle
        from comando.utility import get_latest
        with open(get_latest('*_results.pickle'), 'rb') as f:
            obj_vals, dvs, ovs = pickle.load(f)
        result = obj_vals.values
        expected = [
            # TAC    GWI
            [559.95, 1262.29],
            [660.28, 1155.55],
            [860.67, 1031.97],
            [1021.75, 939.49],
            [1237.34, 865.32],
            [1876.55, 773.03],
            [2311.15, 678.36],
            [2960.84, 581.86]
        ]

        # 5 % error if based on lower bound (≈ 0.05/0.95)
        assert (abs(result - expected)/result).max() <= 5.3e-2


@pytest.mark.skipif(skip_long, reason='Skipping long running tests...')
@pytest.mark.skipif(missing_module("gurobipy"),
                    reason="gurobipy is not installed")
def test_DESTEST_case_study(run_in_tmpdir, clear_components):
    """Execute the DESTEST case study."""
    from examples.DESTEST.case_study import run_destest_case_study
    with run_in_tmpdir:
        run_destest_case_study()


@pytest.mark.skipif(skip_long, reason='Skipping long running tests...')
def test_BDR_case_study(run_in_tmpdir, clear_components):
    """Execute the building demand response case study.

    Takes about 1 s on a 2.8 GHz Quad-Core Intel Core i7
    """
    from examples.BDR.case_study import run_BDR_case_study
    from pyomo.opt.parallel.manager import ActionManagerError
    with run_in_tmpdir:
        try:
            assert run_BDR_case_study() == pytest.approx(24.8975)
        except ActionManagerError as e:
            if '[]' in str(e):
                pytest.xfail('NEOS is unavailable')


@pytest.mark.skipif(skip_long, reason='Skipping long running tests...')
@pytest.mark.skipif(missing_module('symengine'),
                    reason="symengine is not installed")
@pytest.mark.skipif(missing_module('maingopy') and missing_module('pymaingo'),
                    reason="maingopy/pymaingo is not installed")
@pytest.mark.skipif(which("cplex") is None, reason="CPLEX is not installed")
@pytest.mark.skipif(which("baron") is None, reason="BARON is not installed")
def test_ORC_case_study(run_in_tmpdir, clear_components, reset_backend):
    """Execute the organic Rankine cycle case study."""
    from examples.ORC.case_study import run_ORC_case_study
    with run_in_tmpdir:
        res = run_ORC_case_study()
        # assert res[0]['P_net'] == pytest.approx(14.134e6) == res[1]['P_net']
        assert res[2]['P_net'] == pytest.approx(16.485e6, rel=1e-3)
