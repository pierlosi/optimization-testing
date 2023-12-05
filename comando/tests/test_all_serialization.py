"""Tests for the Serialization.

IMPORTANT NOTE:
Using pickle and other serialization modules that result in the execution of
code are inherently unsafe, also see the pickle and dill documentation or
https://www.benfrederickson.com/dont-pickle-your-data/.

In the future we should therefore provide our own serialization, e.g., based on
JSON.
"""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
from contextlib import contextmanager
import pickle
import os

import pytest

import comando

try:
    import dill
    DILL_FOUND = True
except ModuleNotFoundError:
    DILL_FOUND = False
# import json


@contextmanager
def cwd(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


def compare(old, new):
    for attr in ['name', 'timesteps', 'data', 'design', 'operation']:
        comp = getattr(old, attr) == getattr(new, attr)
        try:
            assert comp
        except ValueError:
            assert all(comp)
    assert all(old.data == new.data)
    assert all(old.data == new.data.getter())
    # with pytest.xfail("Pandas doesn't implement __eq__ correctly"):
    #     assert all(old.data.getter() == new.data)
    assert all(old.data.getter() == new.data.getter())

    assert all(old.data != new.data)
    assert all(old.data != new.data.getter())
    # with pytest.xfail("Pandas doesn't implement __eq__ correctly"):
    #     assert all(old.data.getter() != new.data)
    assert all(old.data.getter() != new.data.getter())


def test_pickle_parameter():
    p = comando.Parameter('p', value=3.141)
    s = pickle.dumps(p)
    pl = pickle.loads(s)
    assert pl == p


def test_pickle_variable():
    x = comando.Variable('x', bounds=(0.2, 1.2))
    s = pickle.dumps(x)
    xl = pickle.loads(s)
    assert xl == x


def test_pickle_variable_vector():
    y = comando.VariableVector('y', bounds=(0.2, 1.2))
    s = pickle.dumps(y)
    yl = pickle.loads(s)
    assert yl == y


@pytest.mark.parametrize('scenarios, timesteps',  # for test_problem
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_pickle(test_problem, tmpdir):
    with cwd(tmpdir):
        with open('test.pickle', 'wb') as f:
            pickle.dump(test_problem, f, -1)
        with open('test.pickle', 'rb') as f:
            loaded_problem = pickle.load(f)
    compare(test_problem, loaded_problem)


@pytest.mark.skipif(not DILL_FOUND, reason="dill is not installed")
@pytest.mark.parametrize('scenarios, timesteps',  # for test_problem
                         [(None, (['t1'], 1)),
                          (['s1', 's2'], (['t1'], 1)),
                          (['s1', 's2'], None)])
def test_dill(test_problem, tmpdir):
    with cwd(tmpdir):
        with open('test.dill', 'wb') as f:
            dill.dump(test_problem, f, -1)
        with open('test.dill', 'rb') as f:
            loaded_problem = dill.load(f)
    compare(test_problem, loaded_problem)
    # dill.dump_session("test.sess")
    # dill.load_session("test.sess")

# def test_json():
#     # TypeError: Object of type Problem is not JSON serializable
#     with open('test.json', 'wb') as f:
#         json.dump(P, f)
