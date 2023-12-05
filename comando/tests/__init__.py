"""Automatic tests for the COMANDO Framework to be used with pytest."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import sys
import os
import re

sys.path.insert(0, os.path.abspath('..'))


def cleanup_problem_files():
    """Remove any problem files that may have been created by tests."""
    for f in os.listdir('.'):
        if re.search("tmp.*(pyomo.nl|neos.sol|neos.log)", f):
            os.remove(os.path.join('.', f))


try:
    skip_long = not os.environ["DO_LONG_TESTS"]
except KeyError:
    skip_long = True


def missing_module(module_name):
    """Check whether the module with the given name is missing."""
    import importlib
    try:
        importlib.import_module(module_name)
        return False
    except ModuleNotFoundError:
        return True
