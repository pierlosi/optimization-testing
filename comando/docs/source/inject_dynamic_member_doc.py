"""Code to inject the dynamic classes into the autodoc-generated .rst file."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHOR: Marco Langiu

dynamic_members = """
.. autoclass:: comando.Parameter
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: comando.Variable
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: comando.VariableVector
   :members:
   :undoc-members:
   :show-inheritance:
"""


def inject():
    """Inject the api documentation of dynamic class definitions."""
    import pathlib

    with open(pathlib.Path(__file__).parent / "api/comando.rst", "r") as f1:
        t1 = f1.readlines()

    if any('.. autoclass:: comando.Parameter' in line for line in t1):
        return  # Already inserted!

    lines = iter(t1)

    with open(pathlib.Path(__file__).parent / "api/comando.rst", "w") as f1:
        for line in lines:
            f1.write(line)
            if line.startswith('.. automodule:: comando'):
                f1.write(next(lines))
                f1.write(next(lines))
                f1.write(next(lines))
                f1.write(dynamic_members)
