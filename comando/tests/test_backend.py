"""Tests for different backends."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import comando


def test_backend_setting():
    import sympy
    comando.set_backend('sympy')
    assert comando.get_backend() == sympy
    assert sympy.Symbol in comando.Variable.__mro__

    try:
        import symengine
        comando.set_backend('symengine')
        assert comando.get_backend() == symengine
        assert sympy.Symbol not in comando.Variable.__mro__
        assert symengine.Symbol in comando.Variable.__mro__
        with comando.set_backend('sympy'):  # Temporarily switching to sympy
            assert comando.get_backend() == sympy
            assert sympy.Symbol in comando.Variable.__mro__
            assert symengine.Symbol not in comando.Variable.__mro__
        assert comando.get_backend() == symengine
        assert sympy.Symbol not in comando.Variable.__mro__
        assert symengine.Symbol in comando.Variable.__mro__
    except ModuleNotFoundError:
        pass
    finally:  # cleanup
        comando.set_backend('sympy')
