"""Set up the environment for doctests.

This file is automatically evaluated by pytest. It ensures that we can write
doctests without distracting import statements in the doctest.
"""
import numpy
import pytest
import sympy

import qalgebra


@pytest.fixture(autouse=True)
def set_doctest_env(doctest_namespace):
    """Set up doctest environment."""
    doctest_namespace['numpy'] = numpy
    doctest_namespace['sympy'] = sympy
    doctest_namespace['qalgebra'] = qalgebra
    for attr in qalgebra.__all__:
        doctest_namespace[attr] = getattr(qalgebra, attr)
