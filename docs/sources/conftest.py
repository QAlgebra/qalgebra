"""This file is automatically executed by pytest when testing anything in the
docs folder"""
import pytest

import qalgebra


@pytest.fixture(autouse=True)
def set_doctest_env(doctest_namespace):
    """Inject package itself into doctest namespace.

    This is so we don't need

    .. doctest::

        >>> import qalgebra

    in any doctests
    """
    doctest_namespace['qalgebra'] = qalgebra


def pytest_collectstart(collector):
    """Ignore stderr and javascript output when verifying notebooks.

    Works around a test failure on Travis/AppVeyor (See Travis Build 56).

    https://nbval.readthedocs.io/en/latest/#Skipping-certain-output-types
    """
    if collector.__class__.__name__ == 'IPyNbFile':
        collector.skip_compare += (
            'application/javascript',
            'stderr',
            'application/vnd.jupyter.widget-view+json',
        )
