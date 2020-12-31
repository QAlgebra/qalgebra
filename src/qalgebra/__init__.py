"""Main QAlgebra package

The :mod:`qalgebra` package exposes all of QAlgebra's functionality for easy
interactive or programmative use.

For interactive usage, the package should be initialized as follows::

    >>> import qalgebra
    >>> qalgebra.init_printing()

"""
import importlib

import qalgebra._flat_api_tools
import qalgebra._rules


__doc__ += qalgebra._flat_api_tools.__doc__

__all__ = []  # will be extended by _import_submodules

__known_refs__ = {}

__version__ = '0.2.0'


def _git_version():
    """If installed with 'pip installe -e .' from inside a git repo, the
    current git revision as a string"""

    import os
    import subprocess

    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        FNULL = open(os.devnull, 'w')
        cwd = os.path.dirname(os.path.realpath(__file__))
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=FNULL, env=env, cwd=cwd
        )
        out = proc.communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        return out.strip().decode('ascii')
    except OSError:
        return "unknown"


__git_version__ = _git_version()


def init_algebra(*, default_hs_cls='LocalSpace'):
    """Initialize the algebra system.

    Args:
        default_hs_cls (str): The name of the :class:`.LocalSpace` subclass
            that should be used when implicitly creating Hilbert spaces, e.g.
            in :class:`.OperatorSymbol`

    """
    from qalgebra.core.abstract_quantum_algebra import QuantumExpression
    from qalgebra.core.hilbert_space_algebra import LocalSpace

    default_hs_cls = getattr(
        importlib.import_module('qalgebra'), default_hs_cls
    )
    if issubclass(default_hs_cls, LocalSpace):
        QuantumExpression._default_hs_cls = default_hs_cls
    else:
        raise TypeError("default_hs_cls must be a subclass of LocalSpace")
    # TODO: init_algebra should eventually control which rules QAlgebra uses.


# dynamic initialization

qalgebra._flat_api_tools._import_submodules(__all__, __path__, __name__)

qalgebra._rules._algebraic_rules()
init_algebra()
