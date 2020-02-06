from sympy import S, symbols

from qalgebra.core.hilbert_space_algebra import LocalSpace
from qalgebra.core.operator_algebra import (
    IdentityOperator,
    OperatorPlus,
    OperatorSymbol,
    ZeroOperator,
)


def test_op_plus_scalar():
    """Test the we can add a scalar to an operator"""
    hs = LocalSpace("0")
    A = OperatorSymbol('A', hs=hs)
    alpha = symbols('alpha')
    assert A + 0 == A
    assert OperatorPlus.create(A, 0) == A
    assert 0 + A == A
    assert A + S.Zero == A
    assert ZeroOperator + S.Zero == ZeroOperator
    assert OperatorPlus.create(ZeroOperator, S.Zero) == ZeroOperator
    assert A + S.One == A + 1
    assert A + alpha == OperatorPlus(alpha * IdentityOperator, A)
    assert OperatorPlus.create(alpha, A) == OperatorPlus(
        alpha * IdentityOperator, A
    )
