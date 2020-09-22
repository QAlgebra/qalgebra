"""Test indexed sums over operators"""
import pytest
from sympy import IndexedBase, symbols

from qalgebra import (
    FockIndex,
    IdxSym,
    IndexOverFockSpace,
    IndexOverList,
    KroneckerDelta,
    LocalSigma,
    LocalSpace,
    OperatorIndexedSum,
    OperatorSymbol,
    ScalarTimesOperator,
    StrLabel,
)


def test_pull_constfactor_from_operator_sum():
    """Test that a constant in pulled from a sum of operators."""
    i = IdxSym('i')
    alpha = symbols('alpha')

    def A(i, j):
        return OperatorSymbol(StrLabel(IndexedBase('A')[i, j]), hs=0)

    term = alpha * A(i, i)
    range_i = IndexOverList(i, (1, 2))
    expr = OperatorIndexedSum.create(term, ranges=(range_i,))
    assert isinstance(expr, ScalarTimesOperator)
    assert expr.coeff == alpha
    assert isinstance(expr.term, OperatorIndexedSum)


def test_operator_kronecker_sum():
    """Test that Kronecker delta are eliminiated from indexed sums over
    operators"""
    i = IdxSym('i')
    j = IdxSym('j')
    alpha = symbols('alpha')
    delta_ij = KroneckerDelta(i, j)
    delta_0i = KroneckerDelta(0, i)
    delta_1j = KroneckerDelta(1, j)
    delta_0j = KroneckerDelta(0, j)
    delta_1i = KroneckerDelta(1, i)

    def A(i, j):
        return OperatorSymbol(StrLabel(IndexedBase('A')[i, j]), hs=0)

    term = delta_ij * A(i, j)
    sum = OperatorIndexedSum.create(
        term, ranges=(IndexOverList(i, (1, 2)), IndexOverList(j, (1, 2)))
    )
    assert sum == OperatorIndexedSum.create(
        A(i, i), ranges=(IndexOverList(i, (1, 2)),)
    )
    assert sum.doit() == (
        OperatorSymbol("A_11", hs=0) + OperatorSymbol("A_22", hs=0)
    )

    term = alpha * delta_ij * A(i, j)
    range_i = IndexOverList(i, (1, 2))
    range_j = IndexOverList(j, (1, 2))
    sum = OperatorIndexedSum.create(term, ranges=(range_i, range_j))
    assert isinstance(sum, ScalarTimesOperator)
    expected = alpha * OperatorIndexedSum.create(
        A(i, i), ranges=(IndexOverList(i, (1, 2)),)
    )
    assert sum == expected

    hs = LocalSpace('0', basis=('g', 'e'))
    i_range = IndexOverFockSpace(i, hs)
    j_range = IndexOverFockSpace(j, hs)
    sig_ij = LocalSigma(FockIndex(i), FockIndex(j), hs=hs)
    sig_0j = LocalSigma('g', FockIndex(j), hs=hs)
    sig_i1 = LocalSigma(FockIndex(i), 'e', hs=hs)

    term = delta_0i * delta_1j * sig_ij

    sum = OperatorIndexedSum.create(term, ranges=(i_range,))
    expected = delta_1j * sig_0j
    assert sum == expected

    sum = OperatorIndexedSum.create(term, ranges=(j_range,))
    expected = delta_0i * sig_i1
    assert sum == expected

    term = (delta_0i * delta_1j + delta_0j * delta_1i) * sig_ij
    sum = OperatorIndexedSum.create(term, ranges=(i_range, j_range))
    expected = LocalSigma('g', 'e', hs=hs) + LocalSigma('e', 'g', hs=hs)
    assert sum == expected
