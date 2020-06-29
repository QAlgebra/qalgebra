"""Test for the symbolic_equation package.

This is maintained as an external package, but we want to test that it
integrates well with qalgebra
"""
import pytest
import sympy
from symbolic_equation import Eq
from sympy.core.sympify import SympifyError

from qalgebra import (
    Create,
    Destroy,
    IdentityOperator,
    OperatorSymbol,
    ZeroOperator,
    latex,
)


# These only cover things not already coveraged in the doctest


def test_apply_to_lhs():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq = eq0.apply_to_lhs(lambda expr: expr + E0).tag('new')
    assert eq.lhs == H_0 + E0
    assert eq.rhs == eq0.rhs
    assert eq._tag == 'new'


def test_apply_mtd():
    H_0 = OperatorSymbol('H_0', hs=0)
    H = OperatorSymbol('H', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq = eq0.apply('substitute', {H_0: H, E0: 0}).tag('new')
    assert eq.lhs == H
    assert eq.rhs == ω * Create(hs=0) * Destroy(hs=0)
    assert eq._tag == 'new'


def test_eq_copy():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq = eq0.copy()
    assert eq == eq0
    assert eq is not eq0


def test_eq_add_const():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq = eq0 + E0
    assert eq.lhs == H_0 + E0
    assert eq.rhs == eq0.rhs + E0
    assert eq._tag is None


def test_eq_mult_const():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq = 2 * eq0
    assert eq == eq0 * 2
    assert eq.lhs == 2 * eq0.lhs
    assert eq.rhs == 2 * eq0.rhs
    assert eq._tag is None


def test_eq_div_const():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq = eq0 / 2
    assert eq.lhs == eq0.lhs / 2
    assert eq.rhs == eq0.rhs / 2
    assert eq._tag is None


def test_eq_equals_const():
    H_0 = OperatorSymbol('H_0', hs=0)
    eq0 = Eq(H_0, IdentityOperator)
    assert eq0 - 1 == ZeroOperator


def test_eq_sub_eq():
    ω, E0 = sympy.symbols('omega, E_0')
    H_0 = OperatorSymbol('H_0', hs=0)
    H_1 = OperatorSymbol('H_1', hs=0)
    mu = OperatorSymbol('mu', hs=0)
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq1 = Eq(H_1, mu + E0, tag='1')
    eq = eq0 - eq1
    assert eq.lhs == H_0 - H_1
    assert eq.rhs == ω * Create(hs=0) * Destroy(hs=0) - mu
    assert eq._tag is None


def test_eq_sub_const():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq = eq0 - E0
    assert eq.lhs == H_0 - E0
    assert eq.rhs == ω * Create(hs=0) * Destroy(hs=0)
    assert eq._tag is None


def test_repr_latex():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    Eq.latex_renderer = staticmethod(latex)
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    repr1 = eq0._repr_latex_()
    repr2 = latex(eq0)
    assert repr1 == repr2


def test_eq_str():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    assert str(eq0) == "%s = %s    (0)" % (str(eq0.lhs), str(eq0.rhs))


def test_eq_repr():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    assert repr(eq0) == "%s = %s    (0)" % (repr(eq0.lhs), repr(eq0.rhs))


def test_no_sympify():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    with pytest.raises(SympifyError):
        sympy.sympify(eq0)


def test_eq_substitute():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')
    eq1 = eq0.apply('substitute', {E0: 0}).reset()
    eq2 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0))
    assert eq1 == eq2


def test_unchanged_apply():
    H_0 = OperatorSymbol('H_0', hs=0)
    ω, E0 = sympy.symbols('omega, E_0')
    eq0 = Eq(H_0, ω * Create(hs=0) * Destroy(hs=0) + E0, tag='0')

    assert eq0.apply(lambda s: s.expand()).reset() == eq0
    assert eq0.apply(lambda s: s.expand()) == eq0
    assert eq0.apply(lambda s: s.expand())._lhs is None

    assert eq0.apply('expand').reset() == eq0
    assert eq0.apply('expand') == eq0
    assert eq0.apply('expand')._lhs is None
