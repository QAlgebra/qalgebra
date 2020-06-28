import pytest
from sympy import atan, symbols

from qalgebra.core.operator_algebra import LocalSigma, OperatorSymbol
from qalgebra.core.state_algebra import CoherentStateKet
from qalgebra.printing import (
    ascii,
    configure_printing,
    init_printing,
    latex,
    unicode,
)


def test_sympy_tex_cached():
    """Test that we can use the cache to change how sub-expressions of sympy
    are printed in tex"""
    a = symbols('a')
    A = OperatorSymbol("A", hs=1)
    expr = (a ** 2 / 2) * A

    assert latex(expr) == r'\frac{a^{2}}{2} \hat{A}^{(1)}'

    cache = {a: r'\alpha'}
    assert latex(expr, cache=cache) == r'\frac{\alpha^{2}}{2} \hat{A}^{(1)}'


def test_sympy_setting():
    """Test that we can pass settings to the sympy sub-printer"""
    x = symbols('a')
    A = OperatorSymbol("A", hs=1)
    expr = atan(x) * A
    assert latex(expr) == r'\operatorname{atan}{\left(a \right)} \hat{A}^{(1)}'
    assert (
        latex(expr, inv_trig_style='full')
        == r'\arctan{\left(a \right)} \hat{A}^{(1)}'
    )


def test_custom_options():
    """Test giving options to print routines or using configure_printing"""
    A = OperatorSymbol('A', hs=1)
    CNOT = OperatorSymbol('CNOT', hs=1)
    sig = LocalSigma(0, 1, hs=1)
    ket = CoherentStateKet(symbols('alpha'), hs=1)

    assert ascii(A) == r'A^(1)'
    assert ascii(A, show_hs_label=False) == 'A'
    with pytest.raises(TypeError) as exc_info:
        ascii(A, some_bogus_option=False)
    assert "not a valid setting" in str(exc_info.value)
    assert ascii(sig) == r'|0><1|^(1)'
    assert ascii(ket) == r'|alpha=alpha>^(1)'
    assert unicode(A) == r'Â⁽¹⁾'
    assert unicode(sig) == r'|0⟩⟨1|⁽¹⁾'
    assert unicode(ket) == r'|α=α⟩⁽¹⁾'
    assert latex(A) == r'\hat{A}^{(1)}'
    assert (
        latex(sig)
        == r'\left\lvert 0 \middle\rangle\!\middle\langle 1 \right\rvert^{(1)}'
    )
    assert latex(ket) == r'\left\lvert \alpha=\alpha \right\rangle^{(1)}'

    with configure_printing(
        unicode_op_hats=False, tex_op_macro=r'\Op{{{name}}}'
    ):
        assert unicode(A) == r'A⁽¹⁾'
        assert latex(A) == r'\Op{A}^{(1)}'

    with configure_printing(show_hs_label=False):
        assert ascii(A) == r'A'
        assert ascii(sig) == r'|0><1|'
        assert ascii(ket) == r'|alpha=alpha>'
        assert unicode(A) == r'Â'
        assert unicode(sig) == r'|0⟩⟨1|'
        assert unicode(ket) == r'|α=α⟩'
        assert latex(A) == r'\hat{A}'
        assert latex(A, show_hs_label=True) == r'\hat{A}^{(1)}'
        assert latex(A) == r'\hat{A}'
        assert (
            latex(sig)
            == r'\left\lvert 0 \middle\rangle\!\middle\langle 1 \right\rvert'
        )
        assert latex(ket) == r'\left\lvert \alpha=\alpha \right\rangle'

    assert latex(CNOT) == r'\text{CNOT}^{(1)}'
    with configure_printing(tex_textop_macro=r'\Op{{{name}}}'):
        assert latex(CNOT) == r'\Op{CNOT}^{(1)}'

    init_printing(show_hs_label=False)
    assert unicode(A) == r'Â'
    assert latex(A) == r'\hat{A}'
    with configure_printing(
        unicode_op_hats=False, tex_op_macro=r'\Op{{{name}}}'
    ):
        assert unicode(A) == r'A'
        assert latex(A) == r'\Op{A}'
        with configure_printing(tex_op_macro=r'\op{{{name}}}'):
            assert unicode(A) == r'A'
            assert latex(A) == r'\op{A}'
            with configure_printing(tex_use_braket=True):
                assert latex(sig) == r'\Ket{0}\!\Bra{1}'
    assert unicode(A) == r'Â'
    assert latex(A) == r'\hat{A}'
    init_printing(reset=True)

    assert ascii(A) == r'A^(1)'
    assert ascii(sig) == r'|0><1|^(1)'
    assert ascii(ket) == r'|alpha=alpha>^(1)'
    assert unicode(A) == r'Â⁽¹⁾'
    assert unicode(sig) == r'|0⟩⟨1|⁽¹⁾'
    assert unicode(ket) == r'|α=α⟩⁽¹⁾'
    assert latex(A) == r'\hat{A}^{(1)}'
    assert (
        latex(sig)
        == r'\left\lvert 0 \middle\rangle\!\middle\langle 1 \right\rvert^{(1)}'
    )
    assert latex(ket) == r'\left\lvert \alpha=\alpha \right\rangle^{(1)}'


def test_custom_repr():
    A = OperatorSymbol('A', hs=1)
    assert repr(A) in ['Â⁽¹⁾', 'A^(1)']
    init_printing(repr_format='srepr', reset=True)
    assert repr(A) == "OperatorSymbol('A', hs=LocalSpace('1'))"
    init_printing(reset=True)
    assert repr(A) in ['Â⁽¹⁾', 'A^(1)']
    with configure_printing(repr_format='srepr'):
        assert repr(A) == "OperatorSymbol('A', hs=LocalSpace('1'))"
    assert repr(A) in ['Â⁽¹⁾', 'A^(1)']


def test_exception_teardown():
    """Test that teardown works when breaking out due to an exception"""

    class ConfigurePrintingException(Exception):
        pass

    init_printing(show_hs_label=True, repr_format='ascii')
    try:
        with configure_printing(show_hs_label=False, repr_format='srepr'):
            raise ConfigurePrintingException
    except ConfigurePrintingException:
        A = OperatorSymbol('A', hs=1)
        assert repr(A) == 'A^(1)'
    finally:
        # Even if this failed we don't want to make a mess for other tests
        init_printing(reset=True)
