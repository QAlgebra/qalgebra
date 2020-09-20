from qalgebra.core.abstract_algebra import _apply_rules
from qalgebra.core.exceptions import CannotSimplify
from qalgebra.core.hilbert_space_algebra import LocalSpace
from qalgebra.core.operator_algebra import (
    Operator,
    OperatorPlus,
    OperatorSymbol,
    OperatorTimes,
    ScalarTimesOperator,
)
from qalgebra.pattern_matching import pattern, pattern_head, wc
from qalgebra.printing import srepr
from qalgebra.toolbox.core import temporary_rules


def test_simplify():
    """Test simplification of expr according to manual rules"""
    h1 = LocalSpace("h1")
    a = OperatorSymbol("a", hs=h1)
    b = OperatorSymbol("b", hs=h1)
    c = OperatorSymbol("c", hs=h1)
    d = OperatorSymbol("d", hs=h1)

    expr = 2 * (a * b * c - b * c * a)

    A_ = wc('A', head=Operator)
    B_ = wc('B', head=Operator)
    C_ = wc('C', head=Operator)

    def b_times_c_equal_d(B, C):
        if B.label == 'b' and C.label == 'c':
            return d
        else:
            raise CannotSimplify

    with temporary_rules(OperatorTimes):
        OperatorTimes.add_rule(
            'extra', pattern_head(B_, C_), b_times_c_equal_d
        )
        new_expr = expr.rebuild()

    commutator_rule = (
        pattern(
            OperatorPlus,
            pattern(OperatorTimes, A_, B_),
            pattern(ScalarTimesOperator, -1, pattern(OperatorTimes, B_, A_)),
        ),
        lambda A, B: OperatorSymbol(
            "Commut%s%s" % (A.label.upper(), B.label.upper()), hs=A.space
        ),
    )
    assert commutator_rule[0].match(new_expr.term)

    with temporary_rules(OperatorTimes):
        OperatorTimes.add_rule(
            'extra', pattern_head(B_, C_), b_times_c_equal_d
        )
        new_expr = _apply_rules(expr, [commutator_rule])
    assert (
        srepr(new_expr)
        == "ScalarTimesOperator(ScalarValue(2), OperatorSymbol('CommutAD', "
        "hs=LocalSpace('h1')))"
    )
