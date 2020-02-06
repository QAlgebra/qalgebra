import sympy

import qalgebra.core.operator_algebra
import qalgebra.library.fock_operators
from qalgebra.convert.to_sympy_matrix import convert_to_sympy_matrix
from qalgebra.core.hilbert_space_algebra import LocalSpace


def test_convert_to_sympy_matrix():
    N = 4
    Hil = LocalSpace('full', basis=range(N))

    Hil_q1 = LocalSpace('Q1', basis=range(2))
    Hil_q2 = LocalSpace('Q2', basis=range(2))

    expr = qalgebra.core.operator_algebra.IdentityOperator
    assert convert_to_sympy_matrix(expr, Hil) == sympy.eye(N)

    expr = qalgebra.core.operator_algebra.ZeroOperator
    assert convert_to_sympy_matrix(expr, Hil) == 0

    expr = qalgebra.library.fock_operators.Create(hs=Hil)
    assert convert_to_sympy_matrix(expr, Hil) == sympy.Matrix(
        [
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [0, sympy.sqrt(2), 0, 0],
            [0, 0, sympy.sqrt(3), 0],
        ]
    )

    expr = qalgebra.library.fock_operators.Destroy(hs=Hil)
    assert convert_to_sympy_matrix(expr, Hil) == sympy.Matrix(
        [
            [0, 1, 0, 0],
            [0, 0, sympy.sqrt(2), 0],
            [0, 0, 0, sympy.sqrt(3)],
            [0, 0, 0, 0],
        ]
    )

    phi = sympy.symbols('phi', real=True)
    expr = qalgebra.library.fock_operators.Phase(phi, hs=Hil)
    assert convert_to_sympy_matrix(expr, Hil) == sympy.Matrix(
        [
            [1, 0, 0, 0],
            [0, sympy.exp(sympy.I * phi), 0, 0],
            [0, 0, sympy.exp(2 * sympy.I * phi), 0],
            [0, 0, 0, sympy.exp(3 * sympy.I * phi)],
        ]
    )

    expr = qalgebra.core.operator_algebra.LocalSigma(1, 3, hs=Hil)
    assert convert_to_sympy_matrix(expr, Hil) == sympy.Matrix(
        [[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]]
    )

    expr = qalgebra.library.fock_operators.Create(
        hs=Hil
    ) + qalgebra.library.fock_operators.Destroy(hs=Hil)
    assert convert_to_sympy_matrix(expr, Hil) == sympy.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, sympy.sqrt(2), 0],
            [0, sympy.sqrt(2), 0, sympy.sqrt(3)],
            [0, 0, sympy.sqrt(3), 0],
        ]
    )

    expr = qalgebra.library.fock_operators.Create(
        hs=Hil
    ) * qalgebra.library.fock_operators.Destroy(hs=Hil)
    assert convert_to_sympy_matrix(expr, Hil) == sympy.Matrix(
        [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3]]
    )

    w1, w2 = sympy.symbols('omega_1, omega_2', real=True)
    expr = w1 * (
        qalgebra.library.fock_operators.Create(hs=Hil_q1)
        * qalgebra.library.fock_operators.Destroy(hs=Hil_q1)
    ) + w2 * (
        qalgebra.library.fock_operators.Create(hs=Hil_q2)
        * qalgebra.library.fock_operators.Destroy(hs=Hil_q2)
    )
    assert convert_to_sympy_matrix(expr, expr.space) == sympy.Matrix(
        [[0, 0, 0, 0], [0, w2, 0, 0], [0, 0, w1, 0], [0, 0, 0, w1 + w2]]
    )

    J = sympy.symbols('J', real=True)
    expr = J * (
        qalgebra.library.fock_operators.Create(hs=Hil_q1)
        * qalgebra.library.fock_operators.Destroy(hs=Hil_q2)
        + qalgebra.library.fock_operators.Destroy(hs=Hil_q1)
        * qalgebra.library.fock_operators.Create(hs=Hil_q2)
    )
    assert convert_to_sympy_matrix(expr, expr.space) == sympy.Matrix(
        [[0, 0, 0, 0], [0, 0, J, 0], [0, J, 0, 0], [0, 0, 0, 0]]
    )

    expr = qalgebra.core.operator_algebra.Adjoint(
        qalgebra.library.fock_operators.Create(hs=Hil)
    )
    expr2 = qalgebra.library.fock_operators.Destroy(hs=Hil)
    assert convert_to_sympy_matrix(
        expr, expr.space
    ) == convert_to_sympy_matrix(expr2, expr2.space)
