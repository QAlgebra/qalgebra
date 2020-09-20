"""Conversion of QAlgebra expressions to qutip objects.
"""
import re
from functools import reduce

from numpy import arange
from numpy import cos as np_cos
from numpy import diag as np_diag
from numpy import sin as np_sin
from scipy.linalg import pinv, svd
from scipy.sparse import csr_matrix
from sympy import symbols
from sympy.utilities.lambdify import lambdify

from ..core.abstract_algebra import Operation
from ..core.exceptions import AlgebraError
from ..core.hilbert_space_algebra import TrivialSpace
from ..core.operator_algebra import (
    Adjoint,
    IdentityOperator,
    LocalOperator,
    LocalSigma,
    NullSpaceProjector,
    Operator,
    OperatorPlus,
    OperatorTimes,
    OperatorTrace,
    PseudoInverse,
    ScalarTimesOperator,
    ZeroOperator,
)
from ..core.scalar_algebra import ScalarValue, is_scalar
from ..core.state_algebra import (
    BasisKet,
    BraKet,
    CoherentStateKet,
    KetBra,
    KetPlus,
    OperatorTimesKet,
    ScalarTimesKet,
    State,
    TensorKet,
)
from ..core.super_operator_algebra import (
    IdentitySuperOperator,
    ScalarTimesSuperOperator,
    SPost,
    SPre,
    SuperOperator,
    SuperOperatorPlus,
    SuperOperatorTimes,
    SuperOperatorTimesOperator,
    ZeroSuperOperator,
)
from ..library.fock_operators import Create, Destroy, Displace, Phase, Squeeze
from ..library.spin_algebra import Jminus, Jplus, Jz


try:
    import qutip
except ImportError:
    # we want QAlgebra to be importable even if qutip is not installed. In this
    # case, and exception will occur when any of the qutip conversion routines
    # are called
    pass


DENSE_DIMENSION_LIMIT = 1000

__all__ = ['convert_to_qutip']


def convert_to_qutip(expr, full_space=None, mapping=None):
    """Convert a QAlgebra expression to a qutip object

    Args:
        expr: a QAlgebra expression
        full_space (HilbertSpace): The
            Hilbert space in which `expr` is defined. If not given,
            ``expr.space`` is used. The Hilbert space must have a well-defined
            basis.
        mapping (dict): A mapping of any (sub-)expression to either a
            `quip.Qobj` directly, or to a callable that will convert the
            expression into a `qutip.Qobj`. Useful for e.g. supplying objects
            for symbols
    Raises:
        ValueError: if `expr` is not in `full_space`, or if `expr` cannot be
            converted.
    """
    if full_space is None:
        full_space = expr.space
    if not expr.space.is_tensor_factor_of(full_space):
        raise ValueError(
            "expr '%s' must be in full_space %s" % (expr, full_space)
        )
    if full_space == TrivialSpace:
        raise AlgebraError(
            "Cannot convert object in TrivialSpace to qutip. "
            "You may pass a non-trivial `full_space`"
        )
    if mapping is not None:
        if expr in mapping:
            ret = mapping[expr]
            if isinstance(ret, qutip.Qobj):
                return ret
            else:
                assert callable(ret)
                return ret(expr)
    if expr is IdentityOperator:
        local_spaces = full_space.local_factors
        if len(local_spaces) == 0:
            raise ValueError(
                "full_space %s does not have local factors" % full_space
            )
        else:
            return qutip.tensor(
                *[qutip.qeye(s.dimension) for s in local_spaces]
            )
    elif expr is ZeroOperator:
        return qutip.tensor(
            *[
                qutip.Qobj(csr_matrix((s.dimension, s.dimension)))
                for s in full_space.local_factors
            ]
        )
    elif isinstance(expr, LocalOperator):
        return _convert_local_operator_to_qutip(expr, full_space, mapping)
    elif isinstance(expr, Operator) and isinstance(expr, Operation):
        return _convert_operator_operation_to_qutip(expr, full_space, mapping)
    elif isinstance(expr, OperatorTrace):
        raise NotImplementedError('Cannot convert OperatorTrace to ' 'qutip')
        # actually, this is perfectly doable in principle, but requires a bit
        # of work
    elif isinstance(expr, State):
        return _convert_ket_to_qutip(expr, full_space, mapping)
    elif isinstance(expr, SuperOperator):
        return _convert_superoperator_to_qutip(expr, full_space, mapping)
    elif isinstance(expr, Operation):
        # This is assumed to be an Operation on states, as we have handled all
        # other Operations above. Eventually, a StateOperation should be
        # defined as a common superclass for the Operations in the state
        # algebra
        return _convert_state_operation_to_qutip(expr, full_space, mapping)
    else:
        raise ValueError(
            "Cannot convert '%s' of type %s" % (str(expr), type(expr))
        )


def _convert_local_operator_to_qutip(expr, full_space, mapping):
    """Convert a LocalOperator instance to qutip."""
    n = full_space.dimension
    if full_space != expr.space:
        all_spaces = full_space.local_factors
        own_space_index = all_spaces.index(expr.space)
        return qutip.tensor(
            *(
                [qutip.qeye(s.dimension) for s in all_spaces[:own_space_index]]
                + [convert_to_qutip(expr, expr.space, mapping=mapping)]
                + [
                    qutip.qeye(s.dimension)
                    for s in all_spaces[own_space_index + 1 :]
                ]
            )
        )
    if isinstance(expr, Create):
        return qutip.create(n)
    elif isinstance(expr, Jz):
        return qutip.jmat((expr.space.dimension - 1) / 2.0, "z")
    elif isinstance(expr, Jplus):
        return qutip.jmat((expr.space.dimension - 1) / 2.0, "+")
    elif isinstance(expr, Jminus):
        return qutip.jmat((expr.space.dimension - 1) / 2.0, "-")
    elif isinstance(expr, Destroy):
        return qutip.destroy(n)
    elif isinstance(expr, Phase):
        arg = complex(expr.operands[1]) * arange(n)
        d = np_cos(arg) + 1j * np_sin(arg)
        return qutip.Qobj(np_diag(d))
    elif isinstance(expr, Displace):
        alpha = expr.operands[1]
        return qutip.displace(n, alpha)
    elif isinstance(expr, Squeeze):
        eta = expr.operands[1]
        return qutip.displace(n, eta)
    elif isinstance(expr, LocalSigma):
        j = expr.j
        k = expr.k
        if isinstance(j, str):
            j = expr.space.basis_labels.index(j)
        if isinstance(k, str):
            k = expr.space.basis_labels.index(k)
        ket = qutip.basis(n, j)
        bra = qutip.basis(n, k).dag()
        return ket * bra
    else:
        raise ValueError(
            "Cannot convert '%s' of type %s" % (str(expr), type(expr))
        )


def _convert_operator_operation_to_qutip(expr, full_space, mapping):
    if isinstance(expr, OperatorPlus):
        return sum(
            (
                convert_to_qutip(op, full_space, mapping=mapping)
                for op in expr.operands
            ),
            0,
        )
    elif isinstance(expr, OperatorTimes):
        # if any factor acts non-locally, we need to expand distributively.
        if any(len(op.space) > 1 for op in expr.operands):
            se = expr.expand()
            if se == expr:
                raise ValueError(
                    "Cannot represent as QuTiP object: {!s}".format(expr)
                )
            return convert_to_qutip(se, full_space, mapping=mapping)
        all_spaces = full_space.local_factors
        by_space = []
        ck = 0
        for ls in all_spaces:
            # group factors by associated local space
            ls_ops = [
                convert_to_qutip(o, o.space, mapping=mapping)
                for o in expr.operands
                if o.space == ls
            ]
            if len(ls_ops):
                # compute factor associated with local space
                by_space.append(reduce(lambda a, b: a * b, ls_ops, 1))
                ck += len(ls_ops)
            else:
                # if trivial action, take identity matrix
                by_space.append(qutip.qeye(ls.dimension))
        assert ck == len(expr.operands)
        # combine local factors in tensor product
        return qutip.tensor(*by_space)
    elif isinstance(expr, Adjoint):
        return convert_to_qutip(
            qutip.dag(expr.operands[0]), full_space, mapping=mapping
        )
    elif isinstance(expr, ScalarTimesOperator):
        try:
            coeff = complex(expr.coeff)
        except TypeError:
            raise TypeError(
                "Scalar coefficient '%s' is not numerical" % expr.coeff
            )
        return coeff * convert_to_qutip(
            expr.term, full_space=full_space, mapping=mapping
        )
    elif isinstance(expr, PseudoInverse):
        mo = convert_to_qutip(
            expr.operand, full_space=full_space, mapping=mapping
        )
        if full_space.dimension <= DENSE_DIMENSION_LIMIT:
            arr = mo.data.toarray()
            piarr = pinv(arr)
            pimo = qutip.Qobj(piarr)
            pimo.dims = mo.dims
            pimo.isherm = mo.isherm
            pimo.type = 'oper'
            return pimo
        raise NotImplementedError(
            "Only implemented for smaller state " "spaces"
        )
    elif isinstance(expr, NullSpaceProjector):
        mo = convert_to_qutip(
            expr.operand, full_space=full_space, mapping=mapping
        )
        if full_space.dimension <= DENSE_DIMENSION_LIMIT:
            arr = mo.data.toarray()
            # compute Singular Value Decomposition
            U, s, Vh = svd(arr)
            tol = 1e-8 * s[0]
            zero_svs = s < tol
            Vhzero = Vh[zero_svs, :]
            PKarr = Vhzero.conjugate().transpose().dot(Vhzero)
            PKmo = qutip.Qobj(PKarr)
            PKmo.dims = mo.dims
            PKmo.isherm = True
            PKmo.type = 'oper'
            return PKmo
        raise NotImplementedError(
            "Only implemented for smaller state " "spaces"
        )
    else:
        raise ValueError(
            "Cannot convert '%s' of type %s" % (str(expr), type(expr))
        )


def _convert_state_operation_to_qutip(expr, full_space, mapping):
    if full_space != expr.space:
        all_spaces = full_space.local_factors
        own_space_index = all_spaces.index(expr.space)
        return qutip.tensor(
            *(
                [qutip.qeye(s.dimension) for s in all_spaces[:own_space_index]]
                + convert_to_qutip(expr, expr.space, mapping=mapping)
                + [
                    qutip.qeye(s.dimension)
                    for s in all_spaces[own_space_index + 1 :]
                ]
            )
        )
    if isinstance(expr, BraKet):
        bq = convert_to_qutip(expr.bra, full_space, mapping=mapping)
        kq = convert_to_qutip(expr.ket, full_space, mapping=mapping)
        return bq * kq
    elif isinstance(expr, KetBra):
        bq = convert_to_qutip(expr.bra, full_space, mapping=mapping)
        kq = convert_to_qutip(expr.ket, full_space, mapping=mapping)
        return kq * bq
    else:
        raise ValueError(
            "Cannot convert '%s' of type %s" % (str(expr), type(expr))
        )


def _convert_ket_to_qutip(expr, full_space, mapping):
    n = full_space.dimension
    if full_space != expr.space:
        all_spaces = full_space.local_factors
        own_space_index = all_spaces.index(expr.space)
        factors = (
            [qutip.qeye(s.dimension) for s in all_spaces[:own_space_index]]
            + [convert_to_qutip(expr, expr.space, mapping=mapping)]
            + [
                qutip.qeye(s.dimension)
                for s in all_spaces[own_space_index + 1 :]
            ]
        )
        return qutip.tensor(*factors)
    elif isinstance(expr, KetPlus):
        return sum(
            (
                convert_to_qutip(op, full_space, mapping=mapping)
                for op in expr.operands
            ),
            0,
        )
    elif isinstance(expr, TensorKet):
        if any(len(op.space) > 1 for op in expr.operands):
            se = expr.expand()
            if se == expr:
                raise ValueError(
                    "Cannot represent as QuTiP " "object: {!s}".format(expr)
                )
            return convert_to_qutip(se, full_space, mapping=mapping)
        factors = [
            convert_to_qutip(o, o.space, mapping=mapping)
            for o in expr.operands
        ]
        return qutip.tensor(*factors)
    elif isinstance(expr, ScalarTimesKet):
        return complex(expr.coeff) * convert_to_qutip(
            expr.term, full_space, mapping=mapping
        )
    elif isinstance(expr, OperatorTimesKet):
        return convert_to_qutip(
            expr.coeff, full_space, mapping=mapping
        ) * convert_to_qutip(expr.term, full_space, mapping=mapping)
    elif isinstance(expr, BasisKet):
        return qutip.basis(expr.space.dimension, expr.index)
    else:
        raise ValueError(
            "Cannot convert '%s' of type %s" % (str(expr), type(expr))
        )


def _convert_superoperator_to_qutip(expr, full_space, mapping):
    if full_space != expr.space:
        all_spaces = full_space.local_factors
        own_space_index = all_spaces.index(expr.space)
        return qutip.tensor(
            *(
                [qutip.qeye(s.dimension) for s in all_spaces[:own_space_index]]
                + convert_to_qutip(expr, expr.space, mapping=mapping)
                + [
                    qutip.qeye(s.dimension)
                    for s in all_spaces[own_space_index + 1 :]
                ]
            )
        )
    if isinstance(expr, IdentitySuperOperator):
        return qutip.spre(
            qutip.tensor(
                *[qutip.qeye(s.dimension) for s in full_space.local_factors]
            )
        )
    elif isinstance(expr, SuperOperatorPlus):
        return sum(
            (
                convert_to_qutip(op, full_space, mapping=mapping)
                for op in expr.operands
            ),
            0,
        )
    elif isinstance(expr, SuperOperatorTimes):
        ops_qutip = [
            convert_to_qutip(o, full_space, mapping=mapping)
            for o in expr.operands
        ]
        return reduce(lambda a, b: a * b, ops_qutip, 1)
    elif isinstance(expr, ScalarTimesSuperOperator):
        return complex(expr.coeff) * convert_to_qutip(
            expr.term, full_space, mapping=mapping
        )
    elif isinstance(expr, SPre):
        return qutip.spre(
            convert_to_qutip(expr.operands[0], full_space, mapping)
        )
    elif isinstance(expr, SPost):
        return qutip.spost(
            convert_to_qutip(expr.operands[0], full_space, mapping)
        )
    elif isinstance(expr, SuperOperatorTimesOperator):
        sop, op = expr.operands
        return convert_to_qutip(
            sop, full_space, mapping=mapping
        ) * convert_to_qutip(op, full_space, mapping=mapping)
    elif isinstance(expr, ZeroSuperOperator):
        return qutip.spre(
            convert_to_qutip(ZeroOperator, full_space, mapping=mapping)
        )
    else:
        raise ValueError(
            "Cannot convert '%s' of type %s" % (str(expr), type(expr))
        )


def _time_dependent_to_qutip(
    op,
    full_space=None,
    time_symbol=symbols("t", real=True),
    convert_as='pyfunc',
):
    """Convert a possiblty time-dependent operator into the nested-list
    structure required by QuTiP"""
    if full_space is None:
        full_space = op.space
    if time_symbol in op.free_symbols:
        op = op.expand()
        if isinstance(op, OperatorPlus):
            result = []
            for o in op.operands:
                if time_symbol not in o.free_symbols:
                    if len(result) == 0:
                        result.append(
                            convert_to_qutip(o, full_space=full_space)
                        )
                    else:
                        result[0] += convert_to_qutip(o, full_space=full_space)
            for o in op.operands:
                if time_symbol in o.free_symbols:
                    result.append(
                        _time_dependent_to_qutip(
                            o, full_space, time_symbol, convert_as
                        )
                    )
            return result
        elif isinstance(op, ScalarTimesOperator) and isinstance(
            op.coeff, ScalarValue
        ):
            if convert_as == 'pyfunc':
                func_no_args = lambdify(time_symbol, op.coeff.val)
                if {time_symbol} == op.coeff.free_symbols:

                    def func(t, args):
                        # args are ignored for increased efficiency, since we
                        # know there are no free symbols except t
                        return func_no_args(t)

                else:

                    def func(t, args):
                        return func_no_args(t).subs(args)

                coeff = func
            elif convert_as == 'str':
                # a bit of a hack to replace imaginary unit
                # TODO: we can probably use one of the sympy code generation
                # routines, or lambdify with 'numexpr' to implement this in a
                # more robust way
                coeff = re.sub("I", "(1.0j)", str(op.coeff.val))
            else:
                raise ValueError(
                    (
                        "Invalid value '%s' for `convert_as`, must "
                        "be one of 'str', 'pyfunc'"
                    )
                    % convert_as
                )
            return [convert_to_qutip(op.term, full_space), coeff]
        else:
            raise ValueError(
                "op cannot be expressed in qutip. It must have "
                "the structure op = sum_i f_i(t) * op_i"
            )
    else:
        return convert_to_qutip(op, full_space=full_space)
