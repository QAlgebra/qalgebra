"""Matrices of Expressions."""
import numpy as np
import sympy
from sympy import I, Symbol, sympify

from .abstract_algebra import Expression, substitute
from .abstract_quantum_algebra import QuantumExpression
from .exceptions import NoConjugateMatrix, NonSquareMatrix
from .hilbert_space_algebra import ProductSpace, TrivialSpace
from .operator_algebra import adjoint
from .scalar_algebra import is_scalar


__all__ = [
    'Matrix',
    'block_matrix',
    'diagm',
    'hstackm',
    'identity_matrix',
    'vstackm',
    'zerosm',
]

__private__ = []  # anything not in __all__ must be in __private__


class Matrix:
    """Matrix of Expressions."""

    matrix = None
    _hash = None

    def __init__(self, m):
        if isinstance(m, np.ndarray):
            self.matrix = m
        elif isinstance(m, Matrix):
            self.matrix = np.array(m.matrix)
        else:
            self.matrix = np.array(m)
        if len(self.matrix.shape) < 2:
            self.matrix = self.matrix.reshape((self.matrix.shape[0], 1))
        if len(self.matrix.shape) > 2:
            raise ValueError("Must have a shape of length 2")

    @property
    def shape(self):
        """The shape of the matrix ``(nrows, ncols)``."""
        return self.matrix.shape

    @property
    def block_structure(self):
        """For square matrices this gives the block (-diagonal) structure of
        the matrix as a tuple of integers that sum up to the full dimension.

        :rtype: tuple
        """
        n, m = self.shape
        if n != m:
            raise AttributeError(
                "block_structure only defined for square matrices"
            )
        for k in range(1, n):
            if (self.matrix[:k, k:] == 0).all() and (
                self.matrix[k:, :k] == 0
            ).all():
                return (k,) + self[k:, k:].block_structure
        return (n,)

    def _get_blocks(self, block_structure):
        n, m = self.shape
        if n == m:
            if not sum(block_structure) == n:
                raise ValueError()
            if not len(block_structure):
                return ()
            j = block_structure[0]

            if (self.matrix[:j, j:] == 0).all() and (
                self.matrix[j:, :j] == 0
            ).all():
                return (self[:j, :j],) + self[j:, j:]._get_blocks(
                    block_structure[1:]
                )
            else:
                raise ValueError()
        elif m == 1:
            if not len(block_structure):
                return ()
            else:
                return (self[: block_structure[0], :],) + self[
                    : block_structure[0], :
                ]._get_blocks(block_structure[1:])
        else:
            raise ValueError()

    @property
    def is_zero(self):
        """Are all elements of the matrix zero?"""
        for o in self.matrix.ravel():
            try:
                if not o.is_zero:
                    return False
            except AttributeError:
                if not o == 0:
                    return False
        return True

    def __hash__(self):
        if not self._hash:
            self._hash = hash(
                (tuple(self.matrix.ravel()), self.matrix.shape, Matrix)
            )
        return self._hash

    def __eq__(self, other):
        if isinstance(other, Matrix):
            return np.all(self.matrix == other.matrix)
        else:
            return np.all(self.matrix == other)

    def __add__(self, other):
        if isinstance(other, Matrix):
            return Matrix(self.matrix + other.matrix)
        else:
            return Matrix(self.matrix + other)

    def __radd__(self, other):
        return Matrix(other + self.matrix)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            return Matrix(self.matrix.dot(other.matrix))
        else:
            return Matrix(self.matrix * other)

    def __rmul__(self, other):
        return Matrix(other * self.matrix)

    def __sub__(self, other):
        return self + (-1) * other

    def __rsub__(self, other):
        return (-1) * self + other

    def __neg__(self):
        return (-1) * self

    def __truediv__(self, other):
        if is_scalar(other):
            return self * (sympify(1) / other)
        raise NotImplementedError(
            "Can't divide matrix %s by %s" % (self, other)
        )

    def transpose(self):
        """The transpose matrix"""
        return Matrix(self.matrix.T)

    def conjugate(self):
        """The element-wise conjugate matrix.

        This is defined only if all the entries in the matrix have a defined
        conjugate (i.e., they have a `conjugate` method). This is *not* the
        case for a matrix of operators. In such a case, only an
        elementwise :func:`adjoint` would be applicable, but this is
        mathematically different from a complex conjugate.

        Raises:
            NoConjugateMatrix: if any entries have no `conjugate` method
        """
        try:
            return Matrix(np.conjugate(self.matrix))
        except (AttributeError, TypeError):
            raise NoConjugateMatrix(
                "Matrix %s contains entries that have no defined "
                "conjugate" % str(self)
            )

    @property
    def real(self):
        """Element-wise real part.

        Raises:
            NoConjugateMatrix: if entries have no `conjugate` method and no
                other way to determine the real part

        Note:
            A mathematically equivalent way to obtain a real matrix from a
            complex matrix ``M`` is::

                (M.conjugate() + M) / 2

            However, the result may not be identical to ``M.real``, as the
            latter tries to convert elements of the matrix to real values
            directly, if possible, and only uses the conjugate as a fall-back
        """

        def re(val):
            if hasattr(val, 'real'):
                return val.real
            elif hasattr(val, 'as_real_imag'):
                return val.as_real_imag()[0]
            elif hasattr(val, 'conjugate'):
                return (val.conjugate() + val) / 2
            else:
                raise NoConjugateMatrix(
                    "Matrix entry %s contains has no defined "
                    "conjugate" % str(val)
                )

        # Note: Do NOT use self.matrix.real! This will give wrong results, as
        # numpy thinks of objects (Operators) as real, even if they have no
        # defined real part
        return self.element_wise(re)

    @property
    def imag(self):
        """Element-wise imaginary part.

        Raises:
            NoConjugateMatrix: if entries have no `conjugate` method and no
                other way to determine the imaginary part

        Note:
            A mathematically equivalent way to obtain an imaginary matrix from
            a complex matrix ``M`` is::

                (M.conjugate() - M) / (I * 2)

            with same same caveats as :attr:`real`.
        """

        def im(val):
            if hasattr(val, 'imag'):
                return val.imag
            elif hasattr(val, 'as_real_imag'):
                return val.as_real_imag()[1]
            elif hasattr(val, 'conjugate'):
                return (val.conjugate() - val) / (2 * I)
            else:
                raise NoConjugateMatrix(
                    "Matrix entry %s contains has no defined "
                    "conjugate" % str(val)
                )

        # Note: Do NOT use self.matrix.real! This will give wrong results, as
        # numpy thinks of objects (Operators) as real, even if they have no
        # defined real part
        return self.element_wise(im)

    @property
    def T(self):
        """Alias for :meth:`transpose`."""
        return self.transpose()

    def adjoint(self):
        """Adjoint of the matrix.

        This is the transpose and the Hermitian adjoint of all elements."""
        return self.T.element_wise(adjoint)

    dag = adjoint

    def trace(self):
        if self.shape[0] == self.shape[1]:
            return sum(self.matrix[k, k] for k in range(self.shape[0]))
        raise NonSquareMatrix(repr(self))

    @property
    def H(self):
        """Alias for :meth:`adjoint`."""
        return self.adjoint()

    def __getitem__(self, item_id):
        item = self.matrix.__getitem__(item_id)
        if isinstance(item, np.ndarray):
            return Matrix(item)
        return item

    def element_wise(self, func, *args, **kwargs):
        """Apply a function to each matrix element and return the result in a
        new operator matrix of the same shape.

        Args:
            func (callable): A function to be applied to each element. It
                must take the element as its first argument.
            args: Additional positional arguments to be passed to `func`
            kwargs: Additional keyword arguments to be passed to `func`

        Returns:
            Matrix: Matrix with results of `func`, applied element-wise.
        """
        s = self.shape
        emat = [func(o, *args, **kwargs) for o in self.matrix.ravel()]
        return Matrix(np.array(emat).reshape(s))

    def series_expand(self, param: Symbol, about, order: int):
        """Expand the matrix expression as a truncated power series in a scalar
        parameter.

        Args:
            param: Expansion parameter.
            about (.Scalar): Point about which to expand.
            order: Maximum order of expansion >= 0

        Returns:
            tuple of length (order+1), where the entries are the expansion
            coefficients.
        """
        s = self.shape
        emats = zip(
            *[
                o.series_expand(param, about, order)
                for o in self.matrix.ravel()
            ]
        )
        return tuple((Matrix(np.array(em).reshape(s)) for em in emats))

    def expand(self):
        """Expand each matrix element distributively.

        Returns:
            Matrix: Expanded matrix.
        """
        return self.element_wise(
            lambda o: o.expand() if isinstance(o, QuantumExpression) else o
        )

    def substitute(self, var_map):
        """Perform a substitution in all element of the matrix.

        Equivalent to applying :func:`.substitute` element-wise.

        Returns:
            Matrix: Matrix with substitutions
        """
        if self in var_map:
            return var_map[self]
        else:
            return self.element_wise(substitute, var_map=var_map)

    @property
    def free_symbols(self):
        """Free symbols, across all elements."""
        ret = set()
        for o in self.matrix.ravel():
            try:
                ret = ret | o.free_symbols
            except AttributeError:
                pass
        return ret

    @property
    def space(self):
        """Combined Hilbert space of all matrix elements.

        If none of the elements have an associated hilbert space,
        :obj:`.TrivialSpace`.
        """
        arg_spaces = [
            o.space for o in self.matrix.ravel() if hasattr(o, 'space')
        ]
        if len(arg_spaces) == 0:
            return TrivialSpace
        else:
            return ProductSpace.create(*arg_spaces)

    def simplify_scalar(self, func=sympy.simplify):
        """Simplify all scalar expressions appearing in the Matrix."""

        def element_simplify(v):
            if isinstance(v, sympy.Basic):
                return func(v)
            elif isinstance(v, QuantumExpression):
                return v.simplify_scalar(func=func)
            else:
                return v

        return self.element_wise(element_simplify)

    def _repr_latex_(self):
        from qalgebra import latex

        return "$" + latex(self) + "$"


def hstackm(matrices):
    """Generalizes `numpy.hstack` to :class:`.Matrix` objects."""
    return Matrix(np.hstack(tuple(m.matrix for m in matrices)))


def vstackm(matrices):
    """Generalizes `numpy.vstack` to :class:`.Matrix` objects."""
    arr = np.vstack(tuple(m.matrix for m in matrices))
    #    print(tuple(m.matrix.dtype for m in matrices))
    #    print(arr.dtype)
    return Matrix(arr)


def diagm(v, k=0):
    """Generalizes the diagonal matrix creation capabilities of `numpy.diag` to
    :class:`.Matrix` objects."""
    return Matrix(np.diag(v, k))


def block_matrix(A, B, C, D):
    r"""Generate the operator matrix with quadrants

    .. math::

       \begin{pmatrix} A B \\ C D \end{pmatrix}

    Args:
        A (Matrix): Matrix of shape ``(n, m)``
        B (Matrix): Matrix of shape ``(n, k)``
        C (Matrix): Matrix of shape ``(l, m)``
        D (Matrix): Matrix of shape ``(l, k)``

    Returns:
        Matrix: The combined block matrix ``[[A, B], [C, D]]``.
    """
    return vstackm((hstackm((A, B)), hstackm((C, D))))


def identity_matrix(N):
    """Generate the N-dimensional identity matrix.

    Args:
        N (int): Dimension

    Returns:
        Matrix: Identity matrix in N dimensions

    """
    return diagm(np.ones(N, dtype=int))


def zerosm(shape, *args, **kwargs):
    """Generalizes ``numpy.zeros`` to :class:`.Matrix` objects."""
    return Matrix(np.zeros(shape, *args, **kwargs))
