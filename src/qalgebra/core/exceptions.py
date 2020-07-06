"""Exceptions and Errors raised by QAlgebra."""


__all__ = [
    'AlgebraException',
    'AlgebraError',
    'InfiniteSumError',
    'CannotSimplify',
    'BasisNotSetError',
    'UnequalSpaces',
    'OverlappingSpaces',
    'SpaceTooLargeError',
    'CannotSymbolicallyDiagonalize',
    'NonSquareMatrix',
    'NoConjugateMatrix',
]


class AlgebraException(Exception):
    """Base class for all algebraic exceptions.

    These should only be used for internat signaling, and never be raised to a
    user.
    """


class AlgebraError(ValueError):
    """Base class for all algebraic errors."""


class InfiniteSumError(AlgebraError):
    """Raised when expanding a sum into an infinite number of terms."""


class CannotSimplify(AlgebraException):
    """Raised when a rule cannot further simplify an expression."""


class BasisNotSetError(AlgebraError):
    """Raised if the basis or a Hilbert space dimension is unavailable."""


class UnequalSpaces(AlgebraError):
    """Raised when objects fail to be in the same Hilbert space.

    This happens for example when trying to add two states from different
    Hilbert spaces.
    """


class OverlappingSpaces(AlgebraError):
    """Raised when objects fail to be in separate Hilbert spaces."""


class SpaceTooLargeError(AlgebraError):
    """Raised when objects fail to be have overlapping Hilbert spaces."""


class CannotSymbolicallyDiagonalize(AlgebraException):
    """Matrix cannot be diagonalized analytically.

    Signals that a fallback to numerical diagonalization is required.
    """


class NonSquareMatrix(AlgebraError):
    """Raised when a :class:`.Matrix` fails to be square."""


class NoConjugateMatrix(AlgebraError):
    """Raised when entries of :class:`.Matrix` have no defined conjugate."""
