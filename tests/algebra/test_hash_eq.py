"""Test hash and equality implementation of Expressions."""

from qalgebra.core.hilbert_space_algebra import LocalSpace
from qalgebra.library.fock_operators import Destroy


def test_equal_hash():
    """Test that expressions with and equal hash or not equal, and that they
    can be used as dictionary keys"""
    a = Destroy(hs="0")
    expr1 = -a
    expr2 = -2 * a
    h1 = hash(expr1)
    h2 = hash(expr2)
    # expr1 and expr2 just happen to have a hash collision for the current
    # implementation of the the hash function. This does not mean that they
    # are not distinguishable!
    assert h1 == h2  # this is the point of the test
    assert expr1 != expr2
    d = {}
    d[expr1] = 1
    d[expr2] = 2
    assert d[expr1] == 1
    assert d[expr2] == 2


def test_custom_localspace_identifier_hash():
    """Test hashes for expressions with different local_identifiers for their
    Hilbert spaces have different hashes"""
    hs1 = LocalSpace(1)
    hs1_custom = LocalSpace(1, local_identifiers={'Destroy': 'b'})
    assert hash(hs1) != hash(hs1_custom)
    a = Destroy(hs=hs1)
    b = Destroy(hs=hs1_custom)
    assert hash(a) != hash(b)
