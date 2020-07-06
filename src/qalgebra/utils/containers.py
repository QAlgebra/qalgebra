"""Tools for working with data structures built from native containers.

"""

from collections import OrderedDict
from collections.abc import Container, Iterable, Mapping, Sequence, Sized


__all__ = []

__private__ = [
    "sorted_if_possible",
    "nested_tuple",
]  # anything not in __all__ must be in __private__


def sorted_if_possible(iterable, **kwargs):
    """Create a sorted list of elements of an iterable if they are orderable.

    See `sorted` for details on optional arguments to customize the sorting.

    Args:
        Iterable of a finite number of elements to sort.
    kwargs:
        Keyword arguments are passed on to `sorted`.

    Returns:
        List of elements in `iterable`, sorted if orderable, otherwise kept in
        the order of iteration.

    """
    try:
        return sorted(iterable, **kwargs)
    except TypeError:
        return list(iterable)


def nested_tuple(container):
    """Recursively transform a container structure to a nested tuple.

    The function understands container types inheriting from the selected
    abstract base classes in `collections.abc`, and performs the following
    replacements:

    * `Mapping` to tuple of key-value pair tuples.
    * `Sequence` to tuple containing the same elements in unchanged order.
    * `Collection` to tuple containing the same elements in sorted order if
        orderable and otherwise in the order of iteration.

    The function recurses into these container types to perform the same
    replacement, and leaves objects of other types untouched.

    The returned container is hashable if and only if all the values contained
    in the original data structure are hashable.

    Args:
        container: Data structure to transform into a nested tuple.

    Returns:
        Nested tuple containing the same data as `container`.

    """
    if isinstance(container, OrderedDict):
        return tuple(map(nested_tuple, container.items()))
    if isinstance(container, Mapping):
        return tuple(sorted_if_possible(map(nested_tuple, container.items())))
    if not isinstance(container, (str, bytes)):
        if isinstance(container, Sequence):
            return tuple(map(nested_tuple, container))
        if (
            isinstance(container, Container)
            and isinstance(container, Iterable)
            and isinstance(container, Sized)
        ):
            return tuple(sorted_if_possible(map(nested_tuple, container)))
    return container
