.. _library_structure:

=================
Library Structure
=================


Subpackage Organization
=======================


QAlgebra is organized into the sub-packages outlined below. Each
package may in turn contain several sub-modules.

Every package exports all public symbol from all of its sub-packages/-modules
in a "flat" API. Thus, a user can directly import from the top-level :mod:`qalgebra`
package.

In order from high-level to low-level:

.. autosummary::
    qalgebra
    qalgebra.convert
    qalgebra.printing
    qalgebra.toolbox
    qalgebra.library
    qalgebra.core
    qalgebra.pattern_matching
    qalgebra.utils

See also the full :ref:`modindex`


Class Hierarchy
===============

The following is an inheritance diagram of *all* the classes defined in QAlgebra
(this is best viewed as the full-page SVG):

.. inheritance-diagram:: qalgebra
   :parts: 1
   :top-classes: sympy.core.symbol.Symbol
   :cluster_modules:
