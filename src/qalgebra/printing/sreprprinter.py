"""Provides printers for a full-structured representation"""

from textwrap import dedent, indent

from sympy.core.basic import Basic as SympyBasic

from ..core.abstract_algebra import Expression
from ..utils.singleton import Singleton
from ._render_head_repr import render_head_repr
from .base import QalgebraBasePrinter
from .sympy import SympyReprPrinter


__all__ = []
__private__ = [
    'QalgebraSReprPrinter',
    'IndentedSReprPrinter',
    'IndentedSympyReprPrinter',
]


class QalgebraSReprPrinter(QalgebraBasePrinter):
    """Printer for a string (ASCII) representation."""

    sympy_printer_cls = SympyReprPrinter

    def _render_str(self, string):
        return repr(string)

    def emptyPrinter(self, expr):
        """Fallback printer"""
        return render_head_repr(expr, sub_render=self.doprint)

    def _render_matrix(self, data, prefix='', postfix=''):
        if len(data.shape) == 2:
            rows = []
            for row in data:
                rows.append(
                    '[' + ", ".join([self.doprint(val) for val in row]) + ']'
                )
            return prefix + "[" + ", ".join(rows) + "]" + postfix
        else:
            raise ValueError("Cannot render %r" % data)

    def _print_Matrix(self, expr):
        return self._render_matrix(expr.matrix, prefix='Matrix(', postfix=')')

    def _print_ndarray(self, expr):
        return self._render_matrix(
            expr, prefix='array(', postfix=(', dtype=%s)' % expr.dtype)
        )


class IndentedSympyReprPrinter(SympyReprPrinter):
    """Indented repr printer for Sympy objects."""

    def doprint(self, expr):
        res = super().doprint(expr)
        return "    " * (self._print_level - 1) + res


class IndentedSReprPrinter(QalgebraBasePrinter):
    """Printer for rendering an expression in such a way that the resulting
    string can be evaluated in an appropriate context to re-instantiate an
    identical object, using nested indentation (implementing
    ``srepr(expr, indented=True)``
    """

    sympy_printer_cls = IndentedSympyReprPrinter

    def _get_from_cache(self, expr):
        """Obtain cached result, prepend with the keyname if necessary, and
        indent for the current level"""
        is_cached, res = super()._get_from_cache(expr)
        if is_cached:
            indent_str = "    " * self._print_level
            return True, indent(res, indent_str)
        else:
            return False, None

    def _write_to_cache(self, expr, res):
        """Store the cached result without indentation, and without the
        keyname"""
        res = dedent(res)
        super()._write_to_cache(expr, res)

    def _render_str(self, string):
        return "    " * self._print_level + repr(string)

    def emptyPrinter(self, expr):
        """Fallback printer"""
        indent_str = "    " * (self._print_level - 1)
        lines = []
        if isinstance(expr.__class__, Singleton):
            # We exploit that Singletons override __expr__ to directly return
            # their name
            return indent_str + repr(expr)
        if isinstance(expr, Expression):
            args = expr.args
            keys = expr.minimal_kwargs.keys()
            lines.append(indent_str + expr.__class__.__name__ + "(")
            for arg in args:
                lines.append(self.doprint(arg) + ",")
            for key in keys:
                arg = expr.kwargs[key]
                lines.append(
                    ("    " * self._print_level)
                    + key
                    + '='
                    + self.doprint(arg).lstrip()
                    + ","
                )
            if len(args) > 0 or len(keys) > 0:
                lines[-1] = lines[-1][:-1]  # drop trailing comma for last arg
            lines[-1] += ")"
        elif isinstance(expr, (tuple, list)):
            delims = ("(", ")") if isinstance(expr, tuple) else ("[", "]")
            if len(expr) == 1:
                delims = (delims[0], "," + delims[1])
            lines.append(
                indent_str
                + delims[0]
                + ", ".join([render_head_repr(v) for v in expr])
                + delims[1]
            )
        else:
            lines.append(indent_str + SympyReprPrinter().doprint(expr))
        return "\n".join(lines)

    def _render_matrix(self, data, prefix='', postfix=''):
        indent_str = "    " * (self._print_level - 1)
        if len(data.shape) == 2:
            lines = [
                indent_str + prefix + "[",
            ]
            self._print_level += 1
            for row in data:
                indent_str = "    " * (self._print_level - 1)
                lines.append(indent_str + '[')
                for val in row:
                    lines.append(self.doprint(val) + ",")
                lines[-1] = lines[-1][:-1]
                lines.append(indent_str + '],')
            lines[-1] = lines[-1][:-1] + "]" + postfix
            return "\n".join(lines)
        else:
            raise ValueError("Cannot render %r" % data)

    def _print_Matrix(self, expr):
        return self._render_matrix(expr.matrix, prefix='Matrix(', postfix=')')

    def _print_ndarray(self, expr):
        return self._render_matrix(
            expr, prefix='array(', postfix=(', dtype=%s)' % expr.dtype)
        )
