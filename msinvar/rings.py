r"""
Some base rings in which our invariants take values

EXAMPLES::

    sage: from msinvar.rings import RF
    sage: R=RF('u,v')
    sage: R.inject_variables(verbose=False)
    sage: (u-v).adams(2)/(u+v)
    u - v
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.fraction_field import FractionField_generic
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
from sage.rings.rational_field import QQ
from sage.categories.quotient_fields import QuotientFields
from msinvar.lambda_rings import LambdaRings


class RationalFunctionField(FractionField_generic):
    """
    Field of rational functions in several variables, meaning fractions P/Q, where P, Q
    are polynomials in several variables.

    EXAMPLES::

        sage: from msinvar.rings import RF
        sage: R=RF(); R
        Field of Rational Functions in y
        sage: y=R.gen()
        sage: f=(1-y)**3/(1+y)**2; f
        (-y^3 + 3*y^2 - 3*y + 1)/(y^2 + 2*y + 1)
        sage: factor(f)
        (-1) * (y + 1)^-2 * (y - 1)^3
        sage: f.adams(2)
        (-y^6 + 3*y^4 - 3*y^2 + 1)/(y^4 + 2*y^2 + 1)

        sage: R=RF('x,y'); R
        Field of Rational Functions in x, y
        sage: R.inject_variables(verbose=False)
        sage: f=(1-y**2)**5/(1-y)**5*(x-y)
        sage: f.factor()
        (-1) * (-x + y) * (y + 1)^5
    """

    def __init__(self, vars='y', base=QQ):
        vars = vars.split(',')
        R = MPolynomialRing_libsingular(base, n=len(vars), names=vars)
        # cat=Category.join([QuotientFields(),LambdaRings()])
        super().__init__(R, element_class=RationalFunction)  # , category=cat)
        LambdaRings.add_ring(self)

    def _repr_(self):
        vars = ', '.join(self.variable_names())
        return 'Field of Rational Functions in ' + vars


RF = RationalFunctionField


class RationalFunction(FractionFieldElement):
    """Element class for the parent class :class:`RationalFunctionField`."""

    def root_vars(self, k=2):
        """See :meth:`root_vars`."""
        return root_vars(self, k)

    def simp(f):
        if f == 0:
            return f
        num, den = f.numerator(), f.denominator()
        content = den.content()
        num /= content
        den /= content
        return num / den


def root_vars(f, k=2):
    """
    Substitute every variable x in a fraction f by x^(1/k).
    """
    R = f.parent()
    if R in QuotientFields:
        return R(root_vars(f.numerator(), k) / root_vars(f.denominator(), k))

    def root(e):
        return tuple(i // k for i in e)
    dct = {root(e): c for e, c in f.dict().items()}
    return R(dct)
