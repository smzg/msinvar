r"""
Category of Lambda Rings

EXAMPLE::

    sage: from msinvar.lambda_rings import LambdaRings
    sage: R=PolynomialRing(QQ, 'x')
    sage: LambdaRings.add_ring(R)
    sage: x=R.gen(); (1+x).adams(2)
    x^2 + 1
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.commutative_rings import CommutativeRings
from sage.misc.misc_c import prod


class LambdaRings(Category_singleton):
    r"""
    The category of lambda-rings -- commutative rings with plethystic
    operations.

    To add a parent to the category one needs to call:
    :meth:`LambdaRings.add_ring`.

    EXAMPLE::

        sage: from msinvar.lambda_rings import LambdaRings
        sage: R=PolynomialRing(QQ,2,'x,y')
        sage: LambdaRings.add_ring(R)
        sage: x,y=R.gens()
        sage: (x+y).adams(2)
        x^2 + y^2

    We can add an existing parent to lambda-rings, or we can use the init
    method of a parent. For example, :class:`msinvar.tm_polynomials.TMPoly`
    is automatically equipped with a lambda-ring structure.

    EXAMPLE::

        sage: from msinvar import TMPoly
        sage: R1=TMPoly(R,1,'z'); z=R1.gen()
        sage: (x*z).adams(2)
        x^2*z^2

    The default adams operation is :meth:`default_adams`.
    To override it one should define a new method :meth:`adams` in the parent
    or in the element class.

    For existing parent instances to override the default adams operation
    one can call::

        LambdaRings.add_ring(R, adams)

    where ``adams`` is the new adams operation.
    """

    dct_adams = {}

    @staticmethod
    def add_ring(R, adams=None):
        """Add ``R`` to the category of lambda-rings.

        In particular, equip ``R`` and its elements with the adams operation.
        """
        R._refine_category_(LambdaRings())
        if adams is not None:
            LambdaRings.dct_adams[R] = adams

    def super_categories(self):
        """Return the immediate super categories of ``self``."""
        return [CommutativeRings()]

    class ParentMethods:

        def is_lambda_ring(self):
            return False

        def adams(self, a, n):
            dct = LambdaRings.dct_adams
            if self in dct:
                return dct[self](a, n)
            return default_adams(a, n)

    class ElementMethods:

        def adams(self, n):
            r"""
            Adams operation `\psi_n` applied to ``self``.
            """
            return self.parent().adams(self, n)

        def plet(self, f):
            r"""
            Return plethysm f[self], where f is a symmetric function.

            Note that for symmetric functions the method plethysm(self, a)
            returns self[a]. For this reason we use a different name for our method.
            """
            p = f.parent().realization_of().powersum()

            def g(part):
                return prod(self.adams(n) for n in part)
            return p._apply_module_morphism(p(f), g)


def is_LambdaRingElement(a):
    return a.parent() in LambdaRings


def default_adams(f, n):
    """
    Return the default Adams operation.

    It raises all variables in ``f`` to the ``n``-th power.
    """
    try:
        d = {x: x**n for x in f.parent().gens()}
        return f.subs(d)
    except (AttributeError, TypeError, ValueError):
        return f


def adams(f, n):
    try:
        return f.adams(n)
    except (AttributeError, TypeError, ValueError):
        try:
            return f.parent().adams(f, n)
        except (AttributeError, TypeError, ValueError):
            return default_adams(f, n)
