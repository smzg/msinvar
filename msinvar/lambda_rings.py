r"""
Category of Lambda Rings
::
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
    """
    The category of lambda-rings
    -- commutative rings with plethystic operations.

    To add a parent to the category one needs to call::

        LambdaRings.add_ring(R)

    This can be done with existing parent instances or in the init method
    of a parent.
    The default adams operation is **default_adams**.
    To override it one should define a new method **adams** in the parent
    or in the element class.
    For existing parent instances to override the default adams operation
    one needs to call::

        LambdaRings.add_ring(R, adams)

    where **adams** is a new adams operation. It is stored in the attribute
    _adams of a parent.

    s a dictionary inside
    the class **LambdaRings**.

    """

    @staticmethod
    def add_ring(R, adams=None):
        R._refine_category_(LambdaRings())
        R._adams = adams

    def super_categories(self):
        return [CommutativeRings()]

    class ParentMethods:

        _adams = None

        def is_lambda_ring(self):
            return False

        def adams(self, a, n):
            if self._adams is None:
                return default_adams(a, n)
            return self._adams(a, n)

    class ElementMethods:

        def adams(self, n):
            r"""
            Adams operation $\psi_n(self)$
            """
            return self.parent().adams(self, n)

        def plet(self, f):
            r"""
            Return plethysm f[self], where f is a symmetric function.

            Note that for symmetric functions the method plethysm(self, a)
            returns self[a]. For this reason we use a different name for our method.
            """
            p = f.parent().realization_of().powersum()
            def g(part): return prod(self.adams(n) for n in part)
            return p._apply_module_morphism(p(f), g)


def is_LambdaRingElement(a):
    return a.parent() in LambdaRings


def default_adams(f, n):
    try:
        d = {x: x**n for x in f.parent().gens()}
        return f.subs(d)
    except:
        return f


def adams(f, n):
    try:
        return f.adams(n)
    except:
        try:
            return f.parent().adams(f, n)
        except:
            return default_adams(f, n)
