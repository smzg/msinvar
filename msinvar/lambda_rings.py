r"""
Category of Lambda Rings
::
    sage: from msinvar.lambda_rings import LambdaRings
    sage: R=PolynomialRing(QQ, 'x')
    sage: LambdaRings.add_ring(R)
    sage: x=R.gen(); (1+x).adams(2)
    x^2 + 1
"""

#*****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

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
    The default adams operation is **_default_adams**.
    To override it one should define a new method **adams** in the parent
    or in the element class.
    For existing parent instances to override the default adams operation
    one needs to call::
        LambdaRings.add_ring(R, adams)

    where **adams** is a new adams operation. It is stored is a dictionary inside
    the class **LambdaRings**.

    """
    _dct_adams = {}

    @staticmethod
    def add_ring(R, adams=None):
        R._refine_category_(LambdaRings())
        R.set_adams_operation(adams)

    def _default_adams(a, n):
        try:
            p = a.parent()
        except:
            return a
        try:
            d = {x: x**n for x in p.gens()}
            return a.subs(d)
        except:
            return a

    def super_categories(self):
        return [CommutativeRings()]

    class ParentMethods:

        def is_lambda_ring(self):
            return True

        def set_adams_operation(self, adams=None):
            """
            Define the adams operation for the parent to be equal adams.
            Insert this parent to the common dictionary of the category LambdaRings.
            """
            # set Adams operation psi_n(a) to be a function f=f(a,n)
            if adams is not None:
                LambdaRings._dct_adams[self] = adams

        def get_adams_operation(self):
            # return Adams operation as a function
            dct = LambdaRings._dct_adams
            if self in dct:
                return dct[self]
            return LambdaRings._default_adams

        def adams(self, a, n):
            return self.get_adams_operation()(a, n)

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
