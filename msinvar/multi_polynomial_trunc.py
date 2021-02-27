r"""
Truncated Multivariate polynomials

EXAMPLES::
    sage: from msinvar.multi_polynomial_trunc import MPTRing
    sage: R=MPTRing(QQ,2,'x',prec=(2,2))
    sage: R.inject_variables(verbose=False)
    sage: (x0+x1).Exp()
    1 + x0 + x1 + x0^2 + x0*x1 + x1^2 + x0^2*x1 + x0*x1^2 + x0^2*x1^2
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************

from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict
from sage.rings.infinity import infinity
from sage.arith.misc import moebius
from sage.rings.integer_ring import ZZ
from msinvar.lambda_rings import LambdaRings


def adams(f, n):
    try:
        return f.adams(n)
    except:
        pass
    try:
        return f.parent().adams(f, n)
    except:
        pass
    try:
        gens = tuple(x**n for x in f.parent().gens())
        return f(gens)
    except:
        return f


class MPolynomial_trunc(MPolynomial_polydict):
    """
    The element class for  truncated multivariate polynomials.
    The parent class is **MPolynomialRing_trunc**.

    **Implementation details:**

    1. The parent should have a method **le_prec(v)** to test if degree v should be included.
    2. Truncation takes place in the **__init__** method.
    3. Arithmetic operations are inherited from **MPolynomial_polydict**.
    4. We implement Adams operation **adams**, exp, log and plethystic Exp, Log.
    5. Function for the twisted multilpication is given by self._parent._twist_prod
    """

    def __init__(self, parent, x):
        self._parent = parent
        if isinstance(x, dict):
            d = x
        else:
            super().__init__(parent, x)
            d = self.dict()
        d = {e: c for e, c in d.items() if parent.le_prec(e)}
        super().__init__(parent, d)

    def _new_element(self, d):
        """Construct new poly from a dictionary d"""
        return self.__class__(self.parent(), d)

    def _new_constant_poly(self, x, P):
        # overriden; needed to make 1-x[0] of type MPT
        return self._new_element({self.parent().zero_tuple(): x})

    def le_prec(self, v):  # <= precision
        return self.parent().le_prec(v)

    def prec_num(self):
        return self.parent().prec_num()

    def exp(self):
        N = self.prec_num()
        if N == infinity:
            raise TypeError("Error: series is untruncated")
        s = self.parent().one()
        f, x = 1, self
        for i in range(1, N+1):
            f = x*f/i
            if f == 0:
                break
            s = s+f
        return s

    def log(self):
        N = self.prec_num()
        if N == infinity:
            raise TypeError("Error: series is untruncated")
        s = self.parent().zero()
        f, x = 1, 1-self  # need _new_constasnt_poly to have 1-self of type MPT!
        for i in range(1, N+1):
            f = x*f
            if f == 0:
                break
            s = s-f/i
        return s

    def invert(self):
        """
        For simplicity we invert only poly with constant coefficient 1
        """
        if self.is_one():
            return self
        N = self.prec_num()
        if N == infinity:
            raise TypeError("Error: series is untruncated")
        u = self.constant_coefficient()
        if not u.is_one():
            raise TypeError("The constant coefficient should be 1")
        s = self.parent().one()
        f, x = 1, 1-self
        for i in range(1, N+1):
            f = x*f
            if f == 0:
                break
            s = s + f
        return s

    def _div_(self, right):
        return self*right.invert()

    def Psi(self):
        N = self.prec_num()
        if N == infinity:
            raise TypeError("Error: series is untruncated")
        s = self.parent().zero()
        for i in range(1, N+2):
            f = self.adams(i)
            if f == 0:
                break
            s = s+f/i
        return s

    def IPsi(self):
        N = self.prec_num()
        if N == infinity:
            raise TypeError("Error: series is untruncated")
        s = self.parent().zero()
        for i in range(1, N+2):
            f = self.adams(i)
            if f == 0:
                break
            s = s+f*moebius(i)/i
        return s

    def Exp(self):
        return self.Psi().exp()

    def Log(self):
        return self.log().IPsi()

    def twist(self, f):
        """Multiply every term c*x^v by f(v)"""
        d = {e: c*f(e) for e, c in self.dict().items()}
        return self._new_element(d)

    def _mul_(self, right):
        f = self.parent().twist_prod
        if f is None:
            return super()._mul_(right)
        if len(self.dict()) == 0:   # product is zero anyways
            return self
        if right in self.parent().base_ring():
            return self*self._new_constant_poly(right, self.parent())
        d = {}
        for e0, c0 in self.dict().items():
            for e1, c1 in right.dict().items():
                e = e0.eadd(e1)
                if not self.le_prec(e):
                    continue
                c = c0*c1*f(e0, e1)
                if e in d:
                    d[e] = d[e]+c
                else:
                    d[e] = c
        return self._new_element(d)

    def coeff(self, e):
        # returns element from the base_ring; note that coefficient returns element from the ring
        try:
            return self.element()[list(e)]
        except:
            return 0


class MPolynomialRing_trunc(MPolynomialRing_polydict):
    """
    Truncated Multivariable polynomial ring.

    Multivariate polynomials up to degree prec,
    where prec is a degree vector or an integer (total degree).
    Monomials of degree prec are included.
    Adams operations, plethystic Exp and Log are implemented
    (without connection to the category of lambda rings in sage).

    This is the parent class for truncated multivariate polynomials.
    It has an alias **MPTRing**.
    The element class is **MPolynomial_trunc**.

    PARAMETERS:

        1. prec - precision vector (or integer, or infinity) for truncation
        2. twist_prod - function that maps exponents (d,e) to the base_ring. The new product is x^d*x^e=twist_prod(d,e)*x^(d+e)

    EXAMPLES::

        sage: from msinvar.multi_polynomial_trunc import MPTRing
        sage: R=MPTRing(QQ,2,prec=(2,2)); R
        Multivariate Polynomial Ring in x0, x1 over Rational Field truncated at degree (2, 2)
        sage: x=R.gens(); (x[0]+x[1])**3
        3*x0^2*x1 + 3*x0*x1^2

        sage: QR=Frac(PolynomialRing(QQ,'y,t'))
        sage: S=MPTRing(QR,2,'x',prec=(2,2))
        sage: y,t=QR.gens(); x=S.gens()
        sage: (y*x[0]+x[1]).adams(2)
        y^2*x0^2 + x1^2
        sage: (y*x[0]).Exp()
        1 + y*x0 + y^2*x0^2

        sage: R=MPTRing(QQ,1,prec=2, twist_prod=lambda a,b:2)
        sage: x=R.gen(); x*x
        2*x^2
    """

    Element = MPolynomial_trunc

    def __init__(self, base_ring, n, names='x', order='negdegrevlex',
                 prec=infinity, twist_prod=None):
        """prec is a precision vector for the truncation.
        twist_prod is a function that maps exponents (u,v) to the base_ring"""
        # Parent.__init__(self, category=LambdaRings())#not needed anymore
        self._prec = prec
        self.twist_prod = twist_prod
        super().__init__(base_ring, n, names, order)
        LambdaRings.add_ring(self)

    def adams(self, a, n):
        d = {e.emul(n): adams(c, n) for e, c in a.dict().items()
             if self.le_prec(e.emul(n))}
        return self.Element(self, d)

    def __call__(self, x, check=True):  # needed so that the result is MPT
        x = super().__call__(x, check)
        return self.Element(self, x.dict())

    def _poly_class(self):  # needed so that gens() are of type MPT
        return self.Element

    def _repr_(self):
        s = 'total degree ' if isinstance(self._prec, int) else 'degree '
        return super()._repr_()+' truncated at '+s+str(self._prec)

    def zero_tuple(self):
        return self._zero_tuple

    def prec_num(self):
        if self._prec == infinity:
            return infinity
        if self._prec in ZZ:
            return self._prec
        return sum(self._prec)

    def le_prec(self, v):
        d = self._prec
        if d == infinity:
            return True
        if d in ZZ:
            return sum(v) <= d
        return all(i <= j for i, j in zip(v, d))

    def set_twist_product(self, twist_prod):
        self.twist_prod = twist_prod

    def inject(self):
        self.inject_variables(verbose=False)


MPTRing = MPolynomialRing_trunc
