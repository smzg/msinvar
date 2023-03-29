r"""
Truncated Multivariate polynomials

EXAMPLES::

    sage: from msinvar import TMPoly
    sage: R=TMPoly(QQ,2,'x',prec=(2,2))
    sage: R.inject_variables(verbose=False)
    sage: (x0+x1).Exp()
    1 + x0 + x1 + x0^2 + x0*x1 + x1^2 + x0^2*x1 + x0*x1^2 + x0^2*x1^2
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict
from sage.arith.misc import moebius
from sage.rings.rational_field import QQ
from msinvar.lambda_rings import LambdaRings, adams
from msinvar.utils import isiterable, vec


class TMPolynomial(MPolynomial_polydict):
    """
    The element class for  truncated multivariate polynomials.
    The parent class is :class:`TMPolynomialRing`.

    Implementation details:

    1. The parent should have a method :meth:`le_prec` to test if a given degree should be present.
    2. Truncation is performed in the :meth:`__init__` method.
    3. Arithmetic operations are inherited from :class:`MPolynomial_polydict`.
    4. We implement Adams operation adams, exp, log and plethystic Exp, Log.
    5. Function for the twisted multilpication is given by self.parent().prod_twist
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
        if N is None:
            raise ValueError("Series is untruncated")
        if self.constant_coefficient() != 0:
            raise ValueError("The constant coefficient should be 0")
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
        if N is None:
            raise ValueError("Series is untruncated")
        if self.constant_coefficient() != 1:
            raise ValueError("The constant coefficient should be 1")
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
        N = self.prec_num()
        if N is None:
            raise ValueError("Series is untruncated")
        if self.constant_coefficient() != 1:
            raise ValueError("The constant coefficient should be 1")
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
        if N is None:
            raise ValueError("Error: series is untruncated")
        if self.constant_coefficient() != 0:
            raise ValueError("The constant coefficient should be 0")
        s = self.parent().zero()
        for i in range(1, N+2):
            f = self.adams(i)
            if f == 0:
                break
            s = s+f/i
        return s

    def IPsi(self):
        N = self.prec_num()
        if N is None:
            raise ValueError("Error: series is untruncated")
        if self.constant_coefficient() != 0:
            raise ValueError("The constant coefficient should be 0")
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

    def term_twist(self, f):
        """Multiply every term c*x^v by f(v)"""
        d = {e: c*f(e) for e, c in self.dict().items()}
        return self._new_element(d)

    def _mul_(self, right):
        f = self.parent().prod_twist
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
        """
        Return the coefficient of degree e monomial as an element of
        the base_ring.

        Note that **coefficient** returns an element of the polynomial ring.
        """
        try:
            return self.element()[list(e)]
        except:
            return 0

    def subs_base(self, **kw):
        dct = {}
        for e, c in self.dict().items():
            try:
                c = c.subs(**kw)
            except:
                pass
            dct[e] = c
        return self._new_element(dct)

    def subs(self, **kw):
        f = super().subs(**kw)
        if f != self:
            return f
        return self.subs_base(**kw)

    def restrict(self, z, slope):
        """Restrict polynomial to dimension vectors d such that z(d)=slope,
        where ``z`` is a Stability (or the corresponding vector)."""
        from msinvar.stability import Stability
        z = Stability.check(z)

        dct = {}
        for e, c in self.dict().items():
            if vec.iszero(e) or z.has_slope(e, slope):
                dct[e] = c
        return self._new_element(dct)

    def invar(self):
        from msinvar.invariants import Invariant
        return Invariant(self)

    def pLog(self):
        """
        Plethystic logarithm along every ray.
        """
        return self.invar().pLog().poly()

    def pExp(self):
        """
        Plethystic logarithm along every ray.
        """
        return self.invar().pExp().poly()


class TMPolynomialRing(MPolynomialRing_polydict):
    """
    Truncated Multivariate polynomial ring.

    Multivariate polynomials up to degree ``prec``,
    where ``prec`` is a degree vector or an integer (total degree).
    Monomials of degree ``prec`` are included.
    We implement Adams operations, plethystic Exp and Log
    (we introduce the category of lambda rings
    :class:`msinvar.lambda_rings.LambdaRings` and embed our ring into it).

    This is the parent class for truncated multivariate polynomials.
    It has an alias :class:`TMPoly`.
    The element class is :class:`TMPolynomial`.

    PARAMETERS:

        1. ``base_ring`` -- the base of our polynomial algebra.
        2. ``n`` -- the number of vairables.
        3. ``names`` -- names of variables; can be just 'x' or 'x,y,..'.
        4. ``prec`` - precision vector (or integer, or None) for truncation
        5. ``order`` -- the order of variables.
        6. ``prod_twist`` - function that maps exponents (d,e) to the base_ring. The new product is x^d*x^e=prod_twist(d,e)*x^(d+e)

    EXAMPLES::

        sage: from msinvar import TMPoly
        sage: R=TMPoly(QQ,2,prec=(2,2)); R
        Multivariate Polynomial Ring in x0, x1 over Rational Field truncated at degree (2, 2)
        sage: x=R.gens(); (x[0]+x[1])**3
        3*x0^2*x1 + 3*x0*x1^2

        sage: QR=Frac(PolynomialRing(QQ,'y,t'))
        sage: S=TMPoly(QR,2,'x',prec=(2,2))
        sage: y,t=QR.gens(); x=S.gens()
        sage: (y*x[0]+x[1]).adams(2)
        y^2*x0^2 + x1^2
        sage: (y*x[0]).Exp()
        1 + y*x0 + y^2*x0^2

        sage: R=TMPoly(QQ,1,prec=2, prod_twist=lambda a,b:2)
        sage: x=R.gen(); x*x
        2*x^2
    """

    Element = TMPolynomial

    def __init__(self, base_ring=QQ, n=1, names='x', order='negdegrevlex',
                 prec=None, prod_twist=None):
        # Parent.__init__(self, category=LambdaRings())#not needed anymore
        self._prec = prec
        self.prod_twist = prod_twist
        super().__init__(base_ring, n, names, order)
        LambdaRings.add_ring(self)

    def adams(self, a, n):
        d = {e.emul(n): adams(c, n) for e, c in a.dict().items()
             if self.le_prec(e.emul(n))}
        return self.Element(self, d)

    def __call__(self, x):  # needed so that the result is MPT
        if isinstance(x, list) or isinstance(x, tuple):
            return self.Element(self, {tuple(x): 1})
        x = super().__call__(x)
        return self.Element(self, x.dict())

    def _poly_class(self):  # needed so that gens() are of type MPT
        return self.Element

    def _repr_(self):
        s = 'total degree ' if isinstance(self._prec, int) else 'degree '
        return super()._repr_()+' truncated at '+s+str(self._prec)

    def zero_tuple(self):
        return self._zero_tuple

    def prec(self, d=None):
        if d is not None:
            self._prec = d
            return
        return self._prec

    def prec_num(self):
        d = self._prec
        if d is None:
            return None
        if isiterable(d):
            return sum(d)
        return d

    def le_prec(self, v):
        d = self._prec
        if d is None:
            return True
        if isiterable(d):
            return all(i <= j for i, j in zip(v, d))
        return sum(v) <= d


TMPoly = TMPolynomialRing


def Exp(f):
    return f.Exp()


def Log(f):
    return f.Log()


def Psi(f):
    return f.Psi()


def IPsi(f):
    return f.IPsi()
