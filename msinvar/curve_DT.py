r"""
DT invariants of curves

We compute invariants of moduli stacks of semistable vector bundles on a curve
using the formula of 
`Zagier <https://people.mpim-bonn.mpg.de/zagier/files/mpim/94-5/fulltext.pdf>`_
in the form presented in :arxiv:`1310.4991`.

Then we compute integer DT invariants (by taking plethystic logarithm)
and invariants of stable moduli spaces following :arxiv:`0711.0634`.

EXAMPLES::

    sage: from msinvar.curve_DT import Curve
    sage: C=Curve(g=0,prec=[3,3]); C
    Curve of genus 0
    sage: C.intDT([1,0]).dict() #degree 0 and rank <=3 DT invariants
    {(1, 0): 1}

    sage: C=Curve(g=1,prec=[3,3]); C
    Curve of genus 1
    sage: C.intDT([1,0]).dict() #degree 0 and rank <=3 DT invariants
    {(1, 0): y^2 - 2*y + 1}

    sage: C=Curve(g=2,prec=[2,2]); C
    Curve of genus 2
    sage: C.intDT_val([1,0])
    y^4 - 4*y^3 + 6*y^2 - 4*y + 1   
    sage: C.intDT_val([2,0])
    y^10 - 4*y^9 + 7*y^8 - 8*y^7 + 8*y^6 - 8*y^5 + 8*y^4 - 8*y^3 + 7*y^2 - 4*y + 1
    sage: C.intDT_val([2,1])
    y^10 - 4*y^9 + 7*y^8 - 12*y^7 + 24*y^6 - 32*y^5 + 24*y^4 - 12*y^3 + 7*y^2 - 4*y + 1

    sage: C.stable_val([1,0]) # motivic classes of stable moduli spaces
    y^4 - 4*y^3 + 6*y^2 - 4*y + 1
    sage: C.stable_val([2,0]) # this value is different from int_DT
    y^10 - 4*y^9 + 6*y^8 - 4*y^7 - 4*y^6 + 20*y^5 - 30*y^4 + 20*y^3 - 5*y^2
    sage: C.stable_val([2,1])
    y^10 - 4*y^9 + 7*y^8 - 12*y^7 + 24*y^6 - 32*y^5 + 24*y^4 - 12*y^3 + 7*y^2 - 4*y + 1
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.misc_c import prod
from sage.functions.other import floor
from sage.arith.misc import gcd
from msinvar.iterators import OrderedPartitions_iterator
from msinvar import RF, TMPoly, WCS, Invariant, vec


class Curve(WCS):
    """
    Wall-crossing structure of a curve.

    It is a polynomial algebra in 2 variables (actually a quantum affine plane)
    over the field Q(y,t), where our zeta function of the curve lives.
    The zeta function is
    `(1-yt)^{2g}/(1-t)/(1-y^2t)`.
    """

    def __init__(self, g=0, prec=None):
        self.g = g  # genus
        self.base = RF('y,t')
        self.y = self.base.gen(0)
        self.t = self.base.gen(1)
        self.gm = 1/self.y - self.y
        self.q = self.y**2
        self.rank = 2
        self.R = TMPoly(self.base, self.rank, prec=prec)

    def __repr__(self):
        return "Curve of genus "+str(self.g)

    def eform(self, a, b):
        """
        Euler form between sheaves having Chern characters ``a``, ``b``,
        where a[0] is a rank, a[1] is a degree and similarly for ``b``.
        """
        return a[0]*b[1]-a[1]*b[0]+(1-self.g)*a[0]*b[0]

    def sform(self, a, b):
        """Skew-symmetric Euler form."""
        return 2*(a[0]*b[1]-a[1]*b[0])

    def P(self, **kw):
        """Numerator of the zeta function."""
        y, t, g = self.y, self.t, self.g
        p = (1-y*t)**(2*g)
        if len(kw) == 0:
            return p
        return p(**kw)

    def Z(self, **kw):
        """Zeta function of a curve.
        It has the form `(1-yt)^{2g}/(1-t)/(1-y^2t)`."""
        q, t = self.q, self.t
        z = self.P() / (1-t) / (1-q*t)
        if len(kw) == 0:
            return z
        return z(**kw)

    def Zhat(self, **kw):
        """Twsited version of :meth:`Z`."""
        t, g = self.t, self.g
        z = t**(1-g)*self.Z()
        if len(kw) == 0:
            return z
        return z(**kw)

    def Bun(self, r):
        """Motivic class of the moduli stack of rank ``r`` vector bundles
        and some fixed degree."""
        q = self.q
        return self.P(t=1)/(q-1)*prod(self.Z(t=q**i) for i in range(1, r))

    def Bun_tw(self, r):
        """Twisted version of :meth:`Bun`."""
        return self.Bun(r)*self.twist_T([r, 0])

    def stacky_val(self, r, d=None):
        """Motivic class of the stack of semistable vector bundles having
        rank ``r`` and degree ``d``. If ``d`` is None, then ``r`` should be a
        2-vector.
        """
        if d is None:
            r, d = r[0], r[1]
        q = self.q

        def twist(l):
            p = (r-l[-1])*d
            r1 = 0
            s = 1
            for i in range(len(l)-1):
                r1 += l[i]
                r2 = l[i]+l[i+1]
                p -= r2*floor(r1*d/r)
                s = s*(1-q**r2)
            return q**p/s

        val = 0
        for l in OrderedPartitions_iterator(r):
            val += twist(l)*prod(self.Bun_tw(i) for i in l)
        return val

    def stacky(self, r, d=None):
        """Invariant with values given by 
        motivic classes of stacks of semistable vector bundles along the ray
        of Chern characters proportional to (r,d)."""
        if d is None:
            r, d = r[0], r[1]
        k = gcd(r, d)
        r = r // k
        d = d // k
        if self.prec() is None:
            raise ValueError("prec should be specified.")
        dct = {(0, 0): 1}
        i = 1
        while vec.le([i*r, i*d], self.prec()):
            dct[i*r, i*d] = self.stacky_val(i*r, i*d)
            i += 1
        return Invariant(dct, self)

    def intDT(self, r, d=None):
        """Integer DT invariants along the ray of Chern characters
        proportional to (r,d)."""
        I = self.stk2int(self.stacky(r, d))
        return I.term_twist(lambda a: (-self.y)**(1-self.eform(a, a)))

    def intDT_val(self, r, d=None):
        """Integer DT invariant for the Chern character (r,d)."""
        if d is None:
            r, d = r[0], r[1]
        return self.intDT(r, d)([r, d])

    def stable(self, r, d=None):
        """Invariant with values given by motivic classes of moduli spaces
        of stable vector bundles along the ray of Chern characters proportional
        to (r,d)."""
        I = self.stacky(r, d)
        return self.stable_from_stacky(I)

    def stable_val(self, r, d=None):
        """Motivic class of the moduli space of stable vector bundles 
        having rank ``r`` and degree ``d`` (not necessarily coprime)."""
        if d is None:
            r, d = r[0], r[1]
        return self.stable(r, d)([r, d])
