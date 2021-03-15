r"""
Motivic classes of the moduli spaces of twisted Higgs bundles

EXAMPLES::
    
sage: from msinvar.higgs_bundles import CurveAlgebra
sage: from msinvar.higgs_bundles import twisted_higgs_bundles_invariant as invar
sage: C=CurveAlgebra(g=2)
sage: invar(C,l=2,r=2).factor()
(y - 1)^4 * y^10 * (y^2 + 1) * (y^4 - 4*y + 2)
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.combinat.partition import Partitions
from msinvar.tm_polynomials import TMPoly
from msinvar.rings import RationalFunctionField


class CurveAlgebra(RationalFunctionField):
    """
    The algebra where zeta function of a curve lives.

    For simplicity we define the algebra to be Q(y,t) or Q(u,v,t).
    In the first case the zeta function is
    (1-y*t)^(2*g)/(1-t)/(1-y^2*t)
    """

    def __init__(self, g=0, vars='y'):
        self.g = g  # genus
        vars = vars+',t'
        super().__init__(vars)
        if self.ngens() > 3 or self.ngens() < 2:
            raise ValueError("Wrong number of variables")
        self.u = self.gen(0)
        n = self.ngens()
        self.t = self.gen(n-1)
        self.v = self.gen(n-2)
        self.q = self.u*self.v

    def zeta(self, **kw):
        u, v, t, g = self.u, self.v, self.t, self.g
        z = (1-u*t)**g * (1-v*t)**g / (1-t) / (1-u*v*t)
        if len(kw) == 0:
            return z
        return z(**kw)

    def H(self, part, p):
        """
        Compute the function H from https://arxiv.org/abs/1104.5698 (7)
        """
        t, g, Z = self.t, self.g, self.zeta
        q = self.u * self.v
        s = 1
        for i, j in part.cells():
            a = part.arm_length(i, j)
            l = part.leg_length(i, j)
            h = a+l+1
            s = s*(-t**(a-l)*q**a)**p*t**((1-g)*(2*l+1))*Z(t=t**h*q**a)
        return s


################################

def twisted_higgs_bundles_invariant(C, l, r):
    """
    Invariants of the moduli space of l-twisted rank r Higgs bundles
    over a curve C.

    Expected dimension of the moduli space is 1+l*r^2.
    For l=2g-2 it is 2+l*r^2
    The final power twist is given by arxiv.org/pdf/1901.02439.pdf (6)
    """
    u, v, t, g = C.u, C.v, C.t, C.g
    p = l-(2*g-2)
    R = TMPoly(C, 1, 'T', prec=r)
    T = R.gen()
    F = 0
    for i in range(r+1):
        for part in Partitions(i):
            F += C.H(part, p)*T**i
    H = F.Log().coeff([r])*(1-t)*(1-u*v*t)
    M = H(t=1)*(-1)**(p*r)*(u*v)**((g-1)*r**2+p*r*(r+1)/2)
    if p == 0:
        return M*u*v  # the dimension of the moduli space increases by 1 for p=0
    return M
