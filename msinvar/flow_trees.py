r"""
This module implements 
attractor tree formula (see :arxiv:`1910.03098` and :arxiv:`2101.07636`)
and flow tree formula (see :arxiv:`1804.06928` and :arxiv:`2102.11200`).

EXAMPLES::

    sage: from msinvar import *
    sage: Q=KroneckerQuiver(2)
    sage: W=Q.wcs([2,2])
    sage: z=Stability([1,0])
    sage: OmbAtt=W.ratAtt_default() #rational attractor invariant
    sage: OmbAtt.dict()
    {(0, 1): 1, (0, 2): (-y)/(-2*y^2 - 2), (1, 0): 1, (2, 0): (-y)/(-2*y^2 - 2)}

First we apply attractor tree formula to find rational DT invariants for the
above stability z::

    sage: Omb1=attr_tree_formula(W, z, OmbAtt)
    sage: Omb1.simp().dict()
    {(0, 1): 1,
     (0, 2): 1/2*y/(y^2 + 1),
     (1, 0): 1,
     (1, 1): (-y^2 - 1)/y,
     (1, 2): 1,
     (2, 0): 1/2*y/(y^2 + 1),
     (2, 1): 1,
     (2, 2): (-1/2*y^4 - 1/2)/(y^3 + y)}
    
Next we apply flow tree formula::
    
    sage: Omb2=flow_tree_formula(W, z, OmbAtt)
    sage: Omb2.simp().dict()
    {(0, 1): 1,
     (0, 2): 1/2*y/(y^2 + 1),
     (1, 0): 1,
     (1, 1): (-y^2 - 1)/y,
     (1, 2): 1,
     (2, 0): 1/2*y/(y^2 + 1),
     (2, 1): 1,
     (2, 2): (-1/2*y^4 - 1/2)/(y^3 + y)}
    
Finally, we apply the wall-crossing formula to determine the same invariant
from the total invariant (stacky invariant for the trivial stability)::    
    
    sage: Omb3=W.rat_from_total(z, W.total())
    sage: Omb3.dict()
    {(0, 1): 1,
     (0, 2): 1/2*y/(y^2 + 1),
     (1, 0): 1,
     (1, 1): (-y^2 - 1)/y,
     (1, 2): 1,
     (2, 0): 1/2*y/(y^2 + 1),
     (2, 1): 1,
     (2, 2): (-1/2*y^4 - 1/2)/(y^3 + y)}
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.misc_c import prod
from sage.functions.other import factorial
from sage.combinat.permutation import Permutations
from sage.misc.prandom import random
from sage.rings.rational_field import QQ
from msinvar.utils import vec
from msinvar.invariants import Transform, Invariant
from msinvar.iterators import UnorderedMultiPartitions_iterator
from msinvar.stability import Stability


def attr_tree_formula(W, z, I, t=0):
    r"""Attractor tree formula, conjectured in :arxiv:`1910.03098` and proved
    in :arxiv:`2101.07636`.

    - ``W`` -- wall-crossing structure :class:`WCS`.
    - ``z`` -- stability parameter.
    - ``I`` -- rational attractor invariant `\bar\Omega_*`.
    - ``t`` -- rational number.
    """

    T = attr_tree_transform(W, z, t)
    T.twist = W.twist
    I1 = Invariant(lambda d: I(d)/W.gm, W)
    I2 = T(I1)
    return Invariant(lambda d: I2(d)*W.gm, W)


def attr_tree_transform(W, z, t=0):
    z = Stability.check(z)
    Tz = transform_g(z, t)
    Tstar = transform_g(W.sform, t)
    return Tz.plethysm(Tstar.inverse())


def transform_g(f, t=0):
    """Auxiliary transform (collection) g from :arxiv:`2101.07636`."""
    def F(l):
        if len(l) == 1:
            return 1
        tot = vec.add(l)
        v = [0]*len(tot)
        n1 = 0
        n0 = 0
        for i in range(len(l)-1):
            v = vec.add(v, l[i])
            c = f(v, tot)
            if t == 0 and c < 0:
                return 0
            if t == 1 and c > 0:
                return 0
            if c > 0:
                n1 += 1
            elif c == 0:
                n0 += 1
        nm = len(l)-1-n1-n0
        return QQ(t**nm*(t-1)**n1*(t**(n0+1)-(t-1)**(n0+1))/(n0+1))
    return Transform(F)


def flow_tree_formula(W, z, I, d=None):
    r"""Flow tree formula, conjectured in :arxiv:`1804.06928` and proved
    in :arxiv:`2102.11200`. We use its modification from :arxiv:`2101.07636`.

    - ``W`` -- wall-crossing structure :class:`WCS`.
    - ``z`` -- stability parameter.
    - ``I`` -- rational attractor invariant `\bar\Omega_*`.
    """
    z = Stability.check(z)
    if d is None:
        return Invariant(lambda d: flow_tree_formula(W, z, I, d), W)
    if vec.zero(d):
        return 0
    def kp(m):
        y = W.y
        return (1/y**m-y**m)/(1/y-y)*(-1)**(m-1)

    def mult_factorial(l):
        m = {}
        for i in l:
            i = tuple(i)
            if not i in m:
                m[i] = 0
            m[i] += 1
        return prod(factorial(i) for i in m.values())

    te = z.normalize(d)
    A = 0
    for l in UnorderedMultiPartitions_iterator(d):
        M = [[W.sform(i, j) for j in l] for i in l]
        te1 = [te(i)+random()*1e-5 for i in l]
        te1[0] -= sum(te1)
        s = 0
        for p in Permutations(len(l)):
            p1 = [i-1 for i in p]
            s += permIndex(p1, te1, M, kp)
        A += s*prod(I(e) for e in l)/mult_factorial(l)
    return A


def permIndex(l, te, M, kp):
    """Permutation index used in the flow tree formula."""
    def prod1(M, l1, l2):
        return sum(M[i][j] for i in l1 for j in l2)
    if len(l) == 1:
        return 1
    F = 0
    for i in range(1, len(l)):
        l1 = l[:i]
        l2 = l[i:]
        m = prod1(M, l1, l2)
        c = sum(te[i] for i in l1)
        c = round(c, 10)
        if m <= 0 or c >= 0:
            continue
        te1 = [te[j]-c/m*prod1(M, [j], l) for j in range(len(te))]
        F += kp(m)*permIndex(l1, te1, M, kp)*permIndex(l2, te1, M, kp)
    return F
