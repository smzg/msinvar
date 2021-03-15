r"""
Compute refined invariants for some quivers with potentials.

1. Compute (total refined) invariants for the cyclic
quiver C_n with the potential consisting of one term (the cycle).

2. Compute (total refined) invariants for the McKay quiver of C^3/Z_r with
the induced potential (consisting of 3-cycles), where the action is given
by (1, k, -k-1) for any 0<=k<r. This potential quiver can be obtained as a
translation potential quiver of a cyclic quiver:

EXAMPLES::

    sage: from msinvar.quivers import CyclicQuiver
    ....: from msinvar.potential_quiver_invar import *
    ....: r, k = 3, 1
    ....: CQ = CyclicQuiver(r); CQ
    Quiver with 3 vertices and 3 arrows
    sage: W=CQ.wcs([2]*r); W
    Wall-crossing structure on a lattice of rank 3
    sage: W.intAtt().dict() #Attractor invar for a cyclic quiver without potential
    {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1, (1, 1, 1): -y}

    sage: W.total=cyclic_potential_total(W) #Total invar for a cyclic quiver with po
    ....: tential
    ....: W.intAtt().dict() #Attractor invar for a cyclic quiver with potential 
    {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}

    sage: Q = CQ.translation_PQ(k); Q # translation quiver with potential
    Quiver with 3 vertices, 9 arrows and potential with 6 terms
    sage: W = Q.wcs([2]*r)
    ....: W.total=cyclic_TPQ_total(W, k)
    ....: W.intAtt().simp().dict() #Attractor invar for the translation quiver with
    ....:  potential
    {(1, 0, 0): 1,
     (0, 1, 0): 1,
     (0, 0, 1): 1,
     (1, 1, 1): (-2*y^2 - 1)/y,
     (2, 2, 2): (-2*y^2 - 1)/y}

"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from numpy import array
from msinvar.iterators import Multiplicities_iterator
from msinvar.utils import phi
from msinvar.invariants import Invariant
from msinvar.quivers import Quiver, CyclicQuiver, ChainQuiver


def interaction_invariant(S, si, dim, W):
    r"""
    Return expression of the form:

    :math:`\sum_{m:S\to N}
    (-y)^{-\sum_{i,j\in S} m_im_j\sigma(i,j)}/(1/y^2)_m\cdot
    x^{\sum_{i\in S} m_i \dim(i)}`

    1. W is a WCS, y=W.y, r=W.rank.
    2. S is a list.
    3. :math:`\sigma:S\times S\to Z` is a map (interaction form).
    4. :math:`\dim:S\to Z^r` is a map (dimension vectors).
    5. :math:`(q)_m=\prod_{i\in S}(q)_{m_i},\ (q)_k=(1-q)...(1-q^k)`.
    6. The sum runs over m such that :math:`\sum_i m_i \dim(i)\le` W.prec().
    """
    prec = W.prec()
    y = W.y
    if not isinstance(S, list):
        S = list(S(prec))
    Sdim = array([dim(a) for a in S]).T  # dimensions of indecomposables
    Si = array([[si(i, j) for j in S] for i in S])  # matrix of si-values

    dct = {tuple([0]*len(prec)): 1}
    for m in Multiplicities_iterator(Sdim, prec):
        d = tuple(Sdim@m)
        p = -Si@m@m
        if d not in dct:
            dct[d] = 0
        dct[d] += (-y)**p/phi(1/y**2, m)
    return Invariant(dct, W)


def cyclic_TPQ_total(W, k, prec=None, Q=None):
    """
    Return total stacky invariants (for the trivial stability) for the translation
    potential quiver of a cyclic quiver (with the translation given by the k-shift).

    W is the potential quiver or its wall-crossing structure.
    """
    if isinstance(W, Quiver):
        W = W.wcs()
    W.prec(prec)
    r = W.rank
    if Q is None:
        Q = CyclicQuiver(r)

    def rho(a, b):
        return sum(a[i]*(b[i]-b[(i+1) % r]-b[(i+k) % r]
                         + b[(i-k-1) % r]) for i in range(r))

    def si(a, b):
        return 2*Q.hom(a, b)-2*Q.hom(a, Q.shift_indec(b, k))\
            - rho(Q.dim(a), Q.dim(b))

    return interaction_invariant(Q.indecs, si, Q.dim, W)


def cyclic_potential_total(W, prec=None):
    """
    Return total stacky invariants (for the trivial stability) for a cyclic
    quiver with the cyclic potential.

    W is the cyclic quiver or its wall-crossing structure.
    """
    if isinstance(W, Quiver):
        W = W.wcs()
    W.prec(prec)
    r = W.rank
    Q = ChainQuiver(r)
    S = list(Q.indecs())
    S.remove((0, r-1))  # this representation doesn't satisfy relations

    def si(a, b):
        da, db = Q.dim(a), Q.dim(b)
        return 2*Q.hom(a, b)-2*da[0]*db[-1]-W.eform(da, db)

    return interaction_invariant(S, si, Q.dim, W)
