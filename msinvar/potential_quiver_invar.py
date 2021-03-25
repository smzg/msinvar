r"""
Compute refined invariants for some quivers with potentials.

1. Compute (total refined) invariants for the cyclic
quiver C_n with the potential consisting of one term (the cycle).

2. Compute (total refined) invariants for the McKay quiver of C^3/Z_r with
the induced potential (consisting of 3-cycles), where the action is given
by (1, k, -k-1) for any 0<=k<r. This potential quiver can be obtained as a
translation potential quiver of a cyclic quiver, see :arxiv:`1911.01788`.
Some of the results are presented in :arxiv:`2012.14358`. 

EXAMPLES::

    sage: from msinvar.quivers import CyclicQuiver
    sage: from msinvar.potential_quiver_invar import *
    sage: r, k = 3, 1
    sage: CQ = CyclicQuiver(r); CQ
    Cyclic quiver: Quiver with 3 vertices and 3 arrows
    sage: CQ.prec([2]*r) # precision vector
    sage: CQ.intAtt().dict() #Attractor invar for a cyclic quiver without potential
    {(0, 0, 1): 1, (0, 1, 0): 1, (1, 0, 0): 1, (1, 1, 1): -y}

::
    
    sage: total=cyclic_potential_total(CQ) # Total invar for a cyclic quiver with potential
    sage: CQ.intAtt(total).dict() #Attractor invar for a cyclic quiver with potential 
    {(0, 0, 1): 1, (0, 1, 0): 1, (1, 0, 0): 1}

::
    
    sage: PQ = CQ.translation_PQ(1); PQ # translation quiver with potential
    Translation PQ: Quiver with 3 vertices, 9 arrows and potential with 6 terms
    sage: PQ.prec([2,2,2])
    sage: total=translation_PQ_total(PQ)
    sage: PQ.intAtt(total).simp().dict() #Attractor invar for the translation quiver with potential
    {(0, 0, 1): 1,
     (0, 1, 0): 1,
     (1, 0, 0): 1,
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
from msinvar.quivers import Quiver, ChainQuiver


def interaction_invariant(S, si, dim, Q):
    r"""
    Return expression of the form:

    .. math::
        \sum_{m:S\to N} (-y)^{-\sum_{i,j\in S} m_im_j\sigma(i,j)}
        /(1/y^2)_m\cdot x^{\sum_{i\in S} m_i \dim(i)}

    1. Q is a Quiver, y=Q.y, r=Q.rank.
    2. S is a list.
    3. `\sigma:S\times S\to Z` is a map (interaction form).
    4. `\dim:S\to Z^r` is a map (dimension vectors).
    5. `(q)_m=\prod_{i\in S}(q)_{m_i},\ (q)_k=(1-q)...(1-q^k)`.
    6. The sum runs over m such that `\sum_i m_i \dim(i)\le` Q.prec().
    """
    prec = Q.prec()
    y = Q.y
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
    return Invariant(dct, Q)


def translation_PQ_total(PQ, prec=None):
    """
    Return total invariant (stacky invariant for the trivial stability)
    for the translation potential quiver, assuming that the base quiver has 
    implemented methods :meth:`ind_list`, :meth:`ind_dim`, :meth:`ind_hom`.

    ``PQ`` -- translation potential quiver or its wall-crossing structure.

    See :class:`msinvar.quivers.TranslationPQ`.

    EXAMPLES::

        sage: from msinvar import *
        sage: from msinvar.potential_quiver_invar import *
        sage: Q=CyclicQuiver(3); Q
        Cyclic quiver: Quiver with 3 vertices and 3 arrows
        sage: PQ=Q.translation_PQ(1); PQ
        Translation PQ: Quiver with 3 vertices, 9 arrows and potential with 6 terms
        sage: I=translation_PQ_total(PQ, [2,2,2])
        sage: PQ.intAtt(I).simp().dict()
        {(0, 0, 1): 1,
         (0, 1, 0): 1,
         (1, 0, 0): 1,
         (1, 1, 1): (-2*y^2 - 1)/y,
         (2, 2, 2): (-2*y^2 - 1)/y}
    """
    PQ.prec(prec)
    Q = PQ.base_quiver
    hom = Q.ind_hom
    dim = Q.ind_dim

    def rho(a, b):
        s = PQ.eform(a, b)
        return s + sum([2 * a[PQ.tau(e[1])] * b[e[0]] for e in Q.arrows()])

    def si(a, b):
        return 2*hom(a, b)-2*hom(PQ.ind_tau(a), b) - rho(dim(a), dim(b))

    return interaction_invariant(Q.ind_list, si, dim, PQ)


# cyclic_TPQ_total = translation_PQ_total  # backwards compatibility


def cyclic_potential_total(CQ, prec=None):
    """
    Return total stacky invariants (for the trivial stability) for a cyclic
    quiver with the cyclic potential.

    CQ is the cyclic quiver or its wall-crossing structure.
    """
    if not isinstance(CQ, Quiver):
        CQ = CQ.quiver
    CQ.prec(prec)
    r = CQ.vertex_num()
    Q = ChainQuiver(r)
    S = Q.ind_list()
    S.remove((0, r-1))  # this representation doesn't satisfy relations

    def si(a, b):
        da, db = Q.ind_dim(a), Q.ind_dim(b)
        return 2*Q.ind_hom(a, b)-2*da[0]*db[-1]-CQ.eform(da, db)

    return interaction_invariant(S, si, Q.ind_dim, CQ)


class QuiverExample1(Quiver):
    """Quiver 1->0, 1->2 with methods for indecomposables and a non-trivial
    automorphism.

    EXAMPLES::

        sage: from msinvar.potential_quiver_invar import *
        sage: Q=QuiverExample1(); Q
        Quiver with 3 vertices and 2 arrows
        sage: PQ=Q.translation_PQ(); PQ
        Translation PQ: Quiver with 3 vertices, 7 arrows and potential with 4 terms
        sage: PQ.prec([3,3,3])
        sage: I=PQ.translation_PQ_total()
        sage: PQ.intAtt(I).dict()
        {(0, 0, 1): 1,
         (0, 1, 0): -y,
         (0, 1, 1): 1,
         (1, 0, 0): 1,
         (1, 0, 1): -y,
         (1, 1, 0): 1,
         (1, 1, 1): -y,
         (1, 2, 1): -y}
    """

    def __init__(self):
        super().__init__('1-0,1-2')

    def ind_list(self, *args):
        """
        Return the list of all indecomposable representations, paramertized by
        pairs (i,j) with i<=j.
        """
        n = self.vertex_num()
        l = []
        for i in range(n):
            l += [(i, j) for j in range(i, n)]
        return l

    def ind_dim(self, a):
        n = self.vertex_num()
        i, j = a[0], a[1]
        return [0]*i+[1]*(j-i+1)+[0]*(n-j-1)

    def ind_hom(self, a, b):
        dim = self.ind_dim
        if a in {(0, 2), (0, 0), (2, 2)} or b in {(0, 1), (1, 1), (1, 2)}:
            # a is projective or b is injective
            return self.eform(dim(a), dim(b))
        # if b=(0,0) or (2,2), we get 0 as a is generated at the vertex 1
        # if b=(0,2), we also get 0 -- direct check
        return 0

    def tau(self, a):
        """Non-trivial automorphism of the quiver."""
        if isinstance(a, tuple):  # an arrow
            return (2-a[0], a[1], a[2])
        return 2-a  # a vertex

    def ind_tau(self, a):
        """Bijection on indecomposables induced by :meth:`tau`."""
        return (2-a[1], 2-a[0])

    def translation_PQ(self):
        from msinvar import TranslationPQ
        return TranslationPQ(self, self.tau, self.ind_tau)
