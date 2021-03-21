r"""
Brane tilings and molten crystal counting

This code is based on :arxiv:`0809.0117`.
See also :arXiv:`2012.14358`

Given a brane tiling (a bipartite graph embedded in a torus) we associate with
it a quiver with potential (Q,W).
For any vertex `i\in Q_0`, we consider the set `\Delta_i` of paths starting
at i up to an equivalence relation induced by the potential W.
This set has a poset structure and can be interpreted as a crystal embedded
in `\mathbb R^3`.
For example, a brane tiling for `\mathbb C^3` corresponds to a quiver with one
vertex 0, three loops x,y,z and potential W=x[y,z].
The corresponding Jacobian algebra is `\mathbb C[x,y,z]` which has a basis 
parametrized by `\Delta_0=\mathbb N^3`, our crystal.

We define molten crystals as complements of finite ideals 
`I\subset \Delta_i`, meaning a subset I such that `u\le v` and `v\in I`
implies `u\in I`.
For any path `u\in\Delta_i`, let `t(u)\in Q_0` be the target vertex of u.
Define the dimension vector of I to be 
`\dim I=\sum_{u\in I}e_{t(u)}\in\mathbb N^{Q_0}`.
We compute the partition function
`Z_{\Delta_i}(x)=\sum_{I\subset\Delta_i}x^{\dim I}`
which is closely related to numerical DT invariants of (Q,W).

EXAMPLES::

    sage: from msinvar import *
    sage: from msinvar.brane_tilings import *
    sage: CQ=CyclicQuiver(1)
    sage: PQ=CQ.translation_PQ(); PQ
    Quiver with 1 vertices, 3 arrows and potential with 2 terms
    sage: Q=BTQuiver(potential=PQ._potential)
    sage: Z=partition_func(Q, 0, 8); Z
    1 + x + 3*x^2 + 6*x^3 + 13*x^4 + 24*x^5 + 48*x^6 + 86*x^7 + 160*x^8
    sage: Z.Log()
    x + 2*x^2 + 3*x^3 + 4*x^4 + 5*x^5 + 6*x^6 + 7*x^7 + 8*x^8

::
    
    sage: BTD[1]
    ['Conifold=Y10',
     'Phi[1,2,1]*Phi[2,1,1]*Phi[1,2,2]*Phi[2,1,2]-Phi[1,2,1]*Phi[2,1,2]*Phi[1,2,2]*Phi[2,1,1]']
    sage: Q=BTQuiver(potential=BTD[1][1]); Q
    Quiver with 2 vertices, 4 arrows and potential with 2 terms
    sage: Z=partition_func(Q, 1, 5); Z
    1 + x0 + 2*x0*x1 + 4*x0^2*x1 + x0*x1^2 + 2*x0^3*x1 + 8*x0^2*x1^2 + 14*x0^3*x1^2 + 4*x0^2*x1^3
    sage: Z.Log()
    x0 - x0^2 + 2*x0*x1 + 2*x0^2*x1 + x0*x1^2 + 4*x0^2*x1^2 + 3*x0^3*x1^2 + 2*x0^2*x1^3
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.matrix.constructor import matrix, vector
from msinvar.quivers import Quiver
from msinvar.utils import vec
from msinvar.posets import Poset
from msinvar.tm_polynomials import TMPoly
from sage.rings.rational_field import QQ


class BTQuiver(Quiver):
    """Quiver with potential associated with a brane tiling.

    Contains the method :meth:`wt` such that two paths are equivalent (equal
    in the Jacobian algebra if and only if their weights are equal).
    """

    def __init__(self, **kw):
        super().__init__(**kw)
        A = matrix([self.path2vec(p) for p in self._potential])
        m = len(self._potential)
        A1 = A.augment(vector([1]*m))
        l = matrix(A1.right_kernel().basis()).columns()
        self.wzero = [0]*len(l[0])
        self._wt = {a: l[self.arrow_num(a)] for a in self.arrows()}

    def wt(self, a):
        """Weight of an arrow or a path."""
        if isinstance(a, tuple):  # a is an arrow
            return self._wt[a]
        return vec.add(self._wt[x] for x in a)

    def path2vec(self, p):
        """Multiplicities of all arrows in a path."""
        v = [0]*self.arrow_num()
        for a in p:
            v[self.arrow_num(a)] = +1
        return v

    def ideal_dim(self, I):
        """Dimension vector of an ideal I that consists of atoms (paths)."""
        vert = self.vertices()
        d = {i: 0 for i in vert}
        for u in I:
            d[u.t] += 1
        return tuple(d[vert[i]] for i in range(len(vert)))

    def add_new_arrows(self, P, n):
        """
        Add new arrows to the atoms in the poset ``P``. 

        - ``P`` -- poset of atoms,
        - ``n`` -- index in ``P``, where atoms that require new arrows start.
        """
        n1 = len(P.vert)
        for k in range(n, n1):
            u = P.vert[k]
            for a in self.arrows():
                if u.t == self.s(a):
                    v = u.add_arrow(a)
                    for v1 in P.vert:
                        if v.equal(v1):
                            v = v1
                            break
                    if v != v1:
                        P.vert.append(v)
                    P.rel.append([u, v])
        return n1

    def create_path_poset(self, i, N=10):
        """Creat the poset of atoms (paths) that start at the vertex ``i``
        and have length <= ``N``."""
        u = Atom(self, i)
        P = Poset(rel=[], vert=[u])
        n = 0
        for k in range(N):
            n = self.add_new_arrows(P, n)
        P.vert = set(P.vert)
        return P


class Atom:
    """
    Atom -- equivalence class of a path in the Jacobian algebra.

    - ``Q`` -- a BTQuiver,
    - ``t`` -- a target,
    - ``wt`` -- a weight.
    """

    def __init__(self, Q, t, wt=None):
        self.Q = Q
        self.t = t
        if wt is None:
            wt = Q.wzero
        self.wt = wt

    def add_arrow(self, a):
        wt = vec.add(self.wt, self.Q.wt(a))
        t = self.Q.t(a)
        return Atom(self.Q, t, wt)

    def equal(self, u):
        return self.t == u.t and vec.equal(self.wt, u.wt)

    def __repr__(self):
        return "Atom "+str(self.t)+"-"+str(self.wt)


def partition_func(Q, i, N=5):
    """Partition function of ideals (molten crystals) in the crystal of atoms
    that start at the vertex ``i`` and have length <= ``N``.

    We keep track of the dimension vectors of ideals."""
    R = TMPoly(QQ, Q.vertex_num(), 'x', prec=N)
    P = Q.create_path_poset(i, N)
    dct = {}
    for I in P.ideals(N):
        d = tuple(Q.ideal_dim(I))
        if d not in dct:
            dct[d] = 0
        dct[d] += 1
    return R(dct)


def partition_func1(Q, i, N=5):
    """Partition function of ideals (molten crystals) in the crystal of atoms
    that start at the vertex ``i`` and have length <= ``N``.

    We keep track only of the sizes of ideals."""
    R = TMPoly(QQ, 1, 'x', prec=N)
    P = Q.create_path_poset(i, N)
    dct = {}
    for I in P.ideals(N):
        d = (len(I),)
        if d not in dct:
            dct[d] = 0
        dct[d] += 1
    return R(dct)


# data base of some brane tilings from https://www.lpthe.jussieu.fr/~pioline/computing.html
BTD = [
    ["C^3", 'Phi[1,1,1]*Phi[1,1,2]*Phi[1,1,3]-Phi[1,1,1]*Phi[1,1,3]*Phi[1,1,2]'],
    ["Conifold=Y10", 'Phi[1,2,1]*Phi[2,1,1]*Phi[1,2,2]*Phi[2,1,2]-Phi[1,2,1]*Phi[2,1,2]*Phi[1,2,2]*Phi[2,1,1]'],
    ["C^2xC/2", 'Phi[1,1,1]*Phi[1,2,2]*Phi[2,1,1]+Phi[1,1,1]*Phi[1,2,1]*Phi[2,1,2]-Phi[1,2,2]*Phi[2,1,1]*Phi[2,2,1]+Phi[1,2,1]*Phi[2,1,2]*Phi[2,2,1]'],
    ["C^2xC/3", 'Phi[1,1,1]*Phi[1,2,1]*Phi[2,1,1]+Phi[2,2,1]*Phi[2,3,1]*Phi[3,2,1]+Phi[1,3,1]*Phi[3,1,1]*Phi[3,3,1]-Phi[1,2,1]*Phi[2,1,1]*Phi[2,2,1]+Phi[1,1,1]*Phi[1,3,1]*Phi[3,1,1]+Phi[2,3,1]*Phi[3,2,1]*Phi[3,3,1]'],
    ["PdP6=C^3/2x2", 'Phi[1,3,1]*Phi[2,1,1]*Phi[3,2,1]+Phi[1,2,1]*Phi[2,4,1]*Phi[4,1,1]+Phi[2,3,1]*Phi[3,4,1]*Phi[4,2,1]+Phi[1,4,1]*Phi[3,1,1]*Phi[4,3,1]-Phi[1,2,1]*Phi[2,3,1]*Phi[3,1,1]+Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,4,1]*Phi[2,1,1]*Phi[4,2,1]+Phi[2,4,1]*Phi[3,2,1]*Phi[4,3,1]'],
    ["SPP=L121", 'Phi[1,1,1]*Phi[1,3,1]*Phi[3,1,1]+Phi[1,2,1]*Phi[2,1,1]*Phi[2,3,1]*Phi[3,2,1]-Phi[1,1,1]*Phi[1,2,1]*Phi[2,1,1]+Phi[1,3,1]*Phi[2,3,1]*Phi[3,1,1]*Phi[3,2,1]'],
    ["P2=C^3/(1,1,1)", 'Phi[1,2,2]*Phi[2,3,3]*Phi[3,1,1]+Phi[1,2,3]*Phi[2,3,1]*Phi[3,1,2]+Phi[1,2,1]*Phi[2,3,2]*Phi[3,1,3]-Phi[1,2,3]*Phi[2,3,2]*Phi[3,1,1]+Phi[1,2,1]*Phi[2,3,3]*Phi[3,1,2]+Phi[1,2,2]*Phi[2,3,1]*Phi[3,1,3]'],
    ["F0.1=P1xP1", 'Phi[1,2,2]*Phi[2,4,2]*Phi[4,1,1]+Phi[1,3,1]*Phi[3,4,2]*Phi[4,1,2]+Phi[1,3,2]*Phi[3,4,1]*Phi[4,1,3]+Phi[1,2,1]*Phi[2,4,1]*Phi[4,1,4]-Phi[1,3,2]*Phi[3,4,2]*Phi[4,1,1]+Phi[1,2,2]*Phi[2,4,1]*Phi[4,1,2]+Phi[1,2,1]*Phi[2,4,2]*Phi[4,1,3]+Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,4]'],
    ["F0.2=P1xP1", 'Phi[1,2,1]*Phi[2,3,2]*Phi[3,4,2]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,3,1]*Phi[3,4,2]*Phi[4,1,2]-Phi[1,2,2]*Phi[2,3,2]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,2,2]*Phi[2,3,1]*Phi[3,4,1]*Phi[4,1,2]'],
    ["F1=dP1=Y21=L312", 'Phi[1,2,1]*Phi[2,3,2]*Phi[3,4,2]*Phi[4,1,1]+Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,2]+Phi[2,3,1]*Phi[3,4,3]*Phi[4,2,1]-Phi[1,3,1]*Phi[3,4,3]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,3,1]*Phi[3,4,2]*Phi[4,1,2]+Phi[2,3,2]*Phi[3,4,1]*Phi[4,2,1]'],
    ["F2=C^3/(1,1,2)", 'Phi[1,2,1]*Phi[2,3,2]*Phi[3,1,1]+Phi[1,2,2]*Phi[2,4,1]*Phi[4,1,1]+Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,2]+Phi[2,3,1]*Phi[3,4,2]*Phi[4,2,1]-Phi[1,2,2]*Phi[2,3,1]*Phi[3,1,1]+Phi[1,3,1]*Phi[3,4,2]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,4,1]*Phi[4,1,2]+Phi[2,3,2]*Phi[3,4,1]*Phi[4,2,1]'],
    ["dP2.1", 'Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,4,2]*Phi[4,5,2]*Phi[5,1,1]+Phi[2,4,1]*Phi[4,5,1]*Phi[5,2,1]+Phi[3,4,2]*Phi[4,5,3]*Phi[5,3,1]-Phi[1,2,1]*Phi[2,4,1]*Phi[4,1,1]+Phi[1,3,1]*Phi[3,4,2]*Phi[4,5,1]*Phi[5,1,1]+Phi[2,4,2]*Phi[4,5,3]*Phi[5,2,1]+Phi[3,4,1]*Phi[4,5,2]*Phi[5,3,1]'],
    ["dP2.2", 'Phi[1,2,1]*Phi[2,3,2]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,2,2]*Phi[2,3,1]*Phi[3,5,1]*Phi[5,1,1]+Phi[2,4,1]*Phi[4,5,1]*Phi[5,2,1]-Phi[1,2,2]*Phi[2,4,1]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,3,1]*Phi[3,4,1]*Phi[4,5,1]*Phi[5,1,1]+Phi[2,3,2]*Phi[3,5,1]*Phi[5,2,1]'],
    ["PdP2", 'Phi[1,2,2]*Phi[2,3,1]*Phi[3,1,1]+Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,3,1]*Phi[3,5,1]*Phi[5,1,1]+Phi[2,3,1]*Phi[3,4,1]*Phi[4,5,1]*Phi[5,2,1]-Phi[1,2,1]*Phi[2,3,2]*Phi[3,1,1]+Phi[1,2,2]*Phi[2,4,1]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,4,1]*Phi[4,5,1]*Phi[5,1,1]+Phi[2,3,2]*Phi[3,5,1]*Phi[5,2,1]'],
    ["dP3.1", 'Phi[1,3,1]*Phi[3,5,1]*Phi[5,1,1]+Phi[1,2,1]*Phi[2,3,1]*Phi[3,4,1]*Phi[4,5,1]*Phi[5,6,1]*Phi[6,1,1]+Phi[2,4,1]*Phi[4,6,1]*Phi[6,2,1]-Phi[1,2,1]*Phi[2,4,1]*Phi[4,5,1]*Phi[5,1,1]+Phi[1,3,1]*Phi[3,4,1]*Phi[4,6,1]*Phi[6,1,1]+Phi[2,3,1]*Phi[3,5,1]*Phi[5,6,1]*Phi[6,2,1]'],
    ["dP3.2", 'Phi[1,2,2]*Phi[2,4,1]*Phi[4,1,1]+Phi[1,3,1]*Phi[3,5,1]*Phi[5,1,1]+Phi[1,2,1]*Phi[2,3,1]*Phi[3,4,1]*Phi[4,6,1]*Phi[6,1,1]+Phi[2,5,1]*Phi[5,6,1]*Phi[6,2,1]-Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,5,1]*Phi[5,1,1]+Phi[1,2,2]*Phi[2,3,1]*Phi[3,5,1]*Phi[5,6,1]*Phi[6,1,1]+Phi[2,4,1]*Phi[4,6,1]*Phi[6,2,1]'],
    ["dP3.3", 'Phi[1,3,1]*Phi[3,5,1]*Phi[5,1,1]+Phi[1,2,2]*Phi[2,4,1]*Phi[4,5,1]*Phi[5,1,1]+Phi[1,2,2]*Phi[2,3,1]*Phi[3,6,1]*Phi[6,1,2]+Phi[1,4,1]*Phi[4,6,1]*Phi[6,1,2]-Phi[1,2,1]*Phi[2,3,1]*Phi[3,5,1]*Phi[5,1,2]+Phi[1,4,1]*Phi[4,5,1]*Phi[5,1,2]+Phi[1,3,1]*Phi[3,6,1]*Phi[6,1,1]+Phi[1,2,1]*Phi[2,4,1]*Phi[4,6,1]*Phi[6,1,1]'],
    ["dP3.4", 'Phi[1,2,1]*Phi[2,4,1]*Phi[4,1,1]+Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,2]+Phi[1,2,3]*Phi[2,5,1]*Phi[5,1,1]+Phi[1,3,3]*Phi[3,5,1]*Phi[5,1,2]+Phi[1,2,2]*Phi[2,6,1]*Phi[6,1,1]+Phi[1,3,2]*Phi[3,6,1]*Phi[6,1,2]-Phi[1,3,2]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,2,2]*Phi[2,4,1]*Phi[4,1,2]+Phi[1,3,1]*Phi[3,5,1]*Phi[5,1,1]+Phi[1,2,1]*Phi[2,5,1]*Phi[5,1,2]+Phi[1,3,3]*Phi[3,6,1]*Phi[6,1,1]+Phi[1,2,3]*Phi[2,6,1]*Phi[6,1,2]'],
    ["L131", 'Phi[1,1,1]*Phi[1,2,1]*Phi[2,1,1]+Phi[2,2,1]*Phi[2,3,1]*Phi[3,2,1]+Phi[1,4,1]*Phi[3,4,1]*Phi[4,1,1]*Phi[4,3,1]-Phi[1,2,1]*Phi[2,1,1]*Phi[2,2,1]+Phi[1,1,1]*Phi[1,4,1]*Phi[4,1,1]+Phi[2,3,1]*Phi[3,2,1]*Phi[3,4,1]*Phi[4,3,1]'],
    ["L152", 'Phi[1,2,1]*Phi[2,3,1]*Phi[3,1,1]+Phi[2,3,2]*Phi[3,4,1]*Phi[4,2,1]+Phi[1,4,1]*Phi[4,5,1]*Phi[5,1,1]+Phi[1,2,2]*Phi[2,6,1]*Phi[6,1,1]+Phi[3,4,2]*Phi[4,6,1]*Phi[5,3,1]*Phi[6,5,1]-Phi[1,2,2]*Phi[2,3,2]*Phi[3,1,1]+Phi[2,3,1]*Phi[3,4,2]*Phi[4,2,1]+Phi[3,4,1]*Phi[4,5,1]*Phi[5,3,1]+Phi[1,4,1]*Phi[4,6,1]*Phi[6,1,1]+Phi[1,2,1]*Phi[2,6,1]*Phi[5,1,1]*Phi[6,5,1]'],
    ["C^3/(1,1,3)", 'Phi[1,2,1]*Phi[2,3,2]*Phi[3,1,1]+Phi[2,3,1]*Phi[3,4,2]*Phi[4,2,1]+Phi[1,2,2]*Phi[2,5,1]*Phi[5,1,1]+Phi[1,4,1]*Phi[4,5,1]*Phi[5,1,2]+Phi[3,4,1]*Phi[4,5,2]*Phi[5,3,1]-Phi[1,2,2]*Phi[2,3,1]*Phi[3,1,1]+Phi[2,3,2]*Phi[3,4,1]*Phi[4,2,1]+Phi[1,4,1]*Phi[4,5,2]*Phi[5,1,1]+Phi[1,2,1]*Phi[2,5,1]*Phi[5,1,2]+Phi[3,4,2]*Phi[4,5,1]*Phi[5,3,1]'],
    ["Y23=L153", 'Phi[1,4,1]*Phi[2,1,2]*Phi[3,2,1]*Phi[4,3,1]+Phi[3,5,1]*Phi[4,3,2]*Phi[5,4,1]+Phi[1,6,2]*Phi[2,1,1]*Phi[6,2,1]+Phi[4,6,1]*Phi[5,4,2]*Phi[6,5,1]+Phi[1,6,1]*Phi[5,1,1]*Phi[6,5,2]-Phi[1,4,1]*Phi[2,1,1]*Phi[3,2,1]*Phi[4,3,2]+Phi[3,5,1]*Phi[4,3,1]*Phi[5,4,2]+Phi[1,6,1]*Phi[2,1,2]*Phi[6,2,1]+Phi[1,6,2]*Phi[5,1,1]*Phi[6,5,1]+Phi[4,6,1]*Phi[5,4,1]*Phi[6,5,2]'],
    ["C^3/(1,1,4)", 'Phi[1,2,1]*Phi[2,3,2]*Phi[3,1,1]+Phi[2,3,1]*Phi[3,4,2]*Phi[4,2,1]+Phi[3,4,1]*Phi[4,5,2]*Phi[5,3,1]+Phi[1,2,2]*Phi[2,6,1]*Phi[6,1,1]+Phi[1,5,1]*Phi[5,6,1]*Phi[6,1,2]+Phi[4,5,1]*Phi[5,6,2]*Phi[6,4,1]-Phi[1,2,2]*Phi[2,3,1]*Phi[3,1,1]+Phi[2,3,2]*Phi[3,4,1]*Phi[4,2,1]+Phi[3,4,2]*Phi[4,5,1]*Phi[5,3,1]+Phi[1,5,1]*Phi[5,6,2]*Phi[6,1,1]+Phi[1,2,1]*Phi[2,6,1]*Phi[6,1,2]+Phi[4,5,2]*Phi[5,6,1]*Phi[6,4,1]'],
    ["PdP3a=C^3/(1,2,3)", 'Phi[1,2,1]*Phi[2,4,1]*Phi[4,1,1]+Phi[1,4,1]*Phi[4,5,1]*Phi[5,1,1]+Phi[2,3,1]*Phi[3,5,1]*Phi[5,2,1]+Phi[1,3,1]*Phi[3,6,1]*Phi[6,1,1]+Phi[2,5,1]*Phi[5,6,1]*Phi[6,2,1]+Phi[3,4,1]*Phi[4,6,1]*Phi[6,3,1]-Phi[1,3,1]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,5,1]*Phi[5,1,1]+Phi[2,4,1]*Phi[4,5,1]*Phi[5,2,1]+Phi[1,4,1]*Phi[4,6,1]*Phi[6,1,1]+Phi[2,3,1]*Phi[3,6,1]*Phi[6,2,1]+Phi[3,5,1]*Phi[5,6,1]*Phi[6,3,1]'],
    ["Y30", 'Phi[1,2,2]*Phi[2,3,1]*Phi[3,4,1]*Phi[4,1,1]+Phi[1,2,2]*Phi[2,5,1]*Phi[5,6,1]*Phi[6,1,1]+Phi[3,4,1]*Phi[4,5,1]*Phi[5,6,2]*Phi[6,3,1]-Phi[1,2,1]*Phi[2,3,1]*Phi[3,4,2]*Phi[4,1,1]+Phi[1,2,1]*Phi[2,5,1]*Phi[5,6,2]*Phi[6,1,1]+Phi[3,4,2]*Phi[4,5,1]*Phi[5,6,1]*Phi[6,3,1]'],
    ["Y31", 'Phi[1,2,2]*Phi[2,3,1]*Phi[3,4,1]*Phi[4,1,1]+Phi[3,4,2]*Phi[4,5,1]*Phi[5,6,1]*Phi[6,3,1]+Phi[5,6,2]*Phi[6,1,1]*Phi[1,5,1]+Phi[6,1,2]*Phi[1,2,1]*Phi[2,6,1]-Phi[1,2,1]*Phi[2,3,1]*Phi[3,4,2]*Phi[4,1,1]+Phi[3,4,1]*Phi[4,5,1]*Phi[5,6,2]*Phi[6,3,1]+Phi[5,6,1]*Phi[6,1,2]*Phi[1,5,1]+Phi[6,1,1]*Phi[1,2,2]*Phi[2,6,1]']
]
