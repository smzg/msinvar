r"""
<Short one-line summary that ends with no period>

<Paragraph description>

EXAMPLES::

<Lots and lots of examples>

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
    def __init__(self, **kw):
        super().__init__(**kw)
        A = matrix([self.path2vec(p) for p in self._potential])
        m = len(self._potential)
        A1 = A.augment(vector([1]*m))
        l = matrix(A1.right_kernel().basis()).columns()
        self.wzero = [0]*len(l[0])
        self.wt = {a: l[self.arrow_num(a)] for a in self.arrows()}

    def path2vec(self, p):
        v = [0]*self.arrow_num()
        for a in p:
            v[self.arrow_num(a)] = +1
        return v

    def ideal_dim(self, I):
        """Dimension vector of an ideal I that consists of atoms"""
        vert = self.vertices()
        d = {i: 0 for i in vert}
        for u in I:
            d[u.t] += 1
        return tuple(d[vert[i]] for i in range(len(vert)))


class Atom:
    def __init__(self, Q, t, wt=None):
        """Q is a BTQuiver, t is a target, wt is a weight"""
        self.Q = Q
        self.t = t
        if wt is None:
            wt = Q.wzero
        self.wt = wt

    def add_arrow(self, a):
        wt = vec.add(self.wt, self.Q.wt[a])
        t = self.Q.t(a)
        return Atom(self.Q, t, wt)

    def equal(self, u):
        return self.t == u.t and vec.equal(self.wt, u.wt)

    def __repr__(self):
        return "Atom "+str(self.t)+"-"+str(self.wt)


def add_new_arrows(Q, P, n):
    """
    Q is a quiver.

    P is a poset of atoms.

    n is the index where newly added atoms start.
    """
    n1 = len(P.vert)
    for k in range(n, n1):
        u = P.vert[k]
        for a in Q.arrows():
            if u.t == Q.s(a):
                v = u.add_arrow(a)
                for v1 in P.vert:
                    if v.equal(v1):
                        v = v1
                        break
                if v != v1:
                    P.vert.append(v)
                P.rel.append([u, v])
    return n1


def create_path_poset(Q, i, N=10):
    u = Atom(Q, i)
    P = Poset(rel=[], vert=[u])
    n = 0
    for k in range(N):
        n = add_new_arrows(Q, P, n)
    P.vert = set(P.vert)
    return P


def partition_func(Q, i, N=5):
    R = TMPoly(QQ, Q.vertex_num(), 'x', prec=N)
    P = create_path_poset(Q, i, N)
    dct = {}
    for I in P.ideals(N):
        d = tuple(Q.ideal_dim(I))
        if d not in dct:
            dct[d] = 0
        dct[d] += 1
    return R(dct)


def partition_func1(Q, i, N=5):
    R = TMPoly(QQ, 1, 'x', prec=N)
    P = create_path_poset(Q, i, N)
    dct = {}
    for I in P.ideals(N):
        d = (len(I),)
        if d not in dct:
            dct[d] = 0
        dct[d] += 1
    return R(dct)


# data base of some brane tilings from https://www.lpthe.jussieu.fr/~pioline/computing.html
BTD = [
    ["C^3", 'Phi[1,1,1]*Phi[1,1,2]*Phi[1,1,3]-Phi[1,1,1]*Phi[1,1,2]*Phi[1,1,3]'],
    ["Conifold=Y10", 'Phi[1,2,1]*Phi[1,2,2]*Phi[2,1,1]*Phi[2,1,2]-Phi[1,2,1]*Phi[1,2,2]*Phi[2,1,1]*Phi[2,1,2]'],
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
