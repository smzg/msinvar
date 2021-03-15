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

import re
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from msinvar.utils import vec
from msinvar.wall_crossing import WCS


class Quiver(DiGraph):
    def __init__(self, data=None, potential=None, prec=None,
                 pos=None, loops=None, format=None, weighted=False, data_structure=None, name=None):
        if data is not None and isinstance(data, str):
            data = get_arrows_from_string(data)
        if potential is not None:
            if isinstance(potential, str):
                self._potential = convert_potential_string(potential)
            else:
                self._potential = potential
            if data is None:
                data = get_arrows_from_potential(self._potential)
        else:
            self._potential = None
        super().__init__(data=data, multiedges=True, loops=True)
        self._arrows = []
        vert = self.vertices()
        self.vertex_dict = {vert[i]: i for i in range(self.vertex_num())}
        self.set_quiver_dict()
        self._wcs = WCS(self, prec=prec)

    def wcs(self, prec=None):
        if prec is not None:
            self._wcs.prec(prec)
        return self._wcs

    def prec(self, prec=None):
        return self._wcs.prec(prec)

    def arrows(self, n=None):
        if not n is None:
            return self._arrows[n]
        if len(self._arrows) != self.arrow_num():
            self._arrows = list(self.edge_iterator())
        return self._arrows

    def potential(self): return self._potential

    def arrow_num(self, a=None):
        if a is None:
            return self.size()
        return self.arrows().index(a)

    def vertex_num(self, i=None):
        if i is None:
            return self.order()
        return self.vertex_dict[i]

    def s(self, a):
        """Source of an arrow a; a can be an arrow or its number"""
        if isinstance(a, int):
            a = self.arrows(a)
        return a[0]

    def t(self, a):
        """Target of an arrow a; a can be an arrow or its number"""
        if isinstance(a, int):
            a = self.arrows(a)
        return a[1]

    def set_quiver_dict(self):
        d = {}
        for a in self.arrows():
            i, j = self.vertex_num(a[0]), self.vertex_num(a[1])
            if (i, j) in d:
                d[(i, j)] += 1
            else:
                d[(i, j)] = 1
        self.quiver_dict = d

    def eform(self, a, b):
        """Euler form"""
        s = sum(a[i]*b[i] for i in range(len(a)))
        for (i, j), k in self.quiver_dict.items():
            s -= k * a[i]*b[j]
        return s

    def sform(self, a, b):
        """Skew-symmetrized Euler form"""
        s = 0
        for (i, j), k in self.quiver_dict.items():
            s += k * (a[j]*b[i]-a[i]*b[j])
        return s

    def get_eform_matrix(self):
        n = self.vertex_num()
        M = matrix.identity(n)
        for (i, j), k in self.quiver_dict.items():
            M[(i, j)] -= k
        return M

    def get_sform_matrix(self):
        n = self.vertex_num()
        M = matrix(n)
        for (i, j), k in self.quiver_dict.items():
            M[(i, j)] -= k
            M[(j, i)] += k
        return M

    def _repr_(self):
        sv = "Quiver with "+str(self.vertex_num())+" vertices"
        sa = str(self.arrow_num())+" arrows"
        W = self.potential()
        if W is not None:
            return sv+", "+sa+" and potential with "+str(len(W))+" terms"
        return sv+" and "+sa

    def translation_PQ(self, tau=None):
        return translation_PQ(self, tau)

    def Ginzburg_PQ(self):
        return self.translation_PQ()


def convert_potential_string(s):
    """
    Convert a potential string to a list. We drop the coefficients.
    Potential is of the form

    F[1,2,1]*F[2,3]F[3,1]...-F[1,2,2]F[2,1]...+...
    """
    l = re.split('\\+|\\-', s)
    return list(convert_path_string(t) for t in l)


def convert_path_string(s):
    """
    Convert a path string to a list of arrows.
    A path string is of the form (only the content inside square brackets
    is important)

    F[1,2,1]*F[2,3]a[1,2,2]b[2,2]
    """
    l = re.findall(".*?\\[(.*?)\\]", s)
    # list of strings between [,]; use .*? for a non greeedy search
    return list(convert_arrow_string(t) for t in l)


def convert_arrow_string(s):
    """
    Convert a string to a tuple of integers. Replace pairs (i,j) by triples (i,j,1)
    """
    a = tuple(int(i) for i in s.split(','))
    return (*a, 1) if len(a) == 2 else a


def get_arrows_from_potential(W):
    """List of arrows from a potential list W"""
    arrs = set()
    for p in W:
        arrs.update(p)
    return list(arrs)


def get_arrows_from_string(s):
    arrs = set()
    for p in s.split(','):
        p = p.split('-')
        l = list(int(re.sub(r'\D', '', i)) for i in p)  # remove non-digits
        for i in range(len(l)-1):
            arrs.add((l[i], l[i+1], 1))
    return list(arrs)


def KroneckerQuiver(m=2):
    l = [[1, 2, i] for i in range(m)]
    return Quiver(l)


def JordanQuiver(m=1):
    l = [[1, 1, i] for i in range(m)]
    return Quiver(l)


class ChainQuiver(Quiver):
    def __init__(self, n, prec=None):
        l = [[i, i+1, 1] for i in range(n-1)]
        super().__init__(l, prec=prec)

    def dim(self, a):
        """
        Return the dimension vector of an indecomposable representation
        parametrized by a pair a=(i,j) with i<=j.
        """
        n = self.vertex_num()
        i, j = a[0], a[1]
        return [0]*i+[1]*(j-i+1)+[0]*(n-j-1)

    def hom(self, a, b):
        """
        Return the dimension of the Hom-space between indecomposable
        representations parametrized by pairs a, b.
        """
        i, j, k, l = a[0], a[1], b[0], b[1]
        if k <= i and i <= l and j >= l:
            return 1
        return 0

    def indecs(self):
        """
        Return all indecomposable representations (that is, pairs
        (i,j) representing them).
        """
        n = self.vertex_num()
        for i in range(n):
            for j in range(i, n):
                yield (i, j)


class CyclicQuiver(Quiver):
    def __init__(self, n, prec=None):
        l = [[i, i+1, 1] for i in range(n-1)]+[[n-1, 0, 1]]
        super().__init__(l, prec=prec)

    def dim(self, a):
        """
        Return the dimension vector of an indecomposable representation
        parametrized by a triple a=(i,j,n).
        """
        r = self.vertex_num()
        i, j, n = tuple(a)
        d = [n]*r
        for k in range((j-i) % r + 1):
            d[(i+k) % r] += 1
        return d

    def hom(self, a, b):
        """
        Return the dimension of the Hom-space between indecomposable
        representations parametrized by triples a, b.
        """
        if a[2] == -1 or b[2] == -1:
            return 0
        r = self.vertex_num()
        i, j, m = tuple(a)
        k, l, n = tuple(b)
        d1 = (j-i) % r+r*m
        d2 = (l-k) % r+r*n
        return self.cong_number(max(0, d2-d1), d2, i-k, r)

    def indecs(self, d):
        """
        Return nilpotent indecomposables up to dimension vector d.
        """
        r = self.vertex_num()
        m = max(d)
        for i in range(r):
            for j in range(r):
                for n in range(m):
                    a = (i, j, n)
                    if vec.le(self.dim(a), d):
                        yield a

    def shift_indec(self, a, k):
        """
        Cyclic shift if indecomposable modules: (SM)_i=M_{i+k}.
        """
        r = self.vertex_num()
        return ((a[0]-k) % r, (a[1]-k) % r, a[2])

    def translation(self, k=1):
        """
        Return a quiver automorphism given by the k-cyclic shift.
        """
        r = self.vertex_num()

        def tau(a):
            if isinstance(a, tuple):  # an arrow
                return ((a[0]+k) % r, (a[1]+k) % r, 1)
            return (a+k) % r
        return tau

    def translation_PQ(self, k=1):
        tau = self.translation(k)
        return translation_PQ(self, tau)

    def cong_number(self, a, b, i, r):
        """
        Return the bumber of a<=s<=b with s=i mod r.
        """
        if a > b:
            return 0
        q = (b-a) // r
        if (b-a) % r < (i-a) % r:
            return q
        return q+1


def translation_PQ(Q, tau=None):
    """
    Construct a new quiver with potential based on a quiver Q and its
    automorphism tau. See arXiv:1911.01788.

    1. For any arrow a:i->j we add an arrow tau(j)->i.
    2. For any vertex i we add an arrow i->tau(i).
    3. We consider potential consisting of 3-cycles (for every arrow i->j)
        1. i->j->tau(j)->i
        2. -i->tau(i)->tau(j)->i
    """
    if tau is None:
        def tau(i): return i
    W = []
    for a in Q.arrows():
        i, j = a[0], a[1]
        b = tau(a)
        na = 'a'+str(Q.arrow_num(a))
        nb = 'a'+str(Q.arrow_num(b))
        a = (a[0], a[1], na)
        b = (b[0], b[1], nb)
        a1 = (tau(j), i, na+'*')
        li = (i, tau(i), 'l'+str(i))
        lj = (j, tau(j), 'l'+str(j))
        W += [[a, lj, a1], [li, b, a1]]
    return Quiver(potential=W)
