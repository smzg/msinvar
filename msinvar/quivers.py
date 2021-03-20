r"""
Quivers and quivers with potentials

Define a Quiver class (with or without potental) as well as scpecial
subclasses ChainQuiver, CyclicQuiver and constructors KroneckerQuiver, 
JordanQuiver and translation_PQ (for translation potential quivers of quivers
with automorphisms, see :arxiv:`1911.01788`).

EXAMPLES::

    sage: from msinvar.quivers import *
    sage: Q=KroneckerQuiver(3); Q
    Quiver with 2 vertices and 3 arrows
    sage: Q.vertices()
    [1, 2]
    sage: Q.arrows()
    [(1, 2, 0), (1, 2, 1), (1, 2, 2)]
    
::
    
    sage: Q=Quiver('1-2-3'); Q
    Quiver with 3 vertices and 2 arrows
    sage: Q.arrows()
    [(1, 2, 1), (2, 3, 1)]
    
::
    
    sage: Q=Quiver('1-2,1-3'); Q
    Quiver with 3 vertices and 2 arrows
    sage: Q.arrows()
    [(1, 2, 1), (1, 3, 1)]
    
::
    
    sage: Q=CyclicQuiver(3); Q
    Quiver with 3 vertices and 3 arrows
    sage: Q.arrows()
    [(0, 1, 1), (1, 2, 1), (2, 0, 1)]
    sage: Q2=Q.translation_PQ(1); Q2
    Quiver with 3 vertices, 9 arrows and potential with 6 terms
    sage: Q2.arrows()
    [(0, 1, 'a0'),
     (0, 1, 'a1*'),
     (0, 1, 'l0'),
     (1, 2, 'a1'),
     (1, 2, 'a2*'),
     (1, 2, 'l1'),
     (2, 0, 'a0*'),
     (2, 0, 'a2'),
     (2, 0, 'l2')]
    sage: Q2.potential()
    [[(0, 1, 'a0'), (1, 2, 'l1'), (2, 0, 'a0*')],
     [(0, 1, 'l0'), (1, 2, 'a1'), (2, 0, 'a0*')],
     [(1, 2, 'a1'), (2, 0, 'l2'), (0, 1, 'a1*')],
     [(1, 2, 'l1'), (2, 0, 'a2'), (0, 1, 'a1*')],
     [(2, 0, 'a2'), (0, 1, 'l0'), (1, 2, 'a2*')],
     [(2, 0, 'l2'), (0, 1, 'a0'), (1, 2, 'a2*')]]
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import re
from sage.rings.integer import Integer
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from msinvar.utils import vec
from msinvar.wall_crossing import WCS


class Quiver(DiGraph):
    """
    Create a quiver from a list of arrows or from a potential.

    EXAMPLES::

        sage: from msinvar.quivers import *
        sage: Q=Quiver('1-2-3,3-4'); Q
        Quiver with 4 vertices and 3 arrows
        sage: Q.arrows()
        [(1, 2, 1), (2, 3, 1), (3, 4, 1)]
        
    ::
        
        sage: Q=Quiver(potential='a[1,2]a[2,3]a[3,1]'); Q.arrows()
        [(1, 2, 1), (2, 3, 1), (3, 1, 1)]
        sage: Q.potential()
        [[(1, 2, 1), (2, 3, 1), (3, 1, 1)]]
        
    ::
        
        sage: Q=Quiver(potential=['1-2-3-1',[[1,2,2],[2,1]]]); Q.arrows()
        [(1, 2, 1), (1, 2, 2), (2, 1, 1), (2, 3, 1), (3, 1, 1)]
        sage: Q.potential()
        [[(1, 2, 1), (2, 3, 1), (3, 1, 1)], [(1, 2, 2), (2, 1, 1)]]
    """

    def __init__(self, data=None, potential=None, prec=None,
                 pos=None, loops=None, format=None, weighted=False, data_structure=None, name=None):
        """Init a quiver."""
        if data is not None and isinstance(data, str):
            data = Quiver._get_arrow_set(data)
        if potential is not None:
            self._potential=Quiver._get_potential(potential)
            if data is None:
                data = Quiver._get_arrow_set(self._potential)
        else:
            self._potential = None
        super().__init__(data=data, multiedges=True, loops=True)
        self._arrows = []
        vert = self.vertices()
        self.vertex_dict = {vert[i]: i for i in range(self.vertex_num())}
        self._set_quiver_dict()
        self._wcs = WCS(self, prec=prec)

    def wcs(self, prec=None):
        """Wall-crossing structure of the quiver.
        
        EXAMPLES::
            
            sage: from msinvar.quivers import *
            sage: Q=KroneckerQuiver(2); Q
            Quiver with 2 vertices and 2 arrows
            sage: W=Q.wcs([3,3]); W
            Wall-crossing structure on a lattice of rank 2
            sage: W.intAtt().dict() # integer attractor invariants
            {(0, 1): 1, (1, 0): 1}
        """
        if prec is not None:
            self._wcs.prec(prec)
        return self._wcs

    def prec(self, prec=None):
        """Set (or return) precision vector for the WCS of the quiver."""
        return self._wcs.prec(prec)

    def arrows(self, n=None):
        """Return the ``n``-th arrow or the list of arrows if ``n`` is None."""
        if not n is None:
            return self._arrows[n]
        if len(self._arrows) != self.arrow_num():
            self._arrows = list(self.edges())
        return self._arrows

    def potential(self): return self._potential

    def arrow_num(self, a=None):
        """Return the index of an arrow ``a`` or the number of all arrows 
        if ``a`` is None.
        """
        if a is None:
            return self.size()
        return self.arrows().index(a)

    def vertex_num(self, i=None):
        """Return the index of a vertex ``i`` or the number of all vertices
        if ``i`` is None.
        """
        if i is None:
            return self.order()
        return self.vertex_dict[i]

    def s(self, a):
        """Source of an arrow ``a``, which can be an arrow or an index in the
        list self.arrows().
        
        EXAMPLES::

            sage: from msinvar import *
            sage: Q=Quiver('1-2-3'); Q.arrows()
            [(1, 2, 1), (2, 3, 1)]
            sage: Q.s([1,2,1])
            1
            sage: Q.s(0)
            1
        """
        if isinstance(a, int) or isinstance(a, Integer):
            a = self.arrows(a)
        return a[0]

    def t(self, a):
        """Target of an arrow ``a``, which can be an arrow or an index in the
        list self.arrows().

        EXAMPLES::

            sage: from msinvar import *
            sage: Q=Quiver('1-2-3'); Q.arrows()
            [(1, 2, 1), (2, 3, 1)]
            sage: Q.t([1,2,1])
            2
            sage: Q.t(0)
            2
        """
        if isinstance(a, int) or isinstance(a, Integer):
            a = self.arrows(a)
        return a[1]

    def _set_quiver_dict(self):
        d = {}
        for a in self.arrows():
            i, j = self.vertex_num(a[0]), self.vertex_num(a[1])
            if (i, j) in d:
                d[(i, j)] += 1
            else:
                d[(i, j)] = 1
        self.quiver_dict = d

    def eform(self, a, b):
        """Euler form of the quiver.
        
        EXAMPLES::
            
            sage: from msinvar import *
            sage: Q=KroneckerQuiver(2); Q.arrows()
            [(1, 2, 0), (1, 2, 1)]
            sage: Q.eform([1,0],[0,1])
            -2
            sage: Q.eform([0,1],[1,0])
            0
        """
        s = sum(a[i]*b[i] for i in range(len(a)))
        for (i, j), k in self.quiver_dict.items():
            s -= k * a[i]*b[j]
        return s

    def sform(self, a, b):
        """Skew-symmetrized Euler form of the quiver.
        
        EXAMPLES::
            
            sage: from msinvar import *
            sage: Q=KroneckerQuiver(2); Q.arrows()
            [(1, 2, 0), (1, 2, 1)]
            sage: Q.sform([1,0],[0,1])
            -2
            sage: Q.sform([0,1],[1,0])
            2
        """
        s = 0
        for (i, j), k in self.quiver_dict.items():
            s += k * (a[j]*b[i]-a[i]*b[j])
        return s

    def get_eform_matrix(self):
        """Matrix of the eform."""
        n = self.vertex_num()
        M = matrix.identity(n)
        for (i, j), k in self.quiver_dict.items():
            M[(i, j)] -= k
        return M

    def get_sform_matrix(self):
        """Matrix of the sform."""
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
        """
        See :meth:`translation_PQ`
        """
        return translation_PQ(self, tau)

    def Ginzburg_PQ(self):
        """
        Return the Ginzburg quiver (with potential) for the original quiver.

        This construction corresponds to :meth:`translation_PQ` for the
        trivial automorphism.
        """
        return self.translation_PQ()

    @staticmethod
    def _get_potential(s):
        """
        Convert a string or list potential to a list of cycles and drop the 
        coefficients. 
        
        Potential is of the form:

        - F[1,2,1]*F[2,3]F[3,1]...-F[1,2,2]F[2,1]...+...
        - ['1-2-1','2-3-4-2']
        """
        if isinstance(s, str):
            s = re.split('\\+|\\-', s)
        return [Quiver._get_path(p) for p in s]
    
    @staticmethod
    def _get_arrow_set(s):
        """
        Return the set of arrows that appear in the string or list of paths.
        
        The string can be of the form:
            - '1-2-3,3->4'
            
        The list can be of the form:        
            - [[[1,2],[2,3]],[[3,4]]]
            - ['1-2-3','3-4']
            - ['a[1,2]a[2,3]','a[3,4]']
        """
        if isinstance(s, str):
            s = s.split(',')
        arrs = set()
        for p in s:
            arrs.update(Quiver._get_path(p))
        return list(arrs)


    @staticmethod
    def _get_path(p):
        """
        Convert a path string or list to a list of arrows.

        A path string can be of the form:
            
            - 'a[1,2,1]*a[2,1]a[1,2,2]b[2,3]' (sequence of arrows)
            - '1-2-3' (sequence of vertices)
        
        A path list can be of the form:
            - [[1,2],[2,1],[1,2,2]]
            - ['1,2','2,1','1,2,2']
        
        Only the content inside square brackets is used in the first case.
        """
        if isinstance(p, str):
            # get a list of strings between [,]; use .*? for a non greeedy search
            l = re.findall(".*?\\[(.*?)\\]", p)
            if len(l)!=0:
                return [Quiver._get_arrow(a) for a in l]
            l = [int(re.sub(r'\D', '', i)) for i in p.split('-')] # remove non-digits
            return [(l[i], l[i+1], 1) for i in range(len(l)-1)]                
        return list(Quiver._get_arrow(a) for a in p)

    @staticmethod
    def _get_arrow(a):
        """
        Convert a string or a list to a tuple (arrow).
        
        A string can be of the form:
            - '1,2'
            - '1,2,1' (the last integer is the number of the arrow)
        """
        if isinstance(a, str):
            a = tuple(int(i) for i in a.split(','))
        else:
            a = tuple(a)
        return (*a, 1) if len(a) == 2 else a


def KroneckerQuiver(m=2):
    """Return the Kronecker quiver with m arrows."""
    l = [[1, 2, i] for i in range(m)]
    return Quiver(l)


def JordanQuiver(m=1):
    """Return the quiver with 1 vertex and m loops."""
    l = [[1, 1, i] for i in range(m)]
    return Quiver(l)


class ChainQuiver(Quiver):
    """
    Chain quiver with n vertices.

    PARAMETERS:

    - ``n`` -- number of vertices.
    - ``prec`` -- truncation vector for the quantum affine plane.

    EXAMPLES::

        sage: from msinvar import *
        sage: Q=ChainQuiver(3); Q.arrows()
        [(0, 1, 1), (1, 2, 1)]
        sage: list(Q.indecs())
        [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
        sage: Q.dim([0,2])
        [1, 1, 1]
        sage: Q.dim([2,2])
        [0, 0, 1]
        sage: Q.hom([2,2],[0,2])
        1
    """
    def __init__(self, n, prec=None):
        l = [[i, i+1, 1] for i in range(n-1)]
        super().__init__(l, prec=prec)

    def dim(self, a):
        """
        Return the dimension vector of an indecomposable representation
        parametrized by a pair a=(i,j) with `i<=j`.
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
    """
    Cyclic quiver with n vertices.

    PARAMETERS:

    - ``n`` -- number of vertices.
    - ``prec`` -- truncation vector for the quantum affine plane.

    EXAMPLES::

        sage: from msinvar import *
        sage: Q=CyclicQuiver(3); Q.arrows()
        [(0, 1, 1), (1, 2, 1), (2, 0, 1)]
    """
    def __init__(self, n, prec=None):
        l = [[i, i+1, 1] for i in range(n-1)]
        super().__init__(l, prec=prec)

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
        return _cong_number(max(0, d2-d1), d2, i-k, r)

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
        Cyclic shift for indecomposable modules: (SM)_i=M_{i+k}.

        EXAMPLES:

            sage: from msinvar import *
            sage: Q=CyclicQuiver(3)
            sage: Q.shift_indec((0,1,0),1)
            (2, 0, 0)
        """
        r = self.vertex_num()
        return ((a[0]-k) % r, (a[1]-k) % r, a[2])

    def translation(self, k=1):
        """
        Return a quiver automorphism given by the k-cyclic shift.

        EXAMPLES::

            sage: from msinvar import *
            sage: Q=CyclicQuiver(3)
            sage: f=Q.translation(2)
            sage: {i:f(i) for i in Q.vertices()}
            {0: 2, 1: 0, 2: 1}
            sage: {i:f(i) for i in Q.arrows()}
            {(0, 1, 1): (2, 0, 1), (1, 2, 1): (0, 1, 1), (2, 0, 1): (1, 2, 1)}
        """
        r = self.vertex_num()

        def tau(a):
            if isinstance(a, tuple):  # an arrow
                return ((a[0]+k) % r, (a[1]+k) % r, 1)
            return (a+k) % r
        return tau

    def translation_PQ(self, k=1):
        """
        Return translation potential quiver associated with the cyclic quiver
        and its k-cyclic shift.

        EXAMPLES::

            sage: from msinvar import *
            sage: Q=CyclicQuiver(3)
            sage: Q2=Q.translation_PQ(1); Q2
            Quiver with 3 vertices, 9 arrows and potential with 6 terms
        """
        tau = self.translation(k)
        return translation_PQ(self, tau)

def _cong_number(a, b, i, r):
    """
    Return the number of a<=s<=b with s=i mod r.
    """
    if a > b:
        return 0
    q = (b-a) // r
    if (b-a) % r < (i-a) % r:
        return q
    return q+1


def translation_PQ(Q, tau=None):
    """
    Construct a new quiver with potential based on a quiver ``Q`` and its
    automorphism ``tau``. See :arxiv:`1911.01788`.
    
    If ``tau`` is None, we consider the trivial automorphism.

    1. For any arrow a:i->j we add an arrow tau(j)->i.
    2. For any vertex i we add an arrow i->tau(i).
    3. We consider potential consisting of 3-cycles (for every arrow i->j):
    
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
