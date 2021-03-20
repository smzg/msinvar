r"""
This module contains algorithms related to wall-crossing.

We define the class :class:`WCS` associated to a lattice with a skew-symmetric
form (or to a quiver). This class contains the quantum affine plane where all
invariants live and where wall-crossing takes place. It also contains
algorithms to jump between various types of invariants:
    
    1. Stacky, Rational, Integer invariants associated to a stability parameter.
    2. Total invariant (stacky invariant associated to the trivial stability parameter).
    3. Stacky, Rational, Integer invariants associated to the attractor stability.
    4. Invariant counting stable representations of a quiver for a given stability parameter.

EXAMPLES::

    sage: from msinvar import *
    sage: Q=KroneckerQuiver(2); Q
    Quiver with 2 vertices and 2 arrows
    sage: W=Q.wcs([2,2]); W #fix precision vector
    Wall-crossing structure on a lattice of rank 2
    sage: W.total().dict()
    {(0, 0): 1,
    (0, 1): y/(-y^2 + 1),
    (0, 2): y^2/(y^6 - y^4 - y^2 + 1),
    (1, 0): y/(-y^2 + 1),
    (1, 1): y^4/(y^4 - 2*y^2 + 1),
    (1, 2): y^7/(-y^8 + 2*y^6 - 2*y^2 + 1),
    (2, 0): y^2/(y^6 - y^4 - y^2 + 1),
    (2, 1): y^7/(-y^8 + 2*y^6 - 2*y^2 + 1),
    (2, 2): y^12/(y^12 - 2*y^10 - y^8 + 4*y^6 - y^4 - 2*y^2 + 1)}

::
    
    sage: I=W.stacky([1,0]); I.dict()
    {(0, 0): 1,
     (0, 1): y/(-y^2 + 1),
     (0, 2): y^2/(y^6 - y^4 - y^2 + 1),
     (1, 0): y/(-y^2 + 1),
     (1, 1): (y^2 + 1)/(y^2 - 1),
     (1, 2): y/(-y^2 + 1),
     (2, 0): y^2/(y^6 - y^4 - y^2 + 1),
     (2, 1): y/(-y^2 + 1),
     (2, 2): (y^6 + y^4 + 2*y^2)/(y^6 - y^4 - y^2 + 1)}
    sage: I1=W.stk2int(I); I1.dict()
    {(0, 1): 1, (1, 0): 1, (1, 1): (-y^2 - 1)/y, (1, 2): 1, (2, 1): 1}

::
    
    sage: I=W.stacky([0,1]); I.dict()
    {(0, 0): 1,
     (0, 1): y/(-y^2 + 1),
     (0, 2): y^2/(y^6 - y^4 - y^2 + 1),
     (1, 0): y/(-y^2 + 1),
     (2, 0): y^2/(y^6 - y^4 - y^2 + 1)}
    sage: I1=W.stk2int(I); I1.dict()
    {(0, 1): 1, (1, 0): 1}
    
::

    sage: W.intAtt().dict() #integer attractor invariants
    {(0, 1): 1, (1, 0): 1}

::

    sage: W.stable([1,0], slope=1/2).dict() #invariants of stable moduli spaces
    {(1, 1): y^2 + 1}
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.matrix.constructor import matrix
from msinvar.rings import RF
from msinvar.tm_polynomials import TMPoly
from msinvar.utils import vec, phi, cache
from msinvar.iterators import IntegerVectors_iterator
from msinvar.stability import Stability
from msinvar.invariants import Invariant, Psi_map, IPsi_map
from msinvar.invariants import log_transform,  exp_transform, reineke_transform, joyce_transform
from msinvar.invariants import recursive_inversion


class WallCrossingStructure:
    """
    Wall-crossing structure.

    Contains the quantum affine algebra and methods to compute and transform
    different types of invariants (stacky, rational, integer DT invariants 
    corresponding to different stability parameters).

    INPUT:

    - ``quiver`` -- Quiver (if None, ``rank`` and ``sform`` should be present).
    - ``rank`` -- Rank of the lattice.
    - ``sform`` -- Skew-symmetric form on the lattice (a function).
    - ``prec`` -- precision vector for the quantum affine plane.
    """

    def __init__(self, quiver=None, rank=0, sform=None, prec=None):
        self.base = RF('y,q')
        self.y = self.base.gen(0)
        self.q = self.base.gen(1)
        self.gm = 1/self.y - self.y
        if quiver is not None:
            self.rank = quiver.vertex_num()
            self.sform = quiver.sform
            self.quiver = quiver
        else:
            self.rank = rank
            self.sform = sform
        self.R = TMPoly(self.base, self.rank, prec=prec)
        self._total = None

    def total(self, I=None):
        """Total invariant, meaning the stacky invariant for the trivial stability.
        
        EXAMPLE::
            
            sage: from msinvar import *
            sage: Q=JordanQuiver(1); Q
            Quiver with 1 vertices and 1 arrows
            sage: W=Q.wcs([2]); W
            Wall-crossing structure on a lattice of rank 1
            sage: I=W.total(); I.poly()
            1 + (y^2/(y^2 - 1))*x + (y^6/(y^6 - y^4 - y^2 + 1))*x^2
            sage: I.Log().poly()
            (y^2/(y^2 - 1))*x
        """
        if I is not None:
            self._total = I
            return
        if self._total is None:
            self._total = Invariant(self.total_default, self)
        return self._total

    def __repr__(self):
        return "Wall-crossing structure on a lattice of rank "+str(self.rank)

    def twist(self, a, b=None):
        """
        Product twist for a pair of vectors or a list of vectors.
        """
        if b is not None:
            return (-self.y)**self.sform(a, b)
        if len(a) <= 1:
            return 1
        v = vec.add(a)
        s = 0
        for i in range(len(a)-1):
            v = vec.sub(v, a[i])
            s += self.sform(a[i], v)
        return (-self.y)**s

    def twist_product(self):
        """Twist the product in the quantum affine plane (self.R) using
        self.twist."""
        self.R.prod_twist = self.twist

    def untwist_product(self):
        """Untwist the product in the quantum affine plane 
        (make it commutative)."""
        self.R.prod_twist = None

    def prec(self, d=None):
        """Set or return precision vector for the quantum affine plane."""
        return self.R.prec(d)

    def stacky_from_total(self, z, I=None, algo="fast"):
        """See :meth:`total2stacky_algo1`.
        The value of ``algo`` can be 'fast', 'fast2', 'slow'.
        """
        if I is None:
            I = self.total()
        if algo == "fast":
            return total2stacky_algo1(self, I, z)
        if algo == "fast2":
            return total2stacky_algo2(self, I, z)
        if algo == "slow":
            T = reineke_transform(z)
            T.twist = self.twist
            return T(I)
        raise ValueError("No such algo")

    stacky = stacky_from_total
    Omh = stacky_from_total

    def rat_from_total(self, z, I=None):
        """Return rational invariant for the stability ``z``."""
        Omh = self.stacky(z, I)
        return self.stk2rat(Omh)  # assume that z is generic

    def int_from_total(self, z, I=None):
        """Return integer invariant for the stability ``z``."""
        Omh = self.stacky(z, I)
        return self.stk2int(Omh)  # assume that z is generic

    Omb = rat_from_total
    Om = int_from_total

    def total_from_stacky(self, I, z):
        """See :meth:`stacky2total_algo1`."""
        return stacky2total_algo1(self, I, z)

    def stk2stk(self, I, z, z1, algo="fast"):
        """Calculate stacky invariant for stability ``z1`` assuming that
        stacky invariant for stability ``z`` is ``I``."""
        if algo == "fast":
            I1 = stacky2total_algo1(self, I, z)
            return total2stacky_algo1(self, I1, z1)
        if algo == "slow":
            T = joyce_transform(z, z1)
            T.twist = self.twist
            return T(I)

    def stk2rat(self, I, z=None):
        """Transform stacky invariant ``I`` to the rational invariant.
        By default we assume that ``z`` is a generic stability.
        """
        if z is None:
            I1 = I.plog()
            return Invariant(lambda d: I1(d)*self.gm, self)
        T = log_transform(z)
        T.twist = lambda l: self.twist(l)*self.gm
        return T(I)

    def rat2int(self, I):
        """Transform rational invariant ``I`` to the integer invariant."""
        def I1(d): return I(d)/self.gm
        return Invariant(lambda d: IPsi_map(I1, d)*self.gm, self)

    def stk2int(self, I):
        """Transform stacky invariant ``I`` to the integer invariant."""
        I1 = I.pLog()
        return Invariant(lambda d: I1(d)*self.gm, self)

    def int2rat(self, I):
        def I1(d): return I(d)/self.gm
        return Invariant(lambda d: Psi_map(I1, d)*self.gm, self)

    def rat2stk(self, I, z=None):
        I1 = Invariant(lambda d: I(d)/self.gm, self)
        if z is None:
            return I1.pexp()
        T = exp_transform(z)
        T.twist = self.twist
        return T(I1)

    def int2stk(self, I):
        I1 = I.term_twist(lambda d: 1/self.gm)
        return I1.pExp()

    def series(self, I, stab=None, slope=None):
        return I.poly(self, stab,  slope)

    poly=series
        
    # quiver related methods (requires Euler form)
    def eform(self, a, b):
        return self.quiver.eform(a, b)

    def qform(self, d):
        return self.eform(d, d)

    def total_default(self, d):
        """
        Invariants of a quiver corresponding to the trivial stability
        """
        y = self.y
        return (-y)**(-self.qform(d))/phi(1/y**2, d)

    def twist_T(self, d):
        return (-self.y)**(self.qform(d))

    def twist_TI(self, d):
        return (-self.y)**(-self.qform(d))
    # end of quiver related methods

    def stable_from_stacky(self, I, z, slope=0):
        """
        Count z-stable objects having given slope, assuming that stacky count
        of z-semistable objects is known.

        INPUT:

        I : Invariant counting z-semistable objects.

        z : Stability parameter.

        slope : The slope.
        """
        f = self.series(I, z, slope)
        if f.constant_coefficient() == 0:
            f = 1+f
        self.twist_product()  # quantum torus product
        f1 = f.invert()  # inversion in the quantum torus
        self.untwist_product()  # restore commutative product
        f2 = f1.term_twist(self.twist_TI)  # Exp(a/1-q)
        f3 = f2.Log()*(1-self.y**2)
        return Invariant(f3, self)

    def stable_from_total(self, z, slope=0, I=None):
        I = self.stacky(z, I)
        return self.stable_from_stacky(I, z, slope)

    stable = stable_from_total

    def simple(self, I=None):
        z = [0]*self.rank
        return self.stable_from_total(z, 0, I)

    def self_stab(self, d):
        n = len(d)
        a = [self.sform(vec.basis(i, n), d) for i in range(n)]
        return Stability(a)

    def attr_stab(self, d):
        z = self.self_stab(d)
        while not z.is_generic(d):
            z = z.randomize()
        return z
        # return self.self_stab(d).randomize()

    def stkAtt(self, I=None):
        def f(d):
            z = self.attr_stab(d)
            return self.stacky(z, I)(d)
        return Invariant(f, self)

    def intAtt(self, I=None):
        I = self.stkAtt(I).pLog()
        return Invariant(lambda d: I(d)*self.gm, self)

    def stkAtt2total(self, I):
        return stkAtt2total(self, I)


WCS = WallCrossingStructure


def total2stacky_algo1(W, I, z):
    """
    Calculate stacky invariant for stability ``z``, assuming that the total 
    invariant is ``I``.

    INPUT:

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant,
    - ``z``-- Stability.

    Has the same speed as :meth:`total2stacky_algo2`.
    """
    z = Stability.check(z)

    def matrix_T(d):
        zero = [0]*len(d)
        l = [zero]+list(e for e in IntegerVectors_iterator(d)
                        if z(e) > z(d))+[d]
        n = len(l)
        T = matrix(W.base, n)
        for i in range(n):
            for j in range(i, n):
                a = vec.sub(l[j], l[i])
                if all(x >= 0 for x in a):
                    T[i, j] = I(a) * W.twist(l[i], a)
        return T

    def solve(A, b):
        # A is an upper-triangular nxn matrix (with 1 on the diagonal) and b is an n-vector
        n = len(b)
        x = [0]*n
        for i in range(n-1, -1, -1):
            x[i] = b[i]-sum(A[i, j]*x[j] for j in range(i+1, n))
        return x

    def stacky(d):
        T = matrix_T(d)
        n = T.nrows()
        v = [0]*(n-1)+[1]
        return -solve(T, v)[0]

    return Invariant(stacky, W)


def total2stacky_algo2(W, I, z):
    """
    Calculate stacky invariant for stability ``z``, assuming that the total 
    invariant is ``I``.

    INPUT:

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant,
    - ``z``-- Stability.

    Has the same speed as :meth:`total2stacky_algo1`.
    """
    z = Stability.check(z)

    def stacky(d0):
        if vec.iszero(d0):
            return 1
        te = z.normalize(d0)

        @cache
        def tail(d):
            s = 0
            for e in IntegerVectors_iterator(d):
                d1 = vec.sub(d, e)
                if vec.iszero(d1):
                    s += I(e)
                elif te(d1) < 0:
                    s -= tail(d1)*I(e)*W.twist(e, d1)
            return s
        return tail(d0)
    return Invariant(stacky, W)


def stacky2total_algo1(W, I, z):
    """
    Calculate total invariant, assuming that the stacky invariant for 
    stability ``z`` is ``I``.

    This algorithm is faster than :meth:`stacky2total_algo2`.

    INPUT:

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant,
    - ``z``-- Stability.
    """
    z = Stability.check(z)

    @cache
    def total(d, slope=None):
        if vec.iszero(d):
            return 1
        s = 0
        for e in IntegerVectors_iterator(d):
            if slope is None or z(e) < slope:
                d1 = vec.sub(d, e)
                s += total(d1, z(e))*I(e)*W.twist(e, d1)
        return s
    return Invariant(total, W)


def stacky2total_algo2(W, I, z):
    """
    Calculate total invariant, assuming that the stacky invariant for 
    stability ``z`` is ``I``.

    This is a recursive inversion of total2stacky. It is slower than 
    :meth:`stacky2total_algo1`.

    INPUT:

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant,
    - ``z``-- Stability.
    """
    z = Stability.check(z)

    @cache
    def total(d0):
        if vec.iszero(d0):
            return 1
        te = z.normalize(d0)

        @cache
        def tail(d):
            s = 0
            for e in IntegerVectors_iterator(d):
                d1 = vec.sub(d, e)
                if vec.iszero(d1):
                    s += total(e)
                elif te(d1) < 0:
                    s -= tail(d1)*total(e)*W.twist(e, d1)
            return s
        s = 0
        for e in IntegerVectors_iterator(d0):
            d1 = vec.sub(d0, e)
            if vec.iszero(d1):
                pass
            elif te(d1) < 0:
                s -= tail(d1)*total(e)*W.twist(e, d1)
        return I(d0)-s
    return Invariant(total, W)


def total2stkAtt(W, I):
    """
    Calculate stacky attractor invariant, assuming that the total invariant
    is ``I``.

    INPUT:

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant.
    """
    def f(d):
        z = W.attr_stab(d)
        return total2stacky_algo1(W, z, I)(d)
    return Invariant(f, W)


def total2intAtt(W, I):
    """
    Calculate integer attractor invariant, assuming that the total invariant
    is ``I``.

    INPUT:

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant.
    """
    I = total2stkAtt(W, I).pLog()
    return Invariant(lambda d: I(d)*W.gm, W)


def stkAtt2total(W, I):
    """
    Calculate total invariant, assuming that the stacky attractor invariant
    is ``I``.

    This is a recursive inversion of :meth:`total2stkAtt`.

    INPUT:

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant.
    """
    def T(I): return total2stkAtt(W, I)
    T1 = recursive_inversion(T)
    return T1(I)


def intAtt2total(W, I):
    """
    Calculate total invariant, assuming that the integer attractor invariant
    is ``I``.

    INPUT:

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant.
    """
    I1 = W.int2stk(I)
    return stkAtt2total(W, I1)

# def stkAtt2total_algo(W, I):
#     """
#     Calculate total invariant, assuming that the stacky invariant for
#     stability ``z`` is ``I``.

#     This is a recursive inversion of total2stacky. It is slower than
#     :meth:`stacky2total_algo1`.

#     INPUT:

#     - ``W`` -- Wall-crossing structure,
#     - ``I`` -- Invariant,
#     - ``z``-- Stability.
#     """

#     @cache
#     def total(d0):
#         if vec.iszero(d0):
#             return 1
#         z = W.attr_stab(d0)
#         te = z.normalize(d0)

#         @cache
#         def tail(d):
#             s = 0
#             for e in IntegerVectors_iterator(d):
#                 d1 = vec.sub(d, e)
#                 if vec.iszero(d1):
#                     s += total(e)
#                 elif te(d1) < 0:
#                     s -= tail(d1)*total(e)*W.twist(e, d1)
#             return s
#         s = 0
#         for e in IntegerVectors_iterator(d0):
#             d1 = vec.sub(d0, e)
#             if vec.iszero(d1):
#                 pass
#             elif te(d1) < 0:
#                 s -= tail(d1)*total(e)*W.twist(e, d1)
#         return I(d0)-s
#     return Invariant(total, W)
