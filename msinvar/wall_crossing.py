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
    Kronecker quiver: Quiver with 2 vertices and 2 arrows
    sage: Q.prec([2,2]) #fix precision vector
    sage: Q.total().dict()
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
    
    sage: I=Q.stacky([1,0]); I.dict()
    {(0, 0): 1,
     (0, 1): y/(-y^2 + 1),
     (0, 2): y^2/(y^6 - y^4 - y^2 + 1),
     (1, 0): y/(-y^2 + 1),
     (1, 1): (y^2 + 1)/(y^2 - 1),
     (1, 2): y/(-y^2 + 1),
     (2, 0): y^2/(y^6 - y^4 - y^2 + 1),
     (2, 1): y/(-y^2 + 1),
     (2, 2): (y^6 + y^4 + 2*y^2)/(y^6 - y^4 - y^2 + 1)}
    sage: I1=Q.stk2int(I); I1.dict()
    {(0, 1): 1, (1, 0): 1, (1, 1): (-y^2 - 1)/y, (1, 2): 1, (2, 1): 1}

::
    
    sage: I=Q.stacky([0,1]); I.dict()
    {(0, 0): 1,
     (0, 1): y/(-y^2 + 1),
     (0, 2): y^2/(y^6 - y^4 - y^2 + 1),
     (1, 0): y/(-y^2 + 1),
     (2, 0): y^2/(y^6 - y^4 - y^2 + 1)}
    sage: I1=Q.stk2int(I); I1.dict()
    {(0, 1): 1, (1, 0): 1}
    
::

    sage: Q.intAtt().dict() #integer attractor invariants
    {(0, 1): 1, (1, 0): 1}

::

    sage: Q.stable([1,0], slope=1/2).dict() #invariants of stable moduli spaces
    {(1, 1): y^2 + 1}
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from msinvar.rings import RF
from msinvar.tm_polynomials import TMPoly
from msinvar.utils import vec, cache
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

    - ``rank`` -- Rank of the lattice.
    - ``sform`` -- Skew-symmetric form on the lattice (a function).
    - ``prec`` -- precision vector for the quantum affine plane.
    """

    def __init__(self, rank=0, sform=None, prec=None):
        self._rank = rank
        self._sform = sform
        self.base = RF('y,q')
        self.y = self.base.gen(0)
        self.q = self.base.gen(1)
        self.gm = 1/self.y - self.y
        self.R = TMPoly(self.base, rank, prec=prec)
        self._total = None

    def __repr__(self):
        return f'Wall-crossing structure on a lattice of rank {self.rank()}'

    def total(self, I=None):
        """Total invariant, meaning the stacky invariant for the trivial stability.

        EXAMPLE::

            sage: from msinvar import *
            sage: Q=JordanQuiver(1); Q # Quivers inherit from WCS
            Jordan quiver: Quiver with 1 vertices and 1 arrows
            sage: Q.prec([2]); # set precision vector
            sage: I=Q.total(); I.poly()
            1 + (y^2/(y^2 - 1))*x + (y^6/(y^6 - y^4 - y^2 + 1))*x^2
            sage: I.Log().poly()
            (y^2/(y^2 - 1))*x
        """
        if I is not None:
            self._total = I
            return I
        if self._total is None:
            self._total = self.total_default()
        return self._total

    def rank(self): return self._rank
    def sform(self, a, b): return self._sform(a, b)
    def total_default(self): pass

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
        """Twist the product in the quantum affine plane using
        :meth:`twist`."""
        self.R.prod_twist = self.twist

    def untwist_product(self):
        """Untwist the product in the quantum affine plane 
        (make it commutative)."""
        self.R.prod_twist = None

    def prec(self, d=None):
        """Set or return precision vector for the quantum affine plane."""
        return self.R.prec(d)

    def stacky(self, z, I=None, algo='fast'):
        """
        Stacky invariant for the stability ``z``.

        See :meth:`total2stacky_algo1` for more details.
        The value of ``algo`` can be 'fast', 'fast2', 'slow'.
        If the total invariant ``I`` is None, we consider :meth:`total`.        

        EXAMPLES::

            sage: from msinvar import *
            sage: Q=KroneckerQuiver(2)
            sage: Q.prec([5,5])
            sage: z=Stability([1,0])
            sage: I1=Q.stacky(z,algo='fast')
            sage: I2=Q.stacky(z,algo='fast2')
            sage: I1([1,1])
            (y^2 + 1)/(y^2 - 1)
            sage: I1.poly()-I2.poly()
            0
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

    def stacky_from_total(self, *args, **kw):
        """Alias for :meth:`stacky`."""
        return self.stacky(*args, *kw)

    def Omh(self, *args, **kw):
        """Alias for :meth:`stacky`."""
        return self.stacky(*args, *kw)

    def rat_from_total(self, z, I=None, algo='fast'):
        """Return rational invariant for the stability ``z``."""
        Omh = self.stacky(z, I, algo)
        return self.stk2rat(Omh)  # assume that z is generic

    def Omb(self, *args, **kw):
        """Alias for :meth:`rat_from_total`."""
        return self.rat_from_total(*args, *kw)

    def int_from_total(self, z, I=None, algo='fast'):
        """Return integer invariant for the stability ``z``."""
        Omh = self.stacky(z, I, algo)
        return self.stk2int(Omh)  # assume that z is generic

    def Om(self, *args, **kw):
        """Alias for :meth:`int_from_total`."""
        return self.int_from_total(*args, *kw)

    def stacky2total(self, I, z):
        """See :meth:`stacky2total_algo1`."""
        return stacky2total_algo1(self, I, z)

    def stk2stk(self, I, z, z1, algo="fast"):
        """
        Transform stacky invariant ``I`` for stability ``z`` to the 
        stacky invariant for stability ``z1``.

        The value of ``algo`` can be 'fast' or 'slow'.
        """
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

    def rat2stk(self, I, z=None):
        """Transform rational invariant ``I`` to the stacky invariant.
        By default we assume that ``z`` is a generic stability.
        """
        I1 = Invariant(lambda d: I(d)/self.gm, self)
        if z is None:
            return I1.pexp()
        T = exp_transform(z)
        T.twist = self.twist
        return T(I1)

    def rat2int(self, I):
        """Transform rational invariant ``I`` to the integer invariant."""
        def I1(d): return I(d)/self.gm
        return Invariant(lambda d: IPsi_map(I1, d)*self.gm, self)

    def int2rat(self, I):
        """Transform integer invariant ``I`` to the rational invariant."""
        def I1(d): return I(d)/self.gm
        return Invariant(lambda d: Psi_map(I1, d)*self.gm, self)

    def stk2int(self, I):
        """Transform stacky invariant ``I`` to the integer invariant.
        We assume that the stability parameter is generic.
        """
        I1 = I.pLog()
        return Invariant(lambda d: I1(d)*self.gm, self)

    def int2stk(self, I):
        """Transform integer invariant ``I`` to the stacky invariant.
        We assume that the stability parameter is generic.
        """
        I1 = Invariant(lambda d: I(d)/self.gm, self)
        return I1.pExp()

    def poly(self, I, z=None, slope=None):
        """
        Transform invariant ``I`` to a polynomial,
        considering only degrees having a given ``slope`` with respect to a
        given stability parameter ``z``.
        """
        return I.poly(self, z,  slope)

    def series(self, *args, **kw):
        """Alias for :meth:`poly`."""
        return self.poly(*args, **kw)

    # quiver related methods (requires Euler form)
    def eform(self, a, b):
        """Euler form of a quiver or a curve."""
        pass

    def twist_T(self, d):
        """Auxiliary term twist that depends on :meth:`eform`."""
        return (-self.y)**(self.eform(d, d))

    def twist_TI(self, d):
        """Auxiliary term twist that depends on :meth:`eform`."""
        return (-self.y)**(-self.eform(d, d))
    # end of quiver related methods

    def stable_from_stacky(self, I):
        """
        Count stable objects of a fixed slope, assuming that the stacky
        invariant ``I`` counting semistable objects of that slope is given.
        Based on :arxiv:`0708.1259`.        
        """
        # f = self.series(I, z, slope)
        f = I.poly()
        if f.constant_coefficient() == 0:
            f = 1+f
        self.twist_product()  # quantum torus product
        f1 = f.invert()  # inversion in the quantum torus
        self.untwist_product()  # restore commutative product
        f2 = f1.term_twist(self.twist_TI)  # Exp(a/1-q)
        f3 = f2.Log()*(1-self.y**2)
        return Invariant(f3, self)

    def stable_from_total(self, z, slope=0, I=None):
        """
        Count ``z``-stable representations having a given ``slope``, assuming
        that the total invariant ``I`` is given.
        Based on :arxiv:`0708.1259`.        

        - ``z`` -- Stability parameter.
        - ``slope`` -- Slope value.
        - ``I`` -- Total invariant (if None, we consider :meth:`total`).
        """
        I1 = self.stacky(z, I).restrict(z, slope)
        return self.stable_from_stacky(I1)

    def stable(self, *args, **kw):
        """Alias for :meth:`stable_from_total`."""
        return self.stable_from_total(*args, **kw)

    def simple(self, I=None):
        """Count simple reprentations of a quiver, assuming that the total 
        invariant ``I`` is given. If ``I`` is None, we consider :meth:`total`.
        Based on :arxiv:`0708.1259`.        
        """
        if I is None:
            I = self.total()
        return self.stable_from_stacky(I)

    def self_stab(self, d):
        """Self-stability for the dimension vector ``d``."""
        n = len(d)
        a = [self.sform(vec.basis(i, n), d) for i in range(n)]
        return Stability(a)

    def attr_stab(self, d):
        """Attractor stability for the dimension vector ``d``."""
        z = self.self_stab(d)
        while not z.is_generic(d):
            z = z.randomize()
        return z
        # return self.self_stab(d).randomize()

    def stkAtt(self, I=None):
        """Calculate stacky attractor invariant assuming that the total 
        invariant is ``I``. If ``I`` is None, we consider :meth:`total`."""
        def f(d):
            z = self.attr_stab(d)
            I1 = self.stacky(z, I)
            return I1(d)
        return Invariant(f, self)

    def ratAtt(self, I=None):
        """Calculate rational attractor invariant assuming that the total 
        invariant is ``I``. If ``I`` is None, we consider :meth:`total`."""
        return self.stk2rat(self.stkAtt(I))

    def intAtt(self, I=None):
        """Calculate integer attractor invariant assuming that the total 
        invariant is ``I``. If ``I`` is None, we consider :meth:`total`."""
        return self.stk2int(self.stkAtt(I))

    OmhAtt = stkAtt
    OmbAtt = ratAtt
    OmAtt = intAtt

    def stkAtt2total(self, I=None):
        """
        Calculate total invariant, assuming that the stacky attractor invariant
        is ``I``.
        This is a recursive inversion of :meth:`stkAtt`.
        """
        if I is None:
            I = self.stkAtt_default()
        T = recursive_inversion(self.stkAtt)
        return T(I)

    def ratAtt2total(self, I=None):
        """
        Calculate total invariant, assuming that the rational attractor 
        invariant is ``I``.
        """
        if I is None:
            I = self.ratAtt_default()
        I1 = self.rat2stk(I)
        return self.stkAtt2total(I1)

    def intAtt2total(self, I=None):
        """
        Calculate total invariant, assuming that the integer attractor invariant
        is ``I``.
        """
        if I is None:
            I = self.intAtt_default()
        I1 = self.int2stk(I)
        return self.stkAtt2total(I1)

    def intAtt_default(self):
        def f(d):
            if sum(d) == 1:
                return 1
            return 0
        return Invariant(f, self)

    def ratAtt_default(self):
        return self.int2rat(self.intAtt_default())

    def stkAtt_default(self):
        return self.int2stk(self.intAtt_default())

    def attr_tree_formula(self, z, I, t=0):
        """See :meth:`msinvar.flow_trees.attr_tree_formula`."""
        from msinvar.flow_trees import attr_tree_formula
        return attr_tree_formula(self, z, I, t)

    def flow_tree_formula(self, z, I, **kw):
        """See :meth:`msinvar.flow_trees.flow_tree_formula`."""
        from msinvar.flow_trees import flow_tree_formula
        return flow_tree_formula(self, z, I, **kw)


WCS = WallCrossingStructure


def total2stacky_algo1(W, I, z):
    """
    Calculate stacky invariant for stability ``z``, assuming that the total 
    invariant is ``I``.

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant,
    - ``z``-- Stability.

    Based on :arxiv:`math/0204059` (5.5)
    and implementation by
    `Pieter Belmans <https://github.com/pbelmans/hodge-diamond-cutter>`_.
    Has comparable speed to :meth:`total2stacky_algo2`.
    """
    z = Stability.check(z)

    def stacky(d):
        if vec.iszero(d):
            return I(d)
        te = z.normalize(d)
        zero = [0]*len(d)
        l = [zero]+list(e for e in IntegerVectors_iterator(d)
                        if te(e) > 0)+[d]

        def T(i, j):
            a = vec.sub(l[j], l[i])
            if all(k >= 0 for k in a):
                return I(a) * W.twist(l[i], a)
            return 0

        n = len(l)
        b = [0]*(n-1)+[1]
        x = [0]*n
        for i in range(n-1, -1, -1):
            x[i] = b[i]-sum(T(i, j)*x[j] for j in range(i+1, n))
        return -x[0]

    return Invariant(stacky, W)


def total2stacky_algo2(W, I, z):
    """
    Calculate stacky invariant for stability ``z``, assuming that the total 
    invariant is ``I``.

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant,
    - ``z``-- Stability.

    Has comparable speed to :meth:`total2stacky_algo1`.
    """
    z = Stability.check(z)

    def stacky(d0):
        if vec.iszero(d0):
            return I(d0)
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

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant,
    - ``z``-- Stability.
    """
    z = Stability.check(z)

    @cache
    def total(d, slope=None):
        if vec.iszero(d):
            return I(d)
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

    This is a recursive inversion of :meth:`total2stacky_algo1`.
    It is slower than :meth:`stacky2total_algo1`.

    - ``W`` -- Wall-crossing structure,
    - ``I`` -- Invariant,
    - ``z``-- Stability.
    """
    def T(I): return total2stacky_algo1(W, I, z)
    T1 = recursive_inversion(T)
    return T1(I)


# def total2stkAtt(W, I):
#     """
#     Calculate stacky attractor invariant, assuming that the total invariant
#     is ``I``.

#     INPUT:

#     - ``W`` -- Wall-crossing structure,
#     - ``I`` -- Invariant.
#     """
#     def f(d):
#         z = W.attr_stab(d)
#         return total2stacky_algo1(W, z, I)(d)
#     return Invariant(f, W)


# def total2intAtt(W, I):
#     """
#     Calculate integer attractor invariant, assuming that the total invariant
#     is ``I``.

#     INPUT:

#     - ``W`` -- Wall-crossing structure,
#     - ``I`` -- Invariant.
#     """
#     I = total2stkAtt(W, I).pLog()
#     return Invariant(lambda d: I(d)*W.gm, W)


# def stacky2total_algo2(W, I, z):
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
#     z = Stability.check(z)

#     @cache
#     def total(d0):
#         if vec.iszero(d0):
#             return 1
#         c = z(d0)

#         @cache
#         def tail(d):
#             s = 0
#             for e in IntegerVectors_iterator(d):
#                 d1 = vec.sub(d, e)
#                 if vec.iszero(d1):
#                     s += total(e)
#                 elif z(d1) < c:
#                     s -= tail(d1)*total(e)*W.twist(e, d1)
#             return s
#         s = 0
#         for e in IntegerVectors_iterator(d0):
#             d1 = vec.sub(d0, e)
#             if vec.iszero(d1):
#                 pass
#             elif z(d1) < c:
#                 s -= tail(d1)*total(e)*W.twist(e, d1)
#         return I(d0)-s
#     return Invariant(total, W)


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
