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
    """

    def __init__(self, quiver=None, rank=0, sform=None, prec=None):
        """
        Initialize a wall-crossing structure.

        INPUT:

        quiver : object of type Quiver (if none, rank and sform should be present).

        n : Rank of the lattice.

        sform : Skew-symmetric form on the lattice. A function.

        prec : Degree bound.
        """
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
        if quiver is not None:
            self.total = Invariant(self.total_default, self)

    def __repr__(self):
        return "Wall-crossing structure on a lattice of rank "+str(self.rank)

    def twist(self, a, b=None):
        """
        Product twist (for stacky invariants) for a pair of vectors 
        or a list of vectors
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
        self.R.prod_twist = self.twist

    def untwist_product(self):
        self.R.prod_twist = None

    def prec(self, d=None):
        return self.R.prec(d)

    def stacky_from_total(self, z, algo="fast"):
        if algo == "fast":
            return total2stacky_algo1(self, self.total, z)
        if algo == "fast2":
            return total2stacky_algo2(self, self.total, z)
        if algo == "slow":
            T = reineke_transform(z)
            T.twist = self.twist
            return T(self.total)
        raise ValueError("No such algo")

    def rat_from_total(self, z):
        Omh = total2stacky_algo1(self, self.total, z)
        return self.stk2rat(Omh)  # assume that z is generic

    def int_from_total(self, z):
        Omh = total2stacky_algo1(self, self.total, z)
        return self.stk2int(Omh)  # assume that z is generic

    stacky = stacky_from_total
    Omh = stacky_from_total
    Omb = rat_from_total
    Om = int_from_total

    def total_from_stacky(self, I, z):
        return stacky2total_algo1(self, I, z)

    def stk2stk(self, I, z, z1, algo="fast"):
        if algo == "fast":
            I1 = stacky2total_algo1(self, I, z)
            return total2stacky_algo1(self, I1, z1)
        if algo == "slow":
            T = joyce_transform(z, z1)
            T.twist = self.twist
            return T(I)

    def stk2rat(self, I, z=None):
        if z is None:
            I1 = I.plog()
            return Invariant(lambda d: I1(d)*self.gm, self)
        T = log_transform(z)
        T.twist = lambda l: self.twist(l)*self.gm
        return T(I)

    def rat2int(self, I):
        def I1(d): return I(d)/self.gm
        return Invariant(lambda d: IPsi_map(I1, d)*self.gm, self)

    def stk2int(self, I):
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

    def to_series(self, I, stab=None, slope=None):
        return I.to_series(self, stab,  slope)

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

        Parameters
        ----------
        I : Invariant counting z-semistable objects.

        z : Stability parameter.

        slope : The slope.
        """
        f = self.to_series(I, z, slope)
        if f.constant_coefficient() == 0:
            f = 1+f
        self.twist_product()  # quantum torus product
        f1 = f.invert()  # inversion in the quantum torus
        self.untwist_product()  # restore commutative product
        f2 = f1.term_twist(self.twist_TI)  # Exp(a/1-q)
        f3 = f2.Log()*(1-self.y**2)
        return Invariant(f3, self)

    def stable_from_total(self, z, slope=0):
        I = self.stacky_from_total(z)
        return self.stable_from_stacky(I, z, slope)

    stable = stable_from_total

    def simple(self):
        z = [0]*self.rank
        return self.stable_from_total(z, 0)

    def self_stab(self, d):
        n = len(d)
        a = [self.sform(vec.basis(i, n), d) for i in range(n)]
        return Stability(a)

    def attr_stab(self, d):
        z=self.self_stab(d)
        while not z.is_generic(d):
            z=z.randomize()
        return z
        # return self.self_stab(d).randomize()

    def stkAtt(self):
        def f(d):
            z = self.attr_stab(d)
            return self.stacky(z)(d)
        return Invariant(f, self)

    def intAtt(self):
        I = self.stkAtt().pLog()
        return Invariant(lambda d: I(d)*self.gm, self)

    def stkAtt2total(self, I):
        return stkAtt2total_algo(self, I)


WCS = WallCrossingStructure


def total2stacky_algo1(W, I, z):
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
    Calculate stacky invariants for stability z, assuming that total invariants
    are known

    Parameters
    ----------
    W : Wall-crossing structure

    z : Stability
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
    """Recursive inversion of total2stacky. It is slower than algo1."""
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


def stkAtt2total_algo(W, I):
    @cache
    def total(d0):
        if vec.iszero(d0):
            return 1
        z = W.attr_stab(d0)
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
    def f(d):
        z = W.attr_stab(d)
        return total2stacky_algo1(W, I, z)(d)
    return Invariant(f, W)


def total2intAtt(W, I):
    I = total2stkAtt(W, I).pLog()
    return Invariant(lambda d: I(d)*W.gm, W)


def stkAtt2total(W, I):
    def T(I): return total2stkAtt(W, I)
    T1 = recursive_inversion(T)
    return T1(I)

def intAtt2total(W,I):
    I1=W.int2stk(I)
    return stkAtt2total(W, I1)
