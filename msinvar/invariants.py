r"""
Invariant class that associates a value with every dimension vector

An invariant can be constructed from a TMPolynomial, a dictionary, a function 
or another invariant. Every invariant can be transformed to a TMPolynomial
or a dictionary. We cache the values of an invariant so that they don't have
to be computed repeatedly.

EXAMPLES::

    sage: from msinvar import *
    sage: R=TMPoly(prec=[3]); R
    Multivariate Polynomial Ring in x over Rational Field truncated at degree [3]
    sage: x=R.gen(); x.Exp()
    1 + x + x^2 + x^3
    sage: I=Invariant(x.Exp()); I
    Invariant
    sage: I.dict()
    {(0,): 1, (1,): 1, (2,): 1, (3,): 1}
    sage: I.Log().dict()
    {(1,): 1}
    sage: I.Log().poly()
    x
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from numpy import array
from sage.misc.misc_c import prod
from sage.arith.misc import divisors, gcd, moebius, factor
from sage.functions.other import factorial
from sage.rings.rational_field import QQ
from sage.calculus.functional import expand

from msinvar.tm_polynomials import TMPolynomial, TMPolynomialRing
from msinvar.utils import cache, vec
from msinvar.iterators import IntegerVectors_iterator, ListPartitions_iterator
from msinvar.iterators import OrderedMultiPartitions_iterator, OrderedMultiPartitionsLE_iterator
from msinvar.stability import Stability


class Invariant:
    """
    Invariant class that associates an invariants I(d) with every dimension 
    vector d.

    An argument can be a TMPolynomial, a dictionary, a function or another 
    invariant.
    Objects have methods that return associated TMPolynomials and dictionaries.
    If possible, we save the ring (TMPoly) in order to transform our invariant
    to a polynomial.
    
    - ``f`` -- function, TMPolynomial, dictionary, Invariant.
    - ``R`` -- a ring TMPolynomialRing.
    """

    def __init__(self, f, R=None):
        self.f = f
        if R is not None:
            self.R = self.get_ring(R)
        elif isinstance(f, TMPolynomial):
            self.R = f.parent()
        elif isinstance(f, Invariant):
            self.R = f.R
        else:
            self.R = None

    def __call__(self, d):
        return self.value(d)

    def __repr__(self):
        return 'Invariant'

    @cache
    def value(self, d):
        """The value of the invariant at ``d``."""
        f = self.f
        if isinstance(f, TMPolynomial):
            return f.coeff(d)
        if isinstance(f, dict):
            try:
                return f[tuple(d)]
            except:
                return 0
        return f(d)

    def get_ring(self, I=None):
        """The ring (TMPolynomialRing) saved by the invariant.""" 
        if I is not None:
            if hasattr(I, 'R'):
                return I.R
            elif isinstance(I, TMPolynomialRing):
                return I
            return None
        if self.R is not None:
            return self.R
        raise ValueError("Need a ring")

    def get_prec(self, d=None):
        """Precision vector of the ring."""
        from msinvar.wall_crossing import WCS
        if d is None:
            return self.get_ring().prec()
        if isinstance(d, WCS) or isinstance(d, TMPolynomialRing):
            return d.prec()
        return d

    def dict(self, prec=None, stab=None, slope=0):
        """
        Convert invariant to a dictionary, with keys bounded by self.prec().
        """
        dct = {}
        prec = self.get_prec(prec)

        def add(d):
            c = self.value(d)
            if c != 0:
                dct[tuple(d)] = c

        add([0]*len(prec))
        if stab is not None:
            stab = Stability(stab)
        for d in IntegerVectors_iterator(prec):
            if stab is None or stab(d) == slope:
                add(d)
        return dct

    def series(self, R=None, stab=None, slope=None):
        """
        Convert invariant to a TMPolynomial.
        """

        R = self.get_ring(R)
        prec = R.prec()
        return R(self.dict(prec, stab, slope))

    to_series = series
    poly = series

    def term_twist(self, f):
        """Twist every term using the function ``f``."""
        return Invariant(lambda d: self(d)*f(d), self)

    def restrict(self, z, slope=0):
        """Restrict invariant to dimension vectors d such that z(d)=slope,
        where ``z`` is a Stability (or the corresponding vector)."""
        
        z = Stability.check(z)

        def f(d):
            if vec.iszero(d) or z(d) == slope:
                return self.value(d)
            return 0
        return Invariant(f, self)

    def subs(self, **kw):
        """Apply substitution to every term."""
        def f(d):
            try:
                return self.value(d).subs(**kw)
            except:
                return self.value(d)
        return Invariant(f, self)

    def root_vars(self, k=2):
        """Substitute y^k by y for every variable y."""
        def f(d):
            try:
                return self.value(d).root_vars(k)
            except:
                return self.value(d)
        return Invariant(f, self)

    def simp(self, d=None):
        """Simplify every term."""
        if d is None:
            return Invariant(lambda d: simp(self(d)), self)
        return simp(self(d))

    def plog(self, algo="fast"):
        """Take log along every ray."""
        if algo == "fast":
            return plog_fast(self)
        return plog_map(self)

    def pexp(self, algo="fast"):
        """Take exp along every ray."""
        if algo == "fast":
            return pexp_fast(self)
        return pexp_map(self)

    def exp(self):
        """The usual exponential. This method can be quite slow."""
        return Invariant(self.series().exp())

    def log(self):
        """The usual logarithm. This method can be quite slow."""
        return Invariant(self.series().log())

    def Psi(self):
        """Plethystic map Psi."""
        return Psi_map(self)

    def IPsi(self):
        """Plethystic map Psi inverse."""
        return IPsi_map(self)

    def pLog(self):
        """Take plethystic Log along every ray."""
        return self.plog().IPsi()

    def pExp(self):
        """Take plethystic Exp along every ray."""
        return self.Psi().pexp()

    def Exp(self):
        """Take plethystic Exp."""
        return Invariant(self.series().Exp())

    def Log(self):
        """Take plethystic Log."""
        return Invariant(self.series().Log())


# Useful maps between invariants


def log_map(I, d=None):
    """Usual log."""
    if d is None:
        return Invariant(lambda d: log_map(I, d), I)

    def F(l): return QQ((-1)**(len(l)-1)/len(l))
    val = 0
    for l in OrderedMultiPartitions_iterator(d):
        val += F(l) * prod(I(a) for a in l)
    return val


def exp_map(I, d=None):
    """Usual exp."""
    if d is None:
        return Invariant(lambda d: exp_map(I, d), I)

    def F(l): return QQ(1/factorial(len(l)))
    val = 0
    for l in OrderedMultiPartitions_iterator(d):
        val += F(l) * prod(I(a) for a in l)
    return val


def plog_map(I, d=None):
    """log along every ray."""
    # Note that log_transform() gives the same map, but plog_map is faster
    if d is None:
        return Invariant(lambda d: plog_map(I, d), I)
    if vec.iszero(d):
        return 0
    d = array(d)
    n = gcd(d)
    d0 = d//n
    I1 = Invariant(lambda k: I(k[0]*d0))
    return log_map(I1, d=[n])


def pexp_map(I, d=None):
    """exp along every ray."""
    # Note that exp_transform() gives the same map, but pexp_map is faster
    if d is None:
        return Invariant(lambda d: pexp_map(I, d), I)
    if vec.iszero(d):
        return 1
    d = array(d)
    n = gcd(d)
    d0 = d//n
    I1 = Invariant(lambda k: I(k[0]*d0))
    return exp_map(I1, d=[n])


def Psi_map(I, d=None):
    """Plethystic Psi."""
    if d is None:
        return Invariant(lambda d: Psi_map(I, d), I)
    if vec.iszero(d):
        return 0
    n = gcd(d)
    val = 0
    for k in divisors(n):
        d1 = [i//k for i in d]
        if I(d1) != 0:
            val += I(d1).adams(k)/k
    return val


def IPsi_map(I, d=None):
    """Plethsytic Psi inverse."""
    if d is None:
        return Invariant(lambda d: IPsi_map(I, d), I)
    if vec.iszero(d):
        return 0
    n = gcd(d)
    val = 0
    for k in divisors(n):
        d1 = [i//k for i in d]
        if I(d1) != 0:
            val += moebius(k) * I(d1).adams(k)/k
    return val


def recursive_inversion(T):
    """Invert the transformation (map on invariants) T, assuming that
    T(I)(d) is a sum of I(d) and some expression that depends just on
    I(e) with e<d.

    Note that exp and log are not of this form for d=0.
    """
    def T1(I):
        """The inverse of T"""
        def f(d):
            """The required invariant J=T1(I) such that I=T(J)"""
            def f1(e):  # the fake J
                if vec.equal(d, e):
                    return 0
                return J(e)
            J1 = Invariant(f1, I)
            return I(d)-T(J1)(d)
        J = Invariant(f, I)
        return J
    return T1


class Transform:
    """
    A class to transform invariants. Based on :arxiv:`2101.07636`.

    It is encoded by a map from the set of lists of vectors
    to rational numbers.

    We define actions of transforms on 1-collections (Invariant class),
    plethysm between transforms, inverse transfroms.

    INPUT:
    
    - ``F`` -- transformation map from lists of vectors to Q (or base ring).
    - ``twist`` -- product twist, map from sequences of vectors to the base ring.
    """

    def __init__(self, F, twist=None):
        self.F = F
        if twist is None:
            self.twist = lambda l: 1
        else:
            self.twist = twist

    def __call__(self, I, d=None):
        if isinstance(I, Invariant):
            return self.action(I, d)
        if isinstance(I, Transform):
            return self.plethysm(I, d)
        return self.value(I)

    @cache
    def value(self, l):
        return self.F(l)

    def action(self, I, d):
        """Action of the transformation on an invariant ``I``."""
        if d is None:
            return Invariant(lambda d: self.action(I, d), I)
        val = 0
        for l in OrderedMultiPartitions_iterator(d):
            val += self.F(l) * self.twist(l) * prod(I(a) for a in l)
        return val

    def plethysm(self, T, l=None):
        """
        Plethysm of self and T.
        """
        if l is None:
            return Transform(lambda l: self.plethysm(T, l))
        val = 0
        for part in ListPartitions_iterator(l):
            l1 = [vec.add(p) for p in part]
            val += self(l1) * prod(T(p) for p in part)
        return val

    @cache
    def inverse(self, l=None):
        """Inverse transformation."""
        if l is None:
            return Transform(lambda l: self.inverse(l))
        if len(l) == 1:
            return 1
        # if not isinstance(l, np.ndarray):
        #     l = array(l)
        val = 0
        for part in ListPartitions_iterator(l):
            if len(part) == len(l):
                continue
            l1 = [vec.add(p) for p in part]
            val += self.inverse(l1) * prod(self(p) for p in part)
        return -val

    def dict(self, d):
        """Dictionary of values for all lists of vectors with the sum <=d."""
        def tpl(l):
            return tuple(tuple(x) for x in l)
        dct = {}
        for l in OrderedMultiPartitionsLE_iterator(d):
            c = self.value(l)
            if c != 0:
                dct[tpl(l)] = c
        return dct


# Useful transforms
def proportional(a, b=None):
    """Check if vectors a,b are proportional."""
    if b is not None:
        a, b = array(a), array(b)
        return all(a*sum(b) == b*sum(a))
    return all(proportional(a[0], a[i]) for i in range(1, len(a)))


def exp_transform(z=None):
    if z is None:
        def F(l):
            return QQ(1/factorial(len(l))) if proportional(l) else 0
    else:
        z = Stability.check(z)

        def F(l):
            return QQ(1/factorial(len(l))) if z.equal(l) else 0
    return Transform(F)


def log_transform(z=None):
    if z is None:
        def F(l):
            return QQ((-1)**(len(l)-1)/len(l)) if proportional(l) else 0
    else:
        z = Stability.check(z)

        def F(l):
            return QQ((-1)**(len(l)-1)/len(l)) if z.equal(l) else 0
    return Transform(F)


def reineke_sign(l, z, tot=None):
    """Used in :meth:`reineke_transform`."""
    if tot is None:
        tot = vec.add(l)
    te = z.normalize(tot)
    c = 0
    for i in range(len(l)-1):
        c += te(l[i])
        if c <= 0:
            return 0
    return (-1)**(len(l)-1)


def reineke_transform(z):
    """Based on :arxiv:`math/0204059` and :arxiv:`2101.07636`."""
    z = Stability.check(z)
    return Transform(lambda l: reineke_sign(l, z))


def joyce_sign(l, z, z1, tot=None):
    """Used in :meth:`joyce_transform`."""
    if tot is None:
        tot = vec.add(l)
    v = [0]*len(tot)
    s = 1
    for i in range(len(l)-1):
        v = vec.add(v, l[i])
        c = z(l[i], l[i+1])
        c1 = z1(v, tot)
        if c <= 0 and c1 > 0:
            s = -s
        elif c > 0 and c1 <= 0:
            pass
        else:
            return 0
    return s


def joyce_transform(z, z1):
    """Based on :arxiv:`math/0410268` and :arxiv:`2101.07636`."""
    z = Stability.check(z)
    z1 = Stability.check(z1)
    return Transform(lambda l: joyce_sign(l, z, z1))


def HN_transform(z):
    """See :arxiv:`2101.07636`."""
    z = Stability.check(z)
    z1 = Stability.trivial(z.dim())
    return joyce_transform(z, z1)


def simp(f):
    """Simplify an expression."""
    if f == 0:
        return f
    return expand(factor(f))


def indivisible(d):
    """Return True if ``d`` is an indivisible vector."""
    d = array(d)
    n = gcd(d)
    return d//n


def plog_fast(I):
    """
    Return the logarithm taken separately along each ray.

    We construct a dictionary of such logarithms with keys parametrized
    by indivisible vectors. The logarithm is taken in the ring of truncated
    polynomials. Therefore I.R and its precision vector are required.
    """
    dct = {}

    def f(d):
        if vec.iszero(d):
            return 0
        n = gcd(d)
        d0 = tuple(i//n for i in d)
        if d0 in dct:
            return dct[d0].coeff(d)
        prec = I.get_prec()
        m = min(prec[i]//d0[i] for i in range(len(d0)) if d0[i] != 0)
        dct1 = {}
        for k in range(0, m+1):
            d1 = tuple(k*i for i in d0)
            c = I(d1)
            if c != 0:
                dct1[d1] = c
        f1 = I.R(dct1).log()
        dct[d0] = f1
        return f1.coeff(d)
    return Invariant(f, I)


def pexp_fast(I):
    """See :meth:`plog_fast`."""
    dct = {}

    def f(d):
        if vec.iszero(d):
            return 1
        n = gcd(d)
        d0 = tuple(i//n for i in d)
        if d0 in dct:
            return dct[d0].coeff(d)
        prec = I.get_prec()
        m = min(prec[i]//d0[i] for i in range(len(d0)) if d0[i] != 0)
        dct1 = {}
        for k in range(1, m+1):
            d1 = tuple(k*i for i in d0)
            c = I(d1)
            if c != 0:
                dct1[d1] = c
        f1 = I.R(dct1).exp()
        dct[d0] = f1
        return f1.coeff(d)
    return Invariant(f, I)
