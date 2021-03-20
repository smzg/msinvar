r"""
Stability function

EXAMPLES::

    sage: from msinvar import Stability
    sage: z=Stability([1,0]); z
    Stability function [[1, 0], [1, 1]]
    sage: z([1,2]) #slope
    0.3333333333333333
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import numpy as np
from sage.arith.misc import gcd
from msinvar.iterators import IntegerVectors_iterator


class Stability:
    """Stability function corresponding to the vectors ``a``, ``b``.
    
    EXAMPLES::
    
        sage: from msinvar import Stability
        sage: z=Stability([1,0]); z
        Stability function [[1, 0], [1, 1]]
        sage: z([1,2]) #slope
        0.3333333333333333
        sage: z([1,0],[0,1]) #difference of slopes
        1.0
        sage: z.less([1,0],[0,1])
        False
    """

    def __init__(self, a, b=None):
        if isinstance(a, Stability):
            self.a = a.a
            self.b = a.b
            return
        self.a = np.array(a)
        if b is None:
            self.b = np.ones(len(a), int)
        else:
            self.b = np.array(b)

    def __repr__(self):
        return "Stability function "+str([list(self.a), list(self.b)])

    def __call__(self, d, e=None):
        if e is None:
            return self.slope(d)
        return self.compare(d, e)

    def slope(self, d):
        """Slope of the vector ``d``."""
        return np.dot(self.a, d)/np.dot(self.b, d)

    def compare(self, d, e):
        """Difference of slopes."""
        c = self.slope(d)-self.slope(e)
        return self._rnd(c)

    def less(self, d, e):
        """Return True if slope(d)<slope(e)."""
        return self.slope(d) < self.slope(e)

    def lesseq(self, d, e):
        """Return True if slope(d)<=slope(e)."""
        return self.slope(d) <= self.slope(e)

    def equal(self, d, e=None):
        """Return True if slope(d)=slope(e)."""
        if e is not None:
            return self.slope(d) == self.slope(e)
        if len(d) == 1:
            return True
        te = self.normalize(d[0])
        return all(te(d[i]) == 0 for i in range(1, len(d)))

    def dim(self):  # used in HN_transform
        """Rank of the lattice."""
        return len(self.a)

    def normalize(self, d):
        """Return function f such that f(e)<0 iff slope(e)<slope(d)."""
        c = self.a-self.slope(d)*self.b
        return lambda d: self._rnd(np.dot(c, d))

    def randomize(self):
        """Generic perturbation of self."""
        a, b = self.a, self.b
        return Stability(a+np.random.rand(len(a))*1e-5, b)

    def is_generic(self, prec):
        """Return True if slope(d)=slope(e) for d,e<=prec implies that 
        d, e are proportional."""
        s = set()
        for d in IntegerVectors_iterator(prec):
            if gcd(d) == 1:
                c = self.slope(d)
                if c in s:
                    return False
                s.add(c)
        return True

    @staticmethod
    def trivial(n):
        """Return trivial stability of dimension ``n``."""
        return Stability(np.zeros(n, int))

    @staticmethod
    def check(z):
        """Check if ``z`` is a stability.
        If it is a vector, convert it to a stability."""
        if isinstance(z, Stability):
            return z
        return Stability(z)

    def _rnd(self, c):
        if abs(c) < 1e-8:
            return 0
        return c
