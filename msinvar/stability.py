r"""
Stability function

EXAMPLES::

    sage: from msinvar import Stability
    sage: z=Stability([1,0]); z
    Stability function [[1, 0], [1, 1]]
    sage: z([1,2]) #slope
    0.3333333333
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
from msinvar.utils import vec
from sage.misc.prandom import random


class Stability:
    """Stability function corresponding to the vectors ``a``, ``b``.

    EXAMPLES::

        sage: from msinvar import Stability
        sage: z=Stability([1,0]); z
        Stability function [[1, 0], [1, 1]]
        sage: z([1,2]) #slope
        0.3333333333
        sage: z([1,0],[0,1]) #difference of slopes
        1
        sage: z.less([1,0],[0,1])
        False
    """

    def __init__(self, a, b=None):
        if isinstance(a, Stability):
            self.a = a.a
            self.b = a.b
            return
        self.a = a  # np.array(a)
        if b is None:
            self.b = [1]*len(a)  # np.ones(len(a), int)
        else:
            self.b = b  # np.array(b)

    def __repr__(self):
        return "Stability function "+str([list(self.a), list(self.b)])

    def __call__(self, d, e=None):
        if e is None:
            return self.slope(d)
        return self.compare(d, e)

    def slope(self, d):
        """Slope of the vector ``d``."""
        return np.round(vec.dot(self.a, d)/vec.dot(self.b, d), 10)

    def compare(self, d, e):
        """Difference of slopes."""
        return self.slope(d)-self.slope(e)

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

    def weight(self, d):
        """Return vector w such that w*e<0 iff slope(e)<slope(d)."""
        return vec.sub(self.a, vec.scal(self(d), self.b))

    def normalize(self, d):
        """Return function f such that f(e)<0 iff slope(e)<slope(d)."""
        w = self.weight(d)
        return lambda d: np.round(vec.dot(w, d), 10)

    def randomize(self):
        """Generic perturbation of self."""
        a, b=self.a, self.b
        return Stability([i+random()*1e-5 for i in a], b)

    def is_generic(self, prec):
        """Return True if slope(d)=slope(e) for d,e<=prec implies that
        d, e are proportional."""
        s=set()
        for d in IntegerVectors_iterator(prec):
            if gcd(d) == 1:
                c=self.slope(d)
                if c in s:
                    return False
                s.add(c)
        return True

    @ staticmethod
    def trivial(n):
        """Return trivial stability of dimension ``n``."""
        return Stability([0]*n)

    @ staticmethod
    def check(z):
        """Check if ``z`` is a stability.
        If it is a vector, convert it to a stability."""
        if isinstance(z, Stability):
            return z
        return Stability(z)
