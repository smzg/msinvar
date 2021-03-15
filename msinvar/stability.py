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

import numpy as np
from numpy.random import rand
from sage.arith.misc import gcd

from msinvar.iterators import IntegerVectors_iterator


class Stability:
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
        return np.dot(self.a, d)/np.dot(self.b, d)

    def compare(self, d, e):
        c = self.slope(d)-self.slope(e)
        return rnd(c)

    def less(self, d, e):
        return self.slope(d) < self.slope(e)

    def lesseq(self, d, e):
        return self.slope(d) <= self.slope(e)

    def equal(self, d, e=None):
        if e is not None:
            return self.slope(d) == self.slope(e)
        if len(d) == 1:
            return True
        te = self.normalize(d[0])
        return all(te(d[i]) == 0 for i in range(1, len(d)))

    def dim(self):  # used in HN_transform
        return len(self.a)

    def normalize(self, d):
        c = self.a-self.slope(d)*self.b
        return lambda d: rnd(np.dot(c, d))

    def randomize(self):
        a, b = self.a, self.b
        return Stability(a+rand(len(a))*1e-5, b)

    def is_generic(self, prec):
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
        return Stability(np.zeros(n, int))

    @staticmethod
    def check(z):
        if isinstance(z, Stability):
            return z
        return Stability(z)


def rnd(c):
    if abs(c) < 1e-8:
        return 0
    return c
