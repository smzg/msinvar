r"""
Collection of utilities.
"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from functools import reduce, cmp_to_key
import inspect
from sage.misc.dev_tools import import_statements
from sage.misc.misc_c import prod


def isiterable(i):
    return hasattr(i, '__iter__')


def hashable(x):
    try:
        hash(x)
        return x
    except:
        return tuple(hashable(i) for i in x)


def memoize(f):
    memo = {}

    def helper(*args):
        x = hashable(args)
        if x not in memo:
            memo[x] = f(*args)
        return memo[x]
    return helper


cache = memoize


def sort(l, le, reverse=False):
    """Sort the list l in ascending order according to the function le"""
    def cmp(i, j):
        return -1 if le(i, j) else 1

    return sorted(l, key=cmp_to_key(cmp), reverse=reverse)


sort1 = sort


def prod1(a, b=None):
    if b is not None:
        return a*b
    return reduce(prod1, a, 1)


class vec:
    @staticmethod
    def le(a, b):
        return all(i <= j for i, j in zip(a, b))

    @staticmethod
    def iszero(a):
        return all(i == 0 for i in a)

    zero = iszero

    @staticmethod
    def equal(a, b):
        return all(i == j for i, j in zip(a, b))

    @staticmethod
    def vmin(a, b):
        return [min(i, j) for i, j in zip(a, b)]

    @staticmethod
    def vmax(a, b):
        return [max(i, j) for i, j in zip(a, b)]

    @staticmethod
    def scal(c, v):
        return [c * i for i in v]

    @staticmethod
    def dot(a, b):
        return sum(i * j for i, j in zip(a, b))

    @staticmethod
    def sub(a, b):
        return [i - j for i, j in zip(a, b)]

    @staticmethod
    def add(a, b=None):
        """
        Return the sum of vectors a, b.
        If b=None, return the sum of vectors in the list a.
        """
        if b is not None:
            return [i + j for i, j in zip(a, b)]
        return reduce(vec.add, a)

    @staticmethod
    def basis(i, n):
        """Return the standard ``i``-th basis vector of dimension ``n``.
        Here 0<=i<n."""
        return [0] * i + [1] + [0] * (n - i - 1)
    # def zero(n): return [0]*n


##### Information commands ##########


def info(f):
    if inspect.isclass(f):
        return inspect.getmro(f)
    else:
        print(inspect.getsource(f))


def disp(G):
    from msinvar.quivers import Quiver
    if isinstance(G, Quiver):
        G = G.digraph
    return G.plot(edge_labels=True,
                  layout='circular').matplotlib(figsize=[10, 6])


which = import_statements


def timer(f, n=1000):
    from time import time
    t1 = time()
    for i in range(n):
        f
    return time()-t1


def set_plots():
    from sage.graphs.graph_plot import DEFAULT_PLOT_OPTIONS as DPO
    DPO['layout'] = 'circular'
    DPO['loop_size'] = .3
#####################################


# @cache
def phi(q, n):
    r"""Return `(q)_n=\prod_{i=1}^n(1-q^i)`.

    If ``n`` is a list, apply :meth:`phi` to all entries and take the product.
    """
    if isiterable(n):
        return prod(phi(q, i) for i in n)
    return prod((1-q**i) for i in range(1, n+1))


# @cache
def poch(a, q, n):
    r"""Return the Pochhammer symbol `(a;q)_n=\prod_{i=0}^{n-1}(1-aq^i)`.

    If ``n`` is a list, apply :meth:`poch` to all entries and take the product.
        """
    if isiterable(n):
        return prod(poch(a, q, i) for i in n)
    return prod((1-a*q**i) for i in range(n))
