r"""
Collection of iterators.

EXAMPLES::

    sage: from msinvar.iterators import *
    sage: list(IntegerVectors_iterator([2,2]))
    [[1, 0], [2, 0], [0, 1], [1, 1], [2, 1], [0, 2], [1, 2], [2, 2]]

::
    
    sage: M=[[1,2,1],[3,1,1]]
    sage: list(Multiplicities_iterator(M,[3,4]))
    [[1, 0, 0],
     [0, 1, 0],
     [1, 1, 0],
     [0, 0, 1],
     [1, 0, 1],
     [0, 1, 1],
     [0, 0, 2],
     [0, 0, 3]]

"""

# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import numpy as np
from msinvar.utils import vec


def IntegerVectors_iterator(vect):
    r"""
    Iterator over integer vectors 0 < a <= vect

    It is more efficient than
    :meth:`sage.combinat.vector_partition.IntegerVectorsIterator`.

    ``vect`` -- a list of integers.

    EXAMPLES::

        sage: from msinvar.iterators import *
        sage: list(IntegerVectors_iterator([2,2]))
        [[1, 0], [2, 0], [0, 1], [1, 1], [2, 1], [0, 2], [1, 2], [2, 2]]
    """
    n = len(vect)
    a = list(0 for i in range(n))
    k = 0
    while k < n:
        if a[k] < vect[k]:
            a[k] += 1
            for i in range(k):
                a[i] = 0
            yield a.copy()
            k = -1
        k += 1


def Multiplicities_iterator(M, b):
    r"""
    Iterator over integer vectors a>0 such that M*a<=b,
    where ``M`` is a matrix and ``b`` is a vector.

    - ``M`` -- Matrix of size m x n
    - ``b`` -- Vector of size m

    """
    M = np.array(M)
    b = np.array(b)
    m, n = M.shape
    a = np.zeros(n, int)
    k = 0
    while k < n:
        a[k] += 1
        if all(M[:, k:].dot(a[k:]) <= b):
            a[:k] = 0
            yield list(a)
            k = -1
        else:
            a[k] -= 1
        k += 1


def OrderedMultiPartitionsLE_iterator(vect):
    r"""
    Iterator over collections of vectors
    (a_1,..,a_k) such that a_1+..+a_k <= ``vect`` and a_i>0.
    """
    for b in IntegerVectors_iterator(vect):
        yield [b]
        for a in OrderedMultiPartitionsLE_iterator(vec.sub(vect, b)):
            yield [b]+a


def OrderedMultiPartitions_iterator(vect):
    r"""
    Iterator over collections of vectors
    (a_1,...,a_k) such that a_1+...+a_k = ``vect`` and a_i>0.
    """
    if not vec.zero(vect):
        yield [vect]
    for b in IntegerVectors_iterator(vect):
        for a in OrderedMultiPartitions_iterator(vec.sub(vect, b)):
            yield [b]+a


def OrderedPartitionsLE_iterator(n):
    r"""
    Iterator over collections of positive numbers
    (a_1,..,a_k) such that a_1+..+a_k <= ``n``.
    """
    for b in range(1, n+1):
        yield [b]
        for a in OrderedPartitionsLE_iterator(n-b):
            yield [b]+a


def OrderedPartitions_iterator(n):
    r"""
    Iterator over collections of positive numbers
    (a_1,...,a_k) such that a_1+...+a_k = ``n``.
    """
    if n != 0:
        yield [n]
    for b in range(1, n+1):
        for a in OrderedPartitions_iterator(n-b):
            yield [b]+a


def MultiPartitionsLE_iterator(vect, bound=None):
    r"""
    Iterator over collections of vectors
    ``bound``>=a_1 >=...>=a_k>0 such that a_1+...+a_k <= ``vect``.
    """
    if bound is None:
        bound = vect.copy()
    else:
        bound = vec.vmin(vect, bound)
    for b in IntegerVectors_iterator(bound):
        yield [b]
        for a in MultiPartitionsLE_iterator(vec.sub(vect, b), b):
            yield [b]+a


def MultiPartitions_iterator(vect, bound=None):
    r"""
    Iterator over collections of vectors
    ``bound``>=a_1>=...>=a_k>0 such that a_1+...+a_k = ``vect``.
    """
    if bound is None:
        bound = vect
    if not vec.zero(vect) and vec.le(vect, bound):
        yield [vect]
    bound = vec.vmin(vect, bound)
    for b in IntegerVectors_iterator(bound):
        for a in MultiPartitions_iterator(vec.sub(vect, b), b):
            yield [b]+a


def Subsets_iterator(n):
    r"""
    Iterator over non-empty subsets of the set {0,..,n-1}.
    """
    for i in range(n):  # the maximal element of the subset
        yield [i]
        for s in Subsets_iterator(i):
            yield s+[i]


def ListPartitions_iterator(l):
    r"""
    Iterator over collections of nonempty disjoint lists (l1,..,lk)
    such that l1+...+lk= ``l``.
    """
    for a in OrderedPartitions_iterator(len(l)):
        yield list(l[sum(a[:i]):sum(a[:i+1])] for i in range(len(a)))


def UnorderedMultiPartitions_iterator(vect, bound=None):
    r"""
    Iterator over sets of vectors
    (a_1,...,a_k) such that a_1+...+a_k = ``vect`` and a_i>0.
    We order vectors lexicographically.
    """
    def lex(a, b=None):
        if b is None:
            return True
        for i in range(len(a)):
            if a[i] < b[i]:
                return True
            if a[i] > b[i]:
                return False
        return True

    if not vec.zero(vect) and lex(vect, bound):
        yield [vect]
    for b in IntegerVectors_iterator(vect):
        if lex(b, bound):
            for a in UnorderedMultiPartitions_iterator(vec.sub(vect, b), b):
                yield [b]+a


# def UnorderedPartitions_iterator(d):
#     l=list(IntegerVectors_iterator(d))
#     M=np.array(l).T


# def mult2list(m, l):
#     """
#     Return a list that contains entry l[i] with multiplicity m[i]

#     ``m``, ``l`` are lists of the same length and ``m`` contains integers.
#     """
#     L = []
#     for i in range(len(l)):
#         L += [l[i]]*m[i]
#     return L


# def Multiplicities_iterator1(M,b):
#     """
#     Iterator over integer vectors a>0 such that M*a<=b,
#     where M is an integer matrix and b is an integer vector.

#     Parameters
#     ----------
#     M : Integer matrix of size m x n
#     b : Integer vector of size m

#     Yields
#     -------
#     Integer vector a of size n such that M*a<=b

#     """
#     import numpy as np
#     M=matrix(M)
#     b=matrix(b)
#     n=M.ncols()
#     a=matrix(n,1)
#     k=0
#     while k<n:
#         a[k,0]+=1
#         if M[:,k:]*a[k:]<=b:
#             a[:k]=0
#             yield a
#             k=-1
#         else:
#             a[k,0]-=1
#         k+=1
