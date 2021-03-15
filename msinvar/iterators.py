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
from msinvar.utils import sub_vec, zero_vec, min_vec, le_vec


def IntegerVectors_iterator(vect):
    """
    Iterator over integer vectors 0 < a <= vect

    It is more efficient than
    **sage.combinat.vector_partition.IntegerVectorsIterator**

    Parameters
    ----------
    vect : A list of integers

    Examples::
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
    """
    Iterator over integer vectors a>0 such that M*a<=b,
    where M is a matrix and b is a vector.

    Parameters
    ----------
    M : Matrix of size m x n
    b : Vector of size m

    """
    M = np.array(M)
    b = np.array(b)
    m, n = M.shape
    a = np.zeros(n, int)
    k = 0
    while k < n:
        a[k] += 1
        if all(M[:, k:].dot(a[k:])<= b):
            a[:k] = 0
            yield list(a)
            k = -1
        else:
            a[k] -= 1
        k += 1


def OrderedMultiPartitionsLE_iterator(vect):
    """
    Iterator of collections of vectors
    (a_1,..,a_k) such that a_1+..+a_k<=vect and a_i>0
    """
    for b in IntegerVectors_iterator(vect):
        yield [b]
        for a in OrderedMultiPartitionsLE_iterator(sub_vec(vect, b)):
            yield [b]+a


def OrderedMultiPartitions_iterator(vect):
    """
    Iterator of collections of vectors
    (a_1,..,a_k) such that a_1+..+a_k=vect and a_i>0
    """
    if not zero_vec(vect):
        yield [vect]
    for b in IntegerVectors_iterator(vect):
        for a in OrderedMultiPartitions_iterator(sub_vec(vect, b)):
            yield [b]+a


def OrderedPartitionsLE_iterator(n):
    """
    Iterator of collections of positive numbers
    (a_1,..,a_k) such that a_1+..+a_k<=n
    """
    for b in range(1, n+1):
        yield [b]
        for a in OrderedPartitionsLE_iterator(n-b):
            yield [b]+a


def OrderedPartitions_iterator(n):
    """
    Iterator of collections of positive numbers
    (a_1,..,a_k) such that a_1+..+a_k=n
    """
    if n != 0:
        yield [n]
    for b in range(1, n+1):
        for a in OrderedPartitions_iterator(n-b):
            yield [b]+a


def MultiPartitionsLE_iterator(vect, bound=None):
    """
    Iterator of collections of vectors
    b>=a_1>=..>=a_k>0 such that a_1+..+a_k<=vect
    """
    if bound is None:
        bound = vect.copy()
    else:
        bound = min_vec(vect, bound)
    for b in IntegerVectors_iterator(bound):
        yield [b]
        for a in MultiPartitionsLE_iterator(sub_vec(vect, b), b):
            yield [b]+a


def MultiPartitions_iterator(vect, bound=None):
    """
    Iterator of collections of vectors
    b>=a_1>=..>=a_k>0 such that a_1+..+a_k=vect
    """
    if bound is None:
        bound = vect
    if not zero_vec(vect) and le_vec(vect, bound):
        yield [vect]
    bound = min_vec(vect, bound)
    for b in IntegerVectors_iterator(bound):
        for a in MultiPartitions_iterator(sub_vec(vect, b), b):
            yield [b]+a


def Subsets_iterator(n):
    """
    Non-empty subsets of the set {0,..,n-1}
    """
    for i in range(n):  # the maximal element of the subset
        yield [i]
        for s in Subsets_iterator(i):
            yield s+[i]


def ListPartitions_iterator(l):
    """
    Iterator over collections of nonempty lists (l1,..,lk)
    such that l1+..+lk=l
    """
    for a in OrderedPartitions_iterator(len(l)):
        yield list(l[sum(a[:i]):sum(a[:i+1])] for i in range(len(a)))


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
