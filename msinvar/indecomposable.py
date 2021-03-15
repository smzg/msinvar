# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from sage.misc.misc_c import prod
from msinvar.iterators import MultiPartitionsLE_iterator
from msinvar.utils import vec, phi
from msinvar.invariants import Invariant


def hua_formula(Q):
    """
    Count the number of (absolutely) indecomposable representations of the
    quiver Q.

    We apply the Hua formula as described in arXiv:math/0608321.

    EXAMPLES::

        sage: from msinvar.quivers import Quiver
        ....: from msinvar.indecomposable import hua_formula
        ....: Q=Quiver('1-2') # quiver of type A2
        ....: Q.prec([3,3])
        ....: hua_formula(Q).dict()
        {(1, 0): 1, (0, 1): 1, (1, 1): 1}

        sage: Q=Quiver('1-1') # Jordan quiver
        ....: Q.prec([3])
        ....: hua_formula(Q).dict() # we use q=y^2
        {(1,): y^2, (2,): y^2, (3,): y^2}
    """
    W = Q.wcs()
    q = W.y**2
    dct = hua_formula_dict(Q, q, W.prec())
    return Invariant(W.R(dct).Log()*(q-1))


def hua_formula_dict(Q, q, bound):
    z = tuple([0]*len(bound))
    dct = {z: 1}
    for l in MultiPartitionsLE_iterator(bound):
        d = tuple(vec.add(l))
        p = sum(Q.eform(a, a) for a in l)
        l = l+[z]
        s = q**(-p)/prod(phi(1/q, vec.sub(l[k], l[k+1]))
                         for k in range(len(l)-1))
        if d not in dct:
            dct[d] = 0
        dct[d] += s
    return dct
