r"""
Posets

EXAMPLES::

    sage: from msinvar.posets import Poset
    sage: P=Poset([[1,2],[1,3]]); P
    Poset with vertices
    {1, 2, 3}
    and relations
    [[1, 2], [1, 3]]
    sage: list(P.succ(1))
    [1, 2, 3]
    sage: list(P.ideals())
    [{1, 2, 3}, {1, 2}, {1, 3}, {1}, set()]
"""
# *****************************************************************************
#  Copyright (C) 2021 Sergey Mozgovoy <mozhov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************


class Poset:
    """
    Poset class.

    EXAMPLES::

        sage: from msinvar.posets import Poset
        sage: P=Poset([[1,2],[1,3]]); P
        Poset with vertices
        {1, 2, 3}
        and relations
        [[1, 2], [1, 3]]
        sage: list(P.succ(1))
        [1, 2, 3]
        sage: list(P.ideals())
        [{1, 2, 3}, {1, 2}, {1, 3}, {1}, set()]
    """

    def __init__(self, rel, vert=None):
        if vert is None:
            vert = set().union(*rel)
        self.vert = vert
        self.rel = rel

    def opposite(self):
        return Poset(([r[1], r[0]] for r in self.rel), self.vert)

    def upper_closure(self, vert):
        """
        upper ideal
        """
        vert = set(vert)
        new = vert
        while new:
            yield from new
            new = set(r[1] for r in self.rel if r[0] in new) - vert
            vert.update(new)

    def lower_closure(self, vert):
        """
        lower ideal
        """
        return self.opposite().upper_closure(vert)

    def succ(self, i):
        """
        principal upper ideal
        """
        return self.upper_closure({i})

    def pred(self, i):
        """
        principal lower ideal
        """
        return self.lower_closure({i})

    def ideals(self, size=100):
        """
        List of ideals in a poset having <= ``size`` elements.
        """
        if len(self.vert) == 0 or size == 0:
            yield set()
            return
        nonmin = set(r[1] for r in self.rel)
        i = next(iter(self.vert - nonmin))  # choose a minimal element
        rel = list(r for r in self.rel if r[0] != i)  # relations without i
        P = Poset(rel, self.vert - {i})  # for ideals with i
        for I in P.ideals(size - 1):
            yield I | {i}
        vert = self.vert - set(self.succ(i))
        rel = list(r for r in self.rel if r[1] in vert)  # relations on vert
        P = Poset(rel, vert)  # for ideals without i
        yield from P.ideals(size)

    def __repr__(self):
        return (
            "Poset with vertices "
            + str(self.vert)
            + "\nand relations "
            + str(self.rel)
        )
