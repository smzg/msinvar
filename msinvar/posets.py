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


class Poset():
    def __init__(self, rel, vert=None):
        if vert is None:
            vert = set()
            for r in rel:
                vert.update(r)
        self.vert = vert
        self.rel = rel

    def upper_closure(self, vert):
        vert = set(vert)
        new = vert
        while len(new) > 0:
            for i in new:
                yield i
            new = set(r[1] for r in self.rel if r[0] in new)-vert
            vert.update(new)

    def lower_closure(self, vert):
        vert = set(vert)
        new = vert
        while len(new) > 0:
            for i in new:
                yield i
            new = set(r[0] for r in self.rel if r[1] in new)-vert
            vert.update(new)

    def succ(self, i): return self.upper_closure({i})
    def pred(self, i): return self.lower_closure({i})

    def ideals(self, size=100):
        """
        List if ideals in a poset having <=size elements
        """
        if size<0:
            return
        if len(self.vert) == 0 or size == 0:
            yield set()
            return
        nonmin = set(r[1] for r in self.rel)
        for i in self.vert-nonmin:
            break  # choose a minimal element
        # yield {i}
        rel = list(r for r in self.rel if r[0] != i)  # relations without i
        P = Poset(rel, self.vert-{i})  # for ideals with i
        for I in P.ideals(size-1):
            yield I | {i}
        vert = self.vert-set(self.succ(i))
        rel = list(r for r in self.rel if r[1] in vert)  # relations on vert
        P = Poset(rel, vert)  # for ideals without i
        for I in P.ideals(size):
            yield I

    def __repr__(self):
        return "Poset with vertices\n"+str(self.vert)\
            +"\nand relations\n"+str(self.rel)