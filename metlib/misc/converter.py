#!/usr/bin/env python2.7

# converter.py

from metlib.misc.datatype import Null, isnull

class Converter(object):
    """A class for automatically generating converting functions between each node-pair.
    """
    def __init__(self, relations=[], funcname_sep='_to_'):
        """
        relations: a list of ('nodename1', 'nodename2', func) tuple.
        funcname_sep: sep string between the 2 names.
        """
        self.nodes = set()
        self.relations = dict()
        #self.changed_relations = []
        self.funcname_sep = funcname_sep

        for n1, n2, rel_func in relations:
            self.add_relation(n1, n2, rel_func)

    def add_relation(self, node1, node2, relation):
        self.nodes.add(node1)
        self.nodes.add(node2)
        self.relations[(node1, node2)] = relation
    
    def update(self):
        changed = 1
        while changed:
            changed = self._iter_once()

    def _iter_once(self):
        changed = 0
        for n1, n2 in self.relations.keys():
            n1n2_rel = self.relations.get((n1, n2), None)
            if not isnull(n1n2_rel):
                for n3 in self.nodes - set([n1, n2]):
                    n1n3_rel = self.relations.get((n1, n3))
                    n2n3_rel = self.relations.get((n2, n3))
                    if isnull(n1n3_rel) and not isnull(n2n3_rel):
                        def newfunc(value, *args, **kwargs):
                            return n2n3_rel(n1n2_rel(value, *args, **kwargs), *args, **kwargs)
                        newfunc.func_name = '%s%s%s' % (n1, self.funcname_sep, n3)
                        newfunc.func_doc = "Conversion func between %s and %s .\nGenerated via %s and %s ."  \
                            % (n1, n3, n1n2_rel.func_name, n2n3_rel.func_name)
                        self.relations[(n1, n3)] = newfunc
                        changed += 1
        return changed
    
    def __call__(self, n1, n2):
        return self.relations.get((n1, n2), Null)
        
if __name__ == '__main__':
    def A2B(A):
        return "A(%s) to B" % A

    def B2C(B):
        return "B(%s) to C" % B

    def C2D(C):
        return "C(%s) to D" % C

    def D2A(D):
        return "D(%s) to A" % D
        
    def E2C(E):
        return "E(%s) to C" % E
        
    def A2E(A):
        return "A(%s) to E" % A

    c = Converter( [('A', 'B', A2B),
                    ('B', 'C', B2C),
                    ('C', 'D', C2D),
                    ('D', 'A', D2A)],
                    funcname_sep='2',
                    )

    c.update()

    print c.relations
