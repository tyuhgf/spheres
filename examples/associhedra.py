import time

import numpy as np
import sage.all as sg
from sage.geometry.polyhedron.constructor import Polyhedron

from settings import log_handler
from spheres.gg_cocycle import gg
from spheres.simplicial_complex import Sphere

log_handler.setLevel(100)


def graph_associhedron(g):
    """
    surface of the graph-associhedron as a sphere
    """
    n = len(g.vertices())
    gen_set = [g.subgraph(np.nonzero([(mask >> i) % 2 for i in range(n)])[0]).is_connected() for mask in range(1 << n)]
    k = [0] * (1 << n)
    for mask in range(1 << n):
        for mask2 in range(1, mask + 1):
            if mask & mask2 == mask2:
                k[mask] += gen_set[mask2] * gen_set[mask]
    inequalities = []
    for mask in range(1, 1 << n):
        if gen_set[mask]:
            r = [(mask >> i) % 2 for i in range(n)]
            inequalities.append([-1 * k[mask]] + r)
    s = Polyhedron(Polyhedron(ieqs=inequalities).vertices())

    facets = list(s.facets())
    s_dual = [[i for (i, f) in enumerate(facets) if v in f.vertices()] for v in s.vertices()]

    return Sphere(s_dual)


def run_all_graphs():
    gs = []
    for mask in range(1024):
        edges = []
        for v in range(5):
            for w in range(v):
                if mask % 2:
                    edges.append([w, v])
                mask >>= 1
        g = sg.Graph(edges)
        if len(g.vertices()) < 5 or not g.is_connected():
            continue
        flag = True
        for g_ in gs:
            if g.is_isomorphic(g_):
                flag = False
        if flag:
            gs.append(g)

    for g in gs:
        with open('qqq.txt', 'a') as f:
            f.write(str(list(g.edges(labels=False))))
            f.write('  ')

        s = graph_associhedron(g)

        t = time.time()
        res = gg(s, 15)
        t1 = time.time()
        total = t1 - t

        with open('qqq.txt', 'a') as f:
            f.write(str(res))
            f.write('  ')
            f.write(str(total))
            f.write('\n')


if __name__ == '__main__':
    run_all_graphs()
