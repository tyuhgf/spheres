import time

# noinspection PyUnresolvedReferences
import sage.all as sg
from sage.geometry.polyhedron.library import polytopes

from settings import log_handler
from spheres.gg_cocycle import gg
from spheres.simplicial_complex import Sphere

log_handler.setLevel(100)


def get_cyclic_polytope_as_sphere(d, n):
    cyclic_polytope = polytopes.cyclic_polytope(d, n)
    vertices_dict = {v: i for (i, v) in enumerate(cyclic_polytope.vertices())}

    facets = [[vertices_dict[v] for v in f.vertices()] for f in cyclic_polytope.facets()]
    return Sphere(facets)


if __name__ == '__main__':
    t = time.time()
    n_ = 7
    c_n = get_cyclic_polytope_as_sphere(4, n_)
    print(gg(c_n))
    t1 = time.time()
    print('time', t1 - t)
