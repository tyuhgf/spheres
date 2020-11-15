import time
from itertools import permutations

import numpy as np
# noinspection PyUnresolvedReferences
import sage.all as sg

from settings import log_handler
from spheres.gg_cocycle import gg
from spheres.simplicial_complex import Sphere

log_handler.setLevel(100)


class CompSphere(Sphere):
    """
    n 'bordisms' of 0-sphere yield a n-sphere
    """
    def __init__(self, k):
        self.k = k = np.array(k)
        n = len(k) + 1
        vertices_down = [0] + list(np.cumsum(k[:, 0]))
        vertices_up = [0] + list(np.cumsum(k[:, 1]))
        facets = []

        for p in permutations(range(n)):
            f = []
            for i in range(n-1):
                f.append('s_' + '_'.join([str(v) for v in sorted(p[:i + 1])]))

            facets.append([f'd_{vertices_down[p[0]]}'] + f)
            facets.append([f'u_{vertices_up[p[0]]}'] + f)

            if p[0] < p[1]:
                for v in range(vertices_down[p[0]], vertices_down[p[1]]):
                    facets.append([f'd_{v}', f'd_{v+1}'] + f[1:])
                for v in range(vertices_up[p[0]], vertices_up[p[1]]):
                    facets.append([f'u_{v}', f'u_{v+1}'] + f[1:])

        super().__init__(facets)


if __name__ == '__main__':
    for m_ in range(2, 4):
        for n_ in range(1, 5):
            k_ = [[0, 2], [0, n_], [0, m_]]
            print(k_)

            s = CompSphere(k_)
            s = s.rename_vertices({v: i for (i, v) in enumerate(s.vertices())})

            t = time.time()
            print(gg(s, 15))
            t1 = time.time()
            print('time', t1 - t)
