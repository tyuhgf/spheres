import time

# noinspection PyUnresolvedReferences
import sage.all as sg
from sage.geometry.polyhedron.library import polytopes

from settings import log_handler
from spheres.gg_cocycle import gg
from spheres.simplicial_complex import Sphere, BistellarMove

log_handler.setLevel(100)


def get_cyclic_polytope_as_sphere(d, n):
    cyclic_polytope = polytopes.cyclic_polytope(d, n)
    vertices_dict = {v: i for (i, v) in enumerate(cyclic_polytope.vertices())}

    facets = [[vertices_dict[v] for v in f.vertices()] for f in cyclic_polytope.facets()]
    return Sphere(facets)


def main():
    t = time.time()
    n = 7

    move_simplices = []
    for i in range(n-4):
        for j in range(i + 2, n-3):
            move_simplices.append([i, j])
        if i < n - 5:
            move_simplices.append([i])

    s0 = s = get_cyclic_polytope_as_sphere(4, n)
    moves = []
    for ms in move_simplices:
        moves.append(BistellarMove(s, ms))
        s = moves[-1].t

    res = gg(s0, 3, moves=moves, result='all')
    t1 = time.time()
    print('time', t1 - t)

    move_values = [0 for _ in moves]
    for (i, v, g) in res:
        move_values[i] += g

    with open('res_moves.txt', 'a') as f:
        for i in range(len(move_values)):
            f.write(str(move_simplices[i]))
            f.write('  ')
            f.write(str(move_values[i]))
            f.write('\n')

    with open('res_links.txt', 'a') as f:
        for (i, v, g) in res:
            f.write(str(move_simplices[i]))
            f.write('  ')
            f.write(str(v))
            f.write('  ')
            f.write(str(g))
            f.write('\n')


if __name__ == '__main__':
    main()
