from sage import all as sg
from spheres.simplicial_complex import Sphere, BistellarMove
from spheres.gg_cocycle import to_dual_vertices


def test_to_dual_vetrices():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s = s.rename_vertices()
    s = BistellarMove(s, [1, 2, 6]).t
    sim = sg.Simplex([1])
    chain = to_dual_vertices(s, sim)
    assert chain


if __name__ == '__main__':
    test_to_dual_vetrices()
