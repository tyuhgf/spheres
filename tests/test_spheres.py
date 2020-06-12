from sage import all as sg
from spheres.simplicial_complex import Sphere, BistellarMove


def test_is_sphere():
    circle = sg.SimplicialComplex([[1, 2], [1, 3], [2, 3]])
    assert circle.is_homology_sphere()
    assert not sg.SimplicialComplex([[1, 2], [1, 3], [2, 3], [3, 4]]).is_homology_sphere()

    assert circle.suspension().is_homology_sphere()
    assert not circle.product(circle).is_homology_sphere()
    
    Barnette = sg.SimplicialComplex([[ 1, 2, 3, 4 ], [ 1, 2, 3, 5 ], [ 1, 2, 4, 5 ], [ 1, 3, 4, 6 ],  [ 1, 3, 5, 6 ], [ 1, 4, 5, 7 ], [ 1, 4, 6, 7 ], [ 1, 5, 6, 8 ], [ 1, 5, 7, 8 ], [ 1, 6, 7, 8 ], [ 2, 3, 4, 8 ], [ 2, 3, 5, 8 ], [ 2, 4, 5, 7 ], [ 2, 4, 6, 7 ], [ 2, 4, 6, 8 ], [ 2, 5, 7, 8 ], [ 2, 6, 7, 8 ], [ 3, 4, 6, 8 ], [ 3, 5, 6, 8 ]])
    assert Barnette.is_sphere()
    assert Barnette.is_homology_sphere()
    
    Poincare = sg.SimplicialComplex([[ 1, 2, 4, 9 ], [ 1, 2, 4, 15 ], [ 1, 2, 6, 14 ], [ 1, 2, 6, 15 ],  [ 1, 2, 9, 14 ], [ 1, 3, 4, 12 ], [ 1, 3, 4, 15 ], [ 1, 3, 7, 10 ], [ 1, 3, 7, 12 ], [ 1, 3, 10, 15 ], [ 1, 4, 9, 12 ], [ 1, 5, 6, 13 ], [ 1, 5, 6, 14 ], [ 1, 5, 8, 11 ], [ 1, 5, 8, 13 ], [ 1, 5, 11, 14 ], [ 1, 6, 13, 15 ], [ 1, 7, 8, 10 ], [ 1, 7, 8, 11 ], [ 1, 7, 11, 12 ], [ 1, 8, 10, 13 ], [ 1, 9, 11, 12 ], [ 1, 9, 11, 14 ], [ 1, 10, 13, 15 ], [ 2, 3, 5, 10 ], [ 2, 3, 5, 11 ], [ 2, 3, 7, 10 ], [ 2, 3, 7, 13 ], [ 2, 3, 11, 13 ], [ 2, 4, 9, 13 ], [ 2, 4, 11, 13 ], [ 2, 4, 11, 15 ], [ 2, 5, 8, 11 ], [ 2, 5, 8, 12 ], [ 2, 5, 10, 12 ], [ 2, 6, 10, 12 ], [ 2, 6, 10, 14 ], [ 2, 6, 12, 15 ], [ 2, 7, 9, 13 ], [ 2, 7, 9, 14 ], [ 2, 7, 10, 14 ], [ 2, 8, 11, 15 ], [ 2, 8, 12, 15 ], [ 3, 4, 5, 14 ], [ 3, 4, 5, 15 ], [ 3, 4, 12, 14 ], [ 3, 5, 10, 15 ], [ 3, 5, 11, 14 ], [ 3, 7, 12, 13 ], [ 3, 11, 13, 14 ], [ 3, 12, 13, 14 ], [ 4, 5, 6, 7 ], [ 4, 5, 6, 14 ], [ 4, 5, 7, 15 ], [ 4, 6, 7, 11 ], [ 4, 6, 10, 11 ], [ 4, 6, 10, 14 ], [ 4, 7, 11, 15 ], [ 4, 8, 9, 12 ], [ 4, 8, 9, 13 ], [ 4, 8, 10, 13 ], [ 4, 8, 10, 14 ], [ 4, 8, 12, 14 ], [ 4, 10, 11, 13 ], [ 5, 6, 7, 13 ], [ 5, 7, 9, 13 ], [ 5, 7, 9, 15 ], [ 5, 8, 9, 12 ], [ 5, 8, 9, 13 ], [ 5, 9, 10, 12 ], [ 5, 9, 10, 15 ], [ 6, 7, 11, 12 ], [ 6, 7, 12, 13 ], [ 6, 10, 11, 12 ], [ 6, 12, 13, 15 ], [ 7, 8, 10, 14 ], [ 7, 8, 11, 15 ], [ 7, 8, 14, 15 ], [ 7, 9, 14, 15 ], [ 8, 12, 14, 15 ], [ 9, 10, 11, 12 ], [ 9, 10, 11, 16 ], [ 9, 10, 15, 16 ], [ 9, 11, 14, 16 ], [ 9, 14, 15, 16 ], [ 10, 11, 13, 16 ], [ 10, 13, 15, 16 ], [ 11, 13, 14, 16 ], [ 12, 13, 14, 15 ], [ 13, 14, 15, 16 ]])
    assert not Poincare.is_sphere()
    assert Poincare.is_homology_sphere()
    
def test_rename_vertices():
    g = sg.SimplicialComplex([[1, 2, 3], [1, 4, 5], [2, 4], [2, 6]])
    assert g.automorphism_group().order() == 1

    new_names = ['q', 'w', 'e', 'r', 't', 'y']
    substitution = {a: b for (a, b) in zip(range(1, 7), new_names)}
    f = g.rename_vertices(substitution)
    iso = f.is_isomorphic(g, certificate=True)
    assert iso[0]
    assert iso[1] == {v: k for (k, v) in substitution.items()}


def test_sphere():
    s = Sphere([[1, 2], [1, 3], [2, 3]]).join(sg.SimplicialComplex([[5], [6]])).as_sphere()
    s = s.rename_vertices()

    assert Sphere(s.link_oriented([2])).is_valid()

    s = s.join(sg.SimplicialComplex([[7], [8]])).as_sphere()
    s = s.rename_vertices()

    assert Sphere(s.link_oriented([2, 6])).is_valid()


def test_bistellar():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s0 = s.rename_vertices()

    s1 = BistellarMove(s0, [1, 2, 6]).t

    b1 = BistellarMove(s1, [2, 3])
    s2 = b1.t
    assert set(b1.tau) == {6, 7}
    assert s2.link([6, 7]) == sg.SimplicialComplex([[2], [3]])

    b2 = BistellarMove(s2, [6, 7])
    s3 = b2.t
    assert s3.is_isomorphic(s1)

    isom = b2.is_isomorphic(b1.inverse(), certificate=True)
    assert isom[0]
    assert isom[1] == {v: v for v in b1.vertices()}


def test_link_bistellar():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s0 = s.rename_vertices()
    s1 = BistellarMove(s0, [1, 2, 6]).t

    bm = BistellarMove(s1, [3, 4, 7])

    l_sphere = bm.skew_suspension(new_vertices_names=(10, 11))
    l3_sphere = bm.link([3]).skew_suspension(new_vertices_names=(10, 11))

    assert l_sphere.link([10]).is_isomorphic(s1)[0]
    assert l_sphere.link([11]).is_isomorphic(Sphere(bm.t))[0]

    assert l3_sphere.is_isomorphic(Sphere(l_sphere.link_oriented([3])))


def test_chain_in_spheres():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s0 = s.rename_vertices()
    s1 = BistellarMove(s0, [1, 2, 6]).t
    s2 = s1.join(sg.SimplicialComplex([[9], [10]])).as_sphere().rename_vertices()
    s3 = BistellarMove(s2, [8, 2, 6, 10]).t
    chain2 = s3.d()
    chain3 = chain2.d()
    assert len(chain3.data) == 0


def test_barycentric():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s0 = s.rename_vertices()
    s1 = BistellarMove(s0, [1, 2, 6]).t
    s2 = s1.barycentric_subdivision()
    for f in s1.facets_with_orientation:
        assert s2.check_oriented_facet(sg.Simplex([tuple(sorted(f[:i+1])) for i in range(f.dimension() + 1)]))


if __name__ == '__main__':
    test_is_sphere()
    test_rename_vertices()
    test_sphere()
    test_bistellar()
    test_link_bistellar()
    test_chain_in_spheres()
    test_barycentric()
