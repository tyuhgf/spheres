from sage import all as sg
from spheres.simplicial_complex import Sphere, BistellarMove


def test_is_sphere():
    circle = sg.SimplicialComplex([[1, 2], [1, 3], [2, 3]])
    assert circle.is_homology_sphere()
    assert not sg.SimplicialComplex([[1, 2], [1, 3], [2, 3], [3, 4]]).is_homology_sphere()

    assert circle.suspension().is_homology_sphere()
    assert not circle.product(circle).is_homology_sphere()


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
        assert s2.check_oriented_facet(sg.Simplex([tuple(sorted(f[:i+1])) for i in range(len(f))]))


if __name__ == '__main__':
    test_is_sphere()
    test_rename_vertices()
    test_sphere()
    test_bistellar()
    test_link_bistellar()
    test_chain_in_spheres()
    test_barycentric()
