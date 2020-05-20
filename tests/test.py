from sage import all as sg
from simplicial_complex import Sphere


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
    substitution = {str(v): int(str(v)[1:]) for v in s.vertices()}
    s = s.rename_vertices(substitution)

    assert s.link_oriented([2])  # todo

    s = s.join(sg.SimplicialComplex([[7], [8]])).as_sphere()
    substitution = {str(v): int(str(v)[1:]) for v in s.vertices()}
    s = s.rename_vertices(substitution)

    assert s.link_oriented([2, 6])  # todo


if __name__ == '__main__':
    test_is_sphere()
    test_rename_vertices()
    test_sphere()
