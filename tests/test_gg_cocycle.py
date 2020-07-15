from sage import all as sg
from spheres.simplicial_complex import Sphere, BistellarMove
from spheres.gg_cocycle import GGCocycleHelper, chains_tensor_product


def test_to_dual_vetrices():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s = s.rename_vertices()
    s = BistellarMove(s, [1, 2, 6]).t
    bm = BistellarMove(s, [2, 6])

    ggh = GGCocycleHelper(bm)

    sim = sg.Simplex([1])
    q = ggh.to_dual_vertices(sim)

    print(q)


def test_chains11_to_products():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s = s.rename_vertices()
    bm = BistellarMove(s, [1, 2, 6])

    ggh = GGCocycleHelper(bm)

    chains1 = ggh.base_sphere.n_chains(1, sg.QQ)
    chains2 = ggh.base_sphere.n_chains(2, sg.QQ)
    chains22 = chains_tensor_product(chains2, chains2)
    chains12 = chains_tensor_product(chains1, chains2)

    c = chains1(sg.Simplex([1, 2])) +\
        chains1(sg.Simplex([2, 6])) -\
        chains1(sg.Simplex([1, 6]))

    c11 = c.tensor(c)

    q22 = ggh.cycle_in_chains11_to_product_of_cycles(c11)

    q12 = chains22.boundary_0(q22)
    q11 = chains12.boundary_1(q12)

    assert q11 == c11


if __name__ == '__main__':
    test_to_dual_vetrices()
