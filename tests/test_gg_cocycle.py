from sage import all as sg

from settings import log_handler
from spheres.gg_cocycle import GGCocycleHelper, chains_tensor_product
from spheres.simplicial_complex import Sphere, BistellarMove

log_handler.setLevel(100)


def test_to_dual_vertices():
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


# def test_l2_minimal_chain():
#     s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
#     s = s.rename_vertices()
#     bm = BistellarMove(s, [1, 2, 6])
#
#     ggh = GGCocycleHelper(bm)
#
#     chains0 = ggh.base_sphere.n_chains(0, sg.QQ)
#
#     c = chains0(sg.Simplex([1])) +\
#         chains0(sg.Simplex([2])) -\
#         chains0(sg.Simplex([6])) * 2
#
#     q = ggh.l2_minimal_chain(c, [8])
#
#     assert q.boundary() == c
#     print('\n', q, '\n', q.to_vector().norm())


# def test_dualize_cycle():
#     s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
#     s = s.rename_vertices()
#     bm = BistellarMove(s, [1, 2, 6])
#
#     ggh = GGCocycleHelper(bm)
#
#     chains1 = ggh.base_sphere.n_chains(1, sg.QQ)
#
#     c = chains1(sg.Simplex([1, 2])) +\
#         chains1(sg.Simplex([2, 6])) -\
#         chains1(sg.Simplex([1, 6]))
#
#     q = ggh.dualize_cycle(c)
#     print('\n', q)


def test_calc_eta():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s = s.rename_vertices()
    s = BistellarMove(s, [2, 3, 7]).t
    bm = BistellarMove(s, [3, 4, 6])

    ggh = GGCocycleHelper(bm)

    ksi = ggh.ksi
    assert ggh.chains11.boundary_0(ksi) == 0 == ggh.chains11.boundary_1(ksi)


def test_gg_cocycle():
    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s0 = s.rename_vertices()
    bm0 = BistellarMove(s0, [2, 3, 7])  # 8
    s1 = bm0.t
    bm1 = BistellarMove(s1, [4, 5, 7])  # 8
    s2 = bm1.t
    bm2 = BistellarMove(s2, [8])
    s3 = bm2.t
    bm3 = BistellarMove(s3, [9])

    ggh0 = GGCocycleHelper(bm0)
    ggh1 = GGCocycleHelper(bm1)
    ggh2 = GGCocycleHelper(bm2)
    ggh3 = GGCocycleHelper(bm3)
    res = 0
    for ggh in [ggh0, ggh1, ggh2, ggh3]:
        res += ggh.ggh
    assert res == sg.Rational(1) / 210


def test_cp2():
    f = [[int(d) for d in str(n)] for n in
         [1243, 1237, 1276, 2354, 2376, 3476, 3465, 4576, 2385, 2368,
          5386, 4285, 4875, 4817, 4371, 7165, 1785, 1586, 1682, 1284]]
    s = Sphere(f)
    moves = s.path_to_simplex(10)

    res = []
    for move in moves:
        if len(move.sigma) > 1:
            for v in move.sigma:
                ggh = GGCocycleHelper(move.link([v]))
                res.append(ggh.ggh)
        if len(move.tau) > 1:
            for v in move.tau:
                ggh = GGCocycleHelper(move.link([v]))
                res.append(ggh.ggh)

    assert sum(res) == sg.Rational(1)/3


def test_ggh_rename():
    s = Sphere([[1, 4, 8], [1, 3, 4], [1, 8, 6], [4, 6, 8],
                [1, 6, 2], [6, 4, 3], [1, 2, 3], [3, 2, 6]])

    subst = {1: 5, 2: 7, 3: 4, 4: 3, 6: 1, 8: 8, 9: 9, 10: 10}

    move = BistellarMove(s, [1, 3])
    ggh1 = GGCocycleHelper(move)
    ggh_r = GGCocycleHelper(BistellarMove(s.rename_vertices(subst), [subst[q] for q in [1, 3]]))

    assert ggh1.ggh == ggh_r.ggh


def test_ggh_inverse():
    s = Sphere([[1, 4, 8], [1, 3, 4], [1, 8, 6], [4, 6, 8],
                [1, 6, 2], [6, 4, 3], [1, 2, 3], [3, 2, 6]])

    move = BistellarMove(s, [1, 3])
    ggh1 = GGCocycleHelper(move)
    ggh_r = GGCocycleHelper(move.inverse())

    assert ggh1.ggh == -ggh_r.ggh


if __name__ == '__main__':
    test_to_dual_vertices()
    test_chains11_to_products()
    # test_l2_minimal_chain()
    # test_dualize_cycle()
    test_calc_eta()
    test_cp2()
