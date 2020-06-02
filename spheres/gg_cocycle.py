from sage.combinat.free_module import CombinatorialFreeModule_Tensor

from spheres.simplicial_complex import *

ChainElement = CombinatorialFreeModule_Tensor.Element


def to_dual_vertices(sphere: Sphere, sim: sg.Simplex, subdivided_sphere: Sphere = None) -> ChainElement:
    """Weighted average of vertices of the dual cell to a simplex in sphere."""
    if subdivided_sphere is None:
        subdivided_sphere = sphere.barycentric_subdivision()
    chains = subdivided_sphere.n_chains(0, sg.QQ)

    base_link = sphere.link_oriented(sim)
    dual_vertices = [tuple(sorted(list(sim) + list(s))) for s in base_link]
    res = chains(0)
    for v in dual_vertices:
        if chains.monomial(v).to_vector == 0:
            raise ValueError('Occurred unknown monomial!')
        res += chains.monomial(v) / len(dual_vertices)
    return res


def dualize_cycle(sphere: Sphere, cycle: ChainElement) -> ChainElement:
    """Make a corresponding cycle in barycentric subdivision."""
    pass


def l2_minimal_chain(sphere: Sphere, cycle: ChainElement) -> ChainElement:
    """Chain with given boundary, minimal by L2 norm."""
    pass


def gg_cocycle(bistellar_move: BistellarMove):
    """GG-cocycle calculation."""
    # s, t = bistellar_move.s, bistellar_move.t
    # base_sphere = bistellar_move.skew_suspension()
    # subdivised_base_sphere = base_sphere.barycentric_subdivision()
    #
    # chains0_base = base_sphere.n_chains(0, sg.QQ)
    # chains1_base = base_sphere.n_chains(1, sg.QQ)
    #
    # chains0 = subdivised_base_sphere.n_chains(0, sg.QQ)
    # chains1 = subdivised_base_sphere.n_chains(1, sg.QQ)

    pass
