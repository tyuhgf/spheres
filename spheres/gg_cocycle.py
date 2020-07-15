from spheres.simplicial_complex import *

ChainElement = sg.CombinatorialFreeModule.Element


def boundary_n(chain: ChainElement, n: int, g1, g2):
    """Different boundaries for tensor product of chains."""
    if n not in [0, 1]:
        raise ValueError('Parameter n should be 0 or 1!')

    monomials = chain.monomial_coefficients()

    res = 0
    for (a, b), c in monomials.items():
        a = g1(a)
        b = g2(b)
        if n == 0:
            res += a.boundary().tensor(b) * c
        elif n == 1:
            res += a.tensor(b.boundary()) * c

    return res


def chains_tensor_product(a: sg.CombinatorialFreeModule, b: sg.CombinatorialFreeModule):
    """Tensor product of chain groups equipped with boundary operators."""
    res = a.tensor(b)
    res.boundary_0 = lambda chain: boundary_n(chain, 0, a, b)
    res.boundary_1 = lambda chain: boundary_n(chain, 1, a, b)
    return res


class GGCocycleHelper:
    def __init__(self, bistellar_move: BistellarMove):
        self.bistellar_move = bistellar_move
        self.base_sphere = bistellar_move.skew_suspension()
        self.subdivided_base_sphere = self.base_sphere.barycentric_subdivision()

        self.chains0_base = self.base_sphere.n_chains(0, sg.QQ)
        self.chains1_base = self.base_sphere.n_chains(1, sg.QQ)

        self.chains0 = self.subdivided_base_sphere.n_chains(0, sg.QQ)
        self.chains1 = self.subdivided_base_sphere.n_chains(1, sg.QQ)

        self.chains11 = chains_tensor_product(self.chains1_base, self.chains1_base)

    def to_dual_vertices(self, sim: sg.Simplex) -> ChainElement:
        """Weighted average of vertices of the dual cell to a simplex in base_sphere."""
        base_link = self.base_sphere.link_oriented(sim)
        dual_vertices = [tuple(sorted(list(sim) + list(s))) for s in base_link]
        res = self.chains0.zero()
        for v in dual_vertices:
            res += self.chains0(sg.Simplex([v])) / len(dual_vertices)
        return res

    def cycle_in_chains11_to_product_of_cycles(self, cycle):
        cycle = self.chains11(cycle)
        if not self.chains11.boundary_0(cycle) == 0 == self.chains11.boundary_1(cycle):
            raise ValueError('Chain is not a (delta0+delta1)-cycle!')

        chains0 = self.base_sphere.n_chains(0, sg.QQ)
        chains1 = self.base_sphere.n_chains(1, sg.QQ)
        chains2 = self.base_sphere.n_chains(2, sg.QQ)
        chains3 = self.base_sphere.n_chains(3, sg.QQ)
        chains22 = chains_tensor_product(chains2, chains2)
        chains12 = chains_tensor_product(chains1, chains2)
        chains03 = chains_tensor_product(chains0, chains3)
        chains13 = chains_tensor_product(chains1, chains3)

        d1 = sg.matrix([g.boundary().to_vector() for g in chains1.gens()])
        d2 = sg.matrix([g.boundary().to_vector() for g in chains2.gens()])
        d3 = sg.matrix([g.boundary().to_vector() for g in chains3.gens()])

        cycle11 = dict()  # cycle11 to dict form
        for (a, b), c in cycle.monomial_coefficients().items():
            if a in cycle11:
                cycle11[a] += chains1(b) * c
            else:
                cycle11[a] = chains1(b) * c

        chain12_ = chains12.zero()  # (Id * d2) ^-1 of cycle11
        for a, c in cycle11.items():
            a = chains1(a)
            lift = sg.linear_transformation(d2).lift(c.to_vector())
            lift = chains2.from_vector(lift)
            chain12_ += chains12(a.tensor(lift))

        chain02_ = chains12.boundary_0(chain12_)

        chain02 = dict()  # chain02 to dict form
        for (a, b), c in chain02_.monomial_coefficients().items():
            if a in chain02:
                chain02[a] += chains2(b) * c
            else:
                chain02[a] = chains2(b) * c

        chain03_ = chains03.zero()  # (Id * d3) ^-1 of chain02
        for a, c in chain02.items():
            a = chains0(a)
            lift = sg.linear_transformation(d3).lift(c.to_vector())
            lift = chains3.from_vector(lift)
            chain03_ += chains03(a.tensor(lift))

        chain03 = dict()  # chain03 to dict form
        for (a, b), c in chain03_.monomial_coefficients().items():
            if b in chain03:
                chain03[b] += chains0(a) * c
            else:
                chain03[b] = chains0(a) * c

        chain13_ = chains13.zero()  # (d1 * Id) ^-1 of chain03
        for b, c in chain03.items():
            b = chains3(b)
            lift = sg.linear_transformation(d1).lift(c.to_vector())
            lift = chains1.from_vector(lift)
            chain13_ += chains13(lift.tensor(b))

        chain12_ -= chains13.boundary_1(chain13_)

        chain12 = dict()  # chain12 to dict form
        for (a, b), c in chain12_.monomial_coefficients().items():
            if b in chain12:
                chain12[b] += chains1(a) * c
            else:
                chain12[b] = chains1(a) * c

        chain22_ = chains22.zero()  # (d2 * Id) ^-1 of chain12
        for b, c in chain12.items():
            b = chains2(b)
            lift = sg.linear_transformation(d2).lift(c.to_vector())
            lift = chains2.from_vector(lift)
            chain22_ += chains22(lift.tensor(b))

        return chain22_

    def cycle_in_chains11_to_product_of_cycles_base(self, cycle):
        """Inefficient implementation. Todo: delete method when cycle_in_chains11_to_product_of_cycles tested."""
        cycle = self.chains11(cycle)
        if not self.chains11.boundary_0(cycle) == 0 == self.chains11.boundary_1(cycle):
            raise ValueError('Chain is not a (delta0+delta1)-cycle!')

        chains2_base = self.base_sphere.n_chains(2, sg.QQ)
        chains22 = chains_tensor_product(chains2_base, chains2_base)
        chains21 = chains_tensor_product(chains2_base, self.chains1_base)

        gens = chains22.gens()
        matrix = sg.matrix([chains21.boundary_0(chains22.boundary_1(g)).to_vector() for g in gens])

        # that's super-inefficient!
        res = sg.linear_transformation(matrix).lift(cycle.to_vector())
        return chains22.from_vector(res)

    def dualize_cycle(self, cycle: ChainElement) -> ChainElement:
        """Make a corresponding cycle in barycentric subdivision."""
        pass

    def l2_minimal_chain(self, cycle: ChainElement) -> ChainElement:
        """Chain with given boundary, minimal by L2 norm."""
        pass

    def gg_cocycle(self):
        pass


def gg_cocycle(bistellar_move: BistellarMove):
    """GG-cocycle calculation."""
    return gg_cocycle(bistellar_move).gg_cocycle()
