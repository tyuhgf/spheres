import logging
from concurrent.futures.process import ProcessPoolExecutor

from sage import all as sg
from sage.modules.misc import gram_schmidt

import settings
from spheres.simplicial_complex import Sphere, BistellarMove
from spheres.utils import cochain_monomial_to_list, chains_tensor_product

ChainElement = sg.CombinatorialFreeModule.Element


logger = logging.getLogger(__name__)
logger.setLevel(settings.LOG_LEVEL)
logger.addHandler(settings.log_handler)


class GGCocycleHelper:
    def __init__(self, bistellar_move: BistellarMove):
        self.bistellar_move = bistellar_move
        self.base_sphere = bistellar_move.skew_suspension()

        self.chains0_base = self.base_sphere.n_chains(0, sg.QQ)
        self.chains1_base = self.base_sphere.n_chains(1, sg.QQ)
        self.chains2_base = self.base_sphere.n_chains(2, sg.QQ)
        self.chains3_base = self.base_sphere.n_chains(3, sg.QQ)
        self.cochains2_base = self.base_sphere.n_chains(2, sg.QQ, cochains=True)
        self.cochains3_base = self.base_sphere.n_chains(3, sg.QQ, cochains=True)

        self.chains11 = chains_tensor_product(self.chains1_base, self.chains1_base)

        self._link_spheres = dict()
        self._local_dual_chain = dict()

        self.eta = self.calc_eta()
        self.ksi = self.calc_ksi()

        self.chain22 = self.cycle_in_chains11_to_product_of_cycles(self.ksi)
        self.ggh = self.calc_ggh(self.chain22)

    def calc_ksi(self):
        ksi = self.chains11.zero()
        for a, q in self.eta:
            ksi -= a.tensor(a) * q * 2
        for a1, q1 in self.eta:
            for a2, q2 in self.eta:
                ksi += a1.tensor(a2) * q1 * q2
        return ksi

    def calc_ggh(self, chain22):
        logger.info('start')

        result = 0
        for (a, b), c in chain22.monomial_coefficients().items():
            b_ = self.dualize_cycle(self.chains2_base(b).boundary())
            ab = self.chains2_base(a).to_vector() * b_.to_vector()
            result -= ab * c
        logger.info('finish')
        return result

    def calc_eta(self):
        eta = []
        max_vertex = max([i for i in self.base_sphere.vertices() if isinstance(i, int)])
        a, b = max_vertex - 1, max_vertex

        for v in set(self.bistellar_move.s.vertices()) - set(self.bistellar_move.sigma) - set(self.bistellar_move.tau):
            chain = self.chains1_base(sg.Simplex([v, a])) * -1 + self.chains1_base(sg.Simplex([v, b]))
            coef = 1 - sg.Rational(len(self.bistellar_move.s.link(sg.Simplex([v])).vertices())) / 6
            eta.append((chain, coef))
        for v in self.bistellar_move.sigma:
            if v in self.bistellar_move.t.vertices():
                chain = self.chains1_base(sg.Simplex([v, a])) * -1 + self.chains1_base(sg.Simplex([v, b]))
                coef = 1 - sg.Rational(len(self.bistellar_move.t.link(sg.Simplex([v])).vertices())) / 6
                eta.append((chain, coef))
        for v in self.bistellar_move.tau:
            if v in self.bistellar_move.s.vertices():
                chain = self.chains1_base(sg.Simplex([v, a])) * -1 + self.chains1_base(sg.Simplex([v, b]))
                coef = 1 - sg.Rational(len(self.bistellar_move.s.link(sg.Simplex([v])).vertices())) / 6
                eta.append((chain, coef))

        for v1 in self.bistellar_move.sigma:
            for v2 in self.bistellar_move.tau:
                chain = self.chains1_base(sg.Simplex([v1, a])) * -1 + \
                        self.chains1_base(sg.Simplex([v1, v2])) * (-1) ** (v2 < v1) + \
                        self.chains1_base(sg.Simplex([v2, b]))
                coef = sg.Rational(-1)/12 if len(self.bistellar_move.tau) == 2 else sg.Rational(1)/6
                eta.append((chain, coef))

        return eta

    def to_dual_vertices(self, sim: sg.Simplex) -> ChainElement:
        """Weighted average of vertices of the dual cell to a simplex in base_sphere."""
        base_link = self.base_sphere.link_oriented(sim)
        dual_vertices = [tuple(sorted(list(sim) + list(s))) for s in base_link]
        res = self.cochains3_base.zero()
        for v in dual_vertices:
            res += self.cochains3_base(sg.Simplex(v)) / len(dual_vertices)
        return res

    def cycle_in_chains11_to_product_of_cycles(self, cycle):
        logger.info('start')

        cycle = self.chains11(cycle)
        if not self.chains11.boundary_0(cycle) == 0 == self.chains11.boundary_1(cycle):
            raise ValueError('Chain is not a (delta0+delta1)-cycle!')

        chains0 = self.chains0_base
        chains1 = self.chains1_base
        chains2 = self.chains2_base
        chains3 = self.chains3_base
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

    def dualize_cycle(self, cycle: ChainElement) -> ChainElement:
        """Make a corresponding cycle in barycentric subdivision."""
        res = self.base_sphere.n_chains(2, sg.QQ, cochains=True).zero()
        for (k, q) in cycle.monomial_coefficients().items():
            v1, v2 = k.tuple()
            if (v1, v2) not in self._local_dual_chain.keys():
                c1 = self.to_dual_vertices(sg.Simplex([v1])) - self.to_dual_vertices(sg.Simplex([v1, v2]))
                c2 = self.to_dual_vertices(sg.Simplex([v2])) - self.to_dual_vertices(sg.Simplex([v1, v2]))
                self._local_dual_chain[(v1, v2)] = (self.l2_minimal_chain(c1, local=[v1]) -
                                                    self.l2_minimal_chain(c2, local=[v2]))
            res += self._local_dual_chain[(v1, v2)] * q
            logger.info([v1, v2])
        return res

    def l2_minimal_chain(self, cycle: ChainElement, local=None) -> ChainElement:
        """Chain with given boundary, minimal by L2 norm."""
        if local is None or len(local) != 1:
            raise NotImplemented

        v = local[0]
        sim = sg.Simplex([v])

        if v in self._link_spheres.keys():
            lk_sphere, chains0, chains1, ker, d1 = self._link_spheres[v]
        else:
            lk_sphere = Sphere(self.base_sphere.link_oriented(sim))
            chains1 = lk_sphere.n_chains(1, sg.QQ, cochains=True)
            chains0 = lk_sphere.n_chains(2, sg.QQ, cochains=True)
            d1 = []
            for g in chains1.gens():
                d1.append(chains0.zero())
                for m in g.coboundary().monomials():
                    m = cochain_monomial_to_list(m)
                    n = cochain_monomial_to_list(g.monomials()[0])
                    w = list(set(m) - set(n))[0]
                    coef = (-1)**(m.index(w) + lk_sphere.check_oriented_facet(m))
                    d1[-1] += chains0(sg.Simplex(m)) * coef
                d1[-1] = d1[-1].to_vector()
            d1 = sg.matrix(d1)
            ker = d1.left_kernel().basis()
            ker, _ = gram_schmidt(ker)
            self._link_spheres[v] = lk_sphere, chains0, chains1, ker, d1

        cycle_in_link = chains0.zero()  # transform chains0 to a (dual) 0-cycle in link of v
        for (m, c) in cycle.monomial_coefficients().items():
            m = list(m)
            # c *= (-1)**self.base_sphere.check_oriented_facet(m)
            m.remove(v)
            cycle_in_link += chains0(sg.Simplex(m)) * c

        lift = sg.linear_transformation(d1).lift(cycle_in_link.to_vector())

        for k in ker:
            if k:
                lift -= k * (lift * k) / (k * k)
        chain_in_link = chains1.from_vector(lift)
        res = self.cochains2_base.zero()
        for (m, c) in chain_in_link.monomial_coefficients().items():
            m = list(m)
            res += self.cochains2_base(sg.Simplex(m + [v])) * c * (-1) ** sorted(m + [v]).index(v)

        return res

    def inverse_boundary(self, cycle: ChainElement) -> ChainElement:
        """Find 2-chain with given boundary."""
        chains2 = self.base_sphere.n_chains(2, sg.QQ)

        d2 = sg.matrix([g.boundary().to_vector() for g in chains2.gens()])
        lift = sg.linear_transformation(d2).lift(cycle.to_vector())
        return chains2.from_vector(lift)


def gg_cocycle(bistellar_move: BistellarMove):
    """GG-cocycle calculation."""
    ggh = GGCocycleHelper(bistellar_move)
    return ggh.ggh


def gg(sphere, max_workers=1, path_to_simplex_timeout=20):
    if sphere.is_minus_self():
        return 0
    moves = sphere.path_to_simplex(timeout=path_to_simplex_timeout)
    logger.info(f'n_moves: {len(moves)}')
    logger.info([len(move.vertices()) for move in moves])
    res = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for move in moves:
            if move.s.is_minus_self():
                break
            if len(move.sigma) > 1:
                for v in move.sigma:
                    ggh = executor.submit(GGCocycleHelper, move.link([v]))
                    res.append(ggh)
            if len(move.tau) > 1:
                for v in move.tau:
                    ggh = executor.submit(GGCocycleHelper, move.link([v]))
                    res.append(ggh)

    return sum(ggh.result().ggh for ggh in res)
