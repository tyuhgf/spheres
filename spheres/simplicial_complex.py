import re
from typing import List, Union, Set, Iterable, Tuple

from sage import all as sg


def is_homology_sphere(self: sg.SimplicialComplex):
    if not self.is_pseudomanifold() or not self.is_connected():
        return False
    h = self.homology()
    if not all([h[i].order() == 1 for i in range(1, self.dimension())]):
        return False
    if not h[self.dimension()].short_name() == 'Z':
        return False
    return True


def rename_vertices(self: sg.SimplicialComplex, subst: Union[str, dict, callable] = 'int'):
    if subst == 'int':
        subst = {v: int(''.join(re.findall(r'\d', str(v)))) for v in self.vertices()}
        names_function = (lambda x: subst[x] if x in subst else x)
    elif isinstance(subst, dict):
        names_function = (lambda x: subst[x] if x in subst else x)
    elif callable(subst):
        names_function = subst
    else:
        raise ValueError('Argument subst has unknown type!')
    return sg.SimplicialComplex([[names_function(v) for v in f] for f in self.facets()])


def as_sphere(self: sg.SimplicialComplex):
    return Sphere(self.facets())


def change_orientation(self: Union[sg.Simplex, List]):
    s = self.copy() if isinstance(self, list) else list(self)
    if len(s) >= 2:
        s = s[:-2] + s[-2:][::-1]
    if isinstance(self, list):
        return s
    else:
        return sg.Simplex(s)


basic_is_isomorphic = sg.SimplicialComplex.is_isomorphic


def is_isomorphic(self, other, certificate=False):
    res = basic_is_isomorphic(self, other)
    if not certificate:
        return res, None
    else:
        if not res:
            return res, None
        return basic_is_isomorphic(self, other, certificate)


sg.SimplicialComplex.is_homology_sphere = is_homology_sphere
sg.SimplicialComplex.rename_vertices = rename_vertices
sg.SimplicialComplex.as_sphere = as_sphere
sg.SimplicialComplex.is_isomorphic = is_isomorphic  # that is due to bug in sage

sg.Simplex.change_orientation = change_orientation


class Sphere(sg.SimplicialComplex):
    def __init__(self, data: Union[Iterable[List], sg.SimplicialComplex], **kwargs):
        if isinstance(data, sg.SimplicialComplex):
            data = data.facets()

        super(Sphere, self).__init__(data, **kwargs)
        if not self.is_homology_sphere():
            raise ValueError('Complex is not a sphere!')
        self.set_immutable()

        if not self.facets() == set(sg.Simplex(s) for s in data):
            raise ValueError('Data is not list of facets!')

        self.facets_with_orientation = self.orient(list(data)[0])

    def rename_vertices(self, subst='int'):
        return super(Sphere, self).rename_vertices(subst).as_sphere()

    def is_isomorphic(self, other, certificate=False, with_multiple: bool = False):
        res = super(Sphere, self).is_isomorphic(other, certificate=(certificate or with_multiple))
        if with_multiple:
            if not res[0]:
                if certificate:
                    return False, None, None
                else:
                    return False, None
            f = self.facets_with_orientation[0]
            f = [res[1][v] for v in f]
            oriented = other.check_oriented_facet(f)
            oriented = 1 if oriented else -1
            if certificate:
                return res[0], res[1], oriented
            else:
                res = res[0], oriented
        return res

    def is_valid(self) -> bool:
        if not self.is_homology_sphere():
            raise ValueError('Complex is not a sphere!')
        if not set(self.facets_with_orientation) == self.facets():
            raise ValueError('Facets not equals to calculated facets_with_orientation!')
        for e in self.flip_graph().edges():
            v1 = list(set(e[0]) - set(e[1]))
            v2 = list(set(e[1]) - set(e[0]))
            if not len(v1) == len(v2) == 1:
                raise ValueError
            v1 = v1[0]
            v2 = v2[0]
            e0_oriented = [f for f in self.facets_with_orientation if f == sg.Simplex(e[0])][0]
            neg = [v if v != v1 else v2 for v in e0_oriented]
            pos = change_orientation(neg)
            if not self.check_oriented_facet(sg.Simplex(pos)):
                raise ValueError('Not a valid orientation!')
        return True

    def orient(self, oriented_facet) -> List[sg.Simplex]:
        """Create list of oriented facets by BFS starting from oriented_facet."""
        if not isinstance(oriented_facet, sg.Simplex):
            oriented_facet = sg.Simplex(oriented_facet)
        g = self.flip_graph()
        edges = list(g.breadth_first_search(oriented_facet, edges=True))
        oriented_facets = [oriented_facet]
        for e in edges:
            v1 = list(set(e[0]) - set(e[1]))
            v2 = list(set(e[1]) - set(e[0]))
            if not len(v1) == len(v2) == 1:
                raise ValueError
            v1 = v1[0]
            v2 = v2[0]
            e0_oriented = [f for f in oriented_facets if f == sg.Simplex(e[0])][0]
            neg = [v if v != v1 else v2 for v in e0_oriented]
            pos = change_orientation(neg)
            oriented_facets.append(sg.Simplex(pos))
        return oriented_facets

    def d(self):
        return Chain(Sphere, [(1, Sphere(self.link_oriented([i]))) for i in self.vertices()])

    def link_oriented(self, v: Union[sg.Simplex, List]) -> Set[sg.Simplex]:
        """:return: set of simplices s of link(self, v) with orientation such that s+v is oriented simplex of self"""
        facets = list(self.link(v).facets())
        if not isinstance(v, list):
            v = list(v)
        res = set()
        for s in facets:
            s = list(s)
            if len(s) <= 1 or self.check_oriented_facet(sg.Simplex(s + v)):
                res.add(sg.Simplex(s))
            else:
                res.add(sg.Simplex(s).change_orientation())
        return res

    def check_oriented_facet(self, f: Union[List, sg.Simplex]) -> bool:
        if not isinstance(f, sg.Simplex):
            f = sg.Simplex(f)
        q = [s for s in self.facets_with_orientation if sg.Simplex(s) == f]
        if not len(q) == 1:
            raise ValueError('Not a facet!')
        d = {v: i + 1 for (i, v) in enumerate(list(f))}
        perm = [d[v] for v in list(q[0])]
        if sg.Permutation(perm).is_even():
            return True
        return False

    def path_to_simplex(self):
        # todo
        raise NotImplemented


class BistellarMove(sg.SimplicialComplex):
    def __init__(self, sphere: Sphere, sigma: List, new_vertex_name=None):
        s = sphere.__copy__().as_sphere().facets_with_orientation
        super(BistellarMove, self).__init__(s)

        if new_vertex_name is None:
            new_vertex_name = max([i for i in self.vertices() if isinstance(i, int)]) + 1

        self.dim = self.dimension()  # dimension of a sphere
        self.sigma = sigma

        self.s = s

        v = list(self.star(sigma).vertices())
        if len(v) == self.dimension() + 2:
            pass
        elif len(v) == self.dimension() + 1 == len(sigma) and len(self.vertices()) > self.dimension() + 1:
            v.append(new_vertex_name)
        else:
            raise ValueError('Cannot apply bistellar move!')

        self.tau = list(set(v) - set(sigma))

        self.add_face(v)

        self.t = set(self.s).copy()
        sim = sg.SimplicialComplex([v])
        sim.remove_face(v)
        facet = [facet for facet in self.s if facet in sim.facets()][0]
        sim = sim.as_sphere().orient(oriented_facet=facet)

        for facet in sim:
            if facet in self.t:
                self.t.remove(facet)
            else:
                self.t.update({facet.change_orientation()})

    def is_isomorphic(self, other, certificate=False):
        c_self = Sphere(self.s).cone()\
            .rename_vertices({'R0': '__new_vertex__'})\
            .rename_vertices({'L' + str(v): v for v in self.vertices()})
        r_self = sg.SimplicialComplex(self.facets() + c_self.facets())

        c_other = Sphere(other.s).cone()\
            .rename_vertices({'R0': '__new_vertex__'})\
            .rename_vertices({'L' + str(v): v for v in other.vertices()})
        r_other = sg.SimplicialComplex(other.facets() + c_other.facets())

        res = r_self.is_isomorphic(r_other, certificate=certificate)

        if certificate and res[0]:
            if not res[1]['__new_vertex__'] == '__new_vertex__':
                raise ValueError('Invalid arguments!')
            res[1].pop('__new_vertex__')

        return res

    def inverse(self):
        new_vertex_name = self.sigma[0] if len(self.sigma) == 1 else None
        return BistellarMove(Sphere(self.t), self.tau, new_vertex_name=new_vertex_name)

    def link(self, simplex, is_mutable=True):
        res = Sphere(self.s).link_oriented(simplex)
        sigma = list(set(self.sigma) - set(simplex))
        new_vertex_name = self.sigma[0] if len(self.sigma) == 1 else None
        return BistellarMove(Sphere(res), sigma, new_vertex_name)

    def delta(self):
        return Chain(Sphere, [(1, self.t), (-1, self.s)])


class Chain:
    def __init__(self, meta, data: List[Tuple], ring=sg.Integers()):
        self.meta = meta
        self.ring = ring
        self.data = [(self.ring(c), q) for (c, q) in data]

        self.simplify()

    def simplify(self):
        res = []
        for i, (c1, q1) in enumerate(self.data):
            found = False
            for j, (c2, q2) in enumerate(res):
                iso, subst, orientation = q1.is_isomorphic(q2, certificate=True, with_multiple=True)
                if iso:
                    found = True
                    res[j] = c2 + c1 * orientation, q2
                    if res[j][0] == 0:
                        res.pop(j)
                    break
            if not found:
                res.append((c1, q1))
        self.data = res

    def __add__(self, other):
        return Chain(self.meta, self.data + other.data, self.ring)

    def __rmul__(self, x):
        return Chain(self.meta, [(c * x, q) for (c, q) in self.data], self.ring)

    def __mul__(self, x):
        return Chain(self.meta, [(c * x, q) for (c, q) in self.data], self.ring)

    def __sub__(self, other):
        return self + -1 * other

    def d(self):
        res = []
        for (c, q) in self.data:
            chain = q.d() * c
            res += chain.data
        return Chain(self.meta.meta_d, res, self.ring)

    def apply_homomorphism(self, h, new_ring=None):
        if new_ring is None:
            new_ring = self.ring
        return Chain(self.meta, [(h(c), q) for (c, q) in self.data], new_ring)

    def is_valid(self):
        if not all(isinstance(q, self.meta) for (_, q) in self.data):
            raise ValueError('Unknown objects in the chain!')
        if not all(c == self.ring(c) for (c, _) in self.data):
            raise ValueError('Unknown coefficients!')
        return True


# class SimplicialChain(Chain):
#     base_complex
#     def is_valid(self)
#     def d(self)
#
#
# class Subcomplex(SimplicialComplex):
#     base_complex
#     def barycentric_subdivision(self)
#     def is_valid(self)


# def to_base_cycles(chain_in_bistellar_moves: Chain)
#
# def link_coefficient(a: Subcomplex, b: Subcomplex)
