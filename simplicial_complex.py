from typing import List, Union

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


def rename_vertices(self: sg.SimplicialComplex, names_dict):
    return sg.SimplicialComplex([[names_dict[v] if v in names_dict else v for v in f] for f in self.facets()])


def as_sphere(self: sg.SimplicialComplex):
    return Sphere(self.facets())


sg.SimplicialComplex.is_homology_sphere = is_homology_sphere
sg.SimplicialComplex.rename_vertices = rename_vertices
sg.SimplicialComplex.as_sphere = as_sphere


class Sphere(sg.SimplicialComplex):
    def __init__(self, data: Union[List[List], sg.SimplicialComplex], **kwargs):
        if isinstance(data, sg.SimplicialComplex):
            data = data.facets()

        super(Sphere, self).__init__(data, **kwargs)
        if not self.is_homology_sphere():
            raise ValueError('Complex is not a sphere!')
        self.set_immutable()

        if not self.facets() == set(sg.Simplex(s) for s in data):
            raise ValueError('Data is not list of facets!')

        self.facets_with_orientation = self.orient(data[0])

    def rename_vertices(self, names_dict):
        return Sphere(super(Sphere, self).rename_vertices(names_dict))

    def is_oriented(self):
        raise NotImplemented

    def orient(self, oriented_facet):
        if not isinstance(oriented_facet, sg.Simplex):
            oriented_facet = sg.Simplex(oriented_facet)
        g = self.flip_graph()
        edges = list(g.breadth_first_search(oriented_facet, edges=True))
        oriented_facets = [oriented_facet]
        for e in edges:
            oriented_facets.append(e[1])  # todo orient e[1]
        return oriented_facets

    def d(self):
        return [self.link_oriented(sg.Simplex([i])) for i in self.vertices()]  # todo make it a Chain object

    def link_oriented(self, v: Union[sg.Simplex, List]):
        """:return: set of simplices s of link(self, v) with orientation such that s+v is oriented simplex of self"""
        facets = list(self.link(v).facets())
        res = set()
        for s in facets:
            s = list(s)
            if len(s) <= 1 or self.check_oriented_facet(sg.Simplex(s + v)):
                res.add(sg.Simplex(s))
            else:
                res.add(sg.Simplex(s[:-2] + s[-2:][::-1]))
        return res

    def check_oriented_facet(self, f: sg.Simplex):
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
    def __init__(self, sphere, sigma, new_vertex_name_function=None):
        super(BistellarMove, self).__init__(sphere.copy())

        if new_vertex_name_function is None:
            def new_vertex_name_function():
                return max([i for i in self.vertices() if isinstance(i, int)]) + 1

        self.set_mutable()
        self.dim = self.dimension()  # dimension of a sphere
        self.sigma = sigma

        self.s = self.facets()

        v = list(self.star(sigma).vertices())
        if len(v) == self.dimension() + 1:
            self.add_face(v)
        elif len(v) == self.dimension() == len(sigma):
            self.add_face(v + [new_vertex_name_function()])
        else:
            raise ValueError('Cannot apply bistellar move!')

        self.add_face(v)

        self.t = self.s  # todo

    def is_isomorphic(self, other, certificate=False):
        return super(BistellarMove, self).is_isomorphic(other, certificate)

    # def isomorphic(self, bistellar_move)
    # def link(self, simplex)
    # def inverse(self)
    # def delta(self)
    # def is_valid(self)


class Chain:
    def __init__(self, meta, data: List, ring=sg.ZZ):
        self.meta = meta
        self.ring = ring
        self.data = data

        self.simplify()

    def simplify(self):
        for i, c1, q1 in enumerate(self.data):
            for c2, q2 in self.data[:i]:
                if q1.is_isomorphic(q2):
                    pass  # todo

    def __add__(self, other):
        return Chain(self.meta, self.data + other.data, self.ring)

    def __rmul__(self, x):
        return Chain(self.meta, [(c * x, q) for (c, q) in self.data], self.ring)

    def __sub__(self, other):
        return self + -1 * other

    def change_ring(self, ring):
        # todo
        raise NotImplemented

    # def is_valid(self)


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
