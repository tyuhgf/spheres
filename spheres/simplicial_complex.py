import json
import os
import re
from typing import List, Union, Set, Iterable, Tuple

from sage import all as sg

from settings import TMP_DIR
from spheres.gap_executor import gap_execute_commands


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
    if isinstance(self, sg.Simplex):
        s = sg.Simplex(s)
        if s.dimension() == 0:
            if hasattr(self, 'orientation'):
                s.orientation = self.orientation * -1
            else:
                s.orientation = -1
    return s


basic_is_isomorphic = sg.SimplicialComplex.is_isomorphic


def is_isomorphic(self, other, certificate=False):
    res = basic_is_isomorphic(self, other)
    if not certificate:
        return res, None
    else:
        if not res:
            return res, None
        return basic_is_isomorphic(self, other, certificate)


def decorate_init(base_init):
    def new_init(self, *args, **kwargs):
        base_init(self, *args, **kwargs)
        self.orientation = 1

    return new_init


sg.SimplicialComplex.is_homology_sphere = is_homology_sphere
sg.SimplicialComplex.rename_vertices = rename_vertices
sg.SimplicialComplex.as_sphere = as_sphere
sg.SimplicialComplex.is_isomorphic = is_isomorphic  # that is due to bug in sage

sg.Simplex.change_orientation = change_orientation
sg.Simplex.__init__ = decorate_init(sg.Simplex.__init__)  # now sg.Simplex has orientation property


class Sphere(sg.SimplicialComplex):
    def __init__(self, data: Union[Iterable[List], sg.SimplicialComplex], **kwargs):
        self.__class__.meta_d = self.__class__
        if isinstance(data, sg.SimplicialComplex):
            data = data.facets()

        super(Sphere, self).__init__(data, **kwargs)
        if not self.is_homology_sphere():
            raise ValueError('Complex is not a sphere!')
        self.set_immutable()

        if not self.facets() == set(sg.Simplex(s) for s in data):
            raise ValueError('Data is not list of facets!')

        self.facets_with_orientation = self.orient(list(data)[0])

    def rename_vertices(self, subst: Union[str, dict, callable] = 'int'):
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

    def is_minus_self(self):
        """Checks if a sphere has automorphism changing the orientation."""
        gr = self.automorphism_group()
        for gen in gr.gens():
            subst = gen.dict()
            if not self.check_oriented_facet([subst[v] for v in list(self.facets_with_orientation[0])]):
                return True
        return False

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
            if self.check_oriented_facet(sg.Simplex(s + v)):
                res.add(sg.Simplex(s))
            else:
                res.add(sg.Simplex(s).change_orientation())
        return res

    def barycentric_subdivision(self):
        """Makes a barycentric subdivided sphere with induced orientation."""
        if self.dimension() < 1:
            raise ValueError('For barycentric subdivision a sphere has to be of positive dimension.')
        res = super(Sphere, self).barycentric_subdivision()
        facets = list(res.facets())
        flag = sorted([list(v) for v in facets[0]], key=len)
        vertices = flag[0] + [(set(flag[i]) - set(flag[i - 1])).pop() for i in range(1, len(flag))]
        facets[0] = sg.Simplex([tuple(sorted(vertices[:i + 1])) for i in range(len(vertices))])
        if not self.check_oriented_facet(vertices):
            facets[0] = list(facets[0])
            facets[0][0], facets[0][1] = facets[0][1], facets[0][0]
            facets[0] = sg.Simplex(facets[0])
        return Sphere(facets)

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

    def path_to_simplex(self, timeout=3):
        """Returns list of BistellarMove objects representing the path to trivial sphere."""
        vertices_dict = {b: a+1 for (a, b) in enumerate(self.vertices())}
        vertices_dict_reverse = {v: k for (k, v) in vertices_dict.items()}
        s = self.rename_vertices(vertices_dict)

        facets = [list(f) for f in s.facets_with_orientation]

        log_file, out_file, in_file = f'{TMP_DIR}/BISTELLAR.log', f'{TMP_DIR}/BISTELLAR.out', f'{TMP_DIR}/BISTELLAR.in'

        for f in (log_file, out_file, in_file):
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

        bistellar_path = os.path.join(os.path.dirname(__file__), 'BISTELLAR.gap')
        try:
            status = gap_execute_commands([f'log_file := String("{log_file}");',
                                           f'out_file := String("{out_file}");',
                                           f'in_file := String("{in_file}");',
                                           f'type_of_object := 1;',
                                           f'facets := {facets};',
                                           (f'Read("{bistellar_path}");', timeout)])
        except Exception as e:
            raise e
        if not status[-1] == (['The examined complex is a sphere!!!\n'], 'ok', {}):
            raise Exception(f'Gap returned {status[-1]}')

        moves = json.load(open(log_file))
        moves = moves[:-1]  # the last line in log_file is []]
        for move in moves:
            for v in move[0] + move[1]:
                if v not in vertices_dict_reverse:
                    vertices_dict_reverse[v] = '_gap_' + str(v)
            move[0] = [vertices_dict_reverse[v] for v in move[0]]
            move[1] = [vertices_dict_reverse[v] for v in move[1]]

        bm = []
        sphere = self
        for move in moves:
            if len(move[1]) == 1:
                v = move[1][0]
            else:
                v = None
            bm.append(BistellarMove(sphere, move[0], new_vertex_name=v))
            sphere = bm[-1].t

        return bm


Sphere.meta_d = Sphere


class BistellarMove(sg.SimplicialComplex):
    def __init__(self, sphere: Sphere, sigma: List, new_vertex_name=None):
        s = sphere.__copy__().as_sphere().facets_with_orientation
        super(BistellarMove, self).__init__(s)

        if new_vertex_name is None:
            new_vertex_name = max([i for i in self.vertices() if isinstance(i, int)]) + 1

        self.dim = self.dimension()  # dimension of a sphere
        self.sigma = sigma

        v = list(self.star(sigma).vertices())
        if len(v) == self.dimension() + 2:
            pass
        elif len(v) == self.dimension() + 1 == len(sigma) and len(self.vertices()) > self.dimension() + 1:
            if new_vertex_name in self.vertices():
                raise ValueError('Vertex with this name is already in a complex.')
            v.append(new_vertex_name)
        else:
            raise ValueError('Cannot apply bistellar move!')

        self.tau = list(set(v) - set(sigma))

        self.add_face(v)

        t = set(s).copy()
        sim = sg.SimplicialComplex([v])
        sim.remove_face(v)
        facet = [facet for facet in s if facet in sim.facets()][0]
        sim = sim.as_sphere().orient(oriented_facet=facet)

        for facet in sim:
            if facet in t:
                t.remove(facet)
            else:
                t.update({facet.change_orientation()})

        self.s = Sphere(s)
        self.t = Sphere(t)

    def is_isomorphic(self, other, certificate=False):
        c_self = self.s.cone()\
            .rename_vertices({'R0': '__new_vertex__'})\
            .rename_vertices({'L' + str(v): v for v in self.vertices()})
        r_self = sg.SimplicialComplex(self.facets() + c_self.facets())

        c_other = other.s.cone()\
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
        return BistellarMove(self.t, self.tau, new_vertex_name=new_vertex_name)

    def link(self, simplex, is_mutable=True):
        if not set(simplex).issubset(set(self.sigma)):
            raise ValueError('Invalid simplex for link of bistellar move.')
        res = self.s.link_oriented(simplex)
        sigma = list(set(self.sigma) - set(simplex))
        new_vertex_name = self.tau[0] if len(self.tau) == 1 else None
        return BistellarMove(Sphere(res), sigma, new_vertex_name)

    def delta(self):
        return Chain(Sphere, [(1, self.t), (-1, self.s)])

    def skew_suspension(self, new_vertices_names=None) -> Sphere:
        if new_vertices_names is None:
            max_vertex = max([i for i in self.vertices() if isinstance(i, int)])
            new_vertices_names = max_vertex + 1, max_vertex + 2
        a, b = new_vertices_names
        if a in self.vertices() or b in self.vertices():
            raise ValueError('Vertex with this name is already in a complex.')

        s_facets = [change_orientation(list(f) + [a]) for f in self.s.facets_with_orientation]
        t_facets = [list(f) + [b] for f in self.t.facets_with_orientation]
        extra_facet = self.sigma + self.tau  # not oriented!
        return Sphere(s_facets + t_facets + [extra_facet])


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
                iso, orientation = q1.is_isomorphic(q2, certificate=False, with_multiple=True)
                if iso:
                    found = True
                    res[j] = c2 + c1 * orientation, q2
                    if res[j][0] == 0:
                        res.pop(j)
                    break
            if not found:
                res.append((c1, q1))
        self.data = res
        res = []
        for c, q in self.data:
            if q.is_minus_self():
                if self.ring.characteristic() == 2:
                    res.append((c, q))
                elif self.ring.characteristic() == 0:
                    if self.ring == sg.ZZ and c % 2:
                        res.append((c % 2, q))
            else:
                res.append((c, q))
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
