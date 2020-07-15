from spheres.simplicial_complex import *

ChainElement = sg.CombinatorialFreeModule.Element


def cochain_monomial_to_list(m):
    """Turn a SAGE monomial like \\chi_(1, 2, 6) into list [1, 2, 6] of vertices of a simplex."""
    return list(list(m)[0][0])


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

