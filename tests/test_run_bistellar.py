import pytest
import sage.all as sg

from settings import gap_path
from spheres.gap_executor import gap_execute_commands, ListenerNLines
from spheres.simplicial_complex import Sphere, BistellarMove


@pytest.mark.parametrize('commands', [
    ['1+1; 3+3;'],
    ['2+2;',
     {'lines': ['Factorial(4); ', 'Sleep(2);', ' 2+2;'],
      'timeout': 3,
      'listener': ListenerNLines(2)},
     {'lines': '2+3;', 'listener': ListenerNLines(1)}],
])
def test_gap_system(commands):
    if gap_path is None:
        return

    gap_execute_commands(commands)  # todo add answer validation


def test_path_to_simplex():
    if gap_path is None:
        return

    s = Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]).join(sg.SimplicialComplex([[6], [7]])).as_sphere()
    s0 = s.rename_vertices()
    s1 = BistellarMove(s0, [1, 2, 6]).t

    q = s1.path_to_simplex()
    if len(q):
        assert q[-1].t.is_isomorphic(Sphere([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]))
    else:
        assert s1.is_isomorphic(Sphere([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]))


if __name__ == '__main__':
    test_path_to_simplex()
