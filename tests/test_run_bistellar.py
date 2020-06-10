import os
import pytest
import sage.all as sg

from settings import gap_path, TMP_DIR
from spheres.gap_executor import gap_execute_commands, ListenerNLines
from spheres.simplicial_complex import Sphere, BistellarMove


@pytest.mark.gap
@pytest.mark.parametrize('commands, results', [
    (['1+1; 3+3;'], [['2\n6\n']]),
    (['2+2;',
      {'lines': ['Factorial(4); ', 'Sleep(2);', ' 2+2;'],
       'timeout': 3,
       'listener': ListenerNLines(3)},
      {'lines': '2+3;', 'listener': ListenerNLines(1)}], [['4\n'], ['24\n', '', '4\n'], ['5\n']]),
])
def test_gap_system(commands, results):
    if gap_path is None:
        return

    ans = gap_execute_commands(commands)
    for (res, status, _), true_res in zip(ans, results):
        assert status == 'ok'
        assert res == true_res


@pytest.mark.gap
@pytest.mark.parametrize('timeout, result', [(3, 'timeout'), (120, 'ok')])
def test_bistellar(timeout, result):
    """The test on BISTELLAR execution. 3 seconds is not enough, but 120 seconds is enough."""
    log_file, out_file, in_file = f'{TMP_DIR}/BISTELLAR.log', f'{TMP_DIR}/BISTELLAR.out', f'{TMP_DIR}/BISTELLAR.in'

    for f in (log_file, out_file, in_file):
        try:
            os.remove(f)
        except FileNotFoundError:
            pass

    bistellar_path = os.path.join(os.path.dirname(__file__), '../spheres/BISTELLAR.gap')

    status = gap_execute_commands([f'log_file := String("{log_file}");',
                                   f'out_file := String("{out_file}");',
                                   f'in_file := String("{in_file}");',
                                   f'type_of_object := 2;',
                                   f'facets := [];',
                                   f'polytope_dimension := 6;',
                                   f'number_of_vertices := 20;',
                                   (f'Read("{bistellar_path}");', timeout)])
    if result == 'timeout':
        assert status[-1] == ('', 'ok', {'timeout': True})
    elif result == 'ok':
        assert status[-1] == (['The examined complex is a sphere!!!\n'], 'ok', {})


@pytest.mark.gap
@pytest.mark.parametrize('s, timeout', [
    (BistellarMove(Sphere([[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]])
                   .join(sg.SimplicialComplex([[6], [7]]))
                   .as_sphere().rename_vertices(), [1, 2, 6]).t, 3)
])
def test_path_to_simplex(s, timeout):
    if gap_path is None:
        return

    q = s.path_to_simplex(timeout)
    if len(q):
        assert q[-1].t.is_isomorphic(Sphere(sg.Simplex(s.dimension() + 1).faces()))
    else:
        assert s.is_isomorphic(Sphere(sg.Simplex(s.dimension() + 1).faces()))[0]
