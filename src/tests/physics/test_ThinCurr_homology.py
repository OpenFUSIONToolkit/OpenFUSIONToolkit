#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Tests for ThinCurr hole and closure selection (OFT_ThinCurr_holes)
#
# Meshes are built with GMSH by "create_ThinCurr_meshes.py".
import os
import sys
import time
import multiprocessing
import pytest
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
python_dir = os.path.abspath(os.path.join(test_dir, '..','..','python'))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
sys.path.append(python_dir)
from oft_testing import run_command
from OpenFUSIONToolkit.ThinCurr.meshing import read_ThinCurr_mesh
from OpenFUSIONToolkit.ThinCurr._homology_cli import compute_homology

# Expected number of holes and closures for each mesh
meshes = {
    'homology_torus': (2, 1),
    'homology_torus_port': (2, 0),
    'homology_plate': (4, 0),
    'homology_genus4': (8, 1),
    'homology_two_tori': (4, 2),
    'homology_tongue': (1, 0),
}
# The legacy method caps each opening with a fan of cells from one boundary vertex, which fails when that vertex
# is joined to another vertex of the same boundary by a mesh edge (narrow tongue of material)
legacy_unsupported = ('homology_tongue',)

# Options for each method
methods = {
    'tree': [],
    'tree_opt': ['--optimize_holes'],
    'legacy': ['--use_legacy'],
    'legacy_opt': ['--use_legacy', '--optimize_holes'],
}


def run_script(command):
    '''Run a script, printing its output on failure'''
    outs, errs, errcode = run_command(" ".join(command))
    if errcode != 0:
        print("FAILED: {0} exited with non-zero error code!".format(command[0]))
        print("========== STD OUTPUT ==========")
        print(outs.decode())
        print("========== ERR OUTPUT ==========")
        print(errs.decode())
        print("ERRCODE = {0}".format(errcode))
        print("========== END OUTPUT ==========")
    return outs, errs, errcode


def convert_mesh(mesh_name):
    '''Convert a GMSH mesh to native format, returning the path of the converted mesh'''
    os.chdir(test_dir)
    _, _, errcode = run_script(["OFT_convert_gmsh.py", "--in_file={0}.msh".format(mesh_name)])
    assert errcode == 0
    return "{0}.h5".format(mesh_name)


def find_holes(mesh_name, method):
    '''Convert a GMSH mesh and add holes using `OFT_ThinCurr_holes.py`, returning the output file'''
    in_file = convert_mesh(mesh_name)
    out_file = "{0}-{1}.h5".format(mesh_name, method)
    _, _, errcode = run_script(["OFT_ThinCurr_holes.py", "--in_file={0}".format(in_file), "--out_file={0}".format(out_file)] + methods[method])
    assert errcode == 0
    return out_file


def edge_counts(lc):
    '''Number of cells on each edge'''
    out = {}
    for cell in lc:
        for k in range(3):
            e = tuple(sorted((int(cell[k]), int(cell[(k+1) % 3]))))
            out[e] = out.get(e, 0) + 1
    return out


def edge_set(lc):
    return set(edge_counts(lc))


def homology_rank(lc, loops):
    '''Rank of a set of loops in the first homology group (independent of the code under test)'''
    edges = sorted(edge_set(lc))
    eind = {e: i for i, e in enumerate(edges)}
    B = np.zeros((len(edges), lc.shape[0]))
    for f, cell in enumerate(lc):
        for k in range(3):
            a, b = int(cell[k]), int(cell[(k+1) % 3])
            B[eind[(min(a,b),max(a,b))], f] = 1.0 if a < b else -1.0
    Z = np.zeros((len(edges), len(loops)))
    for j, loop in enumerate(loops):
        for k in range(len(loop)):
            a, b = int(loop[k]), int(loop[(k+1) % len(loop)])
            Z[eind[(min(a,b),max(a,b))], j] += 1.0 if a < b else -1.0
    return np.linalg.matrix_rank(np.hstack((B, Z))) - np.linalg.matrix_rank(B)


def is_simple_loop(loop, edges):
    n = len(loop)
    if len(set(int(v) for v in loop)) != n:
        return False
    return all(tuple(sorted((int(loop[k]), int(loop[(k+1) % n])))) in edges for k in range(n))


def has_chord(loop, edges):
    n = len(loop)
    for i in range(n):
        for j in range(i+2, n):
            if (i == 0) and (j == n-1):
                continue
            if tuple(sorted((int(loop[i]), int(loop[j])))) in edges:
                return True
    return False


def loop_length(r, loops):
    total = 0.0
    for loop in loops:
        pts = r[np.asarray(loop)]
        total += np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1).sum()
    return total


@pytest.mark.parametrize("mesh_name,method", [(mesh_name, method) for mesh_name in meshes for method in methods
                                               if not (method.startswith('legacy') and mesh_name in legacy_unsupported)])
def test_holes(mesh_name, method):
    '''Hole and closure counts, loop validity and independence for each method'''
    nholes, nclosures = meshes[mesh_name]
    mesh = read_ThinCurr_mesh(find_holes(mesh_name, method))
    holes = mesh['thincurr'].get('holes', [])
    assert len(holes) == nholes
    assert len(mesh['thincurr'].get('closures', [])) == nclosures
    edges = edge_set(mesh['lc'])
    for loop in holes:
        assert is_simple_loop(loop, edges)
    assert homology_rank(mesh['lc'], holes) == nholes


@pytest.mark.parametrize("mesh_name", ('homology_torus_port', 'homology_genus4'))
def test_verify(mesh_name):
    '''Verification with either method accepts holes from the other method and rejects a mesh without holes'''
    for method in ('tree', 'legacy'):
        out_file = find_holes(mesh_name, method)
        compute_homology(out_file, verify_only=True)
        compute_homology(out_file, verify_only=True, use_legacy=True)
    for use_legacy in (False, True):
        with pytest.raises(ValueError):
            compute_homology(convert_mesh(mesh_name), verify_only=True, use_legacy=use_legacy)


def test_optimize():
    '''Optimization does not lengthen holes'''
    base = read_ThinCurr_mesh(find_holes('homology_genus4', 'tree'))
    opt = read_ThinCurr_mesh(find_holes('homology_genus4', 'tree_opt'))
    assert loop_length(opt['r'], opt['thincurr']['holes']) <= loop_length(base['r'], base['thincurr']['holes']) + 1.E-12


@pytest.mark.parametrize("mesh_name,message", (('homology_tee', 'manifold'),
                                               ('homology_touching', 'touch only at a vertex'),
                                               ('homology_mobius', 'orientable')))
def test_invalid(mesh_name, message):
    in_file = convert_mesh(mesh_name)
    outs, errs, errcode = run_command("OFT_ThinCurr_holes.py --in_file={0}".format(in_file))
    assert errcode != 0
    assert message in (outs + errs).decode()


def mp_run(target, args, timeout=180):
    if os.environ.get('OFT_DEBUG_TEST', 0):
        timeout *= 4
    os.chdir(test_dir)
    mp_q = multiprocessing.Queue()
    p = multiprocessing.Process(target=target, args=args + (mp_q,))
    p.start()
    start = time.time()
    while time.time() - start <= timeout:
        if not p.is_alive():
            break
        time.sleep(.5)
    else:
        print("Timeout reached")
        p.terminate()
        p.join()
        return None
    try:
        result = mp_q.get(timeout=5)
    except:
        print("Failed to get output")
        return None
    p.join()
    return result


def run_eigs(mesh_file, mp_q):
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        tw_model = ThinCurr(myOFT)
        tw_model.setup_model(mesh_file=mesh_file)
        tw_model.set_eta_values(eta_values=np.r_[1.E-6])
        tw_model.compute_Lmat()
        tw_model.compute_Rmat()
        eig_vals, _ = tw_model.get_eigs(10, direct=True)
        mp_q.put(np.sort(np.asarray(eig_vals)))
    except BaseException as e:
        print(e)
        mp_q.put(None)


def test_eigs_match():
    '''Holes from either method give the same current space, and so the same L/R eigenvalues'''
    eigs = {}
    for method in ('tree', 'legacy'):
        eigs[method] = mp_run(run_eigs, (find_holes('homology_torus_port', method),))
        assert eigs[method] is not None
    assert np.allclose(eigs['tree'], eigs['legacy'], rtol=1.E-5, atol=0.0)


def test_ordered_hole():
    '''A hole loop with chords is used in its given order (it cannot be rebuilt from its sorted vertices)'''
    out_file = find_holes('homology_tongue', 'tree')
    mesh = read_ThinCurr_mesh(out_file)
    loop = mesh['thincurr']['holes'][0]
    assert has_chord(loop, edge_set(mesh['lc']))
    assert mp_run(run_eigs, (out_file,)) is not None
