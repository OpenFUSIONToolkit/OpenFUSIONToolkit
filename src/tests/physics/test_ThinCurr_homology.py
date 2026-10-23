#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Tests for ThinCurr hole and closure selection (OFT_ThinCurr_holes)
import os
import sys
import time
import multiprocessing
import pytest
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
python_dir = os.path.abspath(os.path.join(test_dir, '..','..','python'))
sys.path.append(python_dir)
from OpenFUSIONToolkit.ThinCurr.meshing import write_ThinCurr_mesh, read_ThinCurr_mesh
from OpenFUSIONToolkit.ThinCurr.homology import find_holes, count_holes
from OpenFUSIONToolkit.ThinCurr._homology_cli import compute_homology


def sheet(n1, n2, mapping):
    '''Structured triangular sheet over the unit square mapped by `mapping(u, v)`'''
    u, v = np.meshgrid(np.linspace(0.0, 1.0, n1+1), np.linspace(0.0, 1.0, n2+1), indexing='ij')
    r = mapping(u.ravel(), v.ravel())
    idx = lambda i, j: i*(n2+1) + j
    lc = [[idx(i,j), idx(i+1,j), idx(i+1,j+1)] for i in range(n1) for j in range(n2)]
    lc += [[idx(i,j), idx(i+1,j+1), idx(i,j+1)] for i in range(n1) for j in range(n2)]
    return r, np.array(lc)


def merge(*sheets):
    '''Merge sheets, joining coincident vertices and dropping unused ones'''
    r = np.vstack([s[0] for s in sheets])
    offsets = np.cumsum([0] + [s[0].shape[0] for s in sheets[:-1]])
    lc = np.vstack([s[1] + o for s, o in zip(sheets, offsets)])
    _, first, inv = np.unique(np.round(r, 9), axis=0, return_index=True, return_inverse=True)
    r, lc = r[first], inv.ravel()[lc]
    used, lc = np.unique(lc, return_inverse=True)
    return r[used], lc.reshape(-1,3).astype(np.int32)


def torus_map(R=1.0, a=0.3, center=(0.0,0.0,0.0)):
    return lambda u, v: np.c_[(R+a*np.cos(2*np.pi*v))*np.cos(2*np.pi*u),
                              (R+a*np.cos(2*np.pi*v))*np.sin(2*np.pi*u),
                              a*np.sin(2*np.pi*v)] + np.asarray(center)


def torus(n1=24, n2=16):
    return merge(sheet(n1, n2, torus_map()))


def torus_with_port():
    r, lc = torus(32, 24)
    cen = r[lc].mean(axis=1)
    keep = ~((np.abs(np.arctan2(cen[:,1], cen[:,0])) < 0.4) & (cen[:,2] > 0.15))
    return merge((r, lc[keep]))


def plate_cells(N, m, s):
    '''N x N grid of cells with m x m square openings of s x s cells'''
    keep = np.ones((N,N), dtype=bool)
    gap = (N-m*s)//(m+1)
    for i in range(m):
        for j in range(m):
            keep[gap+i*(s+gap):gap+i*(s+gap)+s, gap+j*(s+gap):gap+j*(s+gap)+s] = False
    return keep


def tri_plate(keep, z):
    N = keep.shape[0]
    idx = lambda i, j: i*(N+1) + j
    I, J = np.meshgrid(np.arange(N+1), np.arange(N+1), indexing='ij')
    r = np.c_[I.ravel()/N, J.ravel()/N, np.full(I.size, z)]
    lc = []
    for i, j in zip(*np.nonzero(keep)):
        lc += [[idx(i,j), idx(i+1,j), idx(i+1,j+1)], [idx(i,j), idx(i+1,j+1), idx(i,j+1)]]
    return r, np.array(lc)


def perforated(N=24, m=2, s=3):
    '''Plate with m*m openings (genus 0, m*m+1 boundary loops)'''
    return merge(tri_plate(plate_cells(N, m, s), 0.0))


def slab(N=24, m=2, s=3, h=0.1):
    '''Closed surface of genus m*m: two perforated plates joined by walls at every opening and the outside'''
    rb, lcb = tri_plate(plate_cells(N, m, s), 0.0)
    parts = [(rb, lcb), (rb+[0.0,0.0,h], lcb)]
    edges = {}
    for cell in lcb:
        for k in range(3):
            e = tuple(sorted((cell[k], cell[(k+1) % 3])))
            edges[e] = edges.get(e, 0) + 1
    for (p, q), count in edges.items():
        if count == 1:
            P, Q = rb[p], rb[q]
            parts.append(sheet(1, 2, lambda u, v, P=P, Q=Q: np.c_[P[0]+(Q[0]-P[0])*u, P[1]+(Q[1]-P[1])*u, h*v]))
    return merge(*parts)


def two_tori():
    return merge(sheet(16, 12, torus_map()), sheet(16, 12, torus_map(center=(3.0,0.0,0.0))))


def edge_set(lc):
    out = set()
    for cell in lc:
        for k in range(3):
            out.add(tuple(sorted((int(cell[k]), int(cell[(k+1) % 3])))))
    return out


def homology_rank(r, lc, loops):
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


def check_loop(lc, loop, edges, chord_free=True):
    '''Loop must be a simple chain of mesh edges, optionally without chords'''
    n = len(loop)
    assert len(set(int(v) for v in loop)) == n
    for k in range(n):
        assert tuple(sorted((int(loop[k]), int(loop[(k+1) % n])))) in edges
    if chord_free:
        for i in range(n):
            for j in range(i+2, n):
                if (i == 0) and (j == n-1):
                    continue
                assert tuple(sorted((int(loop[i]), int(loop[j])))) not in edges


meshes = {
    'torus': (torus, 2, 1),
    'torus_with_port': (torus_with_port, 2, 0),
    'perforated': (perforated, 4, 0),
    'slab': (slab, 8, 1),
    'two_tori': (two_tori, 4, 2),
}


@pytest.mark.parametrize("name", meshes.keys())
@pytest.mark.parametrize("optimize", (False, True))
def test_find_holes(name, optimize):
    builder, nholes, nclosures = meshes[name]
    r, lc = builder()
    assert count_holes(r, lc) == (nholes, nclosures)
    out = find_holes(r, lc, optimize=optimize)
    holes = out['boundary_holes'] + out['internal_holes']
    assert len(holes) == nholes
    assert len(out['closures']) == nclosures
    edges = edge_set(lc)
    for loop in out['boundary_holes']:
        check_loop(lc, loop, edges, chord_free=False)
    for loop in out['internal_holes']:
        check_loop(lc, loop, edges)
    assert homology_rank(r, lc, holes) == nholes


def loop_length(r, loops):
    total = 0.0
    for loop in loops:
        pts = r[np.asarray(loop)]
        total += np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1).sum()
    return total


def test_optimize():
    r, lc = slab()
    base = find_holes(r, lc)
    opt = find_holes(r, lc, optimize=True)
    assert loop_length(r, opt['internal_holes']) <= loop_length(r, base['internal_holes']) + 1.E-12


def test_invalid_meshes():
    # Three cells on one edge
    r = np.array([[0.0,0.0,0.0], [1.0,0.0,0.0], [0.5,1.0,0.0], [0.5,-1.0,0.0], [0.5,0.0,1.0]])
    lc = np.array([[0,1,2], [0,1,3], [0,1,4]])
    with pytest.raises(ValueError, match='manifold'):
        find_holes(r, lc)
    # Two cells touching at a vertex
    r = np.array([[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [-1.0,0.0,0.0], [0.0,-1.0,0.0]])
    lc = np.array([[0,1,2], [0,3,4]])
    with pytest.raises(ValueError, match='vertex'):
        find_holes(r, lc)
    # Mobius strip
    n = 24
    u, v = np.meshgrid(np.linspace(0.0, 2.0*np.pi, n, endpoint=False), np.r_[-0.2, 0.2], indexing='ij')
    r = np.c_[((1.0+v*np.cos(u/2))*np.cos(u)).ravel(), ((1.0+v*np.cos(u/2))*np.sin(u)).ravel(), (v*np.sin(u/2)).ravel()]
    idx = lambda i, j: (i % n)*2 + j
    lc = [[idx(i,0), idx(i+1,0), idx(i+1,1)] for i in range(n-1)] + [[idx(i,0), idx(i+1,1), idx(i,1)] for i in range(n-1)]
    lc += [[idx(n-1,0), idx(0,1), idx(0,0)], [idx(n-1,0), idx(0,0), idx(n-1,1)]] # Joined with a half twist
    with pytest.raises(ValueError, match='orientable'):
        find_holes(r, np.array(lc))


@pytest.mark.parametrize("use_legacy", (False, True))
def test_cli(tmp_path, use_legacy):
    r, lc = torus_with_port()
    in_file = str(tmp_path / 'port.h5')
    out_file = str(tmp_path / 'port-homology.h5')
    write_ThinCurr_mesh(in_file, r, lc, np.ones(lc.shape[0], dtype=np.int32))
    with pytest.raises(ValueError):
        compute_homology(in_file, verify_only=True, use_legacy=use_legacy)
    compute_homology(in_file, out_file=out_file, use_legacy=use_legacy)
    mesh = read_ThinCurr_mesh(out_file)
    assert len(mesh['thincurr']['holes']) == 2
    assert homology_rank(r, lc, mesh['thincurr']['holes']) == 2
    # Either method must accept the result
    compute_homology(out_file, verify_only=True)
    compute_homology(out_file, verify_only=True, use_legacy=True)


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


def test_eigs_match_legacy():
    '''Holes from either method must give the same current space, and so the same L/R eigenvalues'''
    r, lc = torus_with_port()
    write_ThinCurr_mesh('homology_port.h5', r, lc, np.ones(lc.shape[0], dtype=np.int32))
    eigs = {}
    for use_legacy in (False, True):
        out_file = 'homology_port-{0}.h5'.format('legacy' if use_legacy else 'new')
        compute_homology('homology_port.h5', out_file=out_file, use_legacy=use_legacy)
        eigs[use_legacy] = mp_run(run_eigs, (out_file,))
        assert eigs[use_legacy] is not None
    assert np.allclose(eigs[False], eigs[True], rtol=1.E-5, atol=0.0)
