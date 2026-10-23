import os
import sys
import numpy as np
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
python_dir = os.path.abspath(os.path.join(test_dir, '..', '..', 'python'))
sys.path.append(python_dir)
from OpenFUSIONToolkit.ThinCurr.junctions import (split_junctions, junction_table, find_junction_holes,
                                                  current_loop_count, has_junctions)
from OpenFUSIONToolkit.ThinCurr.meshing import write_ThinCurr_mesh, read_ThinCurr_mesh
from OpenFUSIONToolkit.ThinCurr._homology_cli import compute_homology
from OpenFUSIONToolkit._gmsh_cli import read_mesh


def sheet(n1, n2, mapping, region):
    u, v = np.meshgrid(np.linspace(0.0, 1.0, n1 + 1), np.linspace(0.0, 1.0, n2 + 1), indexing='ij')
    r = mapping(u.ravel(), v.ravel())
    idx = lambda i, j: i*(n2 + 1) + j
    lc = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)] for i in range(n1) for j in range(n2)]
    lc += [[idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)] for i in range(n1) for j in range(n2)]
    return r, np.array(lc), np.full(len(lc), region)


def merge(*sheets):
    r = np.vstack([s[0] for s in sheets])
    offsets = np.cumsum([0] + [len(s[0]) for s in sheets[:-1]])
    lc = np.vstack([s[1] + o for s, o in zip(sheets, offsets)])
    reg = np.hstack([s[2] for s in sheets])
    _, first, inv = np.unique(np.round(r, 9), axis=0, return_index=True, return_inverse=True)
    return r[first], inv.ravel()[lc].astype(np.int32), reg.astype(np.int32)


def tplate():
    a = sheet(4, 8, lambda u, v: np.c_[-0.5*u, v, 0*u], 1)
    b = sheet(4, 8, lambda u, v: np.c_[0.5*u, v, 0*u], 2)
    c = sheet(4, 8, lambda u, v: np.c_[0*u, v, 0.5*u], 3)
    return merge(a, b, c)


def ring():
    cyl = lambda z0: (lambda u, v: np.c_[0.5*np.cos(2*np.pi*v), 0.5*np.sin(2*np.pi*v), z0 + 0.5*u])
    up = sheet(4, 24, cyl(0.0), 1)
    low = sheet(4, 24, cyl(-0.5), 2)
    flange = sheet(3, 24, lambda u, v: np.c_[(0.5 + 0.3*u)*np.cos(2*np.pi*v), (0.5 + 0.3*u)*np.sin(2*np.pi*v), 0*u], 3)
    return merge(up, low, flange)


def torus(nsect=3):
    sectors = [sheet(4, 8, lambda u, v, k=k: np.c_[(1.0 + 0.3*np.cos(2*np.pi*v))*np.cos(2*np.pi*(k + u)/nsect),
                                                    (1.0 + 0.3*np.cos(2*np.pi*v))*np.sin(2*np.pi*(k + u)/nsect),
                                                    0.3*np.sin(2*np.pi*v)], k + 1) for k in range(nsect)]
    return merge(*sectors)


def test_tplate():
    r, lc, reg = tplate()
    assert has_junctions(lc)
    rs, lcs, lines = split_junctions(r, lc, reg)
    assert len(lines) == 1 and not lines[0]['closed'] and len(lines[0]['vertices']) == 9
    assert len(rs) == len(r) + 2*9
    kje, lje = junction_table(lines)
    assert len(kje) == 9 and len(lje) == 3*8
    out = find_junction_holes(r, lc, reg)
    assert out['n_unknowns'] == out['n_loops'] == current_loop_count(lc, reg)
    assert len(out['info']['holes']) == 0 and len(out['cuts']) == 0
    gap = find_junction_holes(r, lc, reg, insulated=(3,))
    assert gap['n_unknowns'] == out['n_unknowns'] - 7


@pytest.mark.parametrize("insulated,nholes", (((), 1), ((3,), 2)))
def test_ring(insulated, nholes):
    r, lc, reg = ring()
    out = find_junction_holes(r, lc, reg, insulated=insulated)
    assert out['lines'][0]['closed']
    assert out['n_unknowns'] == out['n_loops']
    assert len(out['info']['holes']) == nholes and len(out['cuts']) == 0


def test_torus_sectors():
    r, lc, reg = torus()
    assert not has_junctions(lc)
    out = find_junction_holes(r, lc, reg)
    assert out['n_unknowns'] == out['n_loops']
    assert len(out['cuts']) == 2 and len(out['info']['closures']) == 1 and len(out['holes']) == 2


def test_homology_file(tmp_path):
    r, lc, reg = ring()
    in_file = str(tmp_path / 'ring.h5')
    write_ThinCurr_mesh(in_file, r, lc, reg)
    compute_homology(in_file)
    mesh = read_ThinCurr_mesh(str(tmp_path / 'ring-homology.h5'))
    junc = mesh['thincurr']['junctions']
    assert mesh['r'].shape[0] == r.shape[0] + 2*24
    assert len(junc['holes']) == 1 and len(junc['kje']) == 25
    assert junc['P'].shape[1] == current_loop_count(lc, reg)
    assert junc['P'].shape[0] == np.count_nonzero(junc['pmap']) + len(junc['holes'])


def corner_map(lc, pmap, holes, nexp):
    nc = len(lc)
    I, J = [], []
    for c in range(nc):
        for i in range(3):
            if pmap[lc[c, i]] > 0:
                I.append(3*c + i)
                J.append(pmap[lc[c, i]] - 1)
    nact = np.count_nonzero(pmap)
    for h, hole in enumerate(holes):
        for c, i in hole:
            I.append(3*c + i)
            J.append(nact + h)
    import scipy.sparse
    return scipy.sparse.csr_matrix((np.ones(len(I)), (I, J)), shape=(3*nc, nexp))


@pytest.mark.parametrize("mesh", ('tplate', 'ring', 'torus'))
def test_expanded_map(mesh):
    r, lc, reg = {'tplate': tplate, 'ring': ring, 'torus': torus}[mesh]()
    out = find_junction_holes(r, lc, reg)
    B = corner_map(out['lc'], out['pmap'], out['holes'], out['P'].shape[0])
    assert abs(B @ out['P'] - out['M']).max() < 1.0E-12


def test_gmsh_physical_regions(tmp_path):
    msh = tmp_path / 'two.msh'
    msh.write_text('\n'.join(['$MeshFormat', '2.2 0 8', '$EndMeshFormat', '$Nodes', '4', '1 0 0 0', '2 1 0 0',
                              '3 1 1 0', '4 0 1 0', '$EndNodes', '$Elements', '2', '1 2 2 7 1 1 2 3',
                              '2 2 3 9 2 0 1 3 4', '$EndElements', '']))
    assert list(read_mesh(str(msh))[3]) == [1, 1]
    assert list(read_mesh(str(msh), physical_regions=True)[3]) == [7, 9]
