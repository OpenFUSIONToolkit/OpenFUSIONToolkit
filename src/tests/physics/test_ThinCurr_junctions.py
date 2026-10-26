import os
import sys
import copy
import time
import multiprocessing
import numpy as np
import scipy.sparse
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
    rs, lcs, jedges, lines, vmap = split_junctions(r, lc, reg)
    assert len(lines) == 1 and not lines[0]['closed'] and len(lines[0]['vertices']) == 9
    assert len(rs) == len(r) + 2*9
    kje, lje = junction_table(lc, lcs, jedges)
    assert len(kje) == 9 and len(lje) == 3*8
    out = find_junction_holes(r, lc, reg)
    assert out['n_unknowns'] == out['n_loops'] == current_loop_count(lc, reg)
    gap = find_junction_holes(r, lc, reg, insulated=(3,))
    assert gap['n_unknowns'] == out['n_unknowns'] - 7


@pytest.mark.parametrize("insulated,nholes", (((), 1), ((3,), 2)))
def test_ring(insulated, nholes):
    r, lc, reg = ring()
    out = find_junction_holes(r, lc, reg, insulated=insulated)
    assert out['lines'][0]['closed']
    assert out['n_unknowns'] == out['n_loops']


def test_torus_sectors():
    r, lc, reg = torus()
    assert not has_junctions(lc)
    out = find_junction_holes(r, lc, reg)
    assert out['n_unknowns'] == out['n_loops']


def test_homology_file(tmp_path):
    r, lc, reg = ring()
    in_file = str(tmp_path / 'ring.h5')
    write_ThinCurr_mesh(in_file, r, lc, reg)
    compute_homology(in_file)
    mesh = read_ThinCurr_mesh(str(tmp_path / 'ring-homology.h5'))
    junc = mesh['thincurr']['junctions']
    assert mesh['r'].shape[0] == r.shape[0] + 2*24
    assert len(junc['kje']) == 25
    assert junc['P'].shape[1] == current_loop_count(lc, reg)
    assert junc['P'].shape[0] == np.count_nonzero(junc['pmap']) + len(junc['holes'])


def corner_map(lc, pmap, holes, nexp):
    import scipy.sparse
    nc = len(lc)
    I, J, S = [], [], []
    for c in range(nc):
        for i in range(3):
            if pmap[lc[c, i]] > 0:
                I.append(3*c + i)
                J.append(pmap[lc[c, i]] - 1)
                S.append(1)
    nact = np.count_nonzero(pmap)
    for h, hole in enumerate(holes):
        for c, v, sgn in hole:
            I.append(3*c + list(lc[c]).index(v))
            J.append(nact + h)
            S.append(sgn)
    return scipy.sparse.csr_matrix((np.array(S, dtype=float), (I, J)), shape=(3*nc, nexp))


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


def test_insulate_file(tmp_path):
    from OpenFUSIONToolkit.ThinCurr.junctions import insulate_junction_file
    r, lc, reg = ring()
    in_file = str(tmp_path / 'ring.h5')
    write_ThinCurr_mesh(in_file, r, lc, reg)
    compute_homology(in_file)
    out = insulate_junction_file(str(tmp_path / 'ring-homology.h5'), (3,), str(tmp_path / 'ring-gap.h5'))
    junc = read_ThinCurr_mesh(str(tmp_path / 'ring-gap.h5'))['thincurr']['junctions']
    assert junc['P'].shape[1] == current_loop_count(lc, reg, insulated=(3,)) == out['n_unknowns']
    assert len(junc['holes']) == len(out['holes'])
    assert np.all(junc['vmap'] == out['vmap'])
    assert junc['insulated'] == [3]
    compute_homology(str(tmp_path / 'ring-gap.h5'), verify_only=True)
    import h5py
    with h5py.File(str(tmp_path / 'ring-gap.h5'), 'r+') as h5_file:
        del h5_file['thincurr/junctions/INSULATED']
    with pytest.raises(ValueError, match='unknowns'):
        compute_homology(str(tmp_path / 'ring-gap.h5'), verify_only=True)


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
        return False
    try:
        test_result = mp_q.get(timeout=5)
    except:
        print("Failed to get output")
        return False
    p.join()
    return test_result


def write_pair(r, lc, reg, prefix, insulated=()):
    """Junction mesh file with its map P and the same mesh with P = identity (expanded reference), with an I-coil."""
    from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_Icoil
    out = find_junction_holes(r, lc, reg, insulated=insulated)
    ident = copy.copy(out)
    ident['P'] = scipy.sparse.identity(out['P'].shape[0], format='csr')
    coil = ThinCurr_Icoil('IC')
    coil.add_circular_subcoil([0.7, 0.3])
    files = []
    for tag, junc in (('junction', out), ('expanded', ident)):
        path = os.path.join(test_dir, '{0}_{1}.h5'.format(prefix, tag))
        write_ThinCurr_mesh(path, out['r'], out['lc'], reg, junctions=junc, eta_surf=1.0E-4*np.ones(reg.max()),
                            coil_sets=[coil])
        files.append(path)
    return files, out['P'].toarray()


def run_reduction(files, P, insulated_runtime, mp_q):
    """Junction model matrices must equal P applied to the expanded model's."""
    result = True
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        mats = []
        for path, ins in ((files[0], insulated_runtime), (files[1], None)):
            tw = ThinCurr(myOFT)
            # The expanded reference (P = identity) is not a valid reduced model, so it is not verified
            tw.setup_model(mesh_file=path, junction_insulated=ins, no_verify_holes=(path == files[1]))
            Mc = np.array(tw.compute_Mcoil())
            tw.compute_Lmat()
            tw.compute_Rmat()
            mats.append((np.array(tw.Lmat), tw.Rmat.toarray(), Mc, np.sort(tw.get_eigs(4, direct=True)[0])))
        (L, R, Mc, eig), (Lx, Rx, Mcx, _) = mats
        if P is not None:
            rel = lambda a, b: np.abs(a - b).max()/np.abs(b).max()
            errs = [rel(L, P.T @ Lx @ P), rel(R, P.T @ Rx @ P), rel(Mc, Mcx @ P)]
            print('L, R, Mcoil errors', errs)
            result = max(errs) < 1.0E-12
        else:
            result = eig
    except BaseException as e:
        print(e)
        result = False
    mp_q.put(result)


@pytest.mark.parametrize("mesh", ('tplate', 'ring', 'torus'))
def test_reduction(mesh):
    r, lc, reg = {'tplate': tplate, 'ring': ring, 'torus': torus}[mesh]()
    files, P = write_pair(r, lc, reg, 'junction_' + mesh)
    assert mp_run(run_reduction, (files, P, None))


def run_cache(path, mp_q):
    """Matrices saved to and loaded from cache files must equal those computed directly."""
    result = False
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        from OpenFUSIONToolkit.ThinCurr.sensor import Mirnov, save_sensors
        save_sensors([Mirnov([0.2, 0.5, 0.1], [0.0, 0.0, 1.0], 'B1'), Mirnov([0.1, 0.3, 0.2], [1.0, 0.0, 0.0], 'B2')],
                     'junction_cache.loc')
        myOFT = OFT_env(nthreads=-1)
        mats = []
        for suffix in (None, '.save', '.save'):
            tw = ThinCurr(myOFT)
            tw.setup_model(mesh_file=path)
            tw.compute_Lmat(cache_file=suffix and 'junction_cache_L' + suffix)
            Ms, Msc, _ = tw.compute_Msensor('junction_cache.loc', cache_file=suffix and 'junction_cache_Ms' + suffix)
            mats.append((np.array(tw.Lmat), np.array(Ms), np.array(Msc)))
        result = all(a.shape == b.shape and np.allclose(a, b, rtol=0.0, atol=1.E-12*np.abs(a).max())
                     for m in mats[1:] for a, b in zip(mats[0], m))
    except BaseException as e:
        print(e)
    mp_q.put(result)


def test_cache():
    for f in ('junction_cache_L.save', 'junction_cache_Ms.save'):
        if os.path.exists(os.path.join(test_dir, f)):
            os.remove(os.path.join(test_dir, f))
    files, _ = write_pair(*tplate(), 'junction_cache')
    assert mp_run(run_cache, (files[0],))


def test_runtime_insulation():
    r, lc, reg = ring()
    files_rt, _ = write_pair(r, lc, reg, 'junction_ring_rt')
    files_ins, _ = write_pair(r, lc, reg, 'junction_ring_ins', insulated=(3,))
    eig_rt = mp_run(run_reduction, (files_rt, None, [3]))
    eig_ins = mp_run(run_reduction, (files_ins, None, None))
    assert eig_rt is not False and eig_ins is not False
    assert np.allclose(eig_rt, eig_ins, rtol=1.0E-10, atol=0.0)


@pytest.mark.parametrize("insulated,nholes", (([1], 2), ([(1, 2)], 1), ([(1, 2), (2, 3)], 2), ([1, 2, 3], 3)))
def test_torus_insulation(insulated, nholes):
    r, lc, reg = torus()
    out = find_junction_holes(r, lc, reg, insulated=insulated)
    assert out['n_unknowns'] == out['n_loops'] == current_loop_count(lc, reg, insulated=insulated)


def test_split_regions_file(tmp_path):
    r, lc, reg = torus()
    in_file = str(tmp_path / 'torus.h5')
    write_ThinCurr_mesh(in_file, r, lc, reg)
    compute_homology(in_file, split_regions=True)
    junc = read_ThinCurr_mesh(str(tmp_path / 'torus-homology.h5'))['thincurr']['junctions']
    assert junc['P'].shape[1] == current_loop_count(lc, reg)
    assert len(junc['kje']) - 1 == 3*8


def test_runtime_pair_insulation():
    r, lc, reg = torus()
    files_rt, _ = write_pair(r, lc, reg, 'junction_torus_rt')
    files_ins, _ = write_pair(r, lc, reg, 'junction_torus_ins', insulated=[(1, 2)])
    eig_rt = mp_run(run_reduction, (files_rt, None, [(1, 2)]))
    eig_ins = mp_run(run_reduction, (files_ins, None, None))
    assert eig_rt is not False and eig_ins is not False
    assert np.allclose(eig_rt, eig_ins, rtol=1.0E-10, atol=0.0)


def test_region_both_sides():
    r, lc, reg = tplate()
    reg = np.where(reg == 2, 1, reg)
    out = find_junction_holes(r, lc, reg)
    assert out['n_unknowns'] == out['n_loops']


def test_point_contact():
    a = sheet(4, 4, lambda u, v: np.c_[-0.5 + 0.5*u, -0.5 + 0.5*v, 0*u], 1)
    b = sheet(4, 4, lambda u, v: np.c_[0.5*u, 0.5*v, 0.3*u*v], 2)
    r, lc, reg = merge(a, b)
    assert has_junctions(lc, reg)
    out = find_junction_holes(r, lc, reg)
    assert out['n_unknowns'] == out['n_loops']
    assert len(out['r']) == len(r) + 1


def test_pair_at_three_sheets_warns():
    r, lc, reg = tplate()
    with pytest.warns(UserWarning):
        out = find_junction_holes(r, lc, reg, insulated=[(1, 2)])
    assert out['n_unknowns'] == find_junction_holes(r, lc, reg)['n_unknowns']


def test_torus_with_rib():
    rib = sheet(2, 8, lambda u, v: np.c_[1.0 + (0.3 + 0.15*u)*np.cos(2*np.pi*v), 0*u, (0.3 + 0.15*u)*np.sin(2*np.pi*v)], 4)
    r, lc, reg = merge(*[torus_sector_sheet(k) for k in range(3)], rib)
    assert has_junctions(lc, reg)
    out = find_junction_holes(r, lc, reg)
    assert out['n_unknowns'] == out['n_loops']


def torus_sector_sheet(k, nsect=3):
    return sheet(4, 8, lambda u, v: np.c_[(1.0 + 0.3*np.cos(2*np.pi*v))*np.cos(2*np.pi*(k + u)/nsect),
                                           (1.0 + 0.3*np.cos(2*np.pi*v))*np.sin(2*np.pi*(k + u)/nsect),
                                           0.3*np.sin(2*np.pi*v)], k + 1)


def periodic_torus_file(path, nfp=3, nphi=7, ntheta=12):
    from OpenFUSIONToolkit.ThinCurr.meshing import ThinCurr_periodic_toroid
    import h5py
    phi = np.linspace(0.0, 2.0*np.pi/nfp, nphi)
    th = np.linspace(0.0, 2.0*np.pi, ntheta, endpoint=False)
    P, T = np.meshgrid(phi, th, indexing='ij')
    r_grid = np.stack(((1.0 + 0.3*np.cos(T))*np.cos(P), (1.0 + 0.3*np.cos(T))*np.sin(P), 0.3*np.sin(T)), axis=-1)
    toroid = ThinCurr_periodic_toroid(r_grid, nfp, ntheta, nphi)
    centers = toroid.r[toroid.lc].mean(axis=1)
    reg = np.where(np.mod(np.arctan2(centers[:, 1], centers[:, 0]), 2.0*np.pi/nfp) < np.pi/nfp, 1, 2).astype(np.int32)
    toroid.write_to_file(path, reg=reg)
    with h5py.File(path, 'r+') as h5:
        h5.create_dataset('thincurr/ETA_SURF', data=1.0E-4*np.ones(2))


def test_periodic_count(tmp_path):
    path = str(tmp_path / 'ptorus.h5')
    periodic_torus_file(path)
    info = read_ThinCurr_mesh(path)['thincurr']
    compute_homology(path, split_regions=True)
    junc = read_ThinCurr_mesh(str(tmp_path / 'ptorus-homology.h5'))['thincurr']['junctions']
    assert junc['P'].shape[1] == (info['pmap'].max() + 1) + len(info['holes'])
    assert 'periodic' in junc


def run_periodic(files, mp_q):
    result = False
    try:
        import scipy.linalg
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        eigs = []
        for path in files:
            tw = ThinCurr(myOFT)
            tw.setup_model(mesh_file=path)
            tw.compute_Lmat()
            tw.compute_Rmat()
            eigs.append(np.sort(scipy.linalg.eigh(np.array(tw.Lmat), tw.Rmat.toarray(), eigvals_only=True))[::-1])
        print('periodic eigenvalue difference', np.abs(eigs[0] - eigs[1]).max()/eigs[1].max())
        result = len(eigs[0]) == len(eigs[1]) and np.abs(eigs[0] - eigs[1]).max() < 1.0E-5*eigs[1].max()
    except BaseException as e:
        print(e)
    mp_q.put(result)


def test_periodic_vs_thincurr():
    path = os.path.join(test_dir, 'junction_ptorus.h5')
    periodic_torus_file(path)
    compute_homology(path, split_regions=True)
    assert mp_run(run_periodic, ([os.path.join(test_dir, 'junction_ptorus-homology.h5'), path],))
