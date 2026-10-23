#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Mesh tools for ThinCurr models with junctions (edges shared by more than two sheets)

Each region is a sheet. The mesh is split along region interfaces so every sheet gets its own copy of
the interface vertices and the split mesh is manifold. Current conservation across each junction edge
(Kirchhoff) then joins the sheets again:

  sum_k s_k chi_k = 0 at each end of each junction edge,

over the sheets on that edge, where chi_k is the stream function seen by sheet k's cell on the edge and
s_k = +1 if that cell lists the edge in the junction's direction, otherwise -1.

Holes are found automatically:
  - outer arcs: connected chains of outer edges (boundary edges that are not junction edges); chi is
    uniform along each arc and Kirchhoff links arc values at the ends of junction lines;
  - cuts: a jump h of chi across a closed junction line, or across a loop of the homology basis of a
    manifold mesh, stored as the (cell, corner) pairs on one side of the cut. Cuts may cross junctions;
  - references: arc values (or one vertex, a closure, on a closed surface) fixed at zero to remove
    constant shifts that carry no current.
A cut is kept only if it adds an unknown without a zero-energy mode. The result is checked against the
number of independent current patterns on the original mesh.

@authors Ian Stewart
@date October 2026
@ingroup doxy_oft_python
'''
import os
import io
import tempfile
import contextlib
from fractions import Fraction
import numpy
import scipy.sparse
import scipy.sparse.csgraph
import scipy.linalg
import scipy.sparse.linalg


def _edge_cells(lc):
    '''! Cells on each edge, {(a, b) with a < b: [cells]}'''
    out = {}
    for c, tri in enumerate(lc):
        for i in range(3):
            a, b = int(tri[i]), int(tri[(i + 1) % 3])
            out.setdefault((min(a, b), max(a, b)), []).append(c)
    return out


def _edge_sign(tri, a, b):
    t = list(tri)
    return 1 if t[(t.index(a) + 1) % 3] == b else -1


def orient_cells(lc):
    '''! Cell orientation as ThinCurr sets it on load (`sync_face_normals`)

    In each connected piece the lowest-numbered cell keeps its order and its neighbours are flipped to
    agree with it. Requires a manifold mesh.

    @param lc Cell list [nc,3] (0-based)
    @result Oriented cell list [nc,3]
    '''
    lc = numpy.array(lc).copy()
    ec = _edge_cells(lc)
    if max(len(v) for v in ec.values()) > 2:
        raise ValueError('orient_cells requires a manifold mesh')
    done = numpy.zeros(len(lc), dtype=bool)
    for seed in range(len(lc)):
        if done[seed]:
            continue
        done[seed] = True
        stack = [seed]
        while stack:
            c = stack.pop()
            tri = lc[c]
            for i in range(3):
                a, b = tri[i], tri[(i + 1) % 3]
                for n in ec[(min(a, b), max(a, b))]:
                    if n == c or done[n]:
                        continue
                    t = list(lc[n])
                    if t[(t.index(a) + 1) % 3] == b:
                        lc[n] = lc[n][::-1]
                    done[n] = True
                    stack.append(n)
    return lc


def _pairs(line):
    n = len(line['vertices'])
    return [(i, i + 1) for i in range(n - 1)] + ([(n - 1, 0)] if line['closed'] else [])


def has_junctions(lc):
    '''! True if any edge is shared by more than two cells'''
    for cells in _edge_cells(lc).values():
        if len(cells) > 2:
            return True
    return False


def split_junctions(r, lc, reg):
    '''! Split a mesh along region interfaces, giving each region its own copy of those vertices

    Junction edges are grouped into lines that run between branch points (vertices where the number of
    junction edges is not two) or close on themselves. A vertex shared by several lines has one copy per
    sheet. Each sheet must touch a junction edge from one side only. Cells keep their index.

    @param r Point list [np,3]
    @param lc Cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @result r_split, lc_split, lines: one dict per junction line with
            'vertices' original vertex indices along the line,
            'closed'   True if the line is a loop,
            'segments' per edge of the line, the set of regions on it,
            'copies'   per vertex, {region: split index} for the regions on the adjacent edges,
            'cells'    per edge of the line, {region: cell on that edge}
    '''
    lc = numpy.asarray(lc)
    reg = numpy.asarray(reg)
    ec = _edge_cells(lc)
    jedges = {}
    for e, cells in ec.items():
        regs = [int(k) for k in reg[cells]]
        if len(set(regs)) > 1:
            if len(regs) != len(set(regs)):
                raise ValueError('A sheet lies on both sides of junction edge {0}'.format(e))
            jedges[e] = dict(zip(regs, cells))
        elif len(cells) > 2:
            raise ValueError('Edge {0} is shared by {1} cells of one region'.format(e, len(cells)))
    adj = {}
    for a, b in jedges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    vreg = {}
    for tri, k in zip(lc, reg):
        for v in tri:
            vreg.setdefault(int(v), set()).add(int(k))
    for v, regions in vreg.items():
        if v not in adj and len(regions) > 1:
            raise ValueError('Sheets touch only at vertex {0}'.format(v))
    ends = sorted(v for v, n in adj.items() if len(n) != 2)
    key = lambda a, b: (min(a, b), max(a, b))
    paths, used = [], set()
    for start in ends + sorted(adj):
        for w in sorted(adj[start]):
            if key(start, w) in used:
                continue
            used.add(key(start, w))
            path = [start, w]
            while path[-1] not in ends and path[-1] != start:
                v = path[-1]
                nxt = next(x for x in adj[v] if key(v, x) not in used)
                used.add(key(v, nxt))
                path.append(nxt)
            closed = path[-1] == start and start not in ends
            if closed:
                path.pop()
            paths.append((path, closed))
    r_split, lc_split, copy = [p for p in r], lc.copy(), {}
    for v in sorted(adj):
        regions = sorted(vreg[v])
        copy[v] = {regions[0]: int(v)}
        for k in regions[1:]:
            r_split.append(numpy.array(r[v]).copy())
            copy[v][k] = len(r_split) - 1
            cells = (reg == k) & numpy.any(lc == v, axis=1)
            lc_split[cells] = numpy.where(lc[cells] == v, copy[v][k], lc_split[cells])
    lines = []
    for path, closed in paths:
        line = {'vertices': [int(v) for v in path], 'closed': closed}
        pairs = _pairs(line)
        line['segments'] = [set(jedges[key(path[i], path[j])]) for i, j in pairs]
        line['cells'] = [jedges[key(path[i], path[j])] for i, j in pairs]
        near = [set() for _ in path]
        for p, (i, j) in enumerate(pairs):
            near[i] |= line['segments'][p]
            near[j] |= line['segments'][p]
        line['copies'] = [{k: copy[v][k] for k in sorted(s)} for v, s in zip(path, near)]
        lines.append(line)
    return numpy.array(r_split), lc_split, lines


def junction_table(lines):
    '''! Junction edge table for the mesh file

    @result kje [n_jedges+1] pointer, lje [n,3] (vertex a, vertex b, cell) per sheet on each junction
    edge, 0-based split indices; a -> b is the junction's direction
    '''
    kje, lje = [0], []
    for line in lines:
        for (i, j), cells in zip(_pairs(line), line['cells']):
            for k in sorted(cells):
                lje.append([line['copies'][i][k], line['copies'][j][k], cells[k]])
            kje.append(len(lje))
    return numpy.array(kje, dtype=numpy.int32), numpy.array(lje, dtype=numpy.int32).reshape(-1, 3)


def current_loop_count(lc, reg, insulated=()):
    '''! Number of independent divergence-free current patterns on the original mesh

    Cycles of the dual graph whose nodes are the cells and the junction edges (Kirchhoff), with no flow
    across outer edges or into insulated sheets.
    '''
    reg = numpy.asarray(reg)
    nf = len(lc)
    I, J, nj = [], [], 0
    for e, cells in _edge_cells(lc).items():
        regs = set(int(k) for k in reg[cells])
        if len(regs) == 1:
            if len(cells) == 2:
                I.append(cells[0])
                J.append(cells[1])
            continue
        joined = [c for c in cells if int(reg[c]) not in insulated]
        if len(joined) >= 2:
            for c in joined:
                I.append(c)
                J.append(nf + nj)
        nj += 1
    n = nf + nj
    G = scipy.sparse.coo_matrix((numpy.ones(len(I)), (I, J)), shape=(n, n))
    ncomp = scipy.sparse.csgraph.connected_components(G, directed=False)[0]
    return len(I) - n + ncomp


def _solve(rows):
    '''! Exact sparse elimination of rows (each {var: coef}) = 0, pivoting on the largest variable
    @result {pivot var: {free var: coef}}'''
    expr, users = {}, {}
    for row in rows:
        r = {}
        for v, c in row.items():
            for f, cf in (expr[v].items() if v in expr else ((v, Fraction(1)),)):
                r[f] = r.get(f, 0) + Fraction(c)*cf
        r = {v: c for v, c in r.items() if c != 0}
        if not r:
            continue
        p = max(r)
        cp = r.pop(p)
        e = {v: -c/cp for v, c in r.items()}
        for q in users.pop(p, set()):
            cq = expr[q].pop(p)
            for f, cf in e.items():
                val = expr[q].get(f, 0) + cq*cf
                if val == 0:
                    expr[q].pop(f, None)
                    users.get(f, set()).discard(q)
                else:
                    expr[q][f] = val
                    users.setdefault(f, set()).add(q)
        expr[p] = e
        for f in e:
            users.setdefault(f, set()).add(p)
    return expr


def _ordered_chain(vertices, edges):
    '''! Order a connected chain of edges; @result (points, closed)'''
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    if max(len(n) for n in adj.values()) > 2:
        raise ValueError('Outer arc through vertex {0} branches'.format(
            next(v for v, n in adj.items() if len(n) > 2)))
    ends = [v for v, n in adj.items() if len(n) == 1]
    start = ends[0] if ends else min(adj)
    chain, prev = [start], None
    while True:
        nxt = [w for w in adj[chain[-1]] if w != prev and (len(chain) < 2 or w != chain[-2])]
        if not nxt or nxt[0] == start:
            break
        prev = chain[-1]
        chain.append(nxt[0])
    return numpy.array(chain, dtype=numpy.int32), not ends


def _local_stiffness(r, lc):
    '''! P1 stiffness of each cell, int grad(lambda_i).grad(lambda_j) dA, [nc,3,3]'''
    p = r[lc]
    e = numpy.stack([p[:, 2] - p[:, 1], p[:, 0] - p[:, 2], p[:, 1] - p[:, 0]], axis=1)
    area = 0.5*numpy.linalg.norm(numpy.cross(e[:, 0], e[:, 1]), axis=1)
    return numpy.einsum('cik,cjk->cij', e, e)/(4.0*area)[:, None, None]


def loop_cut(lc, loop):
    '''! Cut across a closed vertex loop of a manifold mesh: the (cell, corner) pairs on its left side

    @param lc Cell list [nc,3] (0-based, manifold)
    @param loop Vertex indices around the loop
    @result [n,2] (cell, corner) pairs, 0-based
    '''
    lco = orient_cells(lc)
    ec = _edge_cells(lco)
    n = len(loop)
    out = []
    for p, v in enumerate(loop):
        v, prev, nxt = int(v), int(loop[p - 1]), int(loop[(p + 1) % n])
        fan = set(numpy.where(numpy.any(lco == v, axis=1))[0])
        left = next(c for c in ec[(min(v, nxt), max(v, nxt))] if _edge_sign(lco[c], v, nxt) == 1)
        side, stack = {left}, [left]
        while stack:
            c = stack.pop()
            for w in lco[c]:
                if w in (v, prev, nxt):
                    continue
                for c2 in ec[(min(v, w), max(v, w))]:
                    if c2 in fan and c2 not in side:
                        side.add(c2)
                        stack.append(c2)
        right = next(c for c in ec[(min(prev, v), max(prev, v))] if _edge_sign(lco[c], prev, v) == 1)
        if right not in side:
            raise ValueError('Loop turns back on itself at vertex {0}'.format(v))
        out += [(c, list(lc[c]).index(v)) for c in sorted(side)]
    return numpy.array(out, dtype=numpy.int32)


def junction_basis(r_split, lc_split, lines, insulated=(), cuts=()):
    '''! Reduced unknowns of a split mesh: map from unknowns to the stream function at each cell corner

    @param r_split, lc_split, lines From @ref split_junctions
    @param insulated Regions with chi = 0 at every junction (left out of Kirchhoff)
    @param cuts List of [n,2] (cell, corner) arrays, each a jump of chi with its own unknown
    @result M (sparse [3 nc, n_reduced], corner (c, i) is row 3c+i), info dict with
            'arcs' [(points, closed)], 'references' (arc indices fixed at zero), 'holes' (free arc indices),
            'cuts' (indices of independent cuts), 'closures' (cells), 'shifts' (zero-current shifts removed)
    '''
    nv, nc = len(r_split), len(lc_split)
    lco = orient_cells(lc_split)
    corner_cut = {}
    for j, cut in enumerate(cuts):
        for c, i in cut:
            corner_cut.setdefault((int(c), int(i)), []).append(-1 - j)
    joined_edges, rows, crows, jverts = set(), [], [], set()
    for line in lines:
        signs = {}
        for (i, j), cells in zip(_pairs(line), line['cells']):
            for k, c in cells.items():
                a, b = line['copies'][i][k], line['copies'][j][k]
                s = _edge_sign(lco[c], a, b)
                if signs.setdefault(k, s) != s:
                    raise ValueError('Sheet {0} changes orientation along a junction line'.format(k))
                jverts |= {a, b}
        for (i, j), regions, cells in zip(_pairs(line), line['segments'], line['cells']):
            joined = [k for k in regions if k not in insulated]
            if len(joined) < 2:
                continue
            for k in joined:
                a, b = line['copies'][i][k], line['copies'][j][k]
                joined_edges.add((min(a, b), max(a, b)))
            for idx in (i, j):
                row, crow = {}, {}
                for k in joined:
                    v = line['copies'][idx][k]
                    row[v] = signs[k]
                    for h in corner_cut.get((cells[k], list(lc_split[cells[k]]).index(v)), []):
                        crow[h] = crow.get(h, 0) + signs[k]
                rows.append(row)
                crows.append(crow)
    ecs = _edge_cells(lc_split)
    outer = [e for e, cells in ecs.items() if len(cells) == 1 and e not in joined_edges]
    var_of = numpy.arange(nv)
    arcs = []
    if outer:
        G = scipy.sparse.coo_matrix((numpy.ones(len(outer)), ([e[0] for e in outer], [e[1] for e in outer])),
                                    shape=(nv, nv))
        lab = scipy.sparse.csgraph.connected_components(G, directed=False)[1]
        groups = {}
        for e in outer:
            groups.setdefault(lab[e[0]], []).append(e)
        for g in sorted(groups):
            pts = sorted(set(v for e in groups[g] for v in e))
            var_of[pts] = nv + len(arcs)
            arcs.append(_ordered_chain(pts, groups[g]))
    rows = [{int(var_of[v]): c for v, c in row.items()} for row in rows]
    expr = _solve([{**row, **crow} for row, crow in zip(rows, crows)])
    G = scipy.sparse.coo_matrix((numpy.ones(3*nc), (lc_split.ravel(), numpy.repeat(lc_split[:, 0], 3))),
                                shape=(nv, nv))
    comp = scipy.sparse.csgraph.connected_components(G, directed=False)[1]
    comp_of = lambda x: comp[x] if x < nv else comp[arcs[x - nv][0][0]]
    ncomp = comp.max() + 1
    gexpr = _solve([{int(comp_of(x)): c for x, c in row.items()} for row in rows])
    gfree = [c for c in range(ncomp) if c not in gexpr]
    gauge = numpy.array([[1.0 if c == f else float(gexpr.get(c, {}).get(f, 0)) for c in range(ncomp)]
                         for f in gfree]).reshape(len(gfree), ncomp)
    free = [x for x in sorted(set(var_of.tolist())) if x not in expr]
    plain = set(numpy.unique(lc_split)) - jverts
    order = [x for x in free if x >= nv] + [x for x in free if x < nv and x in plain]
    fixed, block = [], numpy.zeros((len(gfree), 0))
    for x in order:
        if len(fixed) == len(gfree):
            break
        trial = numpy.hstack((block, gauge[:, [comp_of(x)]]))
        if numpy.linalg.matrix_rank(trial) > block.shape[1]:
            block, fixed = trial, fixed + [x]
    if len(fixed) != len(gfree):
        raise ValueError('Could not remove every zero-current shift')
    hfree = [-1 - j for j in range(len(cuts)) if -1 - j not in expr]
    cols = [x for x in free if x not in fixed] + hfree
    col = {x: j for j, x in enumerate(cols)}
    terms_of = lambda x: expr.get(x, {x: 1})
    I, J, V = [], [], []
    for c in range(nc):
        for i in range(3):
            terms = dict(terms_of(int(var_of[lc_split[c, i]])))
            for h in corner_cut.get((c, i), []):
                for f, cf in terms_of(h).items():
                    terms[f] = terms.get(f, 0) + cf
            for f, cf in terms.items():
                if f in col and cf != 0:
                    I.append(3*c + i)
                    J.append(col[f])
                    V.append(float(cf))
    M = scipy.sparse.csr_matrix((V, (I, J)), shape=(3*nc, len(cols)))
    closures = []
    for x in fixed:
        if x < nv:
            closures.append(int(next(c for c in range(nc) if x in lc_split[c] and all(v in plain for v in lc_split[c]))))
    pmap = numpy.zeros(nv, dtype=numpy.int32)
    active = [v for v in range(nv) if var_of[v] < nv and v not in fixed]
    pmap[active] = numpy.arange(1, len(active) + 1)
    hole_arcs = [a for a in range(len(arcs)) if nv + a not in fixed]
    elements = [int(var_of[v]) for v in active] + [nv + a for a in hole_arcs] + [-1 - j for j in range(len(cuts))]
    I, J, V = [], [], []
    for e, x in enumerate(elements):
        for f, cf in terms_of(x).items():
            if f in col and cf != 0:
                I.append(e)
                J.append(col[f])
                V.append(float(cf))
    P = scipy.sparse.csr_matrix((V, (I, J)), shape=(len(elements), len(cols)))
    cell_of = {}
    for c in range(nc):
        for i in range(3):
            cell_of.setdefault(int(lc_split[c, i]), []).append((c, i))
    explicit = [numpy.array([ci for v in arcs[a][0] for ci in cell_of[int(v)]], dtype=numpy.int32) for a in hole_arcs]
    explicit += [numpy.asarray(cut, dtype=numpy.int32) for cut in cuts]
    info = {'arcs': arcs, 'references': [x - nv for x in fixed if x >= nv],
            'holes': [x - nv for x in cols if x >= nv], 'cuts': [-1 - h for h in hfree],
            'closures': closures, 'shifts': len(gfree), 'pmap': pmap, 'P': P, 'explicit_holes': explicit}
    return M, info


def _min_energy(M, K):
    '''! Smallest eigenvalue of M^T K M relative to its largest entry'''
    A = (M.T @ K @ M).tocsc()
    scale = abs(A).max()
    if A.shape[0] <= 4000:
        return numpy.linalg.eigvalsh(A.toarray()).min()/scale
    val = scipy.sparse.linalg.eigsh(A, k=1, sigma=-1.0E-6*scale, which='LM', return_eigenvectors=False)
    return val.min()/scale


def homology_loops(r, lc):
    '''! Internal (non-boundary) loops of ThinCurr's homology basis for a manifold mesh'''
    from .meshing import write_ThinCurr_mesh, read_ThinCurr_mesh
    from ._homology_cli import compute_homology
    with tempfile.TemporaryDirectory() as tmp:
        in_file = os.path.join(tmp, 'mesh.h5')
        write_ThinCurr_mesh(in_file, r, lc, numpy.ones(len(lc), dtype=numpy.int32))
        with contextlib.redirect_stdout(io.StringIO()):
            compute_homology(in_file)
        holes = read_ThinCurr_mesh(os.path.join(tmp, 'mesh-homology.h5'))['thincurr'].get('holes', [])
    ec = _edge_cells(lc)
    boundary = lambda loop: all(len(ec[(min(a, b), max(a, b))]) == 1 for a, b in zip(loop, numpy.roll(loop, -1)))
    return [numpy.asarray(h) for h in holes if not boundary(h)]


def find_junction_holes(r, lc, reg, insulated=(), tol=1.0E-10):
    '''! Split a junction mesh and choose its holes, cuts and references

    Cut candidates are the closed junction lines and, if the original mesh is manifold, the internal
    loops of ThinCurr's homology basis. Loops that exist only through a non-manifold junction are not
    searched for; the final count check reports them.

    @param r Point list [np,3]
    @param lc Cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @param insulated Regions with chi = 0 at every junction
    @param tol Relative eigenvalue below which a candidate cut is taken to add a zero-energy mode
    @result Dictionary: 'r', 'lc' (split mesh), 'reg', 'lines', 'kje', 'lje' (junction table with
            the sign of each entry in ThinCurr's orientation), 'holes' (explicit holes: [n,2] (cell,
            corner) arrays, value +1 at each, for the non-reference arcs then the cuts), 'cuts', 'pmap'
            (vertex to expanded unknown, 0 inactive), 'vmap' (original vertex of each split vertex), 'P'
            (sparse [n_expanded, n_reduced]; expanded
            unknowns are the active vertices then the explicit holes), 'n_unknowns', 'n_loops', 'M',
            'info'
    '''
    lc = numpy.asarray(lc, dtype=numpy.int32)
    reg = numpy.asarray(reg, dtype=numpy.int32)
    rs, lcs, lines = split_junctions(r, lc, reg)
    candidates = []
    for line in lines:
        if line['closed']:
            k = min(line['segments'][0])
            cells = numpy.where(numpy.any(numpy.isin(lcs, [c[k] for c in line['copies']]), axis=1) & (reg == k))[0]
            on_line = set(c[k] for c in line['copies'])
            candidates.append(numpy.array([(c, i) for c in cells for i in range(3) if lcs[c, i] in on_line],
                                          dtype=numpy.int32))
    if max(len(v) for v in _edge_cells(lc).values()) <= 2:
        candidates += [loop_cut(lc, loop) for loop in homology_loops(r, lc)]
    Kc = _local_stiffness(rs, lcs)
    nc = len(lcs)
    rows = numpy.repeat(3*numpy.arange(nc)[:, None] + numpy.arange(3)[None, :], 3, axis=1).ravel()
    cols = numpy.tile(3*numpy.arange(nc)[:, None] + numpy.arange(3)[None, :], (1, 3)).ravel()
    K = scipy.sparse.csr_matrix((Kc.ravel(), (rows, cols)), shape=(3*nc, 3*nc))
    target = current_loop_count(lc, reg, insulated)
    chosen = []
    M, info = junction_basis(rs, lcs, lines, insulated, chosen)
    for cut in candidates:
        if M.shape[1] >= target:
            break
        M2, info2 = junction_basis(rs, lcs, lines, insulated, chosen + [cut])
        if M2.shape[1] > M.shape[1] and _min_energy(M2, K) > tol:
            chosen.append(cut)
            M, info = M2, info2
    if M.shape[1] != target:
        raise ValueError('Found {0} unknowns but the mesh carries {1} independent currents'.format(M.shape[1], target))
    if _min_energy(M, K) <= tol:
        raise ValueError('Junction model has a zero-energy mode')
    kje, lje = junction_table(lines)
    lco = orient_cells(lcs)
    sign = numpy.array([_edge_sign(lco[c], a, b) for a, b, c in lje], dtype=numpy.int32)
    vmap = numpy.arange(len(rs), dtype=numpy.int32)
    for line in lines:
        for v, copies in zip(line['vertices'], line['copies']):
            vmap[list(copies.values())] = v
    return {'r': rs, 'lc': lcs, 'reg': reg, 'lines': lines, 'kje': kje, 'lje': numpy.c_[lje, sign], 'vmap': vmap,
            'holes': info['explicit_holes'], 'cuts': chosen, 'pmap': info['pmap'], 'P': info['P'],
            'n_unknowns': M.shape[1], 'n_loops': target, 'M': M, 'info': info}


def insulate_junction_file(mesh_file, insulated, out_file):
    '''! Copy a junction mesh file with the junction data for some sheets insulated

    The unsplit mesh is rebuilt from the split one (`vmap`), holes and P are recomputed with the given
    regions insulated (chi = 0 at every junction), and the copy's `thincurr/junctions` group is replaced.

    @param mesh_file Junction mesh file (from the homology tool or @ref find_junction_holes)
    @param insulated Regions to insulate
    @param out_file Output file
    '''
    import shutil
    import h5py
    from .meshing import read_ThinCurr_mesh, write_junction_group
    mesh = read_ThinCurr_mesh(mesh_file)
    vmap = mesh['thincurr']['junctions']['vmap']
    r = numpy.zeros((vmap.max() + 1, 3))
    r[vmap] = mesh['r']
    out = find_junction_holes(r, vmap[mesh['lc']], mesh['reg'], insulated=insulated)
    if (out['lc'].shape != mesh['lc'].shape) or numpy.any(out['lc'] != mesh['lc']):
        raise ValueError('Split mesh in "{0}" does not match the junction split of its unsplit mesh'.format(mesh_file))
    shutil.copyfile(mesh_file, out_file)
    with h5py.File(out_file, 'r+') as h5_file:
        del h5_file['thincurr/junctions']
        write_junction_group(h5_file, out)
    return out
