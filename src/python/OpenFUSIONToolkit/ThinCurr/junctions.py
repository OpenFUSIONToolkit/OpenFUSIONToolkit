#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Mesh tools for ThinCurr models with junctions (edges shared by more than two sheets)

Each region is a sheet. A junction edge is an edge shared by more than two cells; in a mesh that is split, edges
where two regions meet are also junction edges, even where only two cells share them, so that those regions can be
insulated from each other (with two cells, current conservation there is ordinary continuity). The mesh is split
so that, at every vertex, each group of cells connected through ordinary edges (a "fan") has its own copy of the
vertex; the split mesh is manifold and keeps the cells and their regions.
Current conservation across each junction edge (Kirchhoff) joins the cells again:

  sum_c s_c (chi_c(b) - chi_c(a)) = 0 for each junction edge a -> b,

over the cells on that edge, where chi_c(u) is the stream function of cell c at u and s_c = +1 if c lists
the edge in the junction's direction (in ThinCurr's orientation), otherwise -1; s_c (chi_c(b) - chi_c(a))
is the current leaving cell c across the edge.

The currents the model can carry are the cycles of the dual graph whose nodes are the cells and the junction
edges, linked where current can cross (ordinary edges, and joined cells at junction edges). The vertex
unknowns represent most of them; the rest are added as hole elements built from cycles of that graph (a jump
of chi at the vertex where each cell on the cycle is entered and left). Sheet boundaries (chains of outer
edges, with chi uniform along each) and references (sheet boundaries, or one vertex, fixed at zero to remove
constant shifts that carry no current) complete the set. The number of unknowns equals the number of independent cycles.

@authors Ian Stewart
@date October 2026
@ingroup doxy_oft_python
'''
import warnings
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
    '''! Direction in which a cell traverses one of its edges

    @param tri Cell vertices
    @param a, b Edge vertices
    @result +1 if the cell lists the edge as a -> b, -1 if as b -> a
    '''
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


def _junction_edges(lc, reg, ec=None):
    '''! Junction edges, {(a, b): [cells]}: edges with more than two cells, and edges where two regions meet (so that
    those regions can be insulated from each other)'''
    ec = ec or _edge_cells(lc)
    return {e: cells for e, cells in ec.items() if len(cells) > 2 or len(set(int(reg[c]) for c in cells)) > 1}


def _fans(lc, ec, jedges):
    '''! Groups of cells at each vertex connected through ordinary (non-junction, two-cell) edges
    @result {vertex: [[cells], ...]} for vertices with more than one group'''
    vcells = {}
    for c, tri in enumerate(lc):
        for v in tri:
            vcells.setdefault(int(v), []).append(c)
    out = {}
    for v, cells in vcells.items():
        if len(cells) == 1:
            continue
        parent = {c: c for c in cells}

        def find(c):
            '''! Root of the group of cell c (union-find, with path halving)'''
            while parent[c] != c:
                parent[c] = parent[parent[c]]
                c = parent[c]
            return c
        for c in cells:
            for w in lc[c]:
                w = int(w)
                if w == v:
                    continue
                e = (min(v, w), max(v, w))
                if e in jedges or len(ec[e]) != 2:
                    continue
                parent[find(ec[e][0])] = find(ec[e][1])
        groups = {}
        for c in cells:
            groups.setdefault(find(c), []).append(c)
        if len(groups) > 1:
            out[v] = sorted((sorted(g) for g in groups.values()), key=lambda g: g[0])
    return out


def has_junctions(lc, reg=None):
    '''! True if the mesh needs splitting: an edge shared by more than two cells, or (with `reg`) cells of
    different regions that touch only at a vertex'''
    ec = _edge_cells(lc)
    if max(len(v) for v in ec.values()) > 2:
        return True
    if reg is None:
        return False
    reg = numpy.asarray(reg)
    jedges = {e: cells for e, cells in ec.items() if len(cells) > 2}
    for groups in _fans(lc, ec, jedges).values():
        if len(set(int(reg[g[0]]) for g in groups)) > 1:
            return True
    return False


def split_junctions(r, lc, reg):
    '''! Split a mesh so that each fan of cells at a vertex has its own copy of that vertex

    Cells keep their index and region. Junction edges are also grouped into lines between branch points (for
    reporting).

    @param r Point list [np,3]
    @param lc Cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @result r_split, lc_split, jedges ({(a, b) original: [cells]}), lines (list of dicts with 'vertices'
            (original indices) and 'closed'), vmap (original vertex of each split vertex)
    '''
    lc = numpy.asarray(lc)
    reg = numpy.asarray(reg)
    ec = _edge_cells(lc)
    jedges = _junction_edges(lc, reg, ec)
    r_split, vmap, lc_split = [p for p in r], list(range(len(r))), lc.copy()
    for v, groups in _fans(lc, ec, jedges).items():
        for g in groups[1:]:
            r_split.append(numpy.array(r[v]).copy())
            vmap.append(v)
            for c in g:
                lc_split[c][lc[c] == v] = len(r_split) - 1
    adj = {}
    for a, b in jedges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    ends = sorted(v for v, n in adj.items() if len(n) != 2)
    def key(a, b):
        return (min(a, b), max(a, b))
    lines, used = [], set()
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
            lines.append({'vertices': [int(v) for v in path], 'closed': closed})
    return numpy.array(r_split), lc_split, jedges, lines, numpy.array(vmap, dtype=numpy.int32)


def _copy(lc, lc_split, c, v):
    '''! Split vertex of original vertex v in cell c'''
    return int(lc_split[c][list(lc[c]).index(v)])


def junction_table(lc, lc_split, jedges):
    '''! Junction edge table for the mesh file

    @result kje [n_jedges+1] pointer, lje [n,4] (vertex a, vertex b, cell, sign) per cell on each junction
    edge, 0-based split indices; a -> b is the junction's direction and sign is +1 if the cell lists the edge
    in that direction in ThinCurr's orientation
    '''
    lco = orient_cells(lc_split)
    kje, lje = [0], []
    for (a, b), cells in sorted(jedges.items()):
        for c in cells:
            ca, cb = _copy(lc, lc_split, c, a), _copy(lc, lc_split, c, b)
            lje.append([ca, cb, c, _edge_sign(lco[c], ca, cb)])
        kje.append(len(lje))
    return numpy.array(kje, dtype=numpy.int32), numpy.array(lje, dtype=numpy.int32).reshape(-1, 4)


def _insulation(insulated):
    '''! Split an insulation list into insulated regions and insulated pairs of regions'''
    regions, pairs = set(), set()
    for item in insulated:
        if isinstance(item, (tuple, list, set, frozenset, numpy.ndarray)):
            pair = frozenset(int(k) for k in item)
            if len(pair) != 2:
                raise ValueError('Insulated pairs must name two different regions, got {0}'.format(item))
            pairs.add(pair)
        else:
            regions.add(int(item))
    return regions, pairs


def _joined_cells(cells, reg, insulated):
    '''! Cells joined across a junction edge after insulation

    A pair of regions insulated from each other is only separated where they are the only regions on the
    edge. Where a third sheet meets both on the same edge, current can still pass between them through it,
    so the pair has no effect there (a warning is given).'''
    ins_regions, ins_pairs = insulated
    regions = frozenset(int(reg[c]) for c in cells)
    if regions in ins_pairs:
        return []
    for pair in ins_pairs:
        if pair < regions:
            warnings.warn('Regions {0} meet region(s) {1} on the same junction edge; current can pass between them through '
                          'it, so insulating the pair has no effect there'.format(sorted(pair), sorted(regions - pair)))
    return [c for c in cells if int(reg[c]) not in ins_regions]


def dual_graph(lc, reg, insulated=(), jedges=None):
    '''! Dual graph of the currents: nodes are the cells and the junction edges

    @result number of nodes, links [(node u, node w, edge (a, b), cell)]: an ordinary edge links its two
    cells (cell = the first); a junction edge node is linked to each joined cell (node w is the junction).
    '''
    reg = numpy.asarray(reg)
    insulated = _insulation(insulated)
    ec = _edge_cells(lc)
    jedges = jedges if jedges is not None else _junction_edges(lc, reg, ec)
    nc = len(lc)
    links, nj = [], 0
    for e in sorted(ec):
        cells = ec[e]
        if e in jedges:
            joined = _joined_cells(cells, reg, insulated)
            if len(joined) >= 2:
                for c in joined:
                    links.append((c, nc + nj, e, c))
            nj += 1
        elif len(cells) == 2:
            links.append((cells[0], cells[1], e, cells[0]))
    return nc + nj, links


def current_loop_count(lc, reg, insulated=()):
    '''! Number of independent divergence-free current patterns (cycles of @ref dual_graph)'''
    n, links = dual_graph(lc, reg, insulated)
    G = scipy.sparse.coo_matrix((numpy.ones(len(links)), ([l[0] for l in links], [l[1] for l in links])), shape=(n, n))
    return len(links) - n + scipy.sparse.csgraph.connected_components(G, directed=False)[0]


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


def _ordered_chain(edges):
    '''! Order a connected chain of edges; @result (points, closed)'''
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    if max(len(n) for n in adj.values()) > 2:
        return numpy.array(sorted(adj), dtype=numpy.int32), False
    ends = [v for v, n in adj.items() if len(n) == 1]
    start = ends[0] if ends else min(adj)
    chain = [start]
    while True:
        nxt = [w for w in adj[chain[-1]] if len(chain) < 2 or w != chain[-2]]
        if not nxt or nxt[0] == start:
            break
        chain.append(nxt[0])
    return numpy.array(chain, dtype=numpy.int32), not ends


def _local_stiffness(r, lc):
    '''! P1 stiffness of each cell, int grad(lambda_i).grad(lambda_j) dA, [nc,3,3]'''
    p = r[lc]
    e = numpy.stack([p[:, 2] - p[:, 1], p[:, 0] - p[:, 2], p[:, 1] - p[:, 0]], axis=1)
    area = 0.5*numpy.linalg.norm(numpy.cross(e[:, 0], e[:, 1]), axis=1)
    return numpy.einsum('cik,cjk->cij', e, e)/(4.0*area)[:, None, None]


def junction_basis(lc, lc_split, reg, jedges, insulated=(), holes=()):
    '''! Reduced unknowns of a split mesh: map from unknowns to the stream function at each cell corner

    @param lc Original cell list
    @param lc_split, jedges From @ref split_junctions
    @param reg Region list
    @param insulated Insulation, see @ref find_junction_holes
    @param holes Extra hole elements, each [n,3] (cell, split vertex, sign): chi += sign*h at that corner
    @result M (sparse [3 nc, n_reduced], corner (c, i) is row 3c+i), info dict with 'arcs' (sheet boundaries)
            [(points, closed)], 'references' (sheet boundaries fixed at zero), 'closures' (cells), 'shifts', 'pmap',
            'P', 'explicit_holes'
    '''
    lc_split = numpy.asarray(lc_split)
    nv, nc = int(lc_split.max()) + 1, len(lc_split)
    reg = numpy.asarray(reg)
    insulated = _insulation(insulated)
    lco = orient_cells(lc_split)
    corner_hole = {}
    for j, hole in enumerate(holes):
        for c, v, s in hole:
            corner_hole.setdefault((int(c), int(v)), []).append((-1 - j, int(s)))
    joined_edges, rows, crows, jverts = set(), [], [], set()
    for (a, b), cells in jedges.items():
        for c in cells:
            jverts |= {_copy(lc, lc_split, c, a), _copy(lc, lc_split, c, b)}
        joined = _joined_cells(cells, reg, insulated)
        if len(joined) < 2:
            continue
        sign = {}
        for c in joined:
            ca, cb = _copy(lc, lc_split, c, a), _copy(lc, lc_split, c, b)
            joined_edges.add((min(ca, cb), max(ca, cb)))
            sign[c] = _edge_sign(lco[c], ca, cb)
        row, crow = {}, {}
        for c in joined:
            for u, d in ((b, 1), (a, -1)):
                cu = _copy(lc, lc_split, c, u)
                row[cu] = row.get(cu, 0) + d*sign[c]
                for h, s in corner_hole.get((c, cu), []):
                    crow[h] = crow.get(h, 0) + d*sign[c]*s
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
            arcs.append(_ordered_chain(groups[g]))
    rows = [{int(var_of[v]): c for v, c in row.items() if c != 0} for row in rows]
    expr = _solve([{**row, **crow} for row, crow in zip(rows, crows)])
    G = scipy.sparse.coo_matrix((numpy.ones(3*nc), (lc_split.ravel(), numpy.repeat(lc_split[:, 0], 3))),
                                shape=(nv, nv))
    comp = scipy.sparse.csgraph.connected_components(G, directed=False)[1]
    def comp_of(x):
        return comp[x] if x < nv else comp[arcs[x - nv][0][0]]
    ncomp = comp.max() + 1
    grows = []
    for row in rows:
        g = {}
        for x, c in row.items():
            g[int(comp_of(x))] = g.get(int(comp_of(x)), 0) + c
        grows.append({k: c for k, c in g.items() if c != 0})
    gexpr = _solve(grows)
    gfree = [c for c in range(ncomp) if c not in gexpr]
    gauge = numpy.array([[1.0 if c == f else float(gexpr.get(c, {}).get(f, 0)) for c in range(ncomp)]
                         for f in gfree]).reshape(len(gfree), ncomp)
    free = [x for x in sorted(set(var_of.tolist())) if x not in expr]
    plain = set(range(nv)) - jverts
    order = [x for x in free if x >= nv] + [x for x in free if x < nv and x in plain] + [x for x in free if x < nv and x not in plain]
    fixed, block = [], numpy.zeros((len(gfree), 0))
    for x in order:
        if len(fixed) == len(gfree):
            break
        trial = numpy.hstack((block, gauge[:, [comp_of(x)]]))
        if numpy.linalg.matrix_rank(trial) > block.shape[1]:
            block, fixed = trial, fixed + [x]
    if len(fixed) != len(gfree):
        raise ValueError('Could not remove every zero-current shift')
    hfree = [-1 - j for j in range(len(holes)) if -1 - j not in expr]
    cols = [x for x in free if x not in fixed] + hfree
    col = {x: j for j, x in enumerate(cols)}
    def terms_of(x):
        return expr.get(x, {x: 1})
    I, J, V = [], [], []
    for c in range(nc):
        for i in range(3):
            v = int(lc_split[c, i])
            terms = dict(terms_of(int(var_of[v])))
            for h, s in corner_hole.get((c, v), []):
                for f, cf in terms_of(h).items():
                    terms[f] = terms.get(f, 0) + s*cf
            for f, cf in terms.items():
                if f in col and cf != 0:
                    I.append(3*c + i)
                    J.append(col[f])
                    V.append(float(cf))
    M = scipy.sparse.csr_matrix((V, (I, J)), shape=(3*nc, len(cols)))
    closures = []
    for x in fixed:
        if x < nv:
            closures.append(int(next(c for c in range(nc) if x in lc_split[c])))
    pmap = numpy.zeros(nv, dtype=numpy.int32)
    active = [v for v in range(nv) if var_of[v] < nv and v not in fixed]
    pmap[active] = numpy.arange(1, len(active) + 1)
    hole_arcs = [a for a in range(len(arcs)) if nv + a not in fixed]
    elements = [int(var_of[v]) for v in active] + [nv + a for a in hole_arcs] + [-1 - j for j in range(len(holes))]
    I, J, V = [], [], []
    for e, x in enumerate(elements):
        for f, cf in terms_of(x).items():
            if f in col and cf != 0:
                I.append(e)
                J.append(col[f])
                V.append(float(cf))
    P = scipy.sparse.csr_matrix((V, (I, J)), shape=(len(elements), len(cols)))
    cells_of = {}
    for c in range(nc):
        for v in lc_split[c]:
            cells_of.setdefault(int(v), []).append(c)
    explicit = [numpy.array([(c, v, 1) for v in arcs[a][0] for c in cells_of[int(v)]], dtype=numpy.int32)
                for a in hole_arcs]
    explicit += [numpy.asarray(h, dtype=numpy.int32).reshape(-1, 3) for h in holes]
    info = {'arcs': arcs, 'references': [x - nv for x in fixed if x >= nv],
            'holes': [x - nv for x in cols if x >= nv], 'cuts': [-1 - h for h in hfree],
            'closures': closures, 'shifts': len(gfree), 'pmap': pmap, 'P': P, 'explicit_holes': explicit}
    return M, info


def _flux_rows(lc, lc_split, M, links, lco):
    '''! Current across each link (out of its first node) as a sparse row over the reduced unknowns'''
    M = M.tocsr()
    out = []
    for u, w, (a, b), c in links:
        ca, cb = _copy(lc, lc_split, c, a), _copy(lc, lc_split, c, b)
        t = list(lco[c])
        x, y = (ca, cb) if t[(t.index(ca) + 1) % 3] == cb else (cb, ca)
        iy, ix = list(lc_split[c]).index(y), list(lc_split[c]).index(x)
        row = {}
        for ii, s in ((iy, 1.0), (ix, -1.0)):
            r0, r1 = M.indptr[3*c + ii], M.indptr[3*c + ii + 1]
            for j, val in zip(M.indices[r0:r1], M.data[r0:r1]):
                row[int(j)] = row.get(int(j), 0.0) + s*val
        out.append({j: Fraction(val).limit_denominator(1000000) for j, val in row.items() if abs(val) > 1.0E-12})
    return out


def _cycle_holes(lc, lc_split, reg, jedges, insulated, M):
    '''! Hole elements for the currents the reduced unknowns do not represent

    A spanning forest of @ref dual_graph is built; each link outside it closes one cycle, and the cycles
    span all currents. Links are taken in turn and kept as generators when their current, in terms of the
    reduced unknowns, depends on the links already taken (greedy row elimination); each generator's cycle
    becomes a hole element: a jump of chi at the vertex where each cell on the cycle is entered and left.
    @result list of [n,3] (cell, split vertex, sign) arrays
    '''
    n, links = dual_graph(lc, reg, insulated, jedges)
    lco = orient_cells(lc_split)
    nc = len(lc)
    nlink = len(links)
    I = [l[0] for l in links] + [l[1] for l in links]
    J = [l[1] for l in links] + [l[0] for l in links]
    V = list(range(1, nlink + 1)) + list(range(1, nlink + 1))
    A = scipy.sparse.csr_matrix((V, (I, J)), shape=(n, n))
    order, pred = scipy.sparse.csgraph.breadth_first_order(A, 0, directed=False, return_predecessors=True)
    pred = numpy.array(pred)
    seen = numpy.zeros(n, dtype=bool)
    seen[order] = True
    for s in range(n):
        if not seen[s]:
            o, p = scipy.sparse.csgraph.breadth_first_order(A, s, directed=False, return_predecessors=True)
            pred[o] = p[o]
            seen[o] = True
    tree = set()
    for v in range(n):
        if pred[v] >= 0:
            tree.add(int(A[v, pred[v]]) - 1)

    def path_up(v):
        '''! Path from node v up to the root of its spanning tree'''
        out = [v]
        while pred[out[-1]] >= 0:
            out.append(int(pred[out[-1]]))
        return out
    nontree = [k for k in range(nlink) if k not in tree]
    rows = _flux_rows(lc, lc_split, M, [links[k] for k in nontree], lco)
    basis, generators = {}, []
    for k, row in zip(nontree, rows):
        r = dict(row)
        while r:
            p = max(r)
            if p not in basis:
                break
            f = r[p]
            for j, val in basis[p].items():
                r[j] = r.get(j, 0) - f*val
                if r[j] == 0:
                    del r[j]
        if r:
            p = max(r)
            basis[p] = {j: val/r[p] for j, val in r.items()}
        else:
            generators.append(k)
    holes = []
    for k in generators:
        u, w = links[k][0], links[k][1]
        pu, pw = path_up(u), path_up(w)
        common = set(pu) & set(pw)
        lca = next(x for x in pu if x in common)
        seq = [u] + pw[:pw.index(lca) + 1] + pu[1:pu.index(lca)][::-1]
        m = len(seq)
        edge_of = {}
        for (a_, b_) in zip(seq, seq[1:] + seq[:1]):
            edge_of[(a_, b_)] = links[int(A[a_, b_]) - 1][2]
        entries = []
        for i, x in enumerate(seq):
            if x >= nc:
                continue
            e_in, e_out = edge_of[(seq[i - 1], x)], edge_of[(x, seq[(i + 1) % m])]
            vin = {_copy(lc, lc_split, x, e_in[0]), _copy(lc, lc_split, x, e_in[1])}
            vout = (_copy(lc, lc_split, x, e_out[0]), _copy(lc, lc_split, x, e_out[1]))
            wv = (vin & set(vout)).pop()
            t = list(lco[x])
            a_, b_ = vout
            xx, yy = (a_, b_) if t[(t.index(a_) + 1) % 3] == b_ else (b_, a_)
            entries.append((x, wv, 1 if yy == wv else -1))
        holes.append(numpy.array(entries, dtype=numpy.int32))
    return holes


def find_junction_holes(r, lc, reg, insulated=(), tol=1.0E-10):
    '''! Split a junction mesh and choose its holes and references

    @param r Point list [np,3]
    @param lc Cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @param insulated Insulation: regions (insulated at every junction they touch) and/or pairs of regions
    (insulated from each other where they are the only regions on a junction edge)
    @param tol Relative eigenvalue below which the model is taken to have a zero-energy mode
    @result Dictionary: 'r', 'lc' (split mesh), 'reg', 'lines', 'jedges', 'kje', 'lje' (junction table),
            'vmap', 'holes' (explicit holes, [n,3] (cell, split vertex, sign) arrays: non-reference sheet
            boundaries, then cycle holes), 'cuts' (the cycle holes), 'pmap', 'P' (sparse [n_expanded, n_reduced];
            expanded unknowns are the active vertices then the explicit holes), 'n_unknowns', 'n_loops', 'M', 'info',
            'insulated'
    '''
    lc = numpy.asarray(lc, dtype=numpy.int32)
    reg = numpy.asarray(reg, dtype=numpy.int32)
    rs, lcs, jedges, lines, vmap = split_junctions(r, lc, reg)
    target = current_loop_count(lc, reg, insulated)
    M, info = junction_basis(lc, lcs, reg, jedges, insulated)
    cuts = []
    if M.shape[1] < target:
        cuts = _cycle_holes(lc, lcs, reg, jedges, insulated, M)
        M, info = junction_basis(lc, lcs, reg, jedges, insulated, cuts)
    if M.shape[1] != target:
        raise ValueError('Found {0} unknowns but the mesh carries {1} independent currents'.format(M.shape[1], target))
    Kc = _local_stiffness(rs, lcs)
    nc = len(lcs)
    rows = numpy.repeat(3*numpy.arange(nc)[:, None] + numpy.arange(3)[None, :], 3, axis=1).ravel()
    cols = numpy.tile(3*numpy.arange(nc)[:, None] + numpy.arange(3)[None, :], (1, 3)).ravel()
    K = scipy.sparse.csr_matrix((Kc.ravel(), (rows, cols)), shape=(3*nc, 3*nc))
    if _min_energy(M, K) <= tol:
        raise ValueError('Junction model has a zero-energy mode')
    kje, lje = junction_table(lc, lcs, jedges)
    return {'r': rs, 'lc': lcs, 'reg': reg, 'lines': lines, 'jedges': jedges, 'kje': kje, 'lje': lje, 'vmap': vmap,
            'holes': info['explicit_holes'], 'cuts': cuts, 'pmap': info['pmap'], 'P': info['P'],
            'n_unknowns': M.shape[1], 'n_loops': target, 'M': M, 'info': info, 'insulated': list(insulated)}


def _min_energy(M, K):
    '''! Smallest pivot of a sparse LU factorization of M^T K M relative to the largest (zero for a zero-energy mode)'''
    A = (M.T @ K @ M).tocsc()
    try:
        lu = scipy.sparse.linalg.splu(A, permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.0,
                                      options={'SymmetricMode': True})
    except RuntimeError:
        return 0.0
    d = numpy.abs(lu.U.diagonal())
    return d.min()/d.max()


def split_jumper(lc, lc_split, vmap, jedges, jumper):
    '''! Map a jumper (vertex set of the unsplit mesh) onto the split mesh

    The vertices are ordered into a path along mesh edges and each path edge takes the vertex copies of its
    cells. Where the path crosses a junction line, consecutive edges use different copies of the crossing
    vertex and a new segment starts; the jumper current is the sum over the segments. A segment on a sheet
    oriented opposite to the previous one (both traverse a shared junction edge the same way) is reversed, so
    that every segment measures current crossing the path in the same direction.

    @param lc Unsplit cell list [nc,3] (0-based)
    @param lc_split Split cell list [nc,3]
    @param vmap Unsplit vertex of each split vertex
    @param jedges Junction edges {(a, b): [cells]} of the unsplit mesh
    @param jumper Jumper vertices (unsplit, 0-based)
    @result Ordered split vertices, segments separated by -1 (a closed single segment repeats its first vertex)
    '''
    verts = set(int(v) for v in jumper)
    ec = _edge_cells(lc)
    adj = {v: [] for v in verts}
    for a, b in ec:
        if a in verts and b in verts:
            adj[a].append(b)
            adj[b].append(a)
    ends = sorted(v for v in verts if len(adj[v]) == 1)
    if any(len(n) > 2 for n in adj.values()) or len(ends) not in (0, 2):
        raise ValueError('Jumper is not a simple path or loop of mesh edges')
    path = [ends[0] if ends else min(verts)]
    while True:
        nxt = [w for w in sorted(adj[path[-1]]) if w not in path[-2:-1] and (w not in path or (not ends and w == path[0]))]
        if not nxt or (len(path) > 1 and path[-1] == path[0]):
            break
        path.append(nxt[0])
    if len(set(path)) != len(verts):
        raise ValueError('Jumper is not a single connected path or loop of mesh edges')
    edges = []
    for a, b in zip(path[:-1], path[1:]):
        if (min(a, b), max(a, b)) in jedges:
            raise ValueError('Jumper runs along a junction edge ({0}, {1}), where the sheet it measures is ambiguous'.format(a, b))
        c = ec[(min(a, b), max(a, b))][0]
        copy = {int(vmap[v]): int(v) for v in lc_split[c]}
        edges.append((copy[a], copy[b]))
    breaks = [k for k in range(1, len(edges)) if edges[k][0] != edges[k-1][1]]
    if not ends and breaks and edges[-1][1] == edges[0][0]:
        edges = edges[breaks[0]:] + edges[:breaks[0]]
        breaks = [k for k in range(1, len(edges)) if edges[k][0] != edges[k-1][1]]
    segments = [edges[i:j] for i, j in zip([0] + breaks, breaks + [len(edges)])]
    chains = [[seg[0][0]] + [b for _, b in seg] for seg in segments]
    flip = False
    for k in range(1, len(segments)):
        if _copies_opposed(lc_split, vmap, jedges, segments[k-1][-1][1], segments[k][0][0]):
            flip = not flip
        if flip:
            chains[k] = chains[k][::-1]
    out = chains[0]
    for chain in chains[1:]:
        out = out + [-1] + chain
    return numpy.array(out, dtype=numpy.int32)


def _copies_opposed(lc_split, vmap, jedges, va, vb):
    '''! True if the sheets of two copies of a vertex are oppositely oriented (both traverse a junction edge they
    share in the same direction)'''
    v = int(vmap[va])
    for (a, b), cells in jedges.items():
        if v not in (a, b):
            continue
        w = b if a == v else a
        signs = {}
        for c in cells:
            tri = [int(x) for x in lc_split[c]]
            for copy in (va, vb):
                if copy in tri:
                    signs[copy] = 1 if int(vmap[tri[(tri.index(copy) + 1) % 3]]) == w else -1
        if len(signs) == 2:
            return signs[va] == signs[vb]
    raise ValueError('Jumper crosses between sheets that share no junction edge at vertex {0}'.format(v))


def verify_junction_mesh(r, lc, reg, junctions):
    '''! Check a split mesh and its junction data against the unsplit mesh rebuilt from `vmap`

    The split and the number of junction edges must match, and the number of reduced unknowns must equal the
    number of independent currents, with the insulation stored in the file.

    @param r Split point list [np,3]
    @param lc Split cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @param junctions Junction data from @ref meshing.read_ThinCurr_mesh
    @result Number of junction edges, number of reduced unknowns
    '''
    vmap = junctions['vmap']
    r0 = numpy.zeros((vmap.max() + 1, 3))
    r0[vmap] = r
    lc0 = vmap[lc]
    _, lcs, jedges, _, _ = split_junctions(r0, lc0, reg)
    if (lcs.shape != lc.shape) or numpy.any(lcs != lc):
        raise ValueError('Split mesh does not match the junction split of its unsplit mesh, rerun "OFT_ThinCurr_holes"')
    njedges = len(junctions['kje']) - 1
    if njedges != len(jedges):
        raise ValueError('Incorrect number of junction edges found {0} in file (expected {1})'.format(njedges, len(jedges)))
    nunknowns = junctions['P'].shape[1]
    nloops = current_loop_count(lc0, reg, junctions.get('insulated', []))
    if nunknowns != nloops:
        raise ValueError('Incorrect number of junction model unknowns found {0} in file (expected {1})'.format(nunknowns, nloops))
    return njedges, nunknowns


def insulate_junction_file(mesh_file, insulated, out_file):
    '''! Copy a junction mesh file with the junction data for some sheets insulated

    The unsplit mesh is rebuilt from the split one (`vmap`), holes and P are recomputed with the given
    insulation (see @ref find_junction_holes), and the copy's `thincurr/junctions` group is replaced.

    @param mesh_file Junction mesh file (from the homology tool or @ref find_junction_holes)
    @param insulated Insulation: regions and/or pairs of regions
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
