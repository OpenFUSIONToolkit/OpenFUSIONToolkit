#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Mesh tools for ThinCurr models with junctions (edges shared by more than two sheets)

Each region is a sheet. A junction edge is an edge shared by more than two cells; if requested (`split_regions`),
edges where two regions meet are also junction edges, even where only two cells share them, so that those regions
can be insulated from each other (with two cells, current conservation there is ordinary continuity). The mesh is
split so that, at every vertex, each group of cells connected through ordinary edges (a "fan") has its own copy of
the vertex; the split mesh is manifold and keeps the cells and their regions.
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


def _edge_table(lc):
    '''! Edges of a cell list and the cells on each

    Corner 3c+i of cell c is vertex lc[c,i] and the edge lc[c,i] -> lc[c,(i+1)%3].

    @param lc Cell list [nc,3] (0-based)
    @result Dictionary: 'edges' [ne,2] (a < b, sorted), 'ptr' [ne+1] and 'corners' [3 nc] (corners on each edge,
            in increasing order), 'cells' (cell of each entry of 'corners'), 'count' [ne] (number of cells on each
            edge), 'first' [ne] (first corner on each edge)
    '''
    lc = numpy.asarray(lc, dtype=numpy.int64).reshape(-1, 3)
    a = lc.ravel()
    b = lc[:, [1, 2, 0]].ravel()
    lo = numpy.minimum(a, b)
    hi = numpy.maximum(a, b)
    nv = int(a.max()) + 1 if a.size > 0 else 1
    key = lo*nv + hi
    corners = numpy.argsort(key, kind='stable')
    skey = key[corners]
    new = numpy.ones(len(skey), dtype=bool)
    new[1:] = skey[1:] != skey[:-1]
    starts = numpy.flatnonzero(new)
    ptr = numpy.append(starts, len(skey))
    return {'edges': numpy.stack((lo[corners[starts]], hi[corners[starts]]), axis=1), 'ptr': ptr,
            'corners': corners, 'cells': corners//3, 'count': numpy.diff(ptr), 'first': corners[starts]}


def _edge_dict(tab, edges):
    '''! Cells on some edges of an edge table, {(a, b): [cells]}, in the given order of edges'''
    out = {}
    ptr, cells = tab['ptr'], tab['cells']
    for e, (a, b) in zip(edges.tolist(), tab['edges'][edges].tolist()):
        out[(a, b)] = cells[ptr[e]:ptr[e + 1]].tolist()
    return out


def _edge_cells(lc):
    '''! Cells on each edge, {(a, b) with a < b: [cells]}, in order of first appearance'''
    tab = _edge_table(lc)
    return _edge_dict(tab, numpy.argsort(tab['first']))


def _edge_sign(tri, a, b):
    '''! Direction in which a cell traverses one of its edges

    @param tri Cell vertices
    @param a, b Edge vertices
    @result +1 if the cell lists the edge as a -> b, -1 if as b -> a
    '''
    t = list(tri)
    return 1 if t[(t.index(a) + 1) % 3] == b else -1


def _spanning_forest(A):
    '''! Breadth-first spanning forest of an undirected graph, rooted at the lowest node of each connected piece

    Equivalent to a breadth-first search from node 0 and then from each node not yet reached, in turn, but in
    a single search (from an extra node linked to each root).

    @param A Sparse adjacency matrix [n,n]
    @result order (nodes in the order reached), pred (predecessor of each node, -9999 for roots)
    '''
    n = A.shape[0]
    ncomp, lab = scipy.sparse.csgraph.connected_components(A, directed=False)
    roots = numpy.full(ncomp, n, dtype=numpy.int64)
    numpy.minimum.at(roots, lab, numpy.arange(n))
    A = A.tocoo()
    I = numpy.concatenate((A.row, A.col, roots, numpy.full(ncomp, n)))
    J = numpy.concatenate((A.col, A.row, numpy.full(ncomp, n), roots))
    B = scipy.sparse.csr_matrix((numpy.ones(len(I)), (I, J)), shape=(n + 1, n + 1))
    order, pred = scipy.sparse.csgraph.breadth_first_order(B, n, directed=False, return_predecessors=True)
    pred = numpy.array(pred[:n])
    pred[pred == n] = -9999
    return order[1:], pred


def orient_cells(lc):
    '''! Cell orientation as ThinCurr sets it on load (`sync_face_normals`)

    In each connected piece the lowest-numbered cell keeps its order and its neighbours are flipped to
    agree with it. Requires a manifold mesh.

    @param lc Cell list [nc,3] (0-based)
    @result Oriented cell list [nc,3]
    '''
    lc = numpy.array(lc).copy()
    nc = len(lc)
    if nc == 0:
        return lc
    tab = _edge_table(lc)
    if tab['count'].max() > 2:
        raise ValueError('orient_cells requires a manifold mesh')
    two = numpy.flatnonzero(tab['count'] == 2)
    k1 = tab['corners'][tab['ptr'][two]]
    k2 = tab['corners'][tab['ptr'][two] + 1]
    # Neighbours that traverse their shared edge in the same direction disagree
    flat = lc.ravel()
    flip = (flat[k1] == flat[k2]).astype(numpy.int64)
    c1, c2 = k1//3, k2//3
    keep = c1 != c2
    c1, c2, flip = c1[keep], c2[keep], flip[keep]
    _, uniq = numpy.unique(numpy.minimum(c1, c2)*nc + numpy.maximum(c1, c2), return_index=True)
    c1, c2, flip = c1[uniq], c2[uniq], flip[uniq]
    A = scipy.sparse.csr_matrix((numpy.concatenate((flip, flip)) + 1, (numpy.concatenate((c1, c2)),
                                                                      numpy.concatenate((c2, c1)))), shape=(nc, nc))
    order, pred = _spanning_forest(A)
    rel = numpy.zeros(nc, dtype=numpy.int64)
    nonroot = order[pred[order] >= 0]
    rel[nonroot] = numpy.asarray(A[nonroot, pred[nonroot]]).ravel() - 1
    flipped = [0]*nc
    for c, p, f in zip(order.tolist(), pred[order].tolist(), rel[order].tolist()):
        if p >= 0:
            flipped[c] = flipped[p] ^ f
    flipped = numpy.array(flipped, dtype=bool)
    lc[flipped] = lc[flipped][:, ::-1]
    return lc


def _junction_mask(tab, reg, split_regions):
    '''! Junction edges of an edge table: edges with more than two cells and, if `split_regions`, edges where two
    regions meet (so that those regions can be insulated from each other)'''
    mask = tab['count'] > 2
    if split_regions and len(mask) > 0:
        regc = numpy.asarray(reg)[tab['cells']]
        starts = tab['ptr'][:-1]
        mask |= numpy.minimum.reduceat(regc, starts) != numpy.maximum.reduceat(regc, starts)
    return mask


def _junction_edges(tab, jmask):
    '''! Junction edges, {(a, b): [cells]}, in order of first appearance'''
    jedges = numpy.flatnonzero(jmask)
    return _edge_dict(tab, jedges[numpy.argsort(tab['first'][jedges], kind='stable')])


def _fan_groups(lc, tab, jmask):
    '''! Groups of cells at each vertex connected through ordinary (non-junction, two-cell) edges ("fans")

    @param lc Cell list [nc,3] (0-based)
    @param tab Edge table from @ref _edge_table
    @param jmask Junction edges of `tab`
    @result Fan of each corner [3 nc], vertex of each fan, lowest cell of each fan, number of fans at each vertex
    '''
    lc = numpy.asarray(lc, dtype=numpy.int64)
    flat = lc.ravel()
    # Link the corners of the two cells on each ordinary edge at both of its vertices
    ordinary = numpy.flatnonzero((tab['count'] == 2) & (~jmask))
    k1 = tab['corners'][tab['ptr'][ordinary]]
    k2 = tab['corners'][tab['ptr'][ordinary] + 1]
    n1 = 3*(k1//3) + (k1 % 3 + 1) % 3
    n2 = 3*(k2//3) + (k2 % 3 + 1) % 3
    same = flat[k1] == flat[k2]
    I = numpy.concatenate((k1, n1))
    J = numpy.concatenate((numpy.where(same, k2, n2), numpy.where(same, n2, k2)))
    G = scipy.sparse.coo_matrix((numpy.ones(len(I)), (I, J)), shape=(len(flat), len(flat)))
    nlab, lab = scipy.sparse.csgraph.connected_components(G, directed=False)
    lab_vert = numpy.zeros(nlab, dtype=numpy.int64)
    lab_vert[lab] = flat
    lab_cell = numpy.full(nlab, len(lc), dtype=numpy.int64)
    numpy.minimum.at(lab_cell, lab, numpy.arange(len(flat))//3)
    nfans = numpy.bincount(lab_vert, minlength=int(flat.max()) + 1 if flat.size > 0 else 0)
    return lab, lab_vert, lab_cell, nfans


def has_junctions(lc, reg=None):
    '''! True if the mesh needs splitting: an edge shared by more than two cells, or (with `reg`) cells of
    different regions that touch only at a vertex'''
    tab = _edge_table(lc)
    if tab['count'].max() > 2:
        return True
    if reg is None:
        return False
    reg = numpy.asarray(reg)
    _, lab_vert, lab_cell, nfans = _fan_groups(lc, tab, numpy.zeros(len(tab['count']), dtype=bool))
    multi = nfans[lab_vert] > 1
    verts, regs = lab_vert[multi], reg[lab_cell[multi]].astype(numpy.int64)
    rmin = numpy.full(len(nfans), numpy.iinfo(numpy.int64).max, dtype=numpy.int64)
    rmax = numpy.full(len(nfans), numpy.iinfo(numpy.int64).min, dtype=numpy.int64)
    numpy.minimum.at(rmin, verts, regs)
    numpy.maximum.at(rmax, verts, regs)
    return bool(numpy.any(rmin[verts] != rmax[verts]))


def _split(r, lc, tab, jmask):
    '''! Split a mesh so that each fan of cells at a vertex has its own copy of that vertex

    The fan with the lowest cell keeps the vertex. The other copies are added in order of the first appearance of
    their vertex in the cell list, then of their lowest cell.

    @result r_split, lc_split, vmap
    '''
    lc = numpy.asarray(lc)
    r = numpy.asarray(r)
    nv = len(r)
    flat = lc.ravel().astype(numpy.int64)
    lab, lab_vert, lab_cell, _ = _fan_groups(lc, tab, jmask)
    vmin = numpy.full(nv, len(lc), dtype=numpy.int64)
    numpy.minimum.at(vmin, lab_vert, lab_cell)
    extra = numpy.flatnonzero(lab_cell != vmin[lab_vert])
    first = numpy.full(nv, len(flat), dtype=numpy.int64)
    numpy.minimum.at(first, flat, numpy.arange(len(flat)))
    extra = extra[numpy.lexsort((lab_cell[extra], first[lab_vert[extra]]))]
    new_vert = lab_vert.copy()
    new_vert[extra] = nv + numpy.arange(len(extra))
    lc_split = new_vert[lab].reshape(lc.shape).astype(lc.dtype)
    vmap = numpy.concatenate((numpy.arange(nv), lab_vert[extra])).astype(numpy.int32)
    return numpy.vstack((r, r[lab_vert[extra]])), lc_split, vmap


def _junction_lines(jedges):
    '''! Junction edges grouped into lines between branch points (for reporting)'''
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
    return lines


def split_junctions(r, lc, reg, split_regions=False):
    '''! Split a mesh so that each fan of cells at a vertex has its own copy of that vertex

    Cells keep their index and region. Junction edges are also grouped into lines between branch points (for
    reporting).

    @param r Point list [np,3]
    @param lc Cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @param split_regions Also split at all edges where two regions meet?
    @result r_split, lc_split, jedges ({(a, b) original: [cells]}), lines (list of dicts with 'vertices'
            (original indices) and 'closed'), vmap (original vertex of each split vertex)
    '''
    lc = numpy.asarray(lc)
    tab = _edge_table(lc)
    jmask = _junction_mask(tab, reg, split_regions)
    r_split, lc_split, vmap = _split(r, lc, tab, jmask)
    jedges = _junction_edges(tab, jmask)
    return r_split, lc_split, jedges, _junction_lines(jedges), vmap


def _copy(lc, lc_split, c, v):
    '''! Split vertex of original vertex v in cell c'''
    return int(lc_split[c][list(lc[c]).index(v)])


def junction_table(lc, lc_split, jedges, lco=None):
    '''! Junction edge table for the mesh file

    @param lc Original cell list
    @param lc_split, jedges From @ref split_junctions
    @param lco Split cell list oriented by @ref orient_cells (computed if not given)
    @result kje [n_jedges+1] pointer, lje [n,4] (vertex a, vertex b, cell, sign) per cell on each junction
    edge, 0-based split indices; a -> b is the junction's direction and sign is +1 if the cell lists the edge
    in that direction in ThinCurr's orientation
    '''
    lco = orient_cells(lc_split) if lco is None else lco
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


def _check_insulation(tab, jmask, reg, insulated):
    '''! Warn about insulated regions and pairs of regions that are not separated at any junction edge'''
    ins_regions, ins_pairs = _insulation(insulated)
    if not (ins_regions or ins_pairs):
        return
    ptr, cells = tab['ptr'], tab['cells']
    reg = numpy.asarray(reg)
    found_regions, found_pairs = set(), set()
    for e in numpy.flatnonzero(jmask).tolist():
        regions = frozenset(reg[cells[ptr[e]:ptr[e + 1]]].tolist())
        found_regions |= regions
        found_pairs.add(regions)
    for k in sorted(ins_regions - found_regions):
        warnings.warn('Insulated region {0} is not on any junction edge, so insulating it has no effect (use '
                      '"split_regions" to also split the mesh where regions meet)'.format(k))
    for pair in ins_pairs - found_pairs:
        if not any(pair < regions for regions in found_pairs):
            warnings.warn('Insulated regions {0} do not meet on any junction edge, so insulating them has no effect (use '
                          '"split_regions" to also split the mesh where regions meet)'.format(sorted(pair)))


def _dual_links(lc, reg, insulated, tab, jmask):
    '''! Links of @ref dual_graph as arrays

    @result number of nodes, node u, node w, edge vertex a, edge vertex b, cell of each link
    '''
    reg = numpy.asarray(reg)
    insulated = _insulation(insulated)
    nc = len(lc)
    ptr, cells, count = tab['ptr'], tab['cells'], tab['count']
    ordinary = numpy.flatnonzero((count == 2) & (~jmask))
    junc = numpy.flatnonzero(jmask)
    je, jpos, jcell, jnode = [], [], [], []
    for nj, e in enumerate(junc.tolist()):
        joined = _joined_cells(cells[ptr[e]:ptr[e + 1]].tolist(), reg, insulated)
        if len(joined) >= 2:
            for pos, c in enumerate(joined):
                je.append(e)
                jpos.append(pos)
                jcell.append(c)
                jnode.append(nc + nj)
    e = numpy.concatenate((ordinary, numpy.array(je, dtype=numpy.int64)))
    pos = numpy.concatenate((numpy.zeros(len(ordinary), dtype=numpy.int64), numpy.array(jpos, dtype=numpy.int64)))
    u = numpy.concatenate((cells[ptr[ordinary]], numpy.array(jcell, dtype=numpy.int64)))
    w = numpy.concatenate((cells[ptr[ordinary] + 1], numpy.array(jnode, dtype=numpy.int64)))
    order = numpy.lexsort((pos, e))
    e, u, w = e[order], u[order], w[order]
    return nc + len(junc), u, w, tab['edges'][e, 0], tab['edges'][e, 1], u.copy()


def _loop_count(n, u, w):
    '''! Number of independent cycles of a graph with n nodes and links u-w'''
    G = scipy.sparse.coo_matrix((numpy.ones(len(u)), (u, w)), shape=(n, n))
    return len(u) - n + scipy.sparse.csgraph.connected_components(G, directed=False)[0]


def dual_graph(lc, reg, insulated=(), jedges=None, split_regions=False):
    '''! Dual graph of the currents: nodes are the cells and the junction edges

    @param lc Cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @param insulated Insulation, see @ref find_junction_holes
    @param jedges Junction edges (default: from `split_regions`)
    @param split_regions Also treat all edges where two regions meet as junction edges (only used without `jedges`)?
    @result number of nodes, links [(node u, node w, edge (a, b), cell)]: an ordinary edge links its two
    cells (cell = the first); a junction edge node is linked to each joined cell (node w is the junction).
    '''
    tab = _edge_table(lc)
    if jedges is None:
        jmask = _junction_mask(tab, reg, split_regions)
    else:
        jmask = numpy.array([tuple(e) in jedges for e in tab['edges'].tolist()], dtype=bool)
    n, u, w, a, b, cell = _dual_links(lc, reg, insulated, tab, jmask)
    links = [(uu, ww, (aa, bb), cc) for uu, ww, aa, bb, cc in zip(u.tolist(), w.tolist(), a.tolist(), b.tolist(),
                                                                   cell.tolist())]
    return n, links


def current_loop_count(lc, reg, insulated=(), split_regions=False):
    '''! Number of independent divergence-free current patterns (cycles of @ref dual_graph)'''
    tab = _edge_table(lc)
    n, u, w, _, _, _ = _dual_links(lc, reg, insulated, tab, _junction_mask(tab, reg, split_regions))
    return _loop_count(n, u, w)


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


def _add_independent(r, basis):
    '''! Add a sparse row to an echelon basis if it is independent of the rows already in it

    The row is reduced by the basis rows (exact arithmetic, pivoting on the largest index); anything left is
    normalized and added to the basis.

    @param r Row {index: value} (integer or Fraction values; modified)
    @param basis Echelon basis {pivot index: row}, updated in place
    @result True if the row was independent (and added)
    '''
    while r:
        p = max(r)
        if p not in basis:
            break
        f = r[p]
        for j, val in basis[p].items():
            rj = r.get(j, 0) - f*val
            if rj == 0:
                r.pop(j, None)
            else:
                r[j] = rj
    if not r:
        return False
    p = max(r)
    rp = r[p]
    if rp == 1:
        basis[p] = r
    else:
        basis[p] = {j: Fraction(val)/rp for j, val in r.items()}
    return True


def _vertex_cells(lc_split, nv):
    '''! Cells at each vertex, in increasing order; @result pointer [nv+1], cells'''
    flat = numpy.asarray(lc_split, dtype=numpy.int64).ravel()
    order = numpy.argsort(flat, kind='stable')
    ptr = numpy.zeros(nv + 1, dtype=numpy.int64)
    ptr[1:] = numpy.cumsum(numpy.bincount(flat, minlength=nv))
    return ptr, order//3


def junction_basis(lc, lc_split, reg, jedges, insulated=(), holes=(), lco=None):
    '''! Reduced unknowns of a split mesh: map from unknowns to the stream function at each cell corner

    @param lc Original cell list
    @param lc_split, jedges From @ref split_junctions
    @param reg Region list
    @param insulated Insulation, see @ref find_junction_holes
    @param holes Extra hole elements, each [n,3] (cell, split vertex, sign): chi += sign*h at that corner
    @param lco Split cell list oriented by @ref orient_cells (computed if not given)
    @result M (sparse [3 nc, n_reduced], corner (c, i) is row 3c+i), info dict with 'arcs' (sheet boundaries)
            [(points, closed)], 'references' (sheet boundaries fixed at zero), 'closures' (cells), 'shifts', 'pmap',
            'P', 'explicit_holes'
    '''
    lc = numpy.asarray(lc)
    lc_split = numpy.asarray(lc_split)
    nv, nc = int(lc_split.max()) + 1, len(lc_split)
    reg = numpy.asarray(reg)
    insulated = _insulation(insulated)
    lco = orient_cells(lc_split) if lco is None else lco
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
    # Sheet boundaries: chains of outer edges (one cell, not joined), in order of first appearance
    stab = _edge_table(lc_split)
    single = numpy.flatnonzero(stab['count'] == 1)
    single = single[numpy.argsort(stab['first'][single], kind='stable')]
    outer = stab['edges'][single]
    if joined_edges:
        jkeys = numpy.array([a*nv + b for a, b in joined_edges], dtype=numpy.int64)
        outer = outer[~numpy.isin(outer[:, 0]*nv + outer[:, 1], jkeys)]
    var_of = numpy.arange(nv)
    arcs = []
    if len(outer) > 0:
        G = scipy.sparse.coo_matrix((numpy.ones(len(outer)), (outer[:, 0], outer[:, 1])), shape=(nv, nv))
        lab = scipy.sparse.csgraph.connected_components(G, directed=False)[1][outer[:, 0]]
        order = numpy.argsort(lab, kind='stable')
        bounds = numpy.flatnonzero(numpy.diff(lab[order])) + 1
        for group in numpy.split(order, bounds):
            edges = outer[group]
            var_of[numpy.unique(edges)] = nv + len(arcs)
            arcs.append(_ordered_chain([(a, b) for a, b in edges.tolist()]))
    var_rows = []
    for row in rows:
        vrow = {}
        for v, c in row.items():
            x = int(var_of[v])
            vrow[x] = vrow.get(x, 0) + c
        var_rows.append({x: c for x, c in vrow.items() if c != 0})
    rows = var_rows
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
    gindex = {f: i for i, f in enumerate(gfree)}

    def gauge(c):
        '''! Zero-current shift removed by fixing an unknown of piece c, over the free shifts'''
        if c not in gexpr:
            return {gindex[c]: 1}
        return {gindex[f]: cf for f, cf in gexpr[c].items() if cf != 0}
    free = [x for x in sorted(set(var_of.tolist())) if x not in expr]
    plain = set(range(nv)) - jverts
    order = [x for x in free if x >= nv] + [x for x in free if x < nv and x in plain] + [x for x in free if x < nv and x not in plain]
    # Fix one unknown per zero-current shift; whether an unknown removes a shift depends only on its piece, so each
    # piece is tried once
    fixed, basis, tried = [], {}, set()
    for x in order:
        if len(fixed) == len(gfree):
            break
        cx = int(comp_of(x))
        if cx in tried:
            continue
        tried.add(cx)
        if _add_independent(gauge(cx), basis):
            fixed.append(x)
    if len(fixed) != len(gfree):
        raise ValueError('Could not remove every zero-current shift')
    hfree = [-1 - j for j in range(len(holes)) if -1 - j not in expr]
    cols = [x for x in free if x not in fixed] + hfree
    col = {x: j for j, x in enumerate(cols)}
    def terms_of(x):
        return expr.get(x, {x: 1})
    # Map from all variables (vertices, sheet boundaries, holes) to the reduced unknowns
    narcs = len(arcs)
    def var_index(x):
        return x if x >= 0 else nv + narcs - 1 - x
    I, J, V = [], [], []
    for x in cols:
        if x not in expr:
            I.append(var_index(x))
            J.append(col[x])
            V.append(1.0)
    for x, e in expr.items():
        for f, cf in e.items():
            if f in col and cf != 0:
                I.append(var_index(x))
                J.append(col[f])
                V.append(float(cf))
    T = scipy.sparse.csr_matrix((V, (I, J)), shape=(nv + narcs + len(holes), len(cols)))
    # Corners: the variable of their vertex, plus any hole jumps (summed exactly)
    hole_corners = set()
    I, J, V = [], [], []
    for (c, v), hs in corner_hole.items():
        k = 3*c + list(lc_split[c]).index(v)
        hole_corners.add(k)
        terms = dict(terms_of(int(var_of[v])))
        for h, s in hs:
            for f, cf in terms_of(h).items():
                terms[f] = terms.get(f, 0) + s*cf
        for f, cf in terms.items():
            if f in col and cf != 0:
                I.append(k)
                J.append(col[f])
                V.append(float(cf))
    plain_corners = numpy.ones(3*nc, dtype=bool)
    plain_corners[list(hole_corners)] = False
    plain_corners = numpy.flatnonzero(plain_corners)
    S = scipy.sparse.csr_matrix((numpy.ones(len(plain_corners)), (plain_corners, var_of[lc_split.ravel()[plain_corners]])),
                                shape=(3*nc, T.shape[0]))
    M = (S @ T + scipy.sparse.csr_matrix((V, (I, J)), shape=(3*nc, len(cols)))).tocsr()
    first_cell = numpy.full(nv, nc, dtype=numpy.int64)
    numpy.minimum.at(first_cell, lc_split.ravel(), numpy.repeat(numpy.arange(nc), 3))
    closures = [int(first_cell[x]) for x in fixed if x < nv]
    pmap = numpy.zeros(nv, dtype=numpy.int32)
    is_active = var_of < nv
    is_active[[x for x in fixed if x < nv]] = False
    active = numpy.flatnonzero(is_active)
    pmap[active] = numpy.arange(1, len(active) + 1)
    hole_arcs = [a for a in range(narcs) if nv + a not in fixed]
    elements = numpy.concatenate((active, nv + numpy.array(hole_arcs, dtype=numpy.int64),
                                  nv + narcs + numpy.arange(len(holes)))).astype(numpy.int64)
    P = T[elements]
    vptr, vcells = _vertex_cells(lc_split, nv)
    explicit = []
    for a in hole_arcs:
        entries = [(c, v, 1) for v in arcs[a][0].tolist() for c in vcells[vptr[v]:vptr[v + 1]].tolist()]
        explicit.append(numpy.array(entries, dtype=numpy.int32))
    explicit += [numpy.asarray(h, dtype=numpy.int32).reshape(-1, 3) for h in holes]
    info = {'arcs': arcs, 'references': [x - nv for x in fixed if x >= nv],
            'holes': [x - nv for x in cols if x >= nv], 'cuts': [-1 - h for h in hfree],
            'closures': closures, 'shifts': len(gfree), 'pmap': pmap, 'P': P, 'explicit_holes': explicit}
    return M, info


def _exact(values):
    '''! Exact (integer or rational) values of floats that are ratios of small integers'''
    near = numpy.rint(values)
    is_int = numpy.abs(values - near) < 1.0E-9
    return [int(n) if i else Fraction(v).limit_denominator(1000000)
            for v, n, i in zip(values.tolist(), near.tolist(), is_int.tolist())]


def _flux_rows(lc, lc_split, M, cells, a, b, lco):
    '''! Current across each link (out of its first node) as an exact sparse row over the reduced unknowns

    @param lc, lc_split Original and split cell lists
    @param M Corner map from @ref junction_basis
    @param cells, a, b Cell and original edge vertices of each link
    @param lco Split cell list oriented by @ref orient_cells
    @result list of {unknown: coefficient}
    '''
    lc = numpy.asarray(lc)
    lc_split = numpy.asarray(lc_split)
    nk = len(cells)
    pa = numpy.argmax(lc[cells] == a[:, None], axis=1)
    pb = numpy.argmax(lc[cells] == b[:, None], axis=1)
    flipped = numpy.any(numpy.asarray(lco)[cells] != lc_split[cells], axis=1)
    forward = numpy.where(flipped, pa == (pb + 1) % 3, pb == (pa + 1) % 3)
    iy = numpy.where(forward, pb, pa)
    ix = numpy.where(forward, pa, pb)
    D = scipy.sparse.csr_matrix((numpy.concatenate((numpy.ones(nk), -numpy.ones(nk))),
                                 (numpy.tile(numpy.arange(nk), 2), numpy.concatenate((3*cells + iy, 3*cells + ix)))),
                                shape=(nk, M.shape[0]))
    F = (D @ M.tocsr()).tocoo()
    keep = numpy.abs(F.data) > 1.0E-12
    frow, fcol, fval = F.row[keep], F.col[keep], F.data[keep]
    order = numpy.lexsort((fcol, frow))
    frow, fcol, fval = frow[order], fcol[order].tolist(), _exact(fval[order])
    ptr = numpy.searchsorted(frow, numpy.arange(nk + 1)).tolist()
    return [dict(zip(fcol[ptr[k]:ptr[k + 1]], fval[ptr[k]:ptr[k + 1]])) for k in range(nk)]


def _cycle_holes(lc, lc_split, links, M, lco):
    '''! Hole elements for the currents the reduced unknowns do not represent

    A spanning forest of @ref dual_graph is built; each link outside it closes one cycle, and the cycles
    span all currents. Links are taken in turn and kept as generators when their current, in terms of the
    reduced unknowns, depends on the links already taken (greedy row elimination); each generator's cycle
    becomes a hole element: a jump of chi at the vertex where each cell on the cycle is entered and left.

    @param lc, lc_split Original and split cell lists
    @param links Links of the dual graph from @ref _dual_links
    @param M Corner map from @ref junction_basis
    @param lco Split cell list oriented by @ref orient_cells
    @result list of [n,3] (cell, split vertex, sign) arrays
    '''
    n, lu, lw, la, lb, lcell = links
    nc = len(lc)
    nlink = len(lu)
    link_id = numpy.arange(1, nlink + 1)
    A = scipy.sparse.csr_matrix((numpy.concatenate((link_id, link_id)),
                                 (numpy.concatenate((lu, lw)), numpy.concatenate((lw, lu)))), shape=(n, n))
    _, pred = _spanning_forest(A)
    child = numpy.flatnonzero(pred >= 0)
    in_tree = numpy.zeros(nlink, dtype=bool)
    in_tree[numpy.asarray(A[child, pred[child]]).ravel() - 1] = True
    nontree = numpy.flatnonzero(~in_tree)
    rows = _flux_rows(lc, lc_split, M, lcell[nontree], la[nontree], lb[nontree], lco)
    basis, generators = {}, []
    for k, r in zip(nontree.tolist(), rows):
        if not _add_independent(r, basis):
            generators.append(k)
    pred = pred.tolist()

    def path_up(v):
        '''! Path from node v up to the root of its spanning tree'''
        out = [v]
        while pred[out[-1]] >= 0:
            out.append(pred[out[-1]])
        return out
    holes = []
    for k in generators:
        u, w = int(lu[k]), int(lw[k])
        pu, pw = path_up(u), path_up(w)
        common = set(pu) & set(pw)
        lca = next(x for x in pu if x in common)
        seq = [u] + pw[:pw.index(lca) + 1] + pu[1:pu.index(lca)][::-1]
        m = len(seq)
        nxt = seq[1:] + seq[:1]
        ids = numpy.asarray(A[seq, nxt]).ravel() - 1
        edge_of = {}
        for a_, b_, l in zip(seq, nxt, ids.tolist()):
            edge_of[(a_, b_)] = (int(la[l]), int(lb[l]))
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


def find_junction_holes(r, lc, reg, insulated=(), tol=1.0E-10, split_regions=False):
    '''! Split a junction mesh and choose its holes and references

    @param r Point list [np,3]
    @param lc Cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @param insulated Insulation: regions (insulated at every junction they touch) and/or pairs of regions
    (insulated from each other where they are the only regions on a junction edge)
    @param tol Relative eigenvalue below which the model is taken to have a zero-energy mode
    @param split_regions Also split the mesh at all edges where two regions meet (so that those regions can be
    insulated from each other), not only at edges shared by more than two cells?
    @result Dictionary: 'r', 'lc' (split mesh), 'reg', 'lines', 'jedges', 'kje', 'lje' (junction table),
            'vmap', 'holes' (explicit holes, [n,3] (cell, split vertex, sign) arrays: non-reference sheet
            boundaries, then cycle holes), 'cuts' (the cycle holes), 'pmap', 'P' (sparse [n_expanded, n_reduced];
            expanded unknowns are the active vertices then the explicit holes), 'n_unknowns', 'n_loops', 'M', 'info',
            'insulated', 'split_regions'
    '''
    lc = numpy.asarray(lc, dtype=numpy.int32)
    reg = numpy.asarray(reg, dtype=numpy.int32)
    tab = _edge_table(lc)
    jmask = _junction_mask(tab, reg, split_regions)
    rs, lcs, vmap = _split(r, lc, tab, jmask)
    jedges = _junction_edges(tab, jmask)
    lines = _junction_lines(jedges)
    _check_insulation(tab, jmask, reg, insulated)
    links = _dual_links(lc, reg, insulated, tab, jmask)
    target = _loop_count(*links[:3])
    lco = orient_cells(lcs)
    M, info = junction_basis(lc, lcs, reg, jedges, insulated, lco=lco)
    cuts = []
    if M.shape[1] < target:
        cuts = _cycle_holes(lc, lcs, links, M, lco)
        M, info = junction_basis(lc, lcs, reg, jedges, insulated, cuts, lco=lco)
    if M.shape[1] != target:
        raise ValueError('Found {0} unknowns but the mesh carries {1} independent currents'.format(M.shape[1], target))
    Kc = _local_stiffness(rs, lcs)
    nc = len(lcs)
    rows = numpy.repeat(3*numpy.arange(nc)[:, None] + numpy.arange(3)[None, :], 3, axis=1).ravel()
    cols = numpy.tile(3*numpy.arange(nc)[:, None] + numpy.arange(3)[None, :], (1, 3)).ravel()
    K = scipy.sparse.csr_matrix((Kc.ravel(), (rows, cols)), shape=(3*nc, 3*nc))
    if _min_energy(M, K) <= tol:
        raise ValueError('Junction model has a zero-energy mode')
    kje, lje = junction_table(lc, lcs, jedges, lco=lco)
    return {'r': rs, 'lc': lcs, 'reg': reg, 'lines': lines, 'jedges': jedges, 'kje': kje, 'lje': lje, 'vmap': vmap,
            'holes': info['explicit_holes'], 'cuts': cuts, 'pmap': info['pmap'], 'P': info['P'],
            'n_unknowns': M.shape[1], 'n_loops': target, 'M': M, 'info': info, 'insulated': list(insulated),
            'split_regions': bool(split_regions)}


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


def split_jumper(lc, lc_split, vmap, jedges, jumper, ec=None):
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
    @param ec Cells on each edge of the unsplit mesh from @ref _edge_cells (computed if not given)
    @result Ordered split vertices, segments separated by -1 (a closed single segment repeats its first vertex)
    '''
    verts = set(int(v) for v in jumper)
    ec = _edge_cells(lc) if ec is None else ec
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
    number of independent currents, with the insulation and region splitting stored in the file.

    @param r Split point list [np,3]
    @param lc Split cell list [nc,3] (0-based)
    @param reg Region list [nc]
    @param junctions Junction data from @ref meshing.read_ThinCurr_mesh
    @result Number of junction edges, number of reduced unknowns
    '''
    vmap = junctions['vmap']
    split_regions = junctions.get('split_regions', True)
    r0 = numpy.zeros((vmap.max() + 1, 3))
    r0[vmap] = r
    lc0 = vmap[lc]
    _, lcs, jedges, _, _ = split_junctions(r0, lc0, reg, split_regions=split_regions)
    if (lcs.shape != lc.shape) or numpy.any(lcs != lc):
        raise ValueError('Split mesh does not match the junction split of its unsplit mesh, rerun "OFT_ThinCurr_holes"')
    njedges = len(junctions['kje']) - 1
    if njedges != len(jedges):
        raise ValueError('Incorrect number of junction edges found {0} in file (expected {1})'.format(njedges, len(jedges)))
    nunknowns = junctions['P'].shape[1]
    nloops = current_loop_count(lc0, reg, junctions.get('insulated', []), split_regions=split_regions)
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
    junctions = mesh['thincurr']['junctions']
    vmap = junctions['vmap']
    r = numpy.zeros((vmap.max() + 1, 3))
    r[vmap] = mesh['r']
    out = find_junction_holes(r, vmap[mesh['lc']], mesh['reg'], insulated=insulated,
                              split_regions=junctions.get('split_regions', True))
    if (out['lc'].shape != mesh['lc'].shape) or numpy.any(out['lc'] != mesh['lc']):
        raise ValueError('Split mesh in "{0}" does not match the junction split of its unsplit mesh'.format(mesh_file))
    shutil.copyfile(mesh_file, out_file)
    with h5py.File(out_file, 'r+') as h5_file:
        del h5_file['thincurr/junctions']
        write_junction_group(h5_file, out)
    return out
