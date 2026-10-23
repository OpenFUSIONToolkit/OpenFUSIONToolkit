#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
r'''! Hole and closure selection for ThinCurr surface meshes using a tree-cotree approach

ThinCurr represents the current on each surface by a stream function on the vertices, which cannot carry
net current around the non-contractible loops of the surface. A "hole" element (a closed vertex loop across
which the stream function jumps) is added for each independent such loop and a "closure" (one vertex fixed
at zero) is added for each closed surface. For a connected surface with Euler characteristic
\f$ \chi = n_v - n_e + n_f \f$ the number of holes is \f$ 1 - \chi \f$ if the surface has a boundary and
\f$ 2 - \chi \f$ if it is closed.

Holes are chosen for each connected surface as follows:
 1. All boundary loops except the longest are taken as holes.
 2. A tree-cotree decomposition (a spanning tree of the vertices and a spanning tree of the cells over the
    remaining edges) gives one integer functional per hole on closed edge chains. A loop's values under
    these functionals (its signature) are zero exactly when the loop is the boundary of a set of cells, so
    a set of loops is a valid set of holes exactly when its signatures are linearly independent.
 3. Candidate loops are the fundamental cycles of a shortest-path tree from a base vertex. The candidates are
    taken shortest first and kept when their signature is independent of the loops already kept.
 4. Optionally, the shortest-path search is repeated from the middle of each selected loop and the selection
    is repeated over all candidates found ("optimization").
 5. Chords (mesh edges joining two non-consecutive vertices of a loop) are removed by splitting the loop at
    the chord and keeping the part that leaves the set of holes independent, so each hole is a simple,
    chord-free loop.

The cost is linear in the number of cells for fixed topology, plus a term proportional to the number of
holes times the number of cells crossed by the cuts.

@authors Chris Hansen
@date October 2026
@ingroup doxy_oft_python
'''
from fractions import Fraction
import numpy as np
import scipy.sparse
import scipy.sparse.csgraph
import scipy.sparse.linalg

_block_size = 32 # Number of hole functionals processed together when computing loop signatures
_max_passes = 10 # Maximum number of search passes when optimizing


def _edge_table(lc, nv):
    '''! Build the unique edges of a triangle mesh and the half-edge to edge mapping

    Half-edge `3*i+j` runs from `lc[i,j]` to `lc[i,(j+1)%3]`.

    @param lc Cell list [nf,3] (0-based)
    @param nv Number of vertices
    @result `le` edge list [ne,2] (low, high vertex), `he_edge` edge of each half-edge [3*nf],
    `he_dir` direction of each half-edge (+1 if low -> high) [3*nf], `counts` number of cells on each edge [ne]
    '''
    a = lc.reshape(-1).astype(np.int64)
    b = lc[:,[1,2,0]].reshape(-1).astype(np.int64)
    if np.any(a == b):
        raise ValueError('Degenerate cell(s) with repeated vertices found in mesh')
    keys, he_edge, counts = np.unique(np.minimum(a,b)*nv + np.maximum(a,b), return_inverse=True, return_counts=True)
    le = np.stack((keys // nv, keys % nv), axis=1)
    he_dir = np.where(a < b, 1, -1)
    return le, he_edge.reshape(-1), he_dir, counts


def _edge_cells(he_edge, counts):
    '''! Half-edges on each edge

    @param he_edge Edge of each half-edge
    @param counts Number of half-edges on each edge
    @result Half-edges on each edge [ne,2] (-1 for missing second half-edge)
    '''
    order = np.argsort(he_edge, kind='stable')
    start = np.concatenate(([0], np.cumsum(counts)[:-1]))
    ehe = -np.ones((counts.shape[0],2), dtype=np.int64)
    ehe[:,0] = order[start]
    two = np.nonzero(counts > 1)[0]
    ehe[two,1] = order[start[two]+1]
    return ehe


def _edge_index(keys, nv, u, v):
    '''! Index of the edges joining vertices `u` and `v` (arrays)'''
    return np.searchsorted(keys, np.minimum(u,v).astype(np.int64)*nv + np.maximum(u,v))


def _path_sums(pred, vals):
    '''! Sum values along the tree path from each node to its root by pointer jumping

    @param pred Parent of each node (roots are their own parent)
    @param vals Value of the link from each node to its parent [n,...] (zero for roots)
    @result Sum of `vals` over the path from each node to its root [n,...]
    '''
    P = pred.copy()
    X = vals.copy()
    while True:
        PP = P[P]
        if np.array_equal(PP, P):
            return X
        X += X[P]
        P = PP


def _tree_lca(pred, a, b):
    '''! Lowest common ancestors of node pairs in a tree (binary lifting)

    @param pred Parent of each node (the root is its own parent)
    @param a First node of each pair
    @param b Second node of each pair
    @result Lowest common ancestor of each pair
    '''
    depth = _path_sums(pred, (pred != np.arange(pred.shape[0])).astype(np.int64))
    up = [pred]
    while (1 << len(up)) <= depth.max():
        up.append(up[-1][up[-1]])
    a = a.copy()
    b = b.copy()
    swap = depth[a] < depth[b]
    a[swap], b[swap] = b[swap], a[swap].copy()
    diff = depth[a] - depth[b]
    for k, upk in enumerate(up):
        move = ((diff >> k) & 1) == 1
        a[move] = upk[a[move]]
    for upk in reversed(up):
        move = upk[a] != upk[b]
        a[move] = upk[a[move]]
        b[move] = upk[b[move]]
    return np.where(a == b, a, pred[a])


def _check_mesh(r, lc):
    '''! Check a surface mesh, orient its cells consistently and find its connected surfaces

    @param r Vertex list [nv,3]
    @param lc Cell list [nf,3] (0-based)
    @result Oriented cell list [nf,3], surface index of each cell [nf]
    '''
    nv = r.shape[0]
    nf = lc.shape[0]
    if (lc.ndim != 2) or (lc.shape[1] != 3):
        raise ValueError("Cell list must have shape [nf, 3]")
    if (lc.min() < 0) or (lc.max() >= nv):
        raise ValueError("Cell list references vertices that do not exist")
    le, he_edge, he_dir, counts = _edge_table(lc, nv)
    if counts.max() > 2:
        print('ERROR: Non-manifold surface detected')
        for k in np.nonzero(counts > 2)[0][:20]:
            pt = (r[le[k,0],:] + r[le[k,1],:])/2.0
            print('  {0} cells share edge {1} located at {2} {3} {4}'.format(counts[k], k+1, *pt))
        raise ValueError("ThinCurr requires manifold surfaces")
    # Connected surfaces (through shared vertices)
    G = scipy.sparse.coo_matrix((np.ones(le.shape[0]), (le[:,0], le[:,1])), shape=(nv,nv))
    _, vsurf = scipy.sparse.csgraph.connected_components(G, directed=False)
    fsurf = vsurf[lc[:,0]]
    # Orient cells: neighbors through an edge must traverse it in opposite directions
    ehe = _edge_cells(he_edge, counts)
    inner = np.nonzero(counts == 2)[0]
    h1 = ehe[inner,0]
    h2 = ehe[inner,1]
    f1 = h1 // 3
    f2 = h2 // 3
    same = (he_dir[h1] == he_dir[h2]).astype(np.int64)
    _, fgroup = scipy.sparse.csgraph.connected_components(
        scipy.sparse.coo_matrix((np.ones(inner.shape[0]), (f1, f2)), shape=(nf,nf)), directed=False)
    _, first = np.unique(fgroup, return_index=True)
    if np.unique(vsurf[lc[first,0]]).shape[0] != first.shape[0]:
        raise ValueError("Surfaces that touch only at a vertex are not supported")
    rows = np.concatenate((f1, np.full(first.shape[0], nf)))
    cols = np.concatenate((f2, first))
    Gf = scipy.sparse.coo_matrix((np.ones(rows.shape[0]), (rows, cols)), shape=(nf+1,nf+1)).tocsr()
    _, fpred = scipy.sparse.csgraph.breadth_first_order(Gf, nf, directed=False, return_predecessors=True)
    fpred = fpred[:nf]
    pkeys = np.minimum(f1,f2)*nf + np.maximum(f1,f2)
    psort = np.argsort(pkeys)
    pred = np.arange(nf)
    rel = np.zeros(nf, dtype=np.int64)
    linked = np.nonzero(fpred < nf)[0]
    pred[linked] = fpred[linked]
    lk = np.minimum(linked,pred[linked])*nf + np.maximum(linked,pred[linked])
    rel[linked] = same[psort[np.searchsorted(pkeys[psort], lk)]]
    flip = _path_sums(pred, rel) % 2
    sign = 1 - 2*flip
    if np.any(he_dir[h1]*sign[f1] == he_dir[h2]*sign[f2]):
        raise ValueError("ThinCurr requires orientable surfaces (cells could not be oriented consistently)")
    lco = lc.copy()
    lco[flip == 1,1] = lc[flip == 1,2]
    lco[flip == 1,2] = lc[flip == 1,1]
    return lco, fsurf


class _surface:
    '''! Topology and hole functionals for one connected, consistently oriented surface'''
    def __init__(self, r, lc):
        '''! Analyze surface

        @param r Vertex list [nv,3]
        @param lc Oriented cell list [nf,3] (0-based)
        '''
        self.r = r
        self.lc = lc
        self.nv = r.shape[0]
        self.nf = lc.shape[0]
        self.le, he_edge, self.he_dir, counts = _edge_table(lc, self.nv)
        self.ne = self.le.shape[0]
        self.keys = self.le[:,0]*self.nv + self.le[:,1]
        self.ehe = _edge_cells(he_edge, counts)
        self.elen = np.linalg.norm(r[self.le[:,1],:]-r[self.le[:,0],:], axis=1)
        G = scipy.sparse.coo_matrix((self.elen, (self.le[:,0], self.le[:,1])), shape=(self.nv,self.nv)).tocsr()
        self.graph = (G + G.T).tocsr()
        self.closed = counts.min() > 1
        euler = self.nv - self.ne + self.nf
        self.nholes = (2 - euler) if self.closed else (1 - euler)
        self.boundary_loops = self._boundary_loops(np.nonzero(counts == 1)[0])

    def _boundary_loops(self, bedges):
        '''! Trace boundary loops (oriented with the cells)'''
        hb = self.ehe[bedges,0]
        src = self.lc.reshape(-1)[hb]
        dst = self.lc[:,[1,2,0]].reshape(-1)[hb]
        if np.any(np.bincount(src, minlength=self.nv) > 1):
            raise ValueError("Boundary touches itself at a vertex (pinched boundary), which is not supported")
        nxt = -np.ones(self.nv, dtype=np.int64)
        nxt[src] = dst
        seen = np.zeros(self.nv, dtype=bool)
        loops = []
        for start in src:
            if seen[start]:
                continue
            loop = [int(start)]
            seen[start] = True
            v = nxt[start]
            while v != start:
                loop.append(int(v))
                seen[v] = True
                v = nxt[v]
            loops.append(loop)
        return loops

    def setup_functionals(self, root):
        '''! Build the hole functionals from a tree-cotree decomposition

        A spanning tree of the vertices (shortest paths from `root`) is built and the remaining edges link
        the cells (and, for boundary edges, a "ground" node). Each edge outside a spanning tree of these
        links closes a cycle, which defines one functional: the signed sum of a loop's edges over the cycle.

        @param root Root vertex of the vertex spanning tree
        '''
        nv, nf = self.nv, self.nf
        _, pred = scipy.sparse.csgraph.dijkstra(self.graph, indices=root, return_predecessors=True)
        tree_v = np.nonzero(pred >= 0)[0]
        if tree_v.shape[0] != nv-1:
            raise RuntimeError("Surface vertices are not connected")
        in_tree = np.zeros(self.ne, dtype=bool)
        in_tree[_edge_index(self.keys, nv, tree_v, pred[tree_v])] = True
        # Links between cells (and ground) over edges not in the vertex tree, oriented from the cell that
        # traverses the edge high -> low to the cell that traverses it low -> high
        nt = np.nonzero(~in_tree)[0]
        h1 = self.ehe[nt,0]
        h2 = self.ehe[nt,1]
        ground = nf
        f2 = np.where(h2 >= 0, h2 // 3, ground)
        d1 = self.he_dir[h1]
        tail = np.where(d1 < 0, h1 // 3, f2)
        head = np.where(d1 < 0, f2, h1 // 3)
        nH = nf + (0 if self.closed else 1)
        H = scipy.sparse.coo_matrix((np.ones(nt.shape[0]), (tail, head)), shape=(nH,nH)).tocsr()
        order, hpred = scipy.sparse.csgraph.breadth_first_order(H, 0, directed=False, return_predecessors=True)
        if order.shape[0] != nH:
            raise RuntimeError("Cells are not connected after removing the vertex tree")
        # Link from each node to its parent in the cotree
        tlink = -np.ones(nH, dtype=np.int64)
        for child, parent in ((head, tail), (tail, head)):
            idx = np.nonzero(hpred[child] == parent)[0]
            nodes, first = np.unique(child[idx], return_index=True)
            keep = tlink[nodes] < 0
            tlink[nodes[keep]] = idx[first[keep]]
        if np.count_nonzero(tlink >= 0) != nH-1:
            raise RuntimeError("Failed to build cotree")
        in_cotree = np.zeros(nt.shape[0], dtype=bool)
        in_cotree[tlink[tlink >= 0]] = True
        leftover = np.nonzero(~in_cotree)[0]
        if leftover.shape[0] != self.nholes:
            raise RuntimeError("Found {0} independent loops but expected {1} from the Euler characteristic".format(leftover.shape[0], self.nholes))
        hp = np.where(hpred >= 0, hpred, np.arange(nH))
        depth = _path_sums(hp, (hpred >= 0).astype(np.int64))
        rows, cols, vals = [], [], []
        for k, l in enumerate(leftover):
            cyc = [(nt[l], 1)]
            u, v = head[l], tail[l] # Cycle: tail -> head, head -> ... -> common parent -> ... -> tail
            while u != v:
                if depth[u] >= depth[v]:
                    cyc.append((nt[tlink[u]], 1 if tail[tlink[u]] == u else -1))
                    u = hpred[u]
                else:
                    cyc.append((nt[tlink[v]], 1 if tail[tlink[v]] == hpred[v] else -1))
                    v = hpred[v]
            for e, s in cyc:
                rows.append(k)
                cols.append(e)
                vals.append(s)
        self.Y = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(self.nholes,self.ne), dtype=np.int64)
        self.YT = self.Y.T.tocsr()

    def loop_signature(self, loop):
        '''! Values of the hole functionals for a closed vertex loop (dense integer vector)'''
        loop = np.asarray(loop, dtype=np.int64)
        nxt = np.roll(loop, -1)
        edges = _edge_index(self.keys, self.nv, loop, nxt)
        signs = np.where(loop < nxt, 1, -1)
        return np.asarray(self.YT[edges].multiply(signs[:,None]).sum(axis=0)).reshape(-1)

    def tree(self, root, limit=np.inf):
        '''! Shortest-path tree from a root vertex

        @result distance and parent of each vertex (-9999 if not reached or root)
        '''
        return scipy.sparse.csgraph.dijkstra(self.graph, indices=root, return_predecessors=True, limit=limit)

    def candidates(self, root, limit=np.inf):
        '''! Fundamental cycles of the shortest-path tree from `root` that are not cell boundaries

        Each edge (a,b) outside the tree closes the loop a -> b -> (tree path) -> a. Loops are ranked by their
        length, excluding the part of the tree paths shared by a and b.

        @param root Root vertex
        @param limit Only search vertices closer than this to the root and keep loops shorter than this
        @result Loop lengths [n], edges [n], signatures (sparse [n,nholes])
        '''
        nv, le = self.nv, self.le
        dist, pred = self.tree(root, limit)
        reach = np.isfinite(dist)
        nodes = np.nonzero(reach)[0]
        loc = -np.ones(nv, dtype=np.int64)
        loc[nodes] = np.arange(nodes.shape[0])
        parent = pred[nodes]
        child = np.nonzero(parent >= 0)[0]
        tpred = np.arange(nodes.shape[0])
        tpred[child] = loc[parent[child]]
        tedge = _edge_index(self.keys, nv, nodes[child], parent[child])
        tdir = np.where(parent[child] < nodes[child], 1, -1) # +1 if the step parent -> child runs low -> high
        is_tree = np.zeros(self.ne, dtype=bool)
        is_tree[tedge] = True
        reach_edge = reach[le[:,0]] & reach[le[:,1]]
        ce = np.nonzero(reach_edge & ~is_tree)[0]
        # Functionals that can be nonzero on loops in the reached region
        active = np.unique(self.Y[:,np.nonzero(reach_edge)[0]].nonzero()[0])
        la = loc[le[ce,0]]
        lb = loc[le[ce,1]]
        rows, cols, vals = [], [], []
        for i in range(0, active.shape[0], _block_size):
            block = active[i:i+_block_size]
            YTb = self.YT[:,block]
            X = np.zeros((nodes.shape[0], block.shape[0]), dtype=np.int64)
            X[child,:] = YTb[tedge].toarray()*tdir[:,None]
            X = _path_sums(tpred, X)
            S = X[la,:] - X[lb,:] + YTb[ce].toarray()
            nzr, nzc = np.nonzero(S)
            rows.append(nzr)
            cols.append(block[nzc])
            vals.append(S[nzr,nzc])
        if len(rows) > 0:
            rows = np.concatenate(rows)
            cols = np.concatenate(cols)
            vals = np.concatenate(vals)
        sigs = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(ce.shape[0],self.nholes), dtype=np.int64)
        nonzero = np.nonzero(np.diff(sigs.indptr) > 0)[0]
        ce = ce[nonzero]
        lca = _tree_lca(tpred, la[nonzero], lb[nonzero])
        clen = dist[le[ce,0]] + dist[le[ce,1]] + self.elen[ce] - 2.0*dist[nodes[lca]]
        keep = np.nonzero(clen <= limit)[0]
        return clen[keep], ce[keep], sigs[nonzero[keep]]

    def tree_cycle(self, root, limit, edge):
        '''! Simple vertex loop closed by `edge` in the shortest-path tree from `root` (without the stem to the root)'''
        _, pred = self.tree(root, limit)
        a, b = int(self.le[edge,0]), int(self.le[edge,1])
        pa = [a]
        while pred[pa[-1]] >= 0:
            pa.append(int(pred[pa[-1]]))
        pos = {v: i for i, v in enumerate(pa)}
        pb = [b]
        while pb[-1] not in pos:
            pb.append(int(pred[pb[-1]]))
        return pa[:pos[pb[-1]]+1] + pb[:-1][::-1]

    def tighten(self, loop):
        '''! Shorten a loop by cutting corners: where three consecutive vertices form a cell, the middle
        vertex is removed (moving the loop across the cell, which does not change its signature)'''
        if not hasattr(self, 'face_set'):
            self.face_set = set(map(tuple, np.sort(self.lc, axis=1).tolist()))
        loop = list(loop)
        i = 0
        while (len(loop) > 3) and (i < len(loop)):
            n = len(loop)
            tri = tuple(sorted((loop[i-1], loop[i], loop[(i+1) % n])))
            if tri in self.face_set:
                del loop[i]
                i = max(i-1, 0)
            else:
                i += 1
        return loop

    def find_chord(self, loop):
        '''! Positions (i,j), i < j, of two non-consecutive loop vertices joined by a mesh edge, or None'''
        pos = {v: i for i, v in enumerate(loop)}
        n = len(loop)
        indptr, indices = self.graph.indptr, self.graph.indices
        for i, v in enumerate(loop):
            for u in indices[indptr[v]:indptr[v+1]]:
                j = pos.get(int(u))
                if (j is not None) and (j > i+1) and not ((i == 0) and (j == n-1)):
                    return i, j
        return None


class _independent_set:
    '''! Incremental test for linear independence of sparse integer vectors (exact arithmetic)'''
    def __init__(self):
        self.basis = {}

    def add(self, vec):
        '''! Add `vec` ({index: value}) if it is independent of the vectors already added

        @result True if added
        '''
        r = {k: Fraction(int(v)) for k, v in vec.items() if v != 0}
        basis = self.basis
        while r:
            p = max(r)
            if p not in basis:
                break
            f = r[p]
            for k, v in basis[p].items():
                val = r.get(k, 0) - f*v
                if val == 0:
                    r.pop(k, None)
                else:
                    r[k] = val
        if not r:
            return False
        p = max(r)
        basis[p] = {k: v/r[p] for k, v in r.items()}
        return True


def _dense_to_dict(vec):
    return {int(k): int(vec[k]) for k in np.nonzero(vec)[0]}


def _sparse_row(mat, i):
    r0, r1 = mat.indptr[i], mat.indptr[i+1]
    return dict(zip(mat.indices[r0:r1].tolist(), mat.data[r0:r1].tolist()))


def _select(seeds, pools, nholes):
    '''! Greedily select loops shortest first

    @param seeds Signatures of loops that must be kept (dense)
    @param pools List of (root, limit, lengths, edges, signatures) candidate sets
    @param nholes Required number of independent loops
    @result List of (length, root, limit, edge) for the selected candidates
    '''
    indep = _independent_set()
    for sig in seeds:
        if not indep.add(_dense_to_dict(sig)):
            raise RuntimeError("Boundary loops are not independent")
    nsel = len(seeds)
    lengths = np.concatenate([p[2] for p in pools])
    pool_ind = np.concatenate([np.full(p[2].shape[0], i) for i, p in enumerate(pools)])
    row_ind = np.concatenate([np.arange(p[2].shape[0]) for p in pools])
    selected = []
    for k in np.lexsort((row_ind, pool_ind, lengths)):
        if nsel == nholes:
            break
        root, limit, _, edges, sigs = pools[pool_ind[k]]
        if indep.add(_sparse_row(sigs, row_ind[k])):
            selected.append((lengths[k], root, limit, edges[row_ind[k]]))
            nsel += 1
    if nsel != nholes:
        raise RuntimeError("Only found {0} of {1} independent loops".format(nsel, nholes))
    return selected


def _clean_loops(surf, loops, sigs, first):
    '''! Tighten loops and remove chords while keeping the set of loops independent

    @param surf Surface
    @param loops Loops to process (modified in place)
    @param sigs Signatures of all holes on the surface (dense [nholes,nholes], modified in place)
    @param first Row of `sigs` for the first loop in `loops`
    '''
    lu = None
    for j in range(len(loops)):
        while True:
            loops[j] = surf.tighten(loops[j])
            chord = surf.find_chord(loops[j])
            if chord is None:
                break
            i, k = chord
            loop = loops[j]
            parts = sorted((loop[i:k+1], loop[k:] + loop[:i+1]), key=len)
            if lu is None:
                lu = scipy.sparse.linalg.splu(scipy.sparse.csc_matrix(sigs.T.astype(np.float64)))
            for part in parts:
                sig = surf.loop_signature(part)
                coeff = lu.solve(sig.astype(np.float64))
                if abs(coeff[first+j]) > 1.E-8:
                    loops[j] = part
                    sigs[first+j,:] = sig
                    lu = None
                    break
            else:
                raise RuntimeError("Could not remove chord from hole loop")


def _surface_holes(r, lc, ref_point, optimize, info, indent):
    '''! Select holes for one connected surface

    @result boundary holes, omitted boundary loop (or None), internal holes (local vertex indices)
    '''
    surf = _surface(r, lc)
    nholes = surf.nholes
    bloops = surf.boundary_loops
    omitted = None
    if len(bloops) > 0:
        imax = int(np.argmax([len(loop) for loop in bloops]))
        omitted = bloops[imax]
        bloops = bloops[:imax] + bloops[imax+1:]
    if info:
        print(indent+'  # of vertices = {0}, edges = {1}, cells = {2}'.format(surf.nv, surf.ne, surf.nf))
        print(indent+'  # of boundary loops = {0}'.format(len(surf.boundary_loops)))
        print(indent+'  # of holes = {0} ({1} boundary, {2} internal)'.format(nholes, len(bloops), nholes-len(bloops)))
    if nholes == len(bloops):
        return bloops, omitted, []
    root = int(np.argmin(np.linalg.norm(r-ref_point[np.newaxis,:], axis=1)))
    surf.setup_functionals(root)
    seeds = [surf.loop_signature(loop) for loop in bloops]
    pools = [(root, np.inf) + surf.candidates(root)]
    selected = _select(seeds, pools, nholes)
    # Search again from the middle of each selected loop, until the selection does not change. Loops longer than
    # those already selected are never needed, so each search is limited to that length.
    searched = {root}
    for _ in range(_max_passes if optimize else 1):
        limit = max(s[0] for s in selected)
        roots = set()
        for _, sroot, slimit, edge in selected:
            loop = surf.tighten(surf.tree_cycle(sroot, slimit, edge))
            roots.add(loop[len(loop)//2])
        roots = sorted(roots - searched)
        if len(roots) == 0:
            break
        searched.update(roots)
        if info:
            print(indent+'  Searching from {0} additional base point(s)'.format(len(roots)))
        for new_root in roots:
            pools.append((new_root, limit) + surf.candidates(new_root, limit))
        previous = set((s[1], s[3]) for s in selected)
        selected = _select(seeds, pools, nholes)
        if set((s[1], s[3]) for s in selected) == previous:
            break
    loops = [surf.tree_cycle(sroot, slimit, edge) for _, sroot, slimit, edge in selected]
    sigs = np.array(seeds + [surf.loop_signature(loop) for loop in loops], dtype=np.int64).reshape(nholes, nholes)
    _clean_loops(surf, loops, sigs, len(seeds))
    return bloops, omitted, loops


def count_holes(r, lc):
    '''! Count the holes and closures required for a ThinCurr surface mesh (without selecting them)

    @param r Vertex list [nv,3]
    @param lc Cell list [nf,3] (0-based)
    @result Number of holes, number of closures
    '''
    r = np.asarray(r)
    lc = np.asarray(lc, dtype=np.int64)
    _, fsurf = _check_mesh(r, lc)
    le, _, _, counts = _edge_table(lc, r.shape[0])
    surfs, fsurf = np.unique(fsurf, return_inverse=True)
    vsurf = -np.ones(r.shape[0], dtype=np.int64)
    vsurf[lc.reshape(-1)] = np.repeat(fsurf, 3)
    esurf = vsurf[le[:,0]]
    ns = surfs.shape[0]
    euler = np.bincount(vsurf[vsurf >= 0], minlength=ns) - np.bincount(esurf, minlength=ns) + np.bincount(fsurf, minlength=ns)
    closed = np.bincount(esurf[counts == 1], minlength=ns) == 0
    nholes = int(np.sum(np.where(closed, 2-euler, 1-euler)))
    return nholes, int(np.count_nonzero(closed))


def find_holes(r, lc, ref_point=None, optimize=False, info=False, indent=''):
    '''! Select holes and closures for a ThinCurr surface mesh

    @param r Vertex list [nv,3]
    @param lc Cell list [nf,3] (0-based)
    @param ref_point Reference location for the base point of each surface [x,y,z] (default: [0,0,0])
    @param optimize Search again from the middle of each selected loop to shorten the internal holes?
    @param info Print information on each surface?
    @param indent Indentation for printed information
    @result Dictionary with 'boundary_holes' (boundary loops used as holes), 'omitted' (one boundary loop per
    surface with a boundary, not used), 'internal_holes' (other loops) and 'closures' (cell indices, one per
    closed surface); loops are lists of vertex indices (0-based, not repeated at the end)
    '''
    r = np.asarray(r, dtype=np.float64)
    lc = np.asarray(lc, dtype=np.int64)
    ref_point = np.zeros(3) if ref_point is None else np.asarray(ref_point, dtype=np.float64)
    lco, fsurf = _check_mesh(r, lc)
    order = np.argsort(fsurf, kind='stable')
    surfs, starts = np.unique(fsurf[order], return_index=True)
    starts = np.append(starts, order.shape[0])
    out = {'boundary_holes': [], 'omitted': [], 'internal_holes': [], 'closures': []}
    for i in range(surfs.shape[0]):
        faces = order[starts[i]:starts[i+1]]
        verts, lc_local = np.unique(lco[faces,:], return_inverse=True)
        lc_local = lc_local.reshape(-1,3)
        if info:
            print()
            print(indent+'Analyzing surface {0} of {1}'.format(i+1, surfs.shape[0]))
        bloops, omitted, loops = _surface_holes(r[verts,:], lc_local, ref_point, optimize, info, indent)
        out['boundary_holes'] += [verts[loop] for loop in bloops]
        out['internal_holes'] += [verts[loop] for loop in loops]
        if omitted is None:
            out['closures'].append(int(faces.min()))
        else:
            out['omitted'].append(verts[omitted])
    return out
