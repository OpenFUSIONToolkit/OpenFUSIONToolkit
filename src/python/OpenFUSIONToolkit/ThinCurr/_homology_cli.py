#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python cli for ThinCurr hole construction using homology computations

See @ref thincurr_homology_tool

@authors Chris Hansen
@date July 2026
@ingroup doxy_oft_python
'''
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy
from .meshing import read_ThinCurr_mesh, write_ThinCurr_mesh

tri_ed = np.asarray([[2,1], [0,2], [1,0]]) # Triangle edge list
indent_level = ''

class trimesh:
    '''Triangular mesh class for topology analysis
    '''
    def __init__(self, r, lf, info=True):
        '''Construct mesh and orient

        @param r Vertex list [nv,3]
        @param lf Face list [nf,3]
        '''
        if lf.shape[1] != 3:
            raise ValueError("Face list must have shape [nf, 3]")
        self.r=r
        self.lf=lf
        self.np=r.shape[0]
        self.nf=lf.shape[0]
        self.ne = 0
        # Setup mesh
        self._setup_edges()
        self._setup_neighbors()
        self._setup_boundary()
        if len(self._non_manifold_edges) > 0:
            print('ERROR: Non-manifold surface detected')
            for k, cells in self._non_manifold_edges.items():
                pt = (self.r[self.le[k,0],:] + self.r[self.le[k,1],:])/2.0
                print('  {0} cells share edge {1} located at {2} {3} {4}'.format(2+len(cells), k+1, *pt))
            raise ValueError("ThinCurr requires manifold surfaces")
        self.surf_tag = self._orient_surface()
        self.nsurfs = np.max(self.surf_tag)+1
        if info:
            print(indent_level+"Mesh constructed:")
            print(indent_level+"  Found {0} distinct surface(s)".format(self.nsurfs))
            print(indent_level+"  # of vertices = {0} ({1})".format(self.np, self.nbp))
            print(indent_level+"  # of edges    = {0} ({1})".format(self.ne, self.nbe))
            print(indent_level+"  # of faces    = {0} ({1})".format(self.nf, self.nbf))

    def _setup_edges(self):
        '''Construct edges from mesh, building `klpe`, `llpe`, `le`, and `lfe`
        '''
        nr=np.zeros((self.np+1,), dtype=np.int32) # initialize raw edge counter
        ir=np.zeros((self.np+1,), dtype=np.int32)
        jr=np.zeros((self.nf*3,), dtype=np.int32)
        for i in range(self.nf): # loop over cells & count raw edges
            k=self.lf[i,:] # cell corners
            for j in range(3): # loop over edges
                js=min(k[tri_ed[j,0]],k[tri_ed[j,1]]) # low points
                nr[js]+=1                             # edge counter
        ir[self.np]=3*self.nf
        for i in range(self.np-1,-1,-1): # cumulative raw edge count
            ir[i]=ir[i+1]-nr[i]
        if ir[0] != 0:
            raise ValueError('Bad raw edge count')
        nr=np.zeros((self.np+1,), dtype=np.int32) # reset raw edge counter
        for i in range(self.nf): # loop over cells & index raw edges
            k=self.lf[i,:] # cell corners
            for j in range(3): # loop over edges
                js=min(k[tri_ed[j,0]],k[tri_ed[j,1]]) # low points
                je=max(k[tri_ed[j,0]],k[tri_ed[j,1]]) # high points
                jr[ir[js]+nr[js]]=je                  # high point list
                nr[js]+=1                             # edge counter
        for i in range(self.np): # sort raw high point list
            if(nr[i]>1):
                jr[ir[i]:ir[i]+nr[i]].sort()
        self.klpe=np.zeros((self.np+1,), dtype=np.int32)
        for i in range(self.np): # loop over raw edges & count unique edges
            js=ir[i]
            je=ir[i+1]
            for j in range(js,je):
                if j==js:
                    self.klpe[i]=1
                else:
                    if jr[j]>jr[j-1]:
                        self.klpe[i]=self.klpe[i]+1
        self.ne=np.sum(self.klpe) # total number of unique edges
        self.klpe[self.np]=self.ne
        for i in range(self.np-1,-1,-1): # cumulative unique edge count
            self.klpe[i]=self.klpe[i+1]-self.klpe[i]
        if self.klpe[0]!=0:
            raise ValueError('Bad unique edge count')
        self.llpe=(self.np+1)*np.ones((self.ne,), dtype=np.int32)
        for i in range(self.np): # loop over raw edges & index unique edges
            js=ir[i]
            je=ir[i+1]
            for j in range(js,je):
                if j==js:
                    jn=self.klpe[i]
                    self.llpe[jn]=jr[j]
                else:
                    if jr[j]>jr[j-1]:
                        jn=jn+1
                        self.llpe[jn]=jr[j]
        self.le=np.zeros((self.ne,2), dtype=np.int32)
        for i in range(self.np): # construct global edge list
            self.le[self.klpe[i]:self.klpe[i+1],0]=i
        for i in range(self.ne):
            self.le[i,1]=self.llpe[i]
        self.lfe=np.zeros((self.nf,3), dtype=np.int32)
        for i in range(self.nf): # loop over cells & index cells to edges
            k=self.lf[i,:] # cell corners
            for j in range(3): # loop over edges
                js=min(k[tri_ed[j,0]],k[tri_ed[j,1]]) # low points
                je=max(k[tri_ed[j,0]],k[tri_ed[j,1]]) # high points
                jp=self.klpe[js]                      # pointer into low point list
                jn=self.klpe[js+1]-jp                 # number of shared edges
                self.lfe[i,j]=np.searchsorted(self.llpe[jp:jp+jn],je)
                if self.llpe[jp+self.lfe[i,j]] != je:
                    print(self.llpe[jp:jp+jn],je)
                    raise ValueError('Face edge not found!')
                self.lfe[i,j] += jp+1
                if k[tri_ed[j,1]]-k[tri_ed[j,0]]<0:
                    self.lfe[i,j] *= -1 # apply orientation

    def _setup_neighbors(self):
        '''Build topology neighbor lists `kpf`, `lpf`, `lef`, and `lff`
        '''
        self.kpf=np.zeros((self.np+1,), dtype=np.int32)
        nr=np.zeros((self.np+1,), dtype=np.int32)
        for i in range(self.nf): # loop over cells
            for j in range(3):   # loop over corners
                k=self.lf[i,j]   # corner number
                nr[k]+=1         # count cell to corner
        self.npf=np.sum(nr)
        self.kpf[self.np]=self.npf
        for i in range(self.np-1,-1,-1): # cumulative point to cell count
            self.kpf[i]=self.kpf[i+1]-nr[i]
        if self.kpf[0]!=0:
            raise ValueError('Bad point to cell count')
        self.lpf=np.zeros((self.npf,), dtype=np.int32)
        nr=np.zeros((self.np+1,), dtype=np.int32)
        for i in range(self.nf):              # loop over cells
            for j in range(3):                # loop over corners
                k=self.lf[i,j]                # corner number
                self.lpf[self.kpf[k]+nr[k]]=i # index cell to corner
                nr[k]+=1                      # count cell to corner
        self.lef=-np.ones((self.ne,2), dtype=np.int32)
        self.lff=np.zeros((self.nf,3), dtype=np.int32)
        self._non_manifold_edges = {}
        for i in range(self.nf):      # loop over cells & index to edges
            for j in range(3):        # loop over edges
                k=abs(self.lfe[i,j])-1  # edge numbers
                if self.lef[k,0]<0:
                    self.lef[k,0]=i   # record first cell
                else:
                    if self.lef[k,1] >= 0: # Handle more than 2 connections
                        curr_list = self._non_manifold_edges.get(k,[])
                        if i not in curr_list:
                            self._non_manifold_edges[k] = curr_list + [i]
                    else:
                        self.lef[k,1]=i   # record second cell
        for i in range(self.nf):     # loop over cells & locate neighbors
            for j in range(3):       # loop over edges
                k=abs(self.lfe[i,j])-1 # edge numbers
                self.lff[i,j]=np.sum(self.lef[k,:])-i

    def _setup_boundary(self):
        '''Locate and mark boundary elements
        '''
        self.bp=np.zeros((self.np,), dtype=np.bool_)
        self.be=np.zeros((self.ne,), dtype=np.bool_)
        self.bf=np.zeros((self.nf,), dtype=np.bool_)
        # boundary cells have 1 or more faces on boundary
        self.be=np.any(self.lef<0,axis=1)
        self.bf=np.any(self.lff<0,axis=1)
        self.nbe=np.sum(self.be)
        self.nbf=np.sum(self.bf)
        self.lbe=np.zeros((self.nbe,), dtype=np.int32)
        self.lbf=np.zeros((self.nbf,), dtype=np.int32)
        j=0
        for i in range(self.ne):
            if self.be[i]:
                self.lbe[j]=i
                j+=1
        j=0
        for i in range(self.nf):
            if self.bf[i]:
                self.lbf[j]=i
                j+=1
        for i in range(self.nbe):
            for j in range(2):
                self.bp[self.le[self.lbe[i],j]]=True
        self.nbp=np.sum(self.bp)
        self.lbp=np.zeros((self.nbp,), dtype=np.int32)
        j=0
        for i in range(self.np):
            if self.bp[i]:
                self.lbp[j]=i
                j+=1

    def _orient_surface(self):
        '''Orient surface(s) in mesh to ensure consistency
        '''
        def orient_neighbors(starting_face,oriented):
            stack = [starting_face]
            while len(stack) > 0:
                face = stack.pop()
                for j in range(3):
                    face2 = self.lff[face,j]
                    if face2 < 0:
                        continue
                    if oriented[face2] >= 0:
                        continue
                    ed = self.lf[face,tri_ed[j,:]]
                    # Ensure same sense as neighbor (opposite direction of shared edge)
                    for k in range(3):
                        if(self.lff[face2,k]==face):
                            break
                    else:
                        raise ValueError("Could not find face!!")
                    ed2 = self.lf[face2,tri_ed[k,:]]
                    if (ed[0]==ed2[0]) and (ed[1]==ed2[1]):
                        self.lf[face2,[1,2]] = self.lf[face2,[2,1]]
                        self.lfe[face2,:] = -self.lfe[face2,[0,2,1]]
                        self.lff[face2,:] = self.lff[face2,[0,2,1]]
                    oriented[face2] = surf_id
                    stack.append(face2)
        #
        oriented = [-1 for _ in range(self.nf)]
        surf_id=-1
        for i in range(self.nf):
            if oriented[i]<0:
                surf_id += 1
                oriented[i] = surf_id
                orient_neighbors(i,oriented)
        return np.array(oriented)

    def get_face_edge_bop(self):
        r''' Compute face to edge boundary operator \partial_2
        '''
        I = np.zeros((self.nf*3,),dtype=np.int32)
        I[::3] = np.arange(self.nf)
        I[1::3] = np.arange(self.nf)
        I[2::3] = np.arange(self.nf)

        J = self.lfe.reshape((self.nf*3,))
        V = np.sign(J, dtype=np.int32)
        J = abs(J)-1
        return scipy.sparse.csc_matrix((V, (I, J)), dtype=np.int32)

    def get_loop_edge_vec(self,path):
        ''' Convert vertex chains to edges and compute path length
        '''
        edges = np.zeros((self.ne,), dtype=np.int32)
        distance = 0.0
        for i in range(1,len(path)):
            js=min(path[i],path[i-1]) # low points
            je=max(path[i],path[i-1]) # high points
            jp=self.klpe[js]          # pointer into low point list
            jn=self.klpe[js+1]-jp     # number of shared edges
            edge = np.searchsorted(self.llpe[jp:jp+jn],je)
            if self.llpe[jp+edge] != je:
                print(self.llpe[jp:jp+jn],je)
                raise ValueError('Edge not found!')
            distance += np.linalg.norm(self.r[path[i],:]-self.r[path[i-1],:])
            edge += jp
            if path[i]-path[i-1] < 0:
                edges[edge] = -1
            else:
                edges[edge] = 1
        return edges, distance

    def boundary_cycles(self, info=True):
        '''Identify all distinct boundary vertex chains
        '''
        edge_marker = -np.ones((self.ne,))
        # Find cycles
        cycle_ind = -1
        cycle_lists = [[] for _ in range(np.max(self.surf_tag)+1)]
        for ie in range(self.nbe):
            i = self.lbe[ie]
            if edge_marker[i] < 0:
                cycle_ind += 1
                # Start cycle
                edge_marker[i] = cycle_ind
                last_pt = self.le[i,0]
                next_pt = self.le[i,1]
                cycle_pts = [last_pt,next_pt]
                cycle_edges = [i]
                for _ in range(self.nbe):
                    for k in range(self.kpf[next_pt],self.kpf[next_pt+1]):
                        face = self.lpf[k]
                        if self.bf[face]:
                            for l in range(3):
                                k = abs(self.lfe[face,l])-1
                                if (edge_marker[k] >= 0) or (not self.be[k]):
                                    continue
                                if np.any(self.le[k,:] == next_pt) and np.all(self.le[k,:] != last_pt):
                                    edge_marker[k] = cycle_ind
                                    if self.le[k,0] == next_pt:
                                        last_pt = next_pt
                                        next_pt = self.le[k,1]
                                    else:
                                        last_pt = next_pt
                                        next_pt = self.le[k,0]
                                    cycle_pts.append(next_pt)
                                    cycle_edges.append(k)
                                    break
                            else:
                                continue
                            break
                # Detect orientation
                face = self.lef[i,0]
                for l in range(3):
                    if abs(self.lfe[face,l])-1 == i:
                        if self.lfe[face,l] < 0:
                            flip = False
                        else:
                            flip = True
                if flip:
                    cycle_tmp = np.flip(cycle_pts)
                else:
                    cycle_tmp = np.array(cycle_pts)
                # Check to make sure we start at a suitable point by avoiding boundary points that connect via non-boundary edges
                for i in range(cycle_tmp.shape[0]):
                    next_pt = cycle_tmp[i]
                    for k in range(self.kpf[next_pt],self.kpf[next_pt+1]):
                        face = self.lpf[k]
                        for l in range(3):
                            if not (self.bp[self.lf[face,l]] and (self.lf[face,l] in cycle_tmp)):
                                break
                        else:
                            break # Point lies on an interior corner of the cycle (skip)
                    else:
                        break
                else:
                    raise ValueError("Could not find suitable starting point for boundary cycle")
                if i > 0:
                    if info:
                        print("Shifting boundary cycle to {0}".format(i))
                    cycle_tmp = np.hstack((cycle_tmp[i:-1], cycle_tmp[:i+1]))
                #
                cycle_lists[self.surf_tag[face]].append(cycle_tmp)
        k = 0
        for cycle_list in cycle_lists:
            k += len(cycle_list)
        if info:
            print(indent_level+'  Found {0} boundary cycles'.format(k))
        return cycle_lists

    def merge_cells(self,eflag):
        ''' Merge all possible cells while retaining marked edge features
        '''
        def flag_cells(face,cell_group):
            next_faces = []
            for j in range(3):
                if eflag[abs(self.lfe[face,j])-1] > 0:
                    continue
                face2 = self.lff[face,j]
                if face2 < 0:
                    continue
                if cell_group[face2] >= 0:
                    continue
                cell_group[face2] = group_id
                next_faces.append(face2)
            for face2 in next_faces:
                flag_cells(face2,cell_group)
        #
        cell_group = [-1 for _ in range(self.nf)]
        recur_lim = sys.getrecursionlimit()
        sys.setrecursionlimit(self.nf)
        group_id=-1
        for i in range(self.nf):
            if cell_group[i]<0:
                group_id += 1
                cell_group[i] = group_id
                flag_cells(i,cell_group)
        sys.setrecursionlimit(recur_lim)
        cell_group = np.array(cell_group)
        # return cell_group, np.where(cell_group[self.lef[:,0]]!=cell_group[self.lef[:,1]])[0]
        keep_edges = np.where(cell_group[self.lef[:,0]]!=cell_group[self.lef[:,1]])[0]
        pt_flag = np.zeros((self.np,), dtype=np.int32)
        for i in keep_edges:
            pt_flag[self.le[i,0]] += 1
            pt_flag[self.le[i,1]] += 1
        keep_edges = np.where(np.logical_and(cell_group[self.lef[:,0]]!=cell_group[self.lef[:,1]], np.logical_or(pt_flag[self.le[:,0]]>2, pt_flag[self.le[:,1]]>2)))[0]
        return cell_group, keep_edges


def update_svd_with_row(U, S, V_T, new_row):
    """
    Updates SVD (U, S, V_T) when a new row is appended to the matrix.
    Assumes matrix shape is (m, n) where rows are added.
    """
    # Project new row into the current V space (V_T is k x n)
    q = V_T @ new_row.reshape(-1, 1)

    # Compute orthogonal residual (projection error)
    residual = new_row.reshape(-1, 1) - V_T.T @ q
    p = np.linalg.norm(residual)

    if p > 1e-10:
        P = residual / p
    else:
        P = np.zeros_like(residual)
        p = 0.0

    # Construct the small core matrix K (size: (k+1) x (k+1))
    upper = np.hstack([np.diag(S), np.zeros((len(S), 1))])
    lower = np.hstack([q.T, [[p]]])
    K_mat = np.vstack([upper, lower])

    # Compute SVD of the small core matrix
    # u_k, s_k, vt_k = np.linalg.svd(K_mat, full_matrices=False)
    u_k, s_k, vt_k, info = scipy.linalg.lapack.sgesdd(K_mat, full_matrices=0)
    if info != 0:
        u_k, s_k, vt_k, info = scipy.linalg.lapack.sgesvd(K_mat, full_matrices=0)
        if info != 0:
            raise np.linalg.LinAlgError("SVD did not converge")

    # Expand U to accommodate the new row index
    u_ext = np.pad(U, ((0, 1), (0, 1)), mode='constant')
    u_ext[-1, -1] = 1.0
    U_new = u_ext @ u_k

    # Update V_T
    V_T_new = vt_k @ np.vstack([V_T, P.T])

    return U_new, s_k, V_T_new


def compute_greedy_homotopy_basis(face,vertex,bi,face_sweight=None):
    '''Compute the single-point Homotopy basis using the greedy method
    of Erickson and Whittlesey
    '''
    def get_path(istart, pred, path):
        if path[istart] is None:
            # Build path from vertex i to root bi
            v = istart
            path_tmp = [v]
            while pred[v] >= 0:
                v = pred[v]
                path_tmp.append(v)
            # Populate subpaths for all vertices in path_tmp
            for j, v in enumerate(path_tmp):
                if path[v] is not None:
                    break
                path[v] = path_tmp[j:]
        return path[istart]

    nf = face.shape[0]
    nv = vertex.shape[0]

    # Compute adjacency matrices
    I = face.reshape((nf*3,))
    J = face[:,[1,2,0]].reshape((nf*3,))
    V = np.zeros((nf*3,),dtype=np.int32)
    V[::3] = np.arange(nf)
    V[1::3] = np.arange(nf)
    V[2::3] = np.arange(nf)

    amd = scipy.sparse.csc_matrix((V+1, (I, J)),shape=(3*nf,3*nf),dtype=np.int32)
    am = amd.copy()
    am.data.fill(1)
    am += am.transpose()

    # Compute edge matrices
    I,J,_ = scipy.sparse.find(am)
    ind = I<J
    edge = np.vstack((I[ind],J[ind]))
    amd_bool = amd.astype('bool')
    am_bool = am.astype('bool')
    tmp_mat = amd-((amd_bool>am_bool)+(amd_bool<am_bool))
    _, _, V = scipy.sparse.find(tmp_mat)
    _, _, V2 = scipy.sparse.find(tmp_mat.transpose())
    eif = np.vstack((V2[ind],V[ind])) - 1

    # Compute edge length-weighted graph
    I, J, _ = scipy.sparse.find(am)
    el = np.sqrt(np.linalg.norm(vertex[I,:]-vertex[J,:],axis=1))
    G = scipy.sparse.csr_matrix((el, (I, J)),shape=(nv,nv))
    if face_sweight is not None: # Add extra weight to covering edges
        I = face[face_sweight:,0].flatten()
        J = face[face_sweight:,2].flatten()
        el2 = 1.E2*np.ones((I.shape[0],))
        G += scipy.sparse.csr_matrix((el2, (I, J)),shape=(nv,nv))
    G += G.transpose()

    # Compute tree (T) of shortest paths in G using Dijkstra's algorithm
    distance, pred = scipy.sparse.csgraph.dijkstra(G,indices=bi,return_predecessors=True)

    # Compute dual graph G*
    ind = np.logical_and(eif[0,:]>=0,eif[1,:]>=0)
    eif2 = eif[:,ind]
    Gs = scipy.sparse.csc_matrix((np.ones((eif2.shape[1],)),(eif2[0,:],eif2[1,:])),shape=(nf,nf))
    Gs += Gs.transpose()

    # Form graph (G\T)* as edges of G* that are not in T
    I = np.arange(nv)
    I = np.delete(I, bi)
    J = pred[I]
    # (I2,J2) and (J2,I2) are faces corresponding to edge (I,pred[I])
    I2 = amd[I,J]
    J2 = amd[J,I]
    ind = np.logical_not(np.logical_or(I2==0,J2==0))
    I2 = I2[ind]-1
    J2 = J2[ind]-1
    Gs[I2,J2] = 0
    Gs[J2,I2] = 0

    # Build spanning tree (T*) of (G\T)*
    I,J,_ = scipy.sparse.find(Gs)
    ind1 = np.hstack((eif[0,:],eif[1,:]))
    ind2 = np.hstack((eif[1,:],eif[0,:]))
    vals = np.hstack((edge[0,:],edge[1,:]))
    F2E = scipy.sparse.csc_matrix((vals+1,(ind1,ind2)),shape=(nf,nf))
    ei = np.vstack((F2E[I,J]-1,F2E[J,I]-1))
    dvi = vertex[ei[0,:],:]-vertex[ei[1,:],:]
    V = -((distance[ei[0,:]]+distance[ei[1,:]]).flatten()+np.linalg.norm(dvi[0,:,:],axis=1))
    GTs = scipy.sparse.csc_matrix((V, (I, J)),shape=(nf,nf))
    tree = scipy.sparse.csgraph.minimum_spanning_tree(GTs)

    # Modify graph G, to contain only edges neither in T nor crossed by edges in T*
    # Remove edges in T
    I = np.arange(nv)
    I = np.delete(I, bi)
    J = pred[I]
    G[I,J] = 0
    G[J,I] = 0
    # Remove edges crossed by edges in T*
    I,J,_ = scipy.sparse.find(tree)
    ei = np.vstack((F2E[I,J]-1,F2E[J,I]-1))
    G[ei[0,:],ei[1,:]] = 0
    G[ei[1,:],ei[0,:]] = 0

    # The homotopy basis consists of all loops (e), where e is an edge of G
    G = scipy.sparse.tril(G)
    I,J,_ = scipy.sparse.find(G)
    basis_cycles = []
    path = [None for _ in range(nv)]
    for i in range(len(I)):
        pi = get_path(I[i], pred, path)
        pj = get_path(J[i], pred, path)
        basis_cycles.append(np.hstack((np.flip(pi),pj)))
    return basis_cycles

def fixup_loop(cycle,mesh,boundary_cycles,debug):
    ''' Fix up vertex chain by replacing paths that cross the boundary and performing some simplification
    '''
    # Shrink corners
    while True:
        if len(cycle) < 2:
            raise ValueError("Cycle has shrunk to zero size")
        if cycle[-2] == cycle[1]:
            cycle = np.delete(cycle,[0,-1])
            continue
        for k in range(mesh.kpf[cycle[0]],mesh.kpf[cycle[0]+1]):
            face = mesh.lpf[k]
            if (cycle[1] in mesh.lf[face,:]) and (cycle[-2] in mesh.lf[face,:]):
                if cycle[-2] in mesh.lf[face,:]:
                    cycle = np.delete(cycle,[0,-1])
                    cycle = np.append(cycle,cycle[0])
                    continue
        #
        for i in range(1,cycle.shape[0]-1):
            if cycle[i-1] == cycle[i+1]:
                break
            for k in range(mesh.kpf[cycle[i]],mesh.kpf[cycle[i]+1]):
                face = mesh.lpf[k]
                if (cycle[i-1] in mesh.lf[face,:]) and (cycle[i+1] in mesh.lf[face,:]):
                    break
            else:
                continue
            break
        else:
            break
        cycle = np.delete(cycle,i)
    # Stitch boundary cycle crossings
    cycle_starts = [boundary_cycle[0] for boundary_cycle in boundary_cycles]
    while True:
        for i in range(1,cycle.shape[0]-1):
            if cycle[i] in cycle_starts:
                bCycle = np.where(cycle_starts==cycle[i])[0][0]
                for j in (-1,1):
                    if cycle[i+j] in boundary_cycles[bCycle][2:-2]:
                        sCycle = np.where(boundary_cycles[bCycle]==cycle[i+j])[0][0]
                        if debug:
                            print("  Reconnecting pts {0} -> {1}".format(cycle[i],cycle[i+j]))
                        if sCycle > len(boundary_cycles[bCycle])/2:
                            insert_loop = boundary_cycles[bCycle][sCycle:].tolist()
                        else:
                            insert_loop = np.flip(boundary_cycles[bCycle][:sCycle+1]).tolist()
                        if j < 0:
                            cycle = np.array(cycle[:i+j].tolist() + insert_loop + cycle[i+1:].tolist())
                        else:
                            cycle = np.array(cycle[:i].tolist() + insert_loop + cycle[i+j+1:].tolist())
                        break
                else:
                    continue
                break
        else:
            break
    return cycle


def compute_homology(in_file, out_file=None, plot_final=False, plot_steps=False, show_omitted=False,
                     show_covering=False, debug=False, ref_point=None, optimize_holes=False, verify_only=False):
    '''! Compute holes and closures for ThinCurr meshes using a Greedy Homotopy approach

    @param in_file Input mesh file
    @param out_file Ouput mesh file (default: Input file name with "-homology" appended)
    @param plot_final Show final homology basis?
    @param plot_steps Show intermediate bases for each distinct surface?
    @param show_omitted Show boundary cycles that are omitted?
    @param show_covering Show covering triangles for boundary cycles (only used if `plot_steps`)?
    @param debug Print additional debug information?
    @param ref_point Reference location for base point [x,y,z] (default: [0,0,0])
    @param optimize_holes Sample additional points to attempt to optimize holes?
    @param verify_only Verify the correct number of holes and closures are found, but do not write output file?
    '''
    global indent_level

    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + "-homology.h5"

    if ref_point is not None:
        ref_point = np.asarray(ref_point)
        if ref_point.shape[0] != 3:
            raise ValueError('`ref_point` must be a 3-vector')
    else:
        ref_point = np.r_[0.0,0.0,0.0]

    if optimize_holes and verify_only:
        raise ValueError('Optimizing holes not necessary or permitted when `verify_only` is True')
    print_info = debug or (not verify_only)

    # Load mesh
    input_model = read_ThinCurr_mesh(in_file)
    vertex_full = input_model['r'].copy()
    face_full = input_model['lc'].copy()
    if 'thincurr' in input_model:
        if 'nfp' in input_model['thincurr']:
            if verify_only:
                print("Periodic ThinCurr model detected (nfp > 0), skipping homology calculation")
                return
            else:
                raise ValueError('Homology calculation does not support periodic ThinCurr meshes (nfp > 0) at this time')
    else:
        input_model['thincurr'] = {}


    # Setup full mesh
    full_mesh = trimesh(vertex_full,face_full,info=print_info)
    boundary_cycles = full_mesh.boundary_cycles(info=print_info)
    indent_level = '  '

    internal_holes = []
    holes = []
    skipped_holes = []
    closures = []
    for surf_id in range(np.max(full_mesh.surf_tag)+1):
        if print_info:
            print()
            print("Analyzing surface {0} of {1}".format(surf_id+1,np.max(full_mesh.surf_tag)+1))

        # Isolate surface from full mesh
        face_mask = (full_mesh.surf_tag==surf_id)
        reindex_flag = np.zeros((vertex_full.shape[0]+1,), dtype=np.int32)
        reindex_flag[face_full[face_mask,:].flatten()+1] = 1
        vertex = vertex_full[reindex_flag[1:] == 1,:]
        rindexed_pts = np.cumsum(reindex_flag)
        face = rindexed_pts[face_full[face_mask,:]+1]-1

        mesh = trimesh(vertex,face,info=debug)
        reindex_inv = [0 for _ in range(mesh.np)]
        for i in range(full_mesh.np):
            if reindex_flag[i+1] == 1:
                reindex_inv[rindexed_pts[i+1]-1] = i
        reindex_inv = np.array(reindex_inv,dtype=np.int32)

        # Identify boundary cycles and stitch over each to seal mesh
        new_ne = 0
        new_nf = 0
        if len(boundary_cycles[surf_id]) > 0:
            new_bcycles = []
            new_lc = []
            cycle_max = [0, 0]
            for k, cycle in enumerate(boundary_cycles[surf_id]):
                new_bcycles.append(rindexed_pts[cycle+1]-1)
                cycle_lc = []
                for i in range(1,len(cycle)-2):
                    cycle_lc.append([cycle[i+1], cycle[i], cycle[0]])
                new_lc = new_lc + cycle_lc
                new_ne += (len(cycle_lc)-1)
                new_nf += len(cycle_lc)
                if len(cycle_lc) > cycle_max[0]:
                    cycle_max = [len(cycle_lc), k]
            new_lc = rindexed_pts[np.asarray(new_lc)+1]-1
            face_covered = np.vstack((face,new_lc))
        else:
            new_lc = []
            new_bcycles = []
            face_covered = face
            ind = np.argmax(face_mask)
            if print_info:
                print("  No boundary cycles, adding closure element at face {0}".format(ind))
            closures.append(ind)

        # Compute Homotopy basis from basepoint
        if print_info:
            print("  Computing Homology Basis")
            print('    Euler Characteristic (covered) = {0} ({1})'.format(mesh.np-mesh.ne+mesh.nf,mesh.np-(mesh.ne+new_ne)+(mesh.nf+new_nf)))
            print("    # of boundary cycles = {0}".format(len(boundary_cycles[surf_id])))
        ind = np.linalg.norm(vertex-ref_point[np.newaxis,:],axis=1).argmin()
        hb = compute_greedy_homotopy_basis(face_covered,vertex,ind,face_sweight=mesh.nf)
        if print_info:
            print("    # of internal cycles = {0}".format(len(hb)))
        for i in range(len(hb)):
            hb[i] = fixup_loop(hb[i],mesh,new_bcycles,debug)

        for k, cycle in enumerate(boundary_cycles[surf_id]):
            if k == cycle_max[1]:
                skipped_holes.append(cycle)
            else:
                holes.append(cycle)

        hb_out = hb
        # Optimize over cycles from basepoint homotopy to produce better looking basis
        if optimize_holes and (len(hb) > 0):
            print("  Optimizing internal cycles")
            indent_level += '  '
            mesh_covered = trimesh(vertex,face_covered,info=debug)
            bmat_dense_base = mesh_covered.get_face_edge_bop().todense()

            # Compute several additional basis sets
            for j in range(len(hb)):
                minima_sets = [cycle for cycle in hb_out]

                ind2 = hb[j][int(len(hb[j])/2)]
                hb2 = compute_greedy_homotopy_basis(face_covered,vertex,ind2,face_sweight=mesh.nf)
                for i in range(len(hb2)):
                    hb2[i] = fixup_loop(hb2[i],mesh,new_bcycles,debug)
                minima_sets += hb2

                # Build edge operator for comparison
                minima_counts = []
                he = []
                he_mark = np.zeros((mesh_covered.ne,))
                for i in range(len(minima_sets)):
                    evec, distance = mesh_covered.get_loop_edge_vec(minima_sets[i])
                    he.append(evec)
                    minima_counts.append(distance)
                    he_mark[abs(evec)>0] = 1

                # Shrink graph by grouping cells that don't cross cycles
                cell_flags, keep_edges = mesh_covered.merge_cells(he_mark)
                ncoarse = np.max(cell_flags)+1
                print(indent_level + "[{2}/{3}] Reducing mesh to {0} macro cells with {1} macro edges".format(ncoarse,keep_edges.shape[0],j+1,len(hb)))
                bmat_tmp = bmat_dense_base[:,keep_edges]
                bmat_dense = np.zeros((ncoarse,keep_edges.shape[0]), np.float32)
                for i in range(ncoarse):
                    bmat_dense[i,:] = np.sum(bmat_tmp[cell_flags==i,:],axis=0)

                # Build list of cycles from smallest to largest
                # intial_rank = np.linalg.matrix_rank(bmat_dense)
                # U, S, V_T = np.linalg.svd(bmat_dense, full_matrices=False)
                U, S, V_T, _ = scipy.linalg.lapack.sgesdd(bmat_dense, full_matrices=0)
                intial_rank = np.sum(S > 1e-10)
                hb_out = []
                do_check = False
                for i in np.argsort(minima_counts):
                    bmat_tmp = np.vstack((bmat_dense,he[i][keep_edges]))
                    if (not do_check) and (i == len(hb)): # Only start checking once we are looking at new cycles
                        do_check = True
                        # Update SVD with final set of previously accepted cycles
                        U, S, V_T, _ = scipy.linalg.lapack.sgesdd(bmat_dense, full_matrices=0)
                    if do_check:
                        # aug_rank = np.linalg.matrix_rank(bmat_tmp)
                        try:
                            U, S, V_T = update_svd_with_row(U, S, V_T, he[i][keep_edges])
                        except np.linalg.LinAlgError: # Fall back to full factorization
                            # U, S, V_T = np.linalg.svd(bmat_tmp, full_matrices=False)
                            U, S, V_T, _ = scipy.linalg.lapack.sgesdd(bmat_tmp, full_matrices=0)
                        aug_rank = np.sum(S > 1.e-10)
                    else:
                        aug_rank = intial_rank + 1
                    if aug_rank != intial_rank:
                        if debug:
                            print("Adding cycle {0}".format(i))
                        bmat_dense = bmat_tmp
                        hb_out.append(minima_sets[i])
                        intial_rank = aug_rank
                        if len(hb_out) == len(hb):
                            break
                    else:
                        if debug:
                            print("Skipping cycle {0}".format(i))
            indent_level = indent_level[:-2]

        # Save computed internal cycles to hole list
        for basis_cycle in hb_out:
            internal_holes.append(reindex_inv[basis_cycle])

        # Plot intermediate cycles
        if plot_steps:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection='3d')
            ax.plot_trisurf(vertex[:,0], vertex[:,1], vertex[:,2], triangles=face, color=[0.0, 0.0, 0.0, 0.1])
            if show_covering and (len(new_lc) > 0):
                ax.plot_trisurf(vertex[:,0], vertex[:,1], vertex[:,2], triangles=new_lc, color='g')
            for k, cycle in enumerate(new_bcycles):
                if k == cycle_max[1]:
                    if show_omitted:
                        ax.plot(vertex[cycle,0], vertex[cycle,1], vertex[cycle,2], color='k')
                else:
                    ax.plot(vertex[cycle,0], vertex[cycle,1], vertex[cycle,2], color='tab:blue')
            for basis_cycle in hb_out:
                ax.plot(vertex[basis_cycle,0], vertex[basis_cycle,1], vertex[basis_cycle,2], color='tab:orange')
            ax.plot(vertex[ind,0], vertex[ind,1], vertex[ind,2], 'o', color='tab:orange')
            ax.set_aspect('equal','box')
            plt.show()

    # Display final stats
    all_cycles = holes + internal_holes
    if print_info:
        print()
        print("Final model:")
        print("    # of holes = {0}".format(len(all_cycles)))
        print("    # of closures = {0}".format(len(closures)))

    # Verify if requested
    if verify_only:
        # Check holes
        if 'holes' in input_model['thincurr']:
            if len(input_model['thincurr']['holes']) != len(all_cycles):
                raise ValueError("Incorrect number of holes found {1} in file (expected {0})".format(len(all_cycles), len(input_model['thincurr']['holes'])))
        else:
            if 'nodesets' in input_model:
                if len(input_model['nodesets']) >= len(all_cycles):
                    print("""
  Note: Input file appears to be an older ThinCurr mesh file.
    Found {0} nodesets, expected {1} hole(s).""".format(len(input_model['nodesets']), len(all_cycles)))
                    if len(input_model['nodesets']) > len(all_cycles):
                        print("More nodesets than holes may indicate the presence of jumpers instead of an issue.")
            else:
                if len(all_cycles) > 0:
                    raise ValueError("Incorrect number of holes found 0 in file (expected {0})".format(len(all_cycles)))
        # Check closures
        if 'closures' in input_model['thincurr']:
            if len(input_model['thincurr']['closures']) != len(closures):
                raise ValueError("Incorrect number of closures found {1} in file (expected {0})".format(len(closures), len(input_model['thincurr']['closures'])))
        else:
            if 'sidesets' in input_model:
                if len(input_model['sidesets'][0]) == len(closures):
                    print("""
  Note: Input file appears to be an older ThinCurr mesh file.
    Found {0} entries in first sideset, expected {1} closure(s).""".format(len(input_model['sidesets'][0]), len(closures)))
                else:
                    raise ValueError("No closures found and incorrect number of sidesets found {1} in file (expected {0})".format(len(closures), len(input_model['sidesets'])))
            else:
                if len(closures) > 0:
                    raise ValueError("Incorrect number of closures found 0 in file (expected {0})".format(len(closures)))
        return

    # Copy mesh to new file and replace holes/closures
    holes = []
    if len(all_cycles) > 0:
        for k, cycle in enumerate(all_cycles):
            holes.append(cycle[:-1])
    write_ThinCurr_mesh(out_file, input_model['r'], input_model['lc'], input_model['reg'],
                        holes=holes,
                        jumpers=input_model['thincurr'].get('jumpers'), closures=closures,
                        eta_surf=input_model['thincurr'].get('eta_surf'),
                        eta_vol=input_model['thincurr'].get('eta_vol'),
                        thickness=input_model['thincurr'].get('thickness'),
                        coil_sets=input_model['thincurr'].get('coil_sets'),
                        reg_names=input_model.get('reg_names'), reg_attrs=input_model.get('reg_attrs'))

    # Plot final cycles
    if plot_final:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_trisurf(vertex_full[:,0], vertex_full[:,1], vertex_full[:,2], triangles=face_full, color=[0.0, 0.0, 0.0, 0.1])
        for k, cycle in enumerate(holes):
            ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:blue')
        if show_omitted:
            for k, cycle in enumerate(skipped_holes):
                ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='k')
        for k, cycle in enumerate(internal_holes):
            ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:orange')
        for closure_cell in closures:
            closure = full_mesh.lf[closure_cell,0]
            ax.plot(vertex_full[closure,0], vertex_full[closure,1], vertex_full[closure,2], 'o', color='k')
        ax.set_aspect('equal','box')
        plt.show()


def script_entry():
    '''! Command line interface for computing holes and closures for ThinCurr meshes using a Greedy Homotopy approach

    See @ref thincurr_holes_tool
    '''
    parser = argparse.ArgumentParser()
    parser.description = "Compute holes and closures for ThinCurr meshes using a Greedy Homotopy approach"
    parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
    parser.add_argument("--out_file", type=str, default=None, help='Ouput mesh file (default: Input file name with "-homology" appended)')
    parser.add_argument("--keep_nodeset_start", type=int, default=None, help="REMOVED")
    parser.add_argument("--plot_final", action="store_true", default=False, help="Show final homology basis?")
    parser.add_argument("--plot_steps", action="store_true", default=False, help="Show intermediate bases for each distinct surface?")
    parser.add_argument("--show_omitted", action="store_true", default=False, help="Show boundary cycles that are omitted?")
    parser.add_argument("--show_covering", action="store_true", default=False, help="Show covering triangles for boundary cycles (only used if `--plot_steps`)?")
    parser.add_argument("--debug", action="store_true", default=False, help="Print additional debug information?")
    parser.add_argument("--ref_point", default=None, type=float, nargs='+', help='Reference location for base point [x,y,z] (default: [0,0,0])')
    parser.add_argument("--optimize_holes", action="store_true", default=False, help="Sample additional points to attempt to optimize holes?")
    options = parser.parse_args()

    if options.keep_nodeset_start is not None:
        parser.exit(-1, '"--keep_nodeset_start" has been removed, please use `ThinCurr_mesh_tool.py` to convert non-hole nodsets before running this script.')

    compute_homology(
        in_file=options.in_file,
        out_file=options.out_file,
        plot_final=options.plot_final,
        plot_steps=options.plot_steps,
        show_omitted=options.show_omitted,
        show_covering=options.show_covering,
        debug=options.debug,
        ref_point=options.ref_point,
        optimize_holes=options.optimize_holes
    )


## @page thincurr_holes_tool `OFT_ThinCurr_holes`: ThinCurr hole/closure script
#
# @section OFT_ThinCurr_holes_desc Description and options
# This script identifies holes and closures in a ThinCurr surface mesh using
# homology computations.
#
#```shell
# usage: OFT_ThinCurr_holes [-h] --in_file IN_FILE [--out_file OUT_FILE] [--keep_nodeset_start KEEP_NODESET_START] [--plot_final] [--plot_steps] [--show_omitted]
#                                 [--show_covering] [--debug] [--ref_point REF_POINT [REF_POINT ...]] [--optimize_holes]
#
# options:
#     -h, --help            show this help message and exit
#     --in_file IN_FILE     Input mesh file
#     --out_file OUT_FILE   Ouput mesh file (default: Input file name with "-homology" appended)
#     --keep_nodeset_start KEEP_NODESET_START
#                           REMOVED
#     --plot_final          Show final homology basis?
#     --plot_steps          Show intermediate bases for each distinct surface?
#     --show_omitted        Show boundary cycles that are omitted?
#     --show_covering       Show covering triangles for boundary cycles (only used if `--plot_steps`)?
#     --debug               Print additional debug information?
#     --ref_point REF_POINT [REF_POINT ...]
#                           Reference location for base point [x,y,z] (default: [0,0,0])
#     --optimize_holes      Sample additional points to attempt to optimize holes?
#```