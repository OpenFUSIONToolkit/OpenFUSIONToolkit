import os
import sys
import argparse
import numpy as np
import scipy
import h5py
import matplotlib.pyplot as plt

tri_ed = np.asarray([[2,1], [0,2], [1,0]]) # Triangle edge list

class trimesh:

    def __init__(self, r, lf):
        if lf.shape[1] != 3:
            raise ValueError("Face list must have shape [nf, 3]")
        self.r=r
        self.lf=lf
        self.np=r.shape[0]
        self.nf=lf.shape[0]
        self.ne = 0
        # Setup mesh
        self.setup_edges()
        self.setup_neighbors()
        self.setup_boundary()
        self.surf_tag = self.sync_normals()
        self.nsurfs = np.max(self.surf_tag)+1
        print("Mesh constructed:")
        print("  Found {0} distinct surface(s)".format(self.nsurfs))
        print("  # of vertices = {0} ({1})".format(self.np, self.nbp))
        print("  # of edges    = {0} ({1})".format(self.ne, self.nbe))
        print("  # of faces    = {0} ({1})".format(self.nf, self.nbf))
    
    def setup_edges(self):
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
        if ir[0]!=0 :
            raise ValueError('Bad raw edge count')
        nr=np.zeros((self.np+1,), dtype=np.int32) # reset raw edge counter
        for i in range(self.nf): # loop over cells & index raw edges
            k=self.lf[i,:] # cell corners
            for j in range(3): # loop over edges
                js=min(k[tri_ed[j,0]],k[tri_ed[j,1]]) # low points
                je=max(k[tri_ed[j,0]],k[tri_ed[j,1]]) # high points
                jr[ir[js]+nr[js]]=je                  # high point list
                nr[js]+=1                             # edge counter
        # !$omp parallel do
        for i in range(self.np): # sort raw high point list
            if(nr[i]>1):
                jr[ir[i]:ir[i]+nr[i]].sort()
                # jr[ir[i]:ir[i]+nr[i]]=np.sort(jr[ir[i]:ir[i]+nr[i]])
        # deallocate(nr)
        self.klpe=np.zeros((self.np+1,), dtype=np.int32)
        # !$omp parallel do private(j,js,je)
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
        # !$omp parallel do private(j,js,je,jn)
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
        # deallocate(ir,jr)
        self.le=np.zeros((self.ne,2), dtype=np.int32)
        # !$omp parallel do
        for i in range(self.np): # construct global edge list
            self.le[self.klpe[i]:self.klpe[i+1],0]=i
        # !$omp parallel do
        for i in range(self.ne):
            self.le[i,1]=self.llpe[i]
        self.lfe=np.zeros((self.nf,3), dtype=np.int32)
        # !$omp parallel do private(k,j,js,je,jp,jn)
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
                self.lfe[i,j]+=jp
                if k[tri_ed[j,1]]-k[tri_ed[j,0]]<0:
                    self.lfe[i,j]*=-1 # apply orientation
    
    def setup_neighbors(self):
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
        # deallocate(nr)
        self.lef=-np.ones((self.ne,2), dtype=np.int32)
        self.lff=np.zeros((self.nf,3), dtype=np.int32)
        for i in range(self.nf):      # loop over cells & index to edges
            for j in range(3):        # loop over edges
                k=abs(self.lfe[i,j])  # edge numbers
                if self.lef[k,0]<0:
                    self.lef[k,0]=i   # record first cell
                else:
                    self.lef[k,1]=i   # record second cell
        #!$omp parallel do private(j,k)
        for i in range(self.nf):     # loop over cells & locate neighbors
            for j in range(3):       # loop over edges
                k=abs(self.lfe[i,j]) # edge numbers
                self.lff[i,j]=np.sum(self.lef[k,:])-i
    
    def setup_boundary(self):
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
    
    def boundary_cycles(self):
        edge_marker = -np.ones((self.ne,))
        #
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
                                k = abs(self.lfe[face,l])
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
                itmp = i
                if i == 0:
                    itmp = cycle_edges[1]
                face = self.lef[itmp,0]
                for l in range(3):
                    if abs(self.lfe[face,l]) == itmp:
                        if self.lfe[face,l] < 0:
                            flip = False
                        else:
                            flip = True
                if flip:
                    cycle_lists[self.surf_tag[face]].append(np.flip(cycle_pts))
                else:
                    cycle_lists[self.surf_tag[face]].append(np.array(cycle_pts))
        k = 0
        for cycle_list in cycle_lists:
            k += len(cycle_list)
        print('  Found {0} boundary cycles'.format(k))
        return cycle_lists
    
    def sync_normals(self):
        def orient_neighbors(face,oriented):
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
                ed2 = self.lf[face2,tri_ed[k,:]]
                if (ed[0]==ed2[0]) and (ed[1]==ed2[1]):
                    self.lf[face2,[1,2]] = self.lf[face2,[2,1]]
                    self.lfe[face2,:] = -self.lfe[face2,[0,2,1]]
                    self.lff[face2,:] = self.lff[face2,[0,2,1]]
                oriented[face2] = surf_id
                orient_neighbors(face2,oriented)
        #
        oriented = [-1 for _ in range(self.nf)]
        recur_lim = sys.getrecursionlimit()
        sys.setrecursionlimit(self.nf)
        surf_id=-1
        for i in range(self.nf):
            if oriented[i]<0:
                surf_id+=1
                oriented[i] = surf_id
                orient_neighbors(i,oriented)
        sys.setrecursionlimit(recur_lim)
        return np.array(oriented)

def compute_greedy_homotopy_basis(face,vertex,bi,face_sweight=None):
    '''Compute the single-point Homotopy basis using the greedy method
    of Erickson and Whittlesey ()
    '''
    nf = face.shape[0]
    nv = vertex.shape[0]

    # Compute adjacency matrices
    I = face.reshape((nf*3,))
    J = face[:,[1,2,0]].reshape((nf*3,))
    V = np.zeros((nf*3,),dtype=np.int32); V[::3] = np.arange(nf); V[1::3] = np.arange(nf); V[2::3] = np.arange(nf)
    amd = scipy.sparse.csc_matrix((V+1, (I, J)),dtype=np.int32)
    am = amd.copy()
    am.data.fill(1)
    am = am+am.transpose()

    # Compute edge matrices
    I,J,_ = scipy.sparse.find(am)
    ind = I<J
    edge = (np.vstack((I[ind],J[ind]))).transpose()
    amd_bool = amd.astype('bool')
    am_bool = am.astype('bool')
    tmp_mat = amd-((amd_bool>am_bool)+(amd_bool<am_bool))
    _, _, V = scipy.sparse.find(tmp_mat)
    _, _, V2 = scipy.sparse.find(tmp_mat.transpose())
    eif = (np.vstack((V2[ind],V[ind]))).transpose() - 1
    
    # Compute edge length-weighted graph
    I, J, _ = scipy.sparse.find(am)
    el = np.sqrt(np.linalg.norm(vertex[I,:]-vertex[J,:],axis=1))
    G = scipy.sparse.csr_matrix((el, (I, J)),shape=(nv,nv))
    if face_sweight is not None: # Add extra weight to covering edges
        I = face[face_sweight:,0].flatten()
        J = face[face_sweight:,2].flatten()
        el2 = 1.E2*np.ones((I.shape[0],))
        G += scipy.sparse.csr_matrix((el2, (I, J)),shape=(nv,nv))
    G = G + G.transpose()
    
    # Compute tree (T) of shortest paths in G using the Dijkstra's algorithm
    distance,pred = scipy.sparse.csgraph.dijkstra(G,indices=bi,return_predecessors=True)
    path = []
    for i in range(nv):
        v = i
        path_tmp = [v]
        while pred[v] >= 0:
            v = pred[v]
            path_tmp.append(v)
        path.append(np.flip(path_tmp))
    
    # Compute dual graph
    ind = np.logical_and(eif[:,0]>=0,eif[:,1]>=0)
    eif2 = eif[ind,:]
    amf = scipy.sparse.csc_matrix((np.ones((eif2.shape[0],)),(eif2[:,0],eif2[:,1])),shape=(nf,nf))
    amf = amf + amf.transpose()

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
    amf[I2,J2] = 0
    amf[J2,I2] = 0
    
    # Build spanning tree (T*) of (G\T)*
    I,J,_ = scipy.sparse.find(amf)
    ind = np.logical_not(np.logical_or(eif[:,0]==-2,eif[:,1]==-2))
    eif = eif[ind,:]
    edge = edge[ind,:]
    ind1 = np.hstack((eif[:,0],eif[:,1]))
    ind2 = np.hstack((eif[:,1],eif[:,0]))
    vals = np.hstack((edge[:,0],edge[:,1]))
    F2E = scipy.sparse.csc_matrix((vals+1,(ind1,ind2)),shape=(nf,nf))
    ei = (np.vstack((F2E[I,J]-1,F2E[J,I]-1))).transpose()
    dvi = vertex[ei[:,0],:]-vertex[ei[:,1],:]
    V = -((distance[ei[:,0]]+distance[ei[:,1]]).flatten()+np.linalg.norm(dvi[:,0,:],axis=1))
    amf_w = scipy.sparse.csc_matrix((V, (I, J)),shape=(nf,nf))
    tree = scipy.sparse.csgraph.minimum_spanning_tree(amf_w).transpose()
    
    # G2 is the graph, with edges neither in T nor are crossed by edges in T*
    G2 = G.copy()
    # Remove edges in T
    I = np.arange(nv)
    I = np.delete(I, bi)
    J = pred[I]
    G2[I,J] = 0
    G2[J,I] = 0
    
    # Remove edges crossed by edges in T*
    I,J,_ = scipy.sparse.find(tree)
    ei = (np.vstack((F2E[I,J]-1,F2E[J,I]-1))).transpose()
    G2[ei[:,0],ei[:,1]] = 0
    G2[ei[:,1],ei[:,0]] = 0
    
    # The greedy homotopy basis consists of all loops (e), where e is an edge of G2.
    G2 = scipy.sparse.tril(G2)
    I,J,_ = scipy.sparse.find(G2)
    basis_cycles = []
    for i in range(len(I)):
        pi = path[I[i]]
        pj = path[J[i]]
        basis_cycles.append(np.hstack((pi,np.flip(pj))))
    return basis_cycles

def fixup_loop(cycle,mesh,boundary_cycles):
    # # Shrink corners
    # while True:
    #     # print(len(cycle))
    #     if cycle[-2] == cycle[1]:
    #         cycle = np.delete(cycle,[0,-1])
    #         continue
    #     for k in range(mesh.kpf[cycle[0]],mesh.kpf[cycle[0]+1]):
    #         face = mesh.lpf[k]
    #         if np.sum(mesh.lf[face,:]) == (cycle[-2]+cycle[0]+cycle[1]):
    #             if cycle[-2] in mesh.lf[face,:]:
    #                 cycle = np.delete(cycle,[0,-1])
    #                 cycle = np.append(cycle,cycle[0])
    #                 continue
    #     #
    #     for i in range(1,cycle.shape[0]-1):
    #         if cycle[i-1] == cycle[i+1]:
    #             break
    #         for k in range(mesh.kpf[cycle[i]],mesh.kpf[cycle[i]+1]):
    #             face = mesh.lpf[k]
    #             if (cycle[i-1] in mesh.lf[face,:]) and (cycle[i] in mesh.lf[face,:]) and (cycle[i+1] in mesh.lf[face,:]):
    #                 break
    #         else:
    #             continue
    #         break
    #     else:
    #         break
    #     cycle = np.delete(cycle,i)
    # Stitch boundary cycle crossings
    cycle_starts = [boundary_cycle[0] for boundary_cycle in boundary_cycles]
    while True:
        for i in range(1,cycle.shape[0]-1):
            if cycle[i] in cycle_starts:
                bCycle = np.where(cycle_starts==cycle[i])[0][0]
                for j in (-1,1):
                    if cycle[i+j] in boundary_cycles[bCycle][2:-2]:
                        sCycle = np.where(boundary_cycles[bCycle]==cycle[i+j])[0][0]
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


parser = argparse.ArgumentParser()
parser.description = "Compute holes for ThinCurr meshes"
parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
parser.add_argument("--out_file", type=str, default=None, help="Ouput mesh file")
parser.add_argument("--plot", action="store_true", default=False, help="Show final homology basis")
parser.add_argument("--plot_all", action="store_true", default=False, help="Show intermediate bases for each distinct surface")
parser.add_argument("--ref_point", default=None, type=float, nargs='+', help='Reference location for base point')
options = parser.parse_args()

out_file = options.out_file
if out_file is None:
    out_file = os.path.splitext(options.in_file)[0] + "-homology.h5"

if options.ref_point is not None:
    ref_point = np.asarray(options.ref_point)
    if ref_point.shape[0] != 3:
        parser.exit(-1, '"--ref_point" must have 3 values')
else:
    ref_point = np.r_[0.0,0.0,0.0]

# Load mesh
with h5py.File(options.in_file) as fid:
    vertex_full = np.asarray(fid['mesh/R'])
    face_full = np.asarray(fid['mesh/LC'])-1
    reg_full = np.asarray(fid['mesh/REG'])

full_mesh = trimesh(vertex_full,face_full)
boundary_cycles = full_mesh.boundary_cycles()

internal_holes = []
holes = []
skipped_holes = []
closures = []
for surf_id in range(np.max(full_mesh.surf_tag)+1):
    #
    print()
    print("Analyzing surface {0} of {1}".format(surf_id+1,np.max(full_mesh.surf_tag)+1))
    #
    face_mask = (full_mesh.surf_tag==surf_id)
    reindex_flag = np.zeros((vertex_full.shape[0]+1,), dtype=np.int32)
    reindex_flag[face_full[face_mask,:].flatten()+1] = 1
    vertex = vertex_full[reindex_flag[1:] == 1,:]
    rindexed_pts = np.cumsum(reindex_flag)
    face = rindexed_pts[face_full[face_mask,:]+1]-1
    #
    mesh = trimesh(vertex,face)
    reindex_inv = [0 for _ in range(mesh.np)]
    for i in range(full_mesh.np):
        if reindex_flag[i+1] == 1:
            reindex_inv[rindexed_pts[i+1]-1] = i
    reindex_inv = np.array(reindex_inv,dtype=np.int32)
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
        new_bcycles = []
        face_covered = face
        ind = np.argmax(face_mask)
        print("  No boundary cycles, adding closure element at face {0}".format(ind))
        closures.append(ind)
    
    print("  Computing Homology Basis")
    print("    # of boundary cycles = {0}".format(len(boundary_cycles[surf_id])))
    print('    Euler Characteristic = {0} ({1})'.format(mesh.np-mesh.ne+mesh.nf,mesh.np-(mesh.ne+new_ne)+(mesh.nf+new_nf)))


    # Compute Homotopy basis
    ind = np.linalg.norm(vertex-ref_point[np.newaxis,:],axis=1).argmin()
    hb = compute_greedy_homotopy_basis(face_covered,vertex,ind,face_sweight=mesh.nf)
    print("    # of internal cycles = {0}".format(len(hb)))
    for i in range(len(hb)):
        hb[i] = fixup_loop(hb[i],mesh,new_bcycles)
    
    for k, cycle in enumerate(boundary_cycles[surf_id]):
        if k == cycle_max[1]:
            skipped_holes.append(cycle)
        else:
            holes.append(cycle)
    for basis_cycle in hb:
        internal_holes.append(reindex_inv[basis_cycle])

    # Plot intermediate cycles
    if options.plot_all:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_trisurf(vertex[:,0], vertex[:,1], vertex[:,2], triangles=face, edgecolor=[0.0, 0.0, 0.0, 0.1] , color='none')
        for k, cycle in enumerate(new_bcycles):
            if k == cycle_max[1]:
                ax.plot(vertex[cycle,0], vertex[cycle,1], vertex[cycle,2], color='tab:green')
            else:
                ax.plot(vertex[cycle,0], vertex[cycle,1], vertex[cycle,2], color='tab:red')
        for basis_cycle in hb:
            ax.plot(vertex[basis_cycle,0], vertex[basis_cycle,1], vertex[basis_cycle,2], color='tab:blue')
        ax.plot(vertex[ind,0], vertex[ind,1], vertex[ind,2],'ro')
        ax.set_aspect('equal','box')
        plt.show()

# Plot final cycles
if options.plot:
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_trisurf(vertex_full[:,0], vertex_full[:,1], vertex_full[:,2], triangles=face_full, color=[0.0, 0.0, 0.0, 0.1])
    for k, cycle in enumerate(holes):
        ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:red')
    for k, cycle in enumerate(skipped_holes):
        ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:green')
    for k, cycle in enumerate(internal_holes):
        ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:blue')
    for closure_cell in closures:
        closure = full_mesh.lf[closure_cell,0]
        ax.plot(vertex_full[closure,0], vertex_full[closure,1], vertex_full[closure,2], 'o', color='tab:orange')
    ax.set_aspect('equal','box')

# Save to new mesh file
all_cycles = holes + internal_holes
print()
with h5py.File(out_file,'w') as h5_file:
    h5_file.create_dataset('mesh/R', data=vertex_full, dtype='f8')
    h5_file.create_dataset('mesh/LC', data=face_full+1, dtype='i4')
    h5_file.create_dataset('mesh/REG', data=reg_full, dtype='i4')
    if len(all_cycles) > 0:
        h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(all_cycles),], dtype='i4')
        j=0
        for k, cycle in enumerate(all_cycles):
            j+=1
            h5_file.create_dataset('mesh/NODESET{0:04d}'.format(j), data=cycle[:-1]+1, dtype='i4')
    if len(closures) > 0:
        h5_file.create_dataset('mesh/NUM_SIDESETS', data=[len(closures),], dtype='i4')
        h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(1), data=closures, dtype='i4')

plt.show()