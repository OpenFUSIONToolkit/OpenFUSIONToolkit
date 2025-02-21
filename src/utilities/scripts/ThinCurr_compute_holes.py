import os
import sys
import shutil
import argparse
import numpy as np
import scipy
import h5py
import matplotlib.pyplot as plt

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
        for i in range(self.nf):      # loop over cells & index to edges
            for j in range(3):        # loop over edges
                k=abs(self.lfe[i,j])-1  # edge numbers
                if self.lef[k,0]<0:
                    self.lef[k,0]=i   # record first cell
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
        def orient_neighbors(face,oriented):
            next_faces = []
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
                next_faces.append(face2)
            for face2 in next_faces:
                orient_neighbors(face2,oriented)
        #
        oriented = [-1 for _ in range(self.nf)]
        recur_lim = sys.getrecursionlimit()
        sys.setrecursionlimit(self.nf)
        surf_id=-1
        for i in range(self.nf):
            if oriented[i]<0:
                surf_id += 1
                oriented[i] = surf_id
                orient_neighbors(i,oriented)
        sys.setrecursionlimit(recur_lim)
        return np.array(oriented)
    
    def get_face_edge_bop(self):
        ''' Compute face to edge boundary operator \partial_2
        '''
        I = np.zeros((self.nf*3,),dtype=np.int32); I[::3] = np.arange(self.nf); I[1::3] = np.arange(self.nf); I[2::3] = np.arange(self.nf)
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
    
    def boundary_cycles(self):
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
                        if self.bp[self.lf[face,0]] and self.bp[self.lf[face,1]] and self.bp[self.lf[face,2]]:
                            break
                    else:
                        break
                else:
                    raise ValueError("Could not find suitable starting point for boundary cycle")
                if i > 0:
                    print("Shifting boundary cycle to {0}".format(i))
                    cycle_tmp = np.hstack((cycle_tmp[i:-1], cycle_tmp[:i+1]))
                #
                cycle_lists[self.surf_tag[face]].append(cycle_tmp)
        k = 0
        for cycle_list in cycle_lists:
            k += len(cycle_list)
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
        return cell_group, np.where(cell_group[self.lef[:,0]]!=cell_group[self.lef[:,1]])[0]

def compute_greedy_homotopy_basis(face,vertex,bi,face_sweight=None):
    '''Compute the single-point Homotopy basis using the greedy method
    of Erickson and Whittlesey
    '''
    nf = face.shape[0]
    nv = vertex.shape[0]

    # Compute adjacency matrices
    I = face.reshape((nf*3,))
    J = face[:,[1,2,0]].reshape((nf*3,))
    V = np.zeros((nf*3,),dtype=np.int32); V[::3] = np.arange(nf); V[1::3] = np.arange(nf); V[2::3] = np.arange(nf)
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
    distance,pred = scipy.sparse.csgraph.dijkstra(G,indices=bi,return_predecessors=True)
    path = []
    for i in range(nv):
        v = i
        path_tmp = [v]
        while pred[v] >= 0:
            v = pred[v]
            path_tmp.append(v)
        path.append(np.flip(path_tmp))
    
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
    for i in range(len(I)):
        pi = path[I[i]]
        pj = path[J[i]]
        basis_cycles.append(np.hstack((pi,np.flip(pj))))
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


parser = argparse.ArgumentParser()
parser.description = "Compute holes and closures for ThinCurr meshes using a Greedy Homotopy approach"
parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
parser.add_argument("--out_file", type=str, default=None, help='Ouput mesh file (default: Input file name with "-homology" appended)')
parser.add_argument("--keep_nodeset_start", type=int, default=None, help="Starting index of nodesets to keep from input file")
parser.add_argument("--plot_final", action="store_true", default=False, help="Show final homology basis?")
parser.add_argument("--plot_steps", action="store_true", default=False, help="Show intermediate bases for each distinct surface?")
parser.add_argument("--show_omitted", action="store_true", default=False, help="Show boundary cycles that are omitted?")
parser.add_argument("--show_covering", action="store_true", default=False, help="Show covering triangles for boundary cycles (only used if `--plot_steps`)?")
parser.add_argument("--debug", action="store_true", default=False, help="Print additional debug information?")
parser.add_argument("--ref_point", default=None, type=float, nargs='+', help='Reference location for base point [x,y,z] (default: [0,0,0])')
parser.add_argument("--optimize_holes", action="store_true", default=False, help="Sample additional points to attempt to optimize holes?")
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
    keep_nodesets = []
    if options.keep_nodeset_start is not None:
        if 'mesh/NUM_NODESETS' not in fid:
            parser.exit(-1, '"--keep_nodeset_start" specified but no nodesets available')
        num_nodesets = fid['mesh/NUM_NODESETS'][0]
        if options.keep_nodeset_start < 0:
            if options.keep_nodeset_start < -num_nodesets:
                parser.exit(-1, '"--keep_nodeset_start" exceeds the number of available nodesets')
            options.keep_nodeset_start = num_nodesets + 1 + options.keep_nodeset_start
        if options.keep_nodeset_start > num_nodesets:
            parser.exit(-1, '"--keep_nodeset_start" exceeds the number of available nodesets')
        for j in range(num_nodesets):
            if j+1 >= options.keep_nodeset_start:
                try:
                    keep_nodesets.append(np.asarray(fid['mesh/NODESET{0:04d}'.format(j+1)])-1)
                except:
                    parser.exit(-1, 'Failed to read nodeset {0}'.format(j+1))

# Setup full mesh
full_mesh = trimesh(vertex_full,face_full)
boundary_cycles = full_mesh.boundary_cycles()
indent_level = '  '

internal_holes = []
holes = []
skipped_holes = []
closures = []
for surf_id in range(np.max(full_mesh.surf_tag)+1):
    print()
    print("Analyzing surface {0} of {1}".format(surf_id+1,np.max(full_mesh.surf_tag)+1))

    # Isolate surface from full mesh
    face_mask = (full_mesh.surf_tag==surf_id)
    reindex_flag = np.zeros((vertex_full.shape[0]+1,), dtype=np.int32)
    reindex_flag[face_full[face_mask,:].flatten()+1] = 1
    vertex = vertex_full[reindex_flag[1:] == 1,:]
    rindexed_pts = np.cumsum(reindex_flag)
    face = rindexed_pts[face_full[face_mask,:]+1]-1

    mesh = trimesh(vertex,face,info=options.debug)
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
        print("  No boundary cycles, adding closure element at face {0}".format(ind))
        closures.append(ind)

    # Compute Homotopy basis from basepoint
    print("  Computing Homology Basis")
    print('    Euler Characteristic (covered) = {0} ({1})'.format(mesh.np-mesh.ne+mesh.nf,mesh.np-(mesh.ne+new_ne)+(mesh.nf+new_nf)))
    print("    # of boundary cycles = {0}".format(len(boundary_cycles[surf_id])))
    ind = np.linalg.norm(vertex-ref_point[np.newaxis,:],axis=1).argmin()
    hb = compute_greedy_homotopy_basis(face_covered,vertex,ind,face_sweight=mesh.nf)
    print("    # of internal cycles = {0}".format(len(hb)))
    for i in range(len(hb)):
        hb[i] = fixup_loop(hb[i],mesh,new_bcycles,options.debug)
    
    for k, cycle in enumerate(boundary_cycles[surf_id]):
        if k == cycle_max[1]:
            skipped_holes.append(cycle)
        else:
            holes.append(cycle)

    hb_out = hb
    # Optimize over cycles from basepoint homotopy to produce better looking basis
    if options.optimize_holes and (len(hb) > 0):
        print("  Optimizing internal cycles")
        indent_level += '  '
        mesh_covered = trimesh(vertex,face_covered,info=options.debug)
        bmat_dense_base = mesh_covered.get_face_edge_bop().todense()

        # Compute several additional basis sets
        for j in range(len(hb)):
            minima_sets = [cycle for cycle in hb_out]

            ind2 = hb[j][int(len(hb[j])/2)]
            hb2 = compute_greedy_homotopy_basis(face_covered,vertex,ind2,face_sweight=mesh.nf)
            for i in range(len(hb2)):
                hb2[i] = fixup_loop(hb2[i],mesh,new_bcycles,options.debug)        
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
            print(indent_level + "[{2}/{3}] Reducing mesh to {0} macro cells with {1} edges".format(ncoarse,keep_edges.shape[0],j+1,len(hb)))
            bmat_tmp = bmat_dense_base[:,keep_edges]
            bmat_dense = np.zeros((ncoarse,keep_edges.shape[0]))
            for i in range(ncoarse):
                bmat_dense[i,:] = np.sum(bmat_tmp[cell_flags==i,:],axis=0)
            
            # Build list of cycles from smallest to largest
            intial_rank = np.linalg.matrix_rank(bmat_dense)
            hb_out = []
            for i in np.argsort(minima_counts):
                bmat_tmp = np.vstack((bmat_dense,he[i][keep_edges]))
                aug_rank = np.linalg.matrix_rank(bmat_tmp)
                if aug_rank != intial_rank:
                    if options.debug:
                        print("Adding cycle {0}".format(i))
                    bmat_dense = bmat_tmp
                    hb_out.append(minima_sets[i])
                    intial_rank = aug_rank
                    if len(hb_out) == len(hb):
                        break
                else:
                    if options.debug:
                        print("Skipping cycle {0}".format(i))
        indent_level = indent_level[:-2]
    
    # Save computed internal cycles to hole list
    for basis_cycle in hb_out:
        internal_holes.append(reindex_inv[basis_cycle])

    # Plot intermediate cycles
    if options.plot_steps:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_trisurf(vertex[:,0], vertex[:,1], vertex[:,2], triangles=face, color=[0.0, 0.0, 0.0, 0.1])
        if options.show_covering and (len(new_lc) > 0):
            ax.plot_trisurf(vertex[:,0], vertex[:,1], vertex[:,2], triangles=new_lc, color='g')
        for k, cycle in enumerate(new_bcycles):
            if k == cycle_max[1]:
                if options.show_omitted:
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
print()
print("Final model:")
print("    # of holes = {0}".format(len(all_cycles)))
print("    # of closures = {0}".format(len(closures)))
print("    # of additional nodesets = {0}".format(len(keep_nodesets)))

# Copy mesh to new file and replace holes/closures
shutil.copyfile(options.in_file,out_file)
with h5py.File(out_file,'r+') as h5_file:
    # Replace nodesets
    if 'mesh/NUM_NODESETS' in h5_file:
        for j in range(h5_file['mesh/NUM_NODESETS'][0]):
            del h5_file['mesh/NODESET{0:04d}'.format(j+1)]
        del h5_file['mesh/NUM_NODESETS']
    nodesets = []
    if len(all_cycles) > 0:
        for k, cycle in enumerate(all_cycles):
            nodesets.append(cycle[:-1])
    nodesets = nodesets + keep_nodesets
    if len(nodesets) > 0:
        h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(nodesets),], dtype='i4')
        j=0
        for k, nodeset in enumerate(nodesets):
            j+=1
            h5_file.create_dataset('mesh/NODESET{0:04d}'.format(j), data=nodeset+1, dtype='i4')
    # Replace sidesets
    if 'mesh/NUM_SIDESETS' in h5_file:
        for j in range(h5_file['mesh/NUM_SIDESETS'][0]):
            del h5_file['mesh/SIDESET{0:04d}'.format(j+1)]
        del h5_file['mesh/NUM_SIDESETS']
    if len(closures) > 0:
        closures = [closure+1 for closure in closures]
        h5_file.create_dataset('mesh/NUM_SIDESETS', data=[1,], dtype='i4')
        h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(1), data=closures, dtype='i4')

# Plot final cycles
if options.plot_final:
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_trisurf(vertex_full[:,0], vertex_full[:,1], vertex_full[:,2], triangles=face_full, color=[0.0, 0.0, 0.0, 0.1])
    for k, cycle in enumerate(holes):
        ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:blue')
    if options.show_omitted:
        for k, cycle in enumerate(skipped_holes):
            ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='k')
    for k, cycle in enumerate(internal_holes):
        ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:orange')
    for k, cycle in enumerate(keep_nodesets):
        ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:green')
    for closure_cell in closures:
        closure = full_mesh.lf[closure_cell,0]
        ax.plot(vertex_full[closure,0], vertex_full[closure,1], vertex_full[closure,2], 'o', color='k')
    ax.set_aspect('equal','box')
    plt.show()