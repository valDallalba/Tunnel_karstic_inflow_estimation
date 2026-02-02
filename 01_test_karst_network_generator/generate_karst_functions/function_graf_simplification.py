#Valentin Dall'alba & Alexis Neven

import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from geone import grf as grf_geone
import geone.covModel as gcm
from scipy import interpolate
import math


#################
# GRF class
#################

class GRF_objt:
    '''
    Class GRF to create Gaussian random field simulation.
    '''

    def __init__(self, sillVal=1, rangeVal=[150, 50], alpha=0, typeModel='gaussian'):
        self.model_vario = gcm.CovModel2D(elem=[(typeModel, {'w': sillVal, 'r': rangeVal})],
                                          alpha=alpha, name='model_2D_karst')
        return

    def set_dimension(self, dimension=[2000, 200], spacing=[1, 1], origin=[0, 0]):
        self.dimension = dimension
        self.spacing = spacing
        self.origin = origin
        return

    def run_model_GRF(self, nb_real=1):
        self.nb_real = nb_real
        self.simu = grf_geone.grf2D(self.model_vario, self.dimension, self.spacing,
                                    self.origin, nreal=self.nb_real, printInfo=False)
        return

    def grf_to_cdf(self):
        self.cdf = []
        quantiles = np.arange(0, 1 + 0.01, 0.01)

        for simu in self.simu:
            quantiles_val = np.quantile(simu.flatten(), quantiles)
            cdf_fct = interpolate.interp1d(quantiles_val, quantiles)
            proba_from_cdf = cdf_fct(simu)
            self.cdf.append(proba_from_cdf)
        return

    def cdf_to_diameter(self, path_data='./Size_conduits.csv'):
        self.diameter = []
        data = np.genfromtxt(path_data, delimiter=';')
        quantiles = data[:, 0]
        quantiles_val = np.sqrt(data[:, 1] / math.pi) * 2
        diameter_fct = interpolate.interp1d(quantiles, quantiles_val)

        for cdf in self.cdf:
            diameter = diameter_fct(cdf)
            self.diameter.append(diameter)
        return

##################
# Fonctions divers
##################

def get_position_center(node_str, node_end):
    new_x = np.array((node_str[:, 0] + node_end[:, 0]) / 2)
    new_y = np.array((node_str[:, 1] + node_end[:, 1]) / 2)
    return np.array([new_x, new_y])

def get_sgs_values(sgs, position_center):
    value = []
    for pos in position_center:
        value.append(sgs[np.int(pos[1]), np.int(pos[0])])
    return np.array(value)

def simplify_graph(infile_nodes, infile_edges, bbox, dx, dy, verbose=True, seed=12345678):
    """
    Simplifies a graph from pykasso given in two ASCII files.
    Returns a list of nodes, the new segments with diameters, and the pipe diameter field.
    """
    # Load node and edge data
    tmp = np.loadtxt(infile_nodes)
    seg = np.loadtxt(infile_edges)
    x, y = tmp[:, 1], tmp[:, 2]

    # Step 1: Round node coordinates and assign new node ids
    nodes, nodes_corresp, new_nodes = _round_and_map_nodes(x, y, dx, dy)

    # Check for duplicate nodes after rounding
    _warn_duplicate_nodes(new_nodes, step=1, seed=seed)

    if verbose:
        _print_simplification_summary(x, nodes, step=1)

    # Step 2: Simplify edges and handle diagonals
    new_segs, nodes = _simplify_edges(seg, nodes_corresp, new_nodes, nodes, x, y)

    # Remove duplicate segments
    new_segs = _remove_duplicate_segments(new_segs)

    # Rebuild new_nodes array in case new nodes were added
    new_nodes = _rebuild_new_nodes(nodes)

    # Check for duplicate nodes again
    _warn_duplicate_nodes(new_nodes, step=2, seed=seed)

    if verbose:
        _print_simplification_summary(x, nodes, step=2, seg=seg, new_segs=new_segs, new_nodes=new_nodes, bbox=bbox, dx=dx, dy=dy, y=y)

    # Step 3: Ensure all nodes are connected
    new_segs, seg_nodes = _connect_unlinked_nodes(new_nodes, new_segs, seed, verbose)

    # Step 4: Simulate pipe diameters using a Gaussian random field
    newsegs_with_diam, pipe_diameter = _add_pipe_diameters(new_nodes, new_segs, bbox, seed)

    if verbose:
        print('Pipe diameter is added to newsegs array based on sampled sgs simulation!')

    return new_nodes, newsegs_with_diam, pipe_diameter

def _round_and_map_nodes(x, y, dx, dy):
    """Round node coordinates to the grid and assign new node ids."""
    xa = np.floor(x / dx) * dx + dx / 2
    ya = np.floor(y / dy) * dy + dy / 2
    nodes = {}         # key: (x, y) tuple, value: new node id
    nodes_corresp = {} # key: old node id, value: new node id
    j = 0
    for i in range(len(xa)):
        t = (int(xa[i]), int(ya[i]))
        if t in nodes:
            nodes_corresp[i] = nodes[t]
        else:
            nodes[t] = j
            nodes_corresp[i] = j
            j += 1
    new_nodes = np.zeros((len(nodes), 2))
    for t, i in nodes.items():
        new_nodes[i, 0] = t[0]
        new_nodes[i, 1] = t[1]
    return nodes, nodes_corresp, new_nodes

def _warn_duplicate_nodes(new_nodes, step, seed):
    """Warn if duplicate nodes are found."""
    tmp_labels = [f"{new_nodes[i,0]}-{new_nodes[i,1]}" for i in range(len(new_nodes))]
    nnc = sorted(set(tmp_labels))
    if len(nnc) != len(new_nodes):
        if step == 1:
            print("Warning: Several identical nodes after step 1")
        else:
            print("Warning: Several identical nodes ! Should not happen - stop and debug !!!")
            print(f"The corresponding seed is {seed}!")
            print("==============================")

def _print_simplification_summary(x, nodes, step, seg=None, new_segs=None, new_nodes=None, bbox=None, dx=None, dy=None, y=None):
    """Print summary of the simplification process."""
    print("Summary of graph simplification:")
    print("==============================")
    print("Nb nodes export pykasso:", len(x))
    print(f"Nb nodes after simplification - step {step}:", len(nodes))
    print("==============================")
    if step == 2 and seg is not None and new_segs is not None and new_nodes is not None:
        print("Nb edges export pykasso:", len(seg))
        print("Nb edges after simplification:", len(new_segs))
        print("==============================")
        xa, ya = new_nodes[:, 0], new_nodes[:, 1]
        print("Check nodes coordinates and bounding box")
        print("Domain size:", bbox[0], bbox[1], bbox[2], bbox[3])
        print("Original data range:", np.min(x), np.max(x), np.min(y), np.max(y))
        print("Simplified coordinates range:", np.min(xa), np.max(xa), np.min(ya), np.max(ya))
        print("Theoretical box:", bbox[0] + dx / 2, bbox[1] - dx / 2, bbox[2] + dy / 2, bbox[3] - dy / 2)
        print("==============================")

def _simplify_edges(seg, nodes_corresp, new_nodes, nodes, x, y):
    """Simplify the links (edges) and handle diagonals."""
    new_segs = []
    for i in range(len(seg)):
        j0 = nodes_corresp[seg[i, 0]]
        j1 = nodes_corresp[seg[i, 1]]
        if j0 != j1:
            j0, j1 = min(j0, j1), max(j0, j1)
            if [j0, j1] not in new_segs:
                x0, x1 = new_nodes[j0, 0], new_nodes[j1, 0]
                y0, y1 = new_nodes[j0, 1], new_nodes[j1, 1]
                if (x0 == x1) or (y0 == y1):
                    new_segs.append([j0, j1])
                else:
                    # Diagonal connection: add intermediate node if needed
                    i0, i1 = int(seg[i, 0]), int(seg[i, 1])
                    x0o, x1o = x[i0], x[i1]
                    y0o, y1o = y[i0], y[i1]
                    a = (y0o - y1o) / (x0o - x1o)
                    b = y0o - x0o * a
                    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
                    ymp = a * xm + b
                    ynovel = min(y0, y1) if ymp <= ym else max(y0, y1)
                    xnovel = x1 if ynovel == y0 else x0
                    t = (xnovel, ynovel)
                    if t not in nodes:
                        nnid = len(nodes)
                        nodes[t] = nnid
                    else:
                        nnid = nodes[t]
                    new_segs.append([j0, nnid])
                    new_segs.append([nnid, j1])
    return new_segs, nodes

def _remove_duplicate_segments(new_segs):
    """Remove duplicate segments."""
    test = np.unique(np.sort(np.array(new_segs), axis=1), axis=0, return_index=True)[1]
    return np.array(new_segs)[test]

def _rebuild_new_nodes(nodes):
    """Rebuild the new_nodes array in case new nodes were added."""
    new_nodes = np.zeros((len(nodes), 2))
    for t, i in nodes.items():
        new_nodes[i, 0] = t[0]
        new_nodes[i, 1] = t[1]
    return new_nodes

def _connect_unlinked_nodes(new_nodes, new_segs, seed, verbose):
    """Connect unlinked nodes to the nearest node if needed."""
    seg_nodes = sorted(set([i for seg_pair in new_segs for i in seg_pair]))
    if verbose:
        print("Nb nodes in new segments:", len(seg_nodes))
        print("Nb nodes after simplification:", len(new_nodes))
        print("==============================")
    if len(new_nodes) != len(seg_nodes):
        print("WARNING: some nodes are not connected. Adding a conduit")
        plt.plot(new_nodes[:, 0], new_nodes[:, 1], 'k.')
        for i in range(len(new_nodes)):
            if i not in seg_nodes:
                plt.plot(new_nodes[i, 0], new_nodes[i, 1], 'ro')
                dist = (new_nodes[i, 0] - new_nodes[:, 0]) ** 2 + (new_nodes[i, 1] - new_nodes[:, 1]) ** 2
                dist[i] = np.max(dist)
                ind = np.argmin(dist)
                new_segs = np.append(new_segs, [[i, ind]], axis=0)
                print(f"node {i} not in the list of segments. We added a pipe to {ind}")
        print(f"The corresponding seed is {seed}!")
        print("==============================")
        seg_nodes = sorted(set([i for seg_pair in new_segs for i in seg_pair]))
        if len(new_nodes) != len(seg_nodes):
            print("WARNING: some nodes are not connected. Could not fix it.")
    return new_segs, seg_nodes

def _add_pipe_diameters(new_nodes, new_segs, bbox, seed):
    """Simulate the conduit radius using a Gaussian random field and add to segments."""
    node_str = []
    node_end = []
    newsegs_with_diam = np.copy(new_segs)
    for newseg in newsegs_with_diam:
        node_str.append(new_nodes[newseg[0]])
        node_end.append(new_nodes[newseg[1]])
    node_str = np.array(node_str)
    node_end = np.array(node_end)
    pos = get_position_center(node_str, node_end)
    pos = np.transpose(pos)
    sx, sy = 1, 1
    nx = bbox[1] + sx
    ny = bbox[3] + sy
    np.random.seed(seed)
    simu_grf = GRF_objt(typeModel='exponential')
    simu_grf.set_dimension(dimension=[nx, ny], spacing=[sx, sy])
    simu_grf.run_model_GRF(nb_real=1)
    simu_grf.grf_to_cdf()
    simu_grf.cdf_to_diameter()
    pipe_diameter = simu_grf.diameter[0]
    value_diameter = get_sgs_values(pipe_diameter, pos)
    newsegs_with_diam = np.array(newsegs_with_diam)
    newsegs_with_diam = np.c_[newsegs_with_diam, value_diameter]
    return newsegs_with_diam, pipe_diameter

#################
#Export conduits
#################

def export_conduits(nodes, segs, nodes_fname, edges_fname):
    # Select what to export
    v = segs
    n = nodes
    
    f = open(edges_fname,'w')
    for  i in range(len(v)):
        f.write("{:d} {:d} {}\n".format(np.int(v[i,0]),np.int(v[i,1]),v[i,2]))
    f.close()
    
    f = open(nodes_fname,'w')
    for i in range(len(n)):
        f.write("{} {} {}\n".format(i,n[i,0],n[i,1]))
    f.close()
    
    
##################
#Plot conduits
##################

def plot_conduits(nodes, segs, infile_nodes, infile_edges , scale=1.5, figsize=(20,8), case=1, seed=123):
    plt.figure(figsize=figsize)
    
    if len(infile_nodes)>0:
        tmp = np.loadtxt(infile_nodes)
        seg = np.loadtxt(infile_edges)
        x,y = tmp[:,1],tmp[:,2]
        plt.plot(x,y,'.b')
        
    for i in range(len(segs)):
        i0 = int(segs[i,0])
        i1 = int(segs[i,1])
        x1,x2,y1,y2 = nodes[i0,0], nodes[i1,0], nodes[i0,1], nodes[i1,1]
        plt.plot( [x1,x2], [y1,y2], '-k', linewidth=segs[i,2]*scale )
    plt.title('Karst network case {} seed {}'.format(case, seed), fontsize=20)
    plt.show()