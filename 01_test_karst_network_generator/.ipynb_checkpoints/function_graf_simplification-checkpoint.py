import numpy as np
import matplotlib.pyplot as plt
from function_assign_pipe_diameter import *


###############
#Simplify graph
###############

def simplify_graph(infile_nodes, infile_edges , bbox, dx, dy, verbose = True, seed = 12345678, test_rob=True):
    """
    Simplify a graph from pykasso given in two ascii files
    Returns list of nodes and coordinates
    """
    
    tmp = np.loadtxt(infile_nodes)
    seg = np.loadtxt(infile_edges)
    x,y = tmp[:,1],tmp[:,2]
    #plt.plot(x,y,'.')
    #plt.axis([0,10,0,10])
    #plt.show()
    
    # Round the coordinates of the initial nodes
    xa = np.floor( x/dx )*dx+dx/2
    ya = np.floor( y/dy )*dy+dy/2
    
    # Creates two dictionnaries for retrieving nodes rapidly and creating correspondance table
    nodes = {} # key = tuple with node coordinates, value = new node id
    nodes_corresp = {} # key = old node id, value = new node id
    
    # Loop over the old nodes
    j = 0 # Counters of new nodes id
    
    for i in range(len(xa)):
        t = (int(xa[i]),int(ya[i])) # a tuple representing a new node location
        
        if t in nodes.keys():  # This old node corresponds to a new node already defined
            nodes_corresp[i] = nodes[t]  # we just update the correspondance table
            
        else: # Else, we need to define a novel new node
            nodes[t] = j #  We store it in the dictionnary with a new id
            nodes_corresp[i] = j  # Update correspondance table
            j +=1  # Increase counter of ids
            
    # To facilitate the use of the new nodes, we create an array with nodes coordinates
    new_nodes = np.zeros((len(nodes),2))
    
    for t in nodes.keys():
        i = nodes[t]
        new_nodes[i,0]=t[0]
        new_nodes[i,1]=t[1]
        
    # Check if several nodes have the same positions
    tmp = []
    
    for i in range(len(new_nodes)):
        tmp.append( str(new_nodes[i,0])+"-"+str(new_nodes[i,1]) )
    nnc = sorted(set(tmp))
    
    if len(nnc) != len(new_nodes) :
        print("Warning: Several identical nodes after step 1")
        
    if verbose:
        print("Summary of graph simplification:")
        print("==============================")
        print("Nb nodes export pykasso:",len(xa))
        print("Nb nodes after simplification - step 1:",len(nodes))
        print("==============================")
        
    # Simplification of the links
    new_segs = [] # list of new segments
    
    for i in range(len(seg)): # Loop over the old segments
        j0 = nodes_corresp[ seg[i,0] ] # We get the new nodes id from old ids and correspondance table
        j1 = nodes_corresp[ seg[i,1] ]
        
        if j0 != j1: # If the two new nodes id are different, we should create a new segment
            j0,j1 = min((j0,j1)), max((j0,j1)) # Two avoid adding twice the same segment
            
            if [j0,j1] not in new_segs:
                x0,x1,y0,y1 = new_nodes[j0,0], new_nodes[j1,0], new_nodes[j0,1], new_nodes[j1,1]
                
                if (x0==x1) or (y0==y1): # The two nodes are aligned with the grid we just add segment
                    new_segs.append([j0,j1]) # Add new segments in the list
                    
                else: # The connection is along a diagonal
                    i0 = int(seg[i,0]) # We get the node ids in the original set of nodes
                    i1 = int(seg[i,1])
                    x0o, x1o, y0o, y1o =  x[i0],x[i1],y[i0],y[i1]  # We then get non rounded coordinates (original)
                    a = (y0o-y1o)/(x0o-x1o) # Compute the straight line through the original points
                    b = y0o-x0o*a
                    # Center of the rounded points / corresponds to the edges of the grid cells
                    xm,ym = (x0+x1)/2, (y0+y1)/2
                    # Intersection of the line between original points and egde and mid-point
                    ymp = a*xm+b
                    # Selection of the new point on the that should be connected
                    ynovel = min(y0,y1) if ymp<=ym else max(y0,y1)
                    xnovel = x1 if ynovel==y0 else x0
                    t = (xnovel,ynovel) # tuple with nodes coordinates
                    if t not in nodes.keys():
                        nnid = len(nodes) # New node id
                        nodes[t] = nnid
                    else:
                        nnid = nodes[t]
                    new_segs.append([j0,nnid]) # Add new segments in the list
                    new_segs.append([nnid,j1]) # Add new segments in the list
                    
    test = np.unique(np.sort(np.array(new_segs),axis=1),axis=0, return_index=True)[1]
    new_segs = np.array(new_segs)[test]
    
    # To facilitate the use of the new nodes, we create an array with nodes coordinates
    new_nodes = np.zeros((len(nodes),2))
    for t in nodes.keys():
        i = nodes[t]
        new_nodes[i,0]=t[0]
        new_nodes[i,1]=t[1]
        
    if test_rob:
        # This is for Rob: put the points on the boundary to the bounding box
        for i in range(len(new_nodes)):
            #if new_nodes[i,0] <= dx/2:
            #    new_nodes[i,0] = bbox[0]
            if new_nodes[i,1] > bbox[3]-dy:
                new_nodes[i,1] = bbox[3]

    # Check if several nodes have the same positions
    tmp = []
    
    for i in range(len(new_nodes)):
        tmp.append( str(new_nodes[i,0])+"-"+str(new_nodes[i,1]) )
    nnc = sorted(set(tmp))
    
    if len(nnc) != len(new_nodes) :
        print("Warning: Several identical nodes ! Should not happen - stop and debug !!!")
        print("The corresponding seed is {}!".format(seed))
        print("==============================")
        
    if verbose:
        print("Nb nodes after simplification - step 2:",len(nodes))
        print("==============================")
        print("Nb edges export pykasso:",len(seg))
        print("Nb edges after simplification:",len(new_segs))
        print("==============================")
        xa,ya = new_nodes[:,0],new_nodes[:,1]
        print( "Check nodes coordinates and bounding box")
        print( "Domain size:",bbox[0],bbox[1],bbox[2],bbox[3])
        print( "Original data range:",np.min(x), np.max(x), np.min(y) ,np.max(y) )
        print( "Simplified coordinates range:",np.min(xa), np.max(xa), np.min(ya) ,np.max(ya) )
        print( "Theoretical box:",bbox[0]+dx/2,bbox[1]-dx/2,bbox[2]+dy/2,bbox[3]-dy/2 )# Should be dx/2
        print("==============================")
        
    # Check which nodes are in the list of segments
    nnc = [] # List of nodes belonging to the segments
    
    for i in new_segs: # Take all the segments
        nnc.append(i[0])
        nnc.append(i[1])
    seg_nodes = sorted(set(nnc)) # Remove duplicates (a node may be connected to several segmennts)
    
    if verbose:
        print("Nb nodes in new segments:",len(seg_nodes))
        print("Nb nodes after simplification:",len(new_nodes))
        print("==============================")
        
    if len(new_nodes) != len(seg_nodes) :
        print("WARNING: some nodes are not connected. Adding a conduit")
        plt.plot(new_nodes[:,0],new_nodes[:,1],'k.')
        for i in range(len(new_nodes)):
            if i not in seg_nodes:
                plt.plot(new_nodes[i,0],new_nodes[i,1],'ro')
                dist = (new_nodes[i,0] - new_nodes[:,0])**2 + (new_nodes[i,1] - new_nodes[:,1])**2
                dist[i] = max(dist)
                ind = np.argmin(dist)
                new_segs = np.append(new_segs,[[i,ind]],axis=0)
                print("node",i,"not in the list of segments. We added a pipe to ",ind)
        print("The corresponding seed is {}!".format(seed))
        print("==============================")
        nnc = [] # List of nodes belonging to the segments
        for i in new_segs: # Take all the segments
            nnc.append(i[0])
            nnc.append(i[1])
        seg_nodes = sorted(set(nnc)) # Remove duplicates (a node may be connected to several segmennts)
        if len(new_nodes) != len(seg_nodes) :
            print("WARNING: some nodes are not connected. Could not fix it fuck")
            
    # Simulate the conduit radius
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
    sx = 1
    sy = 1
    nx = bbox[1]+sx
    ny = bbox[3]+sy
    np.random.seed(seed)
    simu_grf = GRF_objt(typeModel='exponential')
    simu_grf.set_dimension(dimension=[nx,ny],spacing=[sx,sy])
    simu_grf.run_model_GRF(nb_real=1)
    simu_grf.grf_to_cdf()
    simu_grf.cdf_to_diameter()
    pipe_diameter = simu_grf.diameter[0] 
    value_diameter = get_sgs_values(pipe_diameter, pos)
    newsegs_with_diam = np.array(newsegs_with_diam)
    newsegs_with_diam = np.c_[newsegs_with_diam,value_diameter]
    
    if verbose:
        print('Pipe diameter is added to newsegs array based on sampled sgs simulation!')
    return new_nodes, newsegs_with_diam, pipe_diameter


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