import numpy as np
import pandas as pd
import pickle 
import os
import shutil


def write_geoheight(newX, dir):
    with open(dir+'/geoheight.txt','w') as f:
        for i in range(len(newX)):
            if i == len(newX)-1:
                f.write('{}\t-1'.format(i+1))
            else:
                f.write('{}\t-1\n'.format(i+1))
                
            
def write_kexch(newX, dir):
    with open(dir+'/K_exch.txt','w') as f:
        for i in range(len(newX)):
            if i == len(newX)-1:
                f.write('{}\t5'.format(i+1))
            else:
                f.write('{}\t5\n'.format(i+1))
                

def write_node_heads(newX, dir, node_out=99):
    with open(dir+'/node_head.txt','w') as f:
        for i in range(len(newX)):
            if i+1 == node_out and node_out!=len(newX):
                f.write('{}\t0\n'.format(i+1))
            elif i+1 == node_out:
                f.write('{}\t0'.format(i+1))
            elif i+1 == len(newX):
                f.write('{}\t-1'.format(i+1))
            else:
                f.write('{}\t-1\n'.format(i+1))

                
def write_pipe_info(newSeg, dir, tortuosity=1, rugosity=0.01, Remin=2000, Remax=3000):
    with open(dir+'/pipe_info.txt','w') as f:
        for i in range(len(newSeg)):
            diam = newSeg[i][2]
            if i+1 != len(newSeg):
                f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(i+1, diam, tortuosity, rugosity, Remin, Remax))
            else:
                f.write('{}\t{}\t{}\t{}\t{}\t{}'.format(i+1, diam, tortuosity, rugosity, Remin, Remax))    
    return


def write_network_info(newX, newsegs, dir, dx, dy):
    with open(dir+'/network_info.txt','w') as f:

        for i in range(len(newX)):      
            X = int(newX[i][0]/dx)+1
            if X <= 0:
                X = 1
            Y =  int(newX[i][1]/dy)+1
            if Y <= 0:
                Y = 1
            id_connect, pipe_connect = find_voisin(newsegs, i)
            
            if i+1 != len(newX):
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(i+1, X, Y, 1, 
                                                            id_connect[0], id_connect[1], id_connect[2], id_connect[3],
                                                            id_connect[4], id_connect[5], pipe_connect[0], pipe_connect[1],
                                                            pipe_connect[2], pipe_connect[3], pipe_connect[4], pipe_connect[5]))
            else:
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(i+1, X, Y, 1, 
                                                            id_connect[0], id_connect[1], id_connect[2], id_connect[3],
                                                            id_connect[4], id_connect[5], pipe_connect[0], pipe_connect[1],
                                                            pipe_connect[2], pipe_connect[3], pipe_connect[4], pipe_connect[5]))
    return


def find_voisin(newsegs, id_node):
    id_connect   = []
    pipe_connect = []
    
    for i, seg in enumerate(newsegs):
        if seg[0] == id_node:
            id_connect.append(int(seg[1]+1))
            pipe_connect.append(i+1)
            
        if seg[1] == id_node:
            id_connect.append(int(seg[0]+1))
            pipe_connect.append(i+1)

    while len(id_connect)!=6:
        id_connect.append(0)
        
    while len(pipe_connect)!=6:
        pipe_connect.append(0)
        
    return id_connect, pipe_connect

