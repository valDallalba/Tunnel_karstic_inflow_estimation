import numpy as np
import pandas as pd

from geone import grf as grf_geone

import pykasso as pk

import flopy as fp
from flopy.utils.gridintersect import GridIntersect

import shapely
from shapely.geometry import LineString


#########################
#########################

def run_all_mf6(kmatrix, nx, ny, sx, sy, sz, cst_bc, silent):
    model = Model_fp_6()
    model.set_dimension(nrow=ny, ncol=nx, delr=sy, delc=sx, zbot=-sz)
    model.set_idomain()
    model.init_model()
    model.set_bdc(bdc=cst_bc)
    model.init_bdc()
    model.set_kmatrix(kmatrix=kmatrix)
    model.init_node_prop()
    model.set_initial_cond()
    model.init_initial_cond()
    model.init_oc()
    model.write_input_file()
    model.run_simulation(silent=silent)
    return    


def assign_k_laminar(kmat_in, conduit_rayon, karst_binaire, cell_size=1, Km=1e-3):
    '''
    km is the k matrix
    conduit is the size of the conduit assign on a matrix shape
    '''
    klaminar = klam(conduit_rayon, d=cell_size, Km=Km)
    kmat_out = np.copy(kmat_in)
    kmat_out[karst_binaire!=0] = klaminar[karst_binaire!=0]
    
    return kmat_out


def klam(r, d=1, Rr=0.02, Km=1e-3):
    ''' 
    Input parameters
    d = 1 # Cell size in [m] (assumed to be cubic)
    r = 0.01 # Conduit radius in [m]
    Rr = 0.02 # Rugosity [dimensionless]
    Km = 1e-3 # Matrix conductivity
    '''
    
    #  Fixed parameters
    rho = 1000 # Fluid density
    mu = 1e-3 # Fluid visocisty
    g = 9.81 # Acceleration of gravity
    tau = 1 # Tortuosity
    
    return np.pi*r**2/d**2 *rho*g*r**2/8/mu/tau + (d**2-np.pi*r**2)/(d**2)*Km


def get_spdis(path):
    """
    Function that returns the specific discharge from the cbcfile
    """
 
    spdobj = fp.utils.CellBudgetFile(path, precision='double')  
    spd    = spdobj.get_data(text="SPDIS")
    
    return spd


def get_velocity(cbc_flux, conduit_rayon, karst_binaire, sx, sy):
    qx = (cbc_flux['qx'].reshape((conduit_rayon.shape[0],conduit_rayon.shape[1])))
    qy = (cbc_flux['qy'].reshape((conduit_rayon.shape[0],conduit_rayon.shape[1])))

    Qx = qx * sx * sy
    Qy = qy * sx * sy
    
    vx = Qx / (np.pi * conduit_rayon**2)
    vy = Qy / (np.pi * conduit_rayon**2)

    vx[karst_binaire==0] = 0
    vy[karst_binaire==0] = 0
    velocity = np.sqrt(vx**2+vy**2)
    
    return velocity


def get_hds(path):
    hds    = fp.utils.HeadFile(path)
    times  = hds.get_times()
    head   = hds.get_data(totim=times[-1])
    
    return head


def get_gradient(heads, karst_binaire, sx=1, sy=1):
    head_m = heads[0]
    
    grad_x = np.zeros(head_m.shape)
    
    for i in range(head_m.shape[1]-2):
        grad_x[:,i+1] = (head_m[:,i+2] - head_m[:,i]) / (2*sx)
    grad_x[:,0]  = (head_m[:,1] - head_m[:,0]) / (sx)
    grad_x[:,-1] = (head_m[:,-1] - head_m[:,-2]) /(sx)
            
    grad_y = np.zeros(head_m.shape)
    
    for j in range(head_m.shape[0]-2):
        grad_y[j+1,:] = (head_m[j+2,:] - head_m[j,:]) / (2*sy)
    grad_y[0,:]  = (head_m[1,:] - head_m[0,:]) /(sy)
    grad_y[-1,:] = (head_m[-1,:] - head_m[-2,:]) /(sy)
    
    grad_x[karst_binaire==0] = 0
    grad_y[karst_binaire==0] = 0
    
    grad_h  = np.sqrt(grad_x**2+grad_y**2)
    
    return grad_h


def get_gradient_y(heads, conduit, karst_binaire, sx=1, sy=1):
    head_m = heads[0]
    grad_x = np.zeros(head_m.shape)
    #for i in range(head_m.shape[1]-2):
    #    grad_x[:,i+1] = (head_m[:,i+2] - head_m[:,i]) / (2*sx)
    #grad_x[:,0]  = (head_m[:,1] - head_m[:,0]) / (sx)
    #grad_x[:,-1] = (head_m[:,-1] - head_m[:,-2]) /(sx)
            
    grad_y = np.zeros(head_m.shape)
    
    for j in range(head_m.shape[0]-2):
        grad_y[j+1,:] = (head_m[j+2,:] - head_m[j,:]) / (2*sy)
    grad_y[0,:]  = (head_m[1,:] - head_m[0,:]) /(sy)
    grad_y[-1,:] = (head_m[-1,:] - head_m[-2,:]) /(sy)
    
    grad_x[karst_binaire==0] = 0
    grad_y[karst_binaire==0] = 0
    
    grad_h  = np.sqrt(grad_x**2+grad_y**2)
    return grad_h


def reynolds(v, diameter, limite_lam=2000):
    rey_mat  = (v * diameter * 1000) / 1e-3
    rey_test = np.ones(v.shape)
    rey_test[rey_mat<limite_lam] = 0
    
    return rey_mat, rey_test


def kturb(conduit_rayon, karst_binaire, cell_size=1, Rr=0.02, Km=1e-5, gradh=1e-3):
    ''' 
    Input parameters
    d = 1 # Cell size in [m] (assumed to be cubic)
    r = 0.01 # Conduit radius in [m]
    Rr = 0.02 # Rugosity [dimensionless]
    Km = 1e-3 # Matrix conductivity
    gradh = 1e-3 # Head gradient (assuming it is prescribed by boundary conditions)
    '''
    
    #  Fixed parameters
    r = conduit_rayon
    d = cell_size
    rho = 1000 # Fluid density
    mu = 1e-3 # Fluid visocisty
    g = 9.81 # Acceleration of gravity
    tau = 1 # Tortuosity
    cst=1e-7
    
    kt = np.pi*r**2/d**2*4*np.sqrt(g*r)*np.log10(3.7/Rr)/np.sqrt(gradh) + cst
    kt[karst_binaire==0] = cst
    kt[np.isnan(kt)] = cst
    kt[kt==np.inf] = cst
    
    return kt


def update_km_old(kmat_old, k_turb, karst_binaire):
    kmat_to_upt = np.copy(kmat_old)
    kmat_to_upt[karst_binaire==1] = k_turb[karst_binaire==1]
    
    return kmat_to_upt


def update_km(kmat_old, k_turb, karst_binaire, div=2):
    k_upd = np.copy(kmat_old)
    k_bin = k_turb-kmat_old
    k_bin[k_bin<0] = -1
    k_bin[k_bin>0] = 1
    
    k_diff = np.abs(k_turb-kmat_old)
    k_new  = kmat_old+(k_diff/div)*k_bin
    k_upd[karst_binaire==1] = k_new[karst_binaire==1]
    
    return k_upd


def get_outflow(position_out):
    out = np.flipud(get_spdis('file_modflow_6/model_mf6.cbc')[0]['qy'].reshape(ny,nx))
    
    return out[position_out[0],position_out[1]]



###########################
###########################

class Model_fp_6():
    def __init__(self, modelName='model_mf6', path_ws='file_modflow_6/', exe_path='wind_bin/mf6.exe', silent_test=True):
        self.model_name = modelName
        self.path_ws    = path_ws
        self.sim        = fp.mf6.MFSimulation(exe_name=exe_path, version='mf6',sim_name=self.model_name, sim_ws=self.path_ws )        
        self.silent_test = silent_test

        if silent_test is False :
            print('*** Model is initialized! 0 ***')
        return 
    
    ########
    ########
    
    def set_dimension(self, nrow=200, ncol=2000, nlay=1, delr=1., delc=1., ztop=0., zbot=-4):
        self.ny = nrow
        self.nx = ncol
        self.nlay = nlay
        self.sy = delr
        self.sx = delc
        self.ztop = ztop
        self.zbot = zbot
        self.lx   = ncol*delr    
        self.ly   = nrow*delc
        self.rch_info = 0
        if self.silent_test is False :
            print('Dimensions are set! 1')   
        return    
    
    #######
    #######
    def set_idomain(self,idomain=1):
        if idomain==1:
            self.idomain = np.ones((self.nx, self.ny), dtype=np.int)
        else:
            self.idomain=idomain

        if self.silent_test is False :
            print('Domain is set! 2')
        return 
        
    def init_model(self):
            
        self.tdis  = fp.mf6.ModflowTdis(self.sim, pname='tdis', time_units='seconds')     
        self.gwf   = fp.mf6.ModflowGwf(self.sim, modelname=self.model_name)
        self.solv  = fp.mf6.ModflowIms(self.sim, print_option='SUMMARY',complexity='complex')                    #Solver
        self.grid  = fp.mf6.ModflowGwfdis(self.gwf, nlay=self.nlay, nrow=self.ny, ncol=self.nx,
                               delr=self.sy, delc=self.sx, top=self.ztop, botm=self.zbot, idomain=self.idomain)
        if self.silent_test is False :
            print('Model is discretized! 3')
        return
    
    ###########
    def set_bdc(self, bdc=None):
        self.bdc = bdc
        if self.silent_test is False :
            print('BDC are set! 4')
        return 
    
    def init_bdc(self):
        if self.bdc is not None:
            for i, bd in enumerate(self.bdc):
                fp.mf6.ModflowGwfchd(self.gwf, pname=str(i), filename='bdc_{}'.format(i), stress_period_data={0:bd}, save_flows=True)
        
        if self.silent_test is False :
            print('BDC are initialized! 4_bis')
        return 
    
    ###########
    def set_kmatrix(self, kmatrix=None):
        if kmatrix is None:
            self.kmatrix = np.ones((self.ny,self.nx)) *1e-3
        else:
            self.kmatrix = kmatrix

        if self.silent_test is False :
            print('k matrix is set! 5')
        return 
            
    def init_node_prop(self):
        self.npf = fp.mf6.ModflowGwfnpf(self.gwf, k=self.kmatrix, save_flows=True, save_specific_discharge=True)
        if self.silent_test is False :
            print('Node properties are initialiazed! 6')
        return 
    ###########
    
    def set_initial_cond(self, ic=1):
        if ic==1:
            self.ic=np.ones((self.ny, self.nx))
        else:
            self.ic=ic
        if self.silent_test is False :
            print('Initial condition are set! 7')
        return
    
    def init_initial_cond(self):
        self.gwic = fp.mf6.ModflowGwfic(self.gwf,strt=self.ic)
        if self.silent_test is False :
            print('Initial condition are initialiazed! 8')
        return 

    def init_recharge(self):
        self.rch = fp.mf6.ModflowGwfrcha(self.gwf, recharge=0, save_flows=True)
    ###########
    
    def init_oc(self):
        self.oc = fp.mf6.ModflowGwfoc(self.gwf, head_filerecord='model_mf6.hds',budget_filerecord='model_mf6.cbc',saverecord=[['budget','all'],['head','all']], printrecord=[['head','all']])
        
        if self.silent_test is False :
            print('Output control is initialiazed! 9')
        return 

    def write_input_file(self):
        self.sim.write_simulation(silent=True)
        if self.silent_test is False :
            print('Input file are written! 10')
        return 

    def run_simulation(self,silent):
        self.sim.run_simulation(silent=silent)
        if self.silent_test is False :
            print('Simulation is completed! 11')
        return 