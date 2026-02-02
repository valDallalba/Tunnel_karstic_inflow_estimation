import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import flopy as fp
import flopy.utils.binaryfile as bf
from flopy.utils.gridintersect import GridIntersect

import mfcfp
import os
import shutil


    
class Model_cfp:
    '''
    Class to run create and run 2D model using modflow 2005 cfp
    '''
    
    #######
    #Init the model
    #######
    
    def __init__(self, modelName='model_mf5cfp', path_ws='.mf5_cfp_files/'):
        self.model_name = modelName
        self.path_ws    = path_ws
        self.model      = fp.modflow.Modflow(modelname=self.model_name, exe_name='mf2005cfp.exe')
        return print('*** Model is initialized! 0 ***')
    
    
    ########
    #Init the dimension of the grid
    ########
    
    def set_dimension(self, nrow=5, ncol=5, nlay=1, delr=100., delc=100., ztop=0., zbot=-4):
        self.nrow = nrow
        self.ncol = ncol
        self.nlay = nlay
        self.delr = delr
        self.delc = delc
        self.ztop = ztop
        self.zbot = zbot
        self.lx   = ncol*delr    
        self.ly   = nrow*delc
        self.rch_info = 0
        return print('Dimensions are set! 1')
    
    def set_time_periode(self, nbPer=1, lenPer=1, nbStp=[1], steady=[True]):
        self.nbPer  = nbPer
        self.lenPer = lenPer
        self.nbStp  = nbStp
        self.steady = steady
        return print('Time period is set! 2')

    def init_discretization(self):
        #lower-left corner of the grid at (0, 0)
        self.dis = fp.modflow.ModflowDis(model=self.model, nlay=self.nlay, nrow=self.nrow, ncol=self.ncol, nper=self.nbPer, 
                                            delc=self.delc, delr=self.delr, top=self.ztop, botm=self.zbot, perlen=self.lenPer,
                                            nstp=self.nbStp, steady=self.steady,  xul=0+self.ncol*self.delc, yul=0+self.nrow*self.delr)
        return print('*** Model is discretized! 3 ***')
    
    
    ########
    #Init the boundary condition and the constance head
    ########
    
    def set_bdc_type(self, matrix_bound=False):
        #ibound < 0 : cst head 
        #ibound > 0 : variable head
        #ibound = 0 : no flow
        if matrix_bound is False:
            self.bdc_type = np.ones((self.nlay, self.nrow, self.ncol),dtype=np.int32)
        else:
            self.bdc_type = matrix_bound                      
        return print('Boundary cell types are set! 4')
        
    def set_bdc_val(self, matrix_ini_bdc=False, val=1):
        if matrix_ini_bdc is False:
            self.bdc_val = np.ones((self.nlay, self.nrow, self.ncol),dtype=np.int32)*val
        else:
            self.bdc_val = val
        return print('Boundary cell values are set! 5')
    
    def init_bdc(self):
        self.bdc = fp.modflow.ModflowBas(model=self.model, ibound=self.bdc_type, strt=self.bdc_val)
        return print('*** Boundary conditions are initialized! 6 ****')
    
    
    ########
    #Init flow properties
    ########    
    
    def set_flow_prop(self, layerType=0, cstAni=1, layVka=0, hK=86.4, vK=86.4, ss=1e-6, sy=0.1):
        #layerType = 0 : confine layer
        #layerType > 0 : unconfine layer/ convertible layer
        #layerType < 0 : convertible unless THICKSTART option is set
        #cstAni > 0 : cstAni is the horizontal anisotropy
        #cstAni <= 0 : must define the hani variable that contain the anisotropy
        #layVka = 0 : vka = vertcial hydraulic conductivity
        #layVka != 0 : vka = ration of horizontal to vertical hydraulic conductivity
        #hKa and vKa are in [m/j] (86.4 [m/j] = 0.001 [m/s])
        self.layerType = layerType
        self.cstAni    = cstAni
        self.layVka    = layVka
        self.hK = hK
        self.vK = vK
        self.ss = ss
        self.sy = sy
        return print('Flow properties are defined! 7')
    
    def init_flow_prop(self):
        self.flow_p = fp.modflow.ModflowLpf(model=self.model, laytyp=self.layerType, chani=self.cstAni,
                                               layvka=self.layVka, hk=self.hK, vka=self.vK, ss=self.ss, 
                                               sy=self.sy, ipakcb=53,laywet=0,hdry=999.,wetdry=0)
        return print('*** Flow properties are initialized! 8 ***')
    
    
    ########
    #Init output control files
    ########
    
    def init_output_control(self):
        spd = {(0,0):['print head', 'save head', 'save budget']} 
        self.oc  = fp.modflow.ModflowOc(self.model, stress_period_data=spd)   
        self.pcg = fp.modflow.ModflowPcg(model=self.model) #convergence criteria
        return print('*** Output control parameters are initialized! 9 ***')
    
    
    ########
    #Init recharge (non activate in our case)
    ########
    
    def set_recharge(self,valueRech=0.002):
        self.nrchop = 1
        #self.inrech = 1
        rech = {}
        rech[0] = np.ones([self.nrow, self.ncol]) * valueRech
        self.valueRech = rech[0]
        return print('Recharge is set!')

    def init_recharge(self):
        self.rch_info = 1
        self.rch = fp.modflow.mfrch.ModflowRch(self.model, nrchop=self.nrchop, rech=self.valueRech)
        return print('*** Recharge is initialized! ****')
    
    
    ########
    #Init cfp node properties
    ########
    
    def set_cfp_node_oc(self, path):
        node_head = np.genfromtxt(path+'node_head.txt')
        self.nb_node   = len(node_head)
        self.list_node = np.arange(1,self.nb_node+1)
        self.n_nts     = 1  

        pipe_info = np.genfromtxt(path+'pipe_info.txt')
        self.nb_pipe   = len(pipe_info)
        self.list_pipe = np.arange(1,self.nb_pipe+1)
        self.t_nts     = 1
        self.path_conduits = path
        return print('Output control nodes are set! 10')
    
    def init_cfp_node_oc(self, nb_node_to_save, node_nb_to_save, nb_pipes_to_save, pipe_nb_to_save):
        self.coc_unit_num = 16
        self.coc =  mfcfp.ModflowCoc(nnodes=nb_node_to_save, node_nums=node_nb_to_save, n_nts=self.n_nts,
                                     npipes=nb_pipes_to_save, pipe_nums=pipe_nb_to_save, t_nts=self.t_nts)
        return print('*** Output control nodes are initialized! 11 ***')
    
    
    ########
    #Init conduit recharge properties (doline like) set to 0 here
    ########
    
    def set_conduit_recharge_prop(self, p_crch=[0.00]):
        self.iflag_crch = [1]
        self.spers      = np.arange(1,self.nbPer+1)
        self.p_crch     = p_crch*self.nb_node
        return print('Conduit recharge properties are set! 12')
    
    def init_conduit_recharge_prop(self):
        self.crch_unit_num = 18  
        self.crch = mfcfp.ModflowCrch(self.list_node, self.spers, self.iflag_crch, self.p_crch)
        return print('*** Conduit recharge properties are initialized 13 ***')
    
    
    ########
    #Init conduit flow properties
    ########
    
    def set_conduit_flow_prop(self, temp=25.):
        #sa_exch : integer flag, 1 = assign conduit wall permeability & let cfp comput pipe surface area
        self.temp    = temp
        self.sa_exch = 1 #to modify I guess if we want to modify pipe area/diameter
        return print('Conduit flow properties are set! 14')
    
    def init_conduit_flow_prop(self):
        self.cfp_unit_num = 17  
        self.cfp = mfcfp.ModflowCfp(nnodes=self.nb_node, npipes=self.nb_pipe, nlay=self.nlay,
                                    network_info_file=self.path_conduits+'network_info.txt', geoheight_file=self.path_conduits+'geoheight.txt',
                                    pipe_info_file=self.path_conduits+'pipe_info.txt', node_head_file=self.path_conduits+'node_head.txt',
                                    K_exch_file=self.path_conduits+'K_exch.txt', mode=1, temp=self.temp, sa_exch=self.sa_exch, p_nr=1)
        return print('*** Conduit flow properties are initialized! 15***')
     
        
    ########
    #Write input files
    ########   
    
    def write_input_files(self):
        self.model.write_input()
        return print('*** Input files are written! 16***')
    
    def update_input_files(self):
        filenames = ['coc', 'cfp', 'crch']     
        dataset_strings = [self.coc, self.cfp, self.crch]     
        mfcfp.cfp_write_input(self.model_name, dataset_strings, filenames)
        mfcfp.update_nam(self.model_name, self.coc_unit_num, self.cfp_unit_num, self.crch_unit_num)
        return print('*** Input files are modified! 17***')
    
    
    ########
    #Run simulation
    ########
    
    def run_simu(self):
        self.model.run_model()        
        return print('*** Model is run! 18***')
        
        
    ########
    #Clean repo 
    ########
    
    def clean_repo(self, path='mf5_cfp_files/'):
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
            
        os.mkdir(path)
        
        all_files = os.listdir()
        for name in all_files:
            if name[:len('model_mf5_cfp')-1]=='model_mf5cfp':
                os.rename(name,'mf5_cfp_files/'+name)
        return print('Repository is clean!')
    
    
    def clean_repo_old(self):
        list_ext = ['.bas', '.cbc', '.cfp','.coc','.crch','.dis','.hds','.list','.lpf','.nam','.oc','.pcg']
        
        if self.rch_info == 1:
            list_ext.append('.rch')
            
        if os.path.isdir('file_modflow'):
            shutil.rmtree('file_modflow/', ignore_errors=True)

        if os.path.isdir('output_node'):
            shutil.rmtree('output_node/', ignore_errors=True)

        if os.path.isdir('output_pipe'):
            shutil.rmtree('output_pipe/', ignore_errors=True)

        os.mkdir('file_modflow')
        os.mkdir('output_node')
        os.mkdir('output_pipe')
            
        for  ext in list_ext:
            os.rename(self.model_name+ext,'file_modflow/'+self.model_name+ext)
        for i in range(self.nb_node):
            os.rename('NODE{:04d}.OUT'.format(i+1),'output_node/NODE{:04d}.OUT'.format(i+1))
        for i in range(self.nb_pipe):
            os.rename('TUBE{:04d}.OUT'.format(i+1),'output_pipe/TUBE{:04d}.OUT'.format(i+1))
        return('*** Depositories are cleaned! 19***')
    
    
    ########
    #Some plots
    ########
    
    def plot_heads(self, path='mf5_cfp_files/'):
        hds    = bf.HeadFile(path+self.model_name+'.hds')
        times  = hds.get_times()
        head   = hds.get_data(totim=times[-1])
        #plt.contour(head[0, :, :])
        #plt.colorbar()
        plt.figure(figsize=(15,15))
        plt.imshow(head[0], origin='bottom')
        plt.title('Hydraulic heads')
        plt.colorbar(shrink=0.1)
        plt.tight_layout()
        plt.axis("scaled")
        plt.show()
        return head

    def plot_contour(self, path='./'):
        hds    = bf.HeadFile(path+self.model_name+'.hds')
        times  = hds.get_times()
        head   = hds.get_data(totim=times[-1])
        plt.contour(head[0, :, :])
        plt.colorbar()
        plt.figure(figsize=(15,15))
        #plt.imshow(head[0], origin='bottom')
        plt.title('Hydraulic heads')
        plt.colorbar(shrink=0.1)
        plt.tight_layout()
        plt.axis("scaled")
        plt.show()
        return head
