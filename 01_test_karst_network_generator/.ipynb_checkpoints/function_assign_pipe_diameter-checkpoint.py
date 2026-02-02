import numpy as np
from geone import grf as grf_geone
import geone.covModel as gcm
import geone.img as img
import geone.imgplot as imgplt
from scipy import interpolate
import math


#################
#GRF class
#################

class GRF_objt:
    '''
    Class GRF to create GRF simulation.
    '''
    def __init__(self, sillVal=1, rangeVal=[150,50], alpha=0, typeModel='gaussian'):
        self.model_vario= gcm.CovModel2D(elem=[(typeModel, {'w':sillVal, 'r':rangeVal})], 
                                         alpha=alpha, name='model_2D_karst')
        return

    
    def set_dimension(self, dimension=[2000,200], spacing=[1,1], origin=[0,0]):
        self.dimension = dimension
        self.spacing   = spacing
        self.origin    = origin
        return
    
    def run_model_GRF(self, nb_real=1):
        self.nb_real = nb_real
        self.simu    = grf_geone.grf2D(self.model_vario, self.dimension, self.spacing, 
                                 self.origin, nreal=self.nb_real, printInfo=False)
        #return print('Simulations are completed!')
        return
        
    def grf_to_cdf(self):
        self.cdf  = []
        quantiles = np.arange(0,1+0.01,0.01)
        
        for simu in self.simu:
            quantiles_val  = np.quantile(simu.flatten(), quantiles)
            cdf_fct        = interpolate.interp1d(quantiles_val, quantiles)
            proba_from_cdf = cdf_fct(simu)
            self.cdf.append(proba_from_cdf)
        return   
    
    
    def cdf_to_diameter(self, path_data='./Size_conduits.csv'):
        self.diameter = []
        data          = np.genfromtxt(path_data, delimiter=';')
        quantiles     = data[:,0]
        quantiles_val = np.sqrt(data[:,1]/math.pi) * 2
        diameter_fct  = interpolate.interp1d(quantiles, quantiles_val)
        
        for cdf in self.cdf:
            diameter = diameter_fct(cdf)
            self.diameter.append(diameter)
        #return print('Simulation from SGS are transformed to conduits diameters!')
        return
    
##################
#Fonctions
##################
    
def get_position_center(node_str, node_end):
    new_x = np.array((node_str[:,0]+node_end[:,0])/2)
    new_y = np.array((node_str[:,1]+node_end[:,1])/2)
    return np.array([new_x, new_y])


def get_sgs_values(sgs, position_center):
    value = []
    for pos in position_center:
        value.append(sgs[np.int(pos[1]), np.int(pos[0])])
    return np.array(value)


#############################
#############################