#Valentin Dall'alba

import os
import numpy as np
import matplotlib.pyplot as plt
import yaml
import pykasso as pk

from geone import grf as grf_geone
import geone.covModel as gcm
import geone.img as img
import geone.imgplot as imgplt

import scipy.stats
from scipy import interpolate
import matplotlib.pyplot as plt


class GRF_objt:
    '''
    Function to create and define GRF object for karst diameter assignation.
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
        return print('Simulations are completed!')
    
        
    def grf_to_cdf(self):
        self.cdf  = []
        quantiles = np.arange(0,1+0.00001,0.00001)
        
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
        return print('Simulation from SGS are transformed to conduits diameters!')


class KarstNetworks:
    """
    Classe pour générer et manipuler des réseaux karstiques
    en utilisant le moteur Pykasso.
    """

    def __init__(self, sx=1, sy=1, nx=2000, ny=200):
        """
        Initialise les dimensions du modèle.

        Parameters:
        - sx, sy : résolution spatiale en X et Y.
        - nx, ny : nombre de cellules en X et Y.
        """
        self.dim_x = nx * sx
        self.dim_y = ny * sy
        self.sx = sx
        self.sy = sy
        self.nx = nx
        self.ny = ny

    def set_karst_input(self, nb_inlets_top=12, nb_inlets_rgt=5, out_position=[[0, 0]]):
        """
        Définit les points d'entrée (inlets) et de sortie (outlets) du réseau karstique.

        Parameters:
        - nb_inlets_top : nombre d'entrées sur la face supérieure.
        - nb_inlets_rgt : nombre d'entrées sur la face droite.
        - out_position  : liste de positions [x, y] des sorties.
        """
        # Création du dossier temporaire pour stocker les fichiers
        os.makedirs('../model_info', exist_ok=True)

        self.nb_inlets_top = nb_inlets_top
        self.nb_inlets_rgt = nb_inlets_rgt
        self.out_position = np.array(out_position)

        # Génération des entrées en haut
        x_top = np.abs(np.random.random(nb_inlets_top) * self.dim_x - self.sx)
        y_top = np.ones(nb_inlets_top) * self.dim_y - self.sy
        self.inputs_top = [x_top, y_top]

        # Génération des entrées à droite
        x_right = np.ones(nb_inlets_rgt) * self.dim_x - self.sx
        y_right = np.abs(np.random.random(nb_inlets_rgt) * self.dim_y - self.sy)
        self.inputs_rgt = [x_right, y_right]

        # Sauvegarde des entrées
        with open('../model_info/points_in.csv', 'w') as f:
            for node in range(nb_inlets_top):
                f.write(f"{x_top[node]} {y_top[node]}\n")
            for node in range(nb_inlets_rgt):
                f.write(f"{x_right[node]} {y_right[node]}\n")

        # Sauvegarde des sorties
        with open('../model_info/points_out.csv', 'w') as f:
            for node in self.out_position:
                f.write(f"{node[0]} {node[1]}\n")

        print("Outlets redefined")

    def set_pykasso_param(self, param_pykasso):
        """
        Enregistre les paramètres YAML pour Pykasso.

        Parameters:
        - param_pykasso : dictionnaire des paramètres Pykasso.
        """
        self.path_param_pykasso = '../model_info/model1.yaml'
        with open(self.path_param_pykasso, 'w') as f:
            yaml.dump(param_pykasso, f)

        print("Pykasso parameters are set!")

    def run_pykasso_sks(self):
        """
        Initialise l’objet SKS de Pykasso avec les paramètres chargés.
        """
        self.sks = pk.SKS(self.path_param_pykasso)
        print("Pykasso SKS is initialized!")

    def run_pykasso_karst(self):
        """
        Lance la génération du réseau karstique.
        """
        self.sks.compute_karst_network()
        print("Pykasso karst generation is done!")

    def plot_fractures(self):
        """
        Affiche la carte des fractures du modèle.
        """
        fractures_data = self.sks.geology.data['fractures']['data']
        plt.figure(figsize=(30, 15))
        plt.imshow(
            fractures_data,
            extent=[self.sks.grid.x[0], self.sks.grid.x[-1], self.sks.grid.y[0], self.sks.grid.y[-1]],
            origin='lower',
            cmap='binary'
        )
        plt.title("Fractures in the model")
        plt.xlabel("X")
        plt.ylabel("Z")
        plt.show()

    def plot_karst(self):
        """
        Affiche la carte du réseau karstique simulé.
        """
        karst_network = self.sks.karst_simulations[-1]
        plt.figure(figsize=(40, 40))
        plt.imshow(
            karst_network.maps['karst'][-1],
            extent=[self.sks.grid.x[0], self.sks.grid.x[-1], self.sks.grid.y[0], self.sks.grid.y[-1]],
            origin='lower',
            cmap='binary'
        )
        plt.title("Karst in the model")
        plt.xlabel("X")
        plt.ylabel("Z")
        plt.show()

    def add_bed_fractures(self, nb_bed=10):
        """
        Ajoute des fractures horizontales ("bed fractures") au modèle.

        Parameters:
        - nb_bed : nombre de lits fracturés à ajouter.
        """
        self.nb_bed = nb_bed
        fractures_data = self.sks.geology.data['fractures']['data']
        rows = np.round(np.random.random(nb_bed) * (fractures_data.shape[0] - 1))

        for row in rows:
            fractures_data[int(row), :] = 1.0

    def update(self, seed):
        """
        Met à jour le modèle avec une nouvelle graine aléatoire.

        Parameters:
        - seed : valeur de la graine.
        """
        self.sks.set_rand_seed(seed)
        self.sks.update_inlets()
        self.sks.update_fractures()

    def export(self, seed, case):
        """
        Exporte les nœuds et les arêtes du réseau karstique simulé.

        Parameters:
        - seed : graine utilisée pour la simulation.
        - case : nom ou identifiant du cas.
        """
        os.makedirs('pykasso_networks', exist_ok=True)

        # Export des nœuds
        with open(f'pykasso_networks/nodes_{case}_{seed}.txt', 'w') as f:
            for key, value in self.sks.karst_simulations[-1].network['nodes'].items():
                f.write(f"{key} {value[0]} {value[1]}\n")

        # Export des arêtes
        with open(f'pykasso_networks/edges_{case}_{seed}.txt', 'w') as f:
            for edge in self.sks.karst_simulations[-1].network['edges']:
                f.write(f"{edge[0]} {edge[1]}\n")
