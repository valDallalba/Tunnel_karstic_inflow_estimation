import numpy as np
import flopy as fp

#########################
# Fonctions utilitaires #
#########################

def run_all_mf6(kmatrix, nx, ny, sx, sy, sz, cst_bc, silent):
    """
    Initialise et exécute un modèle MODFLOW 6 avec les paramètres fournis.

    Paramètres
    ----------
    kmatrix : ndarray
        Matrice de conductivités hydrauliques.
    nx, ny : int
        Dimensions du modèle (nombre de colonnes et de lignes).
    sx, sy, sz : float
        Tailles des mailles selon x, y, z.
    cst_bc : list
        Conditions aux limites.
    silent : bool
        Si True, exécution silencieuse.

    Retour
    ------
    None
    """
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
    """
    Associe une conductivité laminaire dans les cellules karstiques.

    kmat_in : ndarray
        Matrice de conductivités initiale.
    conduit_rayon : ndarray
        Rayon des conduits par cellule.
    karst_binaire : ndarray
        Matrice binaire (1 = conduit, 0 = matrice).
    cell_size : float
        Taille de maille en mètres.
    Km : float
        Conductivité de la matrice.

    Retour
    ------
    kmat_out : ndarray
        Matrice de conductivités mise à jour.
    """
    klaminar = klam(conduit_rayon, d=cell_size, Km=Km)
    kmat_out = np.copy(kmat_in)
    kmat_out[karst_binaire != 0] = klaminar[karst_binaire != 0]
    return kmat_out


def klam(r, d=1, Rr=0.02, Km=1e-3):
    """
    Calcule la conductivité équivalente en régime laminaire.

    Paramètres
    ----------
    r : float ou ndarray
        Rayon du conduit [m].
    d : float
        Taille de cellule [m].
    Rr : float
        Rugosité (sans dimension).
    Km : float
        Conductivité de la matrice.

    Retour
    ------
    keq : float ou ndarray
        Conductivité équivalente.
    """
    # Paramètres fixes
    rho = 1000   # Masse volumique de l'eau [kg/m3]
    mu = 1e-3    # Viscosité dynamique [Pa.s]
    g = 9.81     # Gravité [m/s²]
    tau = 1      # Tortuosité

    # Formule : Darcy + Hagen-Poiseuille
    keq = (np.pi * r**2 / d**2) * (rho * g * r**2 / (8 * mu * tau)) \
          + ((d**2 - np.pi * r**2) / d**2) * Km
    return keq


def get_spdis(path):
    """
    Récupère le débit spécifique (SPDIS) dans un fichier CBC.

    Paramètres
    ----------
    path : str
        Chemin vers le fichier CBC.

    Retour
    ------
    spd : list
        Débits spécifiques.
    """
    spdobj = fp.utils.CellBudgetFile(path, precision='double')
    spd = spdobj.get_data(text="SPDIS")
    return spd


def get_velocity(cbc_flux, conduit_rayon, karst_binaire, sx, sy):
    """
    Calcule les vitesses dans les conduits.

    Paramètres
    ----------
    cbc_flux : dict
        Flux (qx, qy) extraits du CBC.
    conduit_rayon : ndarray
        Rayon des conduits.
    karst_binaire : ndarray
        Masque binaire (1 = conduit).
    sx, sy : float
        Dimensions de maille.

    Retour
    ------
    velocity : ndarray
        Module de la vitesse.
    """
    qx = cbc_flux['qx'].reshape(conduit_rayon.shape)
    qy = cbc_flux['qy'].reshape(conduit_rayon.shape)

    Qx = qx * sx * sy
    Qy = qy * sx * sy

    vx = Qx / (np.pi * conduit_rayon**2)
    vy = Qy / (np.pi * conduit_rayon**2)

    vx[karst_binaire == 0] = 0
    vy[karst_binaire == 0] = 0
    velocity = np.sqrt(vx**2 + vy**2)
    return velocity


def get_hds(path):
    """
    Charge les niveaux piézométriques (heads) d’un fichier HDS.

    Paramètres
    ----------
    path : str
        Chemin du fichier HDS.

    Retour
    ------
    head : ndarray
        Dernier état de charge.
    """
    hds = fp.utils.HeadFile(path)
    times = hds.get_times()
    head = hds.get_data(totim=times[-1])
    return head


def get_gradient(heads, karst_binaire, sx=1, sy=1):
    """
    Calcule le gradient hydraulique.

    Paramètres
    ----------
    heads : ndarray
        Matrice des charges.
    karst_binaire : ndarray
        Masque karstique.
    sx, sy : float
        Dimensions de cellule.

    Retour
    ------
    grad_h : ndarray
        Norme du gradient.
    """
    head_m = heads[0]
    grad_x = np.zeros(head_m.shape)
    grad_y = np.zeros(head_m.shape)

    # Gradient en x
    for i in range(head_m.shape[1] - 2):
        grad_x[:, i+1] = (head_m[:, i+2] - head_m[:, i]) / (2 * sx)
    grad_x[:, 0] = (head_m[:, 1] - head_m[:, 0]) / sx
    grad_x[:, -1] = (head_m[:, -1] - head_m[:, -2]) / sx

    # Gradient en y
    for j in range(head_m.shape[0] - 2):
        grad_y[j+1, :] = (head_m[j+2, :] - head_m[j, :]) / (2 * sy)
    grad_y[0, :] = (head_m[1, :] - head_m[0, :]) / sy
    grad_y[-1, :] = (head_m[-1, :] - head_m[-2, :]) / sy

    # Masquage hors conduits
    grad_x[karst_binaire == 0] = 0
    grad_y[karst_binaire == 0] = 0

    grad_h = np.sqrt(grad_x**2 + grad_y**2)
    return grad_h


###########################
# Classe Modflow 6        #
###########################

class Model_fp_6:
    """
    Classe de configuration et exécution d’un modèle MODFLOW 6.
    """

    def __init__(self, modelName='model_mf6',
                 path_ws='file_modflow_6/',
                 exe_path='wind_bin/mf6.exe',
                 silent_test=True):
        self.model_name = modelName
        self.path_ws = path_ws
        self.sim = fp.mf6.MFSimulation(
            exe_name=exe_path, version='mf6',
            sim_name=self.model_name, sim_ws=self.path_ws
        )
        self.silent_test = silent_test
        if not self.silent_test:
            print('*** Model is initialized! 0 ***')

    def set_dimension(self, nrow=200, ncol=2000,
                      nlay=1, delr=1., delc=1.,
                      ztop=0., zbot=-4):
        """ Définit les dimensions du modèle. """
        self.ny, self.nx, self.nlay = nrow, ncol, nlay
        self.sy, self.sx = delr, delc
        self.ztop, self.zbot = ztop, zbot
        self.lx, self.ly = ncol * delr, nrow * delc
        self.rch_info = 0
        if not self.silent_test:
            print('Dimensions are set! 1')

    def set_idomain(self, idomain=1):
        """ Définit l’IDomain (actif/inactif). """
        if idomain == 1:
            self.idomain = np.ones((self.nx, self.ny), dtype=int)
        else:
            self.idomain = idomain
        if not self.silent_test:
            print('Domain is set! 2')

    def init_model(self):
        """ Initialise la simulation et la discrétisation. """
        self.tdis = fp.mf6.ModflowTdis(self.sim, pname='tdis', time_units='seconds')
        self.gwf = fp.mf6.ModflowGwf(self.sim, modelname=self.model_name)
        self.solv = fp.mf6.ModflowIms(self.sim, print_option='SUMMARY', complexity='complex')
        self.grid = fp.mf6.ModflowGwfdis(
            self.gwf, nlay=self.nlay, nrow=self.ny, ncol=self.nx,
            delr=self.sy, delc=self.sx, top=self.ztop, botm=self.zbot,
            idomain=self.idomain
        )
        if not self.silent_test:
            print('Model is discretized! 3')

