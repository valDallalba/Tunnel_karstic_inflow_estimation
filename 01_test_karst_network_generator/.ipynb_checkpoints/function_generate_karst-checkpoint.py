#functions
class Karst_networks:
    
    def __init__(self, sx=1, sy=1, nx=2000, ny=200):
        self.dim_x = nx * sx
        self.dim_y = ny * sy
        self.sx    = sx
        self.sy    = sy
        self.nx    = nx
        self.ny    = ny
        return
    
    def set_karst_input(self, nb_inlets_top=12, nb_inlets_rgt=5, out_position=[[0,0]]):
        '''
        out_position = [x,y]
        '''
        
        try:
            os.mkdir('model_info')
        except:
            print('')
            
        self.nb_inlets_top = nb_inlets_top
        self.nb_inlets_rgt = nb_inlets_rgt
        self.out_position  = np.array(out_position)
        
        x = np.abs(np.random.random(size=(self.nb_inlets_top,)) * self.dim_x - self.sx)
        y = np.ones((self.nb_inlets_top,)) *self.dim_y -self.sy
        self.inputs_top = [x,y]
        
        x_r = np.ones((self.nb_inlets_rgt,)) *self.dim_x - self.sx
        y_r = np.abs(np.random.random(size=(self.nb_inlets_rgt,)) *self.dim_y - self.sy)
        self.inputs_rgt = [x_r,y_r]
        
        with open('model_info/points_in.csv','w') as f:
            for node in range(self.nb_inlets_top):
                pos = '{} {}\n'.format(x[node], y[node])
                f.write(pos)
            for node in range(self.nb_inlets_rgt):
                pos = '{} {}\n'.format(x_r[node], y_r[node])
                f.write(pos)
        
        with open('model_info/points_out.csv','w') as f:
            for node in self.out_position:
                pos = '{} {}\n'.format(node[0], node[1])
                f.write(pos)
    
        return print('Outlets redefined')
    
    
    
    def set_pykasso_param(self, param_pykasso):
        self.path_param_pykasso = 'model_info/model1.yaml'
        with open(self.path_param_pykasso,'w') as f:
            data = yaml.dump(param_pykasso, f)
        return print('Pykasso parameters are sets!')
    
    def run_pykasso_sks(self):        
        self.sks = pk.SKS(self.path_param_pykasso)
        return print('Pykasso SKS is done!')  
    
    def run_pykasso_karst(self):        
        self.sks.compute_karst_network() 
        return print('Pykasso karst generation is done!')  
    
    def plot_fractures(self):
        d = self.sks.geology.data['fractures']['data']
        fig, ax1 = plt.subplots(1,1,figsize=(30,15))
        cmap='binary'
        im1 = ax1.imshow(d, extent=[self.sks.grid.x[0],self.sks.grid.x[-1],self.sks.grid.y[0],self.sks.grid.y[-1]], origin='lower',cmap=cmap)       
        plt.title('Fractures in the model')
        plt.xlabel('X')
        plt.xlabel('Z')
        plt.show()
        return None 

    def plot_karst(self):
        karst_network = self.sks.karst_simulations[-1]
        fig, ax1 = plt.subplots(1,1,figsize=(40,40))
        cmap='binary'
        im1 = ax1.imshow(karst_network.maps['karst'][-1], extent=[self.sks.grid.x[0],self.sks.grid.x[-1],self.sks.grid.y[0],self.sks.grid.y[-1]], origin='lower', cmap=cmap)

        plt.title('Karst in the model')
        plt.xlabel('X')
        plt.xlabel('Z')

        return None 

    def add_bed_fractures(self, nb_bed=10):
        self.nb_bed = nb_bed
        rows = np.round(np.random.random(nb_bed) * (self.sks.geology.data['fractures']['data'].shape[0] - 1))
        for row in rows:
            self.sks.geology.data['fractures']['data'][int(row)][:] = 1.
        return
    
    def update(self,seed):
        self.sks.set_rand_seed(seed)
        self.sks.update_inlets()
        self.sks.update_fractures()
        return
    
    def export(self,seed,case):
        #export the nodes and the edges files
        with open('pykasso_networks/nodes_'+case+'_'+ str(seed)+'.txt', 'w') as f: 
            for key, value in nagra.sks.karst_simulations[-1].network['nodes'].items(): 
                f.write('%s %s %s\n' % (key, value[0],value[1])) 
        with open('pykasso_networks/edges_'+case+'_'+ str(seed)+'.txt', 'w') as f: 
            for edges in nagra.sks.karst_simulations[-1].network['edges']: 
                    f.write('%s %s\n' % (edges[0],edges[1]))
        return