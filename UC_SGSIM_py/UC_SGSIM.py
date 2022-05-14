import UC_SGSIM_py as UC
from UC_SGSIM_py.Krige.kriging import Kriging
from UC_SGSIM_py.Plot.Plot import Visualize
import numpy as np
from multiprocessing import Pool
import time

class Simulation:
    
    def __init__(self, Y, model, nR, randomseed, krige_method='SimpleKrige'):
        self.Y = Y
        self.model = model
        self.nR = nR
        self.hs=model.hs
        self.bw=model.bw
        self.randomseed = randomseed
        self.krige_method = krige_method
        self.size = len(self.Y)
        self.RandomField = np.empty([self.size, self.nR])

    def compute(self, randomseed, parallel = False):
        
        if parallel == True:
            self.randomseed = randomseed

        initial_seed = self.randomseed

        if self.krige_method == 'SimpleKrige':

            self.krige = UC.SimpleKrige(self.model)

        counts = 0

        start_time = time.time()

        while counts < self.nR // self.n_process:
            
            boundary_constrained = 0 
            unsampled = np.linspace(0, self.size-1, self.size)
            
            if boundary_constrained == 0:
                y_value = np.random.normal(0,1,2).reshape(2,1)
                x_grid = np.array([0, self.size-1]).reshape(2,1)
                Z = np.zeros(self.size)
                Z[0], Z[-1] = y_value[0], y_value[1]
                unsampled = np.delete(unsampled,[0,-1])
                neigh = 0 
            else:
                y_value = np.random.normal(0,1,1).reshape(1,1)
                ridx = np.random.randint(0,self.size-1,1)
                x_grid = np.array([ridx]).reshape(1,1)
                Z = np.zeros(self.size)
                Z[ridx] = y_value[0]
                unsampled = np.delete(unsampled,[ridx])
                neigh = 1

            L = np.hstack([x_grid, y_value])

            np.random.seed(self.randomseed)
            randompath = np.random.choice(unsampled, len(unsampled), replace=False)
            
            for i in range(len(unsampled)):

                Z[int(randompath[i])] = self.krige.compute(L, randompath[i], neigh, self.randomseed)
            
                temp = np.hstack([randompath[i],Z[int(randompath[i])]])
                L = np.vstack([L,temp])
                
                if neigh<6:
                    neigh += 1
                
                self.randomseed += 1

            Z_Gap = abs(Z.max()-Z.min())

            if 2<Z_Gap<=6.5:
                self.RandomField[:,counts] = Z
                counts = counts+1
                #print('Progress = %.2f' % (counts/self.nR*100)+'%', end='\r')

            self.randomseed += 1
            
        #print('Progress = %.2f' % 100+'%\n', end='\r')
        
        end_time = time.time()
        
        #print('Time = %f'%(end_time-start_time),'s\n')
        #print("Last RandomSeed = %d" %(self.randomseed),'\n')
        #print("RandomSeed passed = %d" %(self.randomseed-initial_seed),'\n')
        #print("Theroritical Randomseed = %d" %(initial_seed+(self.end-self.start)*self.nR))
        return self.RandomField

    def compute_async(self, n_process, randomseed):

        pool = Pool(processes = n_process)
        self.n_process = n_process
        self.nR = self.nR * n_process
        self.RandomField = np.empty([self.size, self.nR])

        randomseed=[]
        parallel=[]
        for i in range(n_process):
            s=self.randomseed+int(i)*(self.nR+300)*(self.size)
            randomseed.append(int(s))
            parallel.append(True)

        Z = pool.starmap(self.compute,zip(randomseed,parallel))
        #print("Zshape=",np.shape(Z))
        
        for i in range(n_process):
            for j in range(int(self.nR/n_process)):
                self.RandomField[:,(j+int(i*self.nR/n_process))]=Z[i][:,j]
                
        return self.RandomField

    def variogram_compute(self, n_process = 1):
        
        pool = Pool(processes=n_process)
        model_len = self.size
        
        x=np.linspace(0,self.size-1,model_len).reshape(model_len,1)
        
        L=[]
        for i in range(self.nR):
            L.append(np.hstack([x,self.RandomField[:,i].reshape(model_len,1)]))
            
        self.Variogram = pool.starmap(self.model.Variogram, zip(L))
        self.Variogram = np.array(self.Variogram).T

    def MeanPlot(self,n,mean=0,std=1):
        
        m_plot = Visualize(self.model, self.RandomField)
        m_plot.MeanPlot(n, mean, std)
        
    def VarPlot(self,mean=0,std=1):
    
        s_plot = Visualize(self.model, self.RandomField)
        s_plot.Variance_Plot(mean, std)

    def Cdf_Plot(self):

        c_plot = Visualize(self.model, self.RandomField)
        c_plot.CDF_Plot()
    
    def Hist_Plot(self):

        h_plot = Visualize(self.model, self.RandomField)
        h_plot.HIST()

    def VarioPlot(self):

        v_plot = Visualize(self.model, self.RandomField)
        v_plot.Variogram_Plot(self.Variogram)

    def Savedata(self,path):
                
        for i in range(self.nR):
            
            if i<10:
                number='000'+str(i)
            elif 10<=i<100:
                number='00'+str(i)
            elif 100<=i<1000:
                number='0'+str(i)
            elif i>=1000:
                number=str(i)
            
            with open(path+'Realizations'+number+'.txt', 'w') as f:

                for j in range(0, self.size):

                    print('%.2d' %(j) ,'%10.6f' %(self.RandomField[j,i]), file=f)


