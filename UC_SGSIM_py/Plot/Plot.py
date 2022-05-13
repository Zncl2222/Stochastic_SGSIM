import matplotlib.pyplot as plt
import UC_SGSIM_py as UC
from UC_SGSIM_py.Plot.base import Plot_Base
from UC_SGSIM_py.Cov_Model.model import Gaussian, Spherical
import numpy as np
import time

class Visualize(Plot_Base):

    def __init__(self, model, RandomField):
        super().__init__(model, RandomField)

    def MeanPlot(self,n,mean=0,std=1):
    
        nR = len(self.RandomField[0])
        
        if n=="ALL":
            
            for i in range(nR):
                plt.figure(77879,figsize=self.figsize)
                plt.plot(self.RandomField[:,i]*std+mean)
                plt.title("Realizations: "+self.model_name,fontsize=20)
                plt.xlabel("Distance(m)",fontsize=20)
                plt.axhline(y=mean, color='r', linestyle='--',zorder=1)
                plt.ylabel("Thermal Conductivity",fontsize=20)
            
        else:
            for item in n:
                plt.figure(77879,figsize=self.figsize)
                plt.plot(self.RandomField[:,item]*std+mean)
                plt.title("Realizations: "+self.model_name,fontsize=20)
                plt.xlabel("Distance(m)",fontsize=20)
                plt.axhline(y=mean, color='r', linestyle='--',zorder=1)
                plt.ylabel("Thermal Conductivity",fontsize=20)

        
    def Statistic_Plot(self,mean=0,std=1):

        Zmean=np.zeros(len(self.RandomField[:,0]))

        for i in range(len(self.RandomField[:,0])):

            Zmean[i]=np.mean(self.RandomField[i,:]*std+mean)
        
        plt.figure(5212,figsize=self.figsize)
        plt.plot(Zmean,'-s',color='k',markeredgecolor='k',markerfacecolor='y')
        plt.xlabel("Distance(m)",fontsize=20)
        plt.ylabel("Mean",fontsize=20)
        plt.axhline(y=mean, color='r', linestyle='--',zorder=1)
        plt.xticks(fontsize=17),plt.yticks(fontsize=17)

        Zvar=np.zeros(len(self.RandomField[:,0]))

        for i in range(len(self.RandomField[:,0])):

            Zvar[i]=np.var(self.RandomField[i,:]*std)
            
        plt.figure(52712,figsize=self.figsize)
        plt.plot(Zvar,'-o',color='k',markeredgecolor='k',markerfacecolor='r')
        #plt.title("Variance",fontsize=24)
        plt.xlabel("Distance(m)",fontsize=20)
        plt.ylabel("Variance",fontsize=20)
        plt.axhline(y=std**2, color='b', linestyle='--',zorder=1)
        plt.xticks(fontsize=17),plt.yticks(fontsize=17)
    
    def Variogram_Plot(self):
        
        start_time=time.time()
        
        model_len=self.size
        
        Vario=np.zeros([len(self.hs),self.nR])
        
        x=np.linspace(0,self.size-1,model_len).reshape(model_len,1)
        
        for i in range(self.nR):
            L=np.hstack([x,self.RandomField[:,i].reshape(model_len,1)])
            
            Vario[:,i]=self.model.Variogram(L)
            plt.figure(123456,figsize=(10,6))
            plt.plot(Vario[:,i],alpha=0.1)
            plt.title("Model: "+self.model_name,fontsize=20)
            plt.xlabel("Lag(m)",fontsize=20)
            plt.ylabel("Variogram",fontsize=20)
            plt.xticks(fontsize=17),plt.yticks(fontsize=17)
            print('Progress = %.2f' % (i/self.nR*100)+'%', end='\r')
            
        plt.plot(self.model.Var_compute(self.hs),'o',markeredgecolor='k',markerfacecolor='w')
        
        Vario_mean=np.zeros(len(self.hs))
        for i in range(len(self.hs)):
            
            Vario_mean[i]=np.mean(Vario[i,:])
            
        plt.plot(Vario_mean,'--',color='blue')
        
        print('Progress = %.2f' % 100+'%\n', end='\r')
        
        end_time=time.time()
        
        print('Time = ', end_time-start_time,'s')