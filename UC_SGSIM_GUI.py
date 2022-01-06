
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import numpy as np
import time
from multiprocessing import Pool
import multiprocessing as mp
import tkinter as tk
import ttkthemes as th
import tkinter.ttk as ttk
import threading
import os
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure


class GUI_thread():
    
    def thread(func,*args):
        
        t=threading.Thread(target=func,args=args)
        t.setDaemon(True)
        t.start()

class UC_SGSIM():
    
    def __init__(self,start,end,grid_size,model,hs,bw,nR,a,C0,randomseed,mean,std,check_boundary):
        self.start=0
        self.end=end
        self.grid_size=grid_size
        self.model=model
        self.hs=hs
        self.bw=bw
        self.nR=nR
        self.a=a
        self.C0=C0
        self.randomseed=randomseed
        self.RandomField=np.empty([self.end-self.start+1,self.nR])
        self.mean=mean
        self.std=std
        self.check_boundary=check_boundary
        self.Vario=np.zeros(10)
            
    def Variogram(self,L, hs, bw):

        dist=squareform(pdist(L[:,:1]))

        variogram=[]

        for h in hs:

            Z=[]

            for i in range(len(dist[:,0])):

                for j in range(i+1,len(dist[:,0])):

                    if( dist[i,j] >= h-bw )and( dist[i,j] <= h+bw ):

                        Z.append(np.power(L[i,1]-L[j,1],2))

            if np.sum(Z)>=1e-6:

                variogram.append(np.sum(Z)/(2*len(Z)))

        return np.array(variogram)


    def Gaussian(h,a,C0):
        
        
        x=C0*(1-np.exp(-3*h**2/a**2))
        
        return x


    def Spherical( h, a, C0 ):

        if h <= a:
            return C0*( 1.5*h/a - 0.5*(h/a)**3.0 )
        else:
            return C0
    
    def Exponential(h, a, C0):

        x = C0*(1-np.exp(-3*h/a))
        
        return x


    def Cov_model( Y, model,a,C0=1 ):
        
        Z=np.empty(len(Y))
        
        for i in range(len(Y)):

            Z[i]=C0-model(Y[i],a,C0)

        return Z

    def Var_model(self, Y, model,a,C0=1 ):

        Z=np.empty(len(Y))

        for i in range(len(Y)):

            Z[i]=model(Y[i],a,C0)

        return Z
    
    def SimpleKrige(self, Y, model, hs, bw, u, N ,randomseed,a,C0=1):
        
        
        if N==0:
            np.random.seed(randomseed)
            temp=np.random.normal(0,1,1)
            return temp
        
        dist = abs( Y[:,0]-u )
        dist = dist.reshape(len(dist),1)

        close=0
        for i in range(len(dist)):
            if dist[i]<=a:
                close=close+1

        if close==0:
            np.random.seed(randomseed)
            temp=np.random.normal(0,1,1)
            return temp
        
        Y=np.hstack([Y,dist])

        Y=sorted(Y,key=lambda x:x[2])[:N]
        
        Y=np.array(Y)

        #meanvalue = np.mean( Y[:,1] )
        meanvalue = 0
        
        C_dist = UC_SGSIM.Cov_model(Y[:,2],model,a,C0)

        C_dist = np.matrix( C_dist ).T

        C_data = squareform( pdist( Y[:,:1] ) )

        C_data = UC_SGSIM.Cov_model(C_data.flatten(),model,a)
        
        C_data = np.array( C_data ).reshape(N,N)
        
        try:
            weights = np.linalg.inv(C_data) * C_dist
        except np.linalg.LinAlgError as err:
            return 9999.9999

        residuals = (Y[:,1] - meanvalue)

        estimation=np.dot(weights.T,residuals)+meanvalue
        #estimation=np.dot(weights.T,Y[:,1])+(1-np.sum(weights))*meanvalue
        v=np.dot(weights.T,C_dist)

        var=float(1-v)
        
        if var<0:
            var=0
        
        cstd=np.sqrt(var)

        np.random.seed(randomseed)
        fix=np.random.normal(0,cstd,1)
        

        Simulated=float(estimation+fix)

        return Simulated

    def Shi_Ming(A,B):

        A=A.tolist()
        B=B.tolist()

        for i in range(len(A)):
            if np.float64(A[i][0])>np.float64(B[0]):
                A.insert(i,B)
                return A
        A.append(B)
        return np.array(A)
    
    
    def Simulation_P(start,end,grid_size,model,hs,bw,nR,a,C0,randomseed,constrain_boundary):
        
        import numpy as np
        RandomField=np.empty([int(end-start+1),nR])
        
        counts=0
        warning_count=0
        warning_count2=0

        while counts<nR:
            
            u=np.linspace(start,end,end-start+1)

            if constrain_boundary==False:
                np.random.seed(randomseed)
                first=np.random.randint(start,end,1)
                k=np.array([first]).reshape(1,1)
                t=np.random.normal(0,1,1).reshape(1,1)
                Z = np.zeros(end-start+1)
                Z[first]=t
                u=np.delete(u,[first])
                neigh=0

            elif constrain_boundary==True:
                t=np.random.normal(0,1,2).reshape(2,1)
                k=np.array([start,end]).reshape(2,1)
                Z = np.zeros(end-start+1)
                Z[start]=t[0]
                Z[end]=t[1]
                u=np.delete(u,[start,end])
                neigh=2
                
            L=np.hstack([k,t])
            
            randompath=np.random.choice(u,len(u),replace=False)
            
            
            for i in range(len(u)):

                Z[int(randompath[i])]=UC_SGSIM.SimpleKrige(0,L,model,hs,bw,randompath[i],neigh,randomseed,a,C0)

                if Z[int(randompath[i])]==9999.9999:

                    return 1
                
                temp=np.hstack([randompath[i],Z[int(randompath[i])]])

                L=np.vstack([L,temp])
                
                if neigh<6:
                    neigh=neigh+1
                    
                randomseed=randomseed+1
                
            Z_Gap=abs(Z.max()-Z.min())
                 

            if 2<Z_Gap<8:
                RandomField[:,counts]=Z
                counts=counts+1
                warning_count2=warning_count
 
            randomseed=randomseed+1
            warning_count=warning_count+1

            if abs(warning_count2-warning_count)>100:
                return 0


        return RandomField
        
    def Simulation_Parallel(self,n_thread):
        self.judge=0
        self.nCPUS=n_thread
        # initialize the map list for parallel computation
        pool = Pool(processes=n_thread)

        start=[];end=[];grid_size=[];model=[];hs=[];bw=[];nR=[];a=[];C0=[];randomseed=[];c_boundary=[]

        for i in range(n_thread):
            s=self.randomseed+int(i)*(self.nR/n_thread+300)*(self.end-self.start)
            start.append(self.start)
            end.append(self.end)
            grid_size.append(self.grid_size)
            model.append(self.model)
            hs.append(self.hs)
            bw.append(self.bw)
            nR.append(int(self.nR/n_thread))
            a.append(self.a)
            C0.append(self.C0)
            randomseed.append(int(s))
            c_boundary.append(self.check_boundary)

        Z=pool.starmap(UC_SGSIM.Simulation_P,zip(start,end,grid_size,model,hs,bw,nR,a,C0,randomseed,c_boundary))
        

        if np.max(Z)==0:
            warning=tk.messagebox.showerror(title = 'Error',message="Diverge (Please check the range scale and the model len)")
            self.judge=2
            return 0
        elif np.max(Z)==1:
            warning=tk.messagebox.showerror(title = 'Error',message="Singular Matrix (Please check the range scale and the model len)")
            self.judge=2
            return 0

        for i in range(n_thread):
            for j in range(int(self.nR/n_thread)):
                self.RandomField[:,(j+int(i*self.nR/n_thread))]=Z[i][:,j]
        
        self.judge=1
        self.Statistic_Plot(mean=self.mean,std=self.std)
        self.Plot(mean=self.mean,std=self.std)

        return self.RandomField

    def Plot(self,mean=0,std=1):
        
        fig3.clear()
        f3=fig3.add_subplot(111,title="Realizations",xlabel='Distance (m)',ylabel='Y')
        f3.plot(self.RandomField*std+mean)
        f3.grid(alpha=0.3)
        f3.axhline(y=mean, color='r', linestyle='--',zorder=1)
        canvas3.draw_idle()
        
    
    def Statistic_Plot(self,mean=0,std=1):

        Zmean=np.zeros(len(self.RandomField[:,0]))


        for i in range(len(self.RandomField[:,0])):

            Zmean[i]=np.mean(self.RandomField[i,:]*std+mean)
        
        fig.clear()
        f1=fig.add_subplot(111,title="Mean",xlabel='Distance (m)',ylabel='Mean')#.plot(Zmean,'-o',color='k',markeredgecolor='k',markerfacecolor='r')
        f1.plot(Zmean,'-o',color='k',markeredgecolor='k',markerfacecolor='r')
        f1.grid(alpha=0.3)
        f1.axhline(y=mean, color='r', linestyle='--',zorder=1)
        canvas.draw_idle()

        Zvar=np.zeros(len(self.RandomField[:,0]))

        for i in range(len(self.RandomField[:,0])):

            Zvar[i]=np.var(self.RandomField[i,:]*std)

        fig2.clear()
        f2=fig2.add_subplot(111,title="Variance",xlabel='Distance (m)',ylabel='Variance')
        f2.plot(Zvar,'-o',color='k',markeredgecolor='k',markerfacecolor='b')
        f2.grid(alpha=0.3)
        f2.axhline(y=std**2, color='b', linestyle='--',zorder=1)
        canvas2.draw_idle()
            
        
    def Variogram_Plot_P(self,n_thread):

        self.judge=0
        
        start_time=time.time()

        pool = Pool(processes=n_thread)

        if self.model==UC_SGSIM.Gaussian:
            
            Str_model="Gaussian"
            
        elif self.model==UC_SGSIM.Spherical:
            
            Str_model="Spherical"
                    
        elif self.model==UC_SGSIM.Exponential:
            
            Str_model="Exponential"
        
        model_len=self.end-self.start+1
        
        self.Vario=np.zeros([len(self.hs),self.nR])
        
        x=np.linspace(self.start,self.end,model_len).reshape(model_len,1)
        
        hs=[];bw=[];temp=[];L=[]
        
        for i in range(self.nR):
            temp.append(0)
            L.append(np.hstack([x,self.RandomField[:,i].reshape(model_len,1)]))
            hs.append(self.hs)
            bw.append(self.bw)
            
        Z=pool.starmap(UC_SGSIM.Variogram,zip(temp,L,hs,bw))
        self.Vario=np.array(Z).T
        fig4.clear()
        f4=fig4.add_subplot(111,title="Normalized Variogram",xlabel='Distance (m)',ylabel='Variogram')
        f4.plot(self.Vario,alpha=0.2)
        f4.plot(self.Var_model(self.hs,self.model,self.a),'o',markeredgecolor='k',markerfacecolor='w',label='Theoritical variogram')
        f4.grid(alpha=0.3)
        f4.axhline(y=1, color='r', linestyle='--',zorder=1)
        
        
        Vario_mean=np.zeros(len(self.hs))
        
        for i in range(len(self.hs)):
            
            Vario_mean[i]=np.mean(self.Vario[i,:])
            
        f4.plot(Vario_mean,'--',color='blue',label='Mean variogram')
        f4.legend()
        canvas4.draw_idle()

        self.judge=1

        end_time=time.time()
        print('Time = ', end_time-start_time,'s')

    def Savedata(self,path):

        try:
            os.mkdir(path+"\\UC_SGSIM_Realizations")
        except:
            print("Folder has already existed")

        if np.max(self.Vario)!=0:
            try:
                os.mkdir(path+"\\UC_SGSIM_Realizations\\Variogram")
            except:
                print("Folder has already existed")


        time.sleep(1)
        
        path=path+"\\UC_SGSIM_Realizations\\"
        
        mlen=self.end-self.start+1
   
        with open(path+"ParameterSettings.txt",'w') as ff:

            print(mlen," #Model len", file=ff)
            print(self.a," #Effective range",file=ff)
            if self.model==UC_SGSIM.Gaussian:
                print("Gaussian #Variogram model",file=ff)
            elif self.model==UC_SGSIM.Spherical:
                print("Spherical #Variogram model",file=ff)
            elif self.model==UC_SGSIM.Exponential:
                print("Exponential #Variogram model",file=ff)
            print(self.nR," #Number of Realizations", file=ff)
            print(self.randomseed," #Number of initial seed",file=ff)
            print(self.nCPUS," #Number of CPUs",file=ff)
            print(self.mean," #Mean",file=ff)
            print(self.std," #Standard deviation",file=ff)
            print("Python #Programming language",file=ff)
            print(check_logarithm.get()," #Logrithm data",file=ff)
            if np.max(self.Vario)!=0:
                print("y #Variogram calculation",file=ff)
            else:
                print("n #Variogram calculation",file=ff)
            
        
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

                for j in range(0, mlen):
                    if check_logarithm.get()==False:
                        print('%.2d' %(j) ,'%.12f' %(self.RandomField[j,i]*self.std+self.mean), file=f)
                    else:
                        print('%.2d' %(j) ,'%.12f' %(np.exp(self.RandomField[j,i]*self.std+self.mean)), file=f)

            if np.max(self.Vario)!=0:
                with open(path+'\Variogram\\Variogram'+number+'.txt', 'w') as f:
                    for j in range(0, len(self.hs),self.bw):
                            print('%.2d' %(j) ,'%.12f' %(self.Vario[j,i]), file=f)



    def Progress_start(self):
        import time

        Data_Label["background"]='red'
        
        while self.judge==0:
            Data_Label['text']='Processing.'
            time.sleep(0.2)
            Data_Label['text']='Processing..'
            time.sleep(0.2)
            Data_Label['text']='Processing...'
            time.sleep(0.2)
            Data_Label['text']='Processing....'
            time.sleep(0.2)
            Data_Label['text']='Processing.....'
        if self.judge==2:
            Data_Label["text"]="Error"
        elif self.judge==1:
            Data_Label["text"]="Done"
            Data_Label["background"]='orange'

class Validation():

    def __init__(self,path,n_thread,hs,bw):
        self.path=path
        self.n_thread=n_thread
        self.hs=hs
        self.bw=bw

    def Validate(self):

        self.judge=0
        import pandas as pd
        import os
        tempa=pd.read_table(self.path+"\\ParameterSettings.txt",header=None,sep='#')

        self.mlen=int(tempa.iloc[0,0])
        self.a=float(tempa.iloc[1,0])
    
        model=tempa.iloc[2,0]
 
        if model=="Gaussian ":
            self.model=UC_SGSIM.Gaussian
        elif model=="Spherical ":
            self.model=UC_SGSIM.Spherical
        elif model=="Expenoential ":
            self.model=UC_SGSIM.Exponential
        self.nR=int(tempa.iloc[3,0])
        #self.n_thread=int(tempa.iloc[5,0])
        self.mean=float(tempa.iloc[6,0])
        self.std=float(tempa.iloc[7,0])
        self.logarithm=tempa.iloc[9,0]
        self.vario_check=tempa.iloc[10,0]
          
        self.RandomField=np.zeros([self.mlen,self.nR])

        for i in range(self.nR):

            if i<10:
                number='000'+str(i)
            elif 10<=i<100:
                number='00'+str(i)
            elif 100<=i<1000:
                number='0'+str(i)
            elif i>=1000:
                number=str(i)
            if tempa.iloc[8,0]=="Python ":
                Z=pd.read_table(self.path+"\\Realizations"+number+".txt",header=None,sep=' ')
            elif tempa.iloc[8,0]=="C ":
                Z=pd.read_table(self.path+"\\Realizations"+number+".txt",header=None,sep='\t')
            
            self.RandomField[:,i]=Z.iloc[:,1]
            
        if (self.logarithm).startswith("T")==True:
            self.RandomField=np.log(self.RandomField)
        
        fig3.clear()
        f3=fig3.add_subplot(111,title="Realizations(validation)",xlabel='Distance (m)',ylabel='Y')
        f3.plot(self.RandomField)
        f3.grid(alpha=0.3)
        f3.axhline(y=self.mean, color='r', linestyle='--',zorder=1)
        canvas3.draw_idle()
        
        #--------------------------------------Statistic Plot-----------------------------------------

        Zmean=np.zeros(len(self.RandomField[:,0]))

        for i in range(len(self.RandomField[:,0])):
            
            Zmean[i]=np.mean(self.RandomField[i,:])
        
        fig.clear()
        f1=fig.add_subplot(111,title="Mean(validation)",xlabel='Distance (m)',ylabel='Mean')#.plot(Zmean,'-o',color='k',markeredgecolor='k',markerfacecolor='r')
        f1.plot(Zmean,'-o',color='k',markeredgecolor='k',markerfacecolor='r')
        f1.grid(alpha=0.3)
        f1.axhline(y=self.mean, color='r', linestyle='--',zorder=1)
        canvas.draw_idle()

        Zvar=np.zeros(len(self.RandomField[:,0]))

        for i in range(len(self.RandomField[:,0])):

            Zvar[i]=np.var(self.RandomField[i,:])

        fig2.clear()
        f2=fig2.add_subplot(111,title="Variance(validation)",xlabel='Distance (m)',ylabel='Variance')
        f2.plot(Zvar,'-o',color='k',markeredgecolor='k',markerfacecolor='b')
        f2.grid(alpha=0.3)
        f2.axhline(y=self.std**2, color='b', linestyle='--',zorder=1)
        canvas2.draw_idle()

        #----------------------------------Variogram Plot--------------------------------------------
        if self.vario_check=="y ":

            a=os.listdir(self.path+"\\Variogram\\")
        
            xf=pd.read_table(self.path+"\\Variogram\\"+a[0],header=None,encoding='gb2312',sep='\s+')
            
            self.vario_len=len(xf)
            
            self.Variogram=np.empty([self.vario_len,self.nR])
            
            for i in range(self.nR):
                x = pd.read_table(self.path+"\Variogram\\"+a[i],header=None,encoding='gb2312',sep='\s+')
                self.Variogram[:,i]=x.iloc[:,1]
                
        else:
            pool = Pool(processes=self.n_thread)
            
            self.Variogram=np.zeros([len(self.hs),self.nR])
            
            x=np.linspace(0,self.mlen-1,self.mlen).reshape(self.mlen,1)
            
            hs=[];bw=[];temp=[];L=[]
            
            for i in range(self.nR):
                temp.append(0)
                L.append(np.hstack([x,((self.RandomField[:,i]-self.mean)/self.std).reshape(self.mlen,1)]))
                hs.append(self.hs)
                bw.append(self.bw)
                
            Z=pool.starmap(UC_SGSIM.Variogram,zip(temp,L,hs,bw))
            self.Variogram=np.array(Z).T

        fig4.clear()
        f4=fig4.add_subplot(111,title="Normalized Variogram(validation)",xlabel='Distance (m)',ylabel='Variogram')
        f4.plot(self.Variogram,alpha=0.2)
        f4.plot(UC_SGSIM.Var_model(0,self.hs,self.model,self.a),'o',markeredgecolor='k',markerfacecolor='w',label='Theoritical variogram')
        f4.grid(alpha=0.3)
        f4.axhline(y=1, color='r', linestyle='--',zorder=1)
        
        
        Vario_mean=np.zeros(len(self.hs))
        
        for i in range(len(self.hs)):
            
            Vario_mean[i]=np.mean(self.Variogram[i,:])
            
        f4.plot(Vario_mean,'--',color='blue',label='Mean variogram')
        f4.legend()
        canvas4.draw_idle()
        self.judge=1

    def Progress_start(self):
        import time

        Data_Label["background"]='red'
        
        while self.judge==0:
            Data_Label['text']='Processing.'
            time.sleep(0.2)
            Data_Label['text']='Processing..'
            time.sleep(0.2)
            Data_Label['text']='Processing...'
            time.sleep(0.2)
            Data_Label['text']='Processing....'
            time.sleep(0.2)
            Data_Label['text']='Processing.....'
        if self.judge==2:
            Data_Label["text"]="Error"
        elif self.judge==1:
            Data_Label["text"]="Done"
            Data_Label["background"]='orange'

class Toolbar(NavigationToolbar2Tk):

    def set_message(self, s):
        pass

class Scroll_Bar():

    def __init__(self, container ,root ,*args, **kwargs):

        def _on_mousewheel(event):
            self.scroll_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
                
        self.scroll_canvas=tk.Canvas(container)
        self.scrollbar = ttk.Scrollbar(container, orient="vertical", command=self.scroll_canvas.yview)
        self.scrollbar_h = ttk.Scrollbar(container,orient="horizontal",command=self.scroll_canvas.xview)
  
        self.tab=ttk.Frame(self.scroll_canvas)

        self.tab.bind(
            "<Configure>",
            lambda e: self.scroll_canvas.configure(
                scrollregion=self.scroll_canvas.bbox("all")
            )
        )

        self.scroll_canvas.bind_all("<MouseWheel>", _on_mousewheel)
        self.scroll_canvas.create_window((0, 0), window=self.tab, anchor="nw")
        self.scroll_canvas.configure(yscrollcommand=self.scrollbar.set,xscrollcommand=self.scrollbar_h.set)

        container.grid(row=0,column=0,ipadx=int(root.winfo_screenwidth()/2.5),ipady=int(root.winfo_screenheight()/2.5))
        self.scroll_canvas.grid(row=0,column=0,ipadx=int(root.winfo_screenwidth()/2.5),ipady=int(root.winfo_screenheight()/2.5))

        self.scrollbar.grid(row=0, column=1,  sticky=tk.NS)
        self.scrollbar_h.grid(row=1, column=0, sticky=tk.EW)


if __name__=='__main__':


    from multiprocessing import freeze_support
    freeze_support()


    temp=mp.cpu_count()
    clist=[]
    for i in range(temp):
        clist.append(i+1)

    root = tk.Tk()
    root.wm_title('1D Unconditional Sequential Gaussian Simulation')
    padx= 0
    pady= 0
    root.geometry(f"+{padx}+{pady}")
    #root.maxsize(width=1400,height=950)
    #root.minsize(width=300,height=200)
    #root.geometry(str(int(root.winfo_screenwidth()/1.3))+"x"+str(int(root.winfo_screenheight()/1.3)))
    #root.geometry('1600x900')
    #root.resizable(width=0,height=0)

    style = th.ThemedStyle(root)
    style.set_theme("black")
        
    def drive(pjudge):
        global C1
        if pjudge==0:
            
            end=int(modellen_entry.get())
            bw=int(Lagsteps_entry.get())
            hs=np.arange(0,int(Laglen_entry.get()),bw)
            nR=int(nR_entry.get())
            a=float(a_entry.get())
            C0=1
            randomseed=int(init_randomseed_entry.get())
            n_thread=int(Select.get())
            mean=float(mean_entry.get())
            std=float(std_entry.get())
            check=check_boundary.get()
            warning=0

            if nR%n_thread !=0:
                warning=tk.messagebox.askquestion(title = 'Warning',message="Realizations number can't divide by the CPUs number,"+ 
                                                                            "it might occur errors. \nDo you wish to proceed? ")
            
            if warning=='no':
                return 0

            if model_Select.get()=="Gaussian":
                model=UC_SGSIM.Gaussian

            elif model_Select.get()=="Spherical":
                model=UC_SGSIM.Spherical

            elif model_Select.get()=="Exponential":
                model=UC_SGSIM.Exponential
           

            C1=UC_SGSIM(0,end,1,model,hs,bw,nR,a,C0,randomseed,mean,std,check)

            GUI_thread.thread(C1.Progress_start)
                
            GUI_thread.thread(C1.Simulation_Parallel,n_thread)

        elif pjudge==1:
        
            if int(nR_entry.get())%int(Select.get()) !=0:
                tk.messagebox.showerror(title = 'Error',message="Realization number should be divided evenly by the CPUs number")
                return 0
            try:
                GUI_thread.thread(C1.Progress_start)
                GUI_thread.thread(C1.Variogram_Plot_P,int(Select.get()))
            except NameError as C1:
                tk.messagebox.showerror(title = 'Error',message="Please run the simulation first")

    def save():
        C1.Savedata(savedir_entry.get())  
        #try:
            #C1.Savedata(dir_entry.get())
        #except NameError as C1:
            #tk.messagebox.showerror(title = 'Error',message="Please run the simulation first")

    def validate():

        global C2
        bw=int(Lagsteps_entry.get())
        hs=np.arange(0,int(Laglen_entry.get()),bw)
        C2=Validation(validate_entry.get(),int(Select.get()),hs,bw)
        GUI_thread.thread(C2.Progress_start)
        GUI_thread.thread(C2.Validate)

    tab_main=ttk.Notebook()
    tab_main.grid_columnconfigure(0, weight=1)
    tab_main.grid_rowconfigure(0, weight=1)

    St=Scroll_Bar(tab_main,root)

    tab1=St.tab
    tab_main.add(St.scroll_canvas,text='Stochastic')

    tk.Grid.columnconfigure(tab_main,0,weight=1)
    tk.Grid.rowconfigure(tab_main,0,weight=1)
    #---------------------------------------- Browse function--------------------------------------
    def Save_browse():
        filename = tk.filedialog.askdirectory()
        savedir_entry.delete(0,tk.END)
        savedir_entry.insert(0,filename)
    def Validate_browse():
        filename = tk.filedialog.askdirectory()
        validate_entry.delete(0,tk.END)
        validate_entry.insert(0,filename)
    #---------------------------------------- Simulation GUI----------------------------------------

    Plot1=ttk.Notebook(tab1)
    #Plot1.place(x=550,y=25,width=500,height=400)
    Plot1.grid(row=0,column=3,rowspan=8,padx=30,ipady=100,ipadx=150)

    Plot2=ttk.Notebook(tab1)
    #Plot2.place(x=550,y=450,width=500,height=400)
    Plot2.grid(row=9,column=3,rowspan=8,padx=30,ipady=100,ipadx=150)

    Plot3=ttk.Notebook(tab1)
    #Plot3.place(x=1100,y=25,width=500,height=400)
    Plot3.grid(row=0,column=4,rowspan=8,padx=30,ipady=100,ipadx=150)

    Plot4=ttk.Notebook(tab1)
    #Plot4.place(x=1100,y=450,width=500,height=400)
    Plot4.grid(row=9,column=4,rowspan=8,padx=30,ipady=100,ipadx=150)

    global fig,fig2,fig3,fig4,canvas,canvas2,canvas3,canvas4

    fig = Figure(figsize=(2, 2), dpi=100,tight_layout=True)
    fig2= Figure(figsize=(2, 2), dpi=100,tight_layout=True)
    fig3 = Figure(figsize=(2, 2), dpi=100,tight_layout=True)
    fig4= Figure(figsize=(2, 2), dpi=100,tight_layout=True)

    canvas = FigureCanvasTkAgg(fig, master=Plot1)  
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    toolbar = Toolbar(canvas, Plot1)

    canvas2 = FigureCanvasTkAgg(fig2, master=Plot2)  
    canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    toolbar = Toolbar(canvas2, Plot2)

    canvas3 = FigureCanvasTkAgg(fig3, master=Plot3)  
    canvas3.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    toolbar = Toolbar(canvas3, Plot3)

    canvas4 = FigureCanvasTkAgg(fig4, master=Plot4)  
    canvas4.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    toolbar = Toolbar(canvas4, Plot4)

    fig.add_subplot(111,title="Mean",xlabel='Distance',ylabel='Mean')
    canvas.draw_idle()        
   
    fig2.add_subplot(111,title="Variance",xlabel='Distance',ylabel='Variance')
    canvas2.draw_idle()

    fig3.add_subplot(111,title="Realizations",xlabel='Distance',ylabel='Y')
    canvas3.draw_idle()
    
    fig4.add_subplot(111,title="Variogram",xlabel='Lag dist',ylabel='Variogram')
    canvas4.draw_idle()

    modellen_entry=ttk.Entry(tab1)
    modellen_entry.insert(0,"150")
    #modellen_entry.place(x=250,y=50,width=200)
    modellen_entry.grid(row=1,column=1,ipadx=100,padx=15)

    nR_entry=ttk.Entry(tab1)
    nR_entry.insert(0,"10")
    #nR_entry.place(x=250,y=100,width=200)
    nR_entry.grid(row=2,column=1,ipadx=100,padx=15)

    a_entry=ttk.Entry(tab1)
    a_entry.insert(0,"17.32")
    #a_entry.place(x=250,y=150,width=200)
    a_entry.grid(row=3,column=1,ipadx=100,padx=15)

    init_randomseed_entry=ttk.Entry(tab1)
    init_randomseed_entry.insert(0,"1321")
    #init_randomseed_entry.place(x=250,y=200,width=200)
    init_randomseed_entry.grid(row=4,column=1,ipadx=100,padx=15)

    Laglen_entry=ttk.Entry(tab1)
    Laglen_entry.insert(0,"35")
    #Laglen_entry.place(x=250,y=250,width=90)
    Laglen_entry.grid(row=5,column=1,sticky=tk.W,padx=15,ipadx=10)
    
    Lagsteps_entry=ttk.Entry(tab1)
    Lagsteps_entry.insert(0,"1")
    #Lagsteps_entry.place(x=360,y=250,width=90)
    Lagsteps_entry.grid(row=5,column=1,sticky=tk.E,padx=15,ipadx=10)

    savedir_entry=ttk.Entry(tab1)
    #savedir_entry.place(x=250,y=400,width=200)
    savedir_entry.grid(row=8,column=1,ipadx=100,padx=15)

    mean_entry=ttk.Entry(tab1)
    mean_entry.insert(0,"0")
    #mean_entry.place(x=250,y=450,width=200)
    mean_entry.grid(row=9,column=1,ipadx=100,padx=15)

    std_entry=ttk.Entry(tab1)
    std_entry.insert(0,"1")
    #std_entry.place(x=250,y=500,width=200)
    std_entry.grid(row=10,column=1,ipadx=100,padx=15)

    validate_entry=ttk.Entry(tab1)
    #validate_entry.place(x=250,y=550,width=200)
    validate_entry.grid(row=11,column=1,ipadx=100,padx=15)

    option1_Label=ttk.Label(tab1,text="Simulation Options", foreground='white',font='24')
    #option1_Label.place(x=200,y=15)
    option1_Label.grid(row=0,column=1)

    modellen_Label=ttk.Label(tab1,text="Model len")
    #modellen_Label.place(x=50,y=50)
    modellen_Label.grid(row=1,column=0,padx=10)

    nR_Label=ttk.Label(tab1,text="Number of Realization")
    #nR_Label.place(x=50,y=100)
    nR_Label.grid(row=2,column=0,padx=10)

    a_Label=ttk.Label(tab1,text="Effective range")
    #a_Label.place(x=50,y=150)
    a_Label.grid(row=3,column=0,padx=10)

    seed_Label=ttk.Label(tab1,text="Initial seed")
    #seed_Label.place(x=50,y=200)
    seed_Label.grid(row=4,column=0)

    Lag_Label=ttk.Label(tab1,text="Lag len/steps")
    #Lag_Label.place(x=50,y=250)
    Lag_Label.grid(row=5,column=0)

    model_label=ttk.Label(tab1,text="Variogram model")
    #model_label.place(x=50,y=300)
    model_label.grid(row=6,column=0)

    option2_Label=ttk.Label(tab1,text="  Data Options  ", foreground='white',font='24')
    #option2_Label.place(x=200,y=350)
    option2_Label.grid(row=7,column=1)

    dir_Label=ttk.Label(tab1,text="Path for saving data")
    #dir_Label.place(x=50,y=400)
    dir_Label.grid(row=8,column=0)

    Mean_Label=ttk.Label(tab1,text="Mean Value")
    #Mean_Label.place(x=50,y=450)
    Mean_Label.grid(row=9,column=0)

    std_Label=ttk.Label(tab1,text="Standard deviation")
    #std_Label.place(x=50,y=500)
    std_Label.grid(row=10,column=0)

    validate_Label=ttk.Label(tab1,text="Data path (Validation)")
    #validate_Label.place(x=50,y=550)
    validate_Label.grid(row=11,column=0)

    ncore_Label=ttk.Label(tab1,text="CPUs number")
    #ncore_Label.place(x=50,y=600)
    ncore_Label.grid(row=12,column=0)

    Data_Label=ttk.Label(tab1,text="")
    #Data_Label.place(x=50,y=675)
    Data_Label.grid(row=13,column=0)

    #------------------------------------------Button---------------------------------

    run_button = ttk.Button(tab1, text='Run', command=lambda:drive(0))
    #run_button.place(x=150,y=670,height=30,width=130)
    run_button.grid(row=13,column=1,pady=10)

    var_button = ttk.Button(tab1, text='Variogram plot', command=lambda:drive(1))
    #var_button.place(x=300,y=670,height=30,width=130)
    var_button.grid(row=14,column=1,pady=10)

    save_button = ttk.Button(tab1, text='Save', command=save)
    #save_button.place(x=225,y=720,height=30,width=130)
    save_button.grid(row=15,column=1,pady=10)

    validate_button = ttk.Button(tab1, text='Validate', command=validate)
    #validate_button.place(x=225,y=770,height=30,width=130)
    validate_button.grid(row=16,column=1,pady=10)

    number = tk.StringVar()
    model_Select= ttk.Combobox(tab1, width=12, textvariable=number, state='readonly')
    model_Select['values'] = ["Gaussian","Spherical","Exponential"]
    #model_Select.place(x=250,y=300)
    model_Select.grid(row=6,column=1,sticky=tk.W,padx=15)      
    model_Select.current(0)

    number2 = tk.StringVar()
    Select= ttk.Combobox(tab1, width=12, textvariable=number2, state='readonly')
    Select['values'] = clist
    #Select.place(x=250,y=600)
    Select.grid(row=12,column=1,sticky=tk.W,padx=15,pady=10)
    Select.current(0)

    Save_browse_button=ttk.Button(tab1,text="Browse",command=Save_browse)
    #Save_browse_button.place(x=460,y=398)
    Save_browse_button.grid(row=8,column=2)

    Validate_browse_button=ttk.Button(tab1,text="Browse",command=Validate_browse)
    #Validate_browse_button.place(x=460,y=548)
    Validate_browse_button.grid(row=11,column=2)

    check_boundary=tk.BooleanVar()
    check_boundary.set(True)

    checkbox=ttk.Checkbutton(tab1,text="Constrain boundary",var=check_boundary)
    #checkbox.place(x=50,y=335)
    checkbox.grid(row=6,column=1,sticky=tk.E,padx=50)

    check_logarithm=tk.BooleanVar()
    check_logarithm.set(False)

    checkbox_logarithm=ttk.Checkbutton(tab1,text="Back transformation (natural log)",var=check_logarithm)
    #checkbox_logarithm.place(x=50,y=635)
    checkbox_logarithm.grid(row=12,column=1,sticky=tk.E,padx=10)

    tk.mainloop()
