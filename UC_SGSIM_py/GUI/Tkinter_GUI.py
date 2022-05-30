
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
import ttkthemes as th
import tkinter.ttk as ttk
import os
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

from UC_SGSIM_py.GUI.base import GUI_base

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


class UC_SGSIM_GUI(GUI_base):

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
