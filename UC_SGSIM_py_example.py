import UC_SGSIM_py as UC
import numpy as np
import matplotlib.pyplot as plt
import time

if __name__ == '__main__':
    start=time.time()
    X=range(151)            # Model grid, only 1D case is support now

    bw=1                    # lag step
    hs=np.arange(0.,35,bw)  # lag range
    randomseed=12321        # randomseed for simulation
    a=17.32                 # effective range of covariance model
    C0=1                    # sill of covariance model

    nR=100                  # numbers of realizations in each CPU cores, 
                            # if nR = 1 n_process = 8 
                            # than you will compute total 8 realizations

    Cov_model = UC.Gaussian(hs, bw, randomseed, a, C0)  # Create Covariance model first
    sgsim = UC.Simulation(X, Cov_model, nR, randomseed) # Create simulation and input the Cov model
    
    sgsim.compute_async(n_process=8, randomseed=454)    # Start compute with n CPUs
    mid=time.time()
    
    #sgsim.Savedata(r"The path you want to save data")  # Save data wiht .txt file

    sgsim.MeanPlot("ALL")                 # Plot mean
    sgsim.VarPlot()                       # Plot variance 
    sgsim.Cdf_Plot()                      # CDF
    sgsim.Hist_Plot()                     # Hist
    sgsim.variogram_compute(n_process=8)  # Compute variogram before plotting
    sgsim.VarioPlot()                     # Plot variogram and mean variogram for validation
    
    end=time.time()
    print("SGSIM time =", mid-start)
    print("Plot and variogram time =", end-mid)
    print("total time =", end-start)
    print(np.shape(sgsim.RandomField))

    plt.show()                            # plt.show() to show the matplotlib plot

