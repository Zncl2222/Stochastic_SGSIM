import time

import matplotlib.pyplot as plt
import numpy as np
import uc_sgsim as UC
from uc_sgsim.Cov_Model import Gaussian

if __name__ == '__main__':
    start = time.time()
    X = 151  # Model grid, only 1D case is support now

    bw = 1  # lag step
    hs = np.arange(0.0, 35, bw)  # lag range
    randomseed = 12321  # randomseed for simulation
    a = 17.32  # effective range of covariance model
    C0 = 1  # sill of covariance model

    nR = 10  # numbers of realizations in each CPU cores,
    # if nR = 1 n_process = 8
    # than you will compute total 8 realizations

    # Create Covariance model first
    Cov_model = Gaussian(hs, bw, a, C0)

    # Create simulation and input the Cov model
    # sgsim = UC.Simulation(X, Cov_model, nR)
    sgsim = UC.Simulation_byC(X, Cov_model, nR)

    # Start compute with n CPUs
    # sgsim.compute_async(n_process=1, randomseed=454)
    sgsim.compute_by_dll(n_process=1, randomseed=151)

    mid = time.time()

    # Save data wiht .txt file
    # sgsim.Savedata(r"The path you want to save data")

    sgsim.mean_plot('ALL')  # Plot mean
    sgsim.variance_plot()  # Plot variance
    sgsim.cdf_plot(x_location=10)  # CDF
    sgsim.hist_plot(x_location=10)  # Hist
    sgsim.vario_compute_by_dll(n_process=1)
    # sgsim.variogram_compute(n_process=1)  # Compute variogram before plotting
    # Plot variogram and mean variogram for validation
    sgsim.vario_plot()
    # Save random_field and variogram
    sgsim.save_random_field('randomfiel.csv', save_single=True)
    sgsim.save_variogram('')
    end = time.time()
    print('SGSIM time =', mid - start)
    print('Plot and variogram time =', end - mid)
    print('total time =', end - start)
    print(np.shape(sgsim.random_field))

    # plt.show() to show the matplotlib plot
    plt.show()
