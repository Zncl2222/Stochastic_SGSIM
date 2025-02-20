import time

import matplotlib.pyplot as plt
import uc_sgsim as uc
from uc_sgsim.cov_model import Gaussian

if __name__ == '__main__':
    start = time.time()
    x = 151  # Model grid, only 1D case is support now

    bw_s = 1  # lag step
    bw_l = 35  # lag range
    randomseed = 151  # randomseed for simulation
    k_range = 17.32  # effective range of covariance model
    sill = 1  # sill of covariance model

    nR = 2  # numbers of realizations in each CPU cores,
    # if nR = 1 n_process = 8
    # than you will compute total 8 realizations

    # Create Covariance model first
    cov_model = Gaussian(bw_l, bw_s, k_range, sill)

    # Create simulation and input the Cov model
    # You could also set min_value, max_value and max_neighbor for sgsim by key words
    # sgsim = uc.UCSgsimDLL(x, nR, cov_model, min_value=-6, max_value=6, max_neigh=10)
    # set min_value, max_value and max_neighbor by directly assign
    # sgsim.min_value = -6
    # sgsim.max_value = 6
    # sgsim.max_neigh = 10

    # Create simulation with default min_value, max_value and max_neigh params
    # sgsim = uc.UCSgsim(x, nR, cov_model)
    sgsim_c = uc.UCSgsim(x, nR, cov_model, engine='python')

    # Start compute with n CPUs
    sgsim_c.run(n_processes=2)

    mid = time.time()

    sgsim_c.plot()  # Plot realizations
    sgsim_c.mean_plot()  # Plot mean
    sgsim_c.variance_plot()  # Plot variance
    sgsim_c.cdf_plot(x_location=10)  # CDF
    sgsim_c.hist_plot(x_location=10)  # Hist
    sgsim_c.get_variogram(n_processes=2)  # Compute variogram before plotting
    # Plot variogram and mean variogram for validation
    sgsim_c.variogram_plot()
    # Save random_field and variogram
    sgsim_c.save_random_field('randomfields.csv', save_single=True)
    sgsim_c.save_variogram('variograms.csv', save_single=True)
    end = time.time()
    print('SGSIM time =', mid - start)
    print('Plot and variogram time =', end - mid)
    print('total time =', end - start)

    # show figure
    plt.show()
