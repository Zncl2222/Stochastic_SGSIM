import time

import matplotlib.pyplot as plt
import uc_sgsim as uc
from uc_sgsim.cov_model import Gaussian

if __name__ == '__main__':
    start = time.time()
    x = 151  # Model grid, only 1D case is support now

    bw_s = 1  # lag step
    bw_l = 35  # lag range
    randomseed = 11251  # randomseed for simulation
    k_range = 17.32  # effective range of covariance model
    sill = 1  # sill of covariance model

    nR = 100  # numbers of realizations in each CPU cores,
    # if nR = 1 n_process = 8
    # than you will compute total 8 realizations

    # Create Covariance model first
    cov_model = Gaussian(bw_l, bw_s, k_range, sill)

    # Create simulation and input the Cov model
    sgsim = uc.UCSgsim(x, nR, cov_model)
    # sgsim = uc.UCSgsimDLL(x, nR, cov_model)

    # Start compute with n CPUs
    sgsim.compute(n_process=2, randomseed=randomseed)

    mid = time.time()

    sgsim.plot()  # Plot realizations
    sgsim.mean_plot()  # Plot mean
    sgsim.variance_plot()  # Plot variance
    sgsim.cdf_plot(x_location=10)  # CDF
    sgsim.hist_plot(x_location=10)  # Hist
    sgsim.variogram_compute(n_process=2)  # Compute variogram before plotting
    # Plot variogram and mean variogram for validation
    sgsim.variogram_plot()
    # Save random_field and variogram
    sgsim.save_random_field('randomfields.csv', save_single=True)
    sgsim.save_variogram('variograms.csv', save_single=True)
    end = time.time()
    print('SGSIM time =', mid - start)
    print('Plot and variogram time =', end - mid)
    print('total time =', end - start)

    # show figure
    plt.show()
