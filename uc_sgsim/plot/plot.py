import time

import matplotlib.pyplot as plt
import numpy as np
from uc_sgsim.plot.base import PlotBase


class Visualize(PlotBase):
    def __init__(self, model, random_field):
        super().__init__(model, random_field)

    def mean_plot(self, n, mean=0, std=1):
        realization_number = len(self.random_field[:, 0])

        if n == 'ALL':
            for i in range(realization_number):
                plt.figure(77879, figsize=self.figsize)
                plt.plot(self.random_field[i, :] * std + mean)
                plt.title('Realizations: ' + self.model_name, fontsize=20)
                plt.xlabel('Distance(-)', fontsize=20)
                plt.axhline(y=mean, color='r', linestyle='--', zorder=1)
                plt.ylabel('Y', fontsize=20)

        else:
            for item in n:
                plt.figure(77879, figsize=self.figsize)
                plt.plot(self.random_field[item, i] * std + mean)
                plt.title('Realizations: ' + self.model_name, fontsize=20)
                plt.xlabel('Distance(-)', fontsize=20)
                plt.axhline(y=mean, color='r', linestyle='--', zorder=1)
                plt.ylabel('Y', fontsize=20)

    def variance_plot(self, mean=0, std=1):
        Zmean = np.zeros(len(self.random_field[0, :]))

        for i in range(len(self.random_field[0, :])):
            Zmean[i] = np.mean(self.random_field[:, i] * std + mean)

        plt.figure(5212, figsize=self.figsize)
        plt.plot(
            Zmean,
            '-s',
            color='k',
            markeredgecolor='k',
            markerfacecolor='y',
        )
        plt.xlabel('Distance(-)', fontsize=20)
        plt.ylabel('Mean', fontsize=20)
        plt.axhline(y=mean, color='r', linestyle='--', zorder=1)
        plt.xticks(fontsize=17), plt.yticks(fontsize=17)

        Zvar = np.zeros(len(self.random_field[0, :]))

        for i in range(len(self.random_field[0, :])):
            Zvar[i] = np.var(self.random_field[:, i] * std)

        plt.figure(52712, figsize=self.figsize)
        plt.plot(
            Zvar,
            '-o',
            color='k',
            markeredgecolor='k',
            markerfacecolor='r',
        )
        plt.xlabel('Distance(-)', fontsize=20)
        plt.ylabel('Variance', fontsize=20)
        plt.axhline(y=std**2, color='b', linestyle='--', zorder=1)
        plt.xticks(fontsize=17), plt.yticks(fontsize=17)

    def cdf_plot(self, x_location):

        X = self.random_field[:, x_location]

        mu = np.mean(X)
        sigma = np.std(X)
        n_bins = 50

        _, ax = plt.subplots(figsize=(8, 4))

        _, bins, _ = ax.hist(
            X,
            n_bins,
            density=True,
            histtype='step',
            cumulative=True,
            label='Empirical',
        )

        y = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(
            -0.5 * (1 / sigma * (bins - mu)) ** 2,
        )
        y = y.cumsum()
        y /= y[-1]

        ax.plot(bins, y, 'k--', linewidth=1.5, label='Theoretical')

        ax.grid(True)
        ax.legend(loc='right')
        ax.set_title('Cumulative step histograms, x = ' + str(x_location))
        ax.set_xlabel('Random Variable (mm)')
        ax.set_ylabel('Occurrence')

    def hist_plot(self, x_location):

        X = self.random_field[:, x_location]

        mu = np.mean(X)
        sigma = np.std(X)

        num_bins = 50
        plt.figure(num=1151)
        _, bins, _ = plt.hist(
            X,
            num_bins,
            density=1,
            color='blue',
            alpha=0.5,
            edgecolor='k',
        )

        y = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(
            -0.5 * (1 / sigma * (bins - mu)) ** 2,
        )

        plt.plot(bins, y, '--', color='black')

        plt.xlabel('X-Axis')
        plt.ylabel('Y-Axis')

        plt.title('Histogram, x = ' + str(x_location))

    def variogram_plot(self, Variogram):
        start_time = time.time()
        for i in range(self.realization_number):
            plt.figure(123456, figsize=(10, 6))
            plt.plot(Variogram[i, :], alpha=0.1)
            plt.title('Model: ' + self.model_name, fontsize=20)
            plt.xlabel('Lag(m)', fontsize=20)
            plt.ylabel('Variogram', fontsize=20)
            plt.xticks(fontsize=17), plt.yticks(fontsize=17)
            print('Progress = %.2f' % (i / self.realization_number * 100) + '%', end='\r')

        plt.plot(
            self.model.var_compute(self.bandwidth_step),
            'o',
            markeredgecolor='k',
            markerfacecolor='w',
        )

        Vario_mean = np.zeros(len(self.bandwidth_step))
        for i in range(len(self.bandwidth_step)):
            Vario_mean[i] = np.mean(Variogram[:, i])

        plt.plot(Vario_mean, '--', color='blue')

        print('Progress = %.2f' % 100 + '%\n', end='\r')

        end_time = time.time()

        print('Time = ', end_time - start_time, 's')
