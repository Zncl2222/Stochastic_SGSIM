import matplotlib.pyplot as plt
import numpy as np
from uc_sgsim.exception import VariogramDoesNotCompute
from uc_sgsim.plot.base import PlotBase
from ..cov_model.base import CovModel


class SgsimPlot(PlotBase):
    xlabel = 'Distance(-)'
    curr_fig_num = 0

    def __init__(
        self,
        model: CovModel,
        figsize: tuple = (10, 8),
    ):
        super().__init__(figsize)
        self.__vmodel = model
        self.__vmodel_name = model.model_name
        self.__vbandwidth_step = model.bandwidth_step
        self.__vbandwidth = model.bandwidth
        self.__vk_range = model.k_range
        self.__vsill = model.sill

    @property
    def vmodel(self) -> CovModel:
        return self.__vmodel

    @property
    def vmodel_name(self) -> str:
        return self.__vmodel_name

    @property
    def vbandwidth_step(self) -> int:
        return self.__vbandwidth_step

    @property
    def vbandwidth(self) -> np.array:
        return self.__vbandwidth

    @property
    def vk_range(self) -> float:
        return self.__vk_range

    @property
    def vsill(self) -> float:
        return self.__vsill

    def plot(self, realizations: list[int] = None, mean=0) -> None:
        realization_number = len(self.random_field[:, 0])
        if realizations is None:
            for i in range(realization_number):
                plt.figure(self.curr_fig_num, figsize=self.figsize)
                plt.plot(self.random_field[i, :] + mean)
                plt.title('Realizations: ' + self.vmodel_name, fontsize=20)
                plt.xlabel(self.xlabel, fontsize=20)
                plt.axhline(y=mean, color='r', linestyle='--', zorder=1)
                plt.ylabel('Y', fontsize=20)
        elif type(realizations) == list:
            for item in realizations:
                plt.figure(self.curr_fig_num, figsize=self.figsize)
                plt.plot(self.random_field[:, item] + mean)
                plt.title('Realizations: ' + self.vmodel_name, fontsize=20)
                plt.xlabel(self.xlabel, fontsize=20)
                plt.axhline(y=mean, color='r', linestyle='--', zorder=1)
                plt.ylabel('Y', fontsize=20)
        else:
            raise TypeError("Argument 'realizations' should be list or None")
        self.curr_fig_num += 1

    def mean_plot(self, mean=0) -> None:
        zmean = np.zeros(len(self.random_field[0, :]))
        for i in range(len(self.random_field[0, :])):
            zmean[i] = np.mean(self.random_field[:, i] + mean)

        plt.figure(self.curr_fig_num, figsize=self.figsize)
        plt.plot(
            zmean,
            '-s',
            color='k',
            markeredgecolor='k',
            markerfacecolor='y',
        )
        plt.xlabel(self.xlabel, fontsize=20)
        plt.ylabel('Mean', fontsize=20)
        plt.axhline(y=mean, color='r', linestyle='--', zorder=1)
        plt.xticks(fontsize=17), plt.yticks(fontsize=17)
        self.curr_fig_num += 1

    def variance_plot(self) -> None:
        zvar = np.zeros(len(self.random_field[0, :]))
        for i in range(len(self.random_field[0, :])):
            zvar[i] = np.var(self.random_field[:, i])

        plt.figure(self.curr_fig_num, figsize=self.figsize)
        plt.plot(
            zvar,
            '-o',
            color='k',
            markeredgecolor='k',
            markerfacecolor='r',
        )
        plt.xlabel(self.xlabel, fontsize=20)
        plt.ylabel('Variance', fontsize=20)
        plt.axhline(y=self.vmodel.sill, color='b', linestyle='--', zorder=1)
        plt.xticks(fontsize=17), plt.yticks(fontsize=17)
        self.curr_fig_num += 1

    def cdf_plot(self, x_location: int) -> None:
        x = self.random_field[:, x_location]
        mu = np.mean(x)
        sigma = np.std(x)
        n_bins = 50

        _, ax = plt.subplots(figsize=(8, 4))

        _, bins, _ = ax.hist(
            x,
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
        self.curr_fig_num += 1

    def hist_plot(self, x_location: int) -> None:
        x = self.random_field[:, x_location]
        mu = np.mean(x)
        sigma = np.std(x)
        num_bins = 50
        plt.figure(self.curr_fig_num)
        _, bins, _ = plt.hist(
            x,
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
        self.curr_fig_num += 1

    def variogram_plot(self) -> None:
        self.__variogram_validate()
        for i in range(self.realization_number):
            plt.figure(self.curr_fig_num, figsize=(10, 6))
            plt.plot(self.variogram[i, :], alpha=0.1)
            plt.title('Model: ' + self.vmodel_name, fontsize=20)
            plt.xlabel('Lag(m)', fontsize=20)
            plt.ylabel('Variogram', fontsize=20)
            plt.xticks(fontsize=17), plt.yticks(fontsize=17)

        self.theory_variogram_plot()

        Vario_mean = np.zeros(len(self.vbandwidth))
        for i in range(len(self.vbandwidth)):
            Vario_mean[i] = np.mean(self.variogram[:, i])

        plt.plot(Vario_mean, '--', color='blue')
        self.curr_fig_num += 1

    def theory_variogram_plot(self, fig: int = None) -> None:
        if fig is not None:
            plt.figure(fig, figsize=self.figsize)
        plt.plot(
            self.vmodel.var_compute(self.vbandwidth),
            'o',
            markeredgecolor='k',
            markerfacecolor='w',
        )
        plt.title('Model: ' + self.vmodel_name, fontsize=20)
        plt.xlabel('Lag(m)', fontsize=20)
        plt.ylabel('Variogram', fontsize=20)

    def __variogram_validate(self) -> None:
        if type(self.variogram) == int:
            raise VariogramDoesNotCompute
