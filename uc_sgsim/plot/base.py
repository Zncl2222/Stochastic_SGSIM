import matplotlib.pyplot as plt


class PlotBase:
    def __init__(self, figsize: tuple = (10, 8)):
        self.__figsize = figsize

    @property
    def figsize(self) -> tuple:
        return self.__figsize

    def save_plot(self, figname: str = 'figure.png', **kwargs):
        plt.savefig(figname, **kwargs)
