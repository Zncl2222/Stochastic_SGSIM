import numpy as np
from scipy.spatial.distance import pdist, squareform


class CovModel:
    def __init__(self, bandwidth_step: int, bandwidth: np.array, a: float, sill=1):
        self.__bandwidth_step = bandwidth_step
        self.__bandwidth = bandwidth
        self.__a = a
        self.__sill = sill

    @property
    def bandwidth_step(self) -> int:
        return self.__bandwidth_step

    @property
    def bandwidth(self) -> np.array:
        return self.__bandwidth

    @property
    def a(self) -> float:
        return self.__a

    @property
    def sill(self) -> float:
        return self.__sill

    def cov_compute(self, x: np.array) -> float:
        z = np.empty(len(x))
        for i in range(len(x)):
            z[i] = self.__sill - self.model(x[i])

        return z

    def var_compute(self, x: np.array) -> float:
        z = np.empty(len(x))
        for i in range(len(x)):
            z[i] = self.model(x[i])

        return z

    def variogram(self, x: np.array) -> np.array:
        dist = squareform(pdist(x[:, :1]))
        variogram = []

        for h in self.__bandwidth_step:
            z = []
            for i in range(len(dist[:, 0])):
                for j in range(i + 1, len(dist[:, 0])):
                    if (dist[i, j] >= h - self.__bandwidth) and (
                        dist[i, j] <= h + self.__bandwidth
                    ):
                        z.append(np.power(x[i, 1] - x[j, 1], 2))
            if np.sum(z) >= 1e-7:
                variogram.append(np.sum(z) / (2 * len(z)))

        return np.array(variogram)
