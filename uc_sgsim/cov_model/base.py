import numpy as np
from scipy.spatial.distance import pdist, squareform


class CovModel:
    def __init__(self, bandwidth_step: int, bandwidth: np.array, a: float, C0=1):
        self.__bandwidth_step = bandwidth_step
        self.__bandwidth = bandwidth
        self.__a = a
        self.__C0 = C0

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
    def C0(self) -> float:
        return self.__C0

    def cov_compute(self, Y: np.array) -> float:
        Z = np.empty(len(Y))
        for i in range(len(Y)):
            Z[i] = self.__C0 - self.model(Y[i])

        return Z

    def var_compute(self, Y: np.array) -> float:
        Z = np.empty(len(Y))
        for i in range(len(Y)):
            Z[i] = self.model(Y[i])

        return Z

    def variogram(self, Y: np.array) -> np.array:
        dist = squareform(pdist(Y[:, :1]))
        variogram = []

        for h in self.__bandwidth_step:
            Z = []
            for i in range(len(dist[:, 0])):
                for j in range(i + 1, len(dist[:, 0])):
                    if (dist[i, j] >= h - self.__bandwidth) and (
                        dist[i, j] <= h + self.__bandwidth
                    ):
                        Z.append(np.power(Y[i, 1] - Y[j, 1], 2))
            if np.sum(Z) >= 1e-7:
                variogram.append(np.sum(Z) / (2 * len(Z)))

        return np.array(variogram)
