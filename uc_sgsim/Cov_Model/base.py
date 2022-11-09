import numpy as np
from scipy.spatial.distance import pdist, squareform


class cov_model:

    def __init__(self, hs, bw, a, C0=1):
        self.hs = hs
        self.bw = bw
        self.a = a
        self.C0 = C0

    def Cov_compute(self, Y):
        Z = np.empty(len(Y))
        for i in range(len(Y)):
            Z[i] = self.C0 - self.model(Y[i])

        return Z

    def Var_compute(self, Y):
        Z = np.empty(len(Y))
        for i in range(len(Y)):
            Z[i] = self.model(Y[i])

        return Z

    def Variogram(self, Y):
        dist = squareform(pdist(Y[:, :1]))
        variogram = []

        for h in self.hs:
            Z = []
            for i in range(len(dist[:, 0])):
                for j in range(i+1, len(dist[:, 0])):
                    if ((dist[i, j] >= h - self.bw) and
                            (dist[i, j] <= h + self.bw)):
                        Z.append(np.power(Y[i, 1]-Y[j, 1], 2))
            if np.sum(Z) >= 1e-7:
                variogram.append(np.sum(Z)/(2*len(Z)))

        return np.array(variogram)
