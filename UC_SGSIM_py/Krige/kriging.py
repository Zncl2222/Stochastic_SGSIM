from UC_SGSIM_py.Krige.base import Kriging
from scipy.spatial.distance import pdist, squareform
import numpy as np

class SimpleKrige(Kriging):

    def __init__(self, model):

        super().__init__(model)

    def compute(self, L, u, N, randomseed):

        np.random.seed(randomseed)

        if N == 0:
            return np.random.normal(0,1,1)

        dist = abs(L[:,0] - u)
        dist = dist.reshape(len(dist), 1)

        close_point = 0

        for item in dist:
            if item <= self.a:
                close_point += 1

        if close_point == 0:
            return np.random.normal(0,1,1)

        L = np.hstack([L, dist])
        L = np.array(sorted(L, key = lambda x : x[2])[:N])

        meanvalue = 0

        Cov_dist = np.matrix(self.model.Cov_compute(L[:,2])).T

        Cov_data = squareform(pdist(L[:,:1])).flatten()
        Cov_data = np.array(self.model.Cov_compute(Cov_data))
        Cov_data = Cov_data.reshape(N, N)

        weights = np.linalg.inv(Cov_data) * Cov_dist

        residuals = (L[:,1] - meanvalue)

        estimation = np.dot(weights.T, residuals) + meanvalue

        krige_var = float(1 - np.dot(weights.T, Cov_dist))

        if krige_var < 0:
            krige_var = 0

        krige_std = np.sqrt(krige_var)

        random_fix = np.random.normal(0, krige_std, 1)

        Simulated = float(estimation + random_fix)

        return Simulated

        

            


