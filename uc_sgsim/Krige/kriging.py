import numpy as np
from scipy.spatial.distance import pdist, squareform
from uc_sgsim.krige.base import Kriging


class SimpleKrige(Kriging):
    def __init__(self, model):
        super().__init__(model)

    def prediction(self, sample, unsampled):
        n_sampled = len(sample)
        dist_diff = abs(sample[:, 0] - unsampled)
        dist_diff = dist_diff.reshape(len(dist_diff), 1)

        L = np.hstack([sample, dist_diff])
        meanvalue = 0

        Cov_dist = np.matrix(self.model.cov_compute(L[:, 2])).T
        Cov_data = squareform(pdist(L[:, :1])).flatten()
        Cov_data = np.array(self.model.cov_compute(Cov_data))
        Cov_data = Cov_data.reshape(n_sampled, n_sampled)

        weights = np.linalg.inv(Cov_data) * Cov_dist
        residuals = L[:, 1] - meanvalue
        estimation = np.dot(weights.T, residuals) + meanvalue
        krige_var = float(1 - np.dot(weights.T, Cov_dist))

        if krige_var < 0:
            krige_var = 0

        krige_std = np.sqrt(krige_var)

        return estimation, krige_std

    def simulation(self, x, unsampled, **kwargs):
        neighbor = kwargs.get('neighbor')
        if neighbor:
            dist = abs(x[:, 0] - unsampled)
            dist = dist.reshape(len(dist), 1)
            has_neighbor = self.find_neighbor(dist, neighbor)
            if has_neighbor:
                return has_neighbor
            x = np.hstack([x, dist])
            x = np.array(sorted(x, key=lambda itr: itr[2])[:neighbor])

        estimation, krige_std = self.prediction(x, unsampled)

        random_fix = np.random.normal(0, krige_std, 1)
        return estimation + random_fix

    def find_neighbor(self, dist, neighbor):
        if neighbor == 0:
            return np.random.normal(0, 1, 1)
        close_point = 0
        for item in dist:
            if item <= self.a:
                close_point += 1

        if close_point == 0:
            return np.random.normal(0, 1, 1)
