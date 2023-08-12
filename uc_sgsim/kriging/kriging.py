import numpy as np
from scipy.spatial.distance import pdist, squareform
from uc_sgsim.cov_model.base import CovModel
from uc_sgsim.kriging.base import Kriging


class SimpleKriging(Kriging):
    def __init__(self, model: CovModel):
        super().__init__(model)

    def prediction(self, sample: np.array, unsampled: np.array) -> tuple[float, float]:
        n_sampled = len(sample)
        dist_diff = np.sqrt((sample[:, 0] - unsampled[0]) ** 2 + (sample[:, 1] - unsampled[1]) ** 2)
        dist_diff = dist_diff.reshape(len(dist_diff), 1)

        grid = np.hstack([sample, dist_diff])
        meanvalue = 0

        cov_dist = np.array(self.model.cov_compute(grid[:, 3])).reshape(-1, 1)
        cov_data = squareform(pdist(grid[:, :2])).flatten()
        cov_data = np.array(self.model.cov_compute(cov_data))
        cov_data = cov_data.reshape(n_sampled, n_sampled)
        # cov_data[np.diag_indices_from(cov_data)] += 1e-4

        weights = np.linalg.solve(cov_data, cov_dist)
        residuals = grid[:, 1] - meanvalue
        estimation = np.dot(weights.T, residuals) + meanvalue
        kriging_var = float(self.model.sill - np.dot(weights.T, cov_dist))

        if kriging_var < 0:
            kriging_var = 0

        kriging_std = np.sqrt(kriging_var)

        return estimation, kriging_std

    def simulation(self, x: np.array, unsampled: np.array, **kwargs) -> float:
        neighbor = kwargs.get('neighbor')
        if neighbor is not None:
            dist = np.sqrt((x[:, 0] - unsampled[0]) ** 2 + (x[:, 1] - unsampled[1]) ** 2)
            dist = dist.reshape(len(dist), 1)
            has_neighbor = self.find_neighbor(dist, neighbor)
            if has_neighbor:
                return has_neighbor
            x = np.hstack([x, dist])
            sorted_indices = np.argpartition(x[:, 2], neighbor)[:neighbor]
            x = x[sorted_indices]
            # print(x)
        estimation, kriging_std = self.prediction(x, unsampled)

        random_fix = np.random.normal(0, kriging_std, 1)
        return estimation + random_fix

    def find_neighbor(self, dist: list[float], neighbor: int) -> float:
        if neighbor == 0:
            return np.random.normal(0, self.model.sill**0.5, 1)
        close_point = 0

        criteria = self.k_range * 1.732 if self.model.model_name == 'Gaussian' else self.k_range

        for item in dist:
            if item <= criteria:
                close_point += 1

        if close_point == 0:
            return np.random.normal(0, self.model.sill**0.5, 1)


# class SimpleKriging(Kriging):
#     def __init__(self, model: CovModel):
#         super().__init__(model)

#     def prediction(self, sample: np.array, unsampled: np.array) -> tuple[float, float]:
#         n_sampled = len(sample)
#         dist_diff = abs(sample[:, 0] - unsampled)
#         dist_diff = dist_diff.reshape(len(dist_diff), 1)

#         grid = np.hstack([sample, dist_diff])
#         meanvalue = 0

#         cov_dist = np.array(self.model.cov_compute(grid[:, 2])).reshape(-1, 1)
#         cov_data = squareform(pdist(grid[:, :1])).flatten()
#         cov_data = np.array(self.model.cov_compute(cov_data))
#         cov_data = cov_data.reshape(n_sampled, n_sampled)
#         # Add a small nugget to the diagonal of the covariance matrix for numerical stability
#         cov_data[np.diag_indices_from(cov_data)] += 1e-4

#         weights = np.linalg.solve(cov_data, cov_dist)
#         residuals = grid[:, 1] - meanvalue
#         estimation = np.dot(weights.T, residuals) + meanvalue
#         kriging_var = float(self.model.sill - np.dot(weights.T, cov_dist))

#         if kriging_var < 0:
#             kriging_var = 0

#         kriging_std = np.sqrt(kriging_var)

#         return estimation, kriging_std

#     def simulation(self, x: np.array, unsampled: np.array, **kwargs) -> float:
#         neighbor = kwargs.get('neighbor')
#         if neighbor is not None:
#             dist = abs(x[:, 0] - unsampled)
#             dist = dist.reshape(len(dist), 1)
#             has_neighbor = self.find_neighbor(dist, neighbor)
#             if has_neighbor:
#                 return has_neighbor
#             x = np.hstack([x, dist])
#             sorted_indices = np.argpartition(x[:, 2], neighbor)[:neighbor]
#             x = x[sorted_indices]

#         estimation, kriging_std = self.prediction(x, unsampled)

#         random_fix = np.random.normal(0, kriging_std, 1)
#         return estimation + random_fix

#     def find_neighbor(self, dist: list[float], neighbor: int) -> float:
#         if neighbor == 0:
#             return np.random.normal(0, self.model.sill**0.5, 1)
#         close_point = 0

#         criteria = self.k_range * 1.732 if self.model.model_name == 'Gaussian' else self.k_range

#         for item in dist:
#             if item <= criteria:
#                 close_point += 1

#         if close_point == 0:
#             return np.random.normal(0, self.model.sill**0.5, 1)


class OrdinaryKriging(SimpleKriging):
    def prediction(self, sample: np.array, unsampled: np.array) -> tuple[float, float]:
        n_sampled = len(sample)
        dist_diff = abs(sample[:, 0] - unsampled)
        dist_diff = dist_diff.reshape(len(dist_diff), 1)

        grid = np.hstack([sample, dist_diff])

        cov_dist = np.array(self.model.cov_compute(grid[:, 2])).reshape(-1, 1)
        cov_data = squareform(pdist(grid[:, :1])).flatten()
        cov_data = np.array(self.model.cov_compute(cov_data))
        cov_data = cov_data.reshape(n_sampled, n_sampled)

        # Add a small value to the diagonal of the covariance matrix for numerical stability
        cov_data[np.diag_indices_from(cov_data)] += 1e-4

        cov_data_augmented = self.matrix_agumented(cov_data)
        cov_dist_augmented = np.vstack((cov_dist, [1.0]))
        weights = np.linalg.solve(cov_data_augmented, cov_dist_augmented)[:n_sampled]

        estimation = np.dot(weights.T, grid[:, 1])
        kriging_var = float(self.model.sill - np.dot(weights.T, cov_dist))

        if kriging_var < 0:
            kriging_var = 0

        kriging_std = np.sqrt(kriging_var)

        return estimation, kriging_std

    def matrix_agumented(self, mat: np.array):
        ones_column = np.ones((mat.shape[0], 1))
        cov_data_augmented = np.hstack([mat, ones_column])
        ones_row = np.ones((1, cov_data_augmented.shape[1]))
        ones_row[0][-1] = 0
        cov_data_augmented = np.vstack((cov_data_augmented, ones_row))
        return cov_data_augmented
