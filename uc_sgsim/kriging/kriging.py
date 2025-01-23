from __future__ import annotations

import numpy as np
from scipy.spatial.distance import pdist, squareform
from uc_sgsim.kriging.base import Kriging


class SimpleKriging(Kriging):
    """
    Simple Kriging Class

    This class represents the Simple Kriging interpolation technique, which is used
    to estimate values at unsampled locations based on sampled data and a specified
    covariance model.

    Attributes:
        model (CovModel): The covariance model used for interpolation.
        grid_size (int | list[int, int]): Size of the grid for interpolation.
        cov_cache (bool): Flag indicating whether to use a covariance cache.

    Methods:
        prediction(sample: np.array, unsampled: np.array) -> tuple[float, float]:
            Perform Simple Kriging prediction for an unsampled location.
        simulation(x: np.array, unsampled: np.array, **kwargs) -> float:
            Perform Simple Kriging simulation for unsampled locations.
        _find_neighbor(dist: list[float], neighbor: int) -> float:
            Find a nearby point for simulation based on a neighbor criterion.
    """

    def prediction(
        self,
        unsampled: np.ndarray,
        sampled: np.ndarray,
        dist_diff: np.ndarray = None,
    ) -> tuple[float, float]:
        """
        Perform Simple Kriging prediction for an unsampled location.

        Args:
            unsampled (np.array): Unsamped location for which prediction is made.
            sample (np.array): Sampled data for neighboring locations.

        Returns:
            tuple[float, float]: Estimated value and kriging standard deviation.
        """
        n_sampled = len(sampled)
        dist_diff = (
            dist_diff
            if dist_diff is not None
            else np.linalg.norm(unsampled - sampled[:, [0, 1]], axis=1).flatten()
        )
        sampled[:, 3] = dist_diff

        meanvalue = 0
        if self._cov_cache_flag is True:
            cov_dist = np.array(self.model.cov_compute(sampled[:, 3])).reshape(-1, 1)
            self._cov_cache[f'{unsampled[0]}, {unsampled[1]}'] = cov_dist
            # print(cov_dist)
        elif hasattr(self, '_cov_cache') is True:
            cov_dist = self._cov_cache[f'{unsampled[0]}, {unsampled[1]}']
            # print("HQWIPE", cov_dist)
        else:
            cov_dist = np.array(self.model.cov_compute(sampled[:, 3])).reshape(-1, 1)
        # print(cov_dist)
        cov_data = squareform(pdist(sampled[:, :2])).flatten()
        cov_data = np.array(self.model.cov_compute(cov_data))
        cov_data = cov_data.reshape(n_sampled, n_sampled)
        # Add a small nugget to the diagonal of the covariance matrix for numerical stability
        cov_data[np.diag_indices_from(cov_data)] += 1e-4

        weights = np.linalg.solve(cov_data, cov_dist)
        residuals = sampled[:, 2] - meanvalue
        estimation = np.dot(weights.T, residuals) + meanvalue
        kriging_var = float(self.model.sill - np.dot(weights.T, cov_dist))

        kriging_var = kriging_var if kriging_var > 0 else 0
        kriging_std = np.sqrt(kriging_var)

        return estimation, kriging_std

    def simulation(self, unsampled: np.ndarray, sampled: np.ndarray, **kwargs) -> float:
        """
        Perform Simple Kriging simulation for unsampled locations.

        Args:
            unsampled (np.ndarray): Sampled data for neighboring locations.
            unsampled (np.ndarray): Unsamped location for which simulation is performed.
            neighbor (int): The number of neighbors to consider (optional).

        Returns:
            float: Simulated value for the unsampled location.
        """
        if len(sampled) == 0:
            return np.random.normal(0, self.model.sill**0.5, 1)

        neighbor = kwargs.get('neighbor')
        if neighbor is not None:
            distances = np.linalg.norm(unsampled - sampled[:, [0, 1]], axis=1)

            draw_random_normal = self._find_neighbor(distances, neighbor)
            if draw_random_normal:
                return draw_random_normal

            sorted_indices = np.argsort(distances)
            sampled = np.array(sampled)[sorted_indices][:neighbor]
            distances = distances[sorted_indices][:neighbor]

        dist_diff = distances if neighbor is not None else None
        estimation, kriging_std = self.prediction(unsampled, sampled, dist_diff)

        random_fix = np.random.normal(0, kriging_std, 1)
        return estimation + random_fix

    def _find_neighbor(self, distances: list[float], neighbor: int) -> float | None:
        """
        Find a nearby point for simulation based on a neighbor criterion.

        Args:
            distances (list[float]): Distances from sampled points to the unsampled location.
            neighbor (int): The number of neighbors to consider.

        Returns:
            float: Simulated value based on neighbors or a random value if no neighbors are found.
        """
        if neighbor == 0:
            return np.random.normal(0, self.model.sill**0.5, 1)
        close_point = 0

        criteria = self.k_range * 1.732 if self.model.model_name == 'Gaussian' else self.k_range

        for item in distances:
            if item <= criteria:
                close_point += 1

        if close_point == 0:
            return np.random.normal(0, self.model.sill**0.5, 1)


class OrdinaryKriging(SimpleKriging):
    """
    Ordinary Kriging Class

    This class represents the Ordinary Kriging interpolation technique, which is
    an extension of Simple Kriging. It is used to estimate values at unsampled
    locations based on sampled data and a specified covariance model.

    Attributes:
        model (CovModel): The covariance model used for interpolation.
        grid_size (int | list[int, int]): Size of the grid for interpolation.
        cov_cache (bool): Flag indicating whether to use a covariance cache.

    Methods:
        prediction(sample: np.array, unsampled: np.array) -> tuple[float, float]:
            Perform Ordinary Kriging prediction for an unsampled location.
        matrix_augmented(mat: np.array) -> np.array:
            Augment the covariance matrix for Ordinary Kriging.
    """

    def prediction(
        self,
        unsampled: np.ndarray,
        sampled: np.ndarray,
        dist_diff: np.ndarray = None,
    ) -> tuple[float, float]:
        """
        Perform Ordinary Kriging prediction for an unsampled location.

        Args:
            unsampled (np.array): Unsamped location for which prediction is made.
            sample (np.array): Sampled data for neighboring locations.

        Returns:
            tuple[float, float]: Estimated value and kriging standard deviation.
        """
        n_sampled = len(sampled)
        dist_diff = (
            dist_diff
            if dist_diff is not None
            else np.linalg.norm(unsampled - sampled[:, [0, 1]], axis=1).flatten()
        )
        sampled[:, 3] = dist_diff

        if self._cov_cache_flag:
            cov_dist = np.array(self.model.cov_compute(sampled[:, 3])).reshape(-1, 1)
            self._cov_cache[f'{unsampled[0]}, {unsampled[1]}'] = cov_dist
        elif hasattr(self, '_cov_cache'):
            cov_dist = self._cov_cache[f'{unsampled[0]}, {unsampled[1]}']
        else:
            cov_dist = np.array(self.model.cov_compute(sampled[:, 3])).reshape(-1, 1)

        cov_data = squareform(pdist(sampled[:, :2])).flatten()
        cov_data = np.array(self.model.cov_compute(cov_data))
        cov_data = cov_data.reshape(n_sampled, n_sampled)

        # Add a small value to the diagonal of the covariance matrix for numerical stability
        cov_data[np.diag_indices_from(cov_data)] += 1e-4

        cov_data_augmented = self._matrix_augmented(cov_data)
        cov_dist_augmented = np.vstack((cov_dist, [1.0]))
        weights = np.linalg.solve(cov_data_augmented, cov_dist_augmented)[:n_sampled]

        estimation = np.dot(weights.T, sampled[:, 2])
        kriging_var = float(self.model.sill - np.dot(weights.T, cov_dist))

        kriging_var = kriging_var if kriging_var > 0 else 0
        kriging_std = np.sqrt(kriging_var)

        return estimation, kriging_std

    def _matrix_augmented(self, mat: np.ndarray) -> np.ndarray:
        """
        Augment the covariance matrix to constrain summation of weights = 1 for Ordinary Kriging.

        Args:
            mat (np.array): Covariance matrix.

        Returns:
            np.array: Augmented covariance matrix.
        """
        ones_column = np.ones((mat.shape[0], 1))
        cov_data_augmented = np.hstack([mat, ones_column])
        ones_row = np.ones((1, cov_data_augmented.shape[1]))
        ones_row[0][-1] = 0
        cov_data_augmented = np.vstack((cov_data_augmented, ones_row))
        return cov_data_augmented
