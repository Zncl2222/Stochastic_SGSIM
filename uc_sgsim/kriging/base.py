from __future__ import annotations

import numpy as np
from uc_sgsim.cov_model.base import CovModel


class Kriging:
    """
    Base Class for Kriging Interpolation

    This class serves as the parent class for various kriging interpolation techniques
    such as Simple Kriging and Ordinary Kriging. It define and encapsulates common properties and
    methods used in kriging.

    Attributes:
        model (CovModel): The covariance model used for interpolation.
        bandwidth_step (int): The step size for bandwidth increments.
        bandwidth (np.array): An array of bandwidth values.
        k_range (float): The range parameter for the covariance model.
        sill (float): The sill parameter for the covariance model.
        x_size (int): Size of the x-axis for the interpolation grid.
        y_size (int): Size of the y-axis for the interpolation grid (0 if 1D).
        _cov_cache_flag (bool): Flag indicating whether to use a covariance cache.
        _cov_cache (list | list[list]): Cache for computed covariances (if enabled).

    Methods:
        _create_cov_cache(): Create the covariance cache for faster computations.
    """

    def __init__(
        self,
        model: CovModel,
        grid_size: int | list[int, int],
        cov_cache: bool = False,
    ):
        self._model = model
        self._bandwidth_step = model.bandwidth_step
        self._bandwidth = model.bandwidth
        self._k_range = model.k_range
        self._sill = model.sill
        self.x_size = grid_size if isinstance(grid_size, int) else grid_size[0]
        self.y_size = 0 if isinstance(grid_size, int) else grid_size[1]
        self._cov_cache_flag = cov_cache
        if cov_cache is True:
            self._cov_cache = {}

    @property
    def model(self) -> CovModel:
        return self._model

    @property
    def bandwidth_step(self) -> int:
        return self._bandwidth_step

    @property
    def bandwidth(self) -> np.array:
        return self._bandwidth

    @property
    def k_range(self) -> float:
        return self._k_range

    @property
    def sill(self) -> float:
        return self._sill

    def __repr__(self):
        return f'{self.__class__.__name__}'
