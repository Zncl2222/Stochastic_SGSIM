from __future__ import annotations

import numpy as np
from uc_sgsim.cov_model.base import CovModel


class Kriging:
    def __init__(
        self,
        model: CovModel,
        grid_size: int | list[int, int],
        cov_cache: bool = False,
    ):
        self.__model = model
        self.__bandwidth_step = model.bandwidth_step
        self.__bandwidth = model.bandwidth
        self.__k_range = model.k_range
        self.__sill = model.sill
        self.x_size = grid_size if isinstance(grid_size, int) else grid_size[0]
        self.y_size = 0 if isinstance(grid_size, int) else grid_size[1]
        self._cov_cache_flag = cov_cache
        if cov_cache is True:
            self._create_cov_cache()

    @property
    def model(self) -> CovModel:
        return self.__model

    @property
    def bandwidth_step(self) -> int:
        return self.__bandwidth_step

    @property
    def bandwidth(self) -> np.array:
        return self.__bandwidth

    @property
    def k_range(self) -> float:
        return self.__k_range

    @property
    def sill(self) -> float:
        return self.__sill

    def _create_cov_cache(self):
        if self.y_size != 0:
            self._cov_cache = [[0] * self.x_size for _ in range(self.y_size)]
        else:
            self._cov_cache = [0] * self.x_size
