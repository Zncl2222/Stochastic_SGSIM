from __future__ import annotations

import numpy as np
from uc_sgsim.plotting import SgsimPlot
from uc_sgsim.exception import VariogramDoesNotCompute
from uc_sgsim.kriging import SimpleKriging, OrdinaryKriging, Kriging
from uc_sgsim.utils import save_as_multiple_file, save_as_one_file
from uc_sgsim.cov_model.base import CovModel


class RandomField:
    def __init__(self, grid_size: int | list[int, int], realization_number: int):
        self.__realization_number = realization_number
        self.__grid_size = grid_size
        self._create_grid(grid_size)

    def _create_grid(self, grid_size: int | list[int, int]) -> None:
        self.__x = range(grid_size) if isinstance(grid_size, int) else range(grid_size[0])
        self.__y = 0 if isinstance(grid_size, int) else range(grid_size[1])
        self.__x_size = len(self.__x)
        self.__y_size = 0 if isinstance(grid_size, int) else len(self.__y)
        self.random_field = np.empty([self.__realization_number, self.__x_size])
        self.variogram = 0

    @property
    def grid_size(self) -> int | list[int, int]:
        return self.__grid_size

    @property
    def x(self) -> int:
        return self.__x

    @property
    def x_size(self) -> int:
        return self.__x_size

    @property
    def y(self) -> int:
        return self.__y

    @property
    def y_size(self) -> int:
        return self.__y_size

    @property
    def realization_number(self) -> int:
        return self.__realization_number

    @realization_number.setter
    def realization_number(self, val: int):
        self.__realization_number = val

    def save_random_field(
        self,
        path: str,
        file_type: str = 'csv',
        save_single: bool = False,
    ) -> None:
        digit = int(np.log10(self.realization_number))
        number_head = ''
        for i in range(digit):
            number_head += '0'
        num_val = 1
        if save_single is False:
            for i in range(self.realization_number):
                if i // num_val == 10:
                    num_val *= 10
                    number_head = number_head[:-1]
                number = number_head + str(i)
                save_as_multiple_file(
                    number,
                    self.x_size,
                    self.random_field,
                    file_type,
                    'Realizations',
                )
        else:
            save_as_one_file(path, self.random_field)

    def save_variogram(self, path: str, file_type: str = 'csv', save_single: bool = False) -> None:
        if type(self.variogram) == int:
            raise VariogramDoesNotCompute()
        digit = int(np.log10(self.realization_number))
        number_head = ''
        for i in range(digit):
            number_head += '0'
        num_val = 1
        if save_single is False:
            for i in range(self.realization_number):
                if i // num_val == 10:
                    num_val *= 10
                    number_head = number_head[:-1]
                number = number_head + str(i)
                save_as_multiple_file(
                    number,
                    len(self.bandwidth),
                    self.variogram,
                    file_type,
                    'Variogram',
                )
        else:
            save_as_one_file(path, self.variogram)


class SgsimField(RandomField, SgsimPlot):
    def __init__(
        self,
        grid_size: int | list[int, int],
        realization_number: int,
        model: CovModel,
        kriging: str | Kriging = 'SimpleKriging',
        **kwargs,
    ):
        RandomField.__init__(self, grid_size, realization_number)
        SgsimPlot.__init__(self, model)

        self.__model = model
        self.__bandwidth_step = model.bandwidth_step
        self.__bandwidth = model.bandwidth
        self.__set_kriging_method(kriging, **kwargs)
        self.__set_kwargs(**kwargs)

    def __set_kriging_method(self, kriging, **kwargs) -> None:
        self.constant_path = kwargs.get('constant_path', False)
        cov_cache = kwargs.get('cov_cache', False)

        if self.constant_path is False and cov_cache is True:
            raise ValueError('cov_cache should be False when constant_path is False')

        if kriging == 'SimpleKriging':
            self.kriging = SimpleKriging(self.model, self.grid_size, cov_cache)
        elif kriging == 'OrdinaryKriging':
            self.kriging = OrdinaryKriging(self.model, self.grid_size, cov_cache)
        else:
            if not isinstance(kriging, (SimpleKriging, OrdinaryKriging)):
                raise TypeError('Kriging should be class SimpleKriging or OrdinaryKriging')

    def __set_kwargs(self, **kwargs) -> None:
        self.z_min = kwargs.get('z_min', -(self.model.sill**0.5 * 4))
        self.z_max = kwargs.get('z_max', self.model.sill**0.5 * 4)
        self.max_neighbor = kwargs.get('max_neighbor', 8)

    @property
    def model(self) -> CovModel:
        return self.__model

    @property
    def bandwidth(self) -> np.array:
        return self.__bandwidth

    @property
    def bandwidth_step(self) -> int:
        return self.__bandwidth_step
