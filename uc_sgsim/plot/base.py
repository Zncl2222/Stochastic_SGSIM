import numpy as np
from ..cov_model.base import CovModel


class PlotBase:
    def __init__(
        self,
        model: CovModel,
        figsize: tuple = (10, 8),
    ):
        self.__model = model
        self.__figsize = figsize
        self.__model_name = model.model_name
        self.__bandwidth_step = model.bandwidth_step
        self.__bandwidth = model.bandwidth
        self.__k_range = model.k_range
        self.__sill = model.sill

    @property
    def model(self) -> CovModel:
        return self.__model

    @property
    def figsize(self) -> tuple:
        return self.__figsize

    @property
    def model_name(self) -> str:
        return self.__model_name

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
