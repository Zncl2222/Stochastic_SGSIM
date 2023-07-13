import numpy as np
from ..cov_model.base import CovModel


class PlotBase:
    def __init__(
        self,
        model: CovModel,
        figsize: tuple = (10, 8),
    ):
        self.__vmodel = model
        self.__figsize = figsize
        self.__vmodel_name = model.model_name
        self.__vbandwidth_step = model.bandwidth_step
        self.__vbandwidth = model.bandwidth
        self.__vk_range = model.k_range
        self.__vsill = model.sill

    @property
    def vmodel(self) -> CovModel:
        return self.__vmodel

    @property
    def figsize(self) -> tuple:
        return self.__figsize

    @property
    def vmodel_name(self) -> str:
        return self.__vmodel_name

    @property
    def vbandwidth_step(self) -> int:
        return self.__vbandwidth_step

    @property
    def vbandwidth(self) -> np.array:
        return self.__vbandwidth

    @property
    def vk_range(self) -> float:
        return self.__vk_range

    @property
    def sill(self) -> float:
        return self.__vsill
