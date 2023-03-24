import numpy as np
from uc_sgsim.cov_model.base import CovModel


class Gaussian(CovModel):
    def __init__(
        self,
        bandwidth_step: int,
        bandwidth: np.array,
        k_range: float,
        sill: float,
    ):
        super().__init__(bandwidth_step, bandwidth, k_range, sill)
        self.model_name = 'Gaussian'

    def model(self, h: float) -> float:
        return self.sill * (1 - np.exp(-3 * h**2 / self.k_range**2))


class Spherical(CovModel):
    def __init__(
        self,
        bandwidth_step: int,
        bandwidth: np.array,
        k_range: float,
        sill: float,
    ):
        super().__init__(bandwidth_step, bandwidth, k_range, sill)
        self.model_name = 'Spherical'

    def model(self, h: float) -> float:
        if h <= self.k_range:
            return self.sill * (1.5 * h / self.k_range - 0.5 * (h / self.k_range) ** 3.0)
        else:
            return self.sill


class Exponential(CovModel):
    def __init__(self, bandwidth_step, bandwidth, a, sill):
        super().__init__(bandwidth_step, bandwidth, a, sill)
        self.model_name = 'Exponential'


class Circular(CovModel):
    def __init__(self, bandwidth_step, bandwidth, a, sill):
        super().__init__(bandwidth_step, bandwidth, a, sill)
        self.model_name = 'Circular'


class Linear(CovModel):
    def __init__(self, bandwidth_step, bandwidth, a, sill):
        super().__init__(bandwidth_step, bandwidth, a, sill)
        self.model_name = 'Circular'
