import numpy as np
from uc_sgsim.cov_model.base import cov_model


class Gaussian(cov_model):
    def __init__(self, bandwidth_step, bandwidth, a, C0):
        super().__init__(bandwidth_step, bandwidth, a, C0)
        self.model_name = 'Gaussian'

    def model(self, h):
        return self.C0 * (1 - np.exp(-3 * h**2 / self.a**2))


class Spherical(cov_model):
    def __init__(self, bandwidth_step, bandwidth, a, C0):
        super().__init__(bandwidth_step, bandwidth, a, C0)
        self.model_name = 'Spherical'

    def model(self, h):
        if h <= self.a:
            return self.C0 * (1.5 * h / self.a - 0.5 * (h / self.a) ** 3.0)
        else:
            return self.C0


class Exponential(cov_model):
    def __init__(self, model, bandwidth_step, bandwidth, a, C0):
        super().__init__(bandwidth_step, bandwidth, a, C0)
        self.model_name = 'Exponential'


class Circular(cov_model):
    def __init__(self, model, bandwidth_step, bandwidth, a, C0):
        super().__init__(bandwidth_step, bandwidth, a, C0)
        self.model_name = 'Circular'


class Linear(cov_model):
    def __init__(self, model, bandwidth_step, bandwidth, a, C0):
        super().__init__(bandwidth_step, bandwidth, a, C0)
        self.model_name = 'Circular'
