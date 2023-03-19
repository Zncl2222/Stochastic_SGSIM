import pytest

import numpy as np
from uc_sgsim.krige import SimpleKrige
from uc_sgsim.cov_model import Gaussian


@pytest.mark.krige
class TestSimpleKrige:
    @classmethod
    def setup_class(cls) -> None:
        cls.bw = 1
        cls.hs = np.arange(0.0, 35, cls.bw)
        cls.a = 17.32
        cls.C0 = 1
        cls.gaussian = Gaussian(cls.hs, cls.bw, cls.a, cls.C0)
        cls.krige = SimpleKrige(cls.gaussian)
        cls.x_len = 150
        cls.x_grid = np.linspace(0, cls.x_len, cls.x_len + 1).reshape(cls.x_len + 1, 1)
        cls.sampled = np.array(
            [[0, 10, 29, 74, 85, 115, 146], [-1, 0.7, -0.2, 0.065, 1, 3, 0.5]],
        ).T

    def test_simple_krige_prediction(self):
        res = np.empty(self.x_len)
        for i in range(len(self.sampled[0])):
            res[int(self.sampled.T[0][i])] = self.sampled.T[1][i]

        for i in range(self.x_len):
            if i not in self.sampled[0]:
                res[i], krige_std = self.krige.prediction(self.sampled, i)
                assert krige_std >= 0
