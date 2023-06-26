import pytest

import numpy as np
from uc_sgsim.cov_model import Gaussian, Spherical


@pytest.mark.cov_model
class TestCovModel:
    @classmethod
    def setup_class(cls) -> None:
        cls.bw = 1
        cls.hs = 35
        cls.a = 17.32
        cls.C0 = 1
        cls.gaussian = Gaussian(cls.hs, cls.bw, cls.a, cls.C0)
        cls.spherical = Spherical(cls.hs, cls.bw, cls.a, cls.C0)
        cls.test_array = np.array([1, 2, 3, 4])

    def test_properties(self):
        assert len(self.gaussian.bandwidth), self.hs / self.bw
        assert self.gaussian.bandwidth_step, self.bw
        assert self.gaussian.bandwidth_len, self.hs
        assert self.gaussian.k_range, self.a
        assert self.gaussian.sill, self.C0

    def test_cov_compute(self):
        ans_gaussian = [0.99004925, 0.96078718, 0.91392636, 0.85213579]
        cov_gaussian = self.gaussian.cov_compute(self.test_array)
        assert list(cov_gaussian) == pytest.approx(ans_gaussian)

        ans_spherical = [0.91349115, 0.82755971, 0.74278306, 0.65973862]
        cov_spherical = self.spherical.cov_compute(self.test_array)
        assert list(cov_spherical) == pytest.approx(ans_spherical)

    def test_var_compute(self):
        ans_gaussian = [0.00995075, 0.03921282, 0.08607364, 0.14786421]
        var_gaussian = self.gaussian.var_compute(self.test_array)
        assert list(var_gaussian) == pytest.approx(ans_gaussian)

        ans_spherical = [0.08650885, 0.17244029, 0.25721694, 0.34026138]
        var_spherical = self.spherical.var_compute(self.test_array)
        assert list(var_spherical) == pytest.approx(ans_spherical)

    def test_variogram(self):
        array_col1 = self.test_array.reshape(len(self.test_array), 1)
        array_col2 = np.array([10, 20, 30, 40]).reshape(len(self.test_array), 1)
        array_col2_2 = np.array([7, 6, 23, 99]).reshape(len(self.test_array), 1)
        array2d = np.hstack([array_col1, array_col2])
        ans_gaussian = [30.0, 78.57142857, 166.66666667, 283.33333333, 450.0]
        vario_gaussian = self.gaussian.variogram(array2d)
        assert list(vario_gaussian) == pytest.approx(ans_gaussian)

        array2d = np.hstack([array_col1, array_col2_2])
        ans_spherical = [606.6, 1069.35714286, 1952.91666667, 2894.83333333, 4232.0]
        vario_spherical = self.spherical.variogram(array2d)
        assert list(vario_spherical) == pytest.approx(ans_spherical)

    def test_theory_variogram_plot(self):
        self.gaussian.variogram_plot()
        self.spherical.variogram_plot()
