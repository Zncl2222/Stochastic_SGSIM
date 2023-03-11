import unittest
import pytest

import numpy as np

import uc_sgsim as uc
from uc_sgsim.cov_model import Gaussian, Spherical
from uc_sgsim.exception import VariogramDoesNotCompute


@pytest.mark.sgsim
class UCSgsimTest(unittest.TestCase):
    def setUp(self) -> None:
        self.X = 151
        self.nR = 10
        bw = 1
        hs = np.arange(0.0, 35, bw)
        a = 17.32
        C0 = 1
        self.gaussian = Gaussian(hs, bw, a, C0)
        self.spherical = Spherical(hs, bw, a, C0)

    def sgsim_plot(self, sgsim: uc.UCSgsim) -> None:
        sgsim.mean_plot('ALL')
        sgsim.mean_plot([0, 1, 2])
        sgsim.cdf_plot(x_location=10)
        sgsim.hist_plot(x_location=10)
        sgsim.variogram_compute(n_process=1)
        sgsim.vario_plot()
        sgsim.variance_plot()

    def sgsim_save(self, sgsim: uc.UCSgsim) -> None:
        sgsim.save_random_field('randomfield.csv', save_single=True)
        sgsim.save_variogram('variogram.csv', save_single=True)
        sgsim.save_random_field('randomfield/', save_single=False)
        sgsim.save_variogram('variogram/', save_single=False)

    def test_uc_sgsim_gaussian_py_single_process(self):
        sgsim = uc.UCSgsim(self.X, self.gaussian, self.nR)
        sgsim.compute(randomseed=454)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_py_multi_process(self):
        sgsim = uc.UCSgsim(self.X, self.gaussian, self.nR)
        sgsim.compute_async(n_process=2, randomseed=454)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_spherical_py(self):
        sgsim = uc.UCSgsim(self.X, self.spherical, self.nR)
        sgsim.compute_async(n_process=2, randomseed=454)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_c(self):
        sgsim = uc.UCSgsimDLL(self.X, self.gaussian, self.nR)
        sgsim.compute(n_process=2, randomseed=454)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_varigram_exception(self):
        sgsim = uc.UCSgsimDLL(self.X, self.gaussian, self.nR)
        sgsim.compute(n_process=1, randomseed=454)
        self.assertRaises(VariogramDoesNotCompute, lambda: sgsim.vario_plot())

        sgsim = uc.UCSgsim(self.X, self.gaussian, self.nR)
        sgsim.compute_async(n_process=1, randomseed=454)
        self.assertRaises(VariogramDoesNotCompute, lambda: sgsim.vario_plot())
        self.assertRaises(
            VariogramDoesNotCompute,
            lambda: sgsim.save_variogram('variogram/', save_single=False),
        )
