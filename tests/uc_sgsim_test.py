import pytest

import uc_sgsim as uc
from uc_sgsim.cov_model import Gaussian, Spherical, Exponential
from uc_sgsim.exception import VariogramDoesNotCompute
import os


@pytest.mark.sgsim
class TestUCSgsim:
    @classmethod
    def setup_class(cls) -> None:
        cls.X = 151
        cls.nR = 10
        cls.bw = 1
        cls.hs = 35
        cls.a = 17.32
        cls.C0 = 1
        cls.gaussian = Gaussian(cls.hs, cls.bw, cls.a, cls.C0)
        cls.spherical = Spherical(cls.hs, cls.bw, cls.a, cls.C0)
        cls.exponential = Exponential(cls.hs, cls.bw, cls.a, cls.C0)

    def sgsim_plot(self, sgsim: uc.UCSgsim) -> None:
        sgsim.plot()
        sgsim.plot([0, 1, 2])
        sgsim.save_plot()
        assert os.path.isfile('figure.png') is True
        sgsim.save_plot('test.jpg')
        assert os.path.isfile('test.jpg') is True
        sgsim.cdf_plot(x_location=10)
        sgsim.hist_plot(x_location=10)
        sgsim.variogram_compute(n_process=1)
        sgsim.variogram_plot()
        sgsim.mean_plot()
        sgsim.variance_plot()

    def sgsim_save(self, sgsim: uc.UCSgsim) -> None:
        sgsim.save_random_field('randomfield.csv', save_single=True)
        sgsim.save_variogram('variogram.csv', save_single=True)
        sgsim.save_random_field('randomfield/', save_single=False)
        sgsim.save_variogram('variogram/', save_single=False)

    def test_uc_sgsim_gaussian_py_single_process(self):
        sgsim = uc.UCSgsim(self.X, self.nR, self.gaussian)
        sgsim.compute(n_process=1, randomseed=454)
        sgsim.variogram_compute(n_process=1)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_py_with_params(self):
        sgsim = uc.UCSgsim(
            self.X,
            self.nR,
            self.gaussian,
            z_max=3,
            z_min=-3,
            max_neigh=6,
        )
        sgsim.compute(n_process=1, randomseed=454)
        sgsim.variogram_compute(n_process=1)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_py_multi_process(self):
        sgsim = uc.UCSgsim(self.X, self.nR, self.gaussian)
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_nugget_py(self):
        model = Gaussian(35, 1, 17.32, 1, 0.01)
        sgsim = uc.UCSgsim(self.X, self.nR, model)
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_py_o_kriging_single_process(self):
        sgsim = uc.UCSgsim(self.X, self.nR, self.gaussian, kriging='OrdinaryKriging')
        sgsim.compute(n_process=1, randomseed=454)
        sgsim.variogram_compute(n_process=1)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_py_o_kriging_multi_process(self):
        sgsim = uc.UCSgsim(self.X, self.nR, self.gaussian, kriging='OrdinaryKriging')
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_o_kriging_nugget_py(self):
        model = Gaussian(35, 1, 17.32, 1, 0.01)
        sgsim = uc.UCSgsim(self.X, self.nR, model, kriging='OrdinaryKriging')
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_spherical_py(self):
        sgsim = uc.UCSgsim(self.X, self.nR, self.spherical)
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_spherical_nugget_py(self):
        model = Spherical(35, 1, 17.32, 1, 0.01)
        sgsim = uc.UCSgsim(self.X, self.nR, model)
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_exponential_py(self):
        sgsim = uc.UCSgsim(self.X, self.nR, self.exponential)
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_exponential_nugget_py(self):
        model = Exponential(35, 1, 17.32, 1, 0.01)
        sgsim = uc.UCSgsim(self.X, self.nR, model)
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_c(self):
        sgsim = uc.UCSgsimDLL(self.X, self.nR, self.gaussian)
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_ordinary_kriging_c(self):
        sgsim = uc.UCSgsimDLL(self.X, self.nR, self.gaussian, kriging='OrdinaryKriging')
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_wrong_kriging_method(self):
        with pytest.raises(TypeError):
            sgsim = uc.UCSgsim(  # noqa
                self.X,
                self.nR,
                self.gaussian,
                kriging='wrongKriging',
            )

    def test_uc_sgsim_variogram_exception(self):
        sgsim = uc.UCSgsimDLL(self.X, self.nR, self.gaussian)
        sgsim.compute(n_process=1, randomseed=454)
        with pytest.raises(VariogramDoesNotCompute):
            sgsim.variogram_plot()

        sgsim = uc.UCSgsim(self.X, self.nR, self.gaussian)
        sgsim.compute(n_process=1, randomseed=454)
        with pytest.raises(VariogramDoesNotCompute):
            sgsim.variogram_plot()
        with pytest.raises(VariogramDoesNotCompute):
            sgsim.save_variogram('variogram/', save_single=False)

    def test_sgsim_plot_params(self):
        import numpy as np

        sgsim = uc.UCSgsim(self.X, self.nR, self.gaussian)
        vbandwidth = np.arange(0, self.hs, self.bw)
        np.testing.assert_allclose(sgsim.vbandwidth, vbandwidth)
        assert sgsim.vmodel == self.gaussian
        assert sgsim.vmodel_name == 'Gaussian'
        assert sgsim.vbandwidth_step == 1
        assert sgsim.vsill == 1
        assert sgsim.vk_range == 17.32
