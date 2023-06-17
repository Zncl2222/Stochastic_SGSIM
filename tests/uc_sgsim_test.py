import pytest

import uc_sgsim as uc
from uc_sgsim.cov_model import Gaussian, Spherical, Exponential
from uc_sgsim.exception import VariogramDoesNotCompute


@pytest.mark.sgsim
class TestUCSgsim:
    @classmethod
    def setup_class(cls) -> None:
        cls.X = 151
        cls.nR = 10
        bw = 1
        hs = 35
        a = 17.32
        C0 = 1
        cls.gaussian = Gaussian(hs, bw, a, C0)
        cls.spherical = Spherical(hs, bw, a, C0)
        cls.exponential = Exponential(hs, bw, a, C0)

    def sgsim_plot(self, sgsim: uc.UCSgsim) -> None:
        sgsim.plot()
        sgsim.plot([0, 1, 2])
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
        sgsim = uc.UCSgsim(self.X, self.nR, self.gaussian, krige_method='OrdinaryKrige')
        sgsim.compute(n_process=1, randomseed=454)
        sgsim.variogram_compute(n_process=1)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_py_o_kriging_multi_process(self):
        sgsim = uc.UCSgsim(self.X, self.nR, self.gaussian, krige_method='OrdinaryKrige')
        sgsim.compute(n_process=2, randomseed=454)
        sgsim.variogram_compute(n_process=2)
        self.sgsim_plot(sgsim)
        self.sgsim_save(sgsim)

    def test_uc_sgsim_gaussian_o_kriging_nugget_py(self):
        model = Gaussian(35, 1, 17.32, 1, 0.01)
        sgsim = uc.UCSgsim(self.X, self.nR, model, krige_method='OrdinaryKrige')
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

    def test_uc_sgsim_varigram_exception(self):
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
