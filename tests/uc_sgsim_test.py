import pytest

import uc_sgsim as uc
from uc_sgsim.kriging import SimpleKriging
from uc_sgsim.cov_model import Gaussian, Spherical, Exponential
from uc_sgsim.exception import VariogramDoesNotCompute, IterationError
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

    def test_uc_sgsim_gaussian_py_with_params2(self):
        sgsim = uc.UCSgsim(
            self.X,
            self.nR,
            self.gaussian,
            z_max=3,
            z_min=-3,
            max_neighbor=12,
            constant_path=True,
            cov_cache=True,
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

    def test_uc_sgsim_gaussian_py_o_kriging_single_process_with_kwargs(self):
        sgsim = uc.UCSgsim(
            self.X,
            self.nR,
            self.gaussian,
            kriging='OrdinaryKriging',
            constant_path=True,
            cov_cache=True,
        )
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

    def test_get_all_attributes(self):
        sgsim = uc.UCSgsim(
            self.X,
            self.nR,
            self.gaussian,
            z_max=3,
            z_min=-3,
            max_neighbor=12,
            constant_path=True,
            cov_cache=True,
        )
        attributes = sgsim.get_all_attributes()
        assert attributes['grid_size'] == 151
        assert attributes['realization_number'] == 10
        assert attributes['z_max'] == 3
        assert attributes['z_min'] == -3
        assert attributes['max_neighbor'] == 12
        assert attributes['constant_path'] is True
        assert attributes['cov_cache'] is True
        assert isinstance(attributes['kriging'], SimpleKriging)
        assert attributes['model'] == self.gaussian

    def test_sgsim_repr(self, capsys):
        sgsim = uc.UCSgsimDLL(self.X, self.nR, self.gaussian)
        # Call repr() to get the string representation
        repr_output = repr(sgsim)
        # Print the repr output (optional)
        print(repr_output)

        # Capture the printed output
        captured = capsys.readouterr()

        # Assert the printed output
        assert (
            captured.out.strip()
            == 'SgsimField(151, 10, SimpleKriging, GaussianModel'
            + '(bandwidth_len=35, bandwidth_step=1, k_range=17.32, sill=1, nugget=0))'
        )
        assert captured.err == ''
        assert (
            repr_output
            == 'SgsimField(151, 10, SimpleKriging, GaussianModel'
            + '(bandwidth_len=35, bandwidth_step=1, k_range=17.32, sill=1, nugget=0))'
        )

    def test_sgsim_reach_iteration_limit(self):
        sgsim = uc.UCSgsim(
            self.X,
            100,
            self.gaussian,
            z_max=1,
            z_min=-1,
            iteration_limit=3,
            max_neighbor=12,
        )
        with pytest.raises(IterationError):
            sgsim.compute(n_process=1, randomseed=12312)
