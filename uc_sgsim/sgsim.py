import time
import sys
from pathlib import Path
from ctypes import CDLL, POINTER, c_double, c_int
from multiprocessing import Pool

import numpy as np
from uc_sgsim.exception import VariogramDoesNotCompute
from uc_sgsim.krige import SimpleKrige
from uc_sgsim.random_field import RandomField
from uc_sgsim.plot.plot import Visualize

BASE_DIR = Path(__file__).resolve().parent


class UCSgsim(RandomField):
    def __init__(self, x, model, realization_number, krige_method='SimpleKrige'):
        super().__init__(x, model, realization_number)
        self.krige_method = krige_method

    def compute(self, randomseed=0, parallel=False):
        self.randomseed = randomseed
        if parallel is False:
            self.n_process = 1

        if self.krige_method == 'SimpleKrige':
            self.krige = SimpleKrige(self.model)

        counts = 0

        start_time = time.time()
        np.random.seed(self.randomseed)
        while counts < self.realization_number // self.n_process:
            boundary_constrained = 0
            unsampled = np.linspace(0, self.x_size - 1, self.x_size)

            if boundary_constrained == 0:
                y_value = np.random.normal(0, 1, 2).reshape(2, 1)
                x_grid = np.array([0, self.x_size - 1]).reshape(2, 1)
                z = np.zeros(self.x_size)
                z[0], z[-1] = y_value[0], y_value[1]
                unsampled = np.delete(unsampled, [0, -1])
                neigh = 0
            else:
                y_value = np.random.normal(0, 1, 1).reshape(1, 1)
                ridx = np.random.randint(0, self.x_size - 1, 1)
                x_grid = np.array([ridx]).reshape(1, 1)
                z = np.zeros(self.x_size)
                z[ridx] = y_value[0]
                unsampled = np.delete(unsampled, [ridx])
                neigh = 1

            L = np.hstack([x_grid, y_value])

            randompath = np.random.choice(
                unsampled,
                len(unsampled),
                replace=False,
            )

            for i in range(len(unsampled)):
                z[int(randompath[i])] = self.krige.simulation(
                    L,
                    randompath[i],
                    neighbor=neigh,
                )
                temp = np.hstack([randompath[i], z[int(randompath[i])]])
                L = np.vstack([L, temp])

                if neigh < 8:
                    neigh += 1

            self.randomseed += 1

            z_gap = abs(z.max() - z.min())

            if 2 < z_gap <= 6.5:
                self.random_field[counts, :] = z
                counts = counts + 1
                print('Progress = %.2f' % (counts / self.realization_number * 100) + '%', end='\r')

        print('Progress = %.2f' % 100 + '%\n', end='\r')
        end_time = time.time()
        print('Time = %f' % (end_time - start_time), 's\n')

        return self.random_field

    def compute_async(self, n_process, randomseed):
        pool = Pool(processes=n_process)
        self.n_process = n_process
        self.realization_number = self.realization_number * n_process
        self.random_field = np.empty([self.realization_number, self.x_size])

        rand_list = []
        parallel = []
        for i in range(n_process):
            s = randomseed + int(i)
            rand_list.append(int(s))
            parallel.append(True)

        z = pool.starmap(self.compute, zip(rand_list, parallel))

        for i in range(n_process):
            for j in range(int(self.realization_number / n_process)):
                start = int(i * self.realization_number / n_process)
                self.random_field[j + start, :] = z[i][j, :]

        return self.random_field

    def variogram_compute(self, n_process=1):
        pool = Pool(processes=n_process)
        model_len = self.x_size
        x = np.linspace(0, self.x_size - 1, model_len).reshape(model_len, 1)

        L = []
        for i in range(self.realization_number):
            L.append(
                np.hstack([x, self.random_field[i, :].reshape(model_len, 1)]),
            )

        self.variogram = pool.starmap(self.model.variogram, zip(L))
        self.variogram = np.array(self.variogram)

    def mean_plot(self, n, mean=0, std=1):
        m_plot = Visualize(self.model, self.random_field)
        m_plot.mean_plot(n, mean, std)

    def variance_plot(self, mean=0, std=1):
        s_plot = Visualize(self.model, self.random_field)
        s_plot.variance_plot(mean, std)

    def cdf_plot(self, x_location):
        c_plot = Visualize(self.model, self.random_field)
        c_plot.cdf_plot(x_location)

    def hist_plot(self, x_location):
        h_plot = Visualize(self.model, self.random_field)
        h_plot.hist_plot(x_location)

    def vario_plot(self):
        if type(self.variogram) == int:
            raise VariogramDoesNotCompute()
        v_plot = Visualize(self.model, self.random_field)
        v_plot.variogram_plot(self.variogram)


class UCSgsimDLL(UCSgsim):
    def __init__(self, Y, model, realization_number, krige_method='SimpleKrige'):
        super().__init__(Y, model, realization_number)
        self.krige_method = krige_method

    def lib_read(self):
        if sys.platform.startswith('linux'):
            lib = CDLL(str(BASE_DIR) + r'/c_core/uc_sgsim.so')
        elif sys.platform.startswith('win32'):
            lib = CDLL(str(BASE_DIR) + r'/c_core/uc_sgsim.dll')
        return lib

    def cpdll(self, randomseed):
        lib = self.lib_read()
        sgsim = lib.sgsim_dll
        sgsim.argtypes = (
            POINTER(c_double),
            c_int,
            c_int,
            c_double,
            c_double,
            c_int,
        )
        sgsim.restype = None
        mlen = int(self.x_size)
        realization_number = int(self.realization_number // self.n_process)
        random_field = np.empty([realization_number, self.x_size])
        array = (c_double * (mlen * realization_number))()

        sgsim(array, mlen, realization_number, 17.32, 1, randomseed)

        for i in range(realization_number):
            random_field[i, :] = list(array)[i * mlen : (i + 1) * mlen]
        return random_field

    def compute(self, n_process, randomseed):
        pool = Pool(processes=n_process)
        self.n_process = n_process

        if n_process > 1:
            self.realization_number = self.realization_number * n_process
        else:
            self.realization_number = self.realization_number

        self.random_field = np.empty([self.realization_number, self.x_size])

        rand_list = []
        for i in range(n_process):
            s = randomseed + int(i)
            rand_list.append(int(s))

        z = pool.starmap(self.cpdll, zip(rand_list))

        for i in range(n_process):
            for j in range(int(self.realization_number / n_process)):
                start = int(i * self.realization_number / n_process)
                self.random_field[j + start, :] = z[i][j, :]

        return self.random_field

    def variogram_cpdll(self, cpu_number):
        lib = self.lib_read()
        vario = lib.variogram
        vario.argtypes = (
            POINTER(c_double),
            POINTER(c_double),
            c_int,
            c_int,
            c_int,
        )
        vario.restype = None

        mlen = int(self.x_size)
        realization_number = int(self.realization_number // self.n_process)

        vario_size = len(self.bandwidth_step)

        vario_array = (c_double * (vario_size))()
        random_field_array = (c_double * (mlen))()

        Variogram = np.empty([realization_number, vario_size])

        for i in range(realization_number):
            random_field_array[:] = self.random_field[i + cpu_number * realization_number, :]
            vario(random_field_array, vario_array, mlen, vario_size, 1)
            Variogram[i, :] = list(vario_array)

        return Variogram

    def variogram_compute(self, n_process=1):
        pool = Pool(processes=n_process)
        self.n_process = n_process

        if n_process < 1:
            self.realization_number = self.realization_number * n_process
        else:
            self.realization_number = self.realization_number

        self.variogram = np.empty([self.realization_number, len(self.bandwidth_step)])
        cpu_number = []
        for i in range(self.n_process):
            cpu_number.append(i)

        z = pool.starmap(self.variogram_cpdll, zip(cpu_number))

        for i in range(n_process):
            for j in range(int(self.realization_number / n_process)):
                self.variogram[(j + int(i * self.realization_number / n_process)), :] = z[i][j, :]
        return self.variogram
