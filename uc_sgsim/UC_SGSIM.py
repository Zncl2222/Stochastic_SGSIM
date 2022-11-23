import time
import sys
from pathlib import Path
from ctypes import CDLL, POINTER, c_double, c_int
from multiprocessing import Pool

import numpy as np
import uc_sgsim as UC
from uc_sgsim.Plot.Plot import Visualize

BASE_DIR = Path(__file__).resolve().parent


class Simulation:
    def __init__(self, Y, model, realization_number, randomseed=0, krige_method='SimpleKrige'):
        self.Y = Y
        self.model = model
        self.realization_number = realization_number
        self.bandwidth_step = model.bandwidth_step
        self.bandwidth = model.bandwidth
        self.randomseed = randomseed
        self.krige_method = krige_method
        self.size = len(self.Y)
        self.random_field = np.empty([self.realization_number, self.size])
        self.parallel_times = 0

    def compute(self, randomseed=0, parallel=False):
        if parallel is True:
            self.randomseed = randomseed
        else:
            self.n_process = 1

        initial_seed = self.randomseed

        if self.krige_method == 'SimpleKrige':
            self.krige = UC.SimpleKrige(self.model)

        counts = 0

        start_time = time.time()

        while counts < self.realization_number // self.n_process:
            boundary_constrained = 0
            unsampled = np.linspace(0, self.size - 1, self.size)

            if boundary_constrained == 0:
                y_value = np.random.normal(0, 1, 2).reshape(2, 1)
                x_grid = np.array([0, self.size - 1]).reshape(2, 1)
                Z = np.zeros(self.size)
                Z[0], Z[-1] = y_value[0], y_value[1]
                unsampled = np.delete(unsampled, [0, -1])
                neigh = 0
            else:
                y_value = np.random.normal(0, 1, 1).reshape(1, 1)
                ridx = np.random.randint(0, self.size - 1, 1)
                x_grid = np.array([ridx]).reshape(1, 1)
                Z = np.zeros(self.size)
                Z[ridx] = y_value[0]
                unsampled = np.delete(unsampled, [ridx])
                neigh = 1

            L = np.hstack([x_grid, y_value])

            np.random.seed(self.randomseed)
            randompath = np.random.choice(
                unsampled,
                len(unsampled),
                replace=False,
            )

            for i in range(len(unsampled)):
                Z[int(randompath[i])] = self.krige.compute(
                    L,
                    randompath[i],
                    neigh,
                    self.randomseed,
                )

                temp = np.hstack([randompath[i], Z[int(randompath[i])]])
                L = np.vstack([L, temp])

                if neigh < 6:
                    neigh += 1

                self.randomseed += 1

            Z_Gap = abs(Z.max() - Z.min())

            if 2 < Z_Gap <= 6.5:
                self.random_field[counts, :] = Z
                counts = counts + 1
                print('Progress = %.2f' % (counts / self.realization_number * 100) + '%', end='\r')

            self.randomseed += 1

        print('Progress = %.2f' % 100 + '%\n', end='\r')

        end_time = time.time()

        print('Time = %f' % (end_time - start_time), 's\n')
        print('Last RandomSeed = %d' % (self.randomseed), '\n')
        print('RandomSeed passed = %d' % (self.randomseed - initial_seed), '\n')

        return self.random_field

    def compute_async(self, n_process, randomseed):
        pool = Pool(processes=n_process)
        self.n_process = n_process
        self.realization_number = self.realization_number * n_process
        self.random_field = np.empty([self.realization_number, self.size])

        randomseed = []
        parallel = []
        for i in range(n_process):
            s = self.randomseed + int(i) * (self.realization_number + 300) * (self.size)
            randomseed.append(int(s))
            parallel.append(True)

        Z = pool.starmap(self.compute, zip(randomseed, parallel))

        for i in range(n_process):
            for j in range(int(self.realization_number / n_process)):
                start = int(i * self.realization_number / n_process)
                self.random_field[j + start, :] = Z[i][j, :]

        return self.random_field

    def variogram_compute(self, n_process=1):
        pool = Pool(processes=n_process)
        model_len = self.size
        x = np.linspace(0, self.size - 1, model_len).reshape(model_len, 1)

        L = []
        for i in range(self.realization_number):
            L.append(
                np.hstack([x, self.random_field[i, :].reshape(model_len, 1)]),
            )

        self.Variogram = pool.starmap(self.model.Variogram, zip(L))
        self.Variogram = np.array(self.Variogram).T

    def MeanPlot(self, n, mean=0, std=1):
        m_plot = Visualize(self.model, self.random_field)
        m_plot.MeanPlot(n, mean, std)

    def VarPlot(self, mean=0, std=1):
        s_plot = Visualize(self.model, self.random_field)
        s_plot.Variance_Plot(mean, std)

    def Cdf_Plot(self, x_location):
        c_plot = Visualize(self.model, self.random_field)
        c_plot.CDF_Plot(x_location)

    def Hist_Plot(self, x_location):
        h_plot = Visualize(self.model, self.random_field)
        h_plot.HIST(x_location)

    def VarioPlot(self):
        v_plot = Visualize(self.model, self.random_field)
        v_plot.Variogram_Plot(self.Variogram)

    def Save_random_field(self, path):
        for i in range(self.realization_number):
            if i < 10:
                number = '000' + str(i)
            elif 10 <= i < 100:
                number = '00' + str(i)
            elif 100 <= i < 1000:
                number = '0' + str(i)
            elif i >= 1000:
                number = str(i)

            with open(path + 'Realizations' + number + '.txt', 'w') as f:
                for j in range(0, self.size):
                    print(
                        '%.2d' % (j),
                        '%10.6f' % (self.random_field[j, i]),
                        file=f,
                    )

    def Save_Variogram(self, path):
        for i in range(self.realization_number):
            if i < 10:
                number = '000' + str(i)
            elif 10 <= i < 100:
                number = '00' + str(i)
            elif 100 <= i < 1000:
                number = '0' + str(i)
            elif i >= 1000:
                number = str(i)

            with open(path + 'Variogram' + number + '.txt', 'w') as f:
                for j in range(0, len(self.bandwidth_step)):
                    print(
                        '%.2d' % (j),
                        '%10.6f' % (self.random_field[j, i]),
                        file=f,
                    )


class Simulation_byC(Simulation):
    def __init__(self, Y, model, realization_number, randomseed=0, krige_method='SimpleKrige'):
        super().__init__(Y, model, realization_number, randomseed, krige_method)

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
        mlen = int(self.size)
        realization_number = int(self.realization_number // self.n_process)
        random_field = np.empty([realization_number, self.size])
        array = (c_double * (mlen * realization_number))()

        sgsim(array, mlen, realization_number, 17.32, 1, randomseed)

        for i in range(realization_number):
            random_field[i, :] = list(array)[i * mlen : (i + 1) * mlen]
        return random_field

    def compute_by_dll(self, n_process, randomseed):
        pool = Pool(processes=n_process)
        self.n_process = n_process

        if self.parallel_times < 1:
            self.realization_number = self.realization_number * n_process
            self.parallel_times += 1
        else:
            self.realization_number = self.realization_number

        self.random_field = np.empty([self.realization_number, self.size])

        randomseed = []
        for i in range(n_process):
            s = self.randomseed + int(i) * (self.realization_number + 300) * (self.size)
            randomseed.append(int(s))

        Z = pool.starmap(self.cpdll, zip(randomseed))

        for i in range(n_process):
            for j in range(int(self.realization_number / n_process)):
                start = int(i * self.realization_number / n_process)
                self.random_field[j + start, :] = Z[i][j, :]

        return self.random_field

    def vario_cpdll(self, cpu_number):
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

        mlen = int(self.size)
        realization_number = int(self.realization_number // self.n_process)

        Variogram = np.empty([self.size, realization_number])
        vario_size = len(self.bandwidth_step)

        vario_array = (c_double * (vario_size))()
        random_field_array = (c_double * (mlen))()

        Variogram = np.empty([vario_size, realization_number])

        for i in range(realization_number):
            random_field_array[:] = self.random_field[:, i + cpu_number * realization_number]
            vario(random_field_array, vario_array, mlen, vario_size, 1)
            Variogram[:, i] = list(vario_array)

        return Variogram

    def vario_compute_by_dll(self, n_process=1):
        pool = Pool(processes=n_process)
        self.n_process = n_process

        if self.parallel_times < 1:
            self.realization_number = self.realization_number * n_process
            self.parallel_times += 1
        else:
            self.realization_number = self.realization_number

        self.Variogram = np.empty([len(self.bandwidth_step), self.realization_number])
        cpu_number = []
        for i in range(self.n_process):
            cpu_number.append(i)

        Z = pool.starmap(self.vario_cpdll, zip(cpu_number))

        for i in range(n_process):
            for j in range(int(self.realization_number / n_process)):
                self.Variogram[:, (j + int(i * self.realization_number / n_process))] = Z[i][:, j]
        return self.Variogram
