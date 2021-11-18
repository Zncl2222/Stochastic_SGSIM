## Sequential Gaussian Simulation

This is a stochastic process based on the geostatistc algorithm. This project coded by several languages like C, C++ and pyhton. 

## Features
* Python version (UC_SGSIM_GUI.py)
  * One dimensional unconditional randomfield generation
  * Graphic user interface (GUI)
  * Now support Gaussian, Spherical, Exponential model
  * Enable to use muti-core to run the simulation (mutiprocessing)
  * Update and show the statistic results on the GUI 
 
* Ipython version (UC_SGSIM.ipython)
  * One dimensional unconditional randomfield generation
  * Now support Gaussian, Spherical, Exponential model
  * Enable to use muti-core to run the simulation (ipyparallel)
  * Plot figure on the notebook

* C version (UC_SGSIM.c)
  * One dimensional unconditional randomfield generation
  * Now support Gaussian, Spherical, Exponential model
  * Better computation efficiency


## Efficiency Comparison
<p align="center">
<img src="https://github.com/Zncl2222/Stochastic_SGSIM/blob/main/figure/C_Cpp_py_comparision.png"  width="70%"/>
</p>

```
Parameters for testing:

model len = 150

number of realizations = 1000

Range scale = 17.32

Variogram model = Gaussian model

---------------------------------------------------------------------------------------

Testing platform:

CPU: AMD Ryzen 9 4900 hs

RAM: DDR4 - 3200 40GB (Dual channel 16GB)

Disk: WD SN530
```
