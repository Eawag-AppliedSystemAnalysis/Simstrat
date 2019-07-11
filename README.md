# Simstrat: a one-dimensional physical lake model

Simstrat is a one-dimensional physical lake model for the simulation of stratification and mixing in deep stratified lakes. The model was originally developed by Goudsmit et al. (2002) and has been successfully applied to lakes with different physical properties. A k-ε model is used to model turbulent mixing including energy transfer of internal seiches. River or groundwater inflow can be added at specific depths or as density-dependent intrusions. The newest version of Simstrat (see below) can also simulate ice/snow covers.

## Run Simstrat
Pre-built binaries are available [here](https://github.com/Eawag-AppliedSystemAnalysis/Simstrat/releases).

## Build Simstrat

### Linux or MacOS environment
We suggest to setup the building environment using Docker; a complete step-by-step guide to use a docker container is available
[here](misc/docker_build_env).

Once the building environment is setup, you can build simstrat with the provided build script. More details [here](build).

### Windows environment
Please install the required packages listed below (System requirements) and then from terminal navigate into `build` folder and run 

~~~bash
FoBiS.py build
~~~

**System requirements**

- [Python](https://www.python.org/) 2.7 or later
- [FoBiS.py](https://github.com/szaghi/FoBiS) 2.2.6 or later (available via GitHub or pip using `pip install FoBiS.py`)
- 2 compiler options:
	- [Intel Fortran (Intel Parallel Studio XE 2016)](https://software.intel.com/en-us/parallel-studio-xe/choose-download) (commercial)
	- Gfortran 6.3 or later (free)

In principle, the manual installation is platform independent. Be aware that other programs on your computer might already use some version of Python and thus interfere with any new installation of Python.

## Documentation

The user manual can be found [here](doc).

The developer documentation can be generated with the FORD python module (`pip install ford`).
To generate the documentation, run

~~~bash
FoBiS.py rule -ex makedoc
~~~

The generated code documentation is saved in `doc/developer/ford/ford_doc/index.htlm`

Additionally, a documentation about the numerical scheme of Simstrat can be found [here](doc/developer/dev_manual).