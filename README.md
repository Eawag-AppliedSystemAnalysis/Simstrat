# Simstrat: a one-dimensional physical lake model

Simstrat is a one-dimensional physical lake model for the simulation of stratification and mixing in deep stratified lakes. The model was originally developed by Goudsmit et al. (2002) and has been successfully applied to lakes with different physical properties. A k-ε model is used to model turbulent mixing including energy transfer of internal seiches. River or groundwater inflow can be added at specific depths or as density-dependent intrusions. The newest version of Simstrat (see below) can also simulate ice/snow covers.

## Version 2.0 (2018)
**Architecture**
- Object-oriented, modern Fortran 2003 architecture
- JSON formatted configuration files
- Consistent array indexing in fortran standard (arrays start with index 1)
- Docker container to build the model

**Model update**
- Addition of an ice/snow model based on MyLake
- Two different kinds of inflows:
	- Fixed inflows
	- Inflows which vary with changing water level

**Documentation**
- Updated documentation including numerical scheme

## Run Simstrat
Pre-built binaries are available [here](prebuilt).

## Build Simstrat

### Build using Docker (or in a Linux-based environment)
A complete compile work flow using a docker container is now available
[here](misc/docker_build_env).

## Build manually

**System Requirements**

- [Python](https://www.python.org/) 2.7 or later
- [FoBiS.py](https://github.com/szaghi/FoBiS) 2.2.6 or later (available via GitHub or pip using pip install FoBiS.py)
- 2 compiler options:
	- [Intel Fortran (Intel Parallel Studio XE 2016)](https://software.intel.com/en-us/parallel-studio-xe/choose-download) (not free)
	- Gfortran 6.3 or later (free)

In principle, the manual installation is platform independent. Be aware that other programs on your computer might already use some version of Python and thus interfere with any new installation of Python.
