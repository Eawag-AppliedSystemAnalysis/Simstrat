# Simstrat: a one-dimensional physical lake model coupled to a biogeochemical library

Simstrat is a one-dimensional physical lake model for the simulation of stratification and mixing in deep stratified lakes. The model was originally developed by Goudsmit et al. (2002) and has been successfully applied to lakes with different physical properties. A k-ε model is used to model turbulent mixing including energy transfer of internal seiches. River or groundwater inflow can be added at specific depths or as density-dependent intrusions. The newest version of Simstrat can also simulate ice/snow covers and do biogeochemical simulations by using its coupling with AED2 (https://github.com/AquaticEcoDynamics/libaed2). From version 3.0 onwards, Simstrat always includes a coupling to AED2 but the coupling can be turned off using a switch in the configuration file for purely physical simulations.

## Run Simstrat
Pre-built binaries are available [here](https://github.com/Eawag-AppliedSystemAnalysis/Simstrat/releases).

## Build Simstrat
After cloning the git repository for the first time,  you need to initialize the 3rd-party libraries with the given script
~~~bash
./git_lib_initialize.sh
~~~
If you are not using git, you can manually download the library source files from the addresses listed in the `.gitmodules` file and save them into the `lib` subfolder.


Before building an executable of Simstrat, you need to setup your building environment. We suggest two alternative options:

### 1. Setup building environment using Docker
You can setup the building environment using Docker for Linux, MacOS and Win hosting systems; a complete step-by-step guide to use a docker container is available
[here](misc/docker_build_env).


### 2. Manual setup of the building environment
Please install the following required packages:

**System requirements**

- [Python](https://www.python.org/) 2.7 or later
- [FoBiS.py](https://github.com/szaghi/FoBiS) 3.0.1 or later (available via GitHub or pip using `pip install FoBiS.py`)
- 2 compiler options:
    - [Intel Fortran (Intel Parallel Studio XE 2016)](https://software.intel.com/en-us/parallel-studio-xe/choose-download) (commercial)
    - Gfortran 6.3 or later (free)

In principle, the manual installation is platform independent. Be aware that other programs on your computer might already use some version of Python and thus interfere with any new installation of Python.

### Build
Once the building environment is setup, you can build Simstrat from the `/build` folder with:
~~~bash
FoBiS.py build
~~~

If the AED2 library is not built yet, use make (or mingw32-make in windows) to compile it in the lib/libaed2 folder.

> **N.B.1** MacOS seems to allow only dynamic library linking. To build for this target, use the `release-gnu-dynamic` compiling mode
~~~bash
FoBiS.py build -mode release-gnu-dynamic
~~~

> **N.B.2** you can find more building options [here](build).




## Documentation

The user manual can be found [here](doc).

The developer documentation can be generated with the FORD python module (`pip install ford`).
To generate the documentation, run

~~~bash
FoBiS.py rule -ex makedoc
~~~

The generated code documentation is saved in `doc/developer/ford/ford_doc/index.html`

Additionally, a documentation about the numerical scheme of Simstrat can be found [here](doc/developer/dev_manual).
