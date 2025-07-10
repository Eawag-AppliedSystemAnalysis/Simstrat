# Simstrat: a one-dimensional physical lake model coupled to a biogeochemical library

Simstrat is a one-dimensional physical lake model for the simulation of stratification and mixing in deep stratified lakes. The model was originally developed by Goudsmit et al. (2002) and has been successfully applied to lakes with different physical properties. A k-ε model is used to model turbulent mixing including energy transfer of internal seiches. River or groundwater inflow can be added at specific depths or as density-dependent intrusions. The newest version of Simstrat can also simulate ice/snow covers and do biogeochemical simulations by using its coupling with FABM (https://github.com/fabm-model/fabm). From version 3.0 onwards, Simstrat always includes a coupling to FABM but the coupling can be turned off using a switch in the configuration file for purely physical simulations.

## Run Simstrat

There are a number of options for running Simstrat: 

- Pre-built binaries are available [here](https://github.com/Eawag-AppliedSystemAnalysis/Simstrat/releases).
- Build Simstrat
- [LakeEnsemblR](https://github.com/aemon-j/LakeEnsemblR)
- Docker container

### Pre-build binary

1. Download the appropriate binary for your operating system and version requirements.
2. Make the binary executable `chmod +x simstrat_linux_303`
3. Pass the location of the par file to the binary
~~~bash
/pathtobinary/simstrat_linux_303 /pathtoparfile/test.par
~~~

### Docker container

To use the docker container all input files, including the par file must be in the same folder. The file paths in the 
par file must be relative to the parent folder that is mounted to the docker container.

Download the docker container

`docker pull eawag/simstrat:3.0.3`

Run the docker container for `test.par` where all input files including `test.par` are located in `/pathtoinputfiles`.
```bash
cd /pathtoinputfiles
docker run -v $(pwd):/simstrat/run eawag/simstrat:3.0.3 test.par
```


## Build Simstrat
After cloning the git repository for the first time,  you need to initialize the 3rd-party libraries with the given script
~~~bash
./git_lib_initialize.sh
~~~
If you are not using git, you can manually download the library source files from the addresses listed in the `.gitmodules` file and save them into the `lib` subfolder.


Before building an executable of Simstrat, you need to setup your building environment. We suggest two alternative options:

### 1a. Setup building environment using Docker
You can setup the building environment using Docker for Linux, MacOS and Win hosting systems; a complete step-by-step guide to use a docker container is available
[here](misc/docker_build_env).


### 1b. Manual setup of the building environment
Please install the following required packages:

**System requirements**

- [Python](https://www.python.org/) 2.7 or later
- [FoBiS.py](https://github.com/szaghi/FoBiS) 3.0.1 or later (available via GitHub or pip using `pip install FoBiS.py`)
- 2 compiler options:
    - [Intel Fortran (Intel Parallel Studio XE 2016)](https://software.intel.com/en-us/parallel-studio-xe/choose-download) (commercial)
    - Gfortran 6.3 or later (free)

In principle, the manual installation is platform independent. Be aware that other programs on your computer might already use some version of Python and thus interfere with any new installation of Python.

### 2. Build FABM

Go to `lib/fabm` and run:

~~~bash
cmake -S ./ -B build -DFABM_HOST=simstrat
~~~

to create the new directory `build` and generate the build configuration inside it, with the Simstrat preprocesser definitions for information about the spatial domain.

> **N.B.1** By default FABM includes all biogeochemical models from `lib/fabm/src/models`. You can restrict this to a subset of directories by providing  
> ~~~bash
> -DFABM_INSTITUTES=<INSTITUTE-NAMES>
> ~~~  
> as an argument to the `cmake` command above. To add externally maintained models add the argument  
> ~~~bash
> -DFABM_EXTRA_INSTITUTES=<INSTITUTES-NAMES>
> ~~~  
> to the `cmake` command. Here, `<INSTITUTE-NAMES>` is a semi-colon-separated list of institute names (all in lower case). You may need to enclose this list in quotes to prevent the shell from interpreting the semi-colon. For the externally maintained models, you additionally need to point FABM to the directory with the source code by adding the argument  
> ~~~bash
> -DFABM_<INSTITUTE-NAME>_BASE=<DIR>
> ~~~  
> to the `cmake` command for every `<INSTITUTE-NAME>` (here all in upper case) in `<INSTITUTE-NAMES>`. With `<DIR>` the corresponding directory on your device. A list of available externally maintained models can be found here: [Available Models in FABM Wiki](https://github.com/fabm-model/fabm/wiki/Biogeochemical-models-in-FABM#available-models). For instance to build FABM with [the Selmaprotbas model](https://github.com/jorritmesman/Selmaprotbas) the full cmake command is
> ~~~bash
> cmake -S ./ -B build -DFABM_HOST=simstrat -DFABM_EXTRA_INSTITUTES=selmaprotbas -DFABM_SELMAPROTBAS_BASE=<DIR>
> ~~~ 
> with `<DIR>` the directory of the selmaprotbas repository on your device.

> **N.B.2** Add  
> ~~~bash
> -DCMAKE_BUILD_TYPE=Debug
> ~~~  
> as argument to `cmake` to compile FABM in debug mode.

Then run:

~~~bash
cmake --build build --target install
~~~

to build and install FABM.

### 3. Build Simstrat

Now leave the fabm folder and go to `Simstrat/build` folder with:
~~~bash
cd ../../build
~~~

from where you can build Simstrat with:

~~~bash
FoBiS.py build
~~~

> **N.B.1** macOS seems to allow only dynamic library linking. To build for this target, use the `release-gnu-dynamic` compiling mode  
> ~~~bash
> FoBiS.py build -mode release-gnu-dynamic
> ~~~

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
