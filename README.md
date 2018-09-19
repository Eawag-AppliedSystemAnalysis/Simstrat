# Version 2 (2018)
**Architecture**
- Object-oriented, modern Fortran 2003 Architecture
- JSON formatted configuration files
- Consistent array indexing in fortran standard (arrays start with index 1)
- Docker container to run model

**Model update**
- Addition of an ice/snow model based on MyLake
- Two different kinds of inflows:
	- Fixed inflows
	- Inflows which vary with changing water level

**Documentation**
- Updated documentation including numerical scheme



## Compile Simstrat

**System Requirements**

- [Python](https://www.python.org/) 2.7 or later
- [FoBiS.py](https://github.com/szaghi/FoBiS) 2.2.6 or later (available via GitHub or pip)
- 2 compiler options:
	- [Intel Fortran (Intel Parallel Studio XE 2016)](https://software.intel.com/en-us/parallel-studio-xe/choose-download)
	- Gfortran 6.3 or later


**Building Simstrat v2.0 with FoBiS.py**

From version 2.0 onwards, [FoBiS.py](https://github.com/szaghi/FoBiS) is used as an automated build-system. The FoBiS.py python-package is available on pip:

```
pip install FoBiS.py
```

To compile simstrat 2.0 (beta) navigate to the Simstrat_v2 folder of the repository and run the following command:

- Using Gnu Fortran Compiler

  ```
  FoBiS.py build -mode release-gnu
  ```

- Using Intel Fortran Compiler (run the command in the parallel studio CMD!)

  ```
  FoBiS.py build -mode release-intel
  ```