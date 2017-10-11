# Version 2 (2017)
- Based on the modularized version
- Rebuild with array indexing in fortran standard (1-indexes whenever possible)
- Objectoriented Architecture
- JSON Files for configuration

## Compile Simstrat

From version 2.0 onwards, [FoBiS.py](https://github.com/szaghi/FoBiS) is used as an automated build-system. To get FoBiS.py install python 3.6 and install the FoBiS.py package with pip:

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