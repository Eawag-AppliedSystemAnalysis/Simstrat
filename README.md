# Simstrat

The source code of different simstrat versions are available here:

[Version 1.0](https://github.com/Eawag-AppliedSystemAnalysis/Simstrat/releases/tag/V1.0)  
Original model (Goudsmit, 2002)  
Year: 2000  

[Version 1.2](https://github.com/Eawag-AppliedSystemAnalysis/Simstrat/releases/tag/V1.2)  
Revision of the code, bug fixes  
Year: 2014  

[Version 1.4](https://github.com/Eawag-AppliedSystemAnalysis/Simstrat/releases/tag/V1.4)  
Improved parameterization of heat fluxes, revision of the code (Schmid and Köster, 2016)  
Year: 2016  

## Compile Simstrat

1. Download the simstrat source code of a specific simstrat version using the link above or clone the current master branch:
```
git clone https://github.com/Eawag-AppliedSystemAnalysis/Simstrat
```

2. Navigate to the 'source' folder within the repository.

3. Depending on the compiler on your system run the following command:

   Intel Fortran Compile  
   ```
   ifort simstrat.f90 -O3 -Qipo -Qprec-div -QxHost
   ```

   gfortran
   ```
   gfortran simstrat.f90 -o simstrat.exe -O2 -ffree-line-length-none -g -ffpe-trap=overflow,zero,invalid -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans
   ```