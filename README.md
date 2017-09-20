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

   Intel Fortran Compiler (on Windows)
   ```
   ifort simstrat.f90 /O3 /Qipo /Qprec-div /QxHost
   ```

   Intel Fortran Compiler (on Linux and macOS)
   ```
   ifort simstrat.f90 -O3 -ipo -prec-div -xHost
   ```

   gfortran (on Windows)
   ```
   gfortran simstrat.f90 -o simstrat.exe -O2 -ffree-line-length-none -g -ffpe-trap=overflow,zero,invalid -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -static-libgfortran
   ```

   gfortran (on Linux and macOS)
   ```
   gfortran simstrat.f90 -o simstrat -O2 -ffree-line-length-none -g -ffpe-trap=overflow,zero,invalid -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -static-libgfortran
   ```
   (don't forget to make the file executable, i.e. `chmod +x simstrat`)

## Getting gfortran for Windows

1. To get gfortran on Windows we recommend to download and install [MSYS2](http://www.msys2.org/). Follow the installation instruction on the webpage (update the package database with `pacman -Syu`).

2. Install a toolchain:
   a) for 32-bit:
      ```
      pacman -S mingw-w64-i686-toolchain
      ```
   b) for 64-bit:
      ```
      pacman -S mingw-w64-x86_64-toolchain
      ```
   Select which package to install, default is all

3. To compile simstrat use the 'MinGW-w64 32-bit Shell' or 'MinGW-w64 64-bit Shell' which should be available in your start menu now.

## Getting Intel Fortran Compiler for Windws

Intel supports the next generation of scientists and engineers by ensuring that students and the faculty at qualifying universities have access to the full-featured Intel® Parallel Studio XE.
You can find the free options for Windows, Linux and macOS [here](https://software.intel.com/en-us/parallel-studio-xe/choose-download).
