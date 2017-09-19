%%%%%%%%%%%%%%%%%% README FOR TESTCASES %%%%%%%%%%%%%%%%%%%

Usage: Compile Simstrat code (in source folder, using gfortran) and place 
the .exe file in this folder. Then execute the file with the chosen TestCase 
in the command window: 
simstrat.exe TestCase_x.par


Testcase 1: constant water table

- Data from Lake Zurich (Schmid and Köster, 2016) except for in- output
- Water input is constant in time and uniform over depth range of 0-2m
- Water output is constant in time and uniform over depth range of 0-1m
- Total water input always equals total water output (2 m^3/s)
- Input temperature is 20°C
- Input salinity is 0.1 %
- The simulation lasts for 2 years


Testcase 2: variable water table

- Data from Lake Zurich (Schmid and Köster, 2016) except for in- output
- Water input is variable in time but uniform over depth range of 0-2m
- Water output is constant in time and uniform over depth range of 0-2m
- Mean total water input equals mean total water output (2 m^3/s)
- Input temperature is 20°C
- Input salinity is 0.1 %
- The simulation lasts for 2 years

