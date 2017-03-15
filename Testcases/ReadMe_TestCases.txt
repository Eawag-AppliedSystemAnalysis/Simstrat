%%%%%%%%%%%%%%%%%% README FOR TESTCASES %%%%%%%%%%%%%%%%%%%

Usage: Compile Simstrat code (in source folder, using gfortran) and place 
the .exe file in this folder. Then execute the file with the chosen TestCase 
in the command window: 
simstrat.exe TestCase_x.par


Testcase 1: constant water table

- Data from Lake Zurich (except for in- output)
- Water input is constant in time and uniform over depth range of 0-2m
- Water output is constant in time and uniform over depth range of 0-1m
- Water input equals water output (2 m^3/s)
- Input temperature is 20°C
- Input salinity is 0.1 %
- The simulation lasts for 300 days

Result after 300 days:
- Lake surface:		136.00000 m
- Surface T:		14.49084 °C
- Bottom T:			4.54444 °C

Testcase 2: variable water table

- Data from Lake Zurich (except for in- output)
- Water input is variable in time but uniform over depth range of 0-2m
- Water output is constant in time and uniform over depth range of 0-2m
- Mean water input equals mean water output (2 m^3/s)
- Input temperature is 20°C
- Input salinity is 0.1 %
- The simulation lasts for 300 days

Result after 300 days:
- Lake surface:		136.04271 m (bug in treatment of variable inflow at surface)
- Surface T:		15.14825 °C
- Bottom T:			4.54444 °C
