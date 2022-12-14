# About this test case
This folder contains the test simulation of Lake Zurich with real bathymetry and forcing data. On the other hand, the inflow data and initial conditions are realistic, but fictional. The main goals of this test are i) to determine whether
a newly compiled executable works as intended and ii) to check whether changed to the source files cause changed in the simulation output.
From this folder (`tests`), run:

~~~bash
.\run.bat
~~~

under windows, or run:

~~~bash
./run_testcase.sh
~~~

under linux to start a simulation with the executable "simstrat" and the Simstrat configuration file "TestCase_LakeZurich.par". The AED2 configuration file ("aed2.nml") is called by the Simstrat configuration file and the output of the simulation is stored in ´TestCases_Results´.

To check whether changes in the source code caused changes in simulation output run:

~~~bash
python .\tests.py
~~~

under windows, or run:

~~~bash
python3 ./tests.py
~~~

under linux.

# Expected test case results

The expected results of the test cases are always generated using the latest release. Current Simstrat release is 3.02 and current AED2 release is 1.3.5. In windows, Simstrat and AED2 were compiled using gfortran 8.1 and in Linux using gfortran 7.5.
Fobos settings used for Simstrat calibration can be found in the latest release