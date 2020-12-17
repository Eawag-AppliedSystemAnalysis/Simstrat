# About this test case
This folder contains the test simulation of Lake Zurich with real bathymetry and forcing data. On the other hand, the inflow data and initial conditions are realistic, but fictional. The main goals of this test are i) to determine whether
a newly compiled executable works as intended and ii) to check whether changed to the source files cause changed in the simulation output.
From this folder (`tests`), run:

~~~bash
.\run.bat
~~~

under windows, or run:

~~~bash
.\run_testcase
~~~

under linux to start a simulation with the executable "simstrat.exe" and the Simstrat configuration file "TestCase_LakeZurich.par". The results are stored in ´TestCases_Results´. The configuration file of AED2 "aed2.nml" is
thereby called by the Simstrat configuration file.

To check whether changes in the source code caused changes in simulation output run:

~~~bash
python .\tests.py
~~~
