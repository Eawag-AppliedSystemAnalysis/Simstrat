# How to build Simstrat-FABM with FoBiS.py

Go to `lib/fabm` and run:

~~~bash
cmake -S ./ -B build -DFABM_HOST=simstrat
~~~

to create the new directory `build` and generate the build configuration inside it,
with the Simstrat preprocesser definitions for information about the spatial domain.

Then run:

~~~bash
cmake --build build --target install
~~~

to build and install FABM.

Then, from `Simstrat/build`, run:

~~~bash
FoBiS.py -h
~~~

to get help information about its usage.
To compile Simstrat with default configuration, run:

~~~bash
FoBiS.py build
~~~

## Clean the project
To clean the main project, simply run

~~~bash
FoBiS.py clean
~~~

whilst if you need to clean also all the libraries (e.g. if you change the final target from Linux to Win), you need to run

~~~bash
FoBiS.py rule -ex purge
~~~

## Generate the code documentation
To generate the code documentation, run

~~~bash
FoBiS.py rule -ex makedoc
~~~

The generated code documentation is saved in `doc/developer/ford/ford_doc/index.htlm`


**N.B.** you can use one-line command to call the build procedure (and others) from any folder, e.g. from `tests`
~~~bash
cd ../build; FoBiS.py build; cd -
~~~
