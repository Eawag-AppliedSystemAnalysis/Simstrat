# How to build Simstrat-AED2 with FoBiS.py
In case the AED2 library is not compiled go to ´lib/libaed2´ and run:

~~~bash
export F90=gfortran
mingw32-make
~~~

under windows (note that the export command only works in a bash shell) or

~~~bash
export F90=gfortran
make
~~~

under linux.

Then, from this folder (`build`), run:

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
