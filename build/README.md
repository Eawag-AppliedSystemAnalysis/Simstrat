# How to build Simstrat with FoBiS.py
From this folder (`build`), run:

~~~bash
FoBiS.py -h
~~~

to get help information about its usage.
To compile Simstrat with default configuration, run:

~~~bash
FoBiS.py build
~~~

If you want to generate code documentation, you can run

~~~bash
FoBiS.py rule -ex makedoc
~~~

The generated code documentation is saved in `doc/developer/ford/ford_doc/index.htlm`


> **N.B.** you can use one-line command to call the build procedure from any folder, e.g. from `tests`
~~~bash
cd ../build; FoBiS.py build; cd -
~~~