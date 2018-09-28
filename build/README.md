# How to build Simstrat (via build.sh)
The build script of Simstrat is `build/build.sh`. From this folder (`build`), run

~~~bash
./build.sh -h
~~~

to get help information about its usage:

~~~
Usage: ./build.sh [-m <string>] [-c] [-d]

  -m <string>  compile mode: release (default), debug
  -c           clean previous compiling files before building
  -d           generate code documentation with FORD after building
 in case of no arguments, "./build.sh -m release" is called
~~~

> **N.B.** you can use one-line command to call the script from the root repository folder, e.g.
~~~bash
cd build; ./build.sh; cd -
~~~

## Code documentation with FORD
You can find the generated code documentation (flag `-d`) in `doc/developer/ford/ford_doc/index.htlm`

> **N.B.** you can generate the documentation in any moment from `doc/developer/ford` folder by running
~~~bash
ford ford_projectfile.md
~~~