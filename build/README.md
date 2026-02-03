## How to build Simstrat-FABM with FoBiS.py

Make sure you're inside your Simstrat build environment and go to `Simstrat/lib/fabm`. There run:

~~~bash
cmake -S ./ -B build -DFABM_HOST=simstrat
~~~

to create the new directory `build` and generate the build configuration inside it, with the Simstrat preprocesser definitions for information about the spatial domain.

> **N.B.1** By default FABM includes all biogeochemical models from `lib/fabm/src/models`. You can restrict this to a subset of directories by providing: 
> ~~~bash
> -DFABM_INSTITUTES=<INSTITUTE-NAMES>
> ~~~  
> as an argument to the `cmake` command above. To add externally maintained models add the argument: 
> ~~~bash
> -DFABM_EXTRA_INSTITUTES=<INSTITUTES-NAMES>
> ~~~  
> to the `cmake` command. `<INSTITUTE-NAMES>` is a semi-colon-separated list of institute names (all in lower case). You may need to enclose this list in quotes to prevent the shell from interpreting the semi-colon. For the externally maintained models, you additionally need to clone their repositories and point FABM to the directory with the source code by adding the argument:
> ~~~bash
> -DFABM_<INSTITUTE-NAME>_BASE=<DIR>
> ~~~  
> to the `cmake` command for every `<INSTITUTE-NAME>` (here all in upper case) in `<INSTITUTE-NAMES>`. With `<DIR>` the directory of the corresponding source code on your device. A list of available externally maintained models can be found here: [Available Models in FABM Wiki](https://github.com/fabm-model/fabm/wiki/Biogeochemical-models-in-FABM#available-models). For instance to build FABM with [the Selmaprotbas model](https://github.com/jorritmesman/Selmaprotbas) the full cmake command is:
> ~~~bash
> cmake -S ./ -B build -DFABM_HOST=simstrat -DFABM_EXTRA_INSTITUTES=selmaprotbas -DFABM_SELMAPROTBAS_BASE=pathto/selmaprotbas
> ~~~ 
> with `pathto` the path to the location of the selmaprotbas repository on your device, relative to `lib/fabm`.

> **N.B.2** Add:
> ~~~bash
> -DCMAKE_BUILD_TYPE=Debug
> ~~~  
> as argument to `cmake` to compile FABM in debug mode.

Then run:

~~~bash
cmake --build build --target install
~~~

Now leave `Simstrat/lib/fabm` and go to `Simstrat/build` with:

~~~bash
cd ../../build
~~~

from where you can build Simstrat with:

~~~bash
FoBiS.py build
~~~

> **N.B.1** macOS seems to allow only dynamic library linking. To build for this target, use the `release-gnu-dynamic` compiling mode:
> ~~~bash
> FoBiS.py build -mode release-gnu-dynamic
> ~~~

> **N.B.2** you can find more building options [here](build).

> **N.B.3** you can use one-line command to call the final Simstrat build procedure from any folder, e.g. from `tests`:
~~~bash
cd ../build; FoBiS.py build; cd -
~~~

## Clean the project

To clean the main project, simply run:

~~~bash
FoBiS.py clean
~~~

whilst if you need to clean also all the libraries (e.g. if you change the final target from Linux to Win), you need to run:

~~~bash
FoBiS.py rule -ex purge
~~~

## Generate the code documentation

The developer documentation can be generated with the FORD python module (`pip install ford`).
To generate the code documentation, run:

~~~bash
FoBiS.py rule -ex makedoc
~~~

The generated code documentation is saved in `doc/developer/ford/ford_doc/index.htlm`
