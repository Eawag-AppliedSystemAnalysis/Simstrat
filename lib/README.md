# External (3rd party) Libraries

This directory contains external libraries used in Simstrat.

For now we have:

| Library Name | Version / Date             | License | Linking      |
|--------------|----------------------------|---------|--------------|
| csv_fortran  | 7.1 /19.07.2019			| BSD	  | static       |
| json_fortan  | ??           				| BSD 	  | static       |
| libaed2	   | 1.3 /20.08.2019			| GPL-3.0 | static       |

Please note that the first two libraries are hard copies, while the `libaed2` library is a submodule, added in the first place as
~~~bash
git submodule add remote_git_url_to_external_library path_to_this_directory
~~~

Submodule repositories are not updated automatically. This has to be done manually from time to time, or if needed. Please be sure to run extensive tests to make sure that any changes in the external library does not destroy our project!


# How to clone the Simstrat project with submodules in your local copy

The **first time** you clone locally our repository you have to initialize and update the submodules with the following commands

~~~bash
git submodule init
git submodule update
~~~

> **NB:** for most of the cases, the above two commands is all what you need; what follows is needed if you will modify the libraries.

If new commits (of the submodules) are pushed to our project, everybody should keep their local submodules clones up-to-date with

~~~bash
git submodule update
~~~

Any changes in the file `../.gitmodules` (e.g. changing the branch inside the submodule) should be followed by calling
~~~bash
git submodule sync
~~~

CAUTION, THIS MIGHT HAVE UNPREDICTABLE IMPACT ON OUR PROJECT: To keep the local submodule clones up-to-date with the remote submodule repository (i.e. pulling), you call
~~~bash
git submodule update --init --recursive --remote
~~~