# Convert Output data

Here you can convert the output files to another format than CSV.

## Available conversion functions:

- `csv_to_netcdf(csv_file, netcdf_file, var_units='')` — Convert Output file from CSV to NetCDF

## Use the scripts:

First set the variable names, the name of the file containing the converted output, the path to the output files, as well as the conversion function in `Simstrat/tools/ConvertOutput/convert_output.py`. If input should also be converted and added to the output, set the paths to input files.

To execute the conversion inside the Simstrat build environment make sure you are in `Simstrat/tools/ConvertOutput` and run:

~~~bash
python3 ./convert_output.py
~~~

> **N.B.** The input conversion is quite robust but not 100% to all input formats
> Forcing, Bathymetry and Initial Conditions require header of format: variable [unit]
> For the units Inflow and Absorption require header of format: variable [unit]