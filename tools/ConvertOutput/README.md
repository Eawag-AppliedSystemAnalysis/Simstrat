# Convert Output data

Here you can convert the output files to another format than CSV.

## Available conversion functions:

- `csv_to_netcdf(csv_file, netcdf_file, var_units='')` — Convert Output file from CSV to NetCDF

## Use the scripts:

First set the variables, their units (optionally) and the path to the output files, as well as the conversion function in `Simstrat/tools/ConvertOutput/convert_output.py`.

To execute the conversion inside the Simstrat build environment make sure you are in `Simstrat/tools/ConvertOutput` and run:

~~~bash
python3 ./convert_output.py
~~~