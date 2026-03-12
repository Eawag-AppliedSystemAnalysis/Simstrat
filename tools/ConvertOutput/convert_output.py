# Load libraries

import os
from functions import csv_to_netcdf

# Required parameters

var_names = ['T', 'S'] # Variable names; pass an empty list to convert all variables in the output folder
filename = 'output_converted' # Name of new file
path_to_output = '../../tests/TestCases_Results/' # Path to output folder
conversion_function = 'csv_to_netcdf' # Currently available: csv_to_netcdf

# Optional parameters

var_units = ['degree_Celsius', '1e-3'] # Variable units, assigned to var_names in order
var_exclude = [] # If empty list is passed for var_names, exclude these variables from conversion
var_format = '' # This variable defines the output format; if left empty it is the first variable in var_names

# Execute conversion

if conversion_function == 'csv_to_netcdf':
    if len(var_names) == 0:
        var_names = [f[:-8] for f in os.listdir(path_to_output) if ((f[-8:] == '_out.dat') and not (f[:-8] in var_exclude))]
    if var_format in var_names:
        var_names.remove(var_format)
        var_names.insert(0, var_format)
    csv_to_netcdf(var_names, filename, path_to_output, var_units)
else:
    raise ValueError(f'Conversion {conversion_function} not available.')