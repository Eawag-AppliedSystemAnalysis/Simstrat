# Load libraries

import os
from functions import csv_to_netcdf

# Required parameters

var_names = [] # Variable names; Pass an empty list to convert all variables in the output folder
filename = 'output_converted' # Name of new file
path_to_output = '../../tests/LakeConstance/Results' # Path to output folder
conversion_function = 'csv_to_netcdf' # Currently available: csv_to_netcdf

# Optional parameters

var_units = [] # Variable units, assigned to var_names in order
var_exclude = ['TotalIceH'] # If empty list is passed for var_names, exclude files that are empty or of different format

# Execute conversion

if conversion_function == 'csv_to_netcdf':
    if len(var_names) == 0:
        var_names = [f[:-8] for f in os.listdir(path_to_output) if ((f[-8:] == '_out.dat') and not (f[:-8] in var_exclude))]
    csv_to_netcdf(var_names, filename, path_to_output, var_units)
else:
    raise ValueError(f'Conversion {conversion_function} not available.')