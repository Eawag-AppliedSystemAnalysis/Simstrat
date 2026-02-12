# Load libraries

import os
from functions import csv_to_netcdf

# Set parameters

variables = ['T', 'S']
var_units = ['degree_Celsius', '1e-3']
path_to_output = '../../tests/TestCases_Results/'
conversion_function = 'csv_to_netcdf' # Currently available: csv_to_netcdf

# Execute conversion

if conversion_function == 'csv_to_netcdf':
    for i, variable in enumerate(variables):
        csv_file = os.path.join(path_to_output, f'{variable}_out.dat')
        netcdf_file = os.path.join(path_to_output, f'{variable}_out.nc')
        if len(var_units) >= (i + 1):
            csv_to_netcdf(csv_file, netcdf_file, var_units[i])
        else:
            csv_to_netcdf(csv_file, netcdf_file)
else:
    raise ValueError(f'Conversion {conversion_function} not available.')