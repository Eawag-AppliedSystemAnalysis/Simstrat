# Load libraries

import os
from functions import csv_to_netcdf

# Set parameters

variable = 'T'
var_units = 'degree_Celsius'
path_to_output = '../../tests/TestCases_Results/'
csv_file = os.path.join(path_to_output, f'{variable}_out.dat')
netcdf_file = os.path.join(path_to_output, f'{variable}_out.nc')

# Chose conversion to execute
# Currently available: csv_to_netcdf

csv_to_netcdf(csv_file, netcdf_file, var_units)