# Load libraries

import os
from functions_conversion import csv_to_netcdf

# Set parameters

var_names = [] # Variable names; pass an empty list to convert all variables in the output folder
path_to_output = '../../tests/TestCases_Results' # Path to output folder
paths_to_input = ['../../tests/TestCase_LakeZurich', '../../tests/TestCase_LakeZurich/FABM_initialconditions', '../../tests/TestCase_LakeZurich/FABM_inflow'] # Paths to input folders
eps =  1e-10 # Tolerance for variation among dimension to drop that dimension (optional, default: 1e-10)
skip_initial = False # Whether to skip the first timepoint (optional, default: False)
conversion_function = 'csv_to_netcdf' # Currently available: csv_to_netcdf
inflow_mode = 2 # Inflow mode for deep inflows as set in configuration file
filename = 'output_converted' # Name of new file

# Execute conversion

if conversion_function == 'csv_to_netcdf':
    if len(var_names) == 0:
        var_names = [f[:-8] for f in os.listdir(path_to_output) if (f[-8:] == '_out.dat')]
    if inflow_mode == 1:
        csv_to_netcdf(var_names, filename, path_to_output, paths_to_input, eps, skip_initial, False)
    else:
        csv_to_netcdf(var_names, filename, path_to_output, paths_to_input, eps, skip_initial)
else:
    raise ValueError(f'Conversion {conversion_function} not available.')