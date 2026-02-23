# Load libraries

import os
import pandas as pd
import numpy as np
from netCDF4 import Dataset



# Main function

def csv_to_netcdf(var_names, filename, path_to_output, var_units=[]):
    """
    Convert Output file from CSV to NetCDF
    
    Parameters:
    -----------
    var_names : str
        The names of the variables to convert the output from (e.g. 'T', 'S')
        Output file has to exist
    filename : str
        Name of file containing converted output
        Creates a new file or overwrites the existing one
    path_to_output : str
        Path to folder with output files
    var_units : str, optional
        The units of the variables (e.g., 'degree_Celsius', '1e-3'), default is an empty list
    """

    # Get Depth and Datetime dimensions from first files with values
    for variable in var_names:

        # Create path and check if CSV file exists and if it is not empty
        csv_file = os.path.join(path_to_output, f'{variable}_out.dat')
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f'File {csv_file} does not exist.')
        if os.stat(csv_file).st_size == 0:
            print(f'File {csv_file} is empty.')
            continue

        # Extract the Depth and Datetime columns
        df = pd.read_csv(csv_file)
        depth_vals = df.columns[1:].astype(np.float64)
        datetime_vals = df['Datetime'].values.astype(np.float64)

        # Once is enough
        break

    # Create NetCDF file
    netcdf_file = os.path.join(path_to_output, f'{filename}.nc')
    with Dataset(netcdf_file, 'w', format='NETCDF4') as ncfile:

        # Create Depth and Datetime dimensions
        ncfile.createDimension('Depth', len(depth_vals))
        ncfile.createDimension('Datetime', len(datetime_vals))

        # Create variables for 'Depth' and 'Datetime'
        # Assuming reference date is 1981-01-01 (common for Simstrat)
        depth_var = ncfile.createVariable('Depth', 'f8', ('Depth',))
        depth_var.units = 'm'
        datetime_var = ncfile.createVariable('Datetime', 'f8', ('Datetime',))
        datetime_var.units = 'days since 1981-01-01 00:00:00'
        datetime_var.calendar = 'standard'

        # Assign the Depth and Datetime values to the corresponding variables
        depth_var[:] = depth_vals
        datetime_var[:] = datetime_vals

        # Add data variables (Depth, Datetime)
        for i, variable in enumerate(var_names):

            # Create path and check if CSV file exists and if it is not empty
            csv_file = os.path.join(path_to_output, f'{variable}_out.dat')
            if not os.path.exists(csv_file):
                raise FileNotFoundError(f'File {csv_file} does not exist.')
            if os.stat(csv_file).st_size == 0:
                print(f'File {csv_file} is empty.')
                continue

            # Create the data variable (Depth, Datetime)
            data_var = ncfile.createVariable(var_names[i], 'f8', ('Depth', 'Datetime'))
            
            # Add data variable units, if provided
            if len(var_units) >= (i + 1):
                data_var.units = var_units[i]
            else:
                data_var.units = ''

            # Fill the data variable (the DataFrame values corresponding to Datetime and Depth)
            # The DataFrame starts with 'Datetime' and the rest are the depth columns
            # Transpose to match the data variable dimensions (Depth, Datetime)
            df = pd.read_csv(csv_file)
            data_matrix = df.iloc[:, 1:].values
            data_var[:, :] = np.transpose(data_matrix)