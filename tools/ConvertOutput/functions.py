# Load libraries

import os
import pandas as pd
import numpy as np
from netCDF4 import Dataset



# Main function

def csv_to_netcdf(var_names, filename, path_to_output, var_units=[], dt_freq='h'):
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
    var_unit : str, optional
        The units of the variables (e.g., 'degree_Celsius', '1e-3'), default is an empty list
    dt_freq : str, optional
        The frequency the datetime dimension is round to, must be 's', 'min', 'h' or 'd'
        Default is 'h': round to hours
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
        # Convert Time column (days since reference data) to Datetime
        # Assuming reference Date is 1981-01-01 (common for Simstrat)
        df = pd.read_csv(csv_file)
        depth_vals = df.columns[1:].astype(np.float64)
        time_vals = df['Datetime']
        datetime_vals = (pd.Timestamp('1981-01-01') + pd.to_timedelta(time_vals, unit='D').astype('timedelta64[s]'))
        datetime_vals = datetime_vals.dt.round(freq=dt_freq)
        datetime_vals = datetime_vals.astype(str)
        if dt_freq == 'd':
            datetime_vals = datetime_vals.astype('|S%d' % 10)
        elif dt_freq == 'h':
            datetime_vals = datetime_vals.astype('|S%d' % 13)
        elif dt_freq == 'min':
            datetime_vals = datetime_vals.astype('|S%d' % 16)
        elif dt_freq == 's':
            datetime_vals = datetime_vals.astype('|S%d' % 19)
        else:
            raise ValueError(f'Datetime frequency must be s, min, h or d')

        # Once is enough
        break

    # Create NetCDF file
    netcdf_file = os.path.join(path_to_output, f'{filename}.nc')
    with Dataset(netcdf_file, 'w', format='NETCDF4') as ncfile:

        # Create Depth and Datetime dimensions
        ncfile.createDimension('Depth', len(depth_vals))
        ncfile.createDimension('Datetime', len(datetime_vals))

        # Create variables for 'Depth' and 'Datetime'
        depth_var = ncfile.createVariable('Depth', 'f8', ('Depth',))
        depth_var.units = 'm'
        datetime_var = ncfile.createVariable('Datetime', str, ('Datetime',))

        # Assign the Depth and Datetime values to the corresponding variables
        depth_var[:] = depth_vals
        for i, dt in enumerate(datetime_vals):
            datetime_var[i] = dt

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