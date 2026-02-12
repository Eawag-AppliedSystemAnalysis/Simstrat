# Load libraries

import os
import pandas as pd
import numpy as np
from netCDF4 import Dataset



# Main function

def csv_to_netcdf(csv_file, netcdf_file, var_unit=''):
    """
    Convert Output file from CSV to NetCDF
    
    Parameters:
    -----------
    csv_file : str
        Path to output file in csv format
        File has to exist
    netcdf_file : str
        Path to output file in netcdf format
        Creates a new file or overwrites the existing one
    var_unit : str, optional
        The unit of the variable (e.g., 'degree_Celsius', 'kg m-3'), default is an empty string
    """
    
    # Check if the CSV file exists
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f'File {csv_file} does not exist.')

    # Read CSV into pandas DataFrame
    df = pd.read_csv(csv_file)

    # Check if the file is empty
    if df.empty:
        raise ValueError(f'File {csv_file} is empty.')

    # Extract the Datetime and Depth columns
    datetime_vals = df['Datetime'].values
    depth_vals = df.columns[1:].astype(np.float64)

    # Create NetCDF file
    with Dataset(netcdf_file, 'w', format='NETCDF4') as ncfile:

        # Create Datetime and Depth dimensions
        ncfile.createDimension('Datetime', len(datetime_vals))
        ncfile.createDimension('Depth', len(depth_vals))

        # Create variables for 'Datetime' and 'Depth'
        datetime_var = ncfile.createVariable('Datetime', 'f8', ('Datetime',))
        depth_var = ncfile.createVariable('Depth', 'f8', ('Depth',))
        datetime_var.units = 'd'
        depth_var.units = 'm'

        # Assign the Datetime and Depth values to the corresponding variables
        datetime_var[:] = datetime_vals
        depth_var[:] = depth_vals

        # Create the data variable (Datetime, Depth)
        data_var = ncfile.createVariable('Data', 'f8', ('Datetime', 'Depth'))
        data_var.units = var_unit

        # Fill the data variable (the DataFrame values corresponding to Datetime and Depth)
        # The DataFrame starts with 'Datetime' and the rest are the depth columns
        data_matrix = df.iloc[:, 1:].values
        data_var[:, :] = data_matrix