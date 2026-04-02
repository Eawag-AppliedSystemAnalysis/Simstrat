# Load libraries

import os
import pandas as pd
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from functions_input import read_input_paths, load_absorption_data, load_inflow_data, load_initial_data, load_forcing_data

# Main function

def csv_to_netcdf(var_names, filename, path_to_output, paths_to_input, start_day=0.0, var_units=[], skiprows=0, collapse=True):
    """
    Convert Output file from CSV to NetCDF
    
    Parameters:
    -----------
    var_names : str list
        The names of the variables to convert the output from (e.g. 'T', 'S')
        Output file has to exist
    filename : str
        Name of file containing converted output
        Creates a new file or overwrites the existing one
    path_to_output : str
        Path to folder with output files
    paths_to_input : str list
        Paths to folders with input files 
    var_units : str list, optional
        The units of the variables (e.g., 'degree_Celsius', '1e-3'), default is an empty list
    collapse : bool, optional
        Whether to collapse data matrices with same value at every depth into a time series
        Default is True
    """

    # Initialise empty list of DataArray
    data_list = []

    # Add data variables from output folder
    for i, variable in enumerate(var_names):

        # Create path and check if CSV file exists and if it is not empty
        csv_file = os.path.join(path_to_output, f'{variable}_out.dat')
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f'File {csv_file} does not exist.')
        elif os.stat(csv_file).st_size == 0:
            print(f'Empty file {csv_file} ignored.')
            continue

        # Read values from the CSV file
        # The DataFrame starts with 'Datetime' and the rest are the depth columns
        df = pd.read_csv(csv_file)

        # Drop incomplete rows, save time values and drop time column
        df = df.dropna()
        time_values = df['Datetime']
        df = df.drop(columns='Datetime')

        # Convert to Data Array
        data = xr.DataArray(df.values.astype(float), dims=['Datetime', 'Depth'], coords={'Datetime': time_values.astype(float), 'Depth': df.columns.astype(float)}, name=var_names[i])
        
        # If the value is constant in time, drop that dimension
        if (data.std(dim='Datetime') == 0).all():
            data = data.isel(Datetime=0, drop=True)
        # If the value is constant in depth, drop that dimension (necessary for surface diagnostic variables)
        if (data.std(dim='Depth') == 0).all():
            data = data.isel(Depth=0, drop=True)

        # Append to data list
        data_list.append(data)

    # Add data variables from input folders
    for path_to_input in paths_to_input:

        # Get sorted paths
        absorption_input_paths, inflow_input_paths, initial_input_paths, forcing_input_paths = read_input_paths(path_to_input)

        # Add absorption-type input data
        for i, file_path in enumerate(absorption_input_paths):
            name = os.path.splitext(os.path.basename(file_path))[0]
            df, units = load_absorption_data(file_path, start_day)
            depths = [depth.replace(' ', '').split('(')[0] for depth in df.columns[1:]]
            if len(df.columns) < 2:
                continue
            if len(df.columns) < 3:
                df.columns = ['Datetime', name]
            else:
                df.columns = ['Datetime'] + [(name + '_' + str(j+1)) for j in range(len(df.columns)-1)]
            for j in range(len(df.columns)-1):
                data = df.set_index(['Datetime']).to_xarray()[df.columns[j+1]]
                data.attrs['units'] = units.iloc[0,j]
                data.attrs['Reference Depth'] = depths[j]
                data_list.append(data)

        # Add inflow-type input data
        for i, file_path in enumerate(inflow_input_paths):
            name = os.path.splitext(os.path.basename(file_path))[0]
            deep, surf, units = load_inflow_data(file_path, start_day)
            deep_depths = [depth.replace(' ', '').split('(')[0] for depth in deep.columns[1:]]
            if len(deep.columns) < 2:
                continue
            if len(deep.columns) < 3:
                deep.columns = ['Datetime', name + '_deep']
            else:
                deep.columns = ['Datetime'] + [(name + '_deep_' + str(j+1)) for j in range(len(deep.columns)-1)]
            for j in range(len(deep.columns)-1):
                data = deep.set_index(['Datetime']).to_xarray()[deep.columns[j+1]]
                data.attrs['units'] = units.iloc[0,0]
                data.attrs['Injection Depth'] = deep_depths[j]
                data_list.append(data)
            surf_depths = [depth.replace(' ', '').split('(')[0] for depth in surf.columns[1:]]
            if len(surf.columns) < 2:
                continue
            if len(surf.columns) < 3:
                surf.columns = ['Datetime', name + '_surface_inflow']
            else:
                surf.columns = ['Datetime'] + [(name + '_surface_inflow_' + str(j+1)) for j in range(len(surf.columns)-1)]
            for j in range(len(surf.columns)-1):
                data = surf.set_index(['Datetime']).to_xarray()[surf.columns[j+1]]
                data.attrs['units'] = units.iloc[0,1]
                data.attrs['Injection Depth'] = surf_depths[j]
                data_list.append(data)

        # Add inital-type input data
        for i, file_path in enumerate(initial_input_paths):
            name = os.path.splitext(os.path.basename(file_path))[0]
            df, units = load_initial_data(file_path)
            df.columns = ['Depth'] + [(name + '_' + df.columns[j+1]) for j in range(len(df.columns)-1)]
            for j in range(len(df.columns)-1):
                data = df.set_index(['Depth']).to_xarray()[df.columns[j+1]]
                data.attrs['units'] = units.iloc[0,j]
                data_list.append(data)

        # Add forcing-type input data
        for i, file_path in enumerate(forcing_input_paths):
            name = os.path.splitext(os.path.basename(file_path))[0]
            df, units = load_forcing_data(file_path, start_day)
            df.columns = ['Datetime'] + [(name + '_' + df.columns[j+1]) for j in range(len(df.columns)-1)]
            for j in range(len(df.columns)-1):
                data = df.set_index(['Datetime']).to_xarray()[df.columns[j+1]]
                data.attrs['units'] = units.iloc[0,j]
                data_list.append(data)

    # Assemble all Variables into a Dataset
    full_data = xr.Dataset()
    for data_listed in data_list:
        full_data[data_listed.name] = data_listed

    # Add coordinate units
    full_data['Depth'].attrs['units'] = 'm'
    full_data['Datetime'].attrs['units'] = 'days since 1981-01-01 00:00:00'
    full_data['Datetime'].attrs['calendar'] = 'standard'

    # Create NetCDF file from Dataset
    netcdf_file = os.path.join(path_to_output, f'{filename}.nc')
    full_data.to_netcdf(netcdf_file)