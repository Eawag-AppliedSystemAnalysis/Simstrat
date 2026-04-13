# Load libraries

import os
import pandas as pd
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from functions_input import read_input_paths, load_absorption_data, load_inflow_data, load_initial_data, load_forcing_data

# Main function

def csv_to_netcdf(var_names, filename, path_to_output, paths_to_input, inflow_mode_not_1=True):
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
    inflow_mode_not_1 : bool, optional
        Whether inflow mode in configuration file is not 1 (Manual inflow mode)
        (Default: True)
    """

    # Initialise empty list of DataArray
    data_list = []

    # Get variable attributes from list_variables
    csv_file = os.path.join(path_to_output, f'list_variables.dat')
    if not os.path.exists(csv_file):
        print(f'Variable attribute file {csv_file} does not exist.')
    elif os.stat(csv_file).st_size == 0:
        print(f'Empty variable attribute file {csv_file} ignored.')
    else:
        attributes = pd.read_csv(csv_file, index_col=0, sep = ', ')
    
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

        # Save time values and drop time column
        time_values = df['Datetime']
        df = df.drop(columns='Datetime')

        # Convert to Data Array
        try:
            data = xr.DataArray(df.values.astype(float), dims=['Datetime', 'Depth'], coords={'Datetime': time_values.astype(float), 'Depth': df.columns.astype(float)}, name=var_names[i])
        # Sometimes the last line is not fully read which can cause a ValueError when converting to float
        except ValueError:
            df = df[:len(df)-1]
            time_values = time_values[:len(time_values)-1]
            data = xr.DataArray(df.values.astype(float), dims=['Datetime', 'Depth'], coords={'Datetime': time_values.astype(float), 'Depth': df.columns.astype(float)}, name=var_names[i])
        # If the value is constant in time, drop that dimension
        if np.all(data['Datetime'].values == data['Datetime'].values[0]):
            data = data.isel(Datetime=0, drop=True)
        # If the value is constant in depth, drop that dimension (necessary for FABM surface diagnostic variables)
        if np.all(data['Depth'].values == data['Depth'].values[0]):
            data = data.isel(Depth=0, drop=True)
        # Add attributes
        data.attrs['long_name'] = attributes['Long Name'][var_names[i]]
        data.attrs['units'] = attributes['Units'][var_names[i]]
        # Append to data list
        data_list.append(data)

    # Add data variables from input folders
    for path_to_input in paths_to_input:

        # Get sorted paths
        absorption_input_paths, inflow_input_paths, initial_input_paths, forcing_input_paths = read_input_paths(path_to_input)

        # Add absorption-type input data
        for i, file_path in enumerate(absorption_input_paths):
            # Load data
            df, units = load_absorption_data(file_path)
            # Store name and depths
            name = os.path.splitext(os.path.basename(file_path))[0]
            depths = [depth.replace(' ', '').split('(')[0] for depth in df.columns[1:]]
            # Set column names
            if len(df.columns) == 2:
                df.columns = ['Datetime', name]
            elif len(df.columns) > 2:
                df.columns = ['Datetime'] + [(name + '_' + str(j+1)) for j in range(len(df.columns)-1)]
            # Aggregate duplicate rows
            df = df.groupby('Datetime', as_index=False).sum()
            # Store every column as DataArray
            for j in range(len(df.columns)-1):
                data = df.set_index(['Datetime']).to_xarray()[df.columns[j+1]]
                if len(units) > 0:
                    data.attrs['units'] = units.iloc[0,j]
                data.attrs['Reference Depth'] = depths[j]
                data_list.append(data)
                
        # Add inflow-type input data
        for i, file_path in enumerate(inflow_input_paths):
            # Load deep and surface data
            deep, surf, units = load_inflow_data(file_path)
            # Store name and depths
            name = os.path.splitext(os.path.basename(file_path))[0]
            deep_depths = [depth.replace(' ', '').split('(')[0] for depth in deep.columns[1:]]
            surf_depths = [depth.replace(' ', '').split('(')[0] for depth in surf.columns[1:]]
            # Set column names
            # Deep inflow mode 2
            if inflow_mode_not_1:
                if len(deep.columns) == 2:
                    deep.columns = ['Datetime', name + '_deep']
                elif len(deep.columns) > 2:
                    deep.columns = ['Datetime'] + [(name + '_deep_' + str(j+1)) for j in range(len(deep.columns)-1)]
            # Deep inflow mode 1
            else:
                if len(deep.columns) == 3:
                    deep.columns = ['Datetime', name + '_deep', '']
                elif len(deep.columns) > 3:
                    deep.columns = ['Datetime'] + [str(int((j+1)/2)) if (j % 2) else (name + '_deep_' + str(int(j/2)+1)) for j in range((len(deep.columns)-1))] 
            # Surface inflow (as deep inflow with mode 1)
            if len(surf.columns) == 3:
                surf.columns = ['Datetime', name + '_surface', '']
            elif len(surf.columns) > 3:
                surf.columns = ['Datetime'] + [str(int((j+1)/2)) if (j % 2) else (name + '_surface_' + str(int(j/2)+1)) for j in range((len(surf.columns)-1))]          
            # Aggregate duplicate rows
            deep = deep.groupby('Datetime', as_index=False).sum()
            surf = surf.groupby('Datetime', as_index=False).sum()
            # Store every column as DataArray
            # Deep inflow mode 2
            if inflow_mode_not_1:
                for j in range(len(deep.columns)-1):
                    data = deep.set_index(['Datetime']).to_xarray()[deep.columns[j+1]]
                    if len(units) > 0:
                        data.attrs['units'] = units.iloc[0,0]
                    data.attrs['Injection Depth'] = deep_depths[j]
                    data_list.append(data)
            # Deep inflow mode 1
            else:
                for j in range(len(deep.columns)-1):
                    if not (j % 2):
                        data = deep.set_index(['Datetime']).to_xarray()[deep.columns[j+1]]
                        if len(units) > 0:
                            data.attrs['units'] = units.iloc[0,1]
                        data.attrs['Injection Depth'] = str(deep_depths[j]) + ' m to ' + str(deep_depths[j+1]) + ' m'
                        data_list.append(data)
            # Surface inflow (as deep inflow with mode 1)
            for j in range(len(surf.columns)-1):
                if not (j % 2):
                    data = surf.set_index(['Datetime']).to_xarray()[surf.columns[j+1]]
                    if len(units) > 0:
                        data.attrs['units'] = units.iloc[0,1]
                    data.attrs['Injection Depth'] = str(surf_depths[j]) + ' m to ' + str(surf_depths[j+1]) + ' m'
                    data_list.append(data)
                
        # Add inital-type input data
        for i, file_path in enumerate(initial_input_paths):
            # Load data
            df, units = load_initial_data(file_path)
            # Store name
            name = os.path.splitext(os.path.basename(file_path))[0]
            # Set column names
            df.columns = ['Depth'] + [(name + '_' + df.columns[j+1]) for j in range(len(df.columns)-1)]
            # Aggregate duplicate rows
            df = df.groupby('Depth', as_index=False).sum()
            # Store every column as DataArray
            for j in range(len(df.columns)-1):
                data = df.set_index(['Depth']).to_xarray()[df.columns[j+1]]
                if len(units) > 0:
                    data.attrs['units'] = units.iloc[0,j]
                data_list.append(data)
                
        # Add forcing-type input data
        for i, file_path in enumerate(forcing_input_paths):
            # Load data
            df, units = load_forcing_data(file_path)
            # Store name
            name = os.path.splitext(os.path.basename(file_path))[0]
            # Set column names
            df.columns = ['Datetime'] + [(name + '_' + df.columns[j+1]) for j in range(len(df.columns)-1)]
            # Aggregate duplicate rows
            df = df.groupby('Datetime', as_index=False).sum()
            # Store every column as DataArray
            for j in range(len(df.columns)-1):
                data = df.set_index(['Datetime']).to_xarray()[df.columns[j+1]]
                if len(units) > 0:
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