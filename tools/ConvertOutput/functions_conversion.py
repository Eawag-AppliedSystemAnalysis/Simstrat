# Load libraries

import os
import re
import warnings
import pandas as pd
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from functions_input import read_input_paths, load_absorption_data, load_inflow_data, load_initial_data, load_forcing_data

# Main function

def csv_to_netcdf(var_names, filename, path_to_output, paths_to_input, eps=1e-10, inflow_mode_not_1=True):
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
    eps : float
        Tolerance for variation among dimension to drop that dimension
        (Default: 1e-20)
    inflow_mode_not_1 : bool, optional
        Whether inflow mode in configuration file is not 1 (Manual inflow mode)
        (Default: True)
    """

    # Initialise empty list of DataArray
    data_list = []

    # Get variable attributes from list_variables
    csv_file = os.path.join(path_to_output, f'_variables.dat')
    if not os.path.exists(csv_file):
        print(f'Variable attribute file {csv_file} does not exist.')
    elif os.stat(csv_file).st_size == 0:
        print(f'Empty variable attribute file {csv_file} ignored.')
    else:
        attributes = pd.read_csv(csv_file, skipinitialspace=True, index_col=0)
    
    # Add data variables from output folder
    for i, variable in enumerate(var_names):

        # Create path and check if CSV file exists and if it is not empty
        # Also check whether it exists in list_variables
        csv_file = os.path.join(path_to_output, f'{variable}_out.dat')
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f'File {csv_file} does not exist.')
        elif os.stat(csv_file).st_size == 0:
            print(f'Empty file {csv_file} ignored.')
            continue
        elif not (var_names[i] in attributes.index):
            print(f'Old file {csv_file} ignored.')
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
        # Drop datetime points that have only nan
        data = data.dropna(dim='Datetime', how='all')
        # Drop depths that have only nan
        data = data.dropna(dim='Depth', how='all')
        # Drop constant dimensions
        with warnings.catch_warnings():
            # Ignore warning appearing for empty results
            warnings.filterwarnings('ignore', message='Degrees of freedom <= 0 for slice')
            # If the value is constant in time, drop that dimension
            datetime_mean = data.mean(dim='Datetime', skipna=True).fillna(0.0)
            if ((np.abs(data - datetime_mean) <= eps * np.abs(datetime_mean)) | data.isnull()).all():
                data = datetime_mean
            # If the value is constant in depth, drop that dimension (necessary for FABM surface diagnostic variables)
            depth_mean = data.mean(dim='Depth', skipna=True).fillna(0.0)
            if ((np.abs(data - depth_mean) <= eps * np.abs(depth_mean)) | data.isnull()).all():
                data = depth_mean
        # Add attributes
        if attributes['Long Name'][var_names[i]] != '-':
            data.attrs['long_name'] = attributes['Long Name'][var_names[i]]
        if normalize_units(attributes['Units'][var_names[i]]) != '-':
            data.attrs['units'] = normalize_units(attributes['Units'][var_names[i]])
        if attributes['Minimum'][var_names[i]] != '-':
            data.attrs['valid_min'] = attributes['Minimum'][var_names[i]]
        if attributes['Maximum'][var_names[i]] != '-':
            data.attrs['valid_max'] = attributes['Maximum'][var_names[i]]
        data.attrs['title'] = attributes['Type'][var_names[i]]
        data.attrs['grid_position'] = attributes['Grid Position'][var_names[i]]
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
                if (len(units) > 0) and (units.iloc[0,j] != '-'):
                    data.attrs['units'] = normalize_units(units.iloc[0,j])
                data.attrs['title'] = 'Absorption Input'
                data.attrs['reference_depth'] = depths[j]
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
                    if (len(units) > 0) and (units.iloc[0,0] != '-'):
                        data.attrs['units'] = normalize_units(units.iloc[0,0])
                    data.attrs['title'] = 'Deep Inflow (Density-driven)'
                    data.attrs['injection_depth'] = deep_depths[j]
                    data_list.append(data)
            # Deep inflow mode 1
            else:
                for j in range(len(deep.columns)-1):
                    if not (j % 2):
                        data = deep.set_index(['Datetime']).to_xarray()[deep.columns[j+1]]
                        if (len(units) > 0) and (units.iloc[0,1] != '-'):
                            data.attrs['units'] = normalize_units(units.iloc[0,1])
                        data.attrs['title'] = 'Deep Inflow (Manual)'
                        data.attrs['injection_depth'] = str(deep_depths[j]) + ' m to ' + str(deep_depths[j+1]) + ' m'
                        data_list.append(data)
            # Surface inflow (as deep inflow with mode 1)
            for j in range(len(surf.columns)-1):
                if not (j % 2):
                    data = surf.set_index(['Datetime']).to_xarray()[surf.columns[j+1]]
                    if (len(units) > 0) and (units.iloc[0,1] != '-'):
                        data.attrs['units'] = normalize_units(units.iloc[0,1])
                    data.attrs['title'] = 'Surface Inflow'
                    data.attrs['injection_depth'] = str(surf_depths[j]) + ' m to ' + str(surf_depths[j+1]) + ' m'
                    data_list.append(data)
                
        # Add inital-type input data
        for i, file_path in enumerate(initial_input_paths):
            # Load data
            df, units = load_initial_data(file_path)
            # Store name, adapt for FABM initials
            name = os.path.splitext(os.path.basename(file_path))[0]
            if name[len(name)-8:] == '_initial':
                name = 'InitialConditions'
            # Set column names
            df.columns = ['Depth'] + [(name + '_' + df.columns[j+1]) for j in range(len(df.columns)-1)]
            # Aggregate duplicate rows
            df = df.groupby('Depth', as_index=False).sum()
            # Store every column as DataArray
            for j in range(len(df.columns)-1):
                data = df.set_index(['Depth']).to_xarray()[df.columns[j+1]]
                if (len(units) > 0) and (units.iloc[0,j] != '-'):
                    data.attrs['units'] = normalize_units(units.iloc[0,j])
                data.attrs['title'] = 'Initial Condition'
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
                if (len(units) > 0) and (units.iloc[0,j] != '-'):
                    data.attrs['units'] = normalize_units(units.iloc[0,j])
                data.attrs['title'] = 'Forcing Input'
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

# Helper function

def normalize_units(units):
    """
    Convert units to standard format
    
    Parameters:
    -----------
    units : str
        Units in format m/s2/kg or m s-2 kg-1
    
    Returns:
    -----------
    units_normalized : str
        Units in format m/s2/kg
    """

    # Return - if no units provided
    if not units or units.strip() == '':
        return '-'

    # Split into seperate parts
    parts = units.split()

    # Empty lists for numerator and denominator
    num = []
    den = []

    # Add to numerator or denominator list depending on exponent
    for part in parts:
        # match: base + exponent
        match = re.fullmatch(r'([a-zA-Z]+)([-]?\d+)?', part)
        if not match:
            num.append(part)
        else:
            base = match.group(1)
            exp = int(match.group(2)) if match.group(2) else 1
            if exp > 0:
                num.append(base + (str(exp) if exp > 1 else ''))
            elif exp < 0:
                den.append(base + (str(-exp) if exp < -1 else ''))

    # 1/unit if no numerator
    if not num:
        num_str = '1'
    else:
        num_str = '  '.join(num)

    # Return with denominator if existing
    if not den:
        return num_str
    else:
        return num_str + '/' + '/'.join(den)