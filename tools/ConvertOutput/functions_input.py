# Load libraries

import pandas as pd
import numpy as np
import os

# Read Paths

def read_input_paths(path_to_input):
    """
    Load paths to input data
    
    Parameters:
    -----------
    path_to_input : str
        Path to input
    path_to_output : str
        Path to output      
    
    Returns:
    --------
    absorption_input_paths : str list
        Paths to absorption-type data
    inflow_input_paths : str list
        Paths to inflow-type data
    initial_input_paths : str list
        Paths to initial-type data
    forcing_input_paths : str list
        Paths to input-type data
    """

    # Get paths from input path
    input_paths = [os.path.join(path_to_input, f) for f in os.listdir(path_to_input) if (f[-4:] == '.dat')]

    # Initialise empty lists
    absorption_input_paths = []
    inflow_input_paths = []
    initial_input_paths = []
    forcing_input_paths = []

    # Read paths and sort them according to the first lines
    for i, file_path in enumerate(input_paths):
        # Don't decode header if it's not necessary
        header = b''
        first = []
        second = []
        # Store header and read first two line
        with open(file_path, 'rb') as f:
            i = 0
            for line in f:
                if i == 0:
                    header = line
                if i == 1:
                    line = line.decode('ascii')
                    first = np.fromstring(line, sep=' ')
                if i == 2:
                    line = line.decode('ascii')
                    second = np.fromstring(line, sep=' ')
                i = i + 1
        # Skip one-line or one-column input
        if len(second) < 2:
            continue
        # Absorption-type input (only one type)
        elif len(first) == 1:
            absorption_input_paths.append(file_path)
        elif len(first) == len(second):
            # Inflow-type input (if the above is fulfilled there is only one inflow: deep or surface)
            if (len(first) == 2) and (((first[0] == 1) and (first[1] == 0)) or ((first[0] == 0) and (first[1] == 1))):
                inflow_input_paths.append(file_path)
            # Forcing- and Initial-type need header of certain format
            header = header.decode('ascii', errors='ignore').split()
            if len(first) * 2 != len(header):
                continue
            # Forcing-type input
            elif (header[1] == '[d]'):
                forcing_input_paths.append(file_path)
            # Inital-type input (Bathymetry or InitialConditions)
            elif (header[1] == '[m]'):
                initial_input_paths.append(file_path)
            # Skip if header is of unknown format
            else:
                continue
        # Inflow-type input (first row sum is total amount of inflows)
        elif sum(first) == len(second) - 1:
            inflow_input_paths.append(file_path)
        # Skip if input is of unknown format
        else:
            continue

    return absorption_input_paths, inflow_input_paths, initial_input_paths, forcing_input_paths

# Load data

def load_absorption_data(file_path, load_units=True):
    """
    Load absorption data from input path
    
    Parameters:
    -----------
    file_path : str
        Path to local data file  
    load_units : bool, optional
        Whether to load units into a data frame or not
        (default: True -> load units)        
    
    Returns:
    --------
    data : pd.DataFrame
        DataFrame with absorption data
    units : pd.DataFrame
        DataFrame with units of absorption (empty if load_units is False)
    """
    
    # First load inflow units into a dataframe
    if load_units:
        try:
            data_title = pd.read_csv(file_path, nrows=0).columns[0]
            # Strip surrounding brackets
            units_columns = np.char.strip(np.asarray([s.split(']') for s in data_title.split('[')][1:])[1:,0])
            # Strip .x identifiers added by pandas for duplicate units
            for i in np.arange(1, len(units_columns)+1).astype('str'):
                units_columns = np.array([s[:-len(f'.{i}')] if s.endswith(f'.{i}') else s for s in units_columns])
            variable_columns = ['Depth', 'Absorption']
            units = pd.DataFrame([units_columns], columns = variable_columns, index = ['Unit'])
        except (UnicodeDecodeError, ValueError, IndexError):
            units = pd.DataFrame([])
    else:
        units = pd.DataFrame([])
    
    # Get the amount of absorptions and their depths
    amount_columns = np.asarray(pd.read_csv(file_path, sep=r'\s+', skiprows = 1, nrows=0).columns).astype('int')
    n = amount_columns[0]
    depths = pd.read_csv(file_path, sep=r'\s+', skiprows = 2, nrows=1, names = range(n)).values[0]
    depths = depths.astype('float')
    depths = np.split(depths, [n])

    # Save absorption order
    order = list(range(1,n+1))

    # Load absorptions
    if n > 0:
        data_columns = np.array(['{0} ({1})'.format(a, b) for a, b in zip(depths[0], order)]).astype(np.dtype('<U32'))
        data_columns = np.insert(data_columns, 0, 'Datetime')
        data_usecols = np.arange(n + 1)
        data = pd.read_csv(file_path, sep=r'\s+', skiprows = 3, usecols = data_usecols, names = data_columns).astype('float')
    else:
        data = pd.read_csv(file_path, sep=r'\s+', skiprows = 3, usecols =[0], names = ['Datetime']).astype('float')

    return data, units.iloc[:,1:]

def load_inflow_data(file_path, load_units=True):
    """
    Load inflow data from input path
    
    Parameters:
    -----------
    file_path : str
        Path to local data file
    load_units : bool, optional
        Whether to load units into a data frame or not
        (default: True -> load units)   
    
    Returns:
    --------
    deep : pd.DataFrame
        DataFrame with deep inflow data
    surf : pd.DataFrame
        DataFrame with surface inflow data
    units : pd.DataFrame
        DataFrame with units of inflows (empty if load_units is False)
    """
    
    # First load inflow units into a dataframe
    if load_units:
        try:
            data_title = pd.read_csv(file_path, nrows=0).columns[0]
            # Strip surrounding brackets
            units_columns = np.char.strip(np.asarray([s.split(']') for s in data_title.split('[')][1:])[1:,0])
            # Strip .x identifiers added by pandas for duplicate units
            for i in np.arange(1, len(units_columns)+1).astype('str'):
                units_columns = np.array([s[:-len(f'.{i}')] if s.endswith(f'.{i}') else s for s in units_columns])
            variable_columns = np.char.strip(np.asarray([s.split(']') for s in data_title.split('[')][1:])[:-1,1])
            units = pd.DataFrame([units_columns], columns = variable_columns, index = ['Unit'])
        except (UnicodeDecodeError, ValueError, IndexError):
            units = pd.DataFrame([])
    else:
        units = pd.DataFrame([])
    
    # Get the amount of inflows and their depths
    # Depths contains two arrays: deep inflow depths and surface inflow depths
    amount_columns = np.asarray(pd.read_csv(file_path, sep=r'\s+', skiprows = 1, nrows=0).columns).astype('int')
    deep_n = amount_columns[0]
    surf_n = amount_columns[1]
    depths = pd.read_csv(file_path, sep=r'\s+', skiprows = 2, nrows=1, names = range(deep_n + surf_n)).values[0]
    depths = depths.astype('float')
    depths = np.split(depths, [deep_n])

    # Save deep and surface inflow order
    deep_order = list(range(1,deep_n+1))
    surf_order = list(range(1,surf_n+1))

    # Load deep and surface inflows
    if deep_n > 0:
        deep_columns = np.array(['{0} ({1})'.format(a, b) for a, b in zip(depths[0], deep_order)]).astype(np.dtype('<U32'))
        deep_columns = np.insert(deep_columns, 0, 'Datetime')
        deep_usecols = np.arange(deep_n + 1)
        deep = pd.read_csv(file_path, sep=r'\s+', skiprows = 3, usecols = deep_usecols, names = deep_columns).astype('float')
    else:
        deep = pd.read_csv(file_path, sep=r'\s+', skiprows = 3, usecols =[0], names = ['Datetime']).astype('float')
    if surf_n > 0:
        surf_columns = np.array(['{0} ({1})'.format(a, b) for a, b in zip(depths[1], surf_order)]).astype(np.dtype('<U32'))
        surf_columns = np.insert(surf_columns, 0, 'Datetime')
        surf_usecols = np.concatenate(([0], np.arange(deep_n + 1, deep_n + surf_n + 1)))
        surf = pd.read_csv(file_path, sep=r'\s+', skiprows = 3, usecols = surf_usecols, names = surf_columns).astype('float')
    else:
        surf = pd.read_csv(file_path, sep=r'\s+', skiprows = 3, usecols =[0], names = ['Datetime']).astype('float')

    return deep, surf, units.iloc[:,1:]

def load_initial_data(file_path, load_units=True):
    """
    Load initial data from input path
    
    Parameters:
    -----------
    file_path : str
        Path to local data file
    load_units : bool, optional
        Whether to load units into a data frame or not
        (default: True -> load units)        
    
    Returns:
    --------
    data : pd.DataFrame
        DataFrame with forcing data
    units : pd.DataFrame
        DataFrame with units of variables (empty if load_units is False)
    """

    # First load variables with units
    data_title = pd.read_csv(file_path, sep=r'\s+', nrows=0)

    # Create new header without time and units
    variable_columns = data_title.columns[2::2].tolist()

    # Get units of all variables (except time) into a dataframe
    if load_units:
        try:
            units_columns = np.asarray(data_title.columns[3::2].tolist())
            # Strip .x identifiers added by pandas for duplicate units
            for i in np.arange(1, len(units_columns)+1).astype('str'):
                units_columns = np.array([s[:-len(f'.{i}')] if s.endswith(f'.{i}') else s for s in units_columns])
            # Strip surrounding brackets
            units_columns = np.char.rstrip(units_columns, ']')
            units_columns = np.char.lstrip(units_columns, '[')
            units = pd.DataFrame([units_columns], columns = variable_columns, index = ['Unit'])
        except (UnicodeDecodeError, ValueError, IndexError):
            units = pd.DataFrame([])
    else:
        units = pd.DataFrame([])

    # Add Depth variable
    variable_columns.insert(0, 'Depth')

    # Load data with new variable header
    data = pd.read_csv(file_path, sep=r'\s+', skiprows=1, names = variable_columns).astype('float')

    # Order from bottom to surface
    if data.iloc[0,0] > data.iloc[-1,0]:
        data = data.iloc[::-1].reset_index(drop=True)
    
    return data, units

def load_forcing_data(file_path, load_units=True):
    """
    Load forcing data from input path
    
    Parameters:
    -----------
    file_path : str
        Path to local data file
    load_units : bool, optional
        Whether to load units into a data frame or not
        (default: True -> load units)        
    
    Returns:
    --------
    data : pd.DataFrame
        DataFrame with forcing data
    units : pd.DataFrame
        DataFrame with units of variables (empty if load_units is False)
    """

    # First load variables with units
    data_title = pd.read_csv(file_path, sep=r'\s+', nrows=0)
    

    # Create new header without time and units
    variable_columns = data_title.columns[2::2].tolist()

    # Get units of all variables (except time) into a dataframe
    if load_units:
        try:
            units_columns = np.asarray(data_title.columns[3::2].tolist())
            # Strip .x identifiers added by pandas for duplicate units
            for i in np.arange(1, len(units_columns)+1).astype('str'):
                units_columns = np.array([s[:-len(f'.{i}')] if s.endswith(f'.{i}') else s for s in units_columns])
            # Strip surrounding brackets
            units_columns = np.char.rstrip(units_columns, ']')
            units_columns = np.char.lstrip(units_columns, '[')
            units = pd.DataFrame([units_columns], columns = variable_columns, index = ['Unit'])
        except (UnicodeDecodeError, ValueError, IndexError):
            units = pd.DataFrame([])
    else:
        units = pd.DataFrame([])

    # Add Datetime variable
    variable_columns.insert(0, 'Datetime')

    # Load data with new variable header
    data = pd.read_csv(file_path, sep=r'\s+', skiprows=1, names = variable_columns).astype('float')
    
    return data, units