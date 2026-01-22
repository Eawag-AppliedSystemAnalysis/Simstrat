# Load libraries

import pandas as pd
import numpy as np
import os



# Load data 

def load_forcing_data(path, file, start='file_start', stop='file_stop', load_units=True):
    """
    Smart function to load forcing data from local file.
    
    Parameters:
    -----------
    path : str
        Path to local data file. Can be:
        - DAT file: 'data/geneva/Forcing.dat'
        - Directory: 'data/geneva' (will search for file in that directory)
    file : str
        The name of the file to load
    start : str, optional
        Start datetime in format 'YYYYMMDDHHMM'
        (default: file_start -> start from beginning of file)
    stop : str, optional
        Stop datetime in format 'YYYYMMDDHHMM'
        (default: file_stop -> stop at end of file)
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

    # Check if the file is of .dat format
    # If path is a directory create a path to the file
    file_ext = os.path.splitext(file)[1].lower()
    if file_ext == '.dat':
        if os.path.isdir(path):
            file_path = os.path.join(path, file)
            if os.path.exists(file_path):
                local_file_path = file_path
            else:
                raise ValueError(f"{file} not found in: {path}")
        elif os.path.isfile(path):
            local_file_path = path
        else:
            raise ValueError(f"{file} not found in: {path}")
    else:
        raise ValueError(f"Unsupported file format: {file_ext}. Use .dat")
    
    print(f"Loading data from local file: {local_file_path}")

    # First load variables with units
    data_title = pd.read_csv(local_file_path, sep=r'\s+', nrows=0)

    # Create new header without time and units
    variable_columns = data_title.columns[2::2].tolist()

    # Get units of all variables (except time) into a dataframe
    if load_units:
        units_columns = np.asarray(data_title.columns[3::2].tolist())
        # Strip .x identifiers added by pandas for duplicate units
        for i in np.arange(1, len(units_columns)+1).astype('str'):
            units_columns = np.char.rstrip(units_columns, ' .{}'.format(i))
        # Strip surrounding brackets
        units_columns = np.char.rstrip(units_columns, ']')
        units_columns = np.char.lstrip(units_columns, '[')
        units = pd.DataFrame([units_columns], columns = variable_columns, index = ['Unit'])
    else:
        units = pd.DataFrame([])

    # Add Datetime variable
    variable_columns.insert(0, 'Datetime')

    # Load data with new variable header
    data = pd.read_csv(local_file_path, sep=r'\s+', skiprows=1, names = variable_columns).astype('float')
    
    # Convert first column (days since reference data) to datetime
    # Assuming reference date is 1981-01-01 (common for Simstrat)
    time_values = data['Datetime']
    datetime_values = pd.Timestamp('1981-01-01') + pd.to_timedelta(time_values, unit='D').astype('timedelta64[s]')
    data['Datetime'] = datetime_values
    
    print(f"✓ Loaded {len(data)} rows from {local_file_path}")
    
    # Filter by date range
    if start == 'file_start':
        start_dt = data['Datetime'][0]
    else:
        start_dt = pd.to_datetime(start, format='%Y%m%d%H%M')
    if stop == 'file_stop':
        stop_dt = data['Datetime'][len(data)-1]
    else:
        stop_dt = pd.to_datetime(stop, format='%Y%m%d%H%M')
    data = data[(data['Datetime'] >= start_dt) & (data['Datetime'] <= stop_dt)]
    
    print(f"✓ Filtered to date range: {len(data)} rows between {start_dt.date()} and {stop_dt.date()}")
    
    return data, units

def load_inflow_data(path, file, start='file_start', stop='file_stop', load_units=True):
    """
    Smart function to load inflow data from local file.
    
    Parameters:
    -----------
    path : str
        Path to local data file. Can be:
        - DAT file: 'data/geneva/Qin.dat', 'data/geneva/Sin.dat', 'data/geneva/Tin.dat'
        - Directory: 'data/geneva' (will search for file in that directory)
    file : str
        The name of the file to load
    start : str, optional
        Start datetime in format 'YYYYMMDDHHMM'
        (default: file_start -> start from beginning of file)
    stop : str, optional
        Stop datetime in format 'YYYYMMDDHHMM'
        (default: file_stop -> stop at end of file)
    load_units : bool, optional
        Whether to load units into a data frame or not
        (default: True -> load units)        
    
    Returns:
    --------
    deep : pd.DataFrame
        DataFrame with deep inflow data
    surf : pd.DataFrame
        DataFrame with surface inflow data
    deep_order : list
        The original order of deep inflows
    surf_order : list
        The original order of surface inflows
    units : pd.DataFrame
        DataFrame with units of inflows (empty if load_units is False)
    """

    # Check if the file is of .dat format
    # If path is a directory create a path to the file
    file_ext = os.path.splitext(file)[1].lower()
    if file_ext == '.dat':
        if os.path.isdir(path):
            file_path = os.path.join(path, file)
            if os.path.exists(file_path):
                local_file_path = file_path
            else:
                raise ValueError(f"{file} not found in: {path}")
        elif os.path.isfile(path):
            local_file_path = path
        else:
            raise ValueError(f"{file} not found in: {path}")
    else:
        raise ValueError(f"Unsupported file format: {file_ext}. Use .dat")
    
    print(f"Loading data from local file: {local_file_path}")
    
    # First load inflow units into a dataframe
    if load_units:
        data_title = pd.read_csv(local_file_path, nrows=0).columns[0]
        # Strip surrounding brackets
        units_columns = np.char.strip(np.asarray([s.split(']') for s in data_title.split('[')][1:])[1:,0])
        # Strip .x identifiers added by pandas for duplicate units
        for i in np.arange(1, len(units_columns)+1).astype('str'):
            units_columns = np.char.rstrip(units_columns, ' .{}'.format(i))
        variable_columns = np.char.strip(np.asarray([s.split(']') for s in data_title.split('[')][1:])[:-1,1])
        units = pd.DataFrame([units_columns], columns = variable_columns, index = ['Unit'])
    else:
        units = pd.DataFrame([])
    
    # Get the amount of inflows and their depths
    # Depths contains two arrays: deep inflow depths and surface inflow depths
    amount_columns = np.asarray(pd.read_csv(local_file_path, sep=r'\s+', skiprows = 1, nrows=0).columns).astype('int')
    deep_n = amount_columns[0]
    surf_n = amount_columns[1]
    depths = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 2, nrows=1, names = range(deep_n + surf_n)).values[0]
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
        deep = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols = deep_usecols, names = deep_columns).astype('float')
        print(f"✓ Loaded {deep_n} deep inflows with {len(deep)} rows from {local_file_path}")
    else:
        deep = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols =[0], names = ['Datetime']).astype('float')
        print(f"No deep inflow")
    if surf_n > 0:
        surf_columns = np.array(['{0} ({1})'.format(a, b) for a, b in zip(depths[1], surf_order)]).astype(np.dtype('<U32'))
        surf_columns = np.insert(surf_columns, 0, 'Datetime')
        surf_usecols = np.concatenate(([0], np.arange(deep_n + 1, deep_n + surf_n + 1)))
        surf = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols = surf_usecols, names = surf_columns).astype('float')
        print(f"✓ Loaded {surf_n} surface inflows with {len(surf)} rows from {local_file_path}")
    else:
        surf = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols =[0], names = ['Datetime']).astype('float')
        print(f"No surface inflow")

    # Convert first column (days since reference data) to datetime
    # Assuming reference date is 1981-01-01 (common for Simstrat)
    time_values = deep['Datetime']
    datetime_values = pd.Timestamp('1981-01-01') + pd.to_timedelta(time_values, unit='D').astype('timedelta64[s]')
    deep['Datetime'] = datetime_values
    surf['Datetime'] = datetime_values

    # Filter by date range
    if start == 'file_start':
        start_dt = deep['Datetime'][0]
    else:
        start_dt = pd.to_datetime(start, format='%Y%m%d%H%M')
    if stop == 'file_stop':
        stop_dt = deep['Datetime'][len(deep)-1]
    else:
        stop_dt = pd.to_datetime(stop, format='%Y%m%d%H%M')
    deep = deep[(deep['Datetime'] >= start_dt) & (deep['Datetime'] <= stop_dt)]
    surf = surf[(surf['Datetime'] >= start_dt) & (surf['Datetime'] <= stop_dt)]
    
    print(f"✓ Filtered to date range: {len(deep)} rows between {start_dt.date()} and {stop_dt.date()}")

    return deep, surf, deep_order, surf_order, units



# Save data

def save_forcing_data(data, path, file, copy=True):
    """
    Smart function to save edited forcing data to local file.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame with modified forcing data
    path : str
        Path to local data file. Can be:
        - DAT file: 'data/geneva/Forcing.dat'
        - Directory: 'data/geneva' (will search for file in that directory)
    file : str
        The name of the file to load
    copy : bool, optional
        Whether to create a new file instead of replacing the existing one
        (default: True -> create a new file)
    """

    # Check if the file is of .dat format
    # If path is a directory create a path to the file
    file_ext = os.path.splitext(file)[1].lower()
    if file_ext == '.dat':
        if os.path.isdir(path):
            file_path = os.path.join(path, file)
            if os.path.exists(file_path):
                local_file_path = file_path
            else:
                raise ValueError(f"{file} not found in: {path}")
        elif os.path.isfile(path):
            local_file_path = path
        else:
            raise ValueError(f"{file} not found in: {path}")
    else:
        raise ValueError(f"Unsupported file format: {file_ext}. Use .dat")

    # Load original data with unitless variable header from modified data
    variable_columns = data.columns
    data_original = pd.read_csv(local_file_path, sep=r'\s+', skiprows=1, names = variable_columns).astype('float')
    data_original_title = (pd.read_csv(local_file_path, nrows=0)).columns[0]
    
    # Convert first column (datetime) to days since reference date
    # Assuming reference date is 1981-01-01 (common for Simstrat)
    # Also assuming time values are initially given with 4 decimal place (common for Simstrat)
    datetime_values = data['Datetime']
    time_values = round((datetime_values - pd.Timestamp('1981-01-01')) / pd.Timedelta(days=1), 4)
    start = time_values.min()
    stop = time_values.max()
    
    # Replace original with new data in edited time frame
    for variable in variable_columns[1:]:
        data_original.loc[(data_original['Datetime'] >= start) & (data_original['Datetime'] <= stop), variable] = data[variable]

    # Get path to a copy or the existing file to write to
    if copy == True:
        write_path = '{0}_{1}.{2}'.format(os.path.splitext(file_path)[0], 'modified', 'dat')
    else:
        write_path = local_file_path
        
    print(f"Writing modified data into: {write_path}")

    # Write edited data
    with open(write_path, 'w', encoding='utf-8') as f:
        f.write(data_original_title + '\n')
        for i in range(len(data_original['Datetime'])):
            if any(np.isnan([data_original[c][i] for c in variable_columns])):
                print(f"NaN value in row {i}. Row skipped.")
                continue
            f.write(' '.join(['%10.4f' % data_original[c][i] for c in variable_columns]) + '\n')
            
    print(f"✓ Modified {len(data)} rows")

def save_inflow_data(deep, surf, deep_order, surf_order, path, file, copy=True):
    """
    Smart function to save edited inflow to local file.
    
    Parameters:
    -----------
    deep : pd.DataFrame
        DataFrame with modified deep inflow data
    surf : pd.DataFrame
        DataFrame with modified surface inflow data
    deep_order : list
        The modified order of deep inflows
    surf_order : list
        The modified order of surface inflows
    path : str
        Path to local data file. Can be:
        - DAT file: 'data/geneva/Qin.dat', 'data/geneva/Sin.dat', 'data/geneva/Tin.dat'
        - Directory: 'data/geneva' (will search for file in that directory)
    file : str
        The name of the file to load
    copy : bool, optional
        Whether to create a new file instead of replacing the existing one
        (default: True -> create a new file)
    """

    # Check if the file is of .dat format
    # If path is a directory create a path to the file
    file_ext = os.path.splitext(file)[1].lower()
    if file_ext == '.dat':
        if os.path.isdir(path):
            file_path = os.path.join(path, file)
            if os.path.exists(file_path):
                local_file_path = file_path
            else:
                raise ValueError(f"{file} not found in: {path}")
        elif os.path.isfile(path):
            local_file_path = path
        else:
            raise ValueError(f"{file} not found in: {path}")
    else:
        raise ValueError(f"Unsupported file format: {file_ext}. Use .dat")
    
    # Get new inflow depths for deep inflows
    deep_depths = np.asarray(deep.columns[1:].to_list())
    # Strip indices
    for i in range(1, len(deep_depths) + 1):
        deep_depths = np.char.rstrip(deep_depths, ' ({})'.format(i))
    deep_depths = deep_depths.astype('float')
    
    # Get new inflow depths for surface inflows
    surf_depths = np.asarray(surf.columns[1:].to_list())
    # Strip indices
    for i in range(1, len(surf_depths) + 1):
        surf_depths = np.char.rstrip(surf_depths, ' ({})'.format(i))
    surf_depths = surf_depths.astype('float')
    
    # Load original deep and surface inflows
    amount_columns = np.asarray(pd.read_csv(local_file_path, sep=r'\s+', skiprows = 1, nrows=0).columns).astype('int')
    deep_original_n = amount_columns[0]
    surf_original_n = amount_columns[1]
    if deep_original_n > 0:
        deep_columns = [str(i) for i in range(1, deep_original_n+1)]
        deep_columns.insert(0, 'Datetime')
        deep_usecols = np.arange(deep_original_n + 1)
        deep_original = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols = deep_usecols, names = deep_columns).astype('float')
    else:
        deep_original = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols =[0], names = ['Datetime']).astype('float')
    if surf_original_n > 0:
        surf_original_n = max(surf_order)
        surf_columns = [str(i) for i in range(1, surf_original_n+1)]
        surf_columns.insert(0, 'Datetime')
        surf_usecols = np.concatenate(([0], np.arange(deep_original_n, deep_original_n + surf_original_n) + 1))
        surf_original = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols = surf_usecols, names = surf_columns).astype('float')
    else:
        surf_original = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols =[0], names = ['Datetime']).astype('float')
    
    # Convert first column (datetime) to days since reference date
    # Assuming reference date is 1981-01-01 (common for Simstrat)
    # Also assuming time values are initially given with 4 decimal place (common for Simstrat)
    datetime_values = deep['Datetime']
    time_values = round((datetime_values - pd.Timestamp('1981-01-01')) / pd.Timedelta(days=1), 4)
    start = time_values.min()
    stop = time_values.max()
    
    # Construct new data frame with deep and surface inflows, where:
    # New inflows are zero where not defined in new inflow data
    # Existing inflows get original values where not defined in new inflow data
    data = pd.read_csv(local_file_path, sep=r'\s+', skiprows = 3, usecols =[0], names = ['Datetime']).astype('float')
    for i, col_original in enumerate(deep_order):
        if col_original == 0:
            pre_start = np.zeros(len(deep_original[deep_original['Datetime'] < start]))
            post_stop = np.zeros((len(deep_original[deep_original['Datetime'] > stop])))
            inflow = np.concatenate((pre_start, deep.iloc[:, i+1].values, post_stop))
        else:
            pre_start = deep_original[deep_original['Datetime'] < start][str(col_original)]
            post_stop = deep_original[deep_original['Datetime'] > stop][str(col_original)]
            inflow = np.concatenate((pre_start, deep.iloc[:, i+1].values, post_stop))
        data.insert(i+1, str(i+1), inflow)
    for i, col_original in enumerate(surf_order):
        if col_original == 0:
            pre_start = np.zeros(len(surf_original[surf_original['Datetime'] < start]))
            post_stop = np.zeros((len(surf_original[surf_original['Datetime'] > stop])))
            inflow = np.concatenate((pre_start, surf.iloc[:, i+1].values, post_stop))
        else:
            pre_start = surf_original[surf_original['Datetime'] < start][str(col_original)]
            post_stop = surf_original[surf_original['Datetime'] > stop][str(col_original)]
            inflow = np.concatenate((pre_start, surf.iloc[:, i+1].values, post_stop))
        data.insert(i+len(deep_order)+1, str(i+len(deep_order)+1), inflow)
    
    # Load original title
    data_original_title = (pd.read_csv(local_file_path, nrows=0)).columns[0]
    
    # Get path to a copy or the existing file to write to
    if copy == True:
        write_path = '{0}_{1}.{2}'.format(os.path.splitext(file_path)[0], 'modified', 'dat')
    else:
        write_path = local_file_path
    print(f"Writing modified data into: {write_path}")
    
    # Write edited data
    with open(write_path, 'w', encoding='utf-8') as f:
        f.write(data_original_title + '\n')
        f.write('%10d %10d\n' % (len(deep_depths), len(surf_depths)))
        f.write('-1         ' + ' '.join(['%10.2f' % z for z in deep_depths]) + ' ' + ' '.join(['%10.2f' % z for z in surf_depths]) + '\n')
        for i in range(len(data)):
            if any(np.isnan([data[c][i] for c in data.columns[1:]])):
                print(f"NaN value in row {i}. Row skipped.")
                continue
            f.write('%10.4f ' % data['Datetime'][i])
            f.write(' '.join(['%10.2f' % data[c][i] for c in data.columns[1:]]))
            f.write('\n')
            
    print(f"✓ Modified {len(deep)} rows")



# Modify data

def check_negative(val_mod, val_save, var):
    """
    Function to check for negative values after modification.
    
    Parameters:
    -----------
    val_mod : pd.Series
        Series with modified data
    val_save : pd.Series
        Series with original data
    var: str
        Variable name or inflow index
    
    Returns:
    --------
    val : pd.Series
        Series:
        - val_save if val_mod contains negative values
        - val_mod else
    """

    if np.any(val_mod < 0):
        if var.isdigit():
            print(f"Negative values for inflow ({var}) — no modification")
        else:
            print(f"Negative values for forcing variable '{var}' — no modification")
        val = val_save
    else:
        if var.isdigit():
            print(f"✓ Inflow ({var}) modified")
        else:
            print(f"✓ Forcing variable '{var}' modified")
        val = val_mod

    return val

def change_depth(data, col_index, depth):
    """
    Function to change depth of one inflow.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame with inflow data
    col_index : int
        Index of inflow to edit
    depth : float
        New depth for inflow     
    
    Returns:
    --------
    data_ed : pd.DataFrame
        DataFrame with edited inflow data
    """
    
    # Check if col_index is inside bounds
    if col_index == -1:
        col_index = len(data.columns) - 1
    elif col_index <= 0 or col_index >= len(data.columns):
        raise ValueError(f"index must have an existing inflow")

    # Change depth at col_index
    col_name = data.columns[col_index]
    new_name = str(float(depth)) + ' ({})'.format(col_index)
    data_ed = data.rename(columns={col_name: new_name})

    print(f"✓ Modified depth of inflow ({col_index})")
    
    return data_ed

def add_inflow(data, data_order, col_index=-1, values=0, depth=0, pos=False):
    """
    Function to add one inflow.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame with inflow data
    data_order : list
        The order of inflows
    col_index : int, optional
        Index of inflow to remove (default: -1 -> at last position)
    values : int or 1D-array of length len(data), optional
        Inflow values (default: 0)
    depth : float, optional
        Depth of added inflow (default: 0)
    pos: bool, optional
        Whether the inflow cannot be negative
        (default: False -> negative values are ok)
    
    Returns:
    --------
    data_ed : pd.DataFrame
        DataFrame with edited inflow data
    data_ed_order : list
        The edited order of inflows
    """

    # Check if col_index is inside bounds
    if col_index == -1:
        col_index = len(data.columns)
    elif col_index > len(data.columns):
        raise ValueError(f"col_index cannot be larger than amount of existing inflows + 1 ({len(data.columns)})")
    elif col_index <= 0:
        raise ValueError(f"col_index must be larger than 0, or -1")
    
    # Create the inflow if it is of correct format
    if pos and np.any(values < 0):
        raise ValueError(f"values cannot be negative or contain negative numbers")
    if np.ndim(values) == 0:
        inflow = np.full((len(data)), float(values)) 
    elif np.ndim(values) == 1:
        if len(values) == len(data):
            inflow = values.astype('float')
        else:
            raise ValueError(f"values must be a scalar or 1D array of same length as data")
    else:
        raise ValueError(f"values must be a scalar or a 1D array of same length as data")

    # Insert the created inflow in a copy of data
    data_ed = data.copy(True)
    data_ed.insert(col_index, str(float(depth)) + ' (0)', inflow)

    # Adapt inflow indices in brackets
    data_ed_names = np.asarray(data_ed.columns[1:].to_list())
    for i in range(0, len(data_ed.columns) + 1):
        data_ed_names = np.char.rstrip(data_ed_names, '({})'.format(i))
    data_ed_names = np.array(['{0} ({1})'.format(a, b) for a, b in zip(data_ed_names, np.arange(1, len(data_ed.columns)).astype('str'))]).astype(np.dtype('<U32'))
    data_ed_names = np.insert(data_ed_names, 0, 'Datetime')
    data_ed.columns = data_ed_names

    # Adapt order list
    data_ed_order = data_order.copy()
    data_ed_order.insert(col_index - 1, 0)

    print(f"✓ Added inflow at ({col_index})")

    return data_ed, data_ed_order

def remove_inflow(data, data_order, col_index=-1):
    """
    Function to remove one inflow.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame with inflow data
    data_order : list
        The order of inflows
    col_index : int, optional
        Index of inflow to remove (default: -1 -> at last position)
        
    
    Returns:
    --------
    data_ed : pd.DataFrame
        DataFrame with edited inflow data
    data_ed_order : list
        The edited order of inflows
    """

    # Check if col_index is inside bounds
    if col_index == -1:
        col_index = len(data.columns) - 1
    elif col_index <= 0 or col_index >= len(data.columns):
        raise ValueError(f"index must have an existing inflow")
    
    # Drop inflow
    data_ed = data.drop(data.columns[col_index], axis = 1)
    
    # Adapt inflow indices in brackets
    if len(data_ed.columns) > 1:
        data_ed_names = np.asarray(data_ed.columns[1:].to_list())
        for i in range(1, len(data.columns) + 1):
            data_ed_names = np.char.rstrip(data_ed_names, '({})'.format(i))
        data_ed_names = np.array(['{0} ({1})'.format(a, b) for a, b in zip(data_ed_names, np.arange(1, len(data_ed.columns)).astype('str'))]).astype(np.dtype('<U32'))
        data_ed_names = np.insert(data_ed_names, 0, 'Datetime')
        data_ed.columns = data_ed_names

    # Adapt order list
    data_ed_order = data_order.copy()
    data_ed_order.pop(col_index - 1)

    print(f"✓ Removed inflow ({col_index})")

    return data_ed, data_ed_order

def change_order(data, data_order, new_col_index):
    """
    Function to change order of inflows.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame with inflow data
    data_order : list
        The order of inflows
    new_col_index : tuple with indices of all inflows
        New ordering of inflows
    
    Returns:
    --------
    data_ed : pd.DataFrame
        DataFrame with edited inflow data
    data_ed_order : list
        The edited order of inflows
    """

    # Check if new inflow order is of correct format
    if not len(np.unique(np.asarray(new_col_index))) == len(data.columns) - 1:
        raise ValueError(f"Need one unique index for every inflow. Currently there are {len(data.columns[1:])} inflows.")
    elif not np.all(np.unique(np.asarray(new_col_index)) == np.linspace(1,len(data.columns)-1,len(data.columns)-1)):
        raise ValueError(f"Need one unique index for every inflow. Currently there are {len(data.columns[1:])} inflows.")

    # Change order of inflow
    new_columns = [data.columns[i] for i in new_col_index]
    new_columns.insert(0, 'Datetime')
    data_ed = data.reindex(columns = new_columns)
    
    # Adapt inflow indices in brackets
    if len(data_ed.columns) > 1:
        data_ed_names = np.asarray(data_ed.columns[1:].to_list())
        for i in range(0, len(data_ed.columns) + 1):
            data_ed_names = np.char.rstrip(data_ed_names, '({})'.format(i))
        data_ed_names = np.array(['{0} ({1})'.format(a, b) for a, b in zip(data_ed_names, np.arange(1, len(data_ed.columns)).astype('str'))]).astype(np.dtype('<U32'))
        data_ed_names = np.insert(data_ed_names, 0, 'Datetime')
        data_ed.columns = data_ed_names

    # Adapt order list
    data_ed_order = [data_order[i] for i in np.array(new_col_index) - 1]

    print(f"✓ Changed inflow order")

    return data_ed, data_ed_order
