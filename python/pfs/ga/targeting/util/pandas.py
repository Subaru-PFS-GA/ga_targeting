import re
import pandas as pd
import numpy as np
from collections.abc import Iterable

def pd_append_column(df, name, data, dtype=None):
    """
    Append a column to a DataFrame by allowing to specify the data type.
    """

    if isinstance(data, str) or not isinstance(data, Iterable):
        df[name] = pd.Series([data] * len(df), index=df.index, dtype=dtype)
    elif isinstance(data, pd.DataFrame):
        df[name] = data[data.columns[0]].astype(dtype).values
    elif isinstance(data, pd.Series):
        # Ignore the index of data by converting to an array using .values
        df[name] = data.astype(dtype).values
    else:
        df[name] = pd.Series(data, index=df.index, dtype=dtype)

def pd_update_column(df, loc, name, data, dtype=None):
    """
    Update a column in a DataFrame by allowing to specify the update locations
    and the target data type.
    """

    if isinstance(data, str) or not isinstance(data, Iterable):
        df.loc[loc, name] = pd.Series([data] * len(df), index=df.index, dtype=dtype)
    elif isinstance(data, pd.Series):
        # Ignore the index of data by converting to an array using .values
        df.loc[loc, name] = data.astype(dtype).values
    else:
        df.loc[loc, name] = pd.Series(data, index=df.index, dtype=dtype)

def pd_to_nullable(df: pd.DataFrame, columns=None, in_place=False) -> pd.DataFrame:
    """
    Convert integer columns so that they accept Nan values.
    """

    if not in_place:
        df = df.copy()

    for c in columns if columns is not None else df:
        if df[c].dtype == np.int32:
            df[c] = df[c].astype('Int32')
        elif df[c].dtype == np.int64:
            df[c] = df[c].astype('Int64')
        elif df[c].dtype == np.float32:
            df[c] = df[c].astype('Float32')
        elif df[c].dtype == np.float64:
            df[c] = df[c].astype('Float64')

    return df

def pd_null_to_nan(df: pd.DataFrame, columns=None, in_place=False) -> pd.DataFrame:
    """
    Convert None values of float columns to NaN.
    """

    if not in_place:
        df = df.copy()

    for c in columns if columns is not None else df:
        if df[c].dtype == 'Float32':
            df[c] = df[c].fillna(np.nan).astype('float32')
        if df[c].dtype == 'Float64':
            df[c] = df[c].fillna(np.nan).astype('float64')

    return df

def pd_fillna(df: pd.DataFrame, name, value):
    if df.dtypes[name] == float or df.dtypes[name] == np.float64 or df.dtypes[name] == np.float32:
        df.loc[df[name].isna(), name] = value
    elif df.dtypes[name] == pd.Float64Dtype() or df.dtypes[name] == pd.Float32Dtype():
        df[name] = df[name].fillna(value)
    elif df.dtypes[name] == 'string':
        df[name] = df[name].fillna('')
    elif df.dtypes[name] == object:
        pass
    else:
        pass