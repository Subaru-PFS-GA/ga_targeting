import re
import pandas as pd
import numpy as np
from collections.abc import Iterable

def pd_append_column(df, name, data, dtype=None):
    """
    Append a column to a DataFrame by allowing to specify the data type.
    """

    if isinstance(data, str):
        df[name] = pd.Series([data] * len(df), dtype=dtype)
    elif not isinstance(data, Iterable):
        df[name] = pd.Series([data] * len(df), dtype=dtype)
    elif isinstance(data, pd.Series):
        #df[name] = data.reset_index(drop=True)
        df[name] = data.astype(dtype)
    else:
        df[name] = pd.Series(data, dtype=dtype)

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
