import pandas as pd
import numpy as np
from collections.abc import Iterable

def pd_append_column(df, c, d, dtype=None):
    """
    Append a column to a DataFrame by allowing to specify the data type.
    """

    if isinstance(d, str):
        df[c] = pd.Series([d] * len(df), dtype=dtype)
    elif not isinstance(d, Iterable):
        df[c] = pd.Series([d] * len(df), dtype=dtype)
    elif isinstance(d, pd.Series):
        df[c] = d.reset_index(drop=True)
    else:
        df[c] = pd.Series(d, dtype=dtype)

def pd_to_nullable(df: pd.DataFrame, in_place=False) -> pd.DataFrame:
    """
    Convert integer columns so that they accept Nan values.
    """

    if not in_place:
        df = df.copy()

    for c in df:
        if df[c].dtype == np.int32:
            df[c] = df[c].astype('Int32')
        elif df[c].dtype == np.int64:
            df[c] = df[c].astype('Int64')
        elif df[c].dtype == np.float32:
            df[c] = df[c].astype('Float32')
        elif df[c].dtype == np.float64:
            df[c] = df[c].astype('Float64')

    return df
