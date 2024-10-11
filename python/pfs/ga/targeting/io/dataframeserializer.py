import os
import pandas as pd
from pandas.io.parsers import read_csv
from pandas.io.feather_format import to_feather, read_feather

class DataFrameSerializer():
    """
    Implements functions to better serialize DataFrames into various formats.
    """

    # TODO: add pickle, HDF5 and FITS formats
    
    def __init__(self, format=None, compression=None, version=None, orig=None):
        if not isinstance(orig, DataFrameSerializer):
            self.__format = format
            self.__compression = compression
            self.__version = version
        else:
            self.__compression = orig.__compression
            self.__format = orig.__format
            self.__version = orig.__version

    def __get_format(self) -> str:
        return self.__format
    
    def __set_format(self, format: str):
        self.__format = format

    format = property(__get_format, __set_format)

    def __get_format_functions(self, filename: str, format: str = None):
        format = format if format is not None else self.__format

        if format is None:
            # Get the file extension from the file name
            fn, format = os.path.splitext(filename)

        if format == ".csv":
            read_func = read_csv
            write_func = lambda df, fn, **kwargs: df.to_csv(fn, **kwargs)
        elif format == ".feather":
            read_func = read_feather
            write_func = to_feather
        else:
            raise ValueError(f"Unknown file extension: {format}")
        
        return read_func, write_func

    def read(self, filename: str, **kwargs) -> pd.DataFrame:
        read_func, _ = self.__get_format_functions(filename)
        return read_func(filename, **kwargs)
    
    def write(self, assignments: pd.DataFrame, filename: str, **kwargs):
        _, write_func = self.__get_format_functions(filename)
        write_func(assignments, filename, **kwargs)