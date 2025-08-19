import os
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as feather
import h5py
from astropy.table import Table

from ..util import safe_deep_copy

class DataFrameSerializer():
    """
    Implements functions to better serialize DataFrames into various formats.

    The DataFrameSerializer class is a base class for reading and writing observations.

    Variables
    ---------
    dataset : str
        Dataset, within a container file (eg. HDF5 group or FITS extension)
    columns : list
        List of columns to be read/written.
    column_names : list
        Column names, can be used if column name cannot be inferred, otherwise
        override the column names.
    column_map : dict
        Column name mapping dictionary.
    data_types : dict
        Data types for each column, if cannot be inferred. Otherwise, override
        the data types.
    value_map: dict of dict
        Map values in columns. First key is the column name, second key is the
        value to map, value is the new value. If value is not found in the map,
        the original value in the column is used.
    index : str
        Column to use as the index.
    idcol : str
        When not None, a column with this name will be created as unique integers.
    mask : function, mask or index
        mask function, mask array or index to apply to the DataFrame.
    format : str
        Data file format, if cannot be inferred from file extension. Otherwise,
        override the file format.
    kwargs : dict
        Additional keyword arguments passed to the file format read/write functions.
    compression : str
        Compression algorithm, when supported.
    version : str
        Data format version, reserved for future use.
    """

    # TODO: add pickle, numpy and FITS formats?
    # TODO: implement mappings, etc when writing to a file, not just when reading
    
    def __init__(self,
                 dataset=None,
                 columns=None,
                 column_names=None,
                 column_map=None,
                 data_types=None,
                 value_map=None,
                 index=None,
                 idcol=None,
                 mask=None,
                 format=None,
                 kwargs=None,
                 compression=None,
                 version=None,
                 orig=None):

        if not isinstance(orig, DataFrameSerializer):
            self.__dataset = dataset                            # Dataset, within a container file
            self.__columns = columns                            # List of columns to be read/written
            self.__column_names = column_names                  # Column names, if cannot be inferred
            self.__column_map = column_map                      # Column name mapping dictionary
            self.__data_types = data_types                      # Data types for each column
            self.__value_map = value_map                        # Map values in columns
            self.__index = index                                # Column to use as the index
            self.__idcol = idcol                                # When not None, a column with this name will be created as unique integers
            self.__mask = mask                                  # Filter function, mask or index
            self.__format = format                              # Data file format, can be inferred from file extension
            self.__kwargs = kwargs                              # Additional keyword arguments for read/write functions
            self.__compression = compression                    # Compression algorithm, when supported
            self.__version = version                            # Data format version, reserved for future use
        else:
            self.__dataset = dataset if dataset is not None else orig.__dataset
            self.__columns = columns if columns is not None else safe_deep_copy(orig.__columns)
            self.__column_names = column_names if column_names is not None else safe_deep_copy(orig.__column_names)
            self.__column_map = column_map if column_map is not None else safe_deep_copy(orig.__column_map)
            self.__data_types = data_types if data_types is not None else safe_deep_copy(orig.__data_types)
            self.__value_map = value_map if value_map is not None else safe_deep_copy(orig.__value_map)
            self.__index = index if index is not None else safe_deep_copy(orig.__index)
            self.__idcol = idcol if idcol is not None else orig.__idcol
            self.__mask = mask if mask is not None else orig.__mask
            self.__format = format if format is not None else orig.__format
            self.__kwargs = kwargs if kwargs is not None else safe_deep_copy(orig.__kwargs)
            self.__compression = compression if compression is not None else orig.__compression
            self.__version = version if version is not None else orig.__version

    #region Properties

    def __get_columns(self):
        return self.__columns
    
    def __set_columns(self, columns: list):
        self.__columns = columns

    columns = property(__get_columns, __set_columns)

    def __get_column_names(self):
        return self.__column_names
    
    def __set_column_names(self, column_names: list):
        self.__column_names = column_names

    column_names = property(__get_column_names, __set_column_names)

    def __get_column_map(self):
        return self.__column_map
    
    def __set_column_map(self, column_map: dict):
        self.__column_map = column_map

    column_map = property(__get_column_map, __set_column_map)

    def __get_data_types(self):
        return self.__data_types
    
    def __set_data_types(self, data_types: dict):
        self.__data_types = data_types

    data_types = property(__get_data_types, __set_data_types)

    def __get_value_map(self):
        return self.__value_map
    
    def __set_value_map(self, value_map: dict):
        self.__value_map = value_map

    value_map = property(__get_value_map, __set_value_map)

    def __get_index(self):
        return self.__index
    
    def __set_index(self, index):
        self.__index = index

    index = property(__get_index, __set_index)

    def __get_idcol(self):
        return self.__idcol

    def __set_idcol(self, idcol):
        self.__idcol = idcol

    idcol = property(__get_idcol, __set_idcol)

    def __get_mask(self):
        return self.__mask
    
    def __set_mask(self, mask):
        self.__mask = mask

    mask = property(__get_mask, __set_mask)

    def __get_format(self) -> str:
        return self.__format
    
    def __set_format(self, format: str):
        self.__format = format

    format = property(__get_format, __set_format)

    def __get_kwargs(self):
        return self.__kwargs
    
    def __set_kwargs(self, kwargs: dict):
        self.__kwargs = kwargs

    kwargs = property(__get_kwargs, __set_kwargs)

    def __get_compression(self):
        return self.__compression
    
    def __set_compression(self, compression):
        self.__compression = compression

    compression = property(__get_compression, __set_compression)

    def __get_version(self):
        return self.__version
    
    def __set_version(self, version):
        self.__version = version

    version = property(__get_version, __set_version)

    #endregion

    def __get_format(self, filename: str, format: str = None) -> str:
        """
        Determine the file format based on the file extension. If the format is
        not provided, the file extension is used to determine the format, otherwise
        the provided format is used.
        """

        format = format if format is not None else self.__format

        if format is None:
            # TODO: recognize if end with .gz and return the format
            #       accordingly. Also support other compressions
            _, format = os.path.splitext(filename)
            
        return format

    def __get_format_functions(self, filename: str, format: str = None):
        """
        Return the read and write functions based on the file format.
        """

        format = self.__get_format(filename, format)

        if format == ".csv":
            read_func = self.__read_csv
            write_func = self.__write_csv
        elif format == '.ecsv':
            read_func = self.__read_ecsv
            write_func = self.__write_ecsv
        elif format == ".feather":
            read_func = self.__read_feather
            write_func = self.__write_feather
        elif format == ".h5" or format == ".hdf5":
            read_func = self.__read_hdf5
            write_func = self.__write_hdf5
        elif format == '.fit' or format == '.fits':
            read_func = self.__read_fits
            write_func = self.__write_fits
        else:
            raise ValueError(f"Unknown file extension: {format}")
        
        return read_func, write_func
    
    def __apply_mask(self, df: pd.DataFrame, mask):
        """
        Apply a mask to the DataFrame. The filter can either be a function, mask or index.
        
        Filtering is done on a row-by-row basis, use the `columns` property to filter columns.
        """

        if mask is None:
            return df
        elif callable(mask):
            return df[mask(df)]
        elif isinstance(mask, pd.Series):
            return df[mask]
        elif isinstance(mask, pd.Index):
            # TODO: verify this
            raise NotImplementedError()
            return df.loc[mask]
        elif isinstance(mask, np.ndarray):
            return df[mask]
        elif isinstance(mask, slice):
            return df[mask]
        else:
            raise NotImplementedError()
        
    def __apply_mappings(self, df: pd.DataFrame):
        """
        Apply column mappings and data types to the DataFrame.
        """

        # Rename columns in the order they appear
        if self.__column_names is not None:
            mapping = { old: new  for old, new in zip(df.columns, self.__column_names) }
            df = df.rename(mapping)

        # Rename columns using a mapping dictionary
        if self.__column_map is not None:
            df = df.rename(columns=self.__column_map)

        # Apply the value map
        if self.__value_map is not None:
            df = df.replace(self.__value_map)

        # Convert data types
        if self.__data_types is not None:
            df = df.astype(self.__data_types)

        return df

    def __apply_index(self, df: pd.DataFrame):
        """
        Apply an index to the DataFrame.
        """

        if self.__index is not None:
            df = df.set_index(self.__index)
        else:
            df.reset_index(drop=True, inplace=True)
        
        return df

    def __apply_idcol(self, df: pd.DataFrame):
        """
        Generate a unique integer ID column if `idcol` is set. We intentionally don't use the
        index here because that might be composite or not unique, etc.
        """

        if self.__idcol is not None:
            if self.__idcol in df.columns:
                raise ValueError(f"Column '{self.__idcol}' already exists in the DataFrame.")
            df[self.__idcol] = np.arange(len(df), dtype=int)
        
        return df
    
    def __read_csv(self, filename: str, dataset=None, mask=None, **kwargs) -> pd.DataFrame:
        if self.__column_names is not None:
            kwargs['names'] = self.__column_names
            kwargs['header'] = None

        if self.__columns is not None:
            kwargs['usecols'] = self.__columns

        if self.data_types is not None:
            kwargs['dtype'] = self.__data_types

        df = pd.read_csv(filename, index_col=False, **kwargs)
        df = self.__apply_mask(df, mask)
        return df
    
    def __write_csv(self, df: pd.DataFrame, filename: str, dataset=None, mask=None, **kwargs):
        if self.__columns is not None:
            kwargs['columns'] = self.__columns

        df = self.__apply_mask(df, mask)
        df.to_csv(filename, index=False, **kwargs)

    def __read_feather(self, filename: str, dataset=None, mask=None, **kwargs) -> pd.DataFrame:
        df = pd.read_feather(filename, columns=self.__columns, **kwargs)
        df = self.__apply_mask(df, mask)
        return df
    
    def __write_feather(self, df: pd.DataFrame, filename: str, dataset=None, mask=None, **kwargs):
        df = self.__apply_mask(df, mask)
        df.to_feather(filename, **kwargs)

    def __read_hdf5(self, filename: str, dataset=None, mask=None, **kwargs) -> pd.DataFrame:
        # TODO: apply filter when reading from HDF5 if filter type is supported

        with h5py.File(filename, 'r') as h:
            if dataset is not None:
                g = h[dataset]
            else:
                g = h
            d = { k: g[k][()] for k in g.keys() }
            
        df = pd.DataFrame(d)
        df = self.__apply_mask(df, mask)
        return df
    
    def __write_hdf5(self, df: pd.DataFrame, filename: str, dataset=None, mask=None, **kwargs):
        df = self.__apply_mask(df, mask)
        with h5py.File(filename, 'a') as h:
            # Recreate group if exists
            if dataset is not None and dataset in h:
                del h[dataset]
                g = h.create_group(dataset)
            else:
                g = h

            for col in df.columns:
                if col in g:
                    del g[col]
                
                g.create_dataset(col, data=np.array(df[col]))

    def __read_fits(self, filename: str, dataset=None, mask=None, **kwargs) -> pd.DataFrame:
        df = Table.read(filename, **kwargs).to_pandas()
        df = self.__apply_mask(df, mask)
        return df

    def __write_fits(self, df: pd.DataFrame, filename: str, dataset=None, mask=None, **kwargs):
        raise NotImplementedError()

    def __read_ecsv(self, filename: str, dataset=None, mask=None, **kwargs) -> pd.DataFrame:
        df = Table.read(filename).to_pandas()
        df = self.__apply_mask(df, mask)
        return df

    def __write_ecsv(self, df: pd.DataFrame, filename: str, dataset=None, mask=None, **kwargs):
        raise NotImplementedError()

    def read(self, filename: str, dataset=None, format=None, mask=None, **kwargs) -> pd.DataFrame:
        dataset = dataset if dataset is not None else self.__dataset
        mask = mask if mask is not None else self.__mask
        if self.__kwargs is not None:
            kwargs.update(self.__kwargs)

        read_func, _ = self.__get_format_functions(filename, format)
        df = read_func(filename, dataset=dataset, mask=mask, **kwargs)
        df = self.__apply_mappings(df)
        df = self.__apply_index(df)
        df = self.__apply_idcol(df)
        return df
    
    def write(self, df: pd.DataFrame, filename: str, dataset=None, mask=None, **kwargs):
        dataset = dataset if dataset is not None else self.__dataset
        mask = mask if mask is not None else self.__mask
        if self.__kwargs is not None:
            kwargs.update(self.__kwargs)

        _, write_func = self.__get_format_functions(filename)
        df = self.__apply_mappings(df)
        write_func(df, filename, dataset=dataset, mask=mask, **kwargs)
