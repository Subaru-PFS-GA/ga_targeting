#!/usr/bin/env python3
"""
Script to read all observation log CSV files from 
/home/enk/pfs/spt_ssp_observation/runs/*/obslog/*.csv
and combine them into a single dataframe.
"""

import pandas as pd
from pathlib import Path
from typing import List
import glob


def find_obslog_csvs(base_path: str = "/home/enk/pfs/spt_ssp_observation/runs") -> List[Path]:
    """
    Find all CSV files matching the pattern runs/*/obslog/*.csv
    
    Parameters
    ----------
    base_path : str
        Base directory containing the runs folders
        
    Returns
    -------
    List[Path]
        List of paths to all CSV files found
    """
    pattern = f"{base_path}/*/obslog/*.csv"
    csv_files = glob.glob(pattern)
    csv_files = sorted([Path(f) for f in csv_files])
    print(f"Found {len(csv_files)} CSV files")
    return csv_files


def read_obslog_csv(csv_path: Path) -> pd.DataFrame:
    """
    Read a single observation log CSV file.
    
    Parameters
    ----------
    csv_path : Path
        Path to the CSV file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with the CSV contents and additional metadata columns
    """
    # Read CSV without treating # as comment to preserve header lines starting with #
    df = pd.read_csv(csv_path, comment=None)
    
    # Strip leading # from column names if present
    df.columns = df.columns.str.lstrip('#').str.strip()
    
    # Add source information
    df['source_file'] = csv_path.name
    df['run_id'] = csv_path.parent.parent.name

    return df


def read_all_obslog_csvs(base_path: str = "/home/enk/pfs/spt_ssp_observation/runs") -> pd.DataFrame:
    """
    Read all observation log CSV files and combine them into a single dataframe.
    
    Parameters
    ----------
    base_path : str
        Base directory containing the runs folders
        
    Returns
    -------
    pd.DataFrame
        Combined dataframe with all observation logs
    """
    csv_files = find_obslog_csvs(base_path)
    
    if not csv_files:
        print("No CSV files found!")
        return pd.DataFrame()
    
    # Read all CSV files
    dfs = []
    for csv_file in csv_files:
        try:
            df = read_obslog_csv(csv_file)
            dfs.append(df)
            #print(f"  Read {csv_file.name}: {len(df)} rows")
            #print(f"    Columns: {list(df.columns)}")
        except Exception as e:
            print(f"  Error reading {csv_file}: {e}")
            continue
    
    # Combine all dataframes
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        print(f"\nTotal rows: {len(combined_df)}")
        print(f"Columns: {list(combined_df.columns)}")
        return combined_df
    else:
        print("No data read successfully!")
        return pd.DataFrame()


def main():
    """Main function to demonstrate usage."""
    # Read all CSV files
    df = read_all_obslog_csvs()
    
    # Display summary
    if not df.empty:
        print("\nDataFrame Info:")
        print(df.info())
        print("\nFirst few rows:")
        print(df.head())
        print("\nRun IDs found:")
        print(df['run_id'].unique())
    
    return df


if __name__ == "__main__":
    df = main()
