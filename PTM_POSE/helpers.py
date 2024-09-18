import pandas as pd
import numpy as np


def join_unique_entries(x, sep = ';'):
    """
    For use with groupby, combines all unique entries separated by ';', removing any NaN entries
    """
    #check if only nan entries
    if all(i != i for i in x):
        return np.nan
    if any(sep in i for i in x if i == i): #check if ';' already in entry, if so, split and remove any NaN entries
        split_list = [i.split(sep) for i in x if i == i]
        #split entries in list by ';' and flatten list
        flat_list = [item for sublist in split_list for item in sublist]
        return sep.join(set(flat_list))
    else:
        entry_list = [str(i) for i in x if i == i]
        return sep.join(set(entry_list))

def join_entries(x, sep = ';'):
    """
    For use with groupby, combines all entries separated by ';', removing any NaN entries
    """
    #check if only nan entries
    if all(i != i for i in x):
        return np.nan

    if any(sep in i for i in x if i == i): #check if ';' already in entry, if so, split and remove any NaN entries
        split_list = [i.split(sep) for i in x if i == i]
        #split entries in list by ';' and flatten list
        flat_list = [item for sublist in split_list for item in sublist]
        return sep.join(flat_list)

    else:
        entry_list = [str(i) for i in x if i == i]
        return sep.join(entry_list)

def join_except_self(df, group_col, value_col, new_col, sep = ';'):
    """
    For a given dataframe, combines all entries with the same information except for the current row, adds that to the new_col label, and returns the updated dataframe
    
    Parameters
    ----------
    df: pandas DataFrame
        The dataframe to be updated
    group_col: str
        The column to group the dataframe by
    value_col: str
        The column to be combined
    new_col: str
        The new column to be added to the dataframe with the grouped information (excluding the info from the current row)

    Returns
    -------
    df: pandas DataFrame
        updated dataframe with new col labeled with new_col value
    """
    df = df.copy()
    df[new_col] = df.groupby(group_col)[value_col].transform(join_unique_entries, sep)

    #go through each row and remove the value(s) in the new column that is in the value column
    new_values = []
    for i, row in df.iterrows():
        if row[new_col] == row[new_col] and row[value_col] == row[value_col]:
            new_values.append(';'.join([trans for trans in row[new_col].split(sep) if trans not in row[value_col].split(sep)]))
        elif row[value_col] != row[value_col]:
            new_values.append(row[new_col])
        else:
            new_values.append(np.nan)
    df[new_col] = new_values
    return df