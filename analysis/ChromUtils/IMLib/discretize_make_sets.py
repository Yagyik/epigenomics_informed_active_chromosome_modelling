import numpy as np
import pickle as pkl
import pandas as pd
import os

import json
import matplotlib.pyplot as plt
import itertools
import upsetplot
import matplotlib.pyplot as plt
from upsetplot import from_memberships

import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from scipy.stats import norm, poisson


def discretize_columns(df, column_suffix, n_bins, bin_edges=None):
    """
    Discretizes columns in the dataframe that start with "as" and have the same quantity after the "_".

    Parameters:
    - df: pd.DataFrame - The input dataframe
    - column_suffix: str - The suffix of the columns to be discretized (e.g., "ReturnTime")
    - n_bins: int - The number of bins to use for discretization

    Returns:
    - pd.DataFrame - The dataframe with discretized columns
    """
    # Identify columns to be discretized
    columns_to_discretize = [col for col in df.columns if col.endswith(column_suffix)]
    print(columns_to_discretize)
    # Discretize each column and create a new column with discretized scores
    for col in columns_to_discretize:
        new_col_name = f"{col}_discretized"
        if bin_edges is None:
            df[new_col_name] = pd.cut(df[col], bins=n_bins-1, labels=False, include_lowest=True).astype(int)
        else:
            df[new_col_name] = pd.cut(df[col], bins=bin_edges, labels=False, include_lowest=True).astype(int)
    return df

def discretize_columns_logarithmically(df, column_suffix, n_bins, direction="forward"):
    """
    Discretizes columns in the dataframe logarithmically.

    Parameters:
    - df: pd.DataFrame - The input dataframe
    - column_suffix: str - The suffix of the columns to be discretized (e.g., "ReturnTime")
    - n_bins: int - The number of bins to use for discretization
    - direction: str - "forward" for largest interval last, "backward" for largest interval first

    Returns:
    - pd.DataFrame - The dataframe with discretized columns
    """
    # Identify columns to be discretized
    columns_to_discretize = [col for col in df.columns if col.endswith(column_suffix)]
    
    # Discretize each column and create a new column with discretized scores
    for col in columns_to_discretize:
        new_col_name = f"{col}_discretized"
        min_val = df[col].min()
        max_val = df[col].max()
        shift = 0.01 * (max_val - min_val)
        min_val += shift
        max_val += shift

        print(min_val, max_val)

        if direction == "forward":
            bins = np.logspace(np.log10(shift), np.log10(max_val - min_val), n_bins-1)
        elif direction == "backward":
            bins = np.logspace(np.log10(max_val - min_val), np.log10(shift), n_bins-1)[::-1]
            print(np.logspace(np.log10(max_val - min_val), np.log10(shift), n_bins-1))
            print("bins",bins)
        else:
            raise ValueError("Invalid direction value. Choose from 'forward' or 'backward'.")

        
        # Clip values to ensure they fall within the bin range
        clipped_values = np.clip(df[col] - min_val, shift, max_val - min_val)
        
        print(np.min(clipped_values), np.max(clipped_values))
        # Discretize the values
        print(clipped_values)
        df[new_col_name] = pd.cut(clipped_values, bins=bins, labels=False, include_lowest=True).astype(int)
    
    return df


def filterNan(df,col,filter_to="max"):
    if filter_to == "max":
        df[col].fillna(df[col].max(), inplace=True)
    elif filter_to == "min":
        df[col].fillna(df[col].min(), inplace=True)
    elif filter_to == "zero":
        df[col].fillna(0, inplace=True)
    else:
        raise ValueError("Invalid filter_to value. Choose from 'max', 'min', or 'zero'.")
    return df


def filter_allNan(df):
    discretized_df = df.copy()
    for col in discretized_df.columns:
        if col.startswith("as") and "MaxIntermingling" in col or "InteractionStrength" in col:
            discretized_df = filterNan(discretized_df, col, filter_to="zero")

    return discretized_df


def create_sets_info(df, col1, col2,nbins=4):
    if nbins == 4:
        string_scale = ["1-low", "2-medium", "3-high"]
    if nbins == 6:
        string_scale = ["1-very-low","2-low","3-medium","4-high","5-very-high"]
    if nbins != 4 and nbins != 6:
        print("bad nbins value")
        return None
    print(col1,col2,col1.split('_'),col2.split('_'))
    # in this version, col2 always has the form "InteractionStrength_discretized"
    possible_sets = [col1.split('_')[1] + "-" + string_scale[i] + "," + col2.split('_')[0] + "-" + string_scale[j] for i in range(nbins-1) for j in range(nbins-1)]
    
    print(possible_sets)
    print("col lengths",len(df[col1].values),len(df[col2].values))
    # print("lens uniques",len(string_scale),np.unique(df[col1].values),np.unique(df[col2].values))
    if col1 not in df.columns or col2 not in df.columns:
        print(f"Columns in df: {df.columns}")
        return None
    sets_list = [col1.split('_')[1] + "-" + string_scale[df[col1].values[i]] + "," + col2.split('_')[0] + "-" + string_scale[df[col2].values[i]]
                 for i in range(len(df[col1].values))]
    print("sets list",len(sets_list))
    df['sets_info'] = pd.Categorical(np.array(sets_list), categories=possible_sets, ordered=False)
    for pset in possible_sets:
        if pset not in sets_list:
            new_row = {col1.split('_')[1]: np.nan, col2.split('_')[0]: np.nan, 'sets_info': pset}
            df = df._append(new_row, ignore_index=True)
            # df['sets_info'] = pd.Categorical(df['sets_info'], categories=possible_sets, ordered=True)
    df.sort_values(by=[col1, col2], ascending=[True, True], inplace=True)
    # df = df.sort_values('sets_info')
    df.reset_index(drop=True, inplace=True)
    
    return df,possible_sets

def create_single_set_info(inp_df, col,nbins=4):
    df = inp_df.copy()
    if nbins == 4:
        string_scale = ["1-low", "2-medium", "3-high"]
    if nbins == 6:
        string_scale = ["1-very-low","2-low","3-medium","4-high","5-very-high"]
    if nbins != 4 and nbins != 6:
        print("bad nbins value")
        return None
    str_snip = col.split('_')[1] #.split('_')[0]
    possible_sets = [str_snip + "-" + string_scale[i] for i in range(nbins-1)]
    print(possible_sets)
    
    sets_list = [str_snip+ "-" + string_scale[df[col].values[i]] for i in range(len(df[col].values))]
    print(sets_list)
    df['sis_'+col] = pd.Categorical(np.array(sets_list), categories=possible_sets, ordered=True)
    # df = df.sort_values('sis_'+col)
    df = df.sort_values(col)
    df.reset_index(drop=True, inplace=True)
    
    return df, possible_sets





def create_sets_info_twocol(df, col1, col2,nbins1=3,nbins2=4):
    if nbins1 == 2:
        string_scale1 = ["1-low", "2-high"]
    elif nbins1 == 3:
        string_scale1 = ["1-low", "2-medium", "3-high"]
    elif nbins1 == 5:
        string_scale1 = ["1-very-low", "2-low", "3-medium", "4-high", "5-very-high"]
    else:
        print("bad nbins1 value")
        return None

    if nbins2 == 2:
        string_scale2 = ["1-low", "2-high"]
    elif nbins2 == 3:
        string_scale2 = ["1-low", "2-medium", "3-high"]
    elif nbins2 == 5:
        string_scale2 = ["1-very-low", "2-low", "3-medium", "4-high", "5-very-high"]
    else:
        print("bad nbins2 value")
        return None

    print(col1, col2, col1.split('_'), col2.split('_'))
    possible_sets = [col1.split('_')[1] + "-" + string_scale1[i] + "," + col2.split('_')[1] + "-" + string_scale2[j] for i in range(nbins1) for j in range(nbins2)]
    
    print(possible_sets)
    print("col lengths", len(df[col1].values), len(df[col2].values))
    if col1 not in df.columns or col2 not in df.columns:
        print(f"Columns in df: {df.columns}")
        return None
    sets_list = [col1.split('_')[1] + "-" + string_scale1[int(df[col1].values[i])] + "," + col2.split('_')[1] + "-" + string_scale2[int(df[col2].values[i])]
                 for i in range(len(df[col1].values))]
    print("sets list", len(sets_list))
    df['sets_info'] = pd.Categorical(np.array(sets_list), categories=possible_sets, ordered=False)
    for pset in possible_sets:
        if pset not in sets_list:
            new_row = {col1.split('_')[1]: np.nan, col2.split('_')[0]: np.nan, 'sets_info': pset}
            df = df._append(new_row, ignore_index=True)
    df.sort_values(by=[col1, col2], ascending=[True, True], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    return df, possible_sets

def create_sets_info_twocol_meta(df, col1, col2,nbins1=3,nbins2=4):
    if nbins1 == 2:
        string_scale1 = ["1-low", "2-high"]
    elif nbins1 == 3:
        string_scale1 = ["1-low", "2-medium", "3-high"]
    elif nbins1 == 5:
        string_scale1 = ["1-very-low", "2-low", "3-medium", "4-high", "5-very-high"]
    else:
        print("bad nbins1 value")
        return None

    if nbins2 == 2:
        string_scale2 = ["1-low", "2-high"]
    elif nbins2 == 3:
        string_scale2 = ["1-low", "2-medium", "3-high"]
    elif nbins2 == 5:
        string_scale2 = ["1-very-low", "2-low", "3-medium", "4-high", "5-very-high"]
    else:
        print("bad nbins2 value")
        return None

    print(col1, col2, col1.split('_'), col2.split('_'))
    possible_sets = [col1.split('_')[1] + "-" + string_scale1[i] + "," + col2.split('_')[0] + "-" + string_scale2[j] for i in range(nbins1) for j in range(nbins2)]
    
    print(possible_sets)
    print("col lengths", len(df[col1].values), len(df[col2].values))
    if col1 not in df.columns or col2 not in df.columns:
        print(f"Columns in df: {df.columns}")
        return None
    sets_list = [col1.split('_')[1] + "-" + string_scale1[int(df[col1].values[i])] + "," + col2.split('_')[0] + "-" + string_scale2[int(df[col2].values[i])]
                 for i in range(len(df[col1].values))]
    print("sets list", len(sets_list))
    df['sets_info'] = pd.Categorical(np.array(sets_list), categories=possible_sets, ordered=False)
    for pset in possible_sets:
        if pset not in sets_list:
            new_row = {col1.split('_')[1]: np.nan, col2.split('_')[0]: np.nan, 'sets_info': pset}
            df = df._append(new_row, ignore_index=True)
    df.sort_values(by=[col1, col2], ascending=[True, True], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    return df, possible_sets






def create_possible_sets_single_col(col1, nbins1=3):
    if nbins1 == 2:
        string_scale1 = ["1-low", "2-high"]
    elif nbins1 == 3:
        string_scale1 = ["1-low", "2-medium", "3-high"]
    elif nbins1 == 5:
        string_scale1 = ["1-very-low", "2-low", "3-medium", "4-high", "5-very-high"]
    else:
        print("bad nbins1 value")
        return None

    print(col1, col1.split('_'))
    possible_sets = [col1 + "-" + string_scale1[i] for i in range(nbins1)]
    
    print(possible_sets)
    return possible_sets