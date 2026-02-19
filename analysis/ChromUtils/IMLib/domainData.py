import numpy as np
import pandas as pd
import pickle as pkl
import json
import os
from ChromUtils.conf_reader import parse_file

class DomainData:
    def __init__(self, combined_str_len, args,params,setlist_len = 1):
        self.decomp_series_list = np.empty((combined_str_len, setlist_len), dtype=object)
        self.HiC_series_list = np.empty((combined_str_len, setlist_len), dtype=object)
        ### include parameters of analysis here
        self.args = args
        self.params = params

    def add_decomp_series(self, i, j, decomp_series):
        self.decomp_series_list[i, j] = decomp_series

    def add_HiC_series(self, i, j, HiC_series):
        self.HiC_series_list[i, j] = HiC_series

    def get_decomp_series(self, i, j):
        return self.decomp_series_list[i, j]

    def get_HiC_series(self, i, j):
        return self.HiC_series_list[i, j]
    
    def aggregate_decomp_series_sets(self, i):
        return np.array([self.decomp_series_list[i, j] for j in range(self.decomp_series_list.shape[1])])
    
    def aggregate_HiC_series_sets(self, i):
        return np.array([self.HiC_series_list[i, j] for j in range(self.HiC_series_list.shape[1])])

    def average_decomp_series(self, indices):
        selected_series = [self.decomp_series_list[i, j] for i, j in indices]
        return np.mean(selected_series, axis=0)

    def average_HiC_series(self, indices):
        selected_series = [self.HiC_series_list[i, j] for i, j in indices]
        return np.mean(selected_series, axis=0)
    
    def add_multiple(self,indices,all_paths):
        # Add series to the object
        for i in range(len(indices)):
            path = all_paths[indices[i]]
            print(path)

            # for j, path in enumerate(all_paths):
                # print(path)
            filename = path[1] + "/HiC_dir/pickleDecomp"
            with open(filename, 'rb') as fileObj:
                decomp_series = pkl.load(fileObj)
                print(decomp_series.shape)
                # decomp_series_list[j, i] = decomp_series
            
            filename = path[1] + "/HiC_dir/pickleAllHiC"
            with open(filename, 'rb') as fileObj:
                HiC_series = pkl.load(fileObj)
                print(HiC_series.shape)
                # HiC_series_list[j, i] = HiC_series
            print("loaded")

            self.add_decomp_series(i, 0, decomp_series[:])
            self.add_HiC_series(i, 0, HiC_series[:])

            del HiC_series, decomp_series