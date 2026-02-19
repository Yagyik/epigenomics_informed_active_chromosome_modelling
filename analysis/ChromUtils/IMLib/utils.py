import numpy as np
import pickle as pkl
import json
import itertools

def load_pickled(all_paths,test_decomp_series, test_HiC_series):
    decomp_series_list = np.empty((0, *test_decomp_series.shape))
    HiC_series_list = np.empty((0, *test_HiC_series.shape))


    for path in all_paths:
        print(path)
        filename = path[1] + "/HiC_dir/pickleDecomp"
        fileObj = open(filename, 'rb')
        decomp_series = pkl.load(fileObj)
        decomp_series_list = np.vstack([decomp_series_list, [decomp_series]])
        fileObj.close()
        
        filename = path[1] + "/HiC_dir/pickleAllHiC"
        fileObj = open(filename, 'rb')
        HiC_series = pkl.load(fileObj)
        print(HiC_series.shape)
        HiC_series_list = np.vstack([HiC_series_list, [HiC_series]])
        fileObj.close()
        print("loaded")
    return decomp_series_list, HiC_series_list


def organise_by(combined_str, key_value_pairs):
    # Find indices where key value pairs are present
    print("combined_str",combined_str)
    print("key_value_pairs",key_value_pairs)
    indices = [i for i, d in enumerate(combined_str) if all(d.get(k) == v for k, v in key_value_pairs.items())]
    print("indices",indices)
    # indices = [i for i, d in enumerate(combined_str) if all(d.get(k) == v for k, v in key_value_pairs)]
    
    # Identify variable key value pairs from the subset
    subset = [combined_str[i] for i in indices]
    print("subset",subset)
    print("subset keys",subset[0].keys())  
    print("key subset", [k for k in subset[0].keys() if k not in dict(key_value_pairs)])
    unique_values = {k: list(dict.fromkeys(d[k] for d in subset)) for k in subset[0].keys() if k not in dict(key_value_pairs)}
    print("unique vals",unique_values)
    # Generate combined list of key[i1]value[j1]_key[i2]value[j2]
    combined_list = []
    for combination in itertools.product(*unique_values.values()):
        print(combination)
        combined_str = "_".join(f"{k}{v}" for k, v in zip(unique_values.keys(), combination))
        combined_list.append(combined_str)
    
    return indices, combined_list