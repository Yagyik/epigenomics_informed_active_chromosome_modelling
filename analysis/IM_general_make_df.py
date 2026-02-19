import numpy as np
import pickle as pkl
import pandas as pd
import os
import json
import argparse
import matplotlib.pyplot as plt
import itertools
import upsetplot
from upsetplot import from_memberships
from scipy.stats import norm, poisson, lognorm

import plotly.graph_objects as go
import plotly.express as px
# from ChromUtils.TSeriesLib import *
from ChromUtils import ChromLib, IMLib,utils,conf_reader,Structures
from ChromUtils.IMLib import file_handling
from ChromUtils.IMLib import utils as IMLib_utils
from ChromUtils.IMLib import df_combine as IM_df_combine

from ChromUtils.IMLib.domainData import DomainData


def load_analysis_params(analysis_params_path,mode="rejuv"):
    with open(analysis_params_path, 'r') as file:
        params= json.load(file)
    print(params)

    combined_str = []
    if mode == "rejuv":
        for eij in params['eijslist']:
            for a in params['aslist']:
                for st in params['stlist']:
                    for dsa in params['dsalist']:
                        for val in params['valist']:
                            if eij != 6.0:
                                combined_str.append({
                                    'sga': f'{dsa:.2f}',
                                    'epa': f'{eij:.1f}',
                                    'as': f'{val:.1f}',
                                    'sa': f'{a}',
                                    'tl': f'{st}'
                                })
                            else:
                                combined_str.append({
                                    'sga': f'{dsa:.2f}',
                                    'epa': f'{eij:.2f}',
                                    'as': f'{val:.1f}',
                                    'sa': f'{a}',
                                    'tl': f'{st}'
                                })

    elif mode == "dediff":
        for eij in params['eijslist']:
            for lmn in params['lmnlist']:
                for sgtau in params['sgtaulist']:
                    for dsa in params['dsalist']:
                        formatted_dsa = f'{dsa:.3g}'
                        for val in params['valist']:
                            if eij != 6.0:
                                combined_str.append({
                                    'sga': formatted_dsa,
                                    'epa': f'{eij:.1f}',
                                    'lmn': f'{lmn:.1f}',
                                    'sgat': f'{int(sgtau)}',
                                    'as': f'{val:.2g}' if val != 1.0 and val != 0.0 else f'{val:.1f}'
                                })

                            else:
                                combined_str.append({
                                    'sga': formatted_dsa,
                                    'epa': f'{eij:.2f}',
                                    'lmn': f'{lmn:.1f}',
                                    'sgat': f'{int(sgtau)}',
                                    'as': f'{val:.2g}' if val != 1.0 and val != 0.0 else f'{val:.1f}'
                                    })
        
    print(combined_str)
    return combined_str, params


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--analysis_params',type=str, help='path to analysis params json file')
    parser.add_argument('--mode', type=str, help='mode of analysis')
    parser.add_argument('--set', type=int, help='set to analyze')
    # parser.add_argument('--anaparprefix', type=str, help='prefix for analysis params')
    # parser.add_argument('--simparprefix', type=str, help='prefix for sim params')
    # parser.add_argument('--anaparpostfix', type=str, help='postfix for analysis params')
    # parser.add_argument('--simparpostfix', type=str, help='postfix for sim params')
    args = parser.parse_args()

    print(args.analysis_params)
    print(args.mode)

    
    combined_str,params = load_analysis_params(args.analysis_params,args.mode)
    print("combined_str",combined_str)
    
    for comb in combined_str:
        print(comb)
        for key, value in comb.items():
            print(key,value)
        print()

    inpref = params['inprefix']
    anaparprefix = params['anaparprefix']
    simparprefix = params['simparprefix']
    anaparpostfix = params['anaparpostfix']
    simparpostfix = params['simparpostfix']

    # anaparams_file = f'AnaParamFiles/anaparams_{setlist[1]}_{inprefix}_test_geom.json'
    # simparams_file = f'SimParamFiles/simparams_{setlist[1]}_{inprefix}_test_geom.json'

    # anaparams_filelist = [f'AnaParamFiles/anaparams_reprog_{set}_{params.inprefix}_geom.json' 
    #                       for i,set in enumerate(params.setlist)]
    # simparams_filelist = [f'SimParamFiles/simparams_reprog_{set}_{params.inprefix}_geom.json'
    #                        for i,set in enumerate(params.setlist)]
    # print(anaparams_filelist)
    # test_anaparams = conf_reader.parse_file(anaparams_filelist[0])
    # test_simparams = conf_reader.parse_file(simparams_filelist[0])

    # anaparams = conf_reader.parse_file(
    #     f'AnaParamFiles/anaparams_reprog_{args.set}_{inpref}_geom.json')
    # simparams = conf_reader.parse_file(
    #     f'SimParamFiles/simparams_reprog_{args.set}_{inpref}_geom.json')
    # paths = makepaths(test_anaparams["outprefix"],test_simparams["confdir"],combined_str[0])
    
    anaparams_file = f'{anaparprefix}_{args.set}_{inpref}_{anaparpostfix}.json'
    simparams_file = f'{simparprefix}_{args.set}_{inpref}_{simparpostfix}.json'

    anaparams = conf_reader.parse_file(anaparams_file)
    simparams = conf_reader.parse_file(simparams_file)

    print(anaparams['outprefix'])
    if not anaparams["outprefix"].startswith("/"):
        dirprefix = os.getcwd()
        print("CWD",dirprefix)
        # anaparams["outprefix"] = dirprefix.join(anaparams["outprefix"])
        anaparams["outprefix"] = dirprefix + "/" + anaparams["outprefix"]

    print(anaparams['outprefix'])
    
    # exit(0)
    all_paths = file_handling.make_allpaths(anaparams["outprefix"],
                                                  simparams["confdir"],
                                                  combined_str)
    

    print(all_paths)

    print(len(all_paths),len(combined_str))

    print(all_paths[0])
    # decomp_series = pkl.load(open(all_paths[0][1] + "/HiC_dir/pickleDecomp", 'rb'))
    # HiC_series = pkl.load(open(all_paths[0][1] + "/HiC_dir/pickleAllHiC", 'rb'))

    # print(decomp_series.shape,HiC_series.shape)

    #### this has to be generalised by input
    ## for rejuv
    if args.mode == "rejuv":
        const_dict = {'sga': '0.02', 'as': '0.6', 'sa': '100000', 'epa' : '6.00'}

    ## for dediff
    if args.mode == "dediff":
        # const_dict = {'sgat': '2','lmn': '0.1', 'epa': '4.0', 'sga': '0.001'}
        const_dict = {'sgat': '2','lmn': '0.1', 'sga': '0.001'}
    indices,combined_list = IMLib_utils.organise_by(combined_str,const_dict )
    print(indices)
    print(combined_list)
    print([combined_str[i] for i in indices])

    combined_str_len = len(indices)

    # combined_str_len = 1 ### do onyl one at a time for memory reasons
    # setlist_len = 4

    domain_data = DomainData(combined_str_len,args,params,setlist_len= 1)


    domain_data.add_multiple(indices,all_paths)
    
    start=800000
    step=10000
    stop=800001

    intMat, refconfig_init = file_handling.get_conf(
        'sample_conf/', "as0.6", 1, simparams, start, step, stop)

    refconfig = refconfig_init
    nsurf = simparams["nsurf"]
    nchrom = simparams["nchrom"]
    npatch = simparams["npatch"]
    infileprefix = simparams["infileprefix"]
    nato = nsurf + nchrom + nchrom * 2 * npatch  # 1600 + 8 + 8x2x4 OR 2000 + 46 + 46x2x23

    print(domain_data.aggregate_decomp_series_sets(0).shape, domain_data.aggregate_HiC_series_sets(0).shape)

    # # Aggregate the decomp and HiC series sets
    # decomp_series_sets = domain_data.aggregate_decomp_series_sets(0)
    # HiC_series_sets = domain_data.aggregate_HiC_series_sets(0)
    IM_tavg_range = (0,-1)
    decomp_tavg_range = (0,-1)
    combined_IM_df,combined_decomp_df = IM_df_combine.gen_combined_df(domain_data, indices,combined_list,
                                                    combined_str, refconfig_init,IM_tavg_range,
                                                    decomp_tavg_range)
    # gen_combined_df(domain_data,indices,combined_list,combined_str,refconfig_init)

    aggregate_df = IM_df_combine.aggregate_data(combined_IM_df,combined_list)

    # filename = f"dediff_aggregate_data_sga{const_dict['sga']}_as{const_dict['as']}_sa{const_dict['sa']}_epa{const_dict['epa']}_set{args.set}.csv"
    print(const_dict)
    if args.mode == "rejuv":
        filename = f"rejuv_aggregate_IM_data_sga{const_dict['sga']}_as{const_dict['as']}_sa{const_dict['sa']}_epa{const_dict['epa']}_set{args.set}.csv"
    if args.mode == "dediff":
        # filename = f"flux_dediff_aggregate_IM_data_sga{const_dict['sga']}_sgat{const_dict['sgat']}_set{args.set}.csv"
        filename = f"flux_dediff_aggregate_IM_data.csv"
    filename = os.path.join(anaparams["outprefix"],"HiC_dir",filename)
    print(filename)
    aggregate_df.to_csv(filename, index=False)
    print(f"Data saved to {filename}")

    if args.mode == "rejuv":
        filename = f"rejuv_decomp_data_sga{const_dict['sga']}_as{const_dict['as']}_sa{const_dict['sa']}_epa{const_dict['epa']}_set{args.set}.csv"
    if args.mode == "dediff":
        filename = f"flux_dediff_decomp_data.csv"
    
    os.makedirs(os.path.join(anaparams["outprefix"],"HiC_dir"), exist_ok=True)
    filename = os.path.join(anaparams["outprefix"],"HiC_dir",filename)
    print(filename)

    combined_decomp_df.to_csv(filename, index=False)


    


