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
from scipy.stats import norm, poisson, lognorm


def compute_expected_return_time(intermingling_data, refconfig,tavg_range=None):
    """
    Computes the expected return time for intermingling events based on sign change frequency, persistence time,
    and distance from full transition, aligned with the refconfig structure.

    Parameters:
    - intermingling_data: np.ndarray of shape (n_timepoints, 46, 46, 23, 23, 2)
    - refconfig: Object containing chromosome data and domain mappings

    Returns:
    - df_pairs: DataFrame with return times for each domain pair
    - df_homologs: DataFrame combining homologous chromosome pairs
    """
    n_sets, n_timepoints, nchrom, _, npatch, _, _ = intermingling_data.shape
    if tavg_range is None:
        start_time = 0
        end_time = -1
    else:
        start_time = tavg_range[0]
        end_time = tavg_range[1]
    return_list = []
    homolog_list = []
    domain_list = []
    pairslist = []
    print(n_timepoints, nchrom, npatch)

    for i in range(int(0.5 * nchrom)):
        itype = (i + 1) % 24
        i1 = i
        i2 = 23 + itype - 1
        
        patchlisti1 = refconfig.chromarr[i1].patchlist
        patchpointeri1 = refconfig.chromarr[i1].patchpointer
        patchlisti2 = refconfig.chromarr[i2].patchlist
        patchpointeri2 = refconfig.chromarr[i2].patchpointer
        
        jlist = [j for j in range(i+1,23) if j + 1 != itype]
        # jlist.append(itype - 1)
        # print(i,jlist)
        for jj in range(1, npatch-1):
            j1 = patchpointeri1[jj]
            jg1 = patchlisti1[jj]
            j2 = patchpointeri2[jj]
            jg2 = patchlisti2[jj]

            # print("checking matches", jj, j1, j2, jg1, jg2)
            
            for k in jlist:
                k1 = k
                k2 = 23 + k

                if i1 == k1 or i1 == k2 or i2 == k1 or i2 == k2:
                    continue
                
                patchlistk1 = refconfig.chromarr[k1].patchlist
                patchpointerk1 = refconfig.chromarr[k1].patchpointer
                patchlistk2 = refconfig.chromarr[k2].patchlist
                patchpointerk2 = refconfig.chromarr[k2].patchpointer
                
                for ll in range(1, npatch-1):
                    l1 = patchpointerk1[ll]
                    lg1 = patchlistk1[ll]
                    l2 = patchpointerk2[ll]
                    lg2 = patchlistk2[ll]
                    
                    pairs = [
                        (i1, jj,j1, jg1, k1, ll,l1, lg1),
                        (i1, jj,j1, jg1, k2, ll,l2, lg2),
                        (i2, jj,j2, jg2, k1, ll,l1, lg1),
                        (i2, jj,j2, jg2, k2, ll,l2, lg2)
                    ]
                    # print(pairs)
                    pairslist.append(pairs)

    print(len(pairslist))                
    for pairs in pairslist:
        # print(pairs)
        pair_return_times = []                
        for (chrom1, domain1,d1p2, dg1, chrom2, domain2,d2p2, dg2) in pairs:

            for s in range(n_sets):
                trajectory = intermingling_data[s,start_time:end_time, chrom1, chrom2, domain1, domain2, 1].reshape(-1)


                # print(intermingling_data.shape,trajectory.shape)
                delta_traj = np.diff(trajectory)
                change_points = np.where(np.sign(delta_traj[:-1]) != np.sign(delta_traj[1:]))[0]
                persistence_lengths = np.diff(change_points) if len(change_points) > 1 else []
                distances_from_1 = 1 - trajectory[change_points] if len(change_points) > 0 else []
                
                if len(persistence_lengths) > 1 and len(distances_from_1) > 1:
                    mean_persistence_time = np.mean(persistence_lengths)
                    
                    # if np.std(persistence_lengths) > np.mean(persistence_lengths):
                    #     rate_estimate = poisson.fit(persistence_lengths)[0]
                    if np.std(persistence_lengths) > np.mean(persistence_lengths):
                        rate_estimate = lognorm.fit(persistence_lengths)[0]
                    else:
                        rate_estimate = norm.fit(persistence_lengths)[0]
                    
                    weighted_rate = np.mean(rate_estimate / distances_from_1)
                    return_time = 1 / weighted_rate if weighted_rate > 0 else np.nan
                else:
                    return_time = np.nan
                
                pair_return_times.append(return_time)
            # return_list.append([chrom1, domain1, patchpointeri[domain1], dg1, chrom2, domain2, patchpointerk[domain2], dg2, return_time])
            # return_list.append([chrom2, domain2, chrom1, domain1, return_time])
        [(i1, jj,j1, jg1, k1, ll,l1, lg1),
        (i1, jj,j1, jg1, k2, ll,l2, lg2),
        (i2, jj,j2, jg2, k1, ll,l1, lg1),
        (i2, jj,j2, jg2, k2, ll,l2, lg2)] = pairs
        if j1 != j2 or l1 != l2:
            print("ERROR -- very very wrong",i1,i2,j1,j2,jj,k1,k2,l1,l2,ll)
        homolog_list.append([i1, k1, jj, j1, ll, l1, min(pair_return_times)])
        # domain_list.append([i1, jj, jg1, j1, ll,return_time])
        # domain_list.append([i2, jj, jg2, j2, ll,return_time])
        # domain_list.append([k1, ll, lg1, l1, jj,return_time])
        # domain_list.append([k2, ll, lg2, l2, jj,return_time])
    
    # df_pairs = pd.DataFrame(return_list, columns=["Chrom1", "Domain1", "Domain1-points-to", "Domain1-glob-index", "Chrom2", "Domain2", "Domain2-points-to", "Domain2-glob-index", "ReturnTime"])
    df_homologs = pd.DataFrame(homolog_list, columns=["Chrom1-id", "Chrom2-id", "Domain1", "Domain1-points-to",
                                                       "Domain2","Domain2-points-to","ReturnTime"])
    # df_domains = pd.DataFrame(domain_list, columns=["Chrom", "Domain", "Domain-glob-index", "Domain-points-to",
                                                    #  "Domain2", "ReturnTime"])
    
    return df_homologs #, df_domains


def compute_average_intermingling(intermingling_data, 
                                  refconfig, threshold,tavg_range=None,mode='avg',plotflag=False):
    """
    Computes the average intermingling for each domain pair and applies various filters based on interaction strength.

    Parameters:
    - intermingling_data: np.ndarray of shape (n_timepoints, 46, 46, 23, 23, 2)
    - refconfig: Object containing chromosome data and domain mappings
    - threshold: float - Threshold for classifying interactions

    Returns:
    - df_pairs: DataFrame with intermingling data for each domain pair
    - df_homologs: DataFrame combining homologous chromosome pairs
    - df_domains: DataFrame summarizing intermingling at the domain level
    """
    n_sets, n_timepoints, nchrom, _, npatch, _, _ = intermingling_data.shape
    if tavg_range is None:
        start_time = 0
        end_time = -1
    else:
        start_time = tavg_range[0]
        end_time = tavg_range[1]
    interaction_list = []
    domain_list = []
    homolog_list = []
    pairslist = []
    
    for i in range(int(0.5 * nchrom)):
        itype = (i + 1) % 24
        i1 = i
        i2 = 23 + itype - 1
        
        patchlisti1 = refconfig.chromarr[i1].patchlist
        patchpointeri1 = refconfig.chromarr[i1].patchpointer
        patchlisti2 = refconfig.chromarr[i2].patchlist
        patchpointeri2 = refconfig.chromarr[i2].patchpointer
        
        jlist = [j for j in range(i+1,23) if j + 1 != itype]
        # jlist.append(itype - 1)
        
        for jj in range(1, npatch-1):
            j1 = patchpointeri1[jj]
            jg1 = patchlisti1[jj]
            j2 = patchpointeri2[jj]
            jg2 = patchlisti2[jj]
            
            for k in jlist:
                k1 = k
                k2 = 23 + k
                
                patchlistk1 = refconfig.chromarr[k1].patchlist
                patchpointerk1 = refconfig.chromarr[k1].patchpointer
                patchlistk2 = refconfig.chromarr[k2].patchlist
                patchpointerk2 = refconfig.chromarr[k2].patchpointer
                
                for ll in range(1, npatch-1):
                    l1 = patchpointerk1[ll]
                    lg1 = patchlistk1[ll]
                    l2 = patchpointerk2[ll]
                    lg2 = patchlistk2[ll]
                    
                    pairs = [
                        (i1, jj,j1, jg1, k1, ll,l1, lg1),
                        (i1, jj,j1, jg1, k2, ll,l2, lg2),
                        (i2, jj,j2, jg2, k1, ll,l1, lg1),
                        (i2, jj,j2, jg2, k2, ll,l2, lg2)
                    ]
                    pairslist.append(pairs)
    print(len(pairs))

    if plotflag:
        # nrows = (len(pairslist) + 7) // 8  # Calculate the number of rows needed
        nrows = 4
        fig, axes = plt.subplots(nrows=nrows, ncols=8, figsize=(20, 2.5 * nrows))
        axes = axes.flatten()  # Flatten the axes array for easy indexing
    print("made figures")
    for idx, pairs in enumerate(pairslist):
        pair_intermingling = []
        aggregated_intermingling = []
        # print(pairs)
        for (chrom1, domain1, d1p2, dg1, chrom2, domain2, d2p2, dg2) in pairs:
            
            if mode == 'avg':
                mean_intermingling = np.mean(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 1])
                interaction_strength = np.mean(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 0])
            elif mode == 'max':
                mean_intermingling = np.max(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 1])
                interaction_strength = np.max(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 0])
            elif mode == "median":
                mean_intermingling = np.median(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 1])
                interaction_strength = np.median(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 0])

            elif mode == "sum":
                mean_intermingling = np.sum(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 1])
                interaction_strength = np.sum(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 0])

            interaction_type = "positive" if interaction_strength >= threshold else "negative"

            aggregated_intermingling.extend(intermingling_data[:, start_time:end_time, chrom1, chrom2, domain1, domain2, 1].flatten())

            # interaction_list.append([chrom1, domain1, patchpointeri[domain1], dg1, chrom2, domain2, patchpointerk[domain2], dg2, interaction_strength, mean_intermingling, interaction_type])
            pair_intermingling.append(mean_intermingling)
        # print(len(aggregated_intermingling))
        # max across pairs of homolog domains for this chromdom pair
        if len(pair_intermingling) > 0:
            max_intermingling = max(pair_intermingling)
        else:
            max_intermingling = 0
        
        
    
        [(i1, jj, j1, jg1, k1, ll, l1, lg1),
            (i1, jj, j1, jg1, k2, ll, l2, lg2),
            (i2, jj, j2, jg2, k1, ll, l1, lg1),
            (i2, jj, j2, jg2, k2, ll, l2, lg2)] = pairs

        homolog_list.append([i1, k1, jj, d1p2, ll, d2p2, np.log(0.000001 + max_intermingling),
                             np.log(0.000001 + interaction_strength)])
        
        # domain_list.append([i1, jj, dg1, d1p2, ll, interaction_strength, max_intermingling])
        # domain_list.append([i2, jj, dg2, d2p2, jj, interaction_strength, max_intermingling])
        # domain_list.append([k1, ll, dg1, d1p2, ll, interaction_strength, max_intermingling])
        # domain_list.append([k2, ll, dg2, d2p2, jj, interaction_strength, max_intermingling])
        plotfrom = len(pairslist) - 32
        if plotflag and idx - plotfrom > 0 and idx - plotfrom < 32 :
            axes[idx-plotfrom].hist(aggregated_intermingling, bins=30)
            axes[idx-plotfrom].set_title(f'Pair {idx}- C1{i1}d{jj} - C2{k1}d{ll}')
            axes[idx-plotfrom].set_xlabel('Intermingling')
            axes[idx-plotfrom].set_ylabel('Frequency')
            axes[idx-plotfrom].set_yscale('log')

    if plotflag:
        plt.tight_layout()
        plt.show()
    # df_pairs = pd.DataFrame(interaction_list, columns=["Chrom1", "Domain1", "Domain1-points-to", "Domain1-glob-index", "Chrom2", "Domain2", "Domain2-points-to", "Domain2-glob-index", "InteractionStrength", "MeanIntermingling", "Type"])
    df_homologs = pd.DataFrame(homolog_list, columns=["Chrom1-id", "Chrom2-id", "Domain1", "Domain1-points-to",
                                                       "Domain2", "Domain2-points-to",
                                                         "MaxIntermingling", "InteractionStrength"])
    # df_domains = pd.DataFrame(domain_list, columns=["Chrom", "Domain", "Domain-glob-index", "Domain-points-to", "Domain2",
                                                    #  "InteractionStrength", "MaxIntermingling"])
    
    return df_homologs #, df_domains


    

def construct_baseline_cases(intermingling_data, refconfig, method="random", n_permutations=100, threshold=0.7):
    """
    Constructs baseline cases by applying different randomization strategies to intermingling data and computes
    average intermingling for each domain pair.

    Parameters:
    - intermingling_data: np.ndarray of shape (n_timepoints, 46, 46, 23, 23, 2)
    - refconfig: Object containing chromosome data and domain mappings
    - method: str - Baseline generation method ("random", "bootstrapped", "no-attraction")
    - n_permutations: int - Number of permutations to generate baseline data
    - threshold: float - Threshold for classifying interactions

    Returns:
    - df_pairs: DataFrame with intermingling data for each domain pair
    - df_homologs: DataFrame combining homologous chromosome pairs
    - df_domains: DataFrame summarizing intermingling at the domain level
    """
    baseline_data = np.copy(intermingling_data)
    n_sets, n_timepoints, nchrom, _, npatch, _, _ = intermingling_data.shape
    
    interaction_list = []
    domain_list = []
    homolog_list = []
    pairslist = []
    
    for _ in range(n_permutations):
        if method == "random":
            np.random.shuffle(baseline_data[:,:, :, :, :, :, 1])
        elif method == "bootstrapped":
            for t in range(baseline_data.shape[0]):
                np.random.shuffle(baseline_data[:,t, :, :, :, :, 1])  # Shuffle within each trajectory
        elif method == "no-attraction":
            baseline_data[:, :, :, :, :, 1] = 0  # Remove intermingling effects

        for i in range(int(0.5 * nchrom)):
            itype = (i + 1) % 24
            i1 = i
            i2 = 23 + itype - 1
            
            patchlisti1 = refconfig.chromarr[i1].patchlist
            patchpointeri1 = refconfig.chromarr[i1].patchpointer
            patchlisti2 = refconfig.chromarr[i2].patchlist
            patchpointeri2 = refconfig.chromarr[i2].patchpointer
            
            jlist = [j for j in range(i+1,23) if j + 1 != itype]
            # jlist.append(itype - 1)
            
            for jj in range(1, npatch -1):
                j1 = patchpointeri1[jj]
                jg1 = patchlisti1[jj]
                j2 = patchpointeri2[jj]
                jg2 = patchlisti2[jj]
                
                for k in jlist:
                    k1 = k
                    k2 = 23 + k
                    
                    patchlistk1 = refconfig.chromarr[k1].patchlist
                    patchpointerk1 = refconfig.chromarr[k1].patchpointer
                    patchlistk2 = refconfig.chromarr[k2].patchlist
                    patchpointerk2 = refconfig.chromarr[k2].patchpointer
                    
                    for ll in range(1, npatch - 1):
                        l1 = patchpointerk1[ll]
                        lg1 = patchlistk1[ll]
                        l2 = patchpointerk2[ll]
                        lg2 = patchlistk2[ll]
                        
                        pairs = [
                        (i1, jj,j1, jg1, k1, ll,l1, lg1),
                        (i1, jj,j1, jg1, k2, ll,l2, lg2),
                        (i2, jj,j2, jg2, k1, ll,l1, lg1),
                        (i2, jj,j2, jg2, k2, ll,l2, lg2)
                    ]
                        pairslist.append(pairs)
    
    for pairs in pairslist:
        pair_intermingling = []
        
        for (chrom1, domain1,d1p2, dg1, chrom2, domain2,d2p2, dg2) in pairs:
            
            mean_intermingling = np.mean(baseline_data[:,:, chrom1, chrom2, domain1, domain2, 1])
            interaction_strength = np.mean(baseline_data[:,:, chrom1, chrom2, domain1, domain2, 0])
            interaction_type = "positive" if interaction_strength >= threshold else "negative"
            
            # interaction_list.append([chrom1, domain1, patchpointeri[domain1], dg1, chrom2, domain2, patchpointerk[domain2], dg2, interaction_strength, mean_intermingling, interaction_type])
            pair_intermingling.append(mean_intermingling)
        
        max_intermingling = max(pair_intermingling)
        homolog_list.append([i1, k1, jj, d1p2,dg1, ll, d2p2,dg2, max_intermingling,interaction_strength])
        
        [(i1, jj,j1, jg1, k1, ll,l1, lg1),
        (i1, jj,j1, jg1, k2, ll,l2, lg2),
        (i2, jj,j2, jg2, k1, ll,l1, lg1),
        (i2, jj,j2, jg2, k2, ll,l2, lg2)] = pairs

        domain_list.append([i1, jj, dg1, d1p2, ll, interaction_strength, max_intermingling])
        domain_list.append([k1, ll, dg2, d2p2, jj,interaction_strength, max_intermingling])
    
    # df_pairs = pd.DataFrame(interaction_list, columns=["Chrom1", "Domain1", "Domain1-points-to", "Domain1-glob-index", "Chrom2", "Domain2", "Domain2-points-to", "Domain2-glob-index", "InteractionStrength", "MeanIntermingling", "Type"])
    df_homologs = pd.DataFrame(homolog_list, columns=["Chrom1-id", "Chrom2-id", "Domain1", "Domain1-points-to","Domain1-glob-ind", "Domain2",
                                                       "Domain2-points-to","Domain2-glob-ind", "MaxIntermingling", "InteractionStrength"])
    # df_domains = pd.DataFrame(domain_list, columns=["Chrom", "Domain", "Domain-glob-index", "Domain-points-to", "Domain2",
                                                    #  "InteractionStrength", "MaxIntermingling"])
    
    return  df_homologs #, df_domains

def compute_decomp(decomp_data,
                   refconfig, tavg_range=None,mode='avg',plotflag=False):
    
    print(decomp_data.shape)
    
    n_sets, n_timepoints,nchrom, npatch = decomp_data.shape
    if tavg_range is None:
        start_time = 0
        end_time = -1
    else:
        start_time = tavg_range[0]
        end_time = tavg_range[1]

    decomp_list = []

    
    for i in range(int(0.5 * nchrom)):
        itype = (i + 1) % 24
        i1 = i
        i2 = 23 + itype - 1
        
        patchlisti1 = refconfig.chromarr[i1].patchlist
        patchpointeri1 = refconfig.chromarr[i1].patchpointer
        # patchlisti2 = refconfig.chromarr[i2].patchlist
        # patchpointeri2 = refconfig.chromarr[i2].patchpointer
        
        # jlist = [j for j in range(i+1,23) if j + 1 != itype]
        # jlist.append(itype - 1)
        
        for jj in range(1, npatch-1):
            j1 = patchpointeri1[jj]
            jg1 = patchlisti1[jj]
            # j2 = patchpointeri2[jj]
            # jg2 = patchlisti2[jj]
            gmean = 0.5*(np.mean(decomp_data[:,start_time:end_time,i1,jj]) 
                         + np.mean(decomp_data[:,start_time:end_time,i2,jj]))
            
            

            if j1 == 1107:
                decomp_list.append([i1,jj,j1,0,gmean])
            if j1 == 1108:
                decomp_list.append([i1,jj,j1,2,gmean])
            else:
                # interaction_strength = np.mean(intermingling_data[:, start_time:end_time,
                #                                                i1, chrom2, domain1, domain2, 0])
                decomp_list.append([i1,jj,j1,1,gmean])
    
    # df_decomp = pd.DataFrame(decomp_list, columns=["Chrom1-id", "Domain1", "Domain1-points-to",
    #                                                "Domain1-glob-ind","InteractionStrength",
    #                                                "Type","Decomp"])
    df_decomp = pd.DataFrame(decomp_list, columns=["Chrom1-id", "Domain1", "Domain1-points-to",
                                                   "Specific","Decomp"])
    
    return df_decomp