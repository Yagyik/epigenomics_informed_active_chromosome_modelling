import pandas as pd
import numpy as np
import itertools
from .domainProps import compute_expected_return_time, compute_average_intermingling, compute_decomp

def organise_by(combined_str, key_value_pairs):
    # Find indices where key value pairs are present
    indices = [i for i, d in enumerate(combined_str) if all(d.get(k) == v for k, v in key_value_pairs.items())]
    print(indices)
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

# # indices,combined_list = organise_by(combined_str, {'sga': '0.02', 'as': '0.6', 'sa': '100000', 'tl' : '5000000'})
# indices,combined_list = organise_by(combined_str, {'sga': '0.02', 'as': '0.6', 'sa': '100000', 'epa' : '6.00'})
# print(indices)
# print(combined_list)
# print([combined_str[i] for i in indices])


def gen_combined_df(domain_data,indices,combined_list,combined_str,refconfig_init,
                    IM_tavg_range=None,decomp_tavg_range=None):

    combined_IM_df = pd.DataFrame()
    combined_decomp_df = pd.DataFrame()
    for i in range(len(indices)):
        
        # Aggregate the decomp and HiC series sets
        decomp_series_sets = domain_data.aggregate_decomp_series_sets(i)
        HiC_series_sets = domain_data.aggregate_HiC_series_sets(i)
        if IM_tavg_range == None:
            IM_tavg_range = (HiC_series_sets.shape[1] - 5, HiC_series_sets.shape[1])

        print(IM_tavg_range)

        if decomp_tavg_range == None:
            decomp_tavg_range = (decomp_series_sets.shape[1] - 5, decomp_series_sets.shape[1])
        # tavg_range = (HiC_series_sets.shape[1] // 2, HiC_series_sets.shape[1] // 2 + 5)\
        # tavg_range = (50,55)
        print(decomp_tavg_range)

        vk = combined_list[i]
        print(vk,indices[i],combined_str[indices[i]])

        df_decomp = compute_decomp(decomp_series_sets, refconfig_init, tavg_range=decomp_tavg_range)
        df_decomp = df_decomp.rename(columns=lambda x: f'{vk}_{x}' if x not in ["Chrom1-id", "Domain1","Domain1-points-to",
                                                                                        "Specific"] else x)
        
        
        # Calculate return times
        df_return_homologs = compute_expected_return_time(HiC_series_sets, refconfig_init,
                                                           tavg_range=IM_tavg_range)
        df_return_homologs = df_return_homologs.rename(columns=lambda x: f'{vk}_{x}' if x not in ["Chrom1-id", "Chrom2-id", "Domain1", "Domain1-points-to",
                                                                                        "Domain2","Domain2-points-to",
                                                                                        "InteractionStrength"] else x)
        # # Calculate average intermingling
        df_im_homologs = compute_average_intermingling(HiC_series_sets, refconfig_init,
                                                        threshold=0.5, tavg_range=IM_tavg_range,mode='sum',plotflag=False)
        df_im_homologs = df_im_homologs.rename(columns=lambda x: f'{vk}_{x}' if x not in ["Chrom1-id", "Chrom2-id", "Domain1", "Domain1-points-to",
                                                                                        "Domain2","Domain2-points-to",
                                                                                        "InteractionStrength"] else x)  
        
        # Combine the dataframes

        if i== 0:
            df_return_homologs_data = df_return_homologs.copy()
            df_im_homologs_data = df_im_homologs.drop(columns=["Chrom1-id", "Chrom2-id", "Domain1",
                                                        "Domain1-points-to", "Domain2",
                                                        "Domain2-points-to"])
            df_decomp_data = df_decomp.copy()

        if i > 0:
            df_return_homologs_data = df_return_homologs.drop(columns=["Chrom1-id", "Chrom2-id", "Domain1",
                                                        "Domain1-points-to", "Domain2",
                                                        "Domain2-points-to"])
            
            df_im_homologs_data = df_im_homologs.drop(columns=["Chrom1-id", "Chrom2-id", "Domain1",
                                                        "Domain1-points-to", "Domain2",
                                                        "Domain2-points-to","InteractionStrength"])
            df_decomp_data = df_decomp.drop(columns=["Chrom1-id", "Domain1","Domain1-points-to",
                                                        "Specific"])
        print("returntime data cols",df_return_homologs_data.columns)    
        
        print("int ming data cols",df_im_homologs_data.columns)
        combined_IM_df = pd.concat([combined_IM_df, df_return_homologs_data], axis=1)
        print("combined return time",combined_IM_df.columns)
        combined_IM_df = pd.concat([combined_IM_df, df_im_homologs_data], axis=1)

        combined_decomp_df = pd.concat([combined_decomp_df, df_decomp_data], axis=1)
        print("combined decomp",combined_decomp_df.columns)

        # combined_df = pd.concat([df_return_homologs, df_im_homologs.drop(columns=["Chrom1-id", "Chrom2-id", "Domain1", "Domain1-points-to", "Domain2", "Domain2-points-to"])], axis=1)

        del HiC_series_sets, decomp_series_sets, df_return_homologs, df_im_homologs,
        df_return_homologs_data, df_im_homologs_data

    return combined_IM_df,combined_decomp_df

# print(combined_df)

def aggregate_data(combined_df,combined_list):
    # Group by Chrom1 and Chrom2


    unique_chrom_pairs = combined_df[['Chrom1-id', 'Chrom2-id']].drop_duplicates()

    print("how any chrom pairs??",len(unique_chrom_pairs.index))
    aggregated_df = pd.DataFrame()
    counter = 0
    for _, row in unique_chrom_pairs.iterrows():
        chrom1_id = row['Chrom1-id']
        chrom2_id = row['Chrom2-id']
        # print(chrom1_id,chrom2_id,combined_df['Chrom1-id'] == chrom1_id,combined_df['Chrom2-id'] == chrom2_id)
        subset = combined_df[(combined_df['Chrom1-id'] == chrom1_id) & (combined_df['Chrom2-id'] == chrom2_id)]
        
        interaction_strengths = subset['InteractionStrength'].tolist()
        max_interminglings = subset[[col for col in subset.columns if 'MaxIntermingling' in col]].values.tolist()

        # Check if "Domain1-points-to" = chrom1_id+1 and "Domain2-points-to" = chrom2_id + 1
        # print(np.unique(subset['Domain1-points-to'].values),np.unique(subset['Domain2-points-to'].values),chrom1_id,chrom2_id)
        valid_rows = subset[(subset['Domain1-points-to'] == chrom2_id + 1) & (subset['Domain2-points-to'] == chrom1_id + 1)]
        print("valid rows",len(valid_rows.index))
        if not valid_rows.empty:
            max_interaction_strength = valid_rows['InteractionStrength'].max()
            specific = 1
                
        else:
            max_interaction_strength = subset['InteractionStrength'].mean() ## mean of all interaction strengths for these two chroms
            specific = 0
        row_dict = {
            'Chrom1-id': int(chrom1_id),
            'Chrom2-id': int(chrom2_id),
            'InteractionStrength': max_interaction_strength,
            'Specific': specific
        } 
        print("row_dict",row_dict)

        
        # if not valid_rows.empty:      
        #     target_intermingling = valid_rows[[col for col in valid_rows.columns if 'MaxIntermingling' in col and 'discretized' not in col]].values[0]
        # else:
        #     target_intermingling = subset[[col for col in subset.columns if 'MaxIntermingling' in col and 'discretized' not in col]].mean().values[0]

        # mean_intermingling = subset[[col for col in subset.columns if 'MaxIntermingling' in col and 'discretized' not in col]].mean().values[0]

        # min_return_time = subset[[col for col in subset.columns if 'ReturnTime' in col]].min().values[0]
        # max_max_intermingling = subset[[col for col in subset.columns if 'MaxIntermingling' in col and 'discretized' not in col]].max().values[0]

        # row_dict.update({
        #     'MeanIntermingling': mean_intermingling,
        #     'minReturnTime': min_return_time,
        #     'MaxIntermingling': max_max_intermingling,
        #     'targetIntermingling': target_intermingling
        # })

        # aggregated_df = pd.concat([aggregated_df, pd.DataFrame(row_dict, index=[counter])], axis=0)

        
        
        
        for va in combined_list:
            if not valid_rows.empty:
                target_intermingling = valid_rows[[col for col in valid_rows.columns 
                                                    if col.startswith(va) and 
                                                    'MaxIntermingling' in col and 
                                                    'discretized' not in col]].values[0]
                # print(row_dict,target_intermingling,target_intermingling[0])
                # target_intermingling = target_intermingling[0]

            else:
                target_intermingling = subset[[col for col in subset.columns 
                                                if col.startswith(va) and 
                                                'MaxIntermingling' in col and 
                                                'discretized' not in col]].mean().values[0]
                # target_intermingling = np.nan
                
                print("else",chrom1_id,chrom2_id,target_intermingling)

            # print("subset",subset[[col for col in subset.columns if col.startswith(va) and 'MaxIntermingling' in col]])
            mean_intermingling = subset[[col for col in subset.columns if col.startswith(va) and
                                            'MaxIntermingling' in col and 
                                            'discretized' not in col]].mean().values[0]
            min_return_time = subset[[col for col in subset.columns if col.startswith(va) and 
                                        'ReturnTime' in col]].min().values[0]
            max_max_intermingling = subset[[col for col in subset.columns if col.startswith(va) and
                                                'MaxIntermingling' in col and
                                                'discretized' not in col]].max().values[0]
            print(mean_intermingling,min_return_time,max_max_intermingling,target_intermingling)
            # print(mean_intermingling,min_return_time,max_max_intermingling,target_intermingling)
            row_dict.update({
                f'{va}_MeanIntermingling': mean_intermingling,
                f'{va}_minReturnTime': min_return_time,
                f'{va}_MaxIntermingling': max_max_intermingling,
                f'{va}_targetIntermingling': target_intermingling
            })
        print(row_dict)        
        aggregated_df = pd.concat([aggregated_df, pd.DataFrame(row_dict, index=[counter])], axis=0)
        
        print(f"Chrom1-id: {chrom1_id}, Chrom2-id: {chrom2_id}")
        # print(f"Domain1-points-to: {subset['Domain1-points-to'].tolist()}",len(subset['Domain1-points-to'].tolist()))
        # print(f"Domain2-points-to: {subset['Domain2-points-to'].tolist()}",len(subset['Domain2-points-to'].tolist()))
        # print(f"InteractionStrengths: {[f'{x:.2f}' for x in interaction_strengths]}")
        # print(f"MaxInterminglings: {[f'{x:.2f}' for sublist in max_interminglings for x in sublist]}")
        print(f"Number of rows: {len(subset)}")
        print(f"Length of InteractionStrengths: {len(interaction_strengths)}")
        print(f"Length of MaxInterminglings: {len(max_interminglings)}")
        

        counter += 1
        # if counter == 5:
        #     break
    print(aggregated_df)

    return aggregated_df

# def aggregate_decomp_data(combined_df,combined_list):
#     aggregated_df = pd.DataFrame()
#     counter = 0
    