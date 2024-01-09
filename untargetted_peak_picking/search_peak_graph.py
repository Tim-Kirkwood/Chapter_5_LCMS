# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 23:52:01 2023

@author: u03132tk
"""
import pandas as pd
import time
import networkx as nx
import copy
import numpy as np

#combined file with all peaks
file ='D:/all_mzmatch_data/negative/combined_combined_peakml/combined.txt'
raw_peak_table = pd.read_csv(file, sep = '\t').fillna(0).set_index('id')
min_intensity = 1000
filtered_rows = []
filtered_indicies = []
for index, row in raw_peak_table.iterrows():
    num_peaks = [i for i in row if i > min_intensity]
    if len(num_peaks) >= 3:
        filtered_rows += [row]
        filtered_indicies += [index]
        

#stop doing the filterin ghere.  you are going to cull very similar peaks that are on seperate rows
#instead, collapse the final df for a subgraph, and run the consecutive sample check on  that as a single row.         
        
peak_table = pd.DataFrame(filtered_rows, 
                          index = filtered_indicies, 
                          columns = raw_peak_table.columns)
print ('making graph...')
rt_dict = peak_table['RT'].to_dict()
mz_dict = peak_table['mass'].to_dict()
combined = {}
for peak, rt in rt_dict.items():
    combined[peak] = {'rt' : rt,
                      'mz' : mz_dict[peak]}
edge_table = []
max_rt_diff = 30#3000000000000000
max_mz_diff = 1
start = time.time()
single_nodes = []
for i, (peak, mz_rt) in enumerate(combined.items()):
    found = False
    if i%1000 == 0:
        print (i, ' out of ', len(combined), ' - seconds: ',round(time.time() - start, 2))
    for ii, (other_peak, other_mz_rt) in enumerate(combined.items()):
        if ii<=i:
            #already compared or self
            continue
        rt_diff = abs(mz_rt['rt'] - other_mz_rt['rt'])
        mz_diff = abs(mz_rt['mz'] - other_mz_rt['mz'])
        if rt_diff <= max_rt_diff and mz_diff <= max_mz_diff:
            edge_table += [[peak, other_peak]]
            found = True
    if not found:
        single_nodes += [peak]
            
#make graph from edge list
df = pd.DataFrame(edge_table)
G = nx.from_pandas_edgelist(df, source = 0, target = 1).to_undirected()
G.add_nodes_from(single_nodes)#note this will include later nodes that are in graph already but arent compared to their earlier partner, but if they are already in graph they wont be added again
print ('counting subgraphs...')
num = sum(1 for x in nx.connected_components(G))
print (num, ' subgraphs')
print ('parsing graph...')
min_intensity_cluster = 50000
bgc_map = {}
samples = ['tk24_25_r1_BA8_1_844', 
           'tk24_25_r2_BB1_1_845',
           'tk24_25_r3_BB2_1_846', 
           'm145_10_r1_BB3_1_865', 
           'm145_10_r2_BB4_1_866',
           'm145_10_r3_BB5_1_867', 
           'm145_15_r1_BC1_1_871', 
           'm145_15_r2_BC2_1_872',
           'm145_15_r3_BC3_1_873', 
           'm145_21_r1_BA5_1_859', 
           'm145_21_r2_BA6_1_860',
           'm145_21_r3_BA7_1_861', 
           'm145_25_r1_BC4_1_874', 
           'm145_25_r2_BC5_1_875',
           'm145_25_r3_BC6_1_876', 
           'm1146_25_r1_BB3_1_847',
           'm1146_25_r2_BB4_1_848', 
           'm1146_25_r3_BB5_1_849',
           'm1152_10_r1_BA2_1_856', 
           'm1152_10_r2_BA3_1_857',
           'm1152_10_r3_BA4_1_858', 
           'm1152_21_r1_BB6_1_868',
           'm1152_21_r2_BB7_1_869', 
           'm1152_21_r3_BB8_1_870',
           'tk24_10_r1_BA2_1_838', 
           'tk24_10_r2_BA3_1_839', 
           'tk24_10_r3_BA4_1_840',
           'tk24_15_r1_BA8_1_862', 
           'tk24_15_r2_BB1_1_863', 
           'tk24_15_r3_BB2_1_864',
           'tk24_21_r1_BA5_1_841', 
           'tk24_21_r2_BA6_1_842', 
           'tk24_21_r3_BA7_1_843'
           ]

inconsistent = 0
low_signal = 0
non_specific = 0

for i, subgraph_nodes in enumerate(nx.connected_components(G)):
    if i % 100 ==0:
        print (f'{i} out of {num} subgraphs')
    subgraph_samples = set()
    
    #for this you need to get peaks that are consistent (ie are present in multipl esamples - but you can check this with the first filter step)
    #next you need to find ones that are much higher in the BGC
    keep = False
    for node in subgraph_nodes:
        peak = peak_table.loc[node][samples]
        sample_names = [name for name, score in peak.items() if score > 0]
        subgraph_samples.update(sample_names)
         
        if peak.max() >= min_intensity_cluster:
            keep = True
        
    if not keep:
        low_signal += 1
        continue
    if len(subgraph_samples) < 3:
        inconsistent += 1
        continue
    bgc_set = set([i.split('_')[1] for i in subgraph_samples])
    if len(bgc_set)>1:
        non_specific += 1
        continue
    
    target_bgc = list(bgc_set)[0]
    if target_bgc not in bgc_map:
        bgc_map[target_bgc] =[]
    subgraph_df = peak_table.loc[list(subgraph_nodes)]
    mz_range = (subgraph_df['mass'].max(), subgraph_df['mass'].min())
    mean_mz = subgraph_df['mass'].mean()
    mean_rt = subgraph_df['RT'].mean()
    sample_intensities = subgraph_df[samples].replace(0, np.NaN)
    mean_intensity = sample_intensities.mean().mean()
    bgc_map[target_bgc] += [{'rt' : mean_rt, 
                             'intensity' : mean_intensity, 
                             'mass' : (mean_mz + 1, mean_mz - 1),
                             'data' : copy.deepcopy(subgraph_df)}]
        
print (f'inconsistent: {inconsistent} low_signal: {low_signal} non_specific: {non_specific}')                
    