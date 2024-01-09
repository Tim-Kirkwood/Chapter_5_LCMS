# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 15:29:21 2023

@author: u03132tk
"""

import os
import pandas as pd
import pickle
import copy
import numpy as np
import networkx as nx

def weighted_j_sim(array1, array2):
    return np.minimum(array1, array2).sum()/np.maximum(array1, array2).sum()

folder = r'D:/all_mzmatch_data/negative/combined_txt'

print ('building rt db')


print ('building rt db')
peaks = [i for i in os.listdir(folder) if 'peaks' in i]
raw_peaks = {}
for file in peaks:
    name = file.split('_')
    strain, bgc = name[0].lower(), name[1] 
    with open(f'{folder}/{file}', 'rb') as in_file: 
        data = pickle.load(in_file)
    if strain not in raw_peaks.keys():
        raw_peaks[strain] = {}
    if bgc not in raw_peaks[strain].keys():
        raw_peaks[strain][bgc] = data
        
    else:
        raise ValueError('duplicate bgc')


peaks = [i for i in os.listdir(folder) if 'clustered' in i]
rt_db = {}
peak_db = {}
all_peaks_db = {}
for file in peaks:
    name = file.split('_')
    strain, bgc = name[0].lower(), name[1] 
    with open(f'{folder}/{file}', 'rb') as in_file: 
        data = pickle.load(in_file)
    rts = {peak_id : rt for peak_id, rt in zip(data['ids'], data['RTs'])}
    intensity =  {peak_id : rt for peak_id, rt in zip(data['ids'], data['Int'])}
    if strain not in rt_db.keys():
        rt_db[strain] = {}
        peak_db[strain] = {}
        all_peaks_db[strain] = {}
    if bgc not in rt_db[strain].keys():
        rt_db[strain][bgc] = rts
        peak_db[strain][bgc] = intensity
        all_peaks_db[strain][bgc] = data
    else:
        raise ValueError('duplicate bgc')

print ('building annotations db')
annotations = [i for i in os.listdir(folder) if 'annotations' in i]
db = {}
for file in annotations:
    name = file.split('_')
    strain, bgc = name[0].lower(), name[1] 
    with open(f'{folder}/{file}', 'rb') as in_file: 
        data = pickle.load(in_file)
    if strain not in db.keys():
        db[strain] = {}
    if bgc not in db[strain].keys():
        db[strain][bgc] = data
    else:
        raise ValueError('duplicate bgc')

#known annotations
print ("KNOWN ANNOTATION GRAPH")
known_edge_table = []     
compared = []  

all_predictions = []
log = ''
for strain, bgcs in db.items():
    for bgc, annotations in bgcs.items():
        print (strain, '  ', bgc)
        for other_strain, other_bgcs in db.items():
            for other_bgc, other_annotations in other_bgcs.items():
                #print (other_strain, '  ', other_bgc)
                
                #do not compare self
                if strain == other_strain and bgc == other_bgc:
                    #print ('skipping')
                    continue
                
                #if the target strain and the target BGC has already been compared, do not repeat
                if sorted([f'{strain}{bgc}', f'{other_strain}{other_bgc}']) in compared:
                    print (f'already compared {strain}::{bgc} and {other_strain}::{other_bgc}')
                    #raise ValueError
                    continue
                else:
                    print (f'comparing {strain}::{bgc} and {other_strain}::{other_bgc}')
                    compared += [sorted([f'{strain}{bgc}', f'{other_strain}{other_bgc}'])]
                #count = 0
                #reported = []
                #num_comparisons = len(annotations) * len(other_annotations)
                recorded =[] 
                for index, (annotation, predictions) in enumerate(annotations.items()):
                    #track proportion of annotations that have been compared against every other annotation
                    progress = int(100*index/len(annotations))
                    if progress%20 == 0 and progress not in recorded:
                        print (progress)
                        recorded+= [progress]
                        
                    #try removing the gibbs rule and increasing the requirment for consistency
                    #but as it stands you can find any that are in al thrree samples, confidently 
                    #identified, and are eluting at similar times
                    
                    #do not compare weak annotations
                    #if predictions['post Gibbs'][0] < 0.5:
                    #    continue
                    
                    #do not compare unknown annotations
                    test_prediction = predictions['id'][0]
                    if test_prediction == 'Unknown':
                        continue
                    
                    #keep a list of annotations so you can include singlets in graph
                    all_predictions += [f'{strain}_{bgc}_{annotation}']#for if it is unknown
                    
                    #compare the annotation against the annotations for every other strain/BGC
        
                            
                    #for the annotations in the target BGC
                    for other_annotation, other_predictions in other_annotations.items():
                        #dont compare weak or unknwon BGCs
                        #if other_predictions['post Gibbs'][0] < 0.5:
                        #    continue
                        other_prediction = other_predictions['id'][0]
                        if other_prediction == 'Unknown':
                            continue
                                
                        #draw edge if they are the same 
                        if test_prediction == other_prediction:
                            test_rt = rt_db[strain][bgc][annotation]
                            other_rt = rt_db[other_strain][other_bgc][other_annotation]
                            if abs(test_rt - other_rt) < 30:
                                known_edge_table += [[f'{strain}_{bgc}_{annotation}', f'{other_strain}_{other_bgc}_{other_annotation}']]
df = pd.DataFrame(known_edge_table)
G = nx.from_pandas_edgelist(df, source = 0, target = 1).to_undirected()
G.add_nodes_from(all_predictions)#note this will include later nodes that are in graph already but arent compared to their earlier partner, but if they are already in graph they wont be added again                        

specific_subgraphs = {'10' : [],
                      '15' : [],
                      '21' : [],
                      '25' : [],
                      '47' : []}

for subgraph in nx.connected_components(G):
    if len(subgraph)>1:
        bgcs = [i.split('_')[1] for i in subgraph]
        strains = [i.split('_')[0] for i in subgraph]
        if len(set(bgcs)) == 1:
            
            data = {}
            for i in subgraph:
                strain, bgc, annotation = i.split('_')
                peak_df = all_peaks_db[strain][bgc]
                peak = peak_df[peak_df['ids'] == int(annotation)]
                
                raw_peak_df = raw_peaks[strain][bgc]
                raw_peak = raw_peak_df[raw_peak_df['id'] == int(annotation)] 
                cols = [i for i in raw_peak.columns if i not in {'id', 'mass', 'RT'}]
                if (raw_peak[cols].values != 0).sum()>=2:
                    data[i] = {'annotations' : db[strain][bgc][int(annotation)], 
                               'intensity' : peak_db[strain][bgc][int(annotation)],
                               'peak' : peak,
                               'raw peak' : raw_peak}
                
            if len(data) > 0:
                specific_subgraphs[bgc] += [data]
