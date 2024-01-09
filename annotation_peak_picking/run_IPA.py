# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 15:01:35 2023

@author: u03132tk
"""

import pandas as pd
import numpy as np
from ipaPy2 import PeakMLIO, ipa
import shutil
import gzip
import os
import xml.etree.ElementTree as ET
import pandas
import pickle

def write_pickle(obj, path):
    with open(path, 'wb') as handle:
        pickle.dump(obj, handle)

directory = 'D:/all_mzmatch_data/negative/combined_txt'

ms1_files = [i for i in os.listdir(directory) if i[i.rindex('.'):] == '.txt']
for ms1_data in ms1_files:
    print (ms1_data)
    print ('reading in files...')
    output_file_base = f'{directory}/{ms1_data[0:ms1_data.rindex(".")]}'
    adducts = pd.read_csv('DB/adducts.csv')
    DB = pd.read_csv('DB/IPA_MS1.csv')
    input_data = pd.read_csv(f'{directory}/{ms1_data}', sep = '\t').fillna(0)
    print ('writing peak table')
    write_pickle(input_data, output_file_base + '_peaks.pkl')
    df = ipa.clusterFeatures(input_data)
    print ('writing clustered features')
    write_pickle(df, output_file_base + '_clustered.pkl')
    ipa.map_isotope_patterns(df,ionisation=-1)
    allAddsPos = ipa.compute_all_adducts(adducts, DB, ionisation=-1, ncores=1)
    annotations=ipa.MS1annotation(df,
                                  allAddsPos,
                                  ppm=50,#3,
                                  ncores=1)    
    #note these are different - not sequential - options
    zs = ipa.Gibbs_sampler_add(df,
                               annotations,
                               #noits=1000, try default value
                               #delta_add=0.1, try default value 
                               all_out=True)
    print ('writing annotations...')
    write_pickle(annotations, output_file_base + '_annotations.pkl')
    print ()
