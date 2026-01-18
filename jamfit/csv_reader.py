import numpy as np
#import os
import re
from collections import OrderedDict
import csv
#import matplotlib.pyplot as plt
#import pandas as pd
#import sys
#if sys.version_info[0] < 3: 
#    from StringIO import StringIO
#else:
#    from io import StringIO


def csv_read(file0, silent=True, rmv_lns=False, rmv_col=False, hdr_keys=False, strp_units=False, delim=False):
    #Code to read a csv file and turn it into a dictionary
    
    csv_f = open(file0, mode='r', encoding='utf-8-sig')

    if delim != False:
        data_struct0  = csv.reader(csv_f, delimiter=delim)
    else:
        data_struct0  = csv.reader(csv_f)

    csv_datalst = []
    for row in data_struct0:
        row[:] = ['nan' if x=='' else x for x in row] #Replace all '' entries with nan
        csv_datalst.append(row)
        
    if rmv_lns != False:  #If specified, remove lines from the beginning
        csv_datalst = csv_datalst[rmv_lns:len(csv_datalst)]            
            
    if rmv_col != False:
	    for i in range(len(csv_datalst)):
                csv_datalst[i] = csv_datalst[i][rmv_col:]


    csv_data_arr = np.array(csv_datalst)  #Turn the list of rows into an array for indexing 
    csv_dict = OrderedDict() #Initialize the ordered dictionary 

    #If no header keys are specified, organize the data into a dictionary with columns: "col1, col2, etc..."
    if hdr_keys == False:
        for i in range(np.shape(csv_data_arr)[1]):
            try:
                csv_dict['col'+str(i+1)] = np.asarray(csv_data_arr[:,i], dtype='float')
            except:
                #csv_data_arr[:,i] = [j.strip() for j in csv_data_arr[:,i]]
                csv_data_arr[:,i] = list(map(str.strip, csv_data_arr[:,i])) #strip the blank spaces from the strings
                csv_dict['col'+str(i+1)] = csv_data_arr[:,i]
                
            
    #If hdr_keys == True, use the csv header as the dictionary keys     
    else:
        if strp_units: #If specified remove units in brackets or parathensis 
            for i in range(len(csv_data_arr[0])):
                csv_data_arr[0][i] = re.sub(r'\([^()]*\)', '', csv_data_arr[0][i])
                csv_data_arr[0][i] = re.sub(r'\[[^()]*\]', '', csv_data_arr[0][i])

        key_arr = list(map(str.strip, csv_data_arr[0])) #strip the blank spaces from the strings  
       
        for i in range(np.shape(csv_data_arr)[1]):
            try:
                csv_dict[key_arr[i]] = np.asarray(csv_data_arr[1:,i], dtype='float')
            except:
                csv_data_arr[:,i] = list(map(str.strip, csv_data_arr[:,i])) #strip the blank spaces from the strings
                csv_dict[key_arr[i]] = csv_data_arr[1:,i]


    if silent == False:
        for key in csv_dict:
            print(key, csv_dict[key])
                
    csv_f.close()
    
    return csv_dict


