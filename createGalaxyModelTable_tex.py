#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: createGalaxyModelTable_tex.py, v 1.0 04/28/2018

Create table 2 for the rotation paper -> list of galaxies-target pairs with impact 
parameter, azimuth, velocity, absorption velocity, model velocity ranges

Based on createTargetTable_tex2.py

"""

# from __future__ import division
import optparse
import sys
import os
# import tempfile
import csv
import string
from math import *
import numpy
import getpass
import correlateSingle7 as correlateSingle
from utilities import *

    
################################################################

def main():
    # This function reformats Bart's file
    
    if getpass.getuser() == 'frenchd':
        filename = '/Users/frenchd/Research/correlation/TARGETLIST_10_17_17_TOTAL.csv'
        resultsName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/salt_galaxy_sightlines_cut_plus_ancillary.csv'
        outName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/model_table.txt'
        
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    file = open(filename,'rU')
    reader = csv.DictReader(file)
    
    results = open(resultsName,'rU')
    resultsReader = csv.DictReader(results)
    
    output = open(outName,'wt')
    

    summaryList = []
    summaryList.append((\
    'Galaxy',\
    'Target',\
    '$\\rho$',\
    'Az.',\
    '$V_{sys}$',\
    '$V_{rot}$',\
    '$V_{Ly\\alpha}$',\
    '$W_{Ly\\alpha}$',\
    'Cyl. Model',\
    'NFW Model'))
    
    targetList = ['Target']
    nameList = ['Galaxy']
    rhoList = ['$\\rho$']
    azList = ['Az.']
    v_sysList = ['$V_{sys}$']
    v_rotList = ['$V_{rot}$']
    v_LyaList = ['$V_{Ly\\alpha}$']
    w_LyaList = ['$W_{Ly\\alpha}$']
    cyl_modelList = ['Cyl. Model']
    nfw_modelList = ['NFW Model']

                
    for l in resultsReader:
        name = l['Name']
        target = l['Target']
        Lya_v = l['Lya_v']
        impact = int(round(float(l['impact']),0))
        Lya_w = l['Lya_W']
        Vhel_measured = l['Vhel_measured']
        Vrot_corrected = l['Vrot_corrected']
        model_range = l['model_range']
        NFW_range = l['NFW_range']
        azimuth = int(round(float(l['azimuth']),0))
        
        
        nameList.append(name)
        targetList.append(target)
        rhoList.append(impact)
        azList.append(azimuth)
        v_sysList.append(Vhel_measured)
        v_rotList.append(Vrot_corrected)
        v_LyaList.append(Lya_v)
        w_LyaList.append(Lya_w)
        cyl_modelList.append(model_range)
        nfw_modelList.append(NFW_range)
    
        
        summaryList.append((name,\
        target,\
        str(impact),\
        str(azimuth),\
        Vhel_measured,\
        Vrot_corrected,\
        Lya_v,\
        Lya_w,\
        model_range,\
        NFW_range))
                    

#     padding = 4
#     widths =[\
#     max(len(str(d)) for d in nameList) + padding,\
#     max(len(str(d)) for d in targetList) + padding,\
#     max(len(str(d)) for d in rhoList) + padding,\
#     max(len(str(d)) for d in azList) + padding,\
#     max(len(str(d)) for d in v_sysList) + padding,\
#     max(len(str(d)) for d in v_rotList) + padding,\
#     max(len(str(d)) for d in v_LyaList) + padding,\
#     max(len(str(d)) for d in w_LyaList) + padding,\
#     max(len(str(d)) for d in cyl_modelList) + padding,\
#     max(len(str(d)) for d in nfw_modelList) + padding]
# 
#     for row in summaryList:
#         output.write("".join(str(i + '  &').ljust(width) for i,width in zip(row,widths))+'\\\\\n')

   
    length_d = max(len(str(d)) for d in nameList)
    length_t = max(len(str(d)) for d in targetList)
    length_r = max(len(str(d)) for d in rhoList)
    length_a = max(len(str(d)) for d in azList)
    length_vs = max(len(str(d)) for d in v_sysList)
    length_vr = max(len(str(d)) for d in v_rotList)
    length_vlya = max(len(str(d)) for d in v_LyaList)
    length_wlya = max(len(str(d)) for d in w_LyaList)
    length_cyl = max(len(str(d)) for d in cyl_modelList)
    length_nfw = max(len(str(d)) for d in nfw_modelList)
    
    for d, t, r, a, vs, vr, vlya, wlya, cyl, nfw in zip(nameList, targetList, azList, v_sysList, v_rotList, v_LyaList, w_LyaList, cyl_modelList, nfw_modelList):
        output.write("".join(str(d + ' &').ljust(
    




    file.close()
    results.close()
    output.close()

if __name__=="__main__":
    main()
