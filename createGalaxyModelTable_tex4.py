#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: createGalaxyModelTable_tex4.py, v 4.0 10/26/2018

Create table 2 for the rotation paper -> list of galaxies-target pairs with impact 
parameter, azimuth, velocity, absorption velocity, model velocity ranges

Based on createTargetTable_tex2.py


v4: updated for (hopefully) final version of paper (10/26/18)
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
#         resultsName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/salt_galaxy_sightlines_cut_plus_ancillary_fits.csv'
#         outName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/model_table3.txt'

        resultsName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/salt_galaxy_sightlines_cut_plus_ancillary_fits_newerrs.csv'
        outName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/model_table_redo2.txt'
        
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
    'System #',\
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
    sysNumberList = ['System #']
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
        e_Vhel_measured = l['e_Vhel_measured']
        Vrot_corrected = l['Vrot_corrected']
#         model_range = eval(l['model_range'])
#         NFW_range = eval(l['NFW_range'])
        e_Vrot_corrected = l['e_Vrot_corrected']

        model_range = eval(l['model_range_err'])
        NFW_range = eval(l['NFW_range_err'])
        azimuth = int(round(float(l['azimuth']),0))
        sysNumber = l['number']
        
        vrot = '${0} \pm {1}$'.format(Vrot_corrected, e_Vrot_corrected)
#         vhel = '{0} \pm {1}'.format(Vhel_measured, e_Vhel_measured)

        sysNumberList.append(sysNumber)
        nameList.append(name)
        targetList.append(target)
        rhoList.append(impact)
        azList.append(azimuth)
        v_sysList.append(Vhel_measured)
        v_rotList.append(vrot)
        v_LyaList.append(Lya_v)
        w_LyaList.append(Lya_w)
        cyl_modelList.append(model_range)
        nfw_modelList.append(NFW_range)
    
        
        model_low = int(Vhel_measured) + int(model_range[0])
        model_high = int(Vhel_measured) + int(model_range[1])
        
        NFW_low = int(Vhel_measured) + int(NFW_range[0])
        NFW_high = int(Vhel_measured) + int(NFW_range[1])
  
        summaryList.append((sysNumber,\
        name,\
        target,\
        str(impact),\
        str(azimuth),\
        Vhel_measured,\
        vrot,\
        Lya_v,\
        Lya_w,\
        str(model_low) + ' - ' + str(model_high),\
        str(NFW_low) + ' - ' + str(NFW_high)))
                    

    padding = 4
    widths =[\
    max(len(str(d)) for d in sysNumberList) + padding,\
    max(len(str(d)) for d in nameList) + padding,\
    max(len(str(d)) for d in targetList) + padding,\
    max(len(str(d)) for d in rhoList) + padding,\
    max(len(str(d)) for d in azList) + padding,\
    max(len(str(d)) for d in v_sysList) + padding,\
    max(len(str(d)) for d in v_rotList) + padding,\
    max(len(str(d)) for d in v_LyaList) + padding,\
    max(len(str(d)) for d in w_LyaList) + padding,\
    max(len(str(d)) for d in cyl_modelList) + padding,\
    max(len(str(d)) for d in nfw_modelList) + padding]

    for row in summaryList:
        output.write("".join(str(i + '  &').ljust(width) for i,width in zip(row,widths))+'\\\\\n')

   
#     length_d = max(len(str(d)) for d in nameList)
#     length_t = max(len(str(d)) for d in targetList)
#     length_r = max(len(str(d)) for d in rhoList)
#     length_a = max(len(str(d)) for d in azList)
#     length_vs = max(len(str(d)) for d in v_sysList)
#     length_vr = max(len(str(d)) for d in v_rotList)
#     length_vlya = max(len(str(d)) for d in v_LyaList)
#     length_wlya = max(len(str(d)) for d in w_LyaList)
#     length_cyl = max(len(str(d)) for d in cyl_modelList)
#     length_nfw = max(len(str(d)) for d in nfw_modelList)
    
#     for d, t, r, a, vs, vr, vlya, wlya, cyl, nfw in zip(nameList, targetList, azList, v_sysList, v_rotList, v_LyaList, w_LyaList, cyl_modelList, nfw_modelList):
#         output.write("".join(str(d + ' &').ljust(width) for i, width in zip(row, widths))+'\\\\\n')
    

    file.close()
    results.close()
    output.close()

if __name__=="__main__":
    main()
