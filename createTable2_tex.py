#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: createTargetTable_tex.py, v 1.0 02/22/2018


Based on: createTargetTable_tex.py, v 1.0 07/20/2016

Create table 1 for the pilot paper -> list of targets with ra and dec, z, program ID, 
grating, obs ID, obs date, texp and s/n.

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

from astropy.io import ascii
from astropy import table
    
################################################################

def main():
    # This function reformats Bart's file
    
    if getpass.getuser() == 'frenchd':
        filename = '/Users/frenchd/Research/correlation/TARGETLIST_10_17_17_TOTAL.csv'
#         resultsName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/salt_galaxy_sightlines_cut.csv'
        resultsName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/salt_galaxy_sightlines_cut_plus_ancillary_fits.csv'
        outName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/table2_2.txt'
        outName2 = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/table2_2_alt.txt'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    file = open(filename,'rU')
    reader = csv.DictReader(file)
    
    results = open(resultsName,'rU')
    resultsReader = csv.DictReader(results)
    
    output = open(outName,'wt')
    
    # run through the results file and collect the names of included targets
    targets = {}
#     for l in resultsReader:
#         include = eval(l['include'])
#         
#         if include:
#             targetName = l['AGNname']
#             if not targets.has_key(targetName):
#                 targets[targetName]=1

    for l in resultsReader:
        targetName = l['Target']
        galaxyName = l['Name']
        if not targets.has_key(targetName):
            targets[targetName]=[galaxyName]
        else:
            l = targets[targetName]
            if not galaxyName in l:
                l.append(galaxyName)
                targets[targetName] = l
                
    targetList = targets.keys()
    
    for t in targetList:
        print t
    
    print len(targetList)
    
    summaryList = []
    summaryList.append(('Target',\
    'Galaxy',\
    'R.A.',\
    'Dec.',\
    'z',\
    'Program',\
    'T_exp'))
    
    nameList = ['Target']
    galaxyNameList = ['Galaxy']
    raList = ['R.A.']
    decList = ['Dec.']
    zList = ['z']
    programList = ['Program']
#     gratingList = ['Grating']
#     obsIDList = ['Obs ID']
#     obsDateList = ['Obs Date']
    texpList = ['T_exp [ks]']
#     snList = ['S/N [1238]']
    
    for l in reader:
        target = l['targetName']
        
        if targets.has_key(target):
            galaxies = targets[target]
            
            for galaxy in galaxies:

                ra = l['ra']
                dec = l['dec']
                z = l['z']
#                 sn = l['SNratio']
                programID = l['programID']
#                 programPI = l['programPI']
#                 instrument = l['instrument']
#                 grating = l['grating']
                texp = l['exposureTime']
#                 sn = l['SNexpected']
                
                nameList.append(target)
                galaxyNameList.append(galaxy)
                raList.append(ra)
                decList.append(dec)
                zList.append(z)
                programList.append(programID)
#                 gratingList.append(grating)
#                 obsIDList.append('obsID')
#                 obsDateList.append('Obs Date')
                texpList.append(texp)
#                 snList.append(sn)
                
                summaryList.append((target,\
                galaxy,\
                ra.replace('(','').replace(')','').replace(',',' '),\
                dec.replace('(','').replace(')','').replace(',',' '),\
                z,\
                programID,\
                texp))
                

    padding = 4
    widths =[\
    max(len(str(d)) for d in nameList) + padding,\
    max(len(str(d)) for d in galaxyNameList) + padding,\
    max(len(str(d)) for d in raList) + padding,\
    max(len(str(d)) for d in decList) + padding,\
    max(len(str(d)) for d in zList) + padding,\
    max(len(str(d)) for d in programList) + padding,\
    max(len(str(d)) for d in texpList) + padding]

    for row in summaryList:
        output.write("".join(str(i + '  &').ljust(width) for i,width in zip(row,widths))+'\\\\\n')


#     fieldnames = ('Target',\
#     'Galaxy',\
#     'R.A.',\
#     'Dec.',\
#     'z',\
#     'Program',\
#     'T_exp')
# 
# 
#     d = {'Target':'U29',\
#     'Galaxy':'U29',\
#     'RAdeg':'f8',\
#     'DEdeg':'f8',\
#     'z':'f8',\
#     'Program':'i8',\
#     'T_exp':'i8'}
# 
#     
#     # initiate the table
# #     table = ascii.Table(list(reader),names=fieldnames, dtype=('f4', 'i4', 'S2'))
#     t = table.Table(summaryList, names=fieldnames, dtype=(d['Target'],\
#     d['Galaxy'],\
#     d['RAdeg'],\
#     d['DEdeg'],\
#     d['z'],\
#     d['Program'],\
#     d['T_exp']))
#     
#     
#     ascii.write(t, outName2, format='fixed_width', overwrite=True, delimiter=' & ',\
#     bookend=False,delimiter_pad=None, formats={'Target':'%-29s',\
#     'Galaxy':'%-29s',\
#     'RAdeg':'-.5f',\
#     'DEdeg':'-.5f',\
#     'z':'-.6f',\
#     'Program':'-6d',\
#     'T_exp':'-6d'})


    file.close()
    results.close()
    output.close()

if __name__=="__main__":
    main()
