#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: createGalaxyModelTable_tex4.py, v 4.0 10/26/2018

Create table 2 for the rotation paper -> list of galaxies-target pairs with impact 
parameter, azimuth, velocity, absorption velocity, model velocity ranges

Based on createTargetTable_tex2.py


v4: updated for (hopefully) final version of paper (10/26/18)

v5: updated for referee report (11/11/19)

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
import correlateSingle11_3 as correlateSingle
from utilities3 import *

    
################################################################

def main():
    
    if getpass.getuser() == 'frenchd':
        results_filename = '/Users/frenchd/Research/rotation_paper/salt_galaxy_sightlines_cut_plus_ancillary_fits_newerrs3.csv'
        out_filename = '/Users/frenchd/Research/rotation_paper_data/results_model_table.txt'
        
        models_filename =  '/Users/frenchd/Research/rotation_paper_data/summary_files_4/rotation_model_summary_4.csv'

    else:
        print('Could not determine username. Exiting.')
        sys.exit()
    
    results = open(results_filename,'rU')
    results_reader = csv.DictReader(results)
    
    models = open(models_filename, 'r')
    models_reader = csv.DictReader(models)
    
    output = open(out_filename,'wt')

    summaryList = []
    summaryList.append((\
    'System #',\
    'Galaxy',\
    'Target',\
    '$\\rho$',\
    '$Az.$',\
    '$Inc.$',\
    '$R_{vir}$',\
    '$L_{\**}$',\
    '$v_{sys}$',\
    '$v_{rot}$',\
    '$\Delta v_{Ly\\alpha}$',\
#     '$W_{Ly\\alpha}$',\
    '$b$',\
    '$log N$',\
    '$v_{Steidel}$',\
    '$v_{NFW}$'))
    
    sysNumberList = ['System #']
    targetList = ['Target']
    nameList = ['Galaxy']
    rhoList = ['$\\rho$']
    azList = ['Az.']
    incList = ['Inc.']
    RvirList = ['$R_{vir}$']
    lstarList = ['$L_{\**}$']
    v_sysList = ['$v_{sys}$']
    v_rotList = ['$v_{rot}$']
    lya_v_list = ['$\Delta v_{Ly\\alpha}$']
    lya_b_list = ['$b']
    lya_N_list = ['$log N$']
    steidel_list = ['$v_{Steidel}$']
    nfw_list = ['$v_{NFW}$']


    models_dict = {}
    for m in models_reader:
        Galaxy = m['Galaxy']
        Targetname = m['Targetname']
        Rvir_stocke = m['Rvir_stocke']
        Lstar = m['Lstar']
        inc = m['inc']
        az = m['azimuth']
        vhel = m['Vhel_measured']
        e_vhel = m['e_Vhel_measured']
        vmax = m['vmax']
        e_vmax = m['e_vmax']
        steidel_model = m['steidel_final_vs']
        nfw_model = m['NFW_final_vs']
        
        key = str(Galaxy) + str(Targetname)
        
        models_dict[key] = {'Galaxy':Galaxy,
                            'Targetname':Targetname,
                            'Rvir_stocke':Rvir_stocke,
                            'Lstar':Lstar,
                            'inc':inc,
                            'az':az,
                            'vhel':vhel,
                            'e_vhel':e_vhel,
                            'vmax':vmax,
                            'e_vmax':e_vmax,
                            'steidel_model':steidel_model,
                            'nfw_model':nfw_model}

                
    for l in results_reader:
        include = l['include']
        if include == 'yes' or include == 'maybe':
            # go on with this galaxy
            pass
        else:
            # skip this one
            continue
        
        galaxy = l['Name']
        target = l['Target']
        lya_v = l['fit_v']
        e_lya_v = l['e_fit_v']
        impact = int(round(float(l['impact']),0))
        lya_b = l['fit_b']
        e_lya_b = l['e_fit_b']
        lya_N = l['fit_N']
        e_lya_N = l['e_fit_N']
#         Vhel_measured = l['Vhel_measured']
#         e_Vhel_measured = l['e_Vhel_measured']
#         Vrot_corrected = l['Vrot_corrected']
#         e_Vrot_corrected = l['e_Vrot_corrected']

#         azimuth = int(round(float(l['azimuth']),0))
        sysNumber = l['number']
        impact = l['impact']

        model_key = str(galaxy) + str(target)
        
        model = models_dict[model_key]
        
        Rvir = model['Rvir_stocke']
        Lstar = model['Lstar']
        inc = model['inc']
        az =  model['az']
        vhel = model['vhel']
        e_vhel = model['e_vhel']
        vmax = model['vmax']
        e_vmax = model['e_vmax']
        steidel = eval(model['steidel_model'])
        nfw = eval(model['nfw_model'])
        
        # --- format some values
        az = int(round(float(az),0))
        inc = int(round(float(inc),0))
        vmax = int(round(float(vmax),0))
        e_vmax = int(round(float(e_vmax),0))
        lya_b = round(float(lya_b),1)
        e_lya_b = round(float(e_lya_b),1)
        vhel = int(round(float(vhel),0))
        e_vhel = int(round(float(e_vhel),0))
        lya_N = round(float(lya_N),2)
        e_lya_N = round(float(e_lya_N),2)
        Rvir = int(round(float(Rvir),0))
        impact = int(round(float(impact),0))
        Lstar = round(float(Lstar),1)
        
        steidel_low  = round(float(steidel[0]),1)
        steidel_high  = round(float(steidel[1]),1)
        nfw_low  = round(float(nfw[0]),1)
        nfw_high  = round(float(nfw[1]),1)
        

        steidel = str([steidel_low, steidel_high])
        nfw = str([nfw_low, nfw_high])

        vrot = '${0} \pm {1}$'.format(vmax, e_vmax)
        vsys = '${0} \pm {1}$'.format(vhel, e_vhel)
        b = '${0} \pm {1}$'.format(lya_b, e_lya_b)
        N = '${0} \pm {1}$'.format(lya_N, e_lya_N)
        
        dv = int(round(float(lya_v) - float(vhel),0))
        e_dv = int(round(np.sqrt(float(e_lya_v)**2 + float(e_vhel)**2), 0))
        dv_lya = '${0} \pm {1}$'.format(dv, e_dv)


        # append to the lists
        sysNumberList.append(sysNumber)
        nameList.append(galaxy)
        targetList.append(target)
        rhoList.append(impact)
        azList.append(az)
        incList.append(inc)
        RvirList.append(Rvir)
        lstarList.append(Lstar)
        
        v_sysList.append(vsys)
        v_rotList.append(vrot)
        lya_v_list.append(dv_lya)
        lya_b_list.append(b)
        lya_N_list.append(N)
        steidel_list.append(steidel)
        nfw_list.append(nfw)        

#         model_low = int(Vhel_measured) + int(model_range[0])
#         model_high = int(Vhel_measured) + int(model_range[1])
#         
#         NFW_low = int(Vhel_measured) + int(NFW_range[0])
#         NFW_high = int(Vhel_measured) + int(NFW_range[1])
  
        summaryList.append((sysNumber,\
        galaxy,\
        target,\
        str(impact),\
        str(az),\
        str(inc),\
        str(Rvir),\
        str(Lstar),\
        vsys,\
        vrot,\
        dv_lya,\
        b,\
        N,\
        steidel,
        nfw))
                    

    padding = 4
    widths =[\
    max(len(str(d)) for d in sysNumberList) + padding,\
    max(len(str(d)) for d in nameList) + padding,\
    max(len(str(d)) for d in targetList) + padding,\
    max(len(str(d)) for d in rhoList) + padding,\
    max(len(str(d)) for d in azList) + padding,\
    max(len(str(d)) for d in incList) + padding,\
    max(len(str(d)) for d in RvirList) + padding,\
    max(len(str(d)) for d in lstarList) + padding,\
    max(len(str(d)) for d in v_sysList) + padding,\
    max(len(str(d)) for d in v_rotList) + padding,\
    max(len(str(d)) for d in lya_v_list) + padding,\
    max(len(str(d)) for d in lya_b_list) + padding,\
    max(len(str(d)) for d in lya_N_list) + padding,\
    max(len(str(d)) for d in steidel_list) + padding,\
    max(len(str(d)) for d in nfw_list) + padding]

    for row in summaryList:
    
        stuff = "".join(str(i + '  &').ljust(width) for i,width in zip(row,widths))
        last_amp = stuff.rfind('&')
        
        output.write(stuff[:last_amp] +' \\\\\n')

#         output.write("".join(str(i + '  &').ljust(width) for i,width in zip(row,widths))+' \\\\\n')
        
        
    models.close()
    results.close()
    output.close()

if __name__=="__main__":
    main()
