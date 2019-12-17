#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id: SALT_paper_map4.py, v4 10/25/18

Makes straight on-sky orientation plot as well as cylindrical and NFW model comparison
maps.



Previously: plot_onsky_absorber_vel.py, v2.0 04/04/18

Plot an impact parameter map showing the locations and velocities of each absorber wrt 
the galaxy (2/19/18)

v2: Orient so all the galaxies have approaching side on the left (04/04/18)

v3: Most of the updates (04/30/18)

v3.1: Fixed non-detection labeling issue. THESIS VERSION: (07/15/18)

v4: Final submitted version (hopefully?) 10/25/18

v5: update for the referee report.  (10/14/19)
    Takes in the single file:
    /Research/rotation_paper_data/summary_files/rotation_model_summary.csv
    

'''


import sys
import os
import csv
import time


from pylab import *
# import atpy
# from math import *
from utilities3 import *
from scipy import stats
import getpass
import math
import pickle
import json
import io
import numpy as np
import matplotlib.pyplot as plt
import correlateSingle11_3 as correlateSingle


from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

from scipy.optimize import curve_fit
import matplotlib.ticker as ticker

from matplotlib import rc
fontScale = 14
rc('text', usetex=True)
rc('font', size=fontScale,family='serif',weight='medium')
rc('xtick.major',size=8,width=0.6)
rc('xtick.minor',size=5,width=0.6)
rc('ytick.major',size=8,width=0.6)
rc('ytick.minor',size=5,width=0.6)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
rc('axes',labelsize = fontScale)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
# rc('font', weight = 450)
# rc('axes',labelweight = 'bold')
rc('axes',linewidth = 1,labelweight='normal')
rc('axes',titlesize='small')
rc('xtick',direction='in')
rc('ytick',direction='in')

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# def withinRange(x,range):
#     '''
#         returns True if x is within range
#     '''
#     low = range[0]
#     high = range[1]
#     
#     if x >= low and x <= high:
#         return True
#     else:
#         return False


def withinRange(vel, model_range, error):
    '''
        see if vel +/- error falls within the range of velocities given by
        model_range
        
        returns boolean
    
    '''
    
    
    max_vel = float(vel) + float(error)
    min_vel = float(vel) - float(error)
    
    lower = float(model_range[0])
    upper = float(model_range[1])
    
    answer = False
    if vel >= lower or max_vel >= lower or min_vel >= lower:
        if vel <= upper or max_vel <= upper or min_vel <= upper:
            answer = True
        
    return answer



def get_data(filename):
    # fields in JSON file are:
    #
    # 'name': galaxy name
    # 'vsys_published': Vhel as found on NED
    # 'dvsys_published': error in Vhel as found in NED
    # 'inclination': inclination
    # 'di': inclination error
    # 'centerTrace': center of spectrum
    # 'distance': best distance to galaxy
    # 'vsys_measured': measured Vhel
    # 'vsys_measured_err': measured error in Vhel
    # 'left_vrot_incCorrected_avg': average left wing velocity, corrected for inc
    # 'left_vrot_incCorrected_avg_err': error in average left wing velocity corrected for inc
    # 'right_vrot_incCorrected_avg': average right wing velocity, corrected for inc
    # 'right_vrot_incCorrected_avg_err':  error in average right wing velocity corrected for inc
    # 'left_vrot_avg': average left wing velocity (redshift subtracted)
    # 'left_vrot_avg_err': error in average left wing velocity (redshift subtracted)
    # 'right_vrot_avg': average right wing velocity (redshift subtracted)
    # 'right_vrot_avg_err': error in average right wing velocity (redshift subtracted)
    # 'vrot_vals': observed velocities (redshift but not inc corrected)
    # 'vrot_errs': errors in observed velocities (redshift but not inc corrected)
    # 'vrot_incCorrected_vals': inclination corrected velocities (redshift subtracted)
    # 'vrot_incCorrected_errs': errors in inclination corrected velocities
    # 'vrot_observed': observed velocities (Vhel + rotation)
    # 'agn': include any information about AGN here
    # 'xVals': physical (kpc) x axis along the slit

#     print 'current file: ',filename
#     print
    with open(filename) as data_file:
        data = json.load(data_file)

#         RA_galaxy = data['RAdeg']
#         Dec_galaxy = data['DEdeg']

        Lstar = data['Lstar']
        e_Lstar = data['e_Lstar']
        dist = data['dist']
        right_vrot_incCorrected_avg_err = data['right_vrot_incCorrected_avg_err']
        left_vrot_incCorrected_avg_err = data['left_vrot_incCorrected_avg_err']

        return dist, Lstar, e_Lstar, left_vrot_incCorrected_avg_err, right_vrot_incCorrected_avg_err
        


def main():
    hubbleConstant = 71.0
    
    # where to write to?
#     out_directory = '/Users/frenchd/Research/test/SALT_maps_yes_maybe5/'
#     out_directory = '/Users/frenchd/Research/test/SALT_maps_yes_maybe6/'
#     out_directory = '/Users/frenchd/Research/test/SALT_maps_redo2/'
#     out_directory = '/Users/frenchd/Research/rotation_paper_data/plots/'
#     out_directory = '/Users/frenchd/Research/rotation_paper_data/plots_alt2/'
    out_directory = '/Users/frenchd/Research/rotation_paper_data/plots2/'

#     out_directory = '/Users/frenchd/Research/test/SALT_maps_yes/'

    # only include absorbers that have dv less than or equal to the maximal rotation velocity?
    only_close_velocities = True
    
    # include open circles for sightlines with no absorption detected?
    include_nondetection = True
    
    # what range of Lstar systems to include?
    Lstar_range = [0.0, 100.]
#     Lstar_range = [0.60001, 100.0]
#     Lstar_range = [0.0, 0.6]
#     Lstar_range = [0.0, 0.5]
#     Lstar_range = [0.50001, 100.0]
#     Lstar_range = [0.0, 0.7]
#     Lstar_range = [0.70001, 100.0]
#     Lstar_range = [0.0, 0.8]
#     Lstar_range = [0.80001, 100.0]
#     Lstar_range = [0.0, 0.9]
#     Lstar_range = [0.90001, 100.0]
#     Lstar_range = [0.0, 1.0]
#     Lstar_range = [1.00001, 100.0]

    # if true, makes a L <= 0.6L* note on maps 
    Lstar_note = True

    # azimuth limit for "maybe" trigger
#     az_limit = 85.
    az_limit = 90.

    # size of legend symbols
    legend_size = 12
    
    # size of legend font
    legend_font = 12
    
    # minimum distance to another galaxy
    min_separation = False

    # how far to zoom in for zoom-in plot? Units of R_vir
    zoom_limit = 1.0
    
    # remove systems of lower inclination than this:
#     inc_limit = 75.
    inc_limit = 75.
    
    # --- include e_vhel_measured and Lya v_fit error in the apparent velocity ranges?
    use_apparent_errors = True

    verbose = False
    
    # which plot to make?
    plot_onsky = False
    plot_steidel = False
    plot_NFW = False
    plot_zoom_in = False
    plot_NFW_zoom_in = False
    plot_azmap = True
    
    # include tags to include
#     include_tags = ['yes','maybe', 'no/maybe']
    include_tags = ['yes','maybe']
#     include_tags = ['yes']

##########################################################################################
##########################################################################################
    # the measurement results
#     directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/'
#     filename = '{0}salt_galaxy_sightlines_cut.csv'.format(directory)
#     filename = '{0}salt_galaxy_sightlines_cut_plus_ancillary_fits.csv'.format(directory)  
    directory = '/Users/frenchd/Research/rotation_paper/'

    data_filename = '{0}salt_galaxy_sightlines_cut_plus_ancillary_fits_newerrs3.csv'.format(directory)

    data_file = open(data_filename,'rU')
    data_reader = csv.DictReader(data_file)
    
    # --- get the model results
#     model_directory = '/Users/frenchd/Research/rotation_paper_data/summary_files/'
#     model_filename = '{0}rotation_model_summary.csv'.format(model_directory)

    model_directory = '/Users/frenchd/Research/rotation_paper_data/summary_files_4/'
    model_filename = '{0}rotation_model_summary_4.csv'.format(model_directory)

    model_file = open(model_filename,'rU')
    model_reader = csv.DictReader(model_file)


    # lists to populate
    nameList = []
    targetList = []
    combinedNameList = []
    vList = []
    wList = []
    NList = []
    RA_targetList = []
    Dec_targetList = []
    incList = []
    paList = []
    azList = []
    RvirList = []
    markerColorList = []
    VhelList = []
    markerColorList_steidel = []
    markerColorList_NFW = []
    nameDict = {}
    
    dvList = []
    dv_up_list = []
    dv_down_list = []
    e_dv_list = []
    
    # non-detections/not finished sightlines
    nameList_non = []
    targetList_non = []
    combinedNameList_non = []
    vList_non = []
    wList_non = []
    RA_targetList_non = []
    Dec_targetList_non = []
    incList_non = []
    paList_non = []
    azList_non = []
    RvirList_non = []
    markerColorList_non = []
    VhelList_non = []
    
    apparent_list = []
    
    model_dict = {}
    for m in model_reader:
        Galaxy = m['Galaxy']
        
        Targetname = m['Targetname']
        impact = m['impact']
        azimuth = m['azimuth']
        side = m['side']
        MajDiam = m['MajDiam']
        Rvir = m['Rvir']
        Rvir_stocke = m['Rvir_stocke']
        inc = m['inc']
        e_inc = m['e_inc']
        PA = m['PA']
        e_PA = m['e_PA']
        Vhel_measured = m['Vhel_measured']
        e_Vhel_measured = m['e_Vhel_measured']
        Vhel_published = m['Vhel_published']
        e_Vhel_published = m['e_Vhel_published']
        right_vrot_incCorrected_avg = m['right_vrot_incCorrected_avg']
        left_vrot_incCorrected_avg = m['left_vrot_incCorrected_avg']
        vmax = m['vmax']
        e_vmax = m['e_vmax']
        vmax_incCorrected = m['vmax_incCorrected']
        e_vmax_incCorrected = m['e_vmax_incCorrected']
        V200 = m['V200']
        c = m['c']
        R200 = m['R200']
        V200_min = m['V200_min']
        c_min = m['c_min']
        R200_min = m['R200_min']
        V200_max = m['V200_max']
        c_max = m['c_max']
        R200_max = m['R200_max']
        NFW_final_vs = m['NFW_final_vs']
        NFW_final_ds = m['NFW_final_ds']
        steidel_final_vs = m['steidel_final_vs']
        steidel_final_ds = m['steidel_final_ds']
        min_NFW_vs = m['min_NFW_vs']
        max_NFW_vs = m['max_NFW_vs']
        min_NFW_ds = m['min_NFW_ds']
        max_NFW_ds = m['max_NFW_ds']
        e_min_NFW_vs = m['e_min_NFW_vs']
        e_max_NFW_vs = m['e_max_NFW_vs']
        e_min_NFW_ds = m['e_min_NFW_ds']
        e_max_NFW_ds = m['e_max_NFW_ds']
        min_steidel_vs = m['min_steidel_vs']
        max_steidel_vs = m['max_steidel_vs']
        min_steidel_ds = m['min_steidel_ds']
        max_steidel_ds = m['max_steidel_ds']
        e_min_steidel_vs = m['e_min_steidel_vs']
        e_max_steidel_vs = m['e_max_steidel_vs']
        e_min_steidel_ds = m['e_min_steidel_ds']
        e_max_steidel_ds = m['e_max_steidel_ds']
        zcutoffm = m['zcutoffm']
        rcutoffm = m['rcutoffm']
        spherical_halo = m['spherical_halo']
        z_gradient = m['z_gradient']
        steidel_hv = m['steidel_hv']
        date = m['date']
        
        dictionary_name = Galaxy + Targetname

        model_dict[dictionary_name] = {"Galaxy":Galaxy,
                            "Targetname":Targetname,
                            "impact":impact,
                            "azimuth":azimuth,
                            "side":side,
                            "MajDiam":MajDiam,
                            "Rvir":Rvir,
                            "Rvir_stocke":Rvir_stocke,
                            "inc":inc,
                            "e_inc":e_inc,
                            "PA":PA,
                            "e_PA":e_PA,
                            "Vhel_measured":Vhel_measured,
                            "e_Vhel_measured":e_Vhel_measured,
                            "Vhel_published":Vhel_published,
                            "e_Vhel_published":e_Vhel_published,
                            "right_vrot_incCorrected_avg":right_vrot_incCorrected_avg,
                            "left_vrot_incCorrected_avg":left_vrot_incCorrected_avg,
                            "vmax":vmax,
                            "e_vmax":e_vmax,
                            "vmax_incCorrected":vmax_incCorrected,
                            "e_vmax_incCorrected":e_vmax_incCorrected,
                            "V200":V200,
                            "c":c,
                            "R200":R200,
                            "V200_min":V200_min,
                            "c_min":c_min,
                            "R200_min":R200_min,
                            "V200_max":V200_max,
                            "c_max":c_max,
                            "R200_max":R200_max,
                            "NFW_final_vs":NFW_final_vs,
                            "NFW_final_ds":NFW_final_ds,
                            "steidel_final_vs":steidel_final_vs,
                            "steidel_final_ds":steidel_final_ds,
                            "min_NFW_vs":min_NFW_vs,
                            "max_NFW_vs":max_NFW_vs,
                            "min_NFW_ds":min_NFW_ds,
                            "max_NFW_ds":max_NFW_ds,
                            "e_min_NFW_vs":e_min_NFW_vs,
                            "e_max_NFW_vs":e_max_NFW_vs,
                            "e_min_NFW_ds":e_min_NFW_ds,
                            "e_max_NFW_ds":e_max_NFW_ds,
                            "min_steidel_vs":min_steidel_vs,
                            "max_steidel_vs":max_steidel_vs,
                            "min_steidel_ds":min_steidel_ds,
                            "max_steidel_ds":max_steidel_ds,
                            "e_min_steidel_vs":e_min_steidel_vs,
                            "e_max_steidel_vs":e_max_steidel_vs,
                            "e_min_steidel_ds":e_min_steidel_ds,
                            "e_max_steidel_ds":e_max_steidel_ds,
                            "zcutoffm":zcutoffm,
                            "rcutoffm":rcutoffm,
                            "spherical_halo":spherical_halo,
                            "z_gradient":z_gradient,
                            "steidel_hv":steidel_hv,
                            "date":date}


    for t in data_reader:
        name = t['Name']
        target = t['Target']
        
        dictionary_name = name+target
        
        lyaV = eval(t['Lya_v'])
        lyaW = eval(t['Lya_W'])
        
        # --- from Voigt profile fits
        fit_v   = float(t['fit_v'])
        e_fit_v	= float(t['e_fit_v'])
        fit_b   = float(t['fit_b'])
        e_fit_b = float(t['e_fit_b'])
        fit_N   = float(t['fit_N'])
        e_fit_N = float(t['e_fit_N'])
        
        RA_galaxy = eval(t['RAdeg'])
        Dec_galaxy = eval(t['DEdeg'])
        RA_target = eval(t['RAdeg_target'])
        Dec_target = eval(t['DEdeg_target'])
        vHel = eval(t['Vhel'])
        include_tag = t['include']
        sysNumber = t['number']
        PA_observed = eval(t['PA_observed'])
        PA_adjust = eval(t['PA_adjust'])
        
        print('sysNumber = {}, fit_v = {}'.format(sysNumber, fit_v))
        
        apparent_corotation = t['apparent_corotation']
        
        
        # --- the model velocity ranges for co-roration
        steidel_range   = eval(model_dict[dictionary_name]["steidel_final_vs"])
        NFW_range       = eval(model_dict[dictionary_name]["NFW_final_vs"])


#         gfilename = directory + 'rot_curves/' + name + '-summary4.json'
        gfilename = model_directory + name + '-summary.json'
        dist, Lstar, e_Lstar, left_vrot_incCorrected_avg_err, right_vrot_incCorrected_avg_err = get_data(gfilename)

        Vhel_measured = float(model_dict[dictionary_name]["Vhel_measured"])
        e_Vhel_measured = float(model_dict[dictionary_name]["e_Vhel_measured"])

#         right_vrot_avg = model_dict[name]["right_vrot_avg"]
        right_vrot_incCorrected_avg = float(model_dict[dictionary_name]["right_vrot_incCorrected_avg"])
#         left_vrot_avg = model_dict[name]["left_vrot_avg"]
        left_vrot_incCorrected_avg = float(model_dict[dictionary_name]["left_vrot_incCorrected_avg"])
        
        inc     = float(model_dict[dictionary_name]['inc'])
        PA      = float(model_dict[dictionary_name]['PA'])
        majDiam = float(model_dict[dictionary_name]['MajDiam'])
        impact  = float(model_dict[dictionary_name]['impact'])
        az      = float(model_dict[dictionary_name]['azimuth'])
        Rvir    = float(model_dict[dictionary_name]['Rvir_stocke'])
        
        if verbose:
            print('------------')
            print(name, ' - ', target)
            print('model_dict[dictionary_name]["left_vrot_incCorrected_avg"] = ', model_dict[dictionary_name]["left_vrot_incCorrected_avg"])
            print('left_vrot_incCorrected_avg = ', left_vrot_incCorrected_avg)
            print()
            print(np.sin(inc * np.pi/180.))
            print()

        
        # remove inclination correction to get apparent velocity
        leftVel = left_vrot_incCorrected_avg * np.sin(inc * np.pi/180.)
        rightVel = right_vrot_incCorrected_avg * np.sin(inc * np.pi/180.)
        
        e_leftVel = left_vrot_incCorrected_avg_err * np.sin(inc * np.pi/180.)
        e_rightVel = right_vrot_incCorrected_avg_err * np.sin(inc * np.pi/180.)

        if verbose:
            print('leftVel: ',leftVel)
            print('rightVel: ',rightVel)

        # calculate impact parameter and shit
#         impact = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist)
    
        # RA component of impact parameter - by setting the Dec to be the same for both
        impact_RA = calculateImpactParameter(RA_galaxy, Dec_galaxy, RA_target, Dec_galaxy, dist)
    
        # Dec component of impact parameter - by setting the RA to be the same for both
        impact_Dec = calculateImpactParameter(RA_galaxy, Dec_galaxy, RA_galaxy, Dec_target, dist)
        
        # --- flip around some shit so it's distance in the direction of the AGN     
        if RA_galaxy > RA_target:
            # target is on the 'left' side
            if impact_RA >0:
                impact_RA *= -1
            
        if Dec_galaxy > Dec_target:
            # target is 'below' galaxy
            if impact_Dec >0:
                impact_Dec *= -1
                
        # rotate everything to PA = 90 deg (major axis horizontal)
#         x = 90. - PA
        # rotate everything to PA = 90 deg (major axis horizontal, approaching on left)
#         x = 90. - PA
        x = -PA_adjust

        
        theta = np.radians(x)
        c, s = np.cos(theta), np.sin(theta)
        R = np.matrix('{} {}; {} {}'.format(c, -s, s, c))
        coords = np.array([impact_RA,impact_Dec])
        
        # rotate clockwise
        impact_RA_rot, impact_Dec_rot = np.array(coords.dot(R))[0]
#         print "R: ",R
#         print 'impact_RA_rot: ',impact_RA_rot
#         print 'impact_Dec_rot: ',impact_Dec_rot
        
        if verbose:
            print()
            print('RA_galaxy, Dec_galaxy = ', RA_galaxy, Dec_galaxy)
            print('RA_target, Dec_target = ', RA_target, Dec_target)
            print('impact_RA: ',impact_RA)
            print('impact_Dec: ',impact_Dec)
            print('impact_RA_rot: ',impact_RA_rot)
            print('impact_Dec_rot: ',impact_Dec_rot)
            print()

        # scale to virial radius
        impact_rvir = impact/Rvir
        impact_RA_vir = impact_RA_rot/Rvir
        impact_Dec_vir = impact_Dec_rot/Rvir
                
        # compare to the absorber velocity:
        # negative means the Lya is higher velocity (red)
#         dv = vsys_measured - lyaV

        # switch it around actually, this matches the rotation curve (positive is going
        # away, or higher than systemic velocity gas)
        if fit_v > 0:
            dv = fit_v - Vhel_measured
#             print 'lyaV, dv = ',lyaV,dv
        else:
            print('else. dv = x')
            dv = 'x'


        # check on which 'side' of the galaxy the absorber is found
        if name == 'CGCG039-137':
             # regular - left side of rotation curve is 'left' on sky
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'ESO343-G014':
             # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'IC5325':
             # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'MCG-03-58-009':
             # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC1566':
             # reverse (meaning '<')
             # Changed to regular on 10/30/2019
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC3513':
             # regular
             # Changed to reverse on 10/30/2019
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC3633':
             # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC4536':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC4939':
            # reverse
            # changed to regular on 10/30/2019
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC5364':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC5786':
            # reverse
            # changed to regular on 11/1/2020
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

            print('impact_RA_vir: ',impact_RA_vir)
            print('rot_vel: ',rot_vel)
                
        if name == 'UGC09760':
            # regular
            # changed to reverse on 11/5/2019
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

##########################################################################################
##########################################################################################
          
        if name == 'NGC3198':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC4565':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'UGC04238':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC3351':
            # reverse
            # changed to regular on 11/05/2019
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC4529':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC6140':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC5907':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'UGC06446':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'UGC06399':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC3726':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC3067':
            # regular # changed from reverse on 10/29/20
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC2770':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC3432':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC3666':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC5951':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC7817':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'UGC08146':
            # reverse
            # changed to regular on 10/30/2019
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

        if name == 'NGC3631':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
                e_rot_vel = e_leftVel
            else:
                rot_vel = rightVel
                e_rot_vel = e_rightVel

            
        # now compare to the rotation velocity
        color_yes = '#1b9e77'    # greenish
        color_no = '#d95f02'     # orange
#         color_maybe = '#7570b3'  # blue-ish purple
        color_maybe = 'grey'
        color_nonDetection = 'black'
        markerColor = 'black'
        
        
#         color_yes = '#0652ff'    # blue
#         color_no = '#ff9408'     # orange
        
        color_yes = '#1805db'    # ultramarine blue
        color_yes = '#436bad'      # french blue
        color_no = '#ec2d01'     # tomato red

        
        if isNumber(dv):
            # the absorption and rotation velocity match
#             if (dv > 0 and rot_vel > 0) or (dv < 0 and rot_vel < 0):
#                 markerColor = color_yes
#             
#             # mismatch
#             elif (dv < 0 and rot_vel > 0) or (dv > 0 and rot_vel < 0):
#                 markerColor = color_no
#             
#             else:
#                 markerColor = 'black'

            # velocity difference +/- errors in rotation velocity and Lya line center
#             dv_up = dv + np.sqrt(e_rot_vel**2 + e_fit_v**2)
#             dv_down = dv - np.sqrt(e_rot_vel**2 + e_fit_v**2)

            dv_up = dv + np.sqrt(e_Vhel_measured**2 + e_fit_v**2)
            dv_down = dv - np.sqrt(e_Vhel_measured**2 + e_fit_v**2)

            e_dv = np.sqrt(e_Vhel_measured**2 + e_fit_v**2)
            
            if verbose:
                print('dv_up:  ',dv_up)
                print('dv_down: ',dv_down)
                print('rot_vel: ',rot_vel)
                print()
            
            # --- see if the absorber has the correct sign, taking rotation velocity 
            # --- measurement errors into account
            if use_apparent_errors:
                # take errors into account
                if (dv_up > 0 and rot_vel > 0) or (dv_up < 0 and rot_vel < 0):
                    markerColor = color_yes
                
                elif (dv_down > 0 and rot_vel > 0) or (dv_down < 0 and rot_vel < 0):
                    markerColor = color_yes
            
                # mismatch
                elif (dv_up < 0 and rot_vel > 0) or (dv_up > 0 and rot_vel < 0):
                    markerColor = color_no
                
                elif (dv_down < 0 and rot_vel > 0) or (dv_down > 0 and rot_vel < 0):
                    markerColor = color_no
            
                else:
                    markerColor = 'black'

            # --- see if the absorber has the correct sign, ignoring rotation velocity
            # --- measurement errors
            else:
                # the absorption and rotation velocity match
                if (dv > 0 and rot_vel > 0) or (dv < 0 and rot_vel < 0):
                    markerColor = color_yes
            
                # mismatch
                elif (dv < 0 and rot_vel > 0) or (dv > 0 and rot_vel < 0):
                    markerColor = color_no
            
                else:
                    markerColor = 'black'


            # compare to the models
            model_errors = np.sqrt(e_fit_v**2 + e_Vhel_measured**2)
            steidel_answer = withinRange(dv, steidel_range, model_errors)
            NFW_answer = withinRange(dv, NFW_range, model_errors)
            
            if markerColor == color_yes:
                apparent_answer  = True
            else:
                apparent_answer = False
            
            if verbose:
                print('model_errors: ',model_errors)
                print('steidel_answer: ',steidel_answer)
                print('dv, steidel_range, model_errors: ',dv, steidel_range, model_errors)
                print('----')
                print('NFW_answer: ', NFW_answer)
                print('dv, NFW_range, model_errors: ',dv, NFW_range, model_errors)
                print('-----')
                print('Apparent answer: ',apparent_answer)
                print('apparent_corotation: ',apparent_corotation)
                print()
        
            if steidel_answer:
                markerColor_steidel = color_yes
            else:
                markerColor_steidel = color_no
            
            if NFW_answer:
                markerColor_NFW = color_yes
            else:
                markerColor_NFW = color_no
            
            if az >= az_limit:
                markerColor_steidel = color_maybe
                markerColor_NFW = color_maybe
                
            if verbose:
                print('dv = {0} - markerColor_NFW = {1}'.format(dv, markerColor_NFW))
                print()
            
        else:
            print('dv == x :', dv, name)
            print()
            markerColor = color_nonDetection
        
        
        # if too close to the minor axis, call it uncertain/maybe
        if az >= az_limit:
            markerColor = color_maybe
            markerColor_steidel = color_maybe
            markerColor_NFW = color_maybe
            
        if name == 'NGC3513' or name == 'NGC4536':
            markerColor = color_maybe
            markerColor_steidel = color_maybe
            markerColor_NFW = color_maybe
            

        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = False

#         combinedName = '${\rm '+name + '-' + target+'}$'.format(name,target)
        combinedName = r'$\rm {0} : {1}$'.format(name,target)
        
        if sysNumber not in nameDict:
            if isNumber(sysNumber):
                nameDict[combinedName] = sysNumber
        
#         print '{0} - dv={1} vs model={2} => {3}'.format(combinedName, dv, model_range, model_answer)
#         print

            
        # decide some things
        if withinRange(Lstar, Lstar_range, 0.0):
            add_to_list = True
        else:
            add_to_list = False
            
        if inc <= inc_limit:
            add_to_list = False
            
        if include_tag not in include_tags:
            add_to_list = False
        
        separation = 500
        if min_separation and add_to_list:
            # correlate with environment
            agnSeparation = False
            minVcorr = False
            minSize = False
            correlation = correlateSingle.correlateTarget(name, min_separation, agnSeparation, minVcorr, minSize, slow=False, searchAll=True)
            galaxyInfo = correlation[name]
            
            if verbose:
                print('galaxyInfo: ',galaxyInfo)
                
            for row in galaxyInfo:
                vhel, galaxyRow = row
                separation = galaxyRow['impactParameter (kpc)']
                galaxyVel = galaxyRow['radialVelocity (km/s)']
                
                if verbose:
                    print('separation: ',separation)
                    print('galaxyVel: ',galaxyVel)
                    print('vHel: ',vHel)
                
                    print('withinRange(galaxyVel, [vHel-400, vHel+400], 0.0): ',withinRange(galaxyVel, [vHel-400, vHel+400], 0.0))
                    print()

                if withinRange(galaxyVel, [vHel-400, vHel+400], 0.0) and add_to_list:
                    if separation <= min_separation and separation >0.0:
                        add_to_list = False
                        print('False for {0} - {1}'.format(name, separation))
                        print()
                        break
        
        
        # populate the lists
        if add_to_list:
            # first detections
            if isNumber(dv) and only_close_velocities:
                if abs(dv) <= (abs(rot_vel) + e_rot_vel):

                    nameList.append(name)
                    targetList.append(target)
                    vList.append(fit_v)
                    wList.append(lyaW)
                    NList.append(fit_N)
                    RA_targetList.append(impact_RA_vir)
                    Dec_targetList.append(impact_Dec_vir)
                    incList.append(inc)
                    paList.append(PA)
                    azList.append(az)
                    RvirList.append(Rvir)
                    markerColorList.append(markerColor)
                    combinedNameList.append(combinedName)
                    markerColorList_NFW.append(markerColor_NFW)
                    markerColorList_steidel.append(markerColor_steidel)
                    
                    dvList.append(dv)
                    dv_up_list.append(dv_up)
                    dv_down_list.append(dv_down)
                    e_dv_list.append(e_dv)
                    
                    apparent_list.append(apparent_corotation)

                else:
                    print('too far: ',name,' , dv = ',dv, ' vs rot_vel = ',rot_vel)
                    print()
                
            else:
                if isNumber(dv):
                    nameList.append(name)
                    targetList.append(target)
                    vList.append(fit_v)
                    wList.append(lyaW)
                    NList.append(fit_N)
                    RA_targetList.append(impact_RA_vir)
                    Dec_targetList.append(impact_Dec_vir)
                    incList.append(inc)
                    paList.append(PA)
                    azList.append(az)
                    RvirList.append(Rvir)
                    markerColorList.append(markerColor)
                    combinedNameList.append(combinedName)
                    markerColorList_NFW.append(markerColor_NFW)
                    markerColorList_steidel.append(markerColor_steidel)                

                    dvList.append(dv)
                    dv_up_list.append(dv_up)
                    dv_down_list.append(dv_down)
                    e_dv_list.append(e_dv)

                    apparent_list.append(apparent_corotation)

                
            if include_nondetection and not isNumber(dv):
                # include non-detections?
                nameList_non.append(name)
                targetList_non.append(target)
                vList_non.append(lyaV)
                wList_non.append(lyaW)
                RA_targetList_non.append(impact_RA_vir)
                Dec_targetList_non.append(impact_Dec_vir)
                incList_non.append(inc)
                paList_non.append(PA)
                azList_non.append(az)
                RvirList_non.append(Rvir)
                markerColorList_non.append(markerColor)
                combinedNameList_non.append(combinedName)
        
        else:
            print('outside range: {0} - {1} Lstar, {2} = separation'.format(Lstar, Lstar_range, separation))
            print()

    model_file.close()
    data_file.close()
    
    nos = 0
    for i in apparent_list:
        if i == 'no':
            nos +=1
    
    yess = len(apparent_list) - nos
    
    if verbose:
        print('Len(apparent_list)  : ', len(apparent_list))

        print('yes/total apparent list = ',yess/len(apparent_list))
        print()

##########################################################################################
##########################################################################################
##########################################################################################
    # sort in order of largest to smallest column density
    orderedList = []
    for ra, dec, c, N, name, steidel, NFW, dv in zip(RA_targetList, Dec_targetList, markerColorList, NList, combinedNameList, markerColorList_steidel, markerColorList_NFW, dvList):
        orderedList.append([N,[ra, dec, c, name, steidel, NFW, dv]])
        
    orderedList.sort(reverse=True)
    
    
    RA_targetList2 = []
    Dec_targetList2 = []
    markerColorList2 = []
    NList2 = []
    combinedNameList2 = []
    markerColorList_NFW2 = []
    markerColorList_steidel2 = []
    countList = []
    dvList2 = []
    
    # zoomed in
    RA_targetList_zoom = []
    Dec_targetList_zoom = []
    markerColorList_zoom = []
    NList_zoom = []
    combinedNameList_zoom = []
    markerColorList_NFW_zoom = []
    markerColorList_steidel_zoom = []
    countList_zoom = []
    dvList_zoom = []
    
    
    count = 1
    count_zoom = 1
    count_dict = {}
    count_dict_zoom = {}
    for i in orderedList:
        N, rest = i
        ra, dec, c, name, steidel, NFW, dv = rest
        
        RA_targetList2.append(ra)
        Dec_targetList2.append(dec)
        markerColorList2.append(c)
        NList2.append(N)
        combinedNameList2.append(name)
        markerColorList_NFW2.append(NFW)
        markerColorList_steidel2.append(steidel)
        dvList2.append(dv)
        
        countList.append(nameDict[name])
        print('name = {}, nameDict[name] = {}'.format(name, nameDict[name]))
            
        # separate, 'zoomed-in' set
        if math.sqrt(ra**2 + dec**2) <= zoom_limit:
            RA_targetList_zoom.append(ra)
            Dec_targetList_zoom.append(dec)
            markerColorList_zoom.append(c)
            NList_zoom.append(N)
            combinedNameList_zoom.append(name)
            markerColorList_NFW_zoom.append(NFW)
            markerColorList_steidel_zoom.append(steidel)
            dvList_zoom.append(dv)
        
            # check if this galaxy-QSO pair already has a number
            countList_zoom.append(nameDict[name])
        
    countList_non = []
    countList_non_zoom = []
    RA_targetList_non_zoom = []
    Dec_targetList_non_zoom = []
    markerColorList_non_zoom = []
    combinedNameList_non_zoom = []
    for name, ra, dec, m in zip(combinedNameList_non, RA_targetList_non, Dec_targetList_non, markerColorList_non):
    
        countList_non.append(nameDict[name])
        
        if math.sqrt(ra**2 + dec**2) <= zoom_limit:
            countList_non_zoom.append(nameDict[name])
            
            RA_targetList_non_zoom.append(ra)
            Dec_targetList_non_zoom.append(dec)
            markerColorList_non_zoom.append(m)
            combinedNameList_non_zoom.append(name)

##########################################################################################
##########################################################################################
# simple apparent on-sky velocity plot
#
#
##########################################################################################
##########################################################################################

    if plot_onsky:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
    #     fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)

    ##########################################################################################
        # plot circles
        def xy(r,phi):
          return r*np.cos(phi), r*np.sin(phi)

        phis=np.arange(0,2*np.pi,0.01)
    
        r = 1.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        r = 2.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        r = 3.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
#         r = 4.0
#         ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
#         ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

    
        # plot the rest
        largest_N = max(NList2)
        smallest_N = min(NList2)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for N in NList2:
            newSize = ((float(N) - smallest_N)/(largest_N - smallest_N)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # just plot them all the same
    #     ax.scatter(RA_targetList2,Dec_targetList2, color=markerColorList2,s=newSizeList)

        # make different style markers for different colors
        for i in arange(len(markerColorList2)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList2[i] == color_maybe:
                marker = 'o'
            if markerColorList2[i] == color_no:
                marker = 'X'
                marker_lw = 0.5
            if markerColorList2[i] == color_yes:
                marker = 'D'
        

            ax.scatter(RA_targetList2[i], Dec_targetList2[i], color=markerColorList2[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)


        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)
    #         print 'combinedNameList2[i]: ',combinedNameList2[i], 'dec = ',Dec_targetList2[i]
    #         print 'yTagOffset: ',yTagOffset

#             annotate(countList[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
#             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
                
            tag = nameDict[combinedNameList2[i]]
            if not tag in previousNames:
                print('tag: ',tag)
                annotate(tag,xy=(RA_targetList2[i], Dec_targetList2[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1


    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non)):
            ax.plot(RA_targetList_non[i], Dec_targetList_non[i], color=markerColorList_non[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
#             annotate(countList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        
            tag = nameDict[combinedNameList_non[i]]
        
            if not tag in previousNames:
                print('tag non: ',tag)
                annotate(tag,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1

    ##########################################################################################

        xlabel(r'$\rm x ~[R_{vir}]$')
        ylabel(r'$\rm y ~[R_{vir}]$')

        ax.set_xlim(-3.0, 3.0)
        ax.set_ylim(-3.0, 3.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(2.97, 0.06),\
        xytext=(0.0,0.0),textcoords='offset points',size=9)

        annotate(r'$\rm Apparent~Vel.$', xy=(2.90, 2.72),\
        xytext=(0.0,0.0), textcoords='offset points', size=12, weight='bold')
        
        if only_close_velocities:
            annotate(r'$\rm | \Delta v | \leq v_{\rm rot}$', xy=(2.90, 2.52),\
            xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
        
            if Lstar_note:
                annotate(r'$\rm L \leq 0.5 L^{\**}$', xy=(2.9, 2.32),\
                xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
            


#         x-axis
#         majorLocator   = MultipleLocator(1)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(0.5)
#         ax.xaxis.set_major_locator(majorLocator)
#         ax.xaxis.set_major_formatter(majorFormatter)
#         ax.xaxis.set_minor_locator(minorLocator)
# 
#         y axis
#         majorLocator   = MultipleLocator(1)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(0.5)
#         ax.yaxis.set_major_locator(majorLocator)
#         ax.yaxis.set_major_formatter(majorFormatter)
#         ax.yaxis.set_minor_locator(minorLocator)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
        plt.legend(handles=[corotate, antirotate, nondetection],loc='lower left', 
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALTmap_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_inclim_{4}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation, inc_limit)
    #     savefig("{0}SALT_map1.pdf".format(out_directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        summary_filename = '{0}/{1}_summary.txt'.format(out_directory, save_name)
        summary_file = open(summary_filename,'wt')
    
    #     for k, v in sorted(previousNames.iteritems(), key=lambda (k,v): (v,k)):
    #         summary_file.write('{0}. {1}, \n'.format(v,k))
        
        for num, name in zip(countList, combinedNameList2):
            summary_file.write('{0}. {1}, \n'.format(num,name))
        
        for num, name in zip(countList_non, combinedNameList_non):
            summary_file.write('{0}. {1}, \n'.format(num,name))      
            
        
        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        combinedNameList_corotate = []
        combinedNameList_antirotate = []
        combinedNameList_maybe = []

        combinedNameList_corotate_zoom_in = []
        combinedNameList_antirotate_zoom_in = []
        combinedNameList_maybe_zoom_in = []
        
        combinedNameList_corotate_zoom_out = []
        combinedNameList_antirotate_zoom_out = []
        combinedNameList_maybe_zoom_out = []
        
        for name, ra, dec, c in zip(combinedNameList2, RA_targetList2, Dec_targetList2, markerColorList2):
            if c == color_yes:
                combinedNameList_corotate.append(name)
            if c == color_maybe:
                combinedNameList_maybe.append(name)
            if c == color_no:
                combinedNameList_antirotate.append(name)
                
            if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                if c == color_yes:
                    combinedNameList_corotate_zoom_in.append(name)
                if c == color_maybe:
                    combinedNameList_maybe_zoom_in.append(name)
                if c == color_no:
                    combinedNameList_antirotate_zoom_in.append(name)
                    
            else:
                if c == color_yes:
                    combinedNameList_corotate_zoom_out.append(name)
                if c == color_maybe:
                    combinedNameList_maybe_zoom_out.append(name)
                if c == color_no:
                    combinedNameList_antirotate_zoom_out.append(name)
                
                
        stats_file.write('Total co-rotate = {0} \n'.format(len(combinedNameList_corotate)))
        stats_file.write('Total anti-rotate = {0} \n'.format(len(combinedNameList_antirotate)))
        stats_file.write('Total uncertain = {0} \n'.format(len(combinedNameList_maybe)))
        stats_file.write('\n')
        stats_file.write('Total co-rotate within {0} = {1} \n'.format(zoom_limit, len(combinedNameList_corotate_zoom_in)))
        stats_file.write('Total co-rotate outside {0} = {1} \n'.format(zoom_limit, len(combinedNameList_corotate_zoom_out)))
        stats_file.write('\n')
        stats_file.write('Total anti-rotate within {0} = {1} \n'.format(zoom_limit, len(combinedNameList_antirotate_zoom_in)))
        stats_file.write('Total anti-rotate outside {0} = {1} \n'.format(zoom_limit, len(combinedNameList_antirotate_zoom_out)))
        stats_file.write('\n')
        stats_file.write('Total uncertain within = {0} \n'.format(zoom_limit, len(combinedNameList_maybe_zoom_in)))
        stats_file.write('Total uncertain outside = {0} \n'.format(zoom_limit, len(combinedNameList_maybe_zoom_out)))
        stats_file.write('\n')
        
        stats_file.write('Co-rotating systems: \n')
        for i in combinedNameList_corotate:
            stats_file.write('{0}\n'.format(i))
        
        stats_file.write('\n')
        stats_file.write('Anti-rotating systems: \n')
        for i in combinedNameList_antirotate:
            stats_file.write('{0}\n'.format(i))
    
        stats_file.write('\n')
        stats_file.write('Uncertain systems: \n')
        for i in combinedNameList_maybe:
            stats_file.write('{0}\n'.format(i))
            
            
        stats_file.close()
        summary_file.close()
    
    
##########################################################################################
##########################################################################################
# steidel Model plot
#
#
##########################################################################################
##########################################################################################

    if plot_steidel:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
    #     fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)

    ##########################################################################################
        # plot circles
        def xy(r,phi):
          return r*np.cos(phi), r*np.sin(phi)

        phis=np.arange(0,2*np.pi,0.01)
    
        r = 1.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        r = 2.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        r = 3.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
#         r = 4.0
#         ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
#         ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

    
        # plot the rest
        largest_N = max(NList2)
        smallest_N = min(NList2)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for N in NList2:
            newSize = ((float(N) - smallest_N)/(largest_N - smallest_N)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # make different style markers for different colors
        for i in arange(len(markerColorList_steidel2)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList_steidel2[i] == color_maybe:
                marker = 'o'
            if markerColorList_steidel2[i] == color_no:
                marker = 'X'
                marker_lw = 0.5
            if markerColorList_steidel2[i] == color_yes:
                marker = 'D'
        

            ax.scatter(RA_targetList2[i], Dec_targetList2[i], color=markerColorList_steidel2[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    
    
        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)
    #         print 'combinedNameList2[i]: ',combinedNameList2[i], 'dec = ',Dec_targetList2[i]
    #         print 'yTagOffset: ',yTagOffset
                
#             annotate(countList[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
#             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
                
            tag = nameDict[combinedNameList2[i]]
            if not tag in previousNames:
                annotate(tag,xy=(RA_targetList2[i], Dec_targetList2[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1


    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non)):
            ax.plot(RA_targetList_non[i], Dec_targetList_non[i], color=markerColorList_non[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
#             annotate(countList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        
            tag = nameDict[combinedNameList_non[i]]
        
            if not tag in previousNames:
                annotate(tag,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1

    ##########################################################################################

        xlabel(r'$\rm x ~[R_{vir}]$')
        ylabel(r'$\rm y ~[R_{vir}]$')
    
#         ax.set_xlim(-4.0, 4.0)
#         ax.set_ylim(-4.0, 4.0)
#         ax.invert_xaxis()
#     
#         annotate(r'$\rm Approaching~ Side$',xy=(3.95, 0.06),\
#         xytext=(0.0,0.0),textcoords='offset points',size=9)

        ax.set_xlim(-3.0, 3.0)
        ax.set_ylim(-3.0, 3.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(2.97, 0.06),\
        xytext=(0.0,0.0),textcoords='offset points',size=9)

        annotate(r'$\rm Steidel.~Model$', xy=(2.90, 2.72),\
        xytext=(0.0,0.0), textcoords='offset points', size=12, weight='bold')
        
        if only_close_velocities:
            annotate(r'$\rm | \Delta v | \leq v_{\rm rot}$', xy=(2.90, 2.52),\
            xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
        
            if Lstar_note:
                annotate(r'$\rm L \leq 0.5 L^{\**}$', xy=(2.90, 2.32),\
                xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
            
        # x-axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)

        # y axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
        plt.legend(handles=[corotate, antirotate, nondetection],loc='lower left', 
                                borderpad=0.8, fontsize=legend_font, fancybox=True)


    ##########################################################################################
        
        save_name = 'SALTmap_steidel_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_inclim_{4}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation, inc_limit)    
    
    #     savefig("{0}SALT_map1.pdf".format(out_directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        summary_filename = '{0}/{1}_summary.txt'.format(out_directory, save_name)
        summary_file = open(summary_filename,'wt')
    
    #     for k, v in sorted(previousNames.iteritems(), key=lambda (k,v): (v,k)):
    #         summary_file.write('{0}. {1}, \n'.format(v,k))

        for num, name in zip(countList, combinedNameList2):
            summary_file.write('{0}. {1}, \n'.format(num,name))
        
        for num, name in zip(countList_non, combinedNameList_non):
            summary_file.write('{0}. {1}, \n'.format(num,name))


        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        combinedNameList_corotate = []
        combinedNameList_antirotate = []
        combinedNameList_maybe = []

        combinedNameList_corotate_zoom_in = []
        combinedNameList_antirotate_zoom_in = []
        combinedNameList_maybe_zoom_in = []
        
        combinedNameList_corotate_zoom_out = []
        combinedNameList_antirotate_zoom_out = []
        combinedNameList_maybe_zoom_out = []
        
        for name, ra, dec, c in zip(combinedNameList2, RA_targetList2, Dec_targetList2, markerColorList_steidel2):
            if c == color_yes:
                combinedNameList_corotate.append(name)
            if c == color_maybe:
                combinedNameList_maybe.append(name)
            if c == color_no:
                combinedNameList_antirotate.append(name)
                
            if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                if c == color_yes:
                    combinedNameList_corotate_zoom_in.append(name)
                if c == color_maybe:
                    combinedNameList_maybe_zoom_in.append(name)
                if c == color_no:
                    combinedNameList_antirotate_zoom_in.append(name)
                    
            else:
                if c == color_yes:
                    combinedNameList_corotate_zoom_out.append(name)
                if c == color_maybe:
                    combinedNameList_maybe_zoom_out.append(name)
                if c == color_no:
                    combinedNameList_antirotate_zoom_out.append(name)
                
                
        stats_file.write('Total co-rotate = {0} \n'.format(len(combinedNameList_corotate)))
        stats_file.write('Total anti-rotate = {0} \n'.format(len(combinedNameList_antirotate)))
        stats_file.write('Total uncertain = {0} \n'.format(len(combinedNameList_maybe)))
        stats_file.write('\n')
        stats_file.write('Total co-rotate within {0} = {1} \n'.format(zoom_limit, len(combinedNameList_corotate_zoom_in)))
        stats_file.write('Total co-rotate outside {0} = {1} \n'.format(zoom_limit, len(combinedNameList_corotate_zoom_out)))
        stats_file.write('\n')
        stats_file.write('Total anti-rotate within {0} = {1} \n'.format(zoom_limit, len(combinedNameList_antirotate_zoom_in)))
        stats_file.write('Total anti-rotate outside {0} = {1} \n'.format(zoom_limit, len(combinedNameList_antirotate_zoom_out)))
        stats_file.write('\n')
        stats_file.write('Total uncertain within = {0} \n'.format(zoom_limit, len(combinedNameList_maybe_zoom_in)))
        stats_file.write('Total uncertain outside = {0} \n'.format(zoom_limit, len(combinedNameList_maybe_zoom_out)))
        stats_file.write('\n')
        
        stats_file.write('\n')
        stats_file.write('Co-rotating systems: \n')
        for i in combinedNameList_corotate:
            stats_file.write('{0}\n'.format(i))
        
        stats_file.write('\n')
        stats_file.write('Anti-rotating systems: \n')
        for i in combinedNameList_antirotate:
            stats_file.write('{0}\n'.format(i))
    
        stats_file.write('\n')
        stats_file.write('Uncertain systems: \n')
        for i in combinedNameList_maybe:
            stats_file.write('{0}\n'.format(i))
            
            
        stats_file.close()
        summary_file.close()


##########################################################################################
##########################################################################################
# NFW model plot
#
#
##########################################################################################
##########################################################################################

    if plot_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
    #     fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)

    ##########################################################################################
        # plot circles
        def xy(r,phi):
          return r*np.cos(phi), r*np.sin(phi)

        phis=np.arange(0,2*np.pi,0.01)
    
        r = 1.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        r = 2.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        r = 3.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
#         r = 4.0
#         ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
#         ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

    
        # plot the rest
        largest_N = max(NList2)
        smallest_N = min(NList2)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for N in NList2:
            newSize = ((float(N) - smallest_N)/(largest_N - smallest_N)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # make different style markers for different colors
        for i in arange(len(markerColorList_NFW2)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList_NFW2[i] == color_maybe:
                marker = 'o'
            if markerColorList_NFW2[i] == color_no:
                marker = 'X'
                marker_lw = 0.5
            if markerColorList_NFW2[i] == color_yes:
                marker = 'D'
        

            ax.scatter(RA_targetList2[i], Dec_targetList2[i], color=markerColorList_NFW2[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    
    
        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)
    #         print 'combinedNameList2[i]: ',combinedNameList2[i], 'dec = ',Dec_targetList2[i]
    #         print 'yTagOffset: ',yTagOffset
                
#             annotate(countList[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
#             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
                
            tag = nameDict[combinedNameList2[i]]
            if not tag in previousNames:
                annotate(tag,xy=(RA_targetList2[i], Dec_targetList2[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1


    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non)):
            ax.plot(RA_targetList_non[i], Dec_targetList_non[i], color=markerColorList_non[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
#             annotate(countList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        
            tag = nameDict[combinedNameList_non[i]]
        
            if not tag in previousNames:
                annotate(tag,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1


    ##########################################################################################

        xlabel(r'$\rm x ~[R_{vir}]$')
        ylabel(r'$\rm y ~[R_{vir}]$')
    
#         ax.set_xlim(-4.0, 4.0)
#         ax.set_ylim(-4.0, 4.0)
#         ax.invert_xaxis()
#     
#         annotate(r'$\rm Approaching~ Side$',xy=(3.95, 0.06),\
#         xytext=(0.0,0.0),textcoords='offset points',size=9)

        ax.set_xlim(-3.0, 3.0)
        ax.set_ylim(-3.0, 3.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(2.97, 0.06),\
        xytext=(0.0,0.0),textcoords='offset points',size=9)

        annotate(r'$\rm NFW~Model$', xy=(2.90, 2.72),\
        xytext=(0.0,0.0), textcoords='offset points', size=12, weight='bold')
        
        if only_close_velocities:
            annotate(r'$\rm | \Delta v | \leq v_{\rm rot}$', xy=(2.90, 2.52),\
            xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
        
            if Lstar_note:
                annotate(r'$\rm L \leq 0.5 L^{\**}$', xy=(2.90, 2.32),\
                xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
            
        # x-axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)

        # y axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
        plt.legend(handles=[corotate, antirotate, nondetection],loc='lower left', 
                                borderpad=0.8, fontsize=legend_font, fancybox=True)


    ##########################################################################################
        
        save_name = 'SALTmap_NFW_model_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_inclim_{4}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation, inc_limit)

    #     savefig("{0}SALT_map1.pdf".format(out_directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        summary_filename = '{0}/{1}_summary.txt'.format(out_directory, save_name)
        summary_file = open(summary_filename,'wt')
    
    #     for k, v in sorted(previousNames.iteritems(), key=lambda (k,v): (v,k)):
    #         summary_file.write('{0}. {1}, \n'.format(v,k))

        for num, name in zip(countList, combinedNameList2):
            summary_file.write('{0}. {1}, \n'.format(num,name))
        
        for num, name in zip(countList_non, combinedNameList_non):
            summary_file.write('{0}. {1}, \n'.format(num,name))
    
        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        combinedNameList_corotate = []
        combinedNameList_antirotate = []
        combinedNameList_maybe = []

        combinedNameList_corotate_zoom_in = []
        combinedNameList_antirotate_zoom_in = []
        combinedNameList_maybe_zoom_in = []
        
        combinedNameList_corotate_zoom_out = []
        combinedNameList_antirotate_zoom_out = []
        combinedNameList_maybe_zoom_out = []
        
        for name, ra, dec, c in zip(combinedNameList2, RA_targetList2, Dec_targetList2, markerColorList_NFW2):
            if c == color_yes:
                combinedNameList_corotate.append(name)
            if c == color_maybe:
                combinedNameList_maybe.append(name)
            if c == color_no:
                combinedNameList_antirotate.append(name)
                
            if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                if c == color_yes:
                    combinedNameList_corotate_zoom_in.append(name)
                if c == color_maybe:
                    combinedNameList_maybe_zoom_in.append(name)
                if c == color_no:
                    combinedNameList_antirotate_zoom_in.append(name)
                    
            else:
                if c == color_yes:
                    combinedNameList_corotate_zoom_out.append(name)
                if c == color_maybe:
                    combinedNameList_maybe_zoom_out.append(name)
                if c == color_no:
                    combinedNameList_antirotate_zoom_out.append(name)
                
                
        stats_file.write('Total co-rotate = {0} \n'.format(len(combinedNameList_corotate)))
        stats_file.write('Total anti-rotate = {0} \n'.format(len(combinedNameList_antirotate)))
        stats_file.write('Total uncertain = {0} \n'.format(len(combinedNameList_maybe)))
        stats_file.write('\n')
        stats_file.write('Total co-rotate within {0} = {1} \n'.format(zoom_limit, len(combinedNameList_corotate_zoom_in)))
        stats_file.write('Total co-rotate outside {0} = {1} \n'.format(zoom_limit, len(combinedNameList_corotate_zoom_out)))
        stats_file.write('\n')
        stats_file.write('Total anti-rotate within {0} = {1} \n'.format(zoom_limit, len(combinedNameList_antirotate_zoom_in)))
        stats_file.write('Total anti-rotate outside {0} = {1} \n'.format(zoom_limit, len(combinedNameList_antirotate_zoom_out)))
        stats_file.write('\n')
        stats_file.write('Total uncertain within = {0} \n'.format(zoom_limit, len(combinedNameList_maybe_zoom_in)))
        stats_file.write('Total uncertain outside = {0} \n'.format(zoom_limit, len(combinedNameList_maybe_zoom_out)))
        stats_file.write('\n')
        
        stats_file.write('\n')
        stats_file.write('Co-rotating systems: \n')
        for i in combinedNameList_corotate:
            stats_file.write('{0}\n'.format(i))
        
        stats_file.write('\n')
        stats_file.write('Anti-rotating systems: \n')
        for i in combinedNameList_antirotate:
            stats_file.write('{0}\n'.format(i))
    
        stats_file.write('\n')
        stats_file.write('Uncertain systems: \n')
        for i in combinedNameList_maybe:
            stats_file.write('{0}\n'.format(i))
            
            
        stats_file.close()
        summary_file.close()



##########################################################################################
##########################################################################################
# On-sky apparent plot but zoomed into 1 R_vir radius only
#
#
##########################################################################################
##########################################################################################

    if plot_zoom_in:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
    #     fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)

    ##########################################################################################
        # plot circles
        def xy(r,phi):
          return r*np.cos(phi), r*np.sin(phi)

        phis=np.arange(0,2*np.pi,0.01)
    
        r = 0.5
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        r = 1.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-1,1],c='black',ls='-',lw=0.6)
        ax.plot([-1,1],[0,0],c='black',ls='-',lw=0.6)
    
#         ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

    
        # plot the rest
        largest_N = max(NList_zoom)
        smallest_N = min(NList_zoom)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for N in NList_zoom:
            newSize = ((float(N) - smallest_N)/(largest_N - smallest_N)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # make different style markers for different colors
        for i in arange(len(markerColorList_zoom)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList_zoom[i] == color_maybe:
                marker = 'o'
            if markerColorList_zoom[i] == color_no:
                marker = 'X'
#                 marker_lw = 1.5
                marker_lw = 0.5
            if markerColorList_zoom[i] == color_yes:
                marker = 'D'
        

            ax.scatter(RA_targetList_zoom[i], Dec_targetList_zoom[i], color=markerColorList_zoom[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    
        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList_zoom)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)

#             annotate(countList_zoom[i],xy=(RA_targetList_zoom[i], Dec_targetList_zoom[i]),\
#             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)

            tag = nameDict[combinedNameList_zoom[i]]        
            if not tag in previousNames:
                annotate(tag,xy=(RA_targetList_zoom[i], Dec_targetList_zoom[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1

    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non_zoom)):
            ax.plot(RA_targetList_non_zoom[i], Dec_targetList_non_zoom[i], color=markerColorList_non_zoom[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
#             annotate(countList_non_zoom[i],xy=(RA_targetList_non_zoom[i], Dec_targetList_non_zoom[i]),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
            
            tag = nameDict[combinedNameList_non_zoom[i]]
        
            if not tag in previousNames:
                annotate(tag,xy=(RA_targetList_non_zoom[i], Dec_targetList_non_zoom[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1

    ##########################################################################################

        xlabel(r'$\rm x ~[R_{vir}]$')
        ylabel(r'$\rm y ~[R_{vir}]$')
    
        ax.set_xlim(-zoom_limit, zoom_limit)
        ax.set_ylim(-zoom_limit, zoom_limit)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(zoom_limit-0.06, 0.02),\
        xytext=(0.0,0.0),textcoords='offset points',size=9)
        
        annotate(r'$\rm Apparent~Vel.~(zoom)$', xy=(0.96, 0.91),\
        xytext=(0.0,0.0), textcoords='offset points', size=12, weight='bold')
        
        if only_close_velocities:
            annotate(r'$\rm | \Delta v | \leq v_{\rm rot}$', xy=(0.96, 0.84),\
            xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
        
            if Lstar_note:
                annotate(r'$\rm L \leq 0.5 L^{\**}$', xy=(0.96, 0.77),\
                xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
            
            
        # x-axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)

        # y axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
#         plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower left', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)


    ##########################################################################################
        
        save_name = 'SALTmap_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_zoom_{4}_inclim_{5}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation,zoom_limit, inc_limit)

    #     savefig("{0}SALT_map1.pdf".format(directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        summary_filename = '{0}/{1}_summary.txt'.format(out_directory, save_name)
        summary_file = open(summary_filename,'wt')
    
    #     for k, v in sorted(previousNames.iteritems(), key=lambda (k,v): (v,k)):
    #         summary_file.write('{0}. {1}, \n'.format(v,k))

        for num, name in zip(countList_zoom, combinedNameList_zoom):
            summary_file.write('{0}. {1}, \n'.format(num,name))
        
        for num, name in zip(countList_non_zoom, combinedNameList_non_zoom):
            summary_file.write('{0}. {1}, \n'.format(num,name))
    
    
        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        combinedNameList_corotate = []
        combinedNameList_antirotate = []
        combinedNameList_maybe = []
        for name, ra, dec, c in zip(combinedNameList_zoom, RA_targetList_zoom, Dec_targetList_zoom, markerColorList_zoom):
            if c == color_yes:
                combinedNameList_corotate.append(name)
            if c == color_maybe:
                combinedNameList_maybe.append(name)
            if c == color_no:
                combinedNameList_antirotate.append(name)
                
                
        stats_file.write('Total co-rotate = {0} \n'.format(len(combinedNameList_corotate)))
        stats_file.write('Total anti-rotate = {0} \n'.format(len(combinedNameList_antirotate)))
        stats_file.write('Total uncertain = {0} \n'.format(len(combinedNameList_maybe)))
        
        stats_file.write('\n')
        stats_file.write('Co-rotating systems: \n')
        for i in combinedNameList_corotate:
            stats_file.write('{0}\n'.format(i))
        
        stats_file.write('\n')
        stats_file.write('Anti-rotating systems: \n')
        for i in combinedNameList_antirotate:
            stats_file.write('{0}\n'.format(i))
    
        stats_file.write('\n')
        stats_file.write('Uncertain systems: \n')
        for i in combinedNameList_maybe:
            stats_file.write('{0}\n'.format(i))
            
            
        stats_file.close()
        summary_file.close()



##########################################################################################
##########################################################################################
# NFW model plot but zoomed into 1 R_vir radius only
#
#
##########################################################################################
##########################################################################################

    if plot_NFW_zoom_in:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
    #     fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)

    ##########################################################################################
        # plot circles
        def xy(r,phi):
          return r*np.cos(phi), r*np.sin(phi)

        phis=np.arange(0,2*np.pi,0.01)
    
        r = 0.5
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        r = 1.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-1,1],c='black',ls='-',lw=0.6)
        ax.plot([-1,1],[0,0],c='black',ls='-',lw=0.6)
    
#         ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

    
        # plot the rest
        largest_N = max(NList_zoom)
        smallest_N = min(NList_zoom)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for N in NList_zoom:
            newSize = ((float(N) - smallest_N)/(largest_N - smallest_N)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # make different style markers for different colors
        for i in arange(len(markerColorList_NFW_zoom)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList_NFW_zoom[i] == color_maybe:
                marker = 'o'
            if markerColorList_NFW_zoom[i] == color_no:
                marker = 'X'
#                 marker_lw = 1.5
                marker_lw = 0.5
            if markerColorList_NFW_zoom[i] == color_yes:
                marker = 'D'
        

            ax.scatter(RA_targetList_zoom[i], Dec_targetList_zoom[i], color=markerColorList_NFW_zoom[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    
        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList_zoom)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)

#             annotate(countList_zoom[i],xy=(RA_targetList_zoom[i], Dec_targetList_zoom[i]),\
#             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)

            tag = nameDict[combinedNameList_zoom[i]]        
            if not tag in previousNames:
                annotate(tag,xy=(RA_targetList_zoom[i], Dec_targetList_zoom[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1
                
    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non_zoom)):
            ax.plot(RA_targetList_non_zoom[i], Dec_targetList_non_zoom[i], color=markerColorList_non_zoom[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
#             annotate(countList_non_zoom[i],xy=(RA_targetList_non_zoom[i], Dec_targetList_non_zoom[i]),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
            
            tag = nameDict[combinedNameList_non_zoom[i]]
        
            if not tag in previousNames:
                annotate(tag,xy=(RA_targetList_non_zoom[i], Dec_targetList_non_zoom[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1


    ##########################################################################################

        xlabel(r'$\rm x ~[R_{vir}]$')
        ylabel(r'$\rm y ~[R_{vir}]$')
    
        ax.set_xlim(-zoom_limit, zoom_limit)
        ax.set_ylim(-zoom_limit, zoom_limit)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(zoom_limit-0.06, 0.02),\
        xytext=(0.0,0.0),textcoords='offset points',size=9)
        
        annotate(r'$\rm NFW~model~(zoom)$', xy=(0.96, 0.91),\
        xytext=(0.0,0.0), textcoords='offset points', size=12, weight='bold')
        
        if only_close_velocities:
            annotate(r'$\rm | \Delta v | \leq v_{\rm rot}$', xy=(0.96, 0.84),\
            xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
        
            if Lstar_note:
                annotate(r'$\rm L \leq 0.5 L^{\**}$', xy=(0.96, 0.77),\
                xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
            
        # x-axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)

        # y axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
#         plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower left', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)


    ##########################################################################################
        
        save_name = 'SALTmap_NFW_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_zoom_{4}_inclim_{5}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation,zoom_limit, inc_limit)

    #     savefig("{0}SALT_map1.pdf".format(directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        summary_filename = '{0}/{1}_summary.txt'.format(out_directory, save_name)
        summary_file = open(summary_filename,'wt')
    
    #     for k, v in sorted(previousNames.iteritems(), key=lambda (k,v): (v,k)):
    #         summary_file.write('{0}. {1}, \n'.format(v,k))

        for num, name in zip(countList_zoom, combinedNameList_zoom):
            summary_file.write('{0}. {1}, \n'.format(num,name))
        
        for num, name in zip(countList_non_zoom, combinedNameList_non_zoom):
            summary_file.write('{0}. {1}, \n'.format(num,name))
    
    
        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        combinedNameList_corotate = []
        combinedNameList_antirotate = []
        combinedNameList_maybe = []
        for name, ra, dec, c in zip(combinedNameList_zoom, RA_targetList_zoom, Dec_targetList_zoom, markerColorList_NFW_zoom):
            if c == color_yes:
                combinedNameList_corotate.append(name)
            if c == color_maybe:
                combinedNameList_maybe.append(name)
            if c == color_no:
                combinedNameList_antirotate.append(name)
                
                
        stats_file.write('Total co-rotate = {0} \n'.format(len(combinedNameList_corotate)))
        stats_file.write('Total anti-rotate = {0} \n'.format(len(combinedNameList_antirotate)))
        stats_file.write('Total uncertain = {0} \n'.format(len(combinedNameList_maybe)))
        
        stats_file.write('\n')
        stats_file.write('Co-rotating systems: \n')
        for i in combinedNameList_corotate:
            stats_file.write('{0}\n'.format(i))
        
        stats_file.write('\n')
        stats_file.write('Anti-rotating systems: \n')
        for i in combinedNameList_antirotate:
            stats_file.write('{0}\n'.format(i))
    
        stats_file.write('\n')
        stats_file.write('Uncertain systems: \n')
        for i in combinedNameList_maybe:
            stats_file.write('{0}\n'.format(i))
            
        stats_file.close()
        summary_file.close()


##########################################################################################
##########################################################################################
# Plot folded NFW model plot showing co-rotation vs major, minor axis
#
#
##########################################################################################
##########################################################################################

    if plot_azmap:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
    #     fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)

    ##########################################################################################
        # plot circles
#         def xy(r,phi):
#           return r*np.cos(phi), r*np.sin(phi)
# 
#         phis=np.arange(0,2*np.pi,0.01)
#     
#         r = 1.0
#         ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
#     
#         r = 2.0
#         ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
#     
#         r = 3.0
#         ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
#         
#         ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
#         ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
    ##########################################################################################

        plot_RA_list = np.abs(RA_targetList2)
        plot_Dec_list = np.abs(Dec_targetList2)
        plot_RA_list_non = np.abs(RA_targetList_non)
        plot_Dec_list_non = np.abs(Dec_targetList_non)
    
        # plot the rest
        largest_N = max(NList2)
        smallest_N = min(NList2)
        maxSize = 500
        minSize = 30    
    
        newSizeList = []
        for N in NList2:
            newSize = ((float(N) - smallest_N)/(largest_N - smallest_N)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # make different style markers for different colors
        for i in arange(len(markerColorList2)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList2[i] == color_maybe:
                marker = 'o'
#                 marker = 'x'
            if markerColorList2[i] == color_no:
#                 marker = 'X'
                marker = 'X'
                marker_lw = 0.5
            if markerColorList2[i] == color_yes:
                marker = 'D'
                
#             if abs(dvList[i]) <= 15:
#                 marker = 'o'
#                 color = 'grey'
# 
#                 ax.scatter(plot_RA_list[i], plot_Dec_list[i], color='grey', \
#                 s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
#             else:
#                 ax.scatter(plot_RA_list[i], plot_Dec_list[i], color=markerColorList2[i], \
#                 s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    

            ax.scatter(plot_RA_list[i], plot_Dec_list[i], color=markerColorList2[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    
        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)
    #         print 'combinedNameList2[i]: ',combinedNameList2[i], 'dec = ',Dec_targetList2[i]
    #         print 'yTagOffset: ',yTagOffset
                
#             annotate(countList[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
#             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
                
            tag = nameDict[combinedNameList2[i]]
            print('tag = {0}, name = {1}, fit_v = {2}'.format(tag, combinedNameList2[i], dvList2[i]))
            print()
            if not tag in previousNames:
                annotate(tag,xy=(plot_RA_list[i], plot_Dec_list[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1


    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non)):
            ax.plot(plot_RA_list_non[i], plot_Dec_list_non[i], color=markerColorList_non[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
#             annotate(countList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        
            tag = nameDict[combinedNameList_non[i]]
        
            if not tag in previousNames:
                annotate(tag,xy=(plot_RA_list_non[i], plot_Dec_list_non[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
            
                previousNames[tag] = 1


    ##########################################################################################

        xlabel(r'$\rm Major~Axis ~[R_{vir}]$')
        ylabel(r'$\rm Minor~Axis ~[R_{vir}]$')
    

        ax.set_xlim(-0.1, 2.5)
        ax.set_ylim(-0.1, 2.5)
#         ax.invert_xaxis()
    
#         annotate(r'$\rm Approaching~ Side$',xy=(2.97, 0.06),\
#         xytext=(0.0,0.0),textcoords='offset points',size=9)
# 
#         annotate(r'$\rm NFW~Model$', xy=(2.90, 2.72),\
#         xytext=(0.0,0.0), textcoords='offset points', size=12, weight='bold')
#         
#         if only_close_velocities:
#             annotate(r'$\rm | \Delta v | \leq v_{\rm rot}$', xy=(2.90, 2.52),\
#             xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
#         
#             if Lstar_note:
#                 annotate(r'$\rm L \leq 0.5 L^{\**}$', xy=(2.90, 2.32),\
#                 xytext=(0.0,0.0), textcoords='offset points', size=11, weight='bold')
            
        # x-axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)

        # y axis
    #     majorLocator   = MultipleLocator(0.5)
    #     majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
    #     minorLocator   = MultipleLocator(0.25)
    #     ax.yaxis.set_major_locator(majorLocator)
    #     ax.yaxis.set_major_formatter(majorFormatter)
    #     ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
        plt.legend(handles=[corotate, antirotate, nondetection],loc='upper right', 
                                borderpad=0.8, fontsize=legend_font, fancybox=True)


    ##########################################################################################
        
        save_name = 'SALT_azmap_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_inclim_{4}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation, inc_limit)

    #     savefig("{0}SALT_map1.pdf".format(out_directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
# end   


    
if __name__ == '__main__':
    main()
    