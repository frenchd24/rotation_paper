#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id: SALT_paper_Lstar_fraction4.py, v4.0 10/26/18

Plot co-rotating fraction as a funtion of Lstar

v2: Not sure when this happened... probably mid-June for Alabama WHIM conference (5/10/18)

v4: Hopefully final paper (submitted) version (10/26/18)

v5: Updated for referee report (10/17/2019)

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
from mpl_toolkits.axes_grid1 import make_axes_locatable


from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

from scipy.optimize import curve_fit
import matplotlib.ticker as ticker

from matplotlib import rc
fontScale = 16
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


def equal_hist_edges(x, nbin):
    '''
    Calculates approximately nbin equal size bins of input data x
    
    Inputs:
    -------
    nbin    : number of bins requested (int)
    x       : input data to be binned (1D array)
    
    
    Returns:
    --------
    array of bin edges
    
    '''
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin+1),
                    np.arange(npt),
                    np.sort(x))




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
#     out_directory = '/Users/frenchd/Research/test/SALT_maps_redo/'
#     out_directory = '/Users/frenchd/Research/rotation_paper_data/plots_alt/'
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
    
    # which plot to make?
    plot_b_vs_dv_apparent = False
    plot_b_vs_dv_NFW = False
    plot_b_vs_dv_NFW_Lstar = False
    plot_b_vs_dv_hists_NFW = False
    plot_inc_dv_hists_NFW = False

    plot_b_comparison_NFW = False
    plot_b_comparison = False
    plot_b_hist_NFW_matching = False
    plot_b_multiple_galaxy_NFW = False
    plot_N_multiple_galaxy_NFW = False
    plot_b_azimuth_NFW = False
    plot_dv_inc_NFW = False
    plot_dv_az_NFW = False
    plot_N_inc_NFW = False
    plot_N_az_NFW = True
    plot_az_inc_NFW = False
    
    # --- use Voigt profile fits instead of AOD results?
    use_fits = True
    
    # --- include e_vhel_measured and Lya v_fit error in the apparent velocity ranges?
    use_apparent_errors = True
    
    # --- plot gridlines?
    plotGrid = False
    
    # don't include absorbers with EW above this
    EW_cut = 10000.
    
    # --- remove systems of lower inclination than this:
    inc_limit = 0
    
    # --- remove systems with more nearby companions than this:
    neighbor_limit = 1000
    
    # --- Lstar at which to split results (where applicable)
    L_limit = 0.5
    
    # include tags to include
#     include_tags = ['yes','maybe' 'no/maybe']
    include_tags = ['yes','maybe']
#     include_tags = ['yes']

    # plot properties:
    # markers
    marker_apparent = 's'
    marker_steidel = 'd'
    marker_NFW = 'D'
    
    # line styles
    ls_apparent = '--'
    ls_steidel = '-.'
    ls_NFW = '-'
    
    # --- print detailed information for each galaxy?
    verbose = False
    
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
    bList = []
    e_bList = []
    NList = []
    e_NList = []
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
    LstarList = []
    e_LstarList = []
    impactvirList = []
    impactList = []
    
    apparent_answers = []
    steidel_answers = []
    NFW_answers = []

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
    bList = []
    e_bList = []
    RA_targetList_non = []
    Dec_targetList_non = []
    incList_non = []
    paList_non = []
    azList_non = []
    RvirList_non = []
    markerColorList_non = []
    VhelList_non = []
    impactvirList_non = []
    impactList_non = []


    nameDict = {}
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
        b = eval(t['b'])
        e_b = eval(t['e_b'])
        
        # --- from Voigt profile fits
        fit_v   = float(t['fit_v'])
        e_fit_v = float(t['e_fit_v'])
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
        
        apparent_corotation = t['apparent_corotation']
        n_companions  = int(t['n_companions'])
        
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

        
        # calculate impact parameter and shit
#         impact = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist)
    
        # RA component of impact parameter - by setting the Dec to be the same for both
        impact_RA = calculateImpactParameter(RA_galaxy, Dec_galaxy, RA_target, Dec_galaxy, dist)
    
        # Dec component of impact parameter - by setting the RA to be the same for both
        impact_Dec = calculateImpactParameter(RA_galaxy, Dec_galaxy, RA_galaxy, Dec_target, dist)
        
        # --- flip around some shit        
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
            if verbose:
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
            
        if n_companions > neighbor_limit:
            print('n_companions > limit : ',n_companions)
            add_to_list = False
        
        separation = 500
        if min_separation and add_to_list:
            # correlate with environment
            agnSeparation = False
            minVcorr = False
            minSize = False
            correlation = correlateSingle.correlateTarget(name, min_separation, agnSeparation, minVcorr, minSize, slow=False, searchAll=True)
            galaxyInfo = correlation[name]
            
            print('galaxyInfo: ',galaxyInfo)
                
            for row in galaxyInfo:
                vhel, galaxyRow = row
                separation = galaxyRow['impactParameter (kpc)']
                galaxyVel = galaxyRow['radialVelocity (km/s)']
                
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
                    e_NList.append(e_fit_N)
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
                    
                    steidel_answers.append(steidel_answer)
                    NFW_answers.append(NFW_answer)
                    apparent_answers.append(apparent_answer)
                    
                    dvList.append(dv)
                    dv_up_list.append(dv_up)
                    dv_down_list.append(dv_down)
                    e_dv_list.append(e_dv)
                    
                    if use_fits:
                        bList.append(fit_b)
                        e_bList.append(e_fit_b)
                    else:
                        bList.append(b)
                        e_bList.append(e_b)

                    LstarList.append(Lstar)
                    e_LstarList.append(e_Lstar)
                    impactList.append(impact)
                    impactvirList.append(impact_rvir)
                
                
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
                    e_NList.append(e_fit_N)
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
                    
                    steidel_answers.append(steidel_answer)
                    NFW_answers.append(NFW_answer)
                    apparent_answers.append(apparent_answer)

                    dvList.append(dv)
                    dv_up_list.append(dv_up)
                    dv_down_list.append(dv_down)
                    e_dv_list.append(e_dv)

                    if use_fits:
                        bList.append(fit_b)
                        e_bList.append(e_fit_b)
                    else:
                        bList.append(b)
                        e_bList.append(e_b)

                    LstarList.append(Lstar)
                    e_LstarList.append(e_Lstar)
                    impactList.append(impact)
                    impactvirList.append(impact_rvir)
                    
                
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
                impactList_non.append(impact)
                impactvirList_non.append(impact_rvir)

        else:
            print('{0} excluded. Lstar={1}, separation={2}'.format(name, Lstar, separation))
            print()

        
    model_file.close()
    data_file.close()
        
##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
# Plot b values for co-rotators vs anti-rotators - apparent
#
#
##########################################################################################
##########################################################################################

    if plot_b_comparison:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(311)

#         bins = arange(0,100,10)
        bins = arange(10,90,10)

        alpha_no = 0.55
        alpha_yes = 0.65        

        corotate_b = []
        antirotate_b = []
        
        corotate_b_close = []
        antirotate_b_close = []
        
        corotate_b_far = []
        antirotate_b_far = []
        
        Lstar_high = []
        Lstar_low = []
        
        for ans, b, ra, dec, L in zip(apparent_answers, bList, RA_targetList, Dec_targetList, LstarList):
            if ans:
                corotate_b.append(b)
                
                if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                    corotate_b_close.append(b)
                    
                else:
                    corotate_b_far.append(b)
            
            else:
                antirotate_b.append(b)
                
                if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                    antirotate_b_close.append(b)
                    
                else:
                    antirotate_b_far.append(b)
                    
            if L <= L_limit:
                Lstar_low.append(b)
            else:
                Lstar_high.append(b)
                    
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(antirotate_b,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_no,
            alpha=alpha_no,
            edgecolor='black',
            hatch='/',
            label=r'$\rm Anti-rotators$')

        hist(corotate_b,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_yes,
            alpha=alpha_yes,
            edgecolor='black',
            label=r'$\rm Co-rotators$')
        
        
        ylim(0, 14)
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm b ~ [km~s^{-1}]$')
        ylabel(r'$\rm Number$')


        ax = fig.add_subplot(312)
                    
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(antirotate_b_close,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_no,
            alpha=alpha_no, 
            edgecolor='black',
            hatch='/',
            label=r'$\rm Anti-rotators~(\rho \leq {0} ~R_{{vir}})$'.format(int(zoom_limit)))
        
        hist(corotate_b_close,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_yes,
            alpha=alpha_yes,
            edgecolor='black',
            label=r'$\rm Co-rotators~(\rho \leq {0} ~R_{{vir}})$'.format(int(zoom_limit)))
        
        ylim(0, 4)
        legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
        xlabel(r'$\rm b ~ [km~s^{-1}]$')
        ylabel(r'$\rm Number$')
        
        ax = fig.add_subplot(313)
                    
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(Lstar_high,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_no,
            alpha=alpha_no,
            hatch='/',
            edgecolor='black',
            label=r'$\rm L^{{\**}} > {0})$'.format(L_limit))
            
        hist(Lstar_low,
            bins=bins, 
            histtype='bar', 
            lw=1.5,
            color=color_yes,
            alpha=alpha_yes,
            edgecolor='black',
            label=r'$\rm L^{{\**}} \leq {0})$'.format(L_limit))
        
        
#         ylim(0, 3)
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm b ~ [km~s^{-1}]$')
        ylabel(r'$\rm Number$')
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')


#         import matplotlib.patches as mpatches
#         import matplotlib.lines as mlines
# 
#         corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
# 
#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
#                               
#         antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
#         plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALT_bhist_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')
            
        
        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        avg_b_corotate = mean(corotate_b)
        avg_b_antirotate = mean(antirotate_b)
        
        med_b_corotate = median(corotate_b)
        med_b_antirotate = median(antirotate_b)
        
        std_b_corotate = std(corotate_b)
        std_b_antirotate = std(antirotate_b)

        # now the zoomed in ones
        avg_b_corotate_close = mean(corotate_b_close)
        avg_b_antirotate_close = mean(antirotate_b_close)
        
        med_b_corotate_close = median(corotate_b_close)
        med_b_antirotate_close = median(antirotate_b_close)
        
        std_b_corotate_close = std(corotate_b_close)
        std_b_antirotate_close = std(antirotate_b_close)


        # now the zoomed out ones
        avg_b_corotate_far = mean(corotate_b_far)
        avg_b_antirotate_far = mean(antirotate_b_far)
        
        med_b_corotate_far = median(corotate_b_far)
        med_b_antirotate_far = median(antirotate_b_far)
        
        std_b_corotate_far = std(corotate_b_far)
        std_b_antirotate_far = std(antirotate_b_far)
        
        
        # now the Lstars
        avg_b_Lstar_high = mean(Lstar_high)
        avg_b_Lstar_low = mean(Lstar_low)
        
        med_b_Lstar_high = median(Lstar_high)
        med_b_Lstar_low = median(Lstar_low)
        
        std_b_Lstar_high = std(Lstar_high)
        std_b_Lstar_low = std(Lstar_low)
                
        stats_file.write('Average b co-rotate = {0} \n'.format(avg_b_corotate))
        stats_file.write('Average b anti-rotate = {0} \n'.format(avg_b_antirotate))
        stats_file.write('Median b co-rotate = {0} \n'.format(med_b_corotate))
        stats_file.write('Median b anti-rotate = {0} \n'.format(med_b_antirotate))
        stats_file.write('Std b co-rotate = {0} \n'.format(std_b_corotate))
        stats_file.write('Std b anti-rotate = {0} \n'.format(std_b_antirotate))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('Average b co-rotate zoom_in = {0} \n'.format(avg_b_corotate_close))
        stats_file.write('Average b anti-rotate zoom_in = {0} \n'.format(avg_b_antirotate_close))
        stats_file.write('Median b co-rotate zoom_in = {0} \n'.format(med_b_corotate_close))
        stats_file.write('Median b anti-rotate zoom_in = {0} \n'.format(med_b_antirotate_close))
        stats_file.write('Std b co-rotate zoom_in = {0} \n'.format(std_b_corotate_close))
        stats_file.write('Std b anti-rotate zoom_in = {0} \n'.format(std_b_antirotate_close))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('Average b co-rotate zoom_out = {0} \n'.format(avg_b_corotate_far))
        stats_file.write('Average b anti-rotate zoom_out = {0} \n'.format(avg_b_antirotate_far))
        stats_file.write('Median b co-rotate zoom_out = {0} \n'.format(med_b_corotate_far))
        stats_file.write('Median b anti-rotate zoom_out = {0} \n'.format(med_b_antirotate_far))
        stats_file.write('Std b co-rotate zoom_out = {0} \n'.format(std_b_corotate_far))
        stats_file.write('Std b anti-rotate zoom_out = {0} \n'.format(std_b_antirotate_far))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('Average b Lstar high = {0} \n'.format(avg_b_Lstar_high))
        stats_file.write('Average b Lstar low = {0} \n'.format(avg_b_Lstar_low))
        stats_file.write('Median b Lstar high = {0} \n'.format(med_b_Lstar_high))
        stats_file.write('Median b Lstar low = {0} \n'.format(med_b_Lstar_low))
        stats_file.write('Std b Lstar high = {0} \n'.format(std_b_Lstar_high))
        stats_file.write('Std b Lstar low = {0} \n'.format(std_b_Lstar_low))
        stats_file.write('\n')
        stats_file.write('\n')
        
        
        stats_file.close()    


##########################################################################################
##########################################################################################
# Plot b values for co-rotators vs anti-rotators - NFW
#
#
##########################################################################################
##########################################################################################

    if plot_b_comparison_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(311)

#         bins = arange(0,100,10)
        bins = arange(10,90,10)

        alpha_no = 0.55
        alpha_yes = 0.65        

        corotate_b = []
        antirotate_b = []
        
        corotate_b_close = []
        antirotate_b_close = []
        
        corotate_b_far = []
        antirotate_b_far = []
        
        Lstar_high = []
        Lstar_low = []
        
        for ans, b, ra, dec, L in zip(NFW_answers, bList, RA_targetList, Dec_targetList, LstarList):
            if ans:
                corotate_b.append(b)
                
                if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                    corotate_b_close.append(b)
                    
                else:
                    corotate_b_far.append(b)
            
            else:
                antirotate_b.append(b)
                
                if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                    antirotate_b_close.append(b)
                    
                else:
                    antirotate_b_far.append(b)
                    
            if L <= L_limit:
                Lstar_low.append(b)
            else:
                Lstar_high.append(b)
                    
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(antirotate_b,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_no,
            alpha=alpha_no,
            edgecolor='black',
            hatch='/',
            label=r'$\rm Anti-rotators$')

        hist(corotate_b,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_yes,
            alpha=alpha_yes,
            edgecolor='black',
            label=r'$\rm Co-rotators$')
        
        
        ylim(0, 10)
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm b ~ [km~s^{-1}]$')
        ylabel(r'$\rm Number$')


        ax = fig.add_subplot(312)
                    
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(antirotate_b_close,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_no,
            alpha=alpha_no, 
            edgecolor='black',
            hatch='/',
            label=r'$\rm Anti-rotators~(\rho \leq {0} ~R_{{vir}})$'.format(int(zoom_limit)))
        
        hist(corotate_b_close,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_yes,
            alpha=alpha_yes,
            edgecolor='black',
            label=r'$\rm Co-rotators~(\rho \leq {0} ~R_{{vir}})$'.format(int(zoom_limit)))
        
        ylim(0, 4)
        legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
        xlabel(r'$\rm b ~ [km~s^{-1}]$')
        ylabel(r'$\rm Number$')
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        ax = fig.add_subplot(313)
                    
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(Lstar_high,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_no,
            alpha=alpha_no,
            hatch='/',
            edgecolor='black',
            label=r'$\rm L^{{\**}} > {0})$'.format(L_limit))
            
        hist(Lstar_low,
            bins=bins, 
            histtype='bar', 
            lw=1.5,
            color=color_yes,
            alpha=alpha_yes,
            edgecolor='black',
            label=r'$\rm L^{{\**}} \leq {0})$'.format(L_limit))
        
        
#         ylim(0, 3)
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm b ~ [km~s^{-1}]$')
        ylabel(r'$\rm Number$')
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

#         import matplotlib.patches as mpatches
#         import matplotlib.lines as mlines
# 
#         corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
# 
#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
#                               
#         antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
#         plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALT_bhist_NFW_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')
            
        
        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        avg_b_corotate = mean(corotate_b)
        avg_b_antirotate = mean(antirotate_b)
        
        med_b_corotate = median(corotate_b)
        med_b_antirotate = median(antirotate_b)
        
        std_b_corotate = std(corotate_b)
        std_b_antirotate = std(antirotate_b)

        # now the zoomed in ones
        avg_b_corotate_close = mean(corotate_b_close)
        avg_b_antirotate_close = mean(antirotate_b_close)
        
        med_b_corotate_close = median(corotate_b_close)
        med_b_antirotate_close = median(antirotate_b_close)
        
        std_b_corotate_close = std(corotate_b_close)
        std_b_antirotate_close = std(antirotate_b_close)


        # now the zoomed out ones
        avg_b_corotate_far = mean(corotate_b_far)
        avg_b_antirotate_far = mean(antirotate_b_far)
        
        med_b_corotate_far = median(corotate_b_far)
        med_b_antirotate_far = median(antirotate_b_far)
        
        std_b_corotate_far = std(corotate_b_far)
        std_b_antirotate_far = std(antirotate_b_far)
        
        
        # now the Lstars
        avg_b_Lstar_high = mean(Lstar_high)
        avg_b_Lstar_low = mean(Lstar_low)
        
        med_b_Lstar_high = median(Lstar_high)
        med_b_Lstar_low = median(Lstar_low)
        
        std_b_Lstar_high = std(Lstar_high)
        std_b_Lstar_low = std(Lstar_low)
                
        stats_file.write('Average b co-rotate = {0} \n'.format(avg_b_corotate))
        stats_file.write('Average b anti-rotate = {0} \n'.format(avg_b_antirotate))
        stats_file.write('Median b co-rotate = {0} \n'.format(med_b_corotate))
        stats_file.write('Median b anti-rotate = {0} \n'.format(med_b_antirotate))
        stats_file.write('Std b co-rotate = {0} \n'.format(std_b_corotate))
        stats_file.write('Std b anti-rotate = {0} \n'.format(std_b_antirotate))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('Average b co-rotate zoom_in = {0} \n'.format(avg_b_corotate_close))
        stats_file.write('Average b anti-rotate zoom_in = {0} \n'.format(avg_b_antirotate_close))
        stats_file.write('Median b co-rotate zoom_in = {0} \n'.format(med_b_corotate_close))
        stats_file.write('Median b anti-rotate zoom_in = {0} \n'.format(med_b_antirotate_close))
        stats_file.write('Std b co-rotate zoom_in = {0} \n'.format(std_b_corotate_close))
        stats_file.write('Std b anti-rotate zoom_in = {0} \n'.format(std_b_antirotate_close))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('Average b co-rotate zoom_out = {0} \n'.format(avg_b_corotate_far))
        stats_file.write('Average b anti-rotate zoom_out = {0} \n'.format(avg_b_antirotate_far))
        stats_file.write('Median b co-rotate zoom_out = {0} \n'.format(med_b_corotate_far))
        stats_file.write('Median b anti-rotate zoom_out = {0} \n'.format(med_b_antirotate_far))
        stats_file.write('Std b co-rotate zoom_out = {0} \n'.format(std_b_corotate_far))
        stats_file.write('Std b anti-rotate zoom_out = {0} \n'.format(std_b_antirotate_far))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('Average b Lstar high = {0} \n'.format(avg_b_Lstar_high))
        stats_file.write('Average b Lstar low = {0} \n'.format(avg_b_Lstar_low))
        stats_file.write('Median b Lstar high = {0} \n'.format(med_b_Lstar_high))
        stats_file.write('Median b Lstar low = {0} \n'.format(med_b_Lstar_low))
        stats_file.write('Std b Lstar high = {0} \n'.format(std_b_Lstar_high))
        stats_file.write('Std b Lstar low = {0} \n'.format(std_b_Lstar_low))
        stats_file.write('\n')
        stats_file.write('\n')
        
        
        # Statistical distribution tests
        #
        # corotate_b vs anti-rotate_b
        ks = stats.ks_2samp(corotate_b, antirotate_b)
        ad = stats.anderson_ksamp([corotate_b, antirotate_b])
        rs, p_val = stats.ranksums(corotate_b, antirotate_b)

        stats_file.write('------------------------------------\n')
        stats_file.write('Statistical tests:\n')
        stats_file.write('Co- vs anti-rotating:\n')
        stats_file.write('KS : corotate_b vs antirotate_b: {0}\n'.format(ks))
        stats_file.write('AD : corotate_b vs antirotate_b: {0}\n'.format(ad))
        stats_file.write('Ranksum : corotate_b vs antirotate_b: {0}, pval={1}\n'.format(rs,p_val))
        stats_file.write('\n')
        stats_file.write('\n')
        
        # Lstar_high vs Lstar_low
        ks = stats.ks_2samp(Lstar_high, Lstar_low)
        ad = stats.anderson_ksamp([Lstar_high, Lstar_low])
        rs, p_val = stats.ranksums(Lstar_high, Lstar_low)

        stats_file.write('-----------\n')
        stats_file.write('Lstar_high vs Lstar_low :\n')
        stats_file.write('KS : Lstar_high vs Lstar_low: {0}\n'.format(ks))
        stats_file.write('AD : Lstar_high vs Lstar_low: {0}\n'.format(ad))
        stats_file.write('Ranksum : Lstar_high vs Lstar_low: {0}, pval={1}\n'.format(rs, p_val))
        stats_file.write('\n')
        stats_file.write('\n')
        
        stats_file.close()    



##########################################################################################
##########################################################################################
# Plot b values vs dv - NFW model
#
#
##########################################################################################
##########################################################################################

    if plot_b_vs_dv_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 1.
        alpha_yes = 1.
        
        corotate_b = []
        antirotate_b = []
        e_corotate_b= []
        e_antirotate_b= []

        corotate_dv = []
        e_corotate_dv = []

        antirotate_dv = []
        e_antirotate_dv = []
        
        for nfw, b, e_b, L, dv, e_dv in zip(NFW_answers, bList,  e_bList, LstarList, dvList, e_dv_list):
            print('b +/- e_b, dv +/- e_dv : ',b,'+/-',e_b, ', ',dv,'+/-', e_dv)
            if nfw:
                corotate_b.append(b)
                e_corotate_b.append(e_b)
                corotate_dv.append(dv)
                e_corotate_dv.append(e_dv)
            
            else:
                antirotate_b.append(b)
                e_antirotate_b.append(e_b)
                antirotate_dv.append(dv)
                e_antirotate_dv.append(e_dv)

        errorbar(antirotate_dv,
                antirotate_b,
                xerr=e_antirotate_dv,
                yerr=e_antirotate_b,
                lw=0,
                elinewidth = 0.8,
                marker='X',
                markeredgewidth=0.8,
                markeredgecolor='black',
                ms = 9,
                color=color_no,
                alpha=alpha_no,
                label=r'$\rm Anti-rotators$')

        errorbar(corotate_dv,
                corotate_b,
                xerr=e_corotate_dv,
                yerr=e_corotate_b,
                lw=0,
                elinewidth = 0.8,
                marker='D',
                markeredgewidth=0.8,
                markeredgecolor='black',
                ms = 9,
                color=color_yes,
                alpha=alpha_yes,
                label=r'$\rm Co-rotators$')

        
        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        ylim(0, 140)
        xlim(-250, 250)
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
        ylabel(r'$\rm b ~ [km~s^{-1}]$')
        xlabel(r'$\rm \Delta v ~[km~s^{-1}]$')
        


#         import matplotlib.patches as mpatches
#         import matplotlib.lines as mlines
# 
#         corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
# 
#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
#                               
#         antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
#         plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALT_b_vs_dv_NFW_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
# Plot b values vs dv - NFW model. Include histograms on both edges
#
#
##########################################################################################
##########################################################################################

    if plot_b_vs_dv_hists_NFW:
        # initial figure
        fig = plt.figure(figsize=(9,9))
        subplots_adjust(hspace=0.500)

#         ax = fig.add_subplot(111)


        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
#         spacing = 0.005
        spacing = 0.015

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom + height + spacing, width, 0.2]
        rect_histy = [left + width + spacing, bottom, 0.2, height]

        ax_scatter = plt.axes(rect_scatter)
        ax_scatter.tick_params(direction='in', top=True, right=True)
        
        ax_histx = plt.axes(rect_histx)
        ax_histx.tick_params(direction='in', labelbottom=False)
        
        ax_histy = plt.axes(rect_histy)
        ax_histy.tick_params(direction='in', labelleft=False)

        alpha_no = 0.7
        alpha_yes = 0.7
        
        corotate_b = []
        antirotate_b = []
        e_corotate_b= []
        e_antirotate_b= []

        corotate_dv = []
        e_corotate_dv = []

        antirotate_dv = []
        e_antirotate_dv = []
        
        for nfw, b, e_b, L, dv, e_dv in zip(NFW_answers, bList,  e_bList, LstarList, dvList, e_dv_list):
            print('b +/- e_b, dv +/- e_dv : ',b,'+/-',e_b, ', ',dv,'+/-', e_dv)
            if nfw:
                corotate_b.append(b)
                e_corotate_b.append(e_b)
                corotate_dv.append(dv)
                e_corotate_dv.append(e_dv)
            
            else:
                antirotate_b.append(b)
                e_antirotate_b.append(e_b)
                antirotate_dv.append(dv)
                e_antirotate_dv.append(e_dv)


        ax_scatter.errorbar(corotate_dv,
                        corotate_b,
                        xerr=e_corotate_dv,
                        yerr=e_corotate_b,
                        lw=0,
                        elinewidth = 0.8,
                        marker='D',
                        markeredgewidth=0.8,
                        markeredgecolor='black',
                        ms = 9,
                        color=color_yes,
                        alpha=alpha_yes,
                        label=r'$\rm Co-rotators$')
                        
        ax_scatter.errorbar(antirotate_dv,
                        antirotate_b,
                        xerr=e_antirotate_dv,
                        yerr=e_antirotate_b,
                        lw=0,
                        elinewidth = 0.8,
                        marker='X',
                        markeredgewidth=0.8,
                        markeredgecolor='black',
                        ms = 9,
                        color=color_no,
                        alpha=alpha_no,
                        label=r'$\rm Anti-rotators$')

        nbins = len(ax_scatter.get_xticklabels())
#         ax_scatter.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
#         print('x nbins: ',nbins)
#         print('ax_scatter.get_xticklabels(): ',ax_scatter.get_xticklabels())
        
        # x-axis
        x_majorLocator   = MultipleLocator(50)
        x_majorFormatter = FormatStrFormatter(r'$\rm %d$')
        x_minorLocator   = MultipleLocator(25)
#         ax_scatter.xaxis.set_major_locator(MaxNLocator(nbins=14, prune='upper'))
        ax_scatter.xaxis.set_major_locator(x_majorLocator)

        ax_scatter.xaxis.set_major_formatter(x_majorFormatter)
        ax_scatter.xaxis.set_minor_locator(x_minorLocator)
        
        # y-axis
        y_majorLocator   = MultipleLocator(10)
        y_majorFormatter = FormatStrFormatter(r'$\rm %d$')
        y_minorLocator   = MultipleLocator(5)
        ax_scatter.yaxis.set_major_locator(y_majorLocator)
        ax_scatter.yaxis.set_major_formatter(y_majorFormatter)
        ax_scatter.yaxis.set_minor_locator(y_minorLocator)
        
        ax_scatter.set_ylim(0, 140)
        ax_scatter.set_xlim(-250, 250)
        
        ax_scatter.yaxis.set_ticks_position('both')
        ax_scatter.xaxis.set_ticks_position('both')
        
        ax_scatter.legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
        ax_scatter.set_ylabel(r'$\rm b ~ [km~s^{-1}]$')
        ax_scatter.set_xlabel(r'$\rm \Delta v ~[km~s^{-1}]$')
        
#         nbins = len(ax_scatter.get_yticklabels())
#         print('y nbins: ',nbins)
#         ax_scatter.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
        
        
        bins = np.arange(-250, 300, 50)
        # --- plot them
        ax_histx.hist(corotate_dv,
                    bins=bins,
                    histtype='bar',
                    lw=1.5,
                    color=color_yes,
                    alpha=alpha_yes,
                    edgecolor='black')
            
        ax_histx.hist(antirotate_dv,
                    bins=bins,
                    histtype='bar',
                    lw=1.5,
                    color=color_no,
                    alpha=alpha_no,
                    edgecolor='black',
                    hatch='/')

        ax_histx.set_xlim(ax_scatter.get_xlim())
                    
#         ax.set_yticks(y_count)
#         ax.set_yticklabels(y_vals)
                    
        # --- y axis formatting
#         majorLocator   = MultipleLocator(4)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(2)
#         ax_histx.yaxis.set_major_locator(majorLocator)
#         ax_histx.yaxis.set_major_formatter(majorFormatter)
#         ax_histx.yaxis.set_minor_locator(minorLocator)
        
        bins = np.arange(0, 130, 10)
        # --- plot them
        ax_histy.hist(corotate_b,
                    bins=bins,
                    histtype='bar',
                    lw=1.5,
                    color=color_yes,
                    alpha=alpha_yes,
                    edgecolor='black', orientation='horizontal')
            
        ax_histy.hist(antirotate_b,
                    bins=bins,
                    histtype='bar',
                    lw=1.5,
                    color=color_no,
                    alpha=alpha_no,
                    edgecolor='black',
                    hatch='/', orientation='horizontal')
                    
        ax_histy.set_xlim(0, 10)
        ax_histy.set_ylim(ax_scatter.get_ylim())

#         ax.set_yticks(y_count)
#         ax.set_yticklabels(y_vals)
                    
        # --- y axis formatting
#         majorLocator   = MultipleLocator(10)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(5)
#         ax_histy.yaxis.set_major_locator(majorLocator)
#         ax_histy.yaxis.set_major_formatter(majorFormatter)
#         ax_histy.yaxis.set_minor_locator(minorLocator)


##########################################################################################
        
        save_name = 'SALT_b_vs_dv_hists_NFW_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')



##########################################################################################
##########################################################################################
# Plot histograms of dv at 3 inclinations bins - NFW model
#
#
##########################################################################################
##########################################################################################

    if plot_inc_dv_hists_NFW:
        # initial figure
        fig = plt.figure(figsize=(9,9))
        subplots_adjust(hspace=0.500)

        ax1 = fig.add_subplot(311)

        alpha_no = 0.7
        alpha_yes = 0.7
        
        corotate_dv_inc1 = []
        antirotate_dv_inc1 = []
        corotate_dv_inc2 = []
        antirotate_dv_inc2 = []
        corotate_dv_inc3 = []
        antirotate_dv_inc3 = []

        antirotate_dv = []
        e_antirotate_dv = []
        
        for nfw, b, e_b, L, dv, e_dv, i in zip(NFW_answers, bList,  e_bList, LstarList, dvList, e_dv_list, incList):
            print('b +/- e_b, dv +/- e_dv : ',b,'+/-',e_b, ', ',dv,'+/-', e_dv)
            if nfw:
                if i < 65:
                    corotate_dv_inc1.append(dv)
                elif i >=65 and i <80:
                    corotate_dv_inc2.append(dv)
                else:
                    corotate_dv_inc3.append(dv)
            
            else:
                if i < 65:
                    antirotate_dv_inc1.append(dv)
                elif i >=65 and i <80:
                    antirotate_dv_inc2.append(dv)
                else:
                    antirotate_dv_inc3.append(dv)


        bins = np.arange(-250, 275,25)
        # --- plot them
        ax1.hist(corotate_dv_inc1,
                bins=bins,
                histtype='bar',
                lw=1.5,
                color=color_yes,
                alpha=alpha_yes,
                edgecolor='black')
                
        ax1.hist(antirotate_dv_inc1,
                bins=bins,
                histtype='bar',
                lw=1.5,
                hatch='/',
                color=color_no,
                alpha=alpha_no,
                edgecolor='black')
        
        ax1.set_title(r'$\rm i < 65$')
                
        #####
        ax2 = fig.add_subplot(312)

        ax2.hist(corotate_dv_inc2,
                bins=bins,
                histtype='bar',
                lw=1.5,
                color=color_yes,
                alpha=alpha_yes,
                edgecolor='black')
                
        ax2.hist(antirotate_dv_inc2,
                bins=bins,
                histtype='bar',
                lw=1.5,
                hatch='/',
                color=color_no,
                alpha=alpha_no,
                edgecolor='black')
                
        ax2.set_title(r'$\rm 65 \leq i < 80$')


        #####
        ax3 = fig.add_subplot(313)

        ax3.hist(corotate_dv_inc3,
                bins=bins,
                histtype='bar',
                lw=1.5,
                color=color_yes,
                alpha=alpha_yes,
                edgecolor='black')
                
        ax3.hist(antirotate_dv_inc3,
                bins=bins,
                histtype='bar',
                lw=1.5,
                hatch='/',
                color=color_no,
                alpha=alpha_no,
                edgecolor='black')
                
        ax3.set_title(r'$\rm 80 \leq i \leq 90$')


        # x-axis
        x_majorLocator   = MultipleLocator(50)
        x_majorFormatter = FormatStrFormatter(r'$\rm %d$')
        x_minorLocator   = MultipleLocator(25)
#         ax_scatter.xaxis.set_major_locator(MaxNLocator(nbins=14, prune='upper'))
        ax3.xaxis.set_major_locator(x_majorLocator)
        ax3.xaxis.set_major_formatter(x_majorFormatter)
        ax3.xaxis.set_minor_locator(x_minorLocator)
        
        # y-axis
#         y_majorLocator   = MultipleLocator(2)
#         y_majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         y_minorLocator   = MultipleLocator(1)
#         ax3.yaxis.set_major_locator(y_majorLocator)
#         ax3.yaxis.set_major_formatter(y_majorFormatter)
#         ax3.yaxis.set_minor_locator(y_minorLocator)
        
#         ax.set_ylim(0, 140)
#         ax.set_xlim(-250, 250)
        
        ax3.yaxis.set_ticks_position('both')
        ax3.xaxis.set_ticks_position('both')
        
        ax3.legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
        ax3.set_xlabel(r'$\rm \Delta v ~[km~s^{-1}]$')
        

##########################################################################################
        
        save_name = 'SALT_inc_dv_hists_NFW_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
# Plot b values vs dv 
#
#
##########################################################################################
##########################################################################################

    if plot_b_vs_dv_apparent:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 1.
        alpha_yes = 1.
        
        corotate_b = []
        antirotate_b = []
        e_corotate_b= []
        e_antirotate_b= []

        corotate_dv = []
        e_corotate_dv = []

        antirotate_dv = []
        e_antirotate_dv = []
        
        for ans, b, e_b, L, dv, e_dv in zip(apparent_answers, bList,  e_bList, LstarList, dvList, e_dv_list):
            print('b +/- e_b, dv +/- e_dv : ',b,'+/-',e_b, ', ',dv,'+/-',e_dv)
            if ans:
                corotate_b.append(b)
                e_corotate_b.append(e_b)
                corotate_dv.append(dv)
                e_corotate_dv.append(e_dv)
            
            else:
                antirotate_b.append(b)
                e_antirotate_b.append(e_b)
                antirotate_dv.append(dv)
                e_antirotate_dv.append(e_dv)

        errorbar(antirotate_dv,
                antirotate_b,
                xerr=e_antirotate_dv,
                yerr=e_antirotate_b,
                lw=0,
                elinewidth = 0.8,
                marker='X',
                markeredgewidth=0.8,
                markeredgecolor='black',
                ms = 9,
                color=color_no,
                alpha=alpha_no,
                label=r'$\rm Anti-rotators$')

        errorbar(corotate_dv,
                corotate_b,
                xerr=e_corotate_dv,
                yerr=e_corotate_b,
                lw=0,
                elinewidth = 0.8,
                marker='D',
                markeredgewidth=0.8,
                markeredgecolor='black',
                ms = 9,
                color=color_yes,
                alpha=alpha_yes,
                label=r'$\rm Co-rotators$')

        
        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        ylim(0, 140)
        xlim(-250, 250)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
        ylabel(r'$\rm b ~ [km~s^{-1}]$')
        xlabel(r'$\rm \Delta v ~[km~s^{-1}]$')
        

#         import matplotlib.patches as mpatches
#         import matplotlib.lines as mlines
# 
#         corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
# 
#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
#                               
#         antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
#         plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_b_vs_dv_apparent_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

##########################################################################################
##########################################################################################
# Plot b values vs dv - NFW model, separate based on Lstar, not co-rotation
#
#
##########################################################################################
##########################################################################################

    if plot_b_vs_dv_NFW_Lstar:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 1.
        alpha_yes = 1.
        
        corotate_b = []
        antirotate_b = []
        e_corotate_b= []
        e_antirotate_b= []

        corotate_dv = []
        e_corotate_dv = []

        antirotate_dv = []
        e_antirotate_dv = []
        
        for nfw, b, e_b, L, dv, e_dv in zip(NFW_answers, bList,  e_bList, LstarList, dvList, e_dv_list):
            print('b +/- e_b, dv +/- e_dv : ',b,'+/-',e_b, ', ',dv,'+/-',e_dv)
            if L <= L_limit:
                corotate_b.append(b)
                e_corotate_b.append(e_b)
                corotate_dv.append(dv)
                e_corotate_dv.append(e_dv)
            
            else:
                antirotate_b.append(b)
                e_antirotate_b.append(e_b)
                antirotate_dv.append(dv)
                e_antirotate_dv.append(e_dv)


        errorbar(antirotate_dv,
                antirotate_b,
                xerr=e_antirotate_dv,
                yerr=e_antirotate_b,
                lw=0,
                elinewidth = 0.8,
                marker='X',
                markeredgewidth=0.8,
                markeredgecolor='black',
                ms = 9,
                color=color_no,
                alpha=alpha_no,
                label=r'$\rm L^{{\**}} > {0}$'.format(L_limit))

        errorbar(corotate_dv,
                corotate_b,
                xerr=e_corotate_dv,
                yerr=e_corotate_b,
                lw=0,
                elinewidth = 0.8,
                marker='D',
                markeredgewidth=0.8,
                markeredgecolor='black',
                ms = 9,
                color=color_yes,
                alpha=alpha_yes,
                label=r'$\rm L^{{\**}} \leq {0}$'.format(L_limit))

        
        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        ylim(0, 140)
        xlim(-250, 250)
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        ylabel(r'$\rm b ~ [km~s^{-1}]$')
        xlabel(r'$\rm \Delta v ~[km~s^{-1}]$')
        
        legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)

#         import matplotlib.patches as mpatches
#         import matplotlib.lines as mlines
# 
#         corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
# 
#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
#                               
#         antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
#         plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALT_b_vs_dv_NFW_Lstar_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')




##########################################################################################
##########################################################################################
# Plot histograms of b values for  co-rotation and anti- when there are multiple absorbers
# associated with one galaxy
#
#
##########################################################################################
##########################################################################################

    if plot_b_hist_NFW_matching:
        # initial figure
        fig = plt.figure(figsize=(8,4))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 0.6
        alpha_yes = 0.6
        
        corotate_b = []
        antirotate_b = []
        e_corotate_b= []
        e_antirotate_b= []

        corotate_dv = []
        e_corotate_dv = []

        antirotate_dv = []
        e_antirotate_dv = []
        
        bins = np.arange(0,120,10)
        
        for name, ans, b, e_b in zip(nameList, NFW_answers, bList, e_bList):
            name_indices = np.where(np.array(nameList) == name)[0]
            
            if len(name_indices) >1:
#                 print('Matching name: ',name, len(name_indices))
                print('name = {0}, b= {1} +/- {2}, corotating = {3}'.format(name, b, e_b, ans))

                if ans:
                    corotate_b.append(b)
                    e_corotate_b.append(e_b)
                else:
                    antirotate_b.append(b)
                    e_antirotate_b.append(e_b)

        hist(antirotate_b,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_no,
            alpha=alpha_no,
            edgecolor='black',
            hatch='/',
            label=r'$\rm Anti-rotators$')

        hist(corotate_b,
            bins=bins,
            histtype='bar',
            lw=1.5,
            color=color_yes,
            alpha=alpha_yes,
            edgecolor='black',
            label=r'$\rm Co-rotators$')
            
        

        # perform the K-S and AD tests 
        
        mean_co = mean(corotate_b)
        median_co = median(corotate_b)
        mean_anti = mean(antirotate_b)
        median_anti = median(antirotate_b)
        std_co = std(corotate_b)
        std_anti = std(antirotate_b)
        
        ks_co_anti = stats.ks_2samp(corotate_b, antirotate_b)
        ad_co_anti = stats.anderson_ksamp([corotate_b,antirotate_b])

        print('KS for corotate_b vs antirotate_b : ',ks_co_anti)
        print('AD for corotate_b vs antirotate_b : ',ad_co_anti)
        print()
        
        z_statrb, p_valrb = stats.ranksums(corotate_b, antirotate_b)
        print('ranksum for corotate_b vs antirotate_b : ',z_statrb, p_valrb)

        
#         ylim(0, 10)
#         legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
#         xlabel(r'$b \rm ~[km~s^{-1}]$')
#         ylabel(r'$\rm Number$')
                    
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
#         
#         ylim(0, 140)
#         xlim(-250, 250)
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        ylabel(r'$\rm Number$')
        xlabel(r'$ b \rm~[km~s^{-1}]$')
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')


#         import matplotlib.patches as mpatches
#         import matplotlib.lines as mlines
# 
#         corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
# 
#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
#                               
#         antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
#         plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALT_b_hist_NFW_matching_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        stats_file.write('------------------ \n')
        stats_file.write('\n')
        stats_file.write('Mean b co-rotators = {0} \n'.format(mean_co))
        stats_file.write('Median b co-rotators = {0} \n'.format(median_co))
        stats_file.write('Std b co-rotators = {0} \n'.format(std_co))
        stats_file.write('\n')
        stats_file.write('Mean b anti-rotators = {0} \n'.format(mean_anti))
        stats_file.write('Median b anti-rotators = {0} \n'.format(median_anti))
        stats_file.write('Std b anti-rotators = {0} \n'.format(std_anti))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('---- \n')
        stats_file.write('\n')
        stats_file.write('KS test co vs anti-rotators = {0} \n'.format(ks_co_anti))
        stats_file.write('AD test co vs anti-rotators = {0} \n'.format(ad_co_anti))
        stats_file.write('Ranksum test co vs anti-rotators: z_stat = {0}, p_val = {1} \n'.format(z_statrb, p_valrb))
        stats_file.write('\n')
        stats_file.close()

##########################################################################################
##########################################################################################
# Plot histograms of b values for  co-rotation and anti- when there are multiple absorbers
# associated with one galaxy
#
#
##########################################################################################
##########################################################################################

    if plot_b_multiple_galaxy_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 0.6
        alpha_yes = 0.6
        
        corotate_b = []
        antirotate_b = []
        e_corotate_b= []
        e_antirotate_b= []

        corotate_dv = []
        e_corotate_dv = []

        antirotate_dv = []
        e_antirotate_dv = []
        
        bins = np.arange(0,110,10)
        names_dict = {}
        all_corotate_bs = []
        all_corotate_e_bs = []
        all_antirotate_bs = []
        all_antirotate_e_bs = []
        
        all_corotate_Ns = []
        all_corotate_e_Ns = []
        all_antirotate_Ns = []
        all_antirotate_e_Ns = []
        
        print("NameList: ",nameList)
        print()
        
        # --- first grab the systems that have multiple components
        for name, ans, b, e_b, N, e_N in zip(nameList, NFW_answers, bList, e_bList, NList, e_NList):
            # --- indices of this galaxy name in the data list
            name_indices = np.where(np.array(nameList) == name)[0]
            
            if len(name_indices) >1:
                print('name = {0}, b= {1} +/- {2}, corotating = {3}'.format(name, b, e_b, ans))
                
                if ans:
                    all_corotate_bs.append(b)
                    all_corotate_e_bs.append(e_b)
                    all_corotate_Ns.append(N)
                    all_corotate_e_Ns.append(e_N)
                else:
                    all_antirotate_bs.append(b)
                    all_antirotate_e_bs.append(e_b)
                    all_antirotate_Ns.append(N)
                    all_antirotate_e_Ns.append(e_N)

                # --- for each unique galaxy name add all the measurements to a dict
                if name in names_dict:
                    vals = names_dict[name]
                    vals['b'].append(b)
                    vals['e_b'].append(e_b)
                    vals['N'].append(N)
                    vals['e_N'].append(e_N)
                    vals['ans'].append(ans)
                    
                else:
                    names_dict[name] = {'b':[b], 'e_b':[e_b], 'N':[N], 'e_N':[e_N], 'ans':[ans]}
                    
        
        
        y_vals = list(names_dict.keys())
        y_count = np.arange(len(y_vals))
        
        print('y_vals: ',y_vals)
        print('y_count: ',y_count)
        print()

        count = 0
        # --- loop  through each galaxy in the dictionary
        for n in y_vals:
            
            vals = names_dict[n]
            bs = vals['b']
            e_bs = vals['e_b']
            Ns = vals['N']
            e_Ns = vals['e_N']

            ans = vals['ans']
            colors = []
            symbols = []
            x_vals = []
            
            #  --- loop through each measurement for this galaxy
            for b, e, N, e_N, a in zip(bs, e_bs, Ns, e_Ns, ans):
                if a:
                    symbol = 'D'
                    color = color_yes
                    alpha = alpha_yes+0.2
                else:
                    symbol = 'X'
                    color = color_no
                    alpha = alpha_no+0.2
                
                colors.append(color)
                symbols.append(symbol)
                x_vals.append(count)

                ax.errorbar(b,
                        count,
                        xerr=e,
                        lw=1,
                        elinewidth = 0.8,
                        marker=symbol,
                        markeredgewidth=0.8,
                        markeredgecolor='black',
                        ms = 12,
                        color=color,
                        alpha=alpha)
            
            count+=1


        ax.errorbar(y_count,
                y_count,
                lw=0,
                elinewidth = 0,
                marker='.',
                markeredgewidth=0,
                markeredgecolor='black',
                ms = 0,
                color='white',
                alpha=0)

        # --- x-axis formatting
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
#         ax.set_aspect(5.0)

        
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        ylabel(r'$\rm Galaxy$')
        xlabel(r'$\rm \emph{b} ~[km~s^{-1}]$')

        #  --- make a histogram on top showing the distribution
        divider = make_axes_locatable(ax)
        axHistx = divider.append_axes("top", 1.2, pad=0.09, sharex=ax)
        
        axHistx.xaxis.set_tick_params(labelbottom=False)
        
        # --- plot them
        axHistx.hist(all_corotate_bs,
                    bins=bins,
                    histtype='bar',
                    lw=1.5,
                    color=color_yes,
                    alpha=alpha_yes,
                    edgecolor='black',
                    label=r'$\rm Co-rotators$')
            
        axHistx.hist(all_antirotate_bs,
                    bins=bins,
                    histtype='bar',
                    lw=1.5,
                    color=color_no,
                    alpha=alpha_no,
                    edgecolor='black',
                    hatch='/',
                    label=r'$\rm Anti-rotators$')
            
                    
        ax.set_yticks(y_count)
        ax.set_yticklabels(y_vals)
                    

        # --- y axis formatting
        majorLocator   = MultipleLocator(4)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2)
        axHistx.yaxis.set_major_locator(majorLocator)
        axHistx.yaxis.set_major_formatter(majorFormatter)
        axHistx.yaxis.set_minor_locator(minorLocator)
        
        xlim(0, 100)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')


        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
# 
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,alpha=alpha_yes+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
# 
#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
#                               
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0, alpha=alpha_no+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
        ax.legend(handles=[corotate, antirotate],loc='upper right', 
                                borderpad=0.7, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_b_multipl_NFW_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')
        
        ks_co_anti = stats.ks_2samp(all_corotate_bs, all_antirotate_bs)
        ad_co_anti = stats.anderson_ksamp([all_corotate_bs, all_antirotate_bs])
        z_statrb, p_valrb = stats.ranksums(all_corotate_bs, all_antirotate_bs)
        
        stats_file.write('------------------ \n')
        stats_file.write('\n')
        stats_file.write('Mean b co-rotators = {0} \n'.format(mean(all_corotate_bs)))
        stats_file.write('Median b co-rotators = {0} \n'.format(median(all_corotate_bs)))
        stats_file.write('Std b co-rotators = {0} \n'.format(std(all_corotate_bs)))
        stats_file.write('\n')
        stats_file.write('Mean b anti-rotators = {0} \n'.format(mean(all_antirotate_bs)))
        stats_file.write('Median b anti-rotators = {0} \n'.format(median(all_antirotate_bs)))
        stats_file.write('Std b anti-rotators = {0} \n'.format(std(all_antirotate_bs)))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('---- \n')
        stats_file.write('\n')
        stats_file.write('KS test co vs anti-rotators = {0} \n'.format(ks_co_anti))
        stats_file.write('AD test co vs anti-rotators = {0} \n'.format(ad_co_anti))
        stats_file.write('Ranksum test co vs anti-rotators: z_stat = {0}, p_val = {1} \n'.format(z_statrb, p_valrb))
        stats_file.write('\n')
        stats_file.close()


##########################################################################################
##########################################################################################
# Plot histograms of column densities for  co-rotation and anti- when there are multiple 
# absorbers associated with one galaxy
#
#
##########################################################################################
##########################################################################################

    if plot_N_multiple_galaxy_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 0.6
        alpha_yes = 0.6
        
        corotate_b = []
        antirotate_b = []
        e_corotate_b= []
        e_antirotate_b= []

        corotate_dv = []
        e_corotate_dv = []

        antirotate_dv = []
        e_antirotate_dv = []
        
        bins = np.arange(12,22,1)
        names_dict = {}
        all_corotate_bs = []
        all_corotate_e_bs = []
        all_antirotate_bs = []
        all_antirotate_e_bs = []
        
        all_corotate_Ns = []
        all_corotate_e_Ns = []
        all_antirotate_Ns = []
        all_antirotate_e_Ns = []
        
        print("NameList: ",nameList)
        print()
        
        # --- first grab the systems that have multiple components
        for name, ans, b, e_b, N, e_N in zip(nameList, NFW_answers, bList, e_bList, NList, e_NList):
            # --- indices of this galaxy name in the data list
            name_indices = np.where(np.array(nameList) == name)[0]
            
            if len(name_indices) >1:
                print('name = {0}, b= {1} +/- {2}, corotating = {3}'.format(name, b, e_b, ans))
                
                if ans:
                    all_corotate_bs.append(b)
                    all_corotate_e_bs.append(e_b)
                    all_corotate_Ns.append(N)
                    all_corotate_e_Ns.append(e_N)
                else:
                    all_antirotate_bs.append(b)
                    all_antirotate_e_bs.append(e_b)
                    all_antirotate_Ns.append(N)
                    all_antirotate_e_Ns.append(e_N)

                # --- for each unique galaxy name add all the measurements to a dict
                if name in names_dict:
                    vals = names_dict[name]
                    vals['b'].append(b)
                    vals['e_b'].append(e_b)
                    vals['N'].append(N)
                    vals['e_N'].append(e_N)
                    vals['ans'].append(ans)
                    
                else:
                    names_dict[name] = {'b':[b], 'e_b':[e_b], 'N':[N], 'e_N':[e_N], 'ans':[ans]}
                    
        
        
        y_vals = list(names_dict.keys())
        y_count = np.arange(len(y_vals))
        
        print('y_vals: ',y_vals)
        print('y_count: ',y_count)
        print()

        count = 0
        # --- loop  through each galaxy in the dictionary
        for n in y_vals:
            
            vals = names_dict[n]
            bs = vals['b']
            e_bs = vals['e_b']
            Ns = vals['N']
            e_Ns = vals['e_N']

            ans = vals['ans']
            colors = []
            symbols = []
            x_vals = []
            
            #  --- loop through each measurement for this galaxy
            for b, e, N, e_N, a in zip(bs, e_bs, Ns, e_Ns, ans):
                if a:
                    symbol = 'D'
                    color = color_yes
                    alpha = alpha_yes+0.2
                else:
                    symbol = 'X'
                    color = color_no
                    alpha = alpha_no+0.2
                
                colors.append(color)
                symbols.append(symbol)
                x_vals.append(count)

                ax.errorbar(N,
                        count,
                        xerr=e_N,
                        lw=1,
                        elinewidth = 0.8,
                        marker=symbol,
                        markeredgewidth=0.8,
                        markeredgecolor='black',
                        ms = 12,
                        color=color,
                        alpha=alpha)
            
            count+=1


        ax.errorbar(y_count,
                y_count,
                lw=0,
                elinewidth = 0,
                marker='.',
                markeredgewidth=0,
                markeredgecolor='black',
                ms = 0,
                color='white',
                alpha=0)

        majorLocator   = MultipleLocator(1)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(0.5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)

        
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        ylabel(r'$\rm Galaxy$')
        xlabel(r'$\rm log\emph{N}$')

        #  --- make a histogram on top showing the distribution
        divider = make_axes_locatable(ax)
        axHistx = divider.append_axes("top", 1.2, pad=0.09, sharex=ax)
        
        axHistx.xaxis.set_tick_params(labelbottom=False)
        
        # --- plot them

        axHistx.hist(all_corotate_Ns,
                    bins=bins,
                    histtype='bar',
                    lw=1.5,
                    color=color_yes,
                    alpha=alpha_yes,
                    edgecolor='black',
                    label=r'$\rm Co-rotators$')
            
        axHistx.hist(all_antirotate_Ns,
                    bins=bins,
                    histtype='bar',
                    lw=1.5,
                    color=color_no,
                    alpha=alpha_no,
                    edgecolor='black',
                    hatch='/',
                    label=r'$\rm Anti-rotators$')
            
                    
        ax.set_yticks(y_count)
        ax.set_yticklabels(y_vals)
        xlim(12, 21)

        # --- y axis formatting
        majorLocator   = MultipleLocator(4)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2)
        axHistx.yaxis.set_major_locator(majorLocator)
        axHistx.yaxis.set_major_formatter(majorFormatter)
        axHistx.yaxis.set_minor_locator(minorLocator)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')


        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
# 
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,alpha=alpha_yes+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
# 
#         maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
#                                   markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
#                               
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0, alpha=alpha_no+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
        ax.legend(handles=[corotate, antirotate],loc='upper right', 
                                borderpad=0.7, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_N_multipl_NFW_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')
        
        ks_co_anti = stats.ks_2samp(all_corotate_bs, all_antirotate_bs)
        ad_co_anti = stats.anderson_ksamp([all_corotate_bs, all_antirotate_bs])
        z_statrb, p_valrb = stats.ranksums(all_corotate_bs, all_antirotate_bs)
        
        stats_file.write('------------------ \n')
        stats_file.write('\n')
        stats_file.write('Mean b co-rotators = {0} \n'.format(mean(all_corotate_bs)))
        stats_file.write('Median b co-rotators = {0} \n'.format(median(all_corotate_bs)))
        stats_file.write('Std b co-rotators = {0} \n'.format(std(all_corotate_bs)))
        stats_file.write('\n')
        stats_file.write('Mean b anti-rotators = {0} \n'.format(mean(all_antirotate_bs)))
        stats_file.write('Median b anti-rotators = {0} \n'.format(median(all_antirotate_bs)))
        stats_file.write('Std b anti-rotators = {0} \n'.format(std(all_antirotate_bs)))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('---- \n')
        stats_file.write('\n')
        stats_file.write('KS test co vs anti-rotators = {0} \n'.format(ks_co_anti))
        stats_file.write('AD test co vs anti-rotators = {0} \n'.format(ad_co_anti))
        stats_file.write('Ranksum test co vs anti-rotators: z_stat = {0}, p_val = {1} \n'.format(z_statrb, p_valrb))
        stats_file.write('\n')
        stats_file.close()


##########################################################################################
##########################################################################################
# Plot b values as a function of azimuth angle for anti- and co-rotatorss
#
#
##########################################################################################
##########################################################################################

    if plot_b_azimuth_NFW:
        # initial figure
        fig = plt.figure(figsize=(7.7,5.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)
        
        include_bhist = False

        alpha_no = 0.6
        alpha_yes = 0.6

        min_inc = 0
        
        corotate_b = []
        antirotate_b = []
        corotate_e_b = []
        antirotate_e_b = []

        corotate_az = []
        corotate_e_az = []
        antirotate_az = []
        antirotate_e_az = []
        
        names_dict = {}
        all_corotate_bs = []
        all_corotate_e_bs = []
        all_antirotate_bs = []
        all_antirotate_e_bs = []
        
        # --- separate values based on NFW_answers
        for name, ans, b, e_b, az, inc in zip(nameList, NFW_answers, bList, e_bList, azList, incList):
            if inc >= min_inc:
                if ans:
                    corotate_b.append(b)
                    corotate_e_b.append(e_b)
                    corotate_az.append(az)
                else:
                    antirotate_b.append(b)
                    antirotate_e_b.append(e_b)
                    antirotate_az.append(az)


        ax.errorbar(corotate_az,
                    corotate_b,
                    yerr=corotate_e_b,
                    lw=0,
                    elinewidth = 0.8,
                    marker='D',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_yes,
                    alpha=alpha_yes+0.2)

        ax.errorbar(antirotate_az,
                    antirotate_b,
                    yerr=antirotate_e_b,
                    lw=0,
                    elinewidth = 0.8,
                    marker='X',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_no,
                    alpha=alpha_no+0.2)

        # --- x-axis formatting
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
                
        ylim(0, 130)
        xlim(0, 91)
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        ylabel(r'$\rm \emph{b} ~[km~s^{-1}]$')
        xlabel(r'$\rm Azimuth~[deg]$')
        
        #  --- make a histogram on top showing the distribution
        if include_bhist:
            divider = make_axes_locatable(ax)
            axHistx = divider.append_axes("top", 1.2, pad=0.15, sharex=ax)
        
            axHistx.xaxis.set_tick_params(labelbottom=False)
        
            # --- plot them
        
            bins = np.arange(0,100,15)
            axHistx.hist(corotate_az,
                        bins=bins,
                        histtype='bar',
                        lw=1.5,
    #                     density=True,
                        color=color_yes,
                        alpha=alpha_yes,
                        edgecolor='black',
                        label=r'$\rm Co-rotators$')
            
            axHistx.hist(antirotate_az,
                        bins=bins,
                        histtype='bar',
                        lw=1.5,
    #                     density=True,
                        color=color_no,
                        alpha=alpha_no,
                        edgecolor='black',
                        hatch='/',
                        label=r'$\rm Anti-rotators$')
                                        

            # --- y axis formatting
            majorLocator   = MultipleLocator(2)
            majorFormatter = FormatStrFormatter(r'$\rm %d$')
            minorLocator   = MultipleLocator(1)
            axHistx.yaxis.set_major_locator(majorLocator)
            axHistx.yaxis.set_major_formatter(majorFormatter)
            axHistx.yaxis.set_minor_locator(minorLocator)

#         ax.yaxis.set_ticks_position('both')
#         ax.xaxis.set_ticks_position('both')

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,alpha=alpha_yes+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation ~(NFW)$')

        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0, alpha=alpha_no+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation ~(NFW)$')

        ax.legend(handles=[corotate, antirotate],loc='upper left', 
                                borderpad=0.6, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_b_az_NFW_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_mininc_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, min_inc)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
# Plot dv values as a function of inclination for anti- and co-rotatorss
#
#
##########################################################################################
##########################################################################################

    if plot_dv_inc_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 0.6
        alpha_yes = 0.6

        min_inc =  0
        max_az = 90
        
        corotate_dv = []
        antirotate_dv = []
        corotate_e_dv = []
        antirotate_e_dv = []

        corotate_az = []
        corotate_e_az = []
        antirotate_az = []
        antirotate_e_az = []
        
        corotate_inc = []
        corotate_e_inc = []
        antirotate_inc = []
        antirotate_e_inc = []
        
        names_dict = {}
        all_corotate_bs = []
        all_corotate_e_bs = []
        all_antirotate_bs = []
        all_antirotate_e_bs = []
        
        # --- separate values based on NFW_answers
        for name, ans, b, e_b, az, inc, dv, e_dv in zip(nameList, NFW_answers, bList, e_bList, azList, incList, dvList, e_dv_list):
            if inc >= min_inc and az <= max_az:
                if ans:
                    corotate_dv.append(dv)
                    corotate_e_dv.append(e_dv)
                    corotate_inc.append(inc)
                    corotate_az.append(az)
                else:
                    antirotate_dv.append(dv)
                    antirotate_e_dv.append(e_dv)
                    antirotate_inc.append(inc)
                    antirotate_az.append(az)

        ax.errorbar(corotate_inc,
                    corotate_dv,
                    yerr=corotate_e_dv,
                    lw=0,
                    elinewidth = 0.8,
                    marker='D',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_yes,
                    alpha=alpha_yes+0.2)

        ax.errorbar(antirotate_inc,
                    antirotate_dv,
                    yerr=antirotate_e_dv,
                    lw=0,
                    elinewidth = 0.8,
                    marker='X',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_no,
                    alpha=alpha_no+0.2)

        # --- x-axis formatting
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
                
        ylim(-200, 200)
        xlim(0, 91)
        
        if plotGrid:
            grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        ylabel(r'$\rm \Delta v ~[km~s^{-1}]$')
        xlabel(r'$\rm Inclination~[deg]$')

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,alpha=alpha_yes+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')


        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0, alpha=alpha_no+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')


        ax.legend(handles=[corotate, antirotate],loc='upper left', 
                                borderpad=0.6, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_dv_inc_nfw_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_mininc_{4}_maxaz_{5}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, min_inc, max_az)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
# Plot dv values as a function of azimuth for anti- and co-rotators
#
#
##########################################################################################
##########################################################################################

    if plot_dv_az_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 0.6
        alpha_yes = 0.6

        min_inc =  0
        max_az = 90
        
        corotate_dv = []
        antirotate_dv = []
        corotate_e_dv = []
        antirotate_e_dv = []

        corotate_az = []
        corotate_e_az = []
        antirotate_az = []
        antirotate_e_az = []
        
        corotate_inc = []
        corotate_e_inc = []
        antirotate_inc = []
        antirotate_e_inc = []
        
        names_dict = {}
        all_corotate_bs = []
        all_corotate_e_bs = []
        all_antirotate_bs = []
        all_antirotate_e_bs = []
        
        # --- separate values based on NFW_answers
        for name, ans, b, e_b, az, inc, dv, e_dv in zip(nameList, NFW_answers, bList, e_bList, azList, incList, dvList, e_dv_list):
            if inc >= min_inc and az <= max_az:
                if ans:
                    corotate_dv.append(dv)
                    corotate_e_dv.append(e_dv)
                    corotate_inc.append(inc)
                    corotate_az.append(az)
                else:
                    antirotate_dv.append(dv)
                    antirotate_e_dv.append(e_dv)
                    antirotate_inc.append(inc)
                    antirotate_az.append(az)

        # --- plot them
        ax.errorbar(corotate_az,
                    corotate_dv,
                    yerr=corotate_e_dv,
                    lw=0,
                    elinewidth = 0.8,
                    marker='D',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_yes,
                    alpha=alpha_yes+0.2)

        ax.errorbar(antirotate_az,
                    antirotate_dv,
                    yerr=antirotate_e_dv,
                    lw=0,
                    elinewidth = 0.8,
                    marker='X',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_no,
                    alpha=alpha_no+0.2)

        # --- x-axis formatting
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
                
        ylim(-200, 200)
        xlim(0, 91)
        
        if plotGrid:
            grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        ylabel(r'$\rm \Delta v ~[km~s^{-1}]$')
        xlabel(r'$\rm Azimuth~[deg]$')


        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,alpha=alpha_yes+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')


        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0, alpha=alpha_no+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')


        ax.legend(handles=[corotate, antirotate],loc='upper left', 
                                borderpad=0.6, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_dv_az_nfw_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_mininc_{4}_maxaz_{5}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, min_inc, max_az)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
# Plot N values as a function of inclination for anti- and co-rotators
#
#
##########################################################################################
##########################################################################################

    if plot_N_inc_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 0.6
        alpha_yes = 0.6

        min_inc =  0
        max_az = 90
        
        corotate_N = []
        antirotate_N = []
        corotate_e_N = []
        antirotate_e_N = []

        corotate_az = []
        corotate_e_az = []
        antirotate_az = []
        antirotate_e_az = []
        
        corotate_inc = []
        corotate_e_inc = []
        antirotate_inc = []
        antirotate_e_inc = []
        
        names_dict = {}
        all_corotate_bs = []
        all_corotate_e_bs = []
        all_antirotate_bs = []
        all_antirotate_e_bs = []
        
        # --- separate values based on NFW_answers
        for name, ans, b, e_b, az, inc, dv, e_dv, N, e_N in zip(nameList, NFW_answers, bList, e_bList, azList, incList, dvList, e_dv_list, NList, e_NList):
            if inc >= min_inc and az <= max_az:
                if ans:
                    corotate_N.append(N)
                    corotate_e_N.append(e_N)
                    corotate_inc.append(inc)
                else:
                    antirotate_N.append(N)
                    antirotate_e_N.append(e_N)
                    antirotate_inc.append(inc)


        ax.errorbar(corotate_inc,
                    corotate_N,
                    yerr=corotate_e_N,
                    lw=0,
                    elinewidth = 0.8,
                    marker='D',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_yes,
                    alpha=alpha_yes+0.2)

        ax.errorbar(antirotate_inc,
                    antirotate_N,
                    yerr=antirotate_e_N,
                    lw=0,
                    elinewidth = 0.8,
                    marker='X',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_no,
                    alpha=alpha_no+0.2)

        # --- x-axis formatting
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
#         majorLocator   = MultipleLocator(50)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(25)
#         ax.yaxis.set_major_locator(majorLocator)
#         ax.yaxis.set_major_formatter(majorFormatter)
#         ax.yaxis.set_minor_locator(minorLocator)
                
#         ylim(-200, 200)
        xlim(0, 91)
        
        if plotGrid:
            grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        ylabel(r'$\rm log \emph{N}$')
        xlabel(r'$\rm Inclination~[deg]$')

        #  --- make a histogram on top showing the distribution
#         divider = make_axes_locatable(ax)
#         axHistx = divider.append_axes("top", 1.2, pad=0.15, sharex=ax)
#         
#         axHistx.xaxis.set_tick_params(labelbottom=False)
#         
#         # --- plot them
#         bins = np.arange(0,100,15)
#         axHistx.hist(corotate_az,
#                     bins=bins,
#                     histtype='bar',
#                     lw=1.5,
#                     color=color_yes,
#                     alpha=alpha_yes,
#                     edgecolor='black',
#                     label=r'$\rm Co-rotators$')
#             
#         axHistx.hist(antirotate_az,
#                     bins=bins,
#                     histtype='bar',
#                     lw=1.5,
#                     color=color_no,
#                     alpha=alpha_no,
#                     edgecolor='black',
#                     hatch='/',
#                     label=r'$\rm Anti-rotators$')
#                                         
# 
#         # --- y axis formatting
#         majorLocator   = MultipleLocator(2)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(1)
#         axHistx.yaxis.set_major_locator(majorLocator)
#         axHistx.yaxis.set_major_formatter(majorFormatter)
#         axHistx.yaxis.set_minor_locator(minorLocator)
#         

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,alpha=alpha_yes+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')


        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0, alpha=alpha_no+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')


        ax.legend(handles=[corotate, antirotate],loc='upper left', 
                                borderpad=0.6, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_N_inc_nfw_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_mininc_{4}_maxaz_{5}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, min_inc, max_az)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
# Plot N values as a function of azimuth for anti- and co-rotators
#
#
##########################################################################################
##########################################################################################

    if plot_N_az_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 0.6
        alpha_yes = 0.6

        min_inc =  0
        max_az = 90
        
        corotate_N = []
        antirotate_N = []
        corotate_e_N = []
        antirotate_e_N = []

        corotate_az = []
        corotate_e_az = []
        antirotate_az = []
        antirotate_e_az = []
        
        corotate_inc = []
        corotate_e_inc = []
        antirotate_inc = []
        antirotate_e_inc = []
        
        names_dict = {}
        all_corotate_bs = []
        all_corotate_e_bs = []
        all_antirotate_bs = []
        all_antirotate_e_bs = []
        
        # --- separate values based on NFW_answers
        for name, ans, b, e_b, az, inc, dv, e_dv, N, e_N in zip(nameList, NFW_answers, bList, e_bList, azList, incList, dvList, e_dv_list, NList, e_NList):
            if inc >= min_inc and az <= max_az:
                if ans:
                    corotate_N.append(N)
                    corotate_e_N.append(e_N)
                    corotate_az.append(az)
                else:
                    antirotate_N.append(N)
                    antirotate_e_N.append(e_N)
                    antirotate_az.append(az)


        ax.errorbar(corotate_az,
                    corotate_N,
                    yerr=corotate_e_N,
                    lw=0,
                    elinewidth = 0.8,
                    marker='D',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_yes,
                    alpha=alpha_yes+0.2)

        ax.errorbar(antirotate_az,
                    antirotate_N,
                    yerr=antirotate_e_N,
                    lw=0,
                    elinewidth = 0.8,
                    marker='X',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_no,
                    alpha=alpha_no+0.2)

        # --- x-axis formatting
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
#         majorLocator   = MultipleLocator(50)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(25)
#         ax.yaxis.set_major_locator(majorLocator)
#         ax.yaxis.set_major_formatter(majorFormatter)
#         ax.yaxis.set_minor_locator(minorLocator)
                
#         ylim(-200, 200)
        xlim(0, 91)
        
        if plotGrid:
            grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        ylabel(r'$\rm log \emph{N}$')
        xlabel(r'$\rm Azimuth~[deg]$')

        #  --- make a histogram on top showing the distribution
#         divider = make_axes_locatable(ax)
#         axHistx = divider.append_axes("top", 1.2, pad=0.15, sharex=ax)
#         
#         axHistx.xaxis.set_tick_params(labelbottom=False)
#         
#         # --- plot them
#         bins = np.arange(0,100,15)
#         axHistx.hist(corotate_az,
#                     bins=bins,
#                     histtype='bar',
#                     lw=1.5,
#                     color=color_yes,
#                     alpha=alpha_yes,
#                     edgecolor='black',
#                     label=r'$\rm Co-rotators$')
#             
#         axHistx.hist(antirotate_az,
#                     bins=bins,
#                     histtype='bar',
#                     lw=1.5,
#                     color=color_no,
#                     alpha=alpha_no,
#                     edgecolor='black',
#                     hatch='/',
#                     label=r'$\rm Anti-rotators$')
#                                         
# 
#         # --- y axis formatting
#         majorLocator   = MultipleLocator(2)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(1)
#         axHistx.yaxis.set_major_locator(majorLocator)
#         axHistx.yaxis.set_major_formatter(majorFormatter)
#         axHistx.yaxis.set_minor_locator(minorLocator)
#         

        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,alpha=alpha_yes+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')


        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0, alpha=alpha_no+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')


        ax.legend(handles=[corotate, antirotate],loc='upper left', 
                                borderpad=0.6, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_N_az_nfw_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_mininc_{4}_maxaz_{5}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, min_inc, max_az)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
# Plot az values as a function of inclination for anti- and co-rotators
#
#
##########################################################################################
##########################################################################################

    if plot_az_inc_NFW:
        # initial figure
        fig = plt.figure(figsize=(8,7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

        alpha_no = 0.6
        alpha_yes = 0.6

        min_inc =  0
        max_az = 90
        
        corotate_N = []
        antirotate_N = []
        corotate_e_N = []
        antirotate_e_N = []

        corotate_az = []
        corotate_e_az = []
        antirotate_az = []
        antirotate_e_az = []
        
        corotate_inc = []
        corotate_e_inc = []
        antirotate_inc = []
        antirotate_e_inc = []
        
        names_dict = {}
        all_corotate_bs = []
        all_corotate_e_bs = []
        all_antirotate_bs = []
        all_antirotate_e_bs = []
        
        az_dist = []
        
        # --- separate values based on NFW_answers
        for name, ans, b, e_b, az, inc, dv, e_dv, N, e_N in zip(nameList, NFW_answers, bList, e_bList, azList, incList, dvList, e_dv_list, NList, e_NList):
            if inc >= min_inc and az <= max_az:
                if ans:
                    corotate_az.append(az)
                    corotate_inc.append(inc)
                else:
                    antirotate_az.append(az)
                    antirotate_inc.append(inc)

                if inc >= 70:
                    az_dist.append(az)

        ax.errorbar(corotate_inc,
                    corotate_az,
                    lw=0,
                    elinewidth = 0.8,
                    marker='D',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_yes,
                    alpha=alpha_yes+0.2)

        ax.errorbar(antirotate_inc,
                    antirotate_az,
                    lw=0,
                    elinewidth = 0.8,
                    marker='X',
                    markeredgewidth=0.8,
                    markeredgecolor='black',
                    ms = 12,
                    color=color_no,
                    alpha=alpha_no+0.2)
                    
        print('len(az_dist)  = {}'.format(len(az_dist)))
        
        az_lt_40 = np.where(np.array(az_dist) <=40)
        print('az_dist: ',az_dist)
        print('az_lt_40: ',az_lt_40)
        print('len(az_lt_40) = {}'.format(len(az_lt_40)))
        print()

        # --- x-axis formatting
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
#         majorLocator   = MultipleLocator(50)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(25)
#         ax.yaxis.set_major_locator(majorLocator)
#         ax.yaxis.set_major_formatter(majorFormatter)
#         ax.yaxis.set_minor_locator(minorLocator)
                
#         ylim(-200, 200)
        xlim(0, 91)
        
        if plotGrid:
            grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        ylabel(r'$\rm Azimuth~[deg]$')
        xlabel(r'$\rm Inclination~[deg]$')



        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,alpha=alpha_yes+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')


        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0, alpha=alpha_no+0.2,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')


        ax.legend(handles=[corotate, antirotate],loc='upper left', 
                                borderpad=0.6, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_az_inc_nfw_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_mininc_{4}_maxaz_{5}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, min_inc, max_az)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
    # end

if __name__ == '__main__':
    main()
    