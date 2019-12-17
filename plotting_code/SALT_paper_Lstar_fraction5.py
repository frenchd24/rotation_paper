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
#     Lstar_range = [0.0, 100.]
#     Lstar_range = [0.60001, 100.0]
#     Lstar_range = [0.0, 0.6]
    Lstar_range = [0.0, 0.5]
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
    plot_Lstar_hist = False
    plot_corotate_fraction_old = False
    plot_corotate_fraction_dist_old = False
    plot_corotate_fraction_minimum_old = False

    plot_corotate_fraction_new = False
    plot_corotate_impact_new = False
    plot_corotate_impact_rvir_new = False
    plot_corotate_azimuth = True
    plot_corotate_inc = True
    plot_corotate_N = False
    plot_corotate_vmax = False


    # --- use Voigt profile fits instead of AOD results?
    use_fits = True
    
    # --- include e_vhel_measured and Lya v_fit error in the apparent velocity ranges?
    use_apparent_errors = True
    
    # don't include absorbers with EW above this
    EW_cut = 10000.
    
    # --- remove systems of lower inclination than this:
    inc_limit = 0
    
    # --- remove systems with more nearby companions than this:
    neighbor_limit = 1000
    
    # --- plot gridlines?
    plotGrid = False
    
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
    
    #
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
    vmaxList = []
    

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
                    vmaxList.append(rot_vel)

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
                    vmaxList.append(rot_vel)


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
# Plot Lstar values for co-rotators vs anti-rotators
#
#
##########################################################################################
##########################################################################################

    if plot_Lstar_hist:
        # initial figure
        fig = plt.figure(figsize=(8,8))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(211)

#         bins = arange(0,100,10)
#         bins = arange(0.2, 1.2, 0.2)
        bins = arange(0.25, 5.0, 0.25)

        alpha_no = 0.55
        alpha_yes = 0.65

        L_limit = 0.5

        corotate_Lstar = []
        antirotate_Lstar = []
        
        corotate_Lstar_close = []
        antirotate_Lstar_close = []
        
        corotate_Lstar_far = []
        antirotate_Lstar_far = []
        
        Lstar_high = []
        Lstar_low = []
        
        for m, Lstar, ra, dec, L in zip(markerColorList_NFW, LstarList, RA_targetList, Dec_targetList, LstarList):
            if m == color_yes:
                corotate_Lstar.append(Lstar)
                
                if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                    corotate_Lstar_close.append(Lstar)
                    
                else:
                    corotate_Lstar_far.append(Lstar)
            
            if m == color_no:
                antirotate_Lstar.append(Lstar)
                
                if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                    antirotate_Lstar_close.append(Lstar)
                    
                else:
                    antirotate_Lstar_far.append(Lstar)
                    
            if L <= L_limit:
                Lstar_low.append(Lstar)
            else:
                Lstar_high.append(Lstar)
                    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        hist(corotate_Lstar, bins=bins, histtype='bar', lw=1.5, color = color_yes, alpha=alpha_yes, edgecolor='black', label=r'$\rm Co-rotators$')
        hist(antirotate_Lstar, bins=bins, histtype='bar', lw=1.5, hatch='//', color = color_no, edgecolor='black',alpha=alpha_no, label=r'$\rm Anti-rotators$')

        xlim(0, 4)
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm L^{{\**}}$')
        ylabel(r'$\rm Number$')


        ax = fig.add_subplot(212)
                    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        hist(corotate_Lstar_close, bins=bins, histtype='bar', lw=1.5, color=color_yes, alpha=alpha_yes, edgecolor='black', label=r'$\rm Co-rotators~(\rho \leq {0} ~R_{{vir}})$'.format(zoom_limit))
        hist(antirotate_Lstar_close, bins=bins, histtype='bar', lw=1.5, hatch='//', color=color_no, edgecolor='black', alpha=alpha_no, label=r'$\rm Anti-rotators~(\rho \leq {0} ~R_{{vir}})$'.format(zoom_limit))

        xlim(0, 4)
        legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm L^{{\**}}$')
        ylabel(r'$\rm Number$')
        

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
        
        save_name = 'SALT_Lstar_hist_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}'.format(only_close_velocities, include_nondetection, L_limit, min_separation)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')
            
        
        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(out_directory, save_name)
        stats_file = open(stats_filename,'wt')

        avg_Lstar_corotate = mean(corotate_Lstar)
        avg_Lstar_antirotate = mean(antirotate_Lstar)
        
        med_Lstar_corotate = median(corotate_Lstar)
        med_Lstar_antirotate = median(antirotate_Lstar)
        
        std_Lstar_corotate = std(corotate_Lstar)
        std_Lstar_antirotate = std(antirotate_Lstar)

        # now the zoomed in ones
        avg_Lstar_corotate_close = mean(corotate_Lstar_close)
        avg_Lstar_antirotate_close = mean(antirotate_Lstar_close)
        
        med_Lstar_corotate_close = median(corotate_Lstar_close)
        med_Lstar_antirotate_close = median(antirotate_Lstar_close)
        
        std_Lstar_corotate_close = std(corotate_Lstar_close)
        std_Lstar_antirotate_close = std(antirotate_Lstar_close)


        # now the zoomed out ones
        avg_Lstar_corotate_far = mean(corotate_Lstar_far)
        avg_Lstar_antirotate_far = mean(antirotate_Lstar_far)
        
        med_Lstar_corotate_far = median(corotate_Lstar_far)
        med_Lstar_antirotate_far = median(antirotate_Lstar_far)
        
        std_Lstar_corotate_far = std(corotate_Lstar_far)
        std_Lstar_antirotate_far = std(antirotate_Lstar_far)
                
        stats_file.write('Average Lstar co-rotate = {0} \n'.format(avg_Lstar_corotate))
        stats_file.write('Average Lstar anti-rotate = {0} \n'.format(avg_Lstar_antirotate))
        stats_file.write('Median Lstar co-rotate = {0} \n'.format(med_Lstar_corotate))
        stats_file.write('Median Lstar anti-rotate = {0} \n'.format(med_Lstar_antirotate))
        stats_file.write('Std Lstar co-rotate = {0} \n'.format(std_Lstar_corotate))
        stats_file.write('Std Lstar anti-rotate = {0} \n'.format(std_Lstar_antirotate))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('Average Lstar co-rotate zoom_in = {0} \n'.format(avg_Lstar_corotate_close))
        stats_file.write('Average Lstar anti-rotate zoom_in = {0} \n'.format(avg_Lstar_antirotate_close))
        stats_file.write('Median Lstar co-rotate zoom_in = {0} \n'.format(med_Lstar_corotate_close))
        stats_file.write('Median Lstar anti-rotate zoom_in = {0} \n'.format(med_Lstar_antirotate_close))
        stats_file.write('Std Lstar co-rotate zoom_in = {0} \n'.format(std_Lstar_corotate_close))
        stats_file.write('Std Lstar anti-rotate zoom_in = {0} \n'.format(std_Lstar_antirotate_close))
        stats_file.write('\n')
        stats_file.write('\n')
        stats_file.write('Average Lstar co-rotate zoom_out = {0} \n'.format(avg_Lstar_corotate_far))
        stats_file.write('Average Lstar anti-rotate zoom_out = {0} \n'.format(avg_Lstar_antirotate_far))
        stats_file.write('Median Lstar co-rotate zoom_out = {0} \n'.format(med_Lstar_corotate_far))
        stats_file.write('Median Lstar anti-rotate zoom_out = {0} \n'.format(med_Lstar_antirotate_far))
        stats_file.write('Std Lstar co-rotate zoom_out = {0} \n'.format(std_Lstar_corotate_far))
        stats_file.write('Std Lstar anti-rotate zoom_out = {0} \n'.format(std_Lstar_antirotate_far))
        stats_file.write('\n')
        stats_file.write('\n')
        
        stats_file.close()



##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of MINIMUM Lstar
#
#
##########################################################################################
##########################################################################################

    if plot_corotate_fraction_minimum_old:
        # initial figure
        fig = plt.figure(figsize=(7.7,6.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 1.5
        marker_size = 12

        L_limit = 0.5

        # no model
        corotate_02 = 0
        e_Lstar_02 = 0
        total_02 = 0
        
        corotate_04 = 0
        e_Lstar_04 = 0
        total_04 = 0
        
        corotate_06 = 0
        e_Lstar_06 = 0
        total_06 = 0

        corotate_08 = 0
        e_Lstar_08 = 0
        total_08 = 0

        corotate_1 = 0
        e_Lstar_1 = 0
        total_1 = 0
        
        corotate_12 = 0
        e_Lstar_12 = 0
        total_12 = 0
        
        corotate_14 = 0
        e_Lstar_14 = 0
        total_14 = 0

        corotate_16 = 0
        e_Lstar_16 = 0
        total_16 = 0
        
        
        # Steidel model
        corotate_02_steidel = 0
        e_Lstar_02_steidel = 0
        total_02_steidel = 0
        
        corotate_04_steidel = 0
        e_Lstar_04_steidel = 0
        total_04_steidel = 0
        
        corotate_06_steidel = 0
        e_Lstar_06_steidel = 0
        total_06_steidel = 0

        corotate_08_steidel = 0
        e_Lstar_08_steidel = 0
        total_08_steidel = 0

        corotate_1_steidel = 0
        e_Lstar_1_steidel = 0
        total_1_steidel = 0
        
        corotate_12_steidel = 0
        e_Lstar_12_steidel = 0
        total_12_steidel = 0
        
        corotate_14_steidel = 0
        e_Lstar_14_steidel = 0
        total_14_steidel = 0

        corotate_16_steidel = 0
        e_Lstar_16_steidel = 0
        total_16_steidel = 0
        
        # NFW model
        corotate_02_nfw = 0
        e_Lstar_02_nfw = 0
        total_02_nfw = 0
        
        corotate_04_nfw = 0
        e_Lstar_04_nfw = 0
        total_04_nfw = 0
        
        corotate_06_nfw = 0
        e_Lstar_06_nfw = 0
        total_06_nfw = 0

        corotate_08_nfw = 0
        e_Lstar_08_nfw = 0
        total_08_nfw = 0

        corotate_1_nfw = 0
        e_Lstar_1_nfw = 0
        total_1_nfw = 0
        
        corotate_12_nfw = 0
        e_Lstar_12_nfw = 0
        total_12_nfw = 0
        
        corotate_14_nfw = 0
        e_Lstar_14_nfw = 0
        total_14_nfw = 0

        corotate_16_nfw = 0
        e_Lstar_16_nfw = 0
        total_16_nfw = 0
        
        
        for m, m_nfw, m_steidel, e_L, ra, dec, L in zip(markerColorList, markerColorList_NFW, markerColorList_steidel, e_LstarList, RA_targetList, Dec_targetList, LstarList):
            if withinRange(L, [0., 0.25], 0.0):
                if m == color_yes:
                    corotate_02 +=1.
                    e_Lstar_02 += e_L**2
                    
                if m != color_maybe:
                    total_02 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_02_steidel +=1.
                    e_Lstar_02_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_02_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_02_nfw +=1.
                    e_Lstar_02_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_02_nfw +=1.
            
            if withinRange(L, [0., 0.5], 0.0):
                if m == color_yes:
                    corotate_04 +=1.
                    e_Lstar_04 += e_L**2

                if m != color_maybe:
                    total_04 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_04_steidel +=1.
                    e_Lstar_04_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_04_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_04_nfw +=1.
                    e_Lstar_04_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_04_nfw +=1.
                    

            if withinRange(L, [0., 0.75], 0.0):
                if m == color_yes:
                    corotate_06 +=1.
                    e_Lstar_06 += e_L**2

                if m != color_maybe:
                    total_06 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_06_steidel +=1.
                    e_Lstar_06_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_06_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_06_nfw +=1.
                    e_Lstar_06_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_06_nfw +=1.
                    

            if withinRange(L, [0., 1.0], 0.0):
                if m == color_yes:
                    corotate_08 +=1.
                    e_Lstar_08 += e_L**2

                if m != color_maybe:
                    total_08 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_08_steidel +=1.
                    e_Lstar_08_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_08_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_08_nfw +=1.
                    e_Lstar_08_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_08_nfw +=1.
                    
            
            if withinRange(L, [0., 1.25], 0.0):
                if m == color_yes:
                    corotate_1 +=1.
                    e_Lstar_1 += e_L**2

                if m != color_maybe:
                    total_1 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_1_steidel +=1.
                    e_Lstar_1_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_1_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_1_nfw +=1.
                    e_Lstar_1_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_1_nfw +=1.

            if withinRange(L, [0., 1.5], 0.0):
                if m == color_yes:
                    corotate_12 +=1.
                    e_Lstar_12 += e_L**2
                    print('Lstar: ',L)
                    print('e_Lstar_12: ',e_L)

                if m != color_maybe:
                    total_12 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_12_steidel +=1.
                    e_Lstar_12_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_12_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_12_nfw +=1.
                    e_Lstar_12_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_12_nfw +=1.
                    

            if withinRange(L, [0., 1.75], 0.0):
                if m == color_yes:
                    corotate_14 +=1.
                    e_Lstar_14 += e_L**2
                    print('Lstar: ',L)
                    print('e_Lstar_14: ',e_L)

                if m != color_maybe:
                    total_14 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_14_steidel +=1.
                    e_Lstar_14_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_14_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_14_nfw +=1.
                    e_Lstar_14_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_14_nfw +=1.
                    
                    
            if withinRange(L, [0., 100.0], 0.0):
                if m == color_yes:
                    corotate_16 +=1.
                    e_Lstar_16 += e_L**2
                    print('Lstar: ',L)
                    print('e_Lstar_16: ',e_L)

                if m != color_maybe:
                    total_16 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_16_steidel +=1.
                    e_Lstar_16_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_16_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_16_nfw +=1.
                    e_Lstar_16_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_16_nfw +=1.
                    
                    
                    
        print('total_02 = ',total_02)
        print('total_04 = ',total_04)
        print('total_06 = ',total_06)
        print('total_08 = ',total_08)
        print('total_1 = ',total_1)
        print('total_12 = ',total_12)
        print('total_14 = ',total_14)
        print('total_15 = ',total_16)
        print()
        print('largest Lstar: ',max(LstarList))
        
        
        e_Lstar_02 = np.sqrt(e_Lstar_02)
        e_Lstar_04 = np.sqrt(e_Lstar_04)
        e_Lstar_06 = np.sqrt(e_Lstar_06)
        e_Lstar_08 = np.sqrt(e_Lstar_08)
        e_Lstar_1 = np.sqrt(e_Lstar_1)
        e_Lstar_12 = np.sqrt(e_Lstar_12)
        e_Lstar_12 = np.sqrt(e_Lstar_14)
        e_Lstar_12 = np.sqrt(e_Lstar_16)

        xerr = [e_Lstar_02, e_Lstar_04, e_Lstar_06, e_Lstar_08, e_Lstar_1, e_Lstar_12, e_Lstar_14, e_Lstar_16]
        
        e_02 = np.sqrt(corotate_02)
        e_04 = np.sqrt(corotate_04)
        e_06 = np.sqrt(corotate_06)
        e_08 = np.sqrt(corotate_08)
        e_1 = np.sqrt(corotate_1)
        e_12 = np.sqrt(corotate_12)
        e_14 = np.sqrt(corotate_14)
        e_16 = np.sqrt(corotate_16)
        
        yerr = [e_02, e_04, e_06, e_08, e_1, e_12, e_14, e_16]

        lstars = np.arange(0.25, 2.25, 0.25)
        try:
            fraction = [corotate_02/total_02,
                        corotate_04/total_04,
                        corotate_06/total_06,
                        corotate_08/total_08,
                        corotate_1/total_1,
                        corotate_12/total_12,
                        corotate_14/total_14,
                        corotate_16/total_16]
                        
            fraction_str = [str(int(corotate_02)) + '/' + str(int(total_02)),
                            str(int(corotate_04)) + '/' + str(int(total_04)),
                            str(int(corotate_06)) + '/' + str(int(total_06)),
                            str(int(corotate_08)) + '/' + str(int(total_08)),
                            str(int(corotate_1)) + '/' + str(int(total_1)),
                            str(int(corotate_12)) + '/' + str(int(total_12)),
                            str(int(corotate_14)) + '/' + str(int(total_14)),
                            str(int(corotate_16)) + '/' + str(int(total_16))]

            fraction_steidel = [corotate_02_steidel/total_02_steidel,
                            corotate_04_steidel/total_04_steidel,
                            corotate_06_steidel/total_06_steidel,
                            corotate_08_steidel/total_08_steidel,
                            corotate_1_steidel/total_1_steidel,
                            corotate_12_steidel/total_12_steidel,
                            corotate_14_steidel/total_14_steidel,
                            corotate_16_steidel/total_16_steidel]
                        
            fraction_str_steidel = [str(int(corotate_02_steidel)) + '/' + str(int(total_02_steidel)),
                                str(int(corotate_04_steidel)) + '/' + str(int(total_04_steidel)),
                                str(int(corotate_06_steidel)) + '/' + str(int(total_06_steidel)),
                                str(int(corotate_08_steidel)) + '/' + str(int(total_08_steidel)),
                                str(int(corotate_1_steidel)) + '/' + str(int(total_1_steidel)),
                                str(int(corotate_12_steidel)) + '/' + str(int(total_12_steidel)),
                                str(int(corotate_14_steidel)) + '/' + str(int(total_14_steidel)),
                                str(int(corotate_16_steidel)) + '/' + str(int(total_16_steidel))]


            fraction_NFW = [corotate_02_nfw/total_02_nfw,
                            corotate_04_nfw/total_04_nfw,
                            corotate_06_nfw/total_06_nfw,
                            corotate_08_nfw/total_08_nfw,
                            corotate_1_nfw/total_1_nfw,
                            corotate_12_nfw/total_12_nfw,
                            corotate_14_nfw/total_14_nfw,
                            corotate_16_nfw/total_16_nfw]
                        
            fraction_str_NFW = [str(int(corotate_02_nfw)) + '/' + str(int(total_02_nfw)),
                                str(int(corotate_04_nfw)) + '/' + str(int(total_04_nfw)),
                                str(int(corotate_06_nfw)) + '/' + str(int(total_06_nfw)),
                                str(int(corotate_08_nfw)) + '/' + str(int(total_08_nfw)),
                                str(int(corotate_1_nfw)) + '/' + str(int(total_1_nfw)),
                                str(int(corotate_12_nfw)) + '/' + str(int(total_12_nfw)),
                                str(int(corotate_14_nfw)) + '/' + str(int(total_14_nfw)),
                                str(int(corotate_16_nfw)) + '/' + str(int(total_16_nfw))]

        except Exception as e:
            print('Exception: ',e)
            fraction = [0,
                        corotate_04/total_04,
                        corotate_06/total_06,
                        corotate_08/total_08,
                        corotate_1/total_1,
                        corotate_12/total_12,
                        corotate_14/total_14,
                        corotate_16/total_16]
                        
            sys.exit()

        print('minimum fraction: ',fraction)
        
        
        yerr = np.sqrt(np.array(fraction))*np.array(fraction)
        print('yerr: ',yerr)
        
        # plot apparent
        plot(lstars, fraction, lw=marker_lw, marker=marker_apparent, ls=ls_apparent, color=color_maybe, ms=marker_size, markeredgecolor='black')
        xTagOffset = -11
        yTagOffset = -18
        for l, f, f_str in zip(lstars, fraction, fraction_str):
            annotate(f_str,xy=(l, f),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)     
        
        
        # plot steidel model
        plot(lstars, fraction_steidel, lw=marker_lw, marker=marker_steidel, ls=ls_steidel, color=color_no, ms=marker_size, markeredgecolor='black')
        xTagOffset = -11
        yTagOffset = -18
        for l, f, f_str in zip(lstars, fraction_steidel, fraction_str_steidel):
            annotate(f_str,xy=(l, f),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)
        

        # plot NFW model
        plot(lstars, fraction_NFW, lw=marker_lw, marker=marker_NFW, ls= ls_NFW, color=color_yes, ms=marker_size, markeredgecolor='black')
        xTagOffset = -11
        yTagOffset = -18
        for l, f, f_str in zip(lstars, fraction_NFW, fraction_str_NFW):
            annotate(f_str,xy=(l, f),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)
        
        ylim(0, 1.05)
        xlim(0.125, 2.1)
        
#         hist(antirotate_inc, bins=bins, histtype='bar', lw=1.5, hatch='//', color = color_no, edgecolor='black',alpha=alpha_no, label=r'$\rm Anti-rotators$')
#         errorbar(lstars, fraction, yerr=yerr, lw=2, marker = 'D', color=color_yes, ms=15)

#         ylim(0, 12)
#         legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm L^{\**} ~ CDF $')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(0.25)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.25/2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls= ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls= ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls= ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
        plt.legend(handles=[NFW, steidel, apparent],loc='lower right', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_corotate_vs_Lstar_minimum_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_ALLmodel'.format(only_close_velocities, include_nondetection, L_limit, min_separation)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of binned Lstar, try a new way
#
#
##########################################################################################
##########################################################################################


    if plot_corotate_fraction_new:
        # initial figure
        fig = plt.figure(figsize=(7.7,5.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 2.0
        marker_size = 14
        
        bin_edges = equal_hist_edges(LstarList, 3)
        
        bin_edges = [0, 0.4, 0.7, 1.7, 6.0]
#         bin_edges = [0, 0.7, 6.0]

        right_edge = 2.5
        
        Lstar_bins, Lstar_bin_edges = np.histogram(LstarList, bins=bin_edges)
        
#         Lstar_bin_edges = append(Lstar_bin_edges, 6.0)

        print('Lstar_bins, Lstar_bin_edges : ',Lstar_bins, Lstar_bin_edges )
        print()
        print('bin_edges from equal_hist_edges: ',bin_edges)
        print()
        print('len(LstarList): ',len(LstarList))
        print()
        print('len(apparent_answers) ',len(apparent_answers))
        print('len(NFW_answers) ',len(NFW_answers))
        print('len(steidel_answers) ',len(steidel_answers))
        print('len(e_LstarList) ',len(e_LstarList))
        print('len(LstarList) ',len(LstarList))
        print('len(nameList) ',len(nameList))
        print()
        
        
        x_data = []
        y_apparent = []
        y_NFW = []
        y_steidel = []
        
        y_apparent_num = []
        y_NFW_num = []
        y_steidel_num = []
        
        i = 0
        while i <= len(Lstar_bins)-1:
            y_a = []
            y_n = []
            y_s = []
            x_d = []
            for a, nfw, steidel, e_L, L, name in zip(apparent_answers, NFW_answers, steidel_answers, e_LstarList, LstarList, nameList):
                print('L = ',L, 'Lstar_bin_edges[i] = ',Lstar_bin_edges[i], Lstar_bin_edges[i+1])
                print('name: ',name)
                if float(L) >= Lstar_bin_edges[i] and float(L) < Lstar_bin_edges[i+1]:
                    print('through: ',name)
                    
                    x_d.append(L)

                    if a:
                        y_a.append(1.)
                    else:
                        y_a.append(0)
                    
                    if nfw:
                        y_n.append(1.)
                    else:
                        y_n.append(0)
                    
                    if steidel:
                        y_s.append(1.)
                    else:
                        y_s.append(0)
                    
            print('y_s: ',y_s)
            print('np.mean(y_a): ',np.mean(y_a))
            print()

            y_apparent.append(np.mean(y_a))
            y_NFW.append(np.mean(y_n))
            y_steidel.append(np.mean(y_s))
            x_data.append(np.mean(x_d))
            
            y_apparent_num.append(len(y_a))
            y_NFW_num.append(len(y_n))
            y_steidel_num.append(len(y_s))
            
            i+=1
        
        x_axis = Lstar_bin_edges[1:]
        x_axis[-1:] = right_edge
        print('x_axis: ',x_axis)
        print('x_data: ',x_data)
        
        # --- Apparent results
        plot(x_data,
            y_apparent,
            lw=marker_lw,
            marker=marker_apparent,
            ls=ls_apparent,
            color=color_maybe,
            ms=marker_size,
            markeredgecolor='black')
            
#         xTagOffset = -11
#         yTagOffset = -18
#         for l, y, num in zip(Lstar_bin_edges[1:], y_apparent, y_apparent_num):
#             annotate(num,xy=(l, y),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)
        
        
        #  --- NFW model
        plot(x_data,
            y_NFW,
            lw=marker_lw,
            marker=marker_NFW,
            ls=ls_NFW,
            color=color_yes,
            ms=marker_size,
            markeredgecolor='black')

#         xTagOffset = -11
#         yTagOffset = -18
#         for l, y, num in zip(Lstar_bin_edges[1:], y_NFW, y_NFW_num):
#             annotate(num,xy=(l, y),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)
            
        
        # --- Steidel model
        plot(x_data,
            y_steidel,
            lw=marker_lw,
            marker=marker_steidel,
            ls=ls_steidel,
            color=color_no,
            ms=marker_size,
            markeredgecolor='black')
        
        xTagOffset = -11
        yTagOffset = -18
        for l, y, num in zip(x_data, y_steidel, y_steidel_num):
            annotate(num,xy=(l, y),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=12)
        
        
        ylim(0, 1.05)
        xlim(0, 2.25)
        
        xlabel(r'$\rm L^{\**}$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(0.25)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.25/2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        if plotGrid:
            ax.grid(True, alpha=0.5)
        
        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')

        plt.legend(handles=[NFW, steidel, apparent],loc='lower left', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
#         save_name = 'SALT_corotate_vs_Lstar_newbins'
        save_name = 'SALT_corotate_vs_Lstar_new_velstrict_{0}_non_{1}_minsep_{2}_neighbors_{3}'.format(only_close_velocities, include_nondetection, min_separation, neighbor_limit)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of binned impact parameter, try a new way
#
#
##########################################################################################
##########################################################################################


    if plot_corotate_impact_new:
        # initial figure
        fig = plt.figure(figsize=(7.7,5.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 2.0
        marker_size = 14
        
        bin_edges = equal_hist_edges(impactList, 3)
        
        bin_edges = [0, 75, 200, 600.]
#         bin_edges = [10.92, 121.02, 217.94, 600.]

        right_edge = 2.5
        
        impact_bins, impact_bin_edges = np.histogram(impactList, bins=bin_edges)
        
#         Lstar_bin_edges = append(Lstar_bin_edges, 6.0)

        print('impact_bins, impact_bin_edges : ',impact_bins, impact_bin_edges )
        print()
        print('bin_edges from equal_hist_edges: ',bin_edges)
        print()
        print('len(impactList): ',len(impactList))
        print('impactList: ',impactList)
        print()
        print('len(apparent_answers) ',len(apparent_answers))
        print('len(NFW_answers) ',len(NFW_answers))
        print('len(steidel_answers) ',len(steidel_answers))
        print('len(nameList) ',len(nameList))
        print()
        
        
        x_data = []
        y_apparent = []
        y_NFW = []
        y_steidel = []
        
        y_apparent_num = []
        y_NFW_num = []
        y_steidel_num = []
        
        i = 0
        while i <= len(impact_bins)-1:
            y_a = []
            y_n = []
            y_s = []
            x_d = []
            for a, nfw, steidel, imp, name in zip(apparent_answers, NFW_answers, steidel_answers, impactList, nameList):
                print('name: ',name)
                print('imp = ',imp, 'impact_bin_edges[i] = ',impact_bin_edges[i], impact_bin_edges[i+1])
                if float(imp) >= impact_bin_edges[i] and float(imp) < impact_bin_edges[i+1]:
                    print('through: ',name)
                    
                    x_d.append(imp)

                    if a:
                        y_a.append(1.)
                    else:
                        y_a.append(0)
                    
                    if nfw:
                        y_n.append(1.)
                    else:
                        y_n.append(0)
                    
                    if steidel:
                        y_s.append(1.)
                    else:
                        y_s.append(0)
                    
            print('y_s: ',y_s)
            print('np.mean(y_a): ',np.mean(y_a))
            print()

            y_apparent.append(np.mean(y_a))
            y_NFW.append(np.mean(y_n))
            y_steidel.append(np.mean(y_s))
            x_data.append(np.mean(x_d))
            
            y_apparent_num.append(len(y_a))
            y_NFW_num.append(len(y_n))
            y_steidel_num.append(len(y_s))
            
            i+=1
        
        x_axis = impact_bin_edges[1:]
        x_axis[-1:] = right_edge
        print('x_axis: ',x_axis)
        print('x_data: ',x_data)
        
        # --- Apparent results
        plot(x_data,
            y_apparent,
            lw=marker_lw,
            marker=marker_apparent,
            ls=ls_apparent,
            color=color_maybe,
            ms=marker_size,
            markeredgecolor='black')
            
#         xTagOffset = -11
#         yTagOffset = -18
#         for l, y, num in zip(Lstar_bin_edges[1:], y_apparent, y_apparent_num):
#             annotate(num,xy=(l, y),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)
        
        
        #  --- NFW model
        plot(x_data,
            y_NFW,
            lw=marker_lw,
            marker=marker_NFW,
            ls=ls_NFW,
            color=color_yes,
            ms=marker_size,
            markeredgecolor='black')

#         xTagOffset = -11
#         yTagOffset = -18
#         for l, y, num in zip(Lstar_bin_edges[1:], y_NFW, y_NFW_num):
#             annotate(num,xy=(l, y),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)
            
        
        # --- Steidel model
        plot(x_data,
            y_steidel,
            lw=marker_lw,
            marker=marker_steidel,
            ls=ls_steidel,
            color=color_no,
            ms=marker_size,
            markeredgecolor='black')
        
        xTagOffset = -11
        yTagOffset = -18
        for l, y, num in zip(x_data, y_steidel, y_steidel_num):
            annotate(num,xy=(l, y),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=12)
        
        
        ylim(0, 1.05)
        xlim(0, 350.0)
        
        xlabel(r'$\rm Impact~Parameter~[kpc]$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        if plotGrid:
            ax.grid(True, alpha=0.5)

        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')

        plt.legend(handles=[NFW, steidel, apparent],loc='lower left', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_corotate_vs_impact_new_velstrict_{0}_non_{1}_minsep_{2}_neighbors_{3}'.format(only_close_velocities, include_nondetection, min_separation, neighbor_limit)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of binned impact parameter / Rvir, try a new way
#
#
##########################################################################################
##########################################################################################


    if plot_corotate_impact_rvir_new:
        # initial figure
        fig = plt.figure(figsize=(7.7,5.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 2.0
        marker_size = 14
        
        bin_edges = equal_hist_edges(impactvirList, 3)
        print('bin_edges from equal_hist_edges: ',bin_edges)

        bin_edges = [0, 0.3, 0.7, 1.0, 1.3, 3.0]
#         bin_edges = [0, 0.5, 1.0, 1.5, 3.0]

#         bin_edges = [10.92, 121.02, 217.94, 600.]

        right_edge = 2.5
        
        impact_bins, impact_bin_edges = np.histogram(impactvirList, bins=bin_edges)
        
#         Lstar_bin_edges = append(Lstar_bin_edges, 6.0)

        print('impact_bins, impact_bin_edges : ',impact_bins, impact_bin_edges )
        print()
        print()
        print('len(impactvirList): ',len(impactvirList))
        print('impactvirList: ',impactvirList)
        print()
        print('len(apparent_answers) ',len(apparent_answers))
        print('len(NFW_answers) ',len(NFW_answers))
        print('len(steidel_answers) ',len(steidel_answers))
        print('len(nameList) ',len(nameList))
        print()
        
        
        x_data = []
        y_apparent = []
        y_NFW = []
        y_steidel = []
        
        y_apparent_num = []
        y_NFW_num = []
        y_steidel_num = []
        
        i = 0
        while i <= len(impact_bins)-1:
            y_a = []
            y_n = []
            y_s = []
            x_d = []
            for a, nfw, steidel, imp, name in zip(apparent_answers, NFW_answers, steidel_answers, impactvirList, nameList):
                print('name: ',name)
                print('imp = ',imp, 'impact_bin_edges[i] = ',impact_bin_edges[i], impact_bin_edges[i+1])
                if float(imp) >= impact_bin_edges[i] and float(imp) < impact_bin_edges[i+1]:
                    print('through: ',name)
                    
                    x_d.append(imp)

                    if a:
                        y_a.append(1.)
                    else:
                        y_a.append(0)
                    
                    if nfw:
                        y_n.append(1.)
                    else:
                        y_n.append(0)
                    
                    if steidel:
                        y_s.append(1.)
                    else:
                        y_s.append(0)
                    
            print('y_s: ',y_s)
            print('np.mean(y_a): ',np.mean(y_a))
            print()

            y_apparent.append(np.mean(y_a))
            y_NFW.append(np.mean(y_n))
            y_steidel.append(np.mean(y_s))
            x_data.append(np.mean(x_d))
            
            y_apparent_num.append(len(y_a))
            y_NFW_num.append(len(y_n))
            y_steidel_num.append(len(y_s))
            
            i+=1
        
        x_axis = impact_bin_edges[1:]
        x_axis[-1:] = right_edge
        print('x_axis: ',x_axis)
        print('x_data: ',x_data)
        
        # --- Apparent results
        plot(x_data,
            y_apparent,
            lw=marker_lw,
            marker=marker_apparent,
            ls=ls_apparent,
            color=color_maybe,
            ms=marker_size,
            markeredgecolor='black')
            
#         xTagOffset = -11
#         yTagOffset = -18
#         for l, y, num in zip(Lstar_bin_edges[1:], y_apparent, y_apparent_num):
#             annotate(num,xy=(l, y),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)
        
        
        #  --- NFW model
        plot(x_data,
            y_NFW,
            lw=marker_lw,
            marker=marker_NFW,
            ls=ls_NFW,
            color=color_yes,
            ms=marker_size,
            markeredgecolor='black')

#         xTagOffset = -11
#         yTagOffset = -18
#         for l, y, num in zip(Lstar_bin_edges[1:], y_NFW, y_NFW_num):
#             annotate(num,xy=(l, y),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=9)
            
        
        # --- Steidel model
        plot(x_data,
            y_steidel,
            lw=marker_lw,
            marker=marker_steidel,
            ls=ls_steidel,
            color=color_no,
            ms=marker_size,
            markeredgecolor='black')
        
        xTagOffset = -11
        yTagOffset = -18
        for l, y, num in zip(x_data, y_steidel, y_steidel_num):
            annotate(num,xy=(l, y),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=12)
        
        
        ylim(0, 1.05)
        xlim(0, 2.5)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        if plotGrid:
            ax.grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')

        plt.legend(handles=[NFW, steidel, apparent],loc='lower right', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_corotate_vs_impact_rvir_new_velstrict_{0}_non_{1}_minsep_{2}_neighbors_{3}'.format(only_close_velocities, include_nondetection, min_separation, neighbor_limit)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')



##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of binned azimuth angle
#
#
##########################################################################################
##########################################################################################


    if plot_corotate_azimuth:
        # initial figure
        fig = plt.figure(figsize=(7.7,5.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 2.0
        marker_size = 14
        
        bin_edges = equal_hist_edges(azList, 5)
        print('bin_edges from equal_hist_edges: ',bin_edges)

#         bin_edges = [0, 0.3, 0.7, 1.0, 1.3, 3.0]

        right_edge = 91.
        
        az_bins, az_bin_edges = np.histogram(azList, bins=bin_edges)
        
#         Lstar_bin_edges = append(Lstar_bin_edges, 6.0)

        print('az_bins, az_bin_edges : ',az_bins, az_bin_edges)
        print()
        print()
        print('len(azrList): ',len(azList))
        print('azList: ',azList)
        print()
        print('len(apparent_answers) ',len(apparent_answers))
        print('len(NFW_answers) ',len(NFW_answers))
        print('len(steidel_answers) ',len(steidel_answers))
        print('len(nameList) ',len(nameList))
        print()
        
        
        x_data = []
        y_apparent = []
        y_NFW = []
        y_steidel = []
        
        y_apparent_num = []
        y_NFW_num = []
        y_steidel_num = []
        
        i = 0
        while i <= len(az_bins)-1:
            y_a = []
            y_n = []
            y_s = []
            x_d = []
            for a, nfw, steidel, az, name in zip(apparent_answers, NFW_answers, steidel_answers, azList, nameList):
                print('name: ',name)
                print('az = ',az, 'az_bin_edges[i] = ',az_bin_edges[i], az_bin_edges[i+1])
                if float(az) >= az_bin_edges[i] and float(az) < az_bin_edges[i+1]:
                    print('through: ',name)
                    
                    x_d.append(az)

                    if a:
                        y_a.append(1.)
                    else:
                        y_a.append(0)
                    
                    if nfw:
                        y_n.append(1.)
                    else:
                        y_n.append(0)
                    
                    if steidel:
                        y_s.append(1.)
                    else:
                        y_s.append(0)
                    
            print('y_s: ',y_s)
            print('np.mean(y_a): ',np.mean(y_a))
            print()

            y_apparent.append(np.mean(y_a))
            y_NFW.append(np.mean(y_n))
            y_steidel.append(np.mean(y_s))
            x_data.append(np.mean(x_d))
            
            y_apparent_num.append(len(y_a))
            y_NFW_num.append(len(y_n))
            y_steidel_num.append(len(y_s))
            
            i+=1
        
        x_axis = az_bin_edges[1:]
        x_axis[-1:] = right_edge
        print('x_axis: ',x_axis)
        print('x_data: ',x_data)
        
        # --- Apparent results
        plot(x_data,
            y_apparent,
            lw=marker_lw,
            marker=marker_apparent,
            ls=ls_apparent,
            color=color_maybe,
            ms=marker_size,
            markeredgecolor='black')
        
        
        #  --- NFW model
        plot(x_data,
            y_NFW,
            lw=marker_lw,
            marker=marker_NFW,
            ls=ls_NFW,
            color=color_yes,
            ms=marker_size,
            markeredgecolor='black')
            
        
        # --- Steidel model
        plot(x_data,
            y_steidel,
            lw=marker_lw,
            marker=marker_steidel,
            ls=ls_steidel,
            color=color_no,
            ms=marker_size,
            markeredgecolor='black')
        
        xTagOffset = -11
        yTagOffset = -18
        for l, y, num in zip(x_data, y_steidel, y_steidel_num):
            annotate(num,xy=(l, y),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=12)
        
        
        ylim(0, 1.05)
        xlim(0, 90)
        
        xlabel(r'$\rm Azimuth~[deg]$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        if plotGrid:
            ax.grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model$')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent$')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model$')

        plt.legend(handles=[NFW, steidel, apparent],loc='lower right', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_corotate_vs_azimuth_new_velstrict_{0}_non_{1}_minsep_{2}_neighbors_{3}'.format(only_close_velocities, include_nondetection, min_separation, neighbor_limit)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of binned inclination angle
#
#
##########################################################################################
##########################################################################################


    if plot_corotate_inc:
        # initial figure
        fig = plt.figure(figsize=(7.7,5.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 2.0
        marker_size = 14
        
        max_az = 90
        
        bin_edges = equal_hist_edges(incList, 4)
        print('bin_edges from equal_hist_edges: ',bin_edges)

        bin_edges = np.array([0, 50, 75, 90])

        right_edge = 90.
        
        bin_edges[-1:] = 90
        print('bin edges now: ',bin_edges)
        
        inc_bins, inc_bin_edges = np.histogram(incList, bins=bin_edges)
        
#         Lstar_bin_edges = append(Lstar_bin_edges, 6.0)

        print('inc_bins, inc_bin_edges : ',inc_bins, inc_bin_edges)
        print()
        print()
        print('len(incrList): ',len(incList))
        print('incList: ',incList)
        print()
        print('len(apparent_answers) ',len(apparent_answers))
        print('len(NFW_answers) ',len(NFW_answers))
        print('len(steidel_answers) ',len(steidel_answers))
        print('len(nameList) ',len(nameList))
        print()
        
        
        x_data = []
        y_apparent = []
        y_NFW = []
        y_steidel = []
        
        y_apparent_num = []
        y_NFW_num = []
        y_steidel_num = []
        
        i = 0
        while i <= len(inc_bins)-1:
            y_a = []
            y_n = []
            y_s = []
            x_d = []
            for a, nfw, steidel, inc, name, az in zip(apparent_answers, NFW_answers, steidel_answers, incList, nameList, azList):
                print('name: ',name)
                print('inc = ',inc, 'inc_bin_edges[i] = ',inc_bin_edges[i], inc_bin_edges[i+1])
                if float(inc) >= inc_bin_edges[i] and float(inc) < inc_bin_edges[i+1]:
                    if az <= max_az:
                        print('through: ',name)
                    
                        x_d.append(inc)

                        if a:
                            y_a.append(1.)
                        else:
                            y_a.append(0)
                    
                        if nfw:
                            y_n.append(1.)
                        else:
                            y_n.append(0)
                    
                        if steidel:
                            y_s.append(1.)
                        else:
                            y_s.append(0)
                    
            print('y_s: ',y_s)
            print('np.mean(y_a): ',np.mean(y_a))
            print()

            y_apparent.append(np.mean(y_a))
            y_NFW.append(np.mean(y_n))
            y_steidel.append(np.mean(y_s))
            x_data.append(np.mean(x_d))
            
            y_apparent_num.append(len(y_a))
            y_NFW_num.append(len(y_n))
            y_steidel_num.append(len(y_s))
            
            i+=1
        
        x_axis = inc_bin_edges[1:]
        x_axis[-1:] = right_edge
        print('x_axis: ',x_axis)
        print('x_data: ',x_data)
        
        # --- Apparent results
        plot(x_data,
            y_apparent,
            lw=marker_lw,
            marker=marker_apparent,
            ls=ls_apparent,
            color=color_maybe,
            ms=marker_size,
            markeredgecolor='black')
        
        
        #  --- NFW model
        plot(x_data,
            y_NFW,
            lw=marker_lw,
            marker=marker_NFW,
            ls=ls_NFW,
            color=color_yes,
            ms=marker_size,
            markeredgecolor='black')
            
        
        # --- Steidel model
        plot(x_data,
            y_steidel,
            lw=marker_lw,
            marker=marker_steidel,
            ls=ls_steidel,
            color=color_no,
            ms=marker_size,
            markeredgecolor='black')

        xTagOffset = -11
        yTagOffset = -18
        for l, y, num in zip(x_data, y_steidel, y_steidel_num):
            annotate(num,xy=(l, y),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=12)
        
        
        ylim(0, 1.05)
        xlim(0, 90)
        
        xlabel(r'$\rm Inclination~[deg]$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        if plotGrid:
            ax.grid(True, alpha=0.5)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')

        plt.legend(handles=[NFW, steidel, apparent],loc='lower left', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_corotate_vs_inc_new_velstrict_{0}_non_{1}_minsep_{2}_neighbors_{3}_maxaz_{4}'.format(only_close_velocities, include_nondetection, min_separation, neighbor_limit, max_az)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of binned column density
#
#
##########################################################################################
##########################################################################################


    if plot_corotate_N:
        # initial figure
        fig = plt.figure(figsize=(7.7,5.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 2.0
        marker_size = 14
        
        max_az = 90
        
        bin_edges = equal_hist_edges(NList, 4)
        print('bin_edges from equal_hist_edges: ',bin_edges)

#         bin_edges = np.array([0, 50, 75, 90])

        right_edge = 21.0
        
        bin_edges[-1:] = right_edge
        print('bin edges now: ',bin_edges)
        
        N_bins, N_bin_edges = np.histogram(NList, bins=bin_edges)
        
#         Lstar_bin_edges = append(Lstar_bin_edges, 6.0)

        print('N_bins, N_bin_edges : ',N_bins, N_bin_edges)
        print()
        print()
        print('len(NrList): ',len(NList))
        print('NList: ',NList)
        print()
        print('len(apparent_answers) ',len(apparent_answers))
        print('len(NFW_answers) ',len(NFW_answers))
        print('len(steidel_answers) ',len(steidel_answers))
        print('len(nameList) ',len(nameList))
        print()
        
        
        x_data = []
        y_apparent = []
        y_NFW = []
        y_steidel = []
        
        y_apparent_num = []
        y_NFW_num = []
        y_steidel_num = []
        
        i = 0
        while i <= len(N_bins)-1:
            y_a = []
            y_n = []
            y_s = []
            x_d = []
            for a, nfw, steidel, N, name, az in zip(apparent_answers, NFW_answers, steidel_answers, NList, nameList, azList):
                print('name: ',name)
                print('N = ',N, 'N_bin_edges[i] = ',N_bin_edges[i], N_bin_edges[i+1])
                if float(N) >= N_bin_edges[i] and float(N) < N_bin_edges[i+1]:
                    if az <= max_az:
                        print('through: ',name)
                    
                        x_d.append(N)

                        if a:
                            y_a.append(1.)
                        else:
                            y_a.append(0)
                    
                        if nfw:
                            y_n.append(1.)
                        else:
                            y_n.append(0)
                    
                        if steidel:
                            y_s.append(1.)
                        else:
                            y_s.append(0)
                    
            print('y_s: ',y_s)
            print('np.mean(y_a): ',np.mean(y_a))
            print()

            y_apparent.append(np.mean(y_a))
            y_NFW.append(np.mean(y_n))
            y_steidel.append(np.mean(y_s))
            x_data.append(np.mean(x_d))
            
            y_apparent_num.append(len(y_a))
            y_NFW_num.append(len(y_n))
            y_steidel_num.append(len(y_s))
            
            i+=1
        
        x_axis = N_bin_edges[1:]
        x_axis[-1:] = right_edge
        print('x_axis: ',x_axis)
        print('x_data: ',x_data)
        
        # --- Apparent results
        plot(x_data,
            y_apparent,
            lw=marker_lw,
            marker=marker_apparent,
            ls=ls_apparent,
            color=color_maybe,
            ms=marker_size,
            markeredgecolor='black')
        
        
        #  --- NFW model
        plot(x_data,
            y_NFW,
            lw=marker_lw,
            marker=marker_NFW,
            ls=ls_NFW,
            color=color_yes,
            ms=marker_size,
            markeredgecolor='black')

        xTagOffset = -11
        yTagOffset = -18
        for l, y, num in zip(x_data, y_NFW, y_NFW_num):
            annotate(num,xy=(l, y),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=12)
            
        
        # --- Steidel model
        plot(x_data,
            y_steidel,
            lw=marker_lw,
            marker=marker_steidel,
            ls=ls_steidel,
            color=color_no,
            ms=marker_size,
            markeredgecolor='black')
        
        
        ylim(0, 1.05)
        xlim(12, 16)
        
        xlabel(r'$\rm log \emph{N}$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(1)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(0.5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        if plotGrid:
            ax.grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')

        plt.legend(handles=[NFW, steidel, apparent],loc='lower left', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_corotate_vs_N_new_velstrict_{0}_non_{1}_minsep_{2}_neighbors_{3}_maxaz_{4}'.format(only_close_velocities, include_nondetection, min_separation, neighbor_limit, max_az)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of V_max
#
#
##########################################################################################
##########################################################################################


    if plot_corotate_vmax:
        # initial figure
        fig = plt.figure(figsize=(7.7,5.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 2.0
        marker_size = 14
        
        max_az = 90
        
        vmaxList = abs(np.array(vmaxList))
        
        bin_edges = equal_hist_edges(vmaxList, 6)
        print('bin_edges from equal_hist_edges: ',bin_edges)

#         bin_edges = np.array([0, 50, 75, 90])

        right_edge = 400.0
        
        bin_edges[-1:] = right_edge
        print('bin edges now: ',bin_edges)
        
        vmax_bins, vmax_bin_edges = np.histogram(vmaxList, bins=bin_edges)
        
#         Lstar_bin_edges = append(Lstar_bin_edges, 6.0)

        print('vmax_bins, vmax_bin_edges : ',vmax_bins, vmax_bin_edges)
        print()
        print()
        print('len(vmaxList): ',len(vmaxList))
        print('vmaxList: ',vmaxList)
        print()
        print('len(apparent_answers) ',len(apparent_answers))
        print('len(NFW_answers) ',len(NFW_answers))
        print('len(steidel_answers) ',len(steidel_answers))
        print('len(nameList) ',len(nameList))
        print()
        
        
        x_data = []
        y_apparent = []
        y_NFW = []
        y_steidel = []
        
        y_apparent_num = []
        y_NFW_num = []
        y_steidel_num = []
        
        i = 0
        while i <= len(vmax_bins)-1:
            y_a = []
            y_n = []
            y_s = []
            x_d = []
            for a, nfw, steidel, N, name, az, vmax in zip(apparent_answers, NFW_answers, steidel_answers, NList, nameList, azList, vmaxList):
                print('name: ',name)
                print('vmax = ',vmax, 'vmax_bin_edges[i] = ',vmax_bin_edges[i], vmax_bin_edges[i+1])
                if float(vmax) >= vmax_bin_edges[i] and float(vmax) < vmax_bin_edges[i+1]:
                    if az <= max_az:
                        print('through: ',name)
                    
                        x_d.append(vmax)

                        if a:
                            y_a.append(1.)
                        else:
                            y_a.append(0)
                    
                        if nfw:
                            y_n.append(1.)
                        else:
                            y_n.append(0)
                    
                        if steidel:
                            y_s.append(1.)
                        else:
                            y_s.append(0)
                    
            print('y_s: ',y_s)
            print('np.mean(y_a): ',np.mean(y_a))
            print()

            y_apparent.append(np.mean(y_a))
            y_NFW.append(np.mean(y_n))
            y_steidel.append(np.mean(y_s))
            x_data.append(np.mean(x_d))
            
            y_apparent_num.append(len(y_a))
            y_NFW_num.append(len(y_n))
            y_steidel_num.append(len(y_s))
            
            i+=1
        
        x_axis = N_bin_edges[1:]
        x_axis[-1:] = right_edge
        print('x_axis: ',x_axis)
        print('x_data: ',x_data)
        
        # --- Apparent results
        plot(x_data,
            y_apparent,
            lw=marker_lw,
            marker=marker_apparent,
            ls=ls_apparent,
            color=color_maybe,
            ms=marker_size,
            markeredgecolor='black')
        
        
        #  --- NFW model
        plot(x_data,
            y_NFW,
            lw=marker_lw,
            marker=marker_NFW,
            ls=ls_NFW,
            color=color_yes,
            ms=marker_size,
            markeredgecolor='black')

        xTagOffset = -11
        yTagOffset = -18
        for l, y, num in zip(x_data, y_NFW, y_NFW_num):
            annotate(num,xy=(l, y),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=12)
            
        
        # --- Steidel model
        plot(x_data,
            y_steidel,
            lw=marker_lw,
            marker=marker_steidel,
            ls=ls_steidel,
            color=color_no,
            ms=marker_size,
            markeredgecolor='black')
        
        
        ylim(0, 1.05)
        xlim(0, 300)
        
        xlabel(r'$\rm \emph{v}_{max}~[km~s^{-1}]$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        if plotGrid:
            ax.grid(True, alpha=0.5)
            
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')

        plt.legend(handles=[NFW, steidel, apparent],loc='lower left', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'SALT_corotate_vs_vmax_new_velstrict_{0}_non_{1}_minsep_{2}_neighbors_{3}_maxaz_{4}'.format(only_close_velocities, include_nondetection, min_separation, neighbor_limit, max_az)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of binned Lstar
#
#
##########################################################################################
##########################################################################################

    if plot_corotate_fraction_old:
        # initial figure
        fig = plt.figure(figsize=(7.7,6.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 1.5
        marker_size = 12

        L_limit = 0.5

        # no model
        
        corotate_04 = 0
        e_Lstar_04 = 0
        total_04 = 0

        corotate_08 = 0
        e_Lstar_08 = 0
        total_08 = 0
        
        corotate_12 = 0
        e_Lstar_12 = 0
        total_12 = 0

        corotate_16 = 0
        e_Lstar_16 = 0
        total_16 = 0
        
        
        # steidel model
        corotate_04_steidel = 0
        e_Lstar_04_steidel = 0
        total_04_steidel = 0

        corotate_08_steidel = 0
        e_Lstar_08_steidel = 0
        total_08_steidel = 0

        corotate_12_steidel = 0
        e_Lstar_12_steidel = 0
        total_12_steidel = 0
    
        corotate_16_steidel = 0
        e_Lstar_16_steidel = 0
        total_16_steidel = 0
        
        
        # NFW model
        corotate_04_nfw = 0
        e_Lstar_04_nfw = 0
        total_04_nfw = 0

        corotate_08_nfw = 0
        e_Lstar_08_nfw = 0
        total_08_nfw = 0
        
        corotate_12_nfw = 0
        e_Lstar_12_nfw = 0
        total_12_nfw = 0

        corotate_16_nfw = 0
        e_Lstar_16_nfw = 0
        total_16_nfw = 0
        
        
        for m, m_nfw, m_steidel, e_L, ra, dec, L in zip(markerColorList, markerColorList_NFW, markerColorList_steidel, e_LstarList, RA_targetList, Dec_targetList, LstarList):
            
            if withinRange(L, [0., 0.5], 0.0):
                if m == color_yes:
                    corotate_04 +=1.
                    e_Lstar_04 += e_L**2

                if m != color_maybe:
                    total_04 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_04_steidel +=1.
                    e_Lstar_04_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_04_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_04_nfw +=1.
                    e_Lstar_04_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_04_nfw +=1.
                    

            if withinRange(L, [0.50001, 1.0], 0.0):
                if m == color_yes:
                    corotate_08 +=1.
                    e_Lstar_08 += e_L**2

                if m != color_maybe:
                    total_08 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_08_steidel +=1.
                    e_Lstar_08_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_08_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_08_nfw +=1.
                    e_Lstar_08_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_08_nfw +=1.
                    

            if withinRange(L, [1.00001, 1.7], 0.0):
                if m == color_yes:
                    corotate_12 +=1.
                    e_Lstar_12 += e_L**2
                    print('Lstar: ',L)
                    print('e_Lstar_12: ',e_L)

                if m != color_maybe:
                    total_12 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_12_steidel +=1.
                    e_Lstar_12_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_12_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_12_nfw +=1.
                    e_Lstar_12_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_12_nfw +=1.
                    
                    
            if withinRange(L, [1.70001, 100.0], 0.0):
                if m == color_yes:
                    corotate_16 +=1.
                    e_Lstar_16 += e_L**2
                    print('Lstar: ',L)
                    print('e_Lstar_16: ',e_L)

                if m != color_maybe:
                    total_16 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_16_steidel +=1.
                    e_Lstar_16_steidel += e_L**2
                    
                if m_steidel != color_maybe:
                    total_16_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_16_nfw +=1.
                    e_Lstar_16_nfw += e_L**2
                    
                if m_nfw != color_maybe:
                    total_16_nfw +=1.
                    
                    
                    
        print('total_04 = ',total_04)
        print('total_08 = ',total_08)
        print('total_12 = ',total_12)
        print('total_15 = ',total_16)
        print()
        print('largest Lstar: ',max(LstarList))
        
        
        e_Lstar_04 = np.sqrt(e_Lstar_04)
        e_Lstar_08 = np.sqrt(e_Lstar_08)
        e_Lstar_12 = np.sqrt(e_Lstar_12)
        e_Lstar_12 = np.sqrt(e_Lstar_16)

        xerr = [e_Lstar_04, e_Lstar_08, e_Lstar_12, e_Lstar_16]
        
        e_04 = np.sqrt(corotate_04)
        e_08 = np.sqrt(corotate_08)
        e_12 = np.sqrt(corotate_12)
        e_16 = np.sqrt(corotate_16)
        
        yerr = [e_04, e_08, e_12, e_16]

#         lstars = np.arange(0.4, 2.0, 0.4)
        lstars = np.array([0.5, 1.0, 1.7, 2.4])

        try:
            fraction = [corotate_04/total_04,
                        corotate_08/total_08,
                        corotate_12/total_12,
                        corotate_16/total_16]
                        
            fraction_str = [str(int(corotate_04)) + '/' + str(int(total_04)),
                            str(int(corotate_08)) + '/' + str(int(total_08)),
                            str(int(corotate_12)) + '/' + str(int(total_12)),
                            str(int(corotate_16)) + '/' + str(int(total_16))]

            fraction_steidel = [corotate_04_steidel/total_04_steidel,
                            corotate_08_steidel/total_08_steidel,
                            corotate_12_steidel/total_12_steidel,
                            corotate_16_steidel/total_16_steidel]
                        
            fraction_str_steidel = [str(int(corotate_04_steidel)) + '/' + str(int(total_04_steidel)),
                                str(int(corotate_08_steidel)) + '/' + str(int(total_08_steidel)),
                                str(int(corotate_12_steidel)) + '/' + str(int(total_12_steidel)),
                                str(int(corotate_16_steidel)) + '/' + str(int(total_16_steidel))]


            fraction_NFW = [corotate_04_nfw/total_04_nfw,
                            corotate_08_nfw/total_08_nfw,
                            corotate_12_nfw/total_12_nfw,
                            corotate_16_nfw/total_16_nfw]
                        
            fraction_str_NFW = [str(int(corotate_04_nfw)) + '/' + str(int(total_04_nfw)),
                                str(int(corotate_08_nfw)) + '/' + str(int(total_08_nfw)),
                                str(int(corotate_12_nfw)) + '/' + str(int(total_12_nfw)),
                                str(int(corotate_16_nfw)) + '/' + str(int(total_16_nfw))]

        except Exception as e:
            print('Exception: ',e)
            fraction = [corotate_04/total_04,
                        corotate_08/total_08,
                        corotate_12/total_12,
                        corotate_16/total_16]
                        
            sys.exit()

        print('minimum fraction: ',fraction)
        
        
        yerr = np.sqrt(np.array(fraction))*np.array(fraction)
        print('yerr: ',yerr)
        
        # plot apparent
        plot(lstars, fraction, lw=marker_lw, marker=marker_apparent, ls=ls_apparent, color=color_maybe, ms=marker_size, markeredgecolor='black')
        xTagOffset = -18
        yTagOffset = -14
        for l, f, f_str in zip(lstars, fraction, fraction_str):
            annotate(f_str,xy=(l, f),\
            xytext=(-len(f_str)+xTagOffset, yTagOffset),textcoords='offset points',size=9)     
        
        
        # plot Steidel model
        print('marker_steidel: ',marker_steidel)
        plot(lstars, fraction_steidel, lw=marker_lw, marker=marker_steidel, ls=ls_steidel, color=color_no, ms=marker_size, markeredgecolor='black')
        xTagOffset = -18
        yTagOffset = -14
        for l, f, f_str in zip(lstars, fraction_steidel, fraction_str_steidel):
            annotate(f_str,xy=(l, f),\
            xytext=(-len(f_str)+xTagOffset, yTagOffset),textcoords='offset points',size=9)
        

        # plot NFW model
        plot(lstars, fraction_NFW, lw=marker_lw, ls=ls_NFW, marker=marker_NFW, color=color_yes, ms=marker_size, markeredgecolor='black')
        xTagOffset = -18
        yTagOffset = -14
        for l, f, f_str in zip(lstars, fraction_NFW, fraction_str_NFW):
            print('f_str, len(f_str): ',f_str, len(f_str))
            annotate(f_str,xy=(l, f),\
            xytext=(-len(f_str)+xTagOffset, yTagOffset),textcoords='offset points',size=9)
        
        ylim(0, 1.05)
        xlim(0.125, 2.6)
        
#         hist(antirotate_inc, bins=bins, histtype='bar', lw=1.5, hatch='//', color = color_no, edgecolor='black',alpha=alpha_no, label=r'$\rm Anti-rotators$')
#         errorbar(lstars, fraction, yerr=yerr, lw=2, marker = 'D', color=color_yes, ms=15)

#         ylim(0, 12)
#         legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm L^{\**}$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(0.25)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.25/2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
        plt.legend(handles=[NFW, steidel, apparent],loc='lower right', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALT_corotate_vs_Lstar_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_ALLmodel'.format(only_close_velocities, include_nondetection, L_limit, min_separation)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')



##########################################################################################
##########################################################################################
##########################################################################################
# Plot co-rotation fraction as a function of binned impact / Rvir
#
#
##########################################################################################
##########################################################################################

    if plot_corotate_fraction_dist_old:
        # initial figure
        fig = plt.figure(figsize=(7.7,6.7))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
#         bins = arange(0, 100, 100)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 1.5
        marker_size = 12

        L_limit = 0.5

        # no model
        corotate_05 = 0
        total_05 = 0
        
        corotate_15 = 0
        total_15 = 0
        
        corotate_3 = 0
        total_3 = 0
        
        # Steidel model
        corotate_05_steidel = 0
        total_05_steidel = 0

        corotate_15_steidel = 0
        total_15_steidel = 0

        corotate_3_steidel = 0
        total_3_steidel = 0
        
        # NFW model
        corotate_05_nfw = 0
        total_05_nfw = 0
        
        corotate_15_nfw = 0
        total_15_nfw = 0
        
        corotate_3_nfw = 0
        total_3_nfw = 0
        
        
        for m, m_nfw, m_steidel, imp, ra, dec, L in zip(markerColorList, markerColorList_NFW, markerColorList_steidel, impactvirList, RA_targetList, Dec_targetList, LstarList):
            if withinRange(imp, [0., 0.5], 0.0):
                if m == color_yes:
                    corotate_05 +=1.

                if m != color_maybe:
                    total_05 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_05_steidel +=1.
                    
                if m_steidel != color_maybe:
                    total_05_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_05_nfw +=1.
                    
                if m_nfw != color_maybe:
                    total_05_nfw +=1.
                    

            if withinRange(imp, [0.50001, 1.5], 0.0):
                if m == color_yes:
                    corotate_15 +=1.

                if m != color_maybe:
                    total_15 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_15_steidel +=1.
                    
                if m_steidel != color_maybe:
                    total_15_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_15_nfw +=1.
                    
                if m_nfw != color_maybe:
                    total_15_nfw +=1.

                    
            if withinRange(imp, [1.50001, 3.0], 0.0):
                if m == color_yes:
                    corotate_3 +=1.

                if m != color_maybe:
                    total_3 +=1.
                    
                # steidel model
                if m_steidel == color_yes:
                    corotate_3_steidel +=1.
                    
                if m_steidel != color_maybe:
                    total_3_steidel +=1.
                    
                # NFW model
                if m_nfw == color_yes:
                    corotate_3_nfw +=1.
                    
                if m_nfw != color_maybe:
                    total_3_nfw +=1.
                    
                    
                    
                    
        print('total_05 = ',total_05)
        print('total_15 = ',total_15)
        print('total_3 = ',total_3)
        print()

#         xdata = np.arange(0.5, 3.5, 0.5)
        xdata = np.array([0.5, 1.5, 3.0])

        try:
            fraction = [corotate_05/total_05,
                        corotate_15/total_15,
                        corotate_3/total_3]
                        
            fraction_str = [str(int(corotate_05)) + '/' + str(int(total_05)),
                            str(int(corotate_15)) + '/' + str(int(total_15)),
                            str(int(corotate_3)) + '/' + str(int(total_3))]

            fraction_steidel = [corotate_05_steidel/total_05_steidel,
                            corotate_15_steidel/total_15_steidel,
                            corotate_3_steidel/total_3_steidel]
                        
            fraction_str_steidel = [str(int(corotate_05_steidel)) + '/' + str(int(total_05_steidel)),
                                str(int(corotate_15_steidel)) + '/' + str(int(total_15_steidel)),
                                str(int(corotate_3_steidel)) + '/' + str(int(total_3_steidel))]


            fraction_NFW = [corotate_05_nfw/total_05_nfw,
                            corotate_15_nfw/total_15_nfw,
                            corotate_3_nfw/total_3_nfw]
                        
            fraction_str_NFW = [str(int(corotate_05_nfw)) + '/' + str(int(total_05_nfw)),
                                str(int(corotate_15_nfw)) + '/' + str(int(total_15_nfw)),
                                str(int(corotate_3_nfw)) + '/' + str(int(total_3_nfw))]

        except Exception as e:
            print('Exception: ',e)
            
            sys.exit()
            
        print('minimum fraction: ',fraction)
        
                
        # plot apparent
        plot(xdata, fraction, lw=marker_lw, marker=marker_apparent, ls=ls_apparent, color=color_maybe, ms=marker_size, markeredgecolor='black')
        xTagOffset = -20
        yTagOffset = -14
        for l, f, f_str in zip(xdata, fraction, fraction_str):
            annotate(f_str,xy=(l, f),\
            xytext=(-len(f_str)+xTagOffset, yTagOffset),textcoords='offset points',size=9)     
        
        
        # plot Steidel model
        print('m_steidel: ',m_steidel)
        plot(xdata, fraction_steidel, lw=marker_lw, marker=marker_steidel, ls=ls_steidel, color=color_no, ms=marker_size, markeredgecolor='black')
        xTagOffset = -20
        yTagOffset = -14
        for l, f, f_str in zip(xdata, fraction_steidel, fraction_str_steidel):
            annotate(f_str,xy=(l, f),\
            xytext=(-len(f_str)+xTagOffset, yTagOffset),textcoords='offset points',size=9)
        

        # plot NFW model
        plot(xdata, fraction_NFW, lw=marker_lw, ls=ls_NFW, marker=marker_NFW, color=color_yes, ms=marker_size, markeredgecolor='black')
        xTagOffset = -20
        yTagOffset = -14
        for l, f, f_str in zip(xdata, fraction_NFW, fraction_str_NFW):
            print('f_str, len(f_str): ',f_str, len(f_str))
            annotate(f_str,xy=(l, f),\
            xytext=(-len(f_str)+xTagOffset, yTagOffset),textcoords='offset points',size=9)
        
        ylim(0, 1.05)
        xlim(0.125, 3.1)
        
#         hist(antirotate_inc, bins=bins, histtype='bar', lw=1.5, hatch='//', color = color_no, edgecolor='black',alpha=alpha_no, label=r'$\rm Anti-rotators$')
#         errorbar(lstars, fraction, yerr=yerr, lw=2, marker = 'D', color=color_yes, ms=15)

#         ylim(0, 12)
#         legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Co-rotating ~ Fraction$')


        # x-axis
        majorLocator   = MultipleLocator(0.25)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.25/2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

        tight_layout()
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')

        apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
                              
        steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
                                  markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')
#                               
#         nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
#                                 markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
#                               
        plt.legend(handles=[NFW, steidel, apparent],loc='lower right', labelspacing=1,
                                borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALT_corotate_vs_dist_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_ALLmodel'.format(only_close_velocities, include_nondetection, L_limit, min_separation)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')



if __name__ == '__main__':
    main()
    