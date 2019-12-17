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
    plot_model_visualizer = True


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
    
    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red
    
    
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
    steidel_rangeList = []
    NFW_rangeList = []
    
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
                    steidel_rangeList.append(steidel_range)
                    NFW_rangeList.append(NFW_range)
                    
                    vmaxList.append(rot_vel)

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
                    steidel_rangeList.append(steidel_range)
                    NFW_rangeList.append(NFW_range)
                    
                    vmaxList.append(rot_vel)

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
##########################################################################################
# Plot model result ranges and dv values for each system as a function of impact parameter
#
#
##########################################################################################
##########################################################################################


    if plot_model_visualizer:
        # initial figure
        fig = plt.figure(figsize=(10.0,5.0))
        subplots_adjust(hspace=0.500)

        ax = fig.add_subplot(211)

        alpha_no = 0.55
        alpha_yes = 0.65

        marker_lw = 2.0
        
#         bin_edges = equal_hist_edges(impactList, 3)
#         bin_edges = [0, 75, 200, 600.]
#         right_edge = 2.5
#         impact_bins, impact_bin_edges = np.histogram(impactList, bins=bin_edges)
        
        x_data = []
        y_apparent = []
        y_NFW = []
        y_steidel = []
        
        y_apparent_num = []
        y_NFW_num = []
        y_steidel_num = []
        
        color_apparent = 'grey'
        color_NFW = color_blue
        color_steidel = color_red
        
        symbol_lya = 'D'
        color_lya = 'black'
        alpha_lya = 1.0
        marker_size = 0
        
        max_nfw = max(np.ravel(np.array(NFW_rangeList)))
        min_nfw = min(np.ravel(np.array(NFW_rangeList)))

        max_steidel = max(np.ravel(np.array(steidel_rangeList)))
        min_steidel = min(np.ravel(np.array(steidel_rangeList)))

        for a, a_nfw, nfw, a_steidel, steidel, dv, e_dv, imp, name, inc, az in zip(apparent_answers, NFW_answers, NFW_rangeList, steidel_answers, steidel_rangeList, dvList, e_dv_list, impactList, nameList, incList, azList):
            # apparent:
#             plot([imp, imp],
#                 a,
#                 lw=marker_lw,
#                 marker='_',
#                 ls='-',
#                 color=color_apparent,
#                 ms=marker_size)


            if steidel[0] <0:
                steidel_low = steidel[0]/min_steidel
                steidel_high = steidel[1]/min_steidel
            else:
                steidel_low = steidel[0]/max_steidel
                steidel_high = steidel[1]/max_steidel
                
            dv_norm = dv/max_steidel
            e_dv_norm = e_dv/max_steidel
            
            if az < 100:
                print('nfw: ',nfw)
                print()
            
                # NFW
#                 ax.plot([imp-2, imp-2],
#                     nfw,
#                     lw=3.0,
#                     marker='_',
#                     ls='-',
#                     color=color_NFW,
#                     ms=marker_size)

                # Steidel
                ax.plot([imp+2, imp+2],
                    [steidel_low, steidel_high],
                    lw=3.0,
                    marker='_',
                    ls='-',
                    color=color_steidel,
                    ms=marker_size)
        
                ax.errorbar(imp,
                        dv_norm,
                        yerr=e_dv_norm,
                        lw=1,
                        elinewidth = 0.8,
                        marker=symbol_lya,
                        markeredgewidth=0.8,
                        markeredgecolor='black',
                        ms = 6,
                        color=color_lya,
                        alpha=alpha_lya)

        ax2 = fig.add_subplot(212)

        for a, a_nfw, nfw, a_steidel, steidel, dv, e_dv, imp, name, inc, az in zip(apparent_answers, NFW_answers, NFW_rangeList, steidel_answers, steidel_rangeList, dvList, e_dv_list, impactList, nameList, incList, azList):
            # apparent:
#             plot([imp, imp],
#                 a,
#                 lw=marker_lw,
#                 marker='_',
#                 ls='-',
#                 color=color_apparent,
#                 ms=marker_size)
            
            if nfw[0] <0:
                nfw_low = nfw[0]/min_nfw
                nfw_high = nfw[1]/min_nfw
            else:
                nfw_low = nfw[0]/max_nfw
                nfw_high = nfw[1]/max_nfw
                
            dv_norm = dv/max_nfw
            e_dv_norm = e_dv/max_nfw
            
            if az >= 0:
                print('nfw: ',nfw)
                print()
            
                # NFW
                ax2.plot([imp-2, imp-2],
                    [nfw_low, nfw_high],
                    lw=3.0,
                    marker='_',
                    ls='-',
                    color=color_NFW,
                    ms=marker_size)

                # Steidel
#                 ax2.plot([imp+2, imp+2],
#                     steidel,
#                     lw=3.0,
#                     marker='_',
#                     ls='-',
#                     color=color_steidel,
#                     ms=marker_size)
        
                ax2.errorbar(imp,
                        dv_norm,
                        yerr=e_dv_norm,
                        lw=1,
                        elinewidth = 0.8,
                        marker=symbol_lya,
                        markeredgewidth=0.8,
                        markeredgecolor='black',
                        ms = 6,
                        color=color_lya,
                        alpha=alpha_lya)


        
#         ylim(0, 1.05)
#         xlim(0, 350.0)
        
        xlabel(r'$\rm Impact~Parameter~[kpc]$')
        ylabel(r'$\rm Velocity~[km~s^{-1}]$')


        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax2.xaxis.set_major_locator(majorLocator)
        ax2.xaxis.set_major_formatter(majorFormatter)
        ax2.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax2.yaxis.set_major_locator(majorLocator)
        ax2.yaxis.set_major_formatter(majorFormatter)
        ax2.yaxis.set_minor_locator(minorLocator)
        
        
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax2.yaxis.set_ticks_position('both')
        ax2.xaxis.set_ticks_position('both')

        ax.set_xlim(0, 500.0)
        ax2.set_xlim(0, 500.0)

        
        if plotGrid:
            ax.grid(True, alpha=0.5)

        tight_layout()
        
#         import matplotlib.patches as mpatches
#         import matplotlib.lines as mlines
# 
#         NFW = mlines.Line2D([], [], color=color_yes, marker=marker_NFW, lw=marker_lw, ls=ls_NFW,
#                                   markersize=marker_size, markeredgecolor='black', label=r'$\rm NFW ~Model $')
# 
#         apparent = mlines.Line2D([], [], color=color_maybe, marker=marker_apparent, lw=marker_lw, ls=ls_apparent,
#                                   markersize=marker_size, markeredgecolor='black', label=r'$\rm Apparent $')
#                               
#         steidel = mlines.Line2D([], [], color=color_no, marker=marker_steidel, lw=marker_lw, ls=ls_steidel,
#                                   markersize=marker_size, markeredgecolor='black', label=r'$\rm Steidel~ Model $')
# 
#         plt.legend(handles=[NFW, steidel, apparent],loc='lower left', labelspacing=1,
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)

##########################################################################################
        
        save_name = 'Model_visualizer_velstrict_{0}_nondet_{1}_minsep_{2}_neighbors_{3}'.format(only_close_velocities, include_nondetection, min_separation, neighbor_limit)
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')


if __name__ == '__main__':
    main()
    