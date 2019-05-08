#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id: SALT_paper_map.py, v3.0 04/30/18

Makes straight on-sky orientation plot as well as cylindrical and NFW model comparison
maps.



Previously: plot_onsky_absorber_vel.py, v2.0 04/04/18

Plot an impact parameter map showing the locations and velocities of each absorber wrt 
the galaxy (2/19/18)

v2: Orient so all the galaxies have approaching side on the left (04/04/18)

'''


import sys
import os
import csv
import time


from pylab import *
# import atpy
# from math import *
from utilities import *
from scipy import stats
import getpass
import math
import pickle
import json
import io
import numpy as np
import matplotlib.pyplot as plt
import correlateSingle11 as correlateSingle


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
  
#         vrot_vals = data['vrot_vals']
#         vrot_incCorrected_vals = data['vrot_incCorrected_vals']
        right_vrot_incCorrected_avg = data['right_vrot_incCorrected_avg']
        left_vrot_incCorrected_avg = data['left_vrot_incCorrected_avg']

        right_vrot_avg = data['right_vrot_avg']
        left_vrot_avg = data['left_vrot_avg']
        
#         xVals = data['xVals']
        inc = data['inclination']
        vsys_measured = data['vsys_measured']
#         galaxyName = data['name']
#         RA_galaxy = data['RAdeg']
#         Dec_galaxy = data['DEdeg']
        dist = data['dist']
        majDiam = data['majDiam']
        PA = data['PA']

        Lstar = data['Lstar']
        e_Lstar = data['e_Lstar']
        
#         print "Lstar: ", Lstar, "type(Lstar): ",type(Lstar)
    
#         agn = data['agn']
        
        return vsys_measured, right_vrot_avg, right_vrot_incCorrected_avg, left_vrot_avg, left_vrot_incCorrected_avg, inc, PA, dist, majDiam, Lstar, e_Lstar
        


def main():
    hubbleConstant = 71.0
    
    # only include absorbers that have dv less than or equal to the maximal rotation velocity?
    only_close_velocities = True
    
    # include open circles for sightlines with no absorption detected?
    include_nondetection = True
    
    # what range of Lstar systems to include?
#     Lstar_range = [0.0, 0.6]
#     Lstar_range = [0.60001, 100.0]
    Lstar_range = [0.0, 100.0]
#     Lstar_range = [0.0, 0.8]

    # azimuth limit for "maybe" trigger
    az_limit = 85.
    
    # size of legend symbols
    legend_size = 12
    
    # size of legend font
    legend_font = 12
    
    # minimum distance to another galaxy
    min_separation = False

    # how far to zoom in for zoom-in plot? Units of R_vir
    zoom_limit = 1.0
    
    # which plot to make?
    plot_b_comparison = True

    # use fits vs integrated values?
    use_fits = True
    
##########################################################################################
##########################################################################################
    # collect data
    directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/'
#     filename = '{0}salt_galaxy_sightlines_cut_plus_ancillary.csv'.format(directory)    
    filename = '{0}salt_galaxy_sightlines_cut_plus_ancillary_fits.csv'.format(directory)    

    theFile = open(filename,'rU')
    tableReader = csv.DictReader(theFile)

    # lists to populate
    nameList = []
    targetList = []
    combinedNameList = []
    vList = []
    wList = []
    RA_targetList = []
    Dec_targetList = []
    incList = []
    paList = []
    azList = []
    RvirList = []
    markerColorList = []
    VhelList = []
    markerColorList_model = []
    markerColorList_NFWmodel = []
    bList = []
    LstarList = []
    
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
    
    
    
    for t in tableReader:

        name = t['Name']
        target = t['Target']
        lyaV = eval(t['Lya_v'])
        lyaW = eval(t['Lya_W'])
        b = eval(t['b'])
        fit_b = eval(t['fit_b'])
        RA_galaxy = eval(t['RAdeg'])
        Dec_galaxy = eval(t['DEdeg'])
        RA_target = eval(t['RAdeg_target'])
        Dec_target = eval(t['DEdeg_target'])
        vHel = eval(t['Vhel'])
        
        PA_observed = eval(t['PA_observed'])
        PA_adjust = eval(t['PA_adjust'])
        
        model_range = eval(t['model_range'])
        NFW_range = eval(t['NFW_range'])
        
        gfilename = directory + 'rot_curves/' + name + '-summary4.json'
        vsys_measured, right_vrot_avg, right_vrot_incCorrected_avg, left_vrot_avg, left_vrot_incCorrected_avg, inc, PA, dist, majDiam, Lstar, e_Lstar = get_data(gfilename)
#         vsys_measured, right_vrot_avg, left_vrot_avg, inc, PA, dist, majDiam, Lstar, e_Lstar = get_data(gfilename)
        
        # remove inclination correction to get apparent velocity
        leftVel = left_vrot_incCorrected_avg * np.sin(inc * np.pi/180.)
        rightVel = right_vrot_incCorrected_avg * np.sin(inc * np.pi/180.)
        
        
        # calculate impact parameter and shit
        impact = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist)
    
        # RA component of impact parameter - by setting the Dec to be the same for both
        impact_RA = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_galaxy,dist)
    
        # Dec component of impact parameter - by setting the RA to be the same for both
        impact_Dec = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_galaxy,Dec_target,dist)
    
        # calculate azimuth
        az = calculateAzimuth(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist,PA)
        az2 = calculateAzimuth(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist,PA+10.)
    
        # calculate R_vir
        Rvir = calculateVirialRadius(majDiam)
        
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
        impact_RA_rot,impact_Dec_rot = np.array(coords.dot(R))[0]
#         print "R: ",R
#         print 'impact_RA_rot: ',impact_RA_rot
#         print 'impact_Dec_rot: ',impact_Dec_rot

                
        # scale to virial radius
        impact_rvir = impact/Rvir
        impact_RA_vir = impact_RA_rot/Rvir
        impact_Dec_vir = impact_Dec_rot/Rvir
                
        # compare to the absorber velocity:
        # negative means the Lya is higher velocity (red)
#         dv = vsys_measured - lyaV

        # switch it around actually, this matches the rotation curve (positive is going
        # away, or higher than systemic velocity gas)
        if lyaV != -99:
            dv = lyaV - vsys_measured
#             print 'lyaV, dv = ',lyaV,dv
        else:
            print 'else. dv = x'
            dv = 'x'
            

        # check on which 'side' of the galaxy the absorber is found
        if name == 'CGCG039-137':
             # regular - left side of rotation curve is 'left' on sky
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

        if name == 'ESO343-G014':
             # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'IC5325':
             # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'MCG-03-58-009':
             # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

        if name == 'NGC1566':
             # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC3513':
             # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC3633':
             # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC4536':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC4939':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC5364':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC5786':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
            print 'impact_RA_vir: ',impact_RA_vir
            print 'rot_vel: ',rot_vel
                
        if name == 'UGC09760':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
##########################################################################################
##########################################################################################
          
        if name == 'NGC3198':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC4565':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'UGC04238':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC3351':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC4529':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC6140':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

        if name == 'NGC5907':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

        if name == 'UGC06446':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

        if name == 'UGC06399':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC3726':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

        if name == 'NGC3067':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

        if name == 'NGC2770':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

        if name == 'NGC3432':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC3666':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC5951':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'NGC7817':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel
                
        if name == 'UGC08146':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

            
        # now compare to the rotation velocity
        color_yes = '#1b9e77'    # greenish
        color_no = '#d95f02'     # orange
#         color_maybe = '#7570b3'  # blue-ish purple
        color_maybe = 'grey'
        color_nonDetection = 'grey'
        markerColor = 'black'
        
        if isNumber(dv):
            # the absorption and rotation velocity match
            if (dv > 0 and rot_vel > 0) or (dv < 0 and rot_vel < 0):
                markerColor = color_yes
            
            # mismatch
            elif (dv < 0 and rot_vel > 0) or (dv > 0 and rot_vel < 0):
                markerColor = color_no
            
            else:
                markerColor = 'black'

            
            # compare to the models
            error = 10.
            model_answer = withinRange(dv, model_range, error)
            NFW_model_answer = withinRange(dv, NFW_range, error)
        
            if model_answer:
                markerColor_model = color_yes
            else:
                markerColor_model = color_no
            
            if NFW_model_answer:
                markerColor_NFWmodel = color_yes
            else:
                markerColor_NFWmodel = color_no
            
            if az >= az_limit:
                markerColor_model = color_maybe
                markerColor_NFWmodel = color_maybe

        else:
            print 'dv == x :', dv, name
            print
            markerColor = color_nonDetection
        
        
        # if too close to the minor axis, call it uncertain/maybe
        if az >= az_limit:
            markerColor = color_maybe
            
        if name == 'NGC3513':
            markerColor = color_maybe


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = False

#         combinedName = '${\rm '+name + '-' + target+'}$'.format(name,target)
        combinedName = r'$\rm {0} : {1}$'.format(name,target)
            
        print '{0} - dv={1} vs model={2} => {3}'.format(combinedName, dv, model_range, model_answer)
        print

            
        # decide some things
        if withinRange(Lstar, Lstar_range, 0.0):
            add_to_list = True
        else:
            add_to_list = False
        
        separation = 500.
        if min_separation and add_to_list:
            # correlate with environment
            agnSeparation = False
            minVcorr = False
            minSize = False
            correlation = correlateSingle.correlateTarget(name, min_separation, agnSeparation, minVcorr, minSize, slow=False, searchAll=True)
            galaxyInfo = correlation[name]
            
            print 'galaxyInfo: ',galaxyInfo
                
            for row in galaxyInfo:
                vhel, galaxyRow = row
                separation = galaxyRow['impactParameter (kpc)']
                galaxyVel = galaxyRow['radialVelocity (km/s)']
                
                print 'separation: ',separation
                print 'galaxyVel: ',galaxyVel
                print 'vHel: ',vHel
                
                print 'withinRange(galaxyVel, [vHel-400, vHel+400], 0.0): ',withinRange(galaxyVel, [vHel-400, vHel+400], 0.0)
                print

                if withinRange(galaxyVel, [vHel-400, vHel+400], 0.0) and add_to_list:
                    if separation <= min_separation and separation >0.0:
                        add_to_list = False
                        print 'False for {0} - {1}'.format(name, separation)
                        print
                        break
        
        
        # populate the lists
        if add_to_list:
            # first detections
            if isNumber(dv) and only_close_velocities:
                if abs(dv) <= (abs(rot_vel) + 10):

                    nameList.append(name)
                    targetList.append(target)
                    vList.append(lyaV)
                    wList.append(lyaW)
                    RA_targetList.append(impact_RA_vir)
                    Dec_targetList.append(impact_Dec_vir)
                    incList.append(inc)
                    paList.append(PA)
                    azList.append(az)
                    RvirList.append(Rvir)
                    markerColorList.append(markerColor)
                    combinedNameList.append(combinedName)
                    markerColorList_NFWmodel.append(markerColor_NFWmodel)
                    markerColorList_model.append(markerColor_model)
                    
                    if use_fits:
                        bList.append(fit_b)
                    else:
                        bList.append(b)
                        
                    LstarList.append(Lstar)
                
                else:
                    print 'too far: ',name,' , dv = ',dv, ' vs rot_vel = ',rot_vel
                    print
                
            else:
                if isNumber(dv):
                    nameList.append(name)
                    targetList.append(target)
                    vList.append(lyaV)
                    wList.append(lyaW)
                    RA_targetList.append(impact_RA_vir)
                    Dec_targetList.append(impact_Dec_vir)
                    incList.append(inc)
                    paList.append(PA)
                    azList.append(az)
                    RvirList.append(Rvir)
                    markerColorList.append(markerColor)
                    combinedNameList.append(combinedName)
                    markerColorList_NFWmodel.append(markerColor_NFWmodel)
                    markerColorList_model.append(markerColor_model)
                                  
                    if use_fits:
                        bList.append(fit_b)
                    else:
                        bList.append(b)

                    LstarList.append(Lstar)

                
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
            print 'outside range: {0} - {1} Lstar, {2} = separation'.format(Lstar, Lstar_range,separation)
            print

            
#         print 'added'
#         print 'impact_RA_vir: ',impact_RA_vir
#         print 'impact_Dec_vir: ',impact_Dec_vir
#         print
#         print
        
##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
# Plot b values for co-rotators vs anti-rotators
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

        alpha = 0.6
        
        L_limit = 0.6
        

        corotate_b = []
        antirotate_b = []
        
        corotate_b_close = []
        antirotate_b_close = []
        
        corotate_b_far = []
        antirotate_b_far = []
        
        Lstar_high = []
        Lstar_low = []
        
        for m, b, ra, dec, L in zip(markerColorList_NFWmodel, bList, RA_targetList, Dec_targetList, LstarList):
            if m == color_yes:
                corotate_b.append(b)
                
                if np.sqrt(ra**2 + dec**2) <= zoom_limit:
                    corotate_b_close.append(b)
                    
                else:
                    corotate_b_far.append(b)
            
            if m == color_no:
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
        
        hist(corotate_b, bins=bins, histtype='bar', lw=1.5, color = 'blue', alpha=alpha, label=r'$\rm Co-rotators$')
        hist(antirotate_b, bins=bins, histtype='bar', lw=1.5, color = 'red', alpha=alpha, label=r'$\rm Anti-rotators$')
        
        ylim(0, 14)
        legend(scatterpoints=1,prop={'size':12},loc=2,fancybox=True)
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
        
        hist(corotate_b_close, bins=bins, histtype='bar', lw=1.5, color = 'blue', alpha=alpha, label=r'$\rm Co-rotators~(\rho \leq {0} R_{{vir}})$'.format(zoom_limit))
        hist(antirotate_b_close, bins=bins, histtype='bar', lw=1.5, color = 'red', alpha=alpha, label=r'$\rm Anti-rotators~(\rho \leq {0}R_{{vir}})$'.format(zoom_limit))
        
        ylim(0, 4)
        legend(scatterpoints=1,prop={'size':12},loc=2,fancybox=True)
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
        
        hist(Lstar_low, bins=bins, histtype='bar', lw=1.5, color = 'blue', alpha=alpha, label=r'$\rm Lstar \leq {0})$'.format(L_limit))
        hist(Lstar_high, bins=bins, histtype='bar', lw=1.5, color = 'red', alpha=alpha, label=r'$\rm Lstar > {0})$'.format(L_limit))
        
#         ylim(0, 3)
        legend(scatterpoints=1,prop={'size':12},loc=2,fancybox=True)
        xlabel(r'$\rm b ~ [km~s^{-1}]$')
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
        
        directory = '/Users/frenchd/Research/test/'
        save_name = 'SALT_bhist_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_fits_{4}'.format(only_close_velocities, include_nondetection, L_limit, min_separation, use_fits)
        savefig("{0}/{1}.pdf".format(directory, save_name),bbox_inches='tight')
            
        
        # write out a file breaking down all this shit
        stats_filename = '{0}/{1}_stats.txt'.format(directory, save_name)
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

    
if __name__ == '__main__':
    main()
    