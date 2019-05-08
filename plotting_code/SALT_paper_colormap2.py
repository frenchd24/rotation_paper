#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  $Id: SALT_paper_colormap.py, v3.0 04/30/18

Makes colormap showing the absorption-galaxy velocity difference

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
    
    # where to write to?
    out_directory = '/Users/frenchd/Research/test/SALT_maps_yes_maybe5/'
#     out_directory = '/Users/frenchd/Research/test/SALT_maps_yes/'

    # only include absorbers that have dv less than or equal to the maximal rotation velocity?
    only_close_velocities = False
    
    # include open circles for sightlines with no absorption detected?
    include_nondetection = True
    
    # what range of Lstar systems to include?
#     Lstar_range = [0.0, 0.6]
#     Lstar_range = [0.60001, 100.0]
    Lstar_range = [0.0, 100.0]
#     Lstar_range = [0.0, 0.5]
#     Lstar_range = [0.50001, 100.0]


    # azimuth limit for "maybe" trigger
    az_limit = 85.
    
    # size of legend symbols
    legend_size = 8
    
    # size of legend font
    legend_font = 8
    
    # minimum distance to another galaxy
    min_separation = False

    # how far to zoom in for zoom-in plot? Units of R_vir
    zoom_limit = 1.0
    
    # remove systems of lower inclination than this:
#     inc_limit = 75.
    inc_limit = 0.
    
    # include a +/- 10 km/s error in apparent results also?
    use_apparent_errors = True

    # which plot to make?
    sequential_colormap = True
    sequential_NFW_colormap = True
    colormap_colormap = True
        
    # include tags to include
    include_tags = ['yes','maybe']
#     include_tags = ['yes']

    # colormap
    colmap = cm.RdBu_r
    
    # limit in velocity difference between absorber and galaxy
    vel_limit = 400.
    
    # system number tag font size
    sysNum_size = 8

##########################################################################################
##########################################################################################
    # collect data
    directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/'
#     filename = '{0}salt_galaxy_sightlines_cut.csv'.format(directory)
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
    nameDict = {}
    dvList = []
    NFW_ranges = []
    
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
        RA_galaxy = eval(t['RAdeg'])
        Dec_galaxy = eval(t['DEdeg'])
        RA_target = eval(t['RAdeg_target'])
        Dec_target = eval(t['DEdeg_target'])
        vHel = eval(t['Vhel'])
        
        include_tag = t['include']
        sysNumber = t['number']
        
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

        if name == 'NGC3631':
            # regular
            if impact_RA_vir > 0:
                rot_vel = leftVel
            else:
                rot_vel = rightVel

            
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

        error = 10.
        rot_error = 10
        
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

            dv_up = dv + error
            dv_down = dv - error
            
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
            markerColor_model = color_maybe
            markerColor_NFWmodel = color_maybe
            
        if name == 'NGC3513':
            markerColor = color_maybe
            markerColor_model = color_maybe
            markerColor_NFWmodel = color_maybe
            

        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = False

#         combinedName = '${\rm '+name + '-' + target+'}$'.format(name,target)
        combinedName = r'$\rm {0} : {1}$'.format(name,target)
        
        nameDict[combinedName] = sysNumber
        
        print '{0} - dv={1} vs model={2} => {3}'.format(combinedName, dv, model_range, model_answer)
        print

            
        # decide some things
        if withinRange(Lstar, Lstar_range, 0.0):
            add_to_list = True
        else:
            add_to_list = False
            
        if inc <= inc_limit:
            add_to_list = False
            
        if include_tag not in include_tags:
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
                if abs(dv) <= (abs(rot_vel) + rot_error):

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
                    dvList.append(dv)
                    
                    NFW_ranges.append(NFW_range)
                
                
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
                    dvList.append(dv)

                    NFW_ranges.append(NFW_range)

                
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


        
##########################################################################################
##########################################################################################
##########################################################################################
    # sort in order of largest to smallest equivalent width
    orderedList = []
    for ra,dec,c,w,name,model,NFW,dv,NFW_range in zip(RA_targetList,Dec_targetList,markerColorList, wList,combinedNameList, markerColorList_model, markerColorList_NFWmodel,dvList, NFW_ranges):
        orderedList.append([w,[ra, dec, c, name, model, NFW, dv, NFW_range]])
        
    orderedList.sort(reverse=True)
    RA_targetList2 = []
    Dec_targetList2 = []
    markerColorList2 = []
    wList2 = []
    combinedNameList2 = []
    markerColorList_NFWmodel2 = []
    markerColorList_model2 = []
    countList = []
    dvList2 = []
    
    NFW_ranges2 = []
    
    count = 1
    count_dict = {}
    for i in orderedList:
        w, rest = i
        ra, dec, c, name, model, NFW, dv, NFW_range = rest
        
        print
        print 'dv: {0} - {1}'.format(dv, name)
        print
        
        RA_targetList2.append(ra)
        Dec_targetList2.append(dec)
        markerColorList2.append(c)
        wList2.append(w)
        combinedNameList2.append(name)
        markerColorList_NFWmodel2.append(NFW)
        markerColorList_model2.append(model)
        dvList2.append(dv)
        
        NFW_ranges2.append(NFW_range)
        
        countList.append(nameDict[name])
        
    countList_non = []
    for name in combinedNameList_non:
        countList_non.append(nameDict[name])

##########################################################################################
##########################################################################################
# Color based on velocity difference, gradient scale
#
#
##########################################################################################
##########################################################################################

    if colormap_colormap:
        # initial figure
        fig = plt.figure(figsize=(10,8))
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
    
    #     r = 4.0
    #     ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
        ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

        # colormap stuff
        vmin_val = -400
        vmax_val = +400

    #     vmin_val = min(dvList2)
    #     vmax_val = max(dvList2)
        norm = matplotlib.colors.Normalize(vmin = vmin_val, vmax = vmax_val)
        m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
        # plot the rest
        largestEW = max(wList2)
        smallestEW = min(wList2)
        maxSize = 500
        minSize = 30
    
        # normalize the marker size to reflect EW
        newSizeList = []
        for w in wList2:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)


        # make different style markers for different colors
        for i in arange(len(markerColorList_NFWmodel2)):
            marker = 'D'
            marker_lw = 0.2

    #         if markerColorList_NFWmodel2[i] == color_maybe:
    #             marker = 'o'
    #         if markerColorList_NFWmodel2[i] == color_no:
    #             marker = 'x'
    #             marker_lw = 1.5
    #         if markerColorList_NFWmodel2[i] == color_yes:
    #             marker = 'D'
    
            plot1 = ax.scatter(RA_targetList2[i], 
                                Dec_targetList2[i], 
                                s = newSizeList[i],
                                c = dvList2[i],
                                vmin = vmin_val,
                                vmax = vmax_val,
                                marker = marker,
                                lw = marker_lw,
                                edgecolor = 'black',
                                cmap = colmap)

    #     plot1 = ax.scatter(RA_targetList2, Dec_targetList2, s=0, c=dvList2, 
    #                     vmin=vmin_val, vmax=vmax_val, marker='.', cmap=colmap)

        # now make the colorbar work
        step = 100
        ticks = arange(vmin_val,vmax_val+step, int(step))
        cbar = plt.colorbar(plot1, ticks=ticks, format=r'$\rm %d$', cmap=colmap, orientation='vertical')
        cbar.set_label(r'$\rm \Delta v ~[km ~s^{-1}]$')

        # put labels on each point
        xTagOffset = 2.0
        yTagOffset = 1.
        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
    #                 
    #         annotate(countList[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
    #         xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)

            tag = nameDict[combinedNameList2[i]]
            yTagOffset = 5.0 + (newSizeList[i]/50.)
        
            if not previousNames.has_key(tag):
                annotate(tag,xy=(RA_targetList2[i], Dec_targetList2[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=sysNum_size)
            
                previousNames[tag] = 1

    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(combinedNameList_non)):
            ax.plot(RA_targetList_non[i], Dec_targetList_non[i], color=markerColorList_non[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

    #         yTagOffset = 5.0
    #         annotate(countList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
    #         xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        
            tag = nameDict[combinedNameList_non[i]]
            yTagOffset = 5.0 + (newSizeList[i]/50.)
        
            if not previousNames.has_key(tag):
                annotate(tag,xy=(RA_targetList2[i], Dec_targetList2[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=sysNum_size)
            
                previousNames[tag] = 1
        

    ##########################################################################################

        xlabel(r'$\rm x ~[R_{vir}]$')
        ylabel(r'$\rm y ~[R_{vir}]$')
    
        ax.set_xlim(-3.0, 3.0)
        ax.set_ylim(-3.0, 3.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(2.97, 0.06),\
        xytext=(0.0,0.0),textcoords='offset points',size=9)


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


        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

    #     corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
    #                               markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')
    # 
    #     maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
    #                               markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
    #                               
    #     antirotate = mlines.Line2D([], [], color=color_no, marker='x',lw=0,
    #                               markersize=legend_size, markeredgecolor=color_no, label=r'$\rm Anti-rotation$')
    #                               
    #     nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
    #                             markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
    #                               
    #     plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
    #                             borderpad=0.8, fontsize=legend_font, fancybox=True)



    ##########################################################################################
        
        save_name = 'SALT_colormap_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation)

    #     savefig("{0}SALT_map1.pdf".format(directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        summary_filename = '{0}/{1}_summary.txt'.format(out_directory, save_name)
        summary_file = open(summary_filename,'wt')
    

        for num, name in zip(countList, combinedNameList2):
            summary_file.write('{0}. {1}, \n'.format(num,name))
        
        for num, name in zip(countList_non, combinedNameList_non):
            summary_file.write('{0}. {1}, \n'.format(num,name))
    
        summary_file.close()




##########################################################################################
##########################################################################################
# Color based on velocity difference, sequential scale (discrete color steps)
#
#
##########################################################################################
##########################################################################################

    if sequential_colormap:
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
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
        ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

        # blue colors
        b1 = '#eff3ff'
        b2 = '#bdd7e7'
        b3 = '#6baed6'
        b4 = '#3182bd'
        b5 = '#08519c'
        
        # red colors
        r1 = '#fee5d9'
        r2 = '#fcae91'
        r3 = '#fb6a4a'
        r4 = '#de2d26'
        r5 = '#a50f15'
        
        # blue velocity bins
        bv1 = -50
        bv2 = -100
        bv3 = -150
        bv4 = -200
        bv5 = -400
        
        # red velocity bins
        rv1 = 50
        rv2 = 100
        rv3 = 150
        rv4 = 200
        rv5 = 400
    
        # plot the rest
        largestEW = max(wList2)
        smallestEW = min(wList2)
        maxSize = 500
        minSize = 30

        # normalize the marker size to reflect EW
        newSizeList = []
        for w in wList2:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)


        for i in arange(len(dvList2)):
            marker = 'D'
            marker_lw = 0.2
            
            dv = float(dvList2[i])
            print 'dv: ',dv
            print '[0,bv1]: ',[0,bv1]
            print 'withinRange(dv,[0,bv1],0): ',withinRange(dv,[0,bv1],0)
            print
            print
            
            if withinRange(dv,[bv1,0],0):
                mc = b1
                
            elif withinRange(dv,[bv2,bv1],0):
                mc = b2
                
            elif withinRange(dv,[bv3,bv2],0):
                mc = b3
                
            elif withinRange(dv,[bv4,bv3],0):
                mc = b4
                
            elif withinRange(dv,[bv5,bv4],0):
                mc = b5
                
            elif withinRange(dv,[0,rv1],0):
                mc = r1
                
            elif withinRange(dv,[rv1,rv2],0):
                mc = r2
                
            elif withinRange(dv,[rv2,rv3],0):
                mc = r3
                
            elif withinRange(dv,[rv3,rv4],0):
                mc = r4

            elif withinRange(dv,[rv4,rv5],0):
                mc = r5
                
            else:
                mc = 'black'


            plot1 = ax.scatter(RA_targetList2[i], 
                                Dec_targetList2[i], 
                                s = newSizeList[i],
                                c = mc,
                                marker = marker,
                                lw = marker_lw,
                                edgecolor = 'black')


        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)
                
            tag = nameDict[combinedNameList2[i]]        
            if not previousNames.has_key(tag):
                annotate(tag,xy=(RA_targetList2[i], Dec_targetList2[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=sysNum_size)
            
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
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=sysNum_size)
        
            tag = nameDict[combinedNameList_non[i]]
        
            if not previousNames.has_key(tag):
                annotate(tag,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=sysNum_size)
            
                previousNames[tag] = 1

    ##########################################################################################

        xlabel(r'$\rm x ~[R_{vir}]$')
        ylabel(r'$\rm y ~[R_{vir}]$')
    
        ax.set_xlim(-3.0, 3.0)
        ax.set_ylim(-3.0, 3.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(2.97, 0.06),\
        xytext=(0.0,0.0),textcoords='offset points',size=9)


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


        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        bv1_leg = mlines.Line2D([], [], color=b1, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(0, bv1))
    
        bv2_leg = mlines.Line2D([], [], color=b2, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(bv1,bv2))
                                  
        bv3_leg = mlines.Line2D([], [], color=b3, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(bv2,bv3))
                                  
        bv4_leg = mlines.Line2D([], [], color=b4, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(bv3,bv4))
                                  
        bv5_leg = mlines.Line2D([], [], color=b5, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(bv4,bv5))
                                  
        rv1_leg = mlines.Line2D([], [], color=r1, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(rv1,0))
                                  
        rv2_leg = mlines.Line2D([], [], color=r2, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(rv2,rv1))

        rv3_leg = mlines.Line2D([], [], color=r3, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(rv3,rv2))

        rv4_leg = mlines.Line2D([], [], color=r4, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(rv4,rv3))

        rv5_leg = mlines.Line2D([], [], color=r5, marker=marker,lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm [{0}, {1}]$'.format(rv5,rv4))

        plt.legend(handles=[rv5_leg, rv4_leg, rv3_leg, rv2_leg, rv1_leg, bv1_leg, bv2_leg, bv3_leg, bv4_leg, bv5_leg],
                    loc='lower right', borderpad=0.8, fontsize=legend_font, fancybox=True)

    ##########################################################################################
        
        save_name = 'SALT_seq_colormap_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation)

    #     savefig("{0}SALT_map1.pdf".format(directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        summary_filename = '{0}/{1}_summary.txt'.format(out_directory, save_name)
        summary_file = open(summary_filename,'wt')


        for num, name in zip(countList, combinedNameList2):
            summary_file.write('{0}. {1}, \n'.format(num,name))
        
        for num, name in zip(countList_non, combinedNameList_non):
            summary_file.write('{0}. {1}, \n'.format(num,name))
    
        summary_file.close()
        
        

##########################################################################################
##########################################################################################
# NFW model - Color based on velocity difference, sequential scale (discrete color steps)
# 
# 
##########################################################################################
##########################################################################################

    if sequential_NFW_colormap:
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
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
        ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################
    
        # colors should be:
        # dark blue crosses - below NFW range by >10 km/s
        # blue circles - within 10 km/s of the lower edge of NFW scale
        # green diamonds - anywhere within the accepted range
        # red circles - within 10 km/s of the upper edge of the NFW scale
        # dark red crosses - above NFW range by >10 km/s


        # blue colors
        b1 = '#eff3ff' # vaguely blue-ish
        b2 = '#bdd7e7' # steely-grey light blue
        b3 = '#6baed6' # sky blue
        b4 = '#3182bd' # nice deep, watery blue
        b5 = '#08519c' # dark navy
        
        # red colors
        r1 = '#fee5d9' # vaguely red-ish
        r2 = '#fcae91' # light salmony red
        r3 = '#fb6a4a' # bad sunburn-red 
        r4 = '#de2d26' # stoplight red
        r5 = '#a50f15' # deep red
        
        g1 = '#31a354' # grass-green
        g2 = '#006837' # deep forest-green
        
        # blue velocity bins
        bv1 = -50
        bv2 = -100
        bv3 = -150
        bv4 = -200
        bv5 = -400
        
        # red velocity bins
        rv1 = 50
        rv2 = 100
        rv3 = 150
        rv4 = 200
        rv5 = 400
    
        # plot the rest
        largestEW = max(wList2)
        smallestEW = min(wList2)
        maxSize = 500
        minSize = 30

        # normalize the marker size to reflect EW
        newSizeList = []
        for w in wList2:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)


        for i in arange(len(dvList2)):
            marker = 'D'
            marker_lw = 0.2
            
            dv = float(dvList2[i])
            NFW_range = NFW_ranges2[i]
            
            low_edge = NFW_range[0]
            high_edge = NFW_range[1]
            
            
            
            
            print 'dv: ',dv
            print '[0,bv1]: ',[0,bv1]
            print 'withinRange(dv,[0,bv1],0): ',withinRange(dv,[0,bv1],0)
            print
            print
            
            if withinRange(dv, NFW_range, 0):
                mc = g1
                marker = 'D'
            
            elif withinRange(dv, [low_edge-10, low_edge+10], 0):
                mc = b3
                marker = 'o'
                
            elif withinRange(dv, [high_edge-10, high_edge+10], 0):
                mc = r3
                marker = 'o'
            
            elif dv < low_edge-10.:
                mc = b5
                marker = 'X'

            elif dv > high_edge+10.:
                mc = r5
                marker = 'X'
            
            else:
                print 'wtf? - ',dv, NFW_range
                mc = 'grey'
                marker = 'd'

            plot1 = ax.scatter(RA_targetList2[i], 
                                Dec_targetList2[i], 
                                s = newSizeList[i],
                                c = mc,
                                marker = marker,
                                lw = marker_lw,
                                edgecolor = 'black')

        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)
                
            tag = nameDict[combinedNameList2[i]]        
            if not previousNames.has_key(tag):
                annotate(tag,xy=(RA_targetList2[i], Dec_targetList2[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=sysNum_size)
            
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
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=sysNum_size)
        
            tag = nameDict[combinedNameList_non[i]]
        
            if not previousNames.has_key(tag):
                annotate(tag,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
                xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=sysNum_size)
            
                previousNames[tag] = 1

    ##########################################################################################

        xlabel(r'$\rm x ~[R_{vir}]$')
        ylabel(r'$\rm y ~[R_{vir}]$')
    
        ax.set_xlim(-3.0, 3.0)
        ax.set_ylim(-3.0, 3.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(2.97, 0.06),\
        xytext=(0.0,0.0),textcoords='offset points',size=9)


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


        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines

        below_leg = mlines.Line2D([], [], color=b5, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm \Delta v < NFW~range~\pm 10~[km s^{{-1}}]$')
    
        blue_edge_leg = mlines.Line2D([], [], color=b3, marker='o',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm \Delta v ~within~10~km s^{{-1}}~of~NFW~lower~edge$')
                                  
        within_leg = mlines.Line2D([], [], color=g1, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm \Delta v ~within~NFW~range$')
                                  
        red_edge_leg = mlines.Line2D([], [], color=r3, marker='o',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm \Delta v ~within~10~km s^{{-1}}~of~NFW~upper~edge$')
                                  
        above_leg = mlines.Line2D([], [], color=r5, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm \Delta v > NFW~range~\pm 10~[km s^{{-1}}]$')

        plt.legend(handles=[below_leg, blue_edge_leg, within_leg, red_edge_leg, above_leg],
                    loc='lower right', borderpad=0.8, fontsize=legend_font, fancybox=True)


    ##########################################################################################
        
        save_name = 'SALT_seq_NFW_colormap_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation)

    #     savefig("{0}SALT_map1.pdf".format(directory),dpi=400,bbox_inches='tight')
        savefig("{0}/{1}.pdf".format(out_directory, save_name),bbox_inches='tight')

        summary_filename = '{0}/{1}_summary.txt'.format(out_directory, save_name)
        summary_file = open(summary_filename,'wt')
    

        for num, name in zip(countList, combinedNameList2):
            summary_file.write('{0}. {1}, \n'.format(num,name))
        
        for num, name in zip(countList_non, combinedNameList_non):
            summary_file.write('{0}. {1}, \n'.format(num,name))
    
        summary_file.close()


    
if __name__ == '__main__':
    main()