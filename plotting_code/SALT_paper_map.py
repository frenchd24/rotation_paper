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
    
    # where to write to?
    out_directory = '/Users/frenchd/Research/test/SALT_maps/'
    
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
    legend_size = 12
    
    # size of legend font
    legend_font = 12
    
    # minimum distance to another galaxy
    min_separation = False

    # how far to zoom in for zoom-in plot? Units of R_vir
    zoom_limit = 1.0
    
    # remove systems of lower inclination than this:
#     inc_limit = 75.
    inc_limit = 0.

    # which plot to make?
    plot_onsky = True
    plot_cyl = True
    plot_NFW = True
    plot_zoom_in = True
    
    # include tags to include
    include_tags = ['yes','maybe']

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
        color_nonDetection = 'black'
        markerColor = 'black'
        
        
#         color_yes = '#0652ff'    # blue
#         color_no = '#ff9408'     # orange
        
        color_yes = '#1805db'    # ultramarine blue
        color_yes = '#436bad'      # french blue
        color_no = '#ec2d01'     # tomato red

        
        
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
    # sort in order of largest to smallest equivalent width
    orderedList = []
    for ra,dec,c,w,name,model,NFW in zip(RA_targetList,Dec_targetList,markerColorList, wList,combinedNameList, markerColorList_model, markerColorList_NFWmodel):
        orderedList.append([w,[ra, dec, c, name, model, NFW]])
        
    orderedList.sort(reverse=True)
    
    
    RA_targetList2 = []
    Dec_targetList2 = []
    markerColorList2 = []
    wList2 = []
    combinedNameList2 = []
    markerColorList_NFWmodel2 = []
    markerColorList_model2 = []
    countList = []
    
    # zoomed in
    RA_targetList_zoom = []
    Dec_targetList_zoom = []
    markerColorList_zoom = []
    wList_zoom = []
    combinedNameList_zoom = []
    markerColorList_NFWmodel_zoom = []
    markerColorList_model_zoom = []
    countList_zoom = []
    
    
    count = 1
    count_zoom = 1
    count_dict = {}
    count_dict_zoom = {}
    for i in orderedList:
        w, rest = i
        ra, dec, c, name, model, NFW = rest
        
        RA_targetList2.append(ra)
        Dec_targetList2.append(dec)
        markerColorList2.append(c)
        wList2.append(w)
        combinedNameList2.append(name)
        markerColorList_NFWmodel2.append(NFW)
        markerColorList_model2.append(model)
        
        # check if this galaxy-QSO pair already has a number
        if count_dict.has_key(name):
            system_count = count_dict[name]
#             countList.append(system_count)
            countList.append('')
        
        else:
            count_dict[name] = count
            countList.append(count)
            count +=1
            
            
        # separate, 'zoomed-in' set
        if math.sqrt(ra**2 + dec**2) <= zoom_limit:
            RA_targetList_zoom.append(ra)
            Dec_targetList_zoom.append(dec)
            markerColorList_zoom.append(c)
            wList_zoom.append(w)
            combinedNameList_zoom.append(name)
            markerColorList_NFWmodel_zoom.append(NFW)
            markerColorList_model_zoom.append(model)
        
            # check if this galaxy-QSO pair already has a number
            if count_dict_zoom.has_key(name):
                system_count = count_dict_zoom[name]
                countList_zoom.append('')
        
            else:
                count_dict_zoom[name] = count_zoom
                countList_zoom.append(count_zoom)
                count_zoom +=1
        
    countList_non = []
    countList_non_zoom = []
    RA_targetList_non_zoom = []
    Dec_targetList_non_zoom = []
    markerColorList_non_zoom = []
    combinedNameList_non_zoom = []
    for name, ra, dec, m in zip(combinedNameList_non, RA_targetList_non, Dec_targetList_non, markerColorList_non):
        countList_non.append(count)
        count +=1
        
        if math.sqrt(ra**2 + dec**2) <= zoom_limit:
            countList_non_zoom.append(count_zoom)
            
            RA_targetList_non_zoom.append(ra)
            Dec_targetList_non_zoom.append(dec)
            markerColorList_non_zoom.append(m)
            combinedNameList_non_zoom.append(name)
            
            count_zoom +=1
            

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
    
        r = 4.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
        ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

    
        # plot the rest
        largestEW = max(wList2)
        smallestEW = min(wList2)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for w in wList2:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
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
                
            annotate(countList[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
            xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
                
    #         if not previousNames.has_key(combinedNameList2[i]):
    #             annotate(counter,xy=(RA_targetList2[i], Dec_targetList2[i]),\
    #             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
    # 
    #             previousNames[combinedNameList2[i]] = counter
    #             counter +=1


    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non)):
            ax.plot(RA_targetList_non[i], Dec_targetList_non[i], color=markerColorList_non[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
            annotate(countList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        

    #         if not previousNames.has_key(combinedNameList_non[i]):
    #             annotate(counter,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
    #             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
    # 
    #             previousNames[combinedNameList_non[i]] = counter
    #             counter +=1


    ##########################################################################################

        xlabel(r'$\rm R.A. ~[R_{vir}]$')
        ylabel(r'$\rm Dec. ~[R_{vir}]$')
    
        ax.set_xlim(-4.0, 4.0)
        ax.set_ylim(-4.0, 4.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(3.95, 0.06),\
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
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

        maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
        plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
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
# Cylindrical Model plot
#
#
##########################################################################################
##########################################################################################

    if plot_cyl:
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
    
        r = 4.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
        ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

    
        # plot the rest
        largestEW = max(wList2)
        smallestEW = min(wList2)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for w in wList2:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # make different style markers for different colors
        for i in arange(len(markerColorList_model2)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList_model2[i] == color_maybe:
                marker = 'o'
            if markerColorList_model2[i] == color_no:
                marker = 'X'
                marker_lw = 0.5
            if markerColorList_model2[i] == color_yes:
                marker = 'D'
        

            ax.scatter(RA_targetList2[i], Dec_targetList2[i], color=markerColorList_model2[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    
    
        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)
    #         print 'combinedNameList2[i]: ',combinedNameList2[i], 'dec = ',Dec_targetList2[i]
    #         print 'yTagOffset: ',yTagOffset
                
            annotate(countList[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
            xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
                
    #         if not previousNames.has_key(combinedNameList2[i]):
    #             annotate(counter,xy=(RA_targetList2[i], Dec_targetList2[i]),\
    #             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
    # 
    #             previousNames[combinedNameList2[i]] = counter
    #             counter +=1


    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non)):
            ax.plot(RA_targetList_non[i], Dec_targetList_non[i], color=markerColorList_non[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
            annotate(countList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        

    #         if not previousNames.has_key(combinedNameList_non[i]):
    #             annotate(counter,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
    #             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
    # 
    #             previousNames[combinedNameList_non[i]] = counter
    #             counter +=1


    ##########################################################################################

        xlabel(r'$\rm R.A. ~[R_{vir}]$')
        ylabel(r'$\rm Dec. ~[R_{vir}]$')
    
        ax.set_xlim(-4.0, 4.0)
        ax.set_ylim(-4.0, 4.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(3.95, 0.06),\
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
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

        maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
        plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
                                borderpad=0.8, fontsize=legend_font, fancybox=True)


    ##########################################################################################
        
        save_name = 'SALTmap_cyl_model_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_inclim_{4}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation, inc_limit)    
    
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
        
        for name, ra, dec, c in zip(combinedNameList2, RA_targetList2, Dec_targetList2, markerColorList_model2):
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
    
        r = 4.0
        ax.plot(*xy(r,phis), c='black',ls='-',lw=0.6)
    
        ax.plot([0,0],[-3,3],c='black',ls='-',lw=0.6)
        ax.plot([-3,3],[0,0],c='black',ls='-',lw=0.6)
    
        ax.scatter(0,0,c='black',marker='*',s=25)
    ##########################################################################################

    
        # plot the rest
        largestEW = max(wList2)
        smallestEW = min(wList2)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for w in wList2:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # make different style markers for different colors
        for i in arange(len(markerColorList_NFWmodel2)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList_NFWmodel2[i] == color_maybe:
                marker = 'o'
            if markerColorList_NFWmodel2[i] == color_no:
                marker = 'X'
                marker_lw = 0.5
            if markerColorList_NFWmodel2[i] == color_yes:
                marker = 'D'
        

            ax.scatter(RA_targetList2[i], Dec_targetList2[i], color=markerColorList_NFWmodel2[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    
    
        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList2)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)
    #         print 'combinedNameList2[i]: ',combinedNameList2[i], 'dec = ',Dec_targetList2[i]
    #         print 'yTagOffset: ',yTagOffset
                
            annotate(countList[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
            xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
                
    #         if not previousNames.has_key(combinedNameList2[i]):
    #             annotate(counter,xy=(RA_targetList2[i], Dec_targetList2[i]),\
    #             xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)
    # 
    #             previousNames[combinedNameList2[i]] = counter
    #             counter +=1


    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non)):
            ax.plot(RA_targetList_non[i], Dec_targetList_non[i], color=markerColorList_non[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
            annotate(countList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        

    #         if not previousNames.has_key(combinedNameList_non[i]):
    #             annotate(counter,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
    #             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
    # 
    #             previousNames[combinedNameList_non[i]] = counter
    #             counter +=1


    ##########################################################################################

        xlabel(r'$\rm R.A. ~[R_{vir}]$')
        ylabel(r'$\rm Dec. ~[R_{vir}]$')
    
        ax.set_xlim(-4.0, 4.0)
        ax.set_ylim(-4.0, 4.0)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(3.95, 0.06),\
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
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

        maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
        plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
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
        
        for name, ra, dec, c in zip(combinedNameList2, RA_targetList2, Dec_targetList2, markerColorList_NFWmodel2):
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
# NFW model plot but zoomed into 1 R_vir radius only
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
    
        ax.scatter(0,0,c='black',marker='*',s=25)
    
    ##########################################################################################

    
        # plot the rest
        largestEW = max(wList_zoom)
        smallestEW = min(wList_zoom)
        maxSize = 500
        minSize = 30
    
        newSizeList = []
        for w in wList_zoom:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
            newSizeList.append(newSize)

        # make different style markers for different colors
        for i in arange(len(markerColorList_NFWmodel_zoom)):
            marker = '*'
            marker_lw = 0.6

            if markerColorList_NFWmodel_zoom[i] == color_maybe:
                marker = 'o'
            if markerColorList_NFWmodel_zoom[i] == color_no:
                marker = 'X'
#                 marker_lw = 1.5
                marker_lw = 0.5
            if markerColorList_NFWmodel_zoom[i] == color_yes:
                marker = 'D'
        

            ax.scatter(RA_targetList_zoom[i], Dec_targetList_zoom[i], color=markerColorList_NFWmodel_zoom[i], \
            s=newSizeList[i], marker=marker, edgecolor='black', lw=marker_lw)
    
        xTagOffset = 2.0
        yTagOffset = 1.

        previousNames = {}
        counter = 1
        for i in arange(len(combinedNameList_zoom)):
        
            yTagOffset = 5.0 + (newSizeList[i]/50.)

            annotate(countList_zoom[i],xy=(RA_targetList_zoom[i], Dec_targetList_zoom[i]),\
            xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)


    ##########################################################################################
        # now the non-detections

        non_size = 10
        non_marker = 'o'
        for i in arange(len(markerColorList_non_zoom)):
            ax.plot(RA_targetList_non_zoom[i], Dec_targetList_non_zoom[i], color=markerColorList_non_zoom[i], \
            ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

            yTagOffset = 5.0
            annotate(countList_non_zoom[i],xy=(RA_targetList_non_zoom[i], Dec_targetList_non_zoom[i]),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
        

    #         if not previousNames.has_key(combinedNameList_non[i]):
    #             annotate(counter,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
    #             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
    # 
    #             previousNames[combinedNameList_non[i]] = counter
    #             counter +=1


    ##########################################################################################

        xlabel(r'$\rm R.A. ~[R_{vir}]$')
        ylabel(r'$\rm Dec. ~[R_{vir}]$')
    
        ax.set_xlim(-zoom_limit, zoom_limit)
        ax.set_ylim(-zoom_limit, zoom_limit)
        ax.invert_xaxis()
    
        annotate(r'$\rm Approaching~ Side$',xy=(zoom_limit-0.04, 0.06),\
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
    #     yellow_line = mlines.Line2D([], [], color='blue', marker='o',lw=0,
    #                               markersize=15, label=r'$\rm \Delta v \leq 50 ~km s^{-1}$')
        corotate = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Co-rotation$')

        maybe = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label='Within Uncertainties')
                              
        antirotate = mlines.Line2D([], [], color=color_no, marker='X',lw=0,
                                  markersize=legend_size, markeredgecolor='black', label=r'$\rm Anti-rotation$')
                              
        nondetection = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0,
                                markeredgecolor='grey', markersize=legend_size, markerfacecolor = 'none', label='Non-detection')
                              
        plt.legend(handles=[corotate, maybe, antirotate, nondetection],loc='lower right', 
                                borderpad=0.8, fontsize=legend_font, fancybox=True)


    ##########################################################################################
        
        save_name = 'SALTmap_NFW_model_velstrict_{0}_non_{1}_Lstar_{2}_minsep_{3}_zoom_{4}_inclim_{5}'.format(only_close_velocities, include_nondetection, Lstar_range, min_separation,zoom_limit, inc_limit)

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
        for name, ra, dec, c in zip(combinedNameList_zoom, RA_targetList_zoom, Dec_targetList_zoom, markerColorList_NFWmodel_zoom):
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

    
if __name__ == '__main__':
    main()
    