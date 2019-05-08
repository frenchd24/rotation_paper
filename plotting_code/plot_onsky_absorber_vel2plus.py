#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_onsky_absorber_vel.py, v2.0 04/04/18

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
#         right_vrot_incCorrected_avg = data['right_vrot_incCorrected_avg']
#         left_vrot_incCorrected_avg = data['left_vrot_incCorrected_avg']

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
#         agn = data['agn']
        
        return vsys_measured, right_vrot_avg, left_vrot_avg, inc, PA, dist, majDiam
        


def main():
    hubbleConstant = 71.0


##########################################################################################
    # get the data
##########################################################################################
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

#     directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/rot_curves/'
#     filename = 'CGCG039-137-summary4.json'
#     filename = 'RFGC3781-summary4.json'

    
#     with open(directory+filename) as data_file:
#         data = json.load(data_file)    
#         
#         vrot_vals = data['vrot_vals']
#         vrot_incCorrected_vals = data['vrot_incCorrected_vals']
#         right_vrot_incCorrected_avg = data['right_vrot_incCorrected_avg']
#         left_vrot_incCorrected_avg = data['left_vrot_incCorrected_avg']
#         
#         xVals = data['xVals']
#         inc = data['inclination']
#         vsys_measured = data['vsys_measured']
#         galaxyName = data['name']
#         RA_galaxy = data['RAdeg']
#         Dec_galaxy = data['DEdeg']
#         dist = data['dist']
#         majDiam = data['majDiam']
#         inc = data['inclination']
#         PA = data['PA']
#         agn = data['agn']

        # which agn do you want to target?

        # CGCG039-137
#         agnName = 'RX_J1121.2+0326'
#         RA_target = agn[agnName]['RAdeg']
#         Dec_target = agn[agnName]['DEdeg']
        
        # IC5325
#         agnName = 'RBS2000'
#         RA_target = agn[agnName]['RAdeg']
#         Dec_target = agn[agnName]['DEdeg']
        
        # RFGC3781
#         agnName = 'RBS1768'
#         RA_target = agn[agnName]['RAdeg']
#         Dec_target = agn[agnName]['DEdeg']
        
    # this one is for 0 centered
#     xvalStart = xvals[0]
#     xvalEnd = xvals[-1]
#     step = 5
#  
#     lowMean = mean(vels[:6])
#     highMean = mean(vels[-6:])
# 
#     vels2 = vels
#     xvals2 = xvals

##########################################################################################
##########################################################################################
    # collect data
    directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/'
#     filename = '{0}salt_galaxy_sightlines_cut.csv'.format(directory)
    filename = '{0}salt_galaxy_sightlines_cut_plus_ancillary.csv'.format(directory)    
    
    theFile = open(filename,'rU')
    tableReader = csv.DictReader(theFile)

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
        
        PA_observed = eval(t['PA_observed'])
        PA_adjust = eval(t['PA_adjust'])
        
        gfilename = directory + 'rot_curves/' + name + '-summary4.json'
        vsys_measured, right_vrot_avg, left_vrot_avg, inc, PA, dist, majDiam = get_data(gfilename)
        
        
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
                
        # check which 'side' of the galaxy the absorber is found
        if name == 'CGCG039-137':
             # regular - left side of rotation curve is 'left' on sky
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg

        if name == 'ESO343-G014':
             # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'IC5325':
             # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'MCG-03-58-009':
             # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg

        if name == 'NGC1566':
             # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC3513':
             # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC3633':
             # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC4536':
            # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC4939':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC5364':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC5786':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
            print 'impact_RA_vir: ',impact_RA_vir
            print 'rot_vel: ',rot_vel
                
        if name == 'UGC09760':
            # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
##########################################################################################
##########################################################################################
          
        if name == 'NGC3198':
            # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC4565':
            # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'UGC04238':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC3351':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC4529':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC6140':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg

        if name == 'NGC5907':
            # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg

        if name == 'UGC06446':
            # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg

        if name == 'UGC06399':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg
                
        if name == 'NGC3726':
            # regular
            if impact_RA_vir > 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg

        if name == 'NGC3067':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg

        if name == 'NGC2770':
            # reverse
            if impact_RA_vir < 0:
                rot_vel = left_vrot_avg
            else:
                rot_vel = right_vrot_avg


#         if az2 > az and impact_Dec >0:
#             # if az increases when increasing PA, then the target is on the "left" side of
#             # the galaxy
#             rot_vel = left_vrot_avg
#             print 'az2 < az, using left: ',rot_vel
#         elif az2 < az and impact_Dec < 0:
#             rot_vel = left_vrot_avg
#             print 'az2 < az, using left: ',rot_vel
# 
#         else:
#             rot_vel = right_vrot_avg
#             print 'az2 > az, using right: ',rot_vel

            
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
        else:
            print 'dv == x :', dv, name
            print
            markerColor = color_nonDetection
            
#         if abs(dv - rot_vel) <= 50 or az >= 85.:
#             markerColor = color_maybe
   
        if az >= 85.:
            markerColor = color_maybe
            
        if name == 'NGC3513':
            markerColor = color_maybe

#         name = name.replace('-','\-')
#         name = name.replace('_','\_')
#         target = name.replace('-','\-')
#         target = name.replace('_','\_')
#         target = name.replace('+','\+')

        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = False

#         combinedName = '${\rm '+name + '-' + target+'}$'.format(name,target)
        combinedName = r'$\rm {0} : {1}$'.format(name,target)
            
        # populate the lists
        
        # first detections
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
        else:
            # then non-detections
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
    for ra,dec,c,w,name in zip(RA_targetList,Dec_targetList,markerColorList, wList,combinedNameList):
        orderedList.append([w,[ra,dec,c,name]])
        
    orderedList.sort(reverse=True)
    RA_targetList2 = []
    Dec_targetList2 = []
    markerColorList2 = []
    wList2 = []
    combinedNameList2 = []
    for i in orderedList:
        w, rest = i
        ra, dec, c, name = rest
        RA_targetList2.append(ra)
        Dec_targetList2.append(dec)
        markerColorList2.append(c)
        wList2.append(w)
        combinedNameList2.append(name)
#         print 'name - markerColor',name,' - ',c
    
    
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

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
        if markerColorList2[i] == color_maybe:
            marker = 'o'
        if markerColorList2[i] == color_no:
            marker = 'x'
        if markerColorList2[i] == color_yes:
            marker = 'D'
            
        ax.scatter(RA_targetList2[i], Dec_targetList2[i], color=markerColorList2[i], \
        s=newSizeList[i], marker=marker, edgecolor='black', lw=0.6)
    
    
#     xTagOffset = -35.0
#     yTagOffset = 0.1
    xTagOffset = 2.0
    yTagOffset = 1.

    previousNames = {}
    counter = 1
    for i in arange(len(combinedNameList2)):
#         print 'newSizeList[i]/50. : ', newSizeList[i]/65.
#         print '(newSizeList[i]/50.)**-1 : ',(newSizeList[i]/65.)**-1
#         yTagOffset = 5.0 + (newSizeList[i]/65.)**-1
        
        yTagOffset = 5.0 + (newSizeList[i]/50.)
        print 'combinedNameList2[i]: ',combinedNameList2[i], 'dec = ',Dec_targetList2[i]
        print 'yTagOffset: ',yTagOffset
                
        if not previousNames.has_key(combinedNameList2[i]):
#             annotate(combinedNameList2[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
            annotate(counter,xy=(RA_targetList2[i], Dec_targetList2[i]),\
            xytext=(xTagOffset, yTagOffset),textcoords='offset points',size=7)

            previousNames[combinedNameList2[i]] = counter
            counter +=1
        else:
#             previousNames[combinedNameList2[i]] += 1.
            pass

##########################################################################################
    # now the non-detections

    non_size = 10
    non_marker = 'o'
    for i in arange(len(markerColorList_non)):
        ax.plot(RA_targetList_non[i], Dec_targetList_non[i], color=markerColorList_non[i], \
        ms=non_size, marker=non_marker, markeredgecolor='grey', lw=0.8, markerfacecolor='none')

        yTagOffset = 5.0

        if not previousNames.has_key(combinedNameList_non[i]):
#             annotate(combinedNameList_non[i],xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
#             xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)
            annotate(counter,xy=(RA_targetList_non[i], Dec_targetList_non[i]),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)

            previousNames[combinedNameList_non[i]] = counter
            counter +=1
        else:
#             previousNames[combinedNameList_non[i]] += 1.
            pass

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
    yellow_line = mlines.Line2D([], [], color=color_maybe, marker='o',lw=0,
                              markersize=15, label='Within Uncertainties')
                              
    red_line = mlines.Line2D([], [], color=color_no, marker='x',lw=0,
                              markersize=15, label=r'$\rm Anti-rotation$')
                              
    blue_line = mlines.Line2D([], [], color=color_yes, marker='D',lw=0,
                              markersize=15, label=r'$\rm Co-rotation$')
                              
    non_legend = mlines.Line2D([], [], color=color_nonDetection, marker='o',lw=0, markeredgecolor='grey',
                              markersize=15, markerfacecolor = 'none', label='Non-detection')
                              
    plt.legend(handles=[yellow_line,red_line,blue_line, non_legend],loc='lower right')


        
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
    
#     show()
    
    directory = '/Users/frenchd/Research/test/'
#     savefig("{0}SALT_map1.pdf".format(directory),dpi=400,bbox_inches='tight')
    savefig("{0}SALT_map_plus2.pdf".format(directory),bbox_inches='tight')

    summary_filename = '{0}SALT_map_plus2_summary.txt'.format(directory)
    summary_file = open(summary_filename,'wt')

#     keys = previousNames.keys()
#     values = previousNames.values()
    
#     reverse_name_dict = {}
#     for k,v in zip(keys,values):
#         reverse_name_dict[v] = k

#     for k,v in zip(keys,values):
#         summary_file.write('{0} : {1}\n'.format(v,k))
    
    for k, v in sorted(previousNames.iteritems(), key=lambda (k,v): (v,k)):
        summary_file.write('{0}. {1}, \n'.format(v,k))
    
    summary_file.close()
    
    
if __name__ == '__main__':
    main()
    