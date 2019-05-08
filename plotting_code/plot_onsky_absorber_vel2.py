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

    print 'current file: ',filename
    print
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
    filename = '{0}salt_galaxy_sightlines_cut.csv'.format(directory)
    theFile = open(filename,'rU')
    tableReader = csv.DictReader(theFile)

    nameList = []
    targetList = []
    combinedNameList = []
    vList = []
    wList = []
#     RA_galaxyList = []
#     Dec_galaxyList = []
    RA_targetList = []
    Dec_targetList = []
    incList = []
    paList = []
    azList = []
    RvirList = []
    markerColorList = []
    VhelList = []
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

        
#         theta = math.atan(impact_Dec/impact_RA) - (x * math.pi/180.)
#         
#         impact_RA2 = impact * math.cos(theta)
#         impact_Dec2 = impact * math.sin(theta)
#         
#         print
#         print 'x: ',x
#         print 'theta: ',theta
#         print 'impact_RA2: ',impact_RA2
#         print 'impact_Dec2: ',impact_Dec2
        
#         print 'Name: ',name
#         print 'right_vrot_avg: ',right_vrot_avg
#         print 'left_vrot_avg: ',left_vrot_avg
#         print 'impact_RA: ',impact_RA
#         print 'impact_Dec: ',impact_Dec
#         print 'RA_galaxy vs RA_target: ',RA_galaxy,', ',RA_target
#         print 'Dec_galaxy vs Dec_target: ',Dec_galaxy,', ',Dec_target
#         print 'PA: ',PA
#         print 'Vhel: ',vsys_measured
#         print 'az, az2: ',az, ', ',az2
        
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
        dv = lyaV - vsys_measured
        print 'lyaV, dv = ',lyaV,dv
                
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
        color_maybe = '#7570b3'  # blue-ish purple
        
        markerColor = 'black'
        
        # the absorption and rotation velocity match
        if (dv > 0 and rot_vel > 0) or (dv < 0 and rot_vel < 0):
            markerColor = color_yes
            
        # mismatch
        elif (dv < 0 and rot_vel > 0) or (dv > 0 and rot_vel < 0):
            markerColor = color_no
            
        else:
            markerColor = 'black'
            
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

#         if bfind(combinedName,'_') and not bfind(combinedName,'\_'):
#             combinedName = combinedName.replace('_','\_')
            
            
        # populate the lists
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

#         print 'added'
#         print 'impact_RA_vir: ',impact_RA_vir
#         print 'impact_Dec_vir: ',impact_Dec_vir
#         print
#         print
        
##########################################################################################
    # impact = impact parameter
    # az = azimuth
    # inc = inclination (adjustedInc)

#     CGCG039-137
#     target = 'RX_J1121.2+0326'
#     RA_target = agn[target]['RAdeg']
#     Dec_target = agn[target]['DEdeg']
#     
#     RA_galaxy = 170.36231
#     Dec_galaxy = 3.44491
#     dist = 101.21
#     majDiam = 26.35
#     
#     impact = 98.9
#     R_vir = 166.09
#     az = 71.
#     inc = 63.
#     PA = 157.
#     
#     ESO343-G014
#     target = 'RBS1768'
#     RA_target = 324.7079167
#     Dec_target = -38.47777778
#     
#     RA_galaxy = 324.43825
#     Dec_galaxy = -38.49256
#     dist = 126.07
#     majDiam = 45.23
#     
#     impact = 465.6
#     az = 74.
#     inc = 89.9
#     PA = 158.
#     
#     IC5325
#     target = 'RBS2000'
#     RA_target = 351.18625
#     Dec_target = -40.68027778
#     
#     RA_galaxy = 352.18096
#     Dec_galaxy = -41.33347
#     dist = 18.1
#     majDiam = 20.45
#     
#     impact = 314.335827
#     az = 64.1
#     inc = 25.
#     PA = 15.
#     
#     MCG-03-58-009
#     target = 'MRC2251-178'
#     RA_target = 343.5245833
#     Dec_target = -17.58194444
#     
#     RA_galaxy = 343.42021
#     Dec_galaxy = -17.47889
#     dist = 142.
#     majDiam = 75.31
#     
#     impact = 355.0699641
#     az = 71.
#     inc = 49.
#     PA = 27.
#     
#     NGC1566
#     target = '1H0419-577'
#     RA_target = 66.50291667
#     Dec_target = -57.20055556
#     
#     RA_galaxy = 65.00175
#     Dec_galaxy = -54.93781
#     dist = 7.19
#     majDiam = 15.34
#     
#     impact = 302.77587
#     az = 9.8
#     inc = 48.
#     PA = 170.
# 
#     NGC1566
#     target = 'HE0429-5343'
#     RA_target = 67.66666667
#     Dec_target = -53.61555556
#     
#     RA_galaxy = 65.00175
#     Dec_galaxy = -54.93781
#     dist = 7.19
#     majDiam = 15.34
#     
#     impact = 256.2063291
#     az = 60.1
#     inc = 48.
#     PA = 170.
# 
#     NGC1566
#     target = 'RBS567'
#     RA_target = 69.91125
#     Dec_target = -53.19194444
#     
#     RA_galaxy = 65.00175
#     Dec_galaxy = -54.93781
#     dist = 7.19
#     majDiam = 15.34
#     
#     impact = 422.6192722
#     az = 69.3
#     inc = 48.
#     PA = 170.
    

##########################################################################################
    # Define the galaxy plane



    
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
        print 'name - markerColor',name,' - ',c
    
    
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

    # initial figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1,1,1)
#     fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)


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

    
#     print 'full RA_targetList2: ',RA_targetList2
#     print
#     print
#     print 'full Dec_targetList2: ',Dec_targetList2
#     print
    
    ax.scatter(0,0,c='black',marker='*',s=25)
    
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
        s=newSizeList[i], marker=marker, edgecolor='black', lw=0.8)
    
    
    xTagOffset = -35.0
    yTagOffset = 0.1
    
    # put some labels on it
#     for i in arange(len(combinedNameList2)):
#         annotate(combinedNameList2[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
#         xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=8)

    previousNames = {}
    for i in arange(len(combinedNameList2)):
        print 'newSizeList[i]/50. : ', newSizeList[i]/65.
        print '(newSizeList[i]/50.)**-1 : ',(newSizeList[i]/65.)**-1
        print
        yTagOffset = 4.0 + (newSizeList[i]/65.)**-1
        
        
#         xTagOffset = newSizeList[i]/50.
        
        if not previousNames.has_key(combinedNameList2[i]):
            annotate(combinedNameList2[i],xy=(RA_targetList2[i], Dec_targetList2[i]),\
            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=7)

            previousNames[combinedNameList2[i]] = 1.
        else:
            previousNames[combinedNameList2[i]] += 1.

#         print 'previousName: ',previousName
#         print 'currentName: ',combinedNameList2[i]
#         print
#         previousName = combinedNameList2[i]

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
                              
    plt.legend(handles=[yellow_line,red_line,blue_line],loc='lower right')


        
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
    
#     show()
    
    directory = '/Users/frenchd/Research/test/'
#     savefig("{0}SALT_map1.pdf".format(directory),dpi=400,bbox_inches='tight')
    savefig("{0}SALT_map6.pdf".format(directory),bbox_inches='tight')


    
if __name__ == '__main__':
    main()
    