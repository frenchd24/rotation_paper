#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  rotation_model2.py, v2.0 1/30/18

v2: Calculate AGN position from coordinates. (1/30/18)

v3: Plot a cylinder instead of plane (2/07/18)

v4: Make a movie of the sightline piercing the halo (2/09/18)



previous versions (from rotation_model.py) :
a sandbox for testing the trig based rotation model

v3: Try to include inclination as well as azimuth (12/5/17)

v5: switch to trying to use a vector method with plane/vector intersept

v6: safeguard, because v5 might actually be working

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



'''
========================================================
'''



def adjust(impact,inc,az):
    a = math.sin(az*math.pi/180.) * impact
    
    ds = a * math.tan((90-inc)*math.pi/180.)
    
    return ds
    
    
    
def plot_cylinder(p0,p1,R):
    #   I totally jacked this code from here:
    #
    #   '''
    #     Created on Sun Oct  2 18:33:10 2016
    # 
    #     Modified from https://stackoverflow.com/questions/38076682/how-to-add-colors-to-each-individual-face-of-a-cylinder-using-matplotlib
    #     to add "end caps" and to undo fancy coloring.
    # 
    #     @author: astrokeat
    #   '''

    #axis and radius
#     p0 = np.array([1, 3, 2]) #point at one end
#     p1 = np.array([8, 5, 9]) #point at other end
#     R = 5

    #vector in direction of axis
    v = p1 - p0

    #find magnitude of vector
    mag = norm(v)

    #unit vector in direction of axis
    v = v / mag

    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])

    #make vector perpendicular to v
    n1 = np.cross(v, not_v)
    #normalize n1
    n1 /= norm(n1)

    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)

    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, 2)
    theta = np.linspace(0, 2 * np.pi, 100)
    rsample = np.linspace(0, R, 2)

    #use meshgrid to make 2d arrays
    t, theta2 = np.meshgrid(t, theta)

    rsample,theta = np.meshgrid(rsample, theta)

    #generate coordinates for surface
    # "Tube"
    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta2) * n1[i] + R * np.cos(theta2) *       n2[i] for i in [0, 1, 2]]
    # "Bottom"
    X2, Y2, Z2 = [p0[i] + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    # "Top"
    X3, Y3, Z3 = [p0[i] + v[i]*mag + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in [0, 1, 2]]

    return (X, Y, Z), (X2, Y2, Z2), (X3, Y3, Z3)
    
    
    
    
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
        
        vrot_vals = data['vrot_vals']
        vrot_incCorrected_vals = data['vrot_incCorrected_vals']
        xVals = data['xVals']
        inclination = data['inclination']
        vsys_measured = data['vsys_measured']
        


def find_intersect(planeNormal,planePoint,rayDirection,rayPoint):
    epsilon=1e-6

    ndotu = planeNormal.dot(rayDirection)

    if abs(ndotu) < epsilon:
        print "no intersection or line is within plane"
        return False
    
    else:
        w = rayPoint - planePoint
        si = -planeNormal.dot(w) / ndotu
        Psi = w + si * rayDirection + planePoint
    
        return Psi
        


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

    directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/rot_curves/'
#     filename = 'CGCG039-137-summary4.json'
#     filename = 'ESO343-G014-summary4.json'
#     filename = 'RFGC3781-summary4.json'
#     filename = 'IC5325-summary4.json'
#     filename = 'MCG-03-58-009-summary4.json'
#     filename = 'NGC1566-summary4.json'
#     filename = 'NGC3513-summary4.json'
#     filename = 'NGC3633-summary4.json'
#     filename = 'NGC4536-summary4.json'
#     filename = 'NGC4939-summary4.json'
    filename = 'NGC5364-summary4.json'



    
    with open(directory+filename) as data_file:
        data = json.load(data_file)    
        
        vrot_vals = data['vrot_vals']
        vrot_incCorrected_vals = data['vrot_incCorrected_vals']
        right_vrot_incCorrected_avg = data['right_vrot_incCorrected_avg']
        left_vrot_incCorrected_avg = data['left_vrot_incCorrected_avg']
        
        xVals = data['xVals']
        inc = data['inclination']
        vsys_measured = data['vsys_measured']
        galaxyName = data['name']
        RA_galaxy = data['RAdeg']
        Dec_galaxy = data['DEdeg']
        dist = data['dist']
        majDiam = data['majDiam']
        inc = data['inclination']
        PA = data['PA']
        agn = data['agn']

        # which agn do you want to target?

        # CGCG039-137
#         agnName = 'RX_J1121.2+0326'
#         agnName = 'SDSSJ112224.10+031802.0'
        
        # IC5325
#         agnName = 'RBS2000'
        
        # RFGC3781 or ESO343-G014
#         agnName = 'RBS1768'
        
        # MCG-03-58-009
#         agnName = 'MRC2251-178'
        
        # NGC1566
#         flipInclination = True
#         agnName = 'HE0439-5254'
#         agnName = 'HE0435-5304'
#         agnName = 'RBS567'
#         agnName = 'HE0429-5343'
#         agnName = '1H0419-577'

        # NGC3513
#         flipInclination = True
#         agnName = 'H1101-232'

        
        # NGC3633
#         flipInclination = False
#         agnName = 'SDSSJ112005.00+041323.0'
#         agnName = 'RX_J1121.2+0326'
#         agnName = 'SDSSJ112224.10+031802.0'


        # NGC4536
#         flipInclination = False
#         agnName = '3C273.0'
#         agnName = 'HE1228+0131'

        # NGC4939
#         flipInclination = True
#         reverse = True
#         agnName = 'PG1302-102'


        # NGC5364
        flipInclination = True
        reverse = True
        agnName = 'SDSSJ135726.27+043541.4'


        # grab the coordinates for this target
        RA_target = agn[agnName]['RAdeg']
        Dec_target = agn[agnName]['DEdeg']
        

##########################################################################################
##########################################################################################
    # rename arrays something stupid
    vels = vrot_vals
    xvals = xVals
    vsys = vsys_measured

    # this works for all positive xvals
    # xvalStart = xvals[0]
    # xvalEnd = xvals[-1]
    # step = 5
    #  
    # lowMean = mean(vels[:4])
    # highMean = mean(vels[-4:])
    # 
    # vels2 = vels
    # xvals2 = xvals
    # 
    # 
    # for i in range(10):
    #     vels2.insert(0,lowMean)
    #     vels2.append(highMean)
    # 
    #     xvalStart -=step
    #     xvalEnd +=step
    #     xvals2.insert(0,xvalStart)
    #     xvals2.append(xvalEnd)


    # this one is for 0 centered
    xvalStart = xvals[0]
    xvalEnd = xvals[-1]
    step = 5
 
    lowMean = mean(vels[:6])
    highMean = mean(vels[-6:])
    
    print 'lowMean: ',lowMean
    print 'highMean: ',highMean
    print
    print 'right_vrot_incCorrected_avg: ',right_vrot_incCorrected_avg
    print 'left_vrot_incCorrected_avg: ',left_vrot_incCorrected_avg
    print

    vels2 = vels
    xvals2 = xvals
    

    for i in range(100):
#         vels2.insert(0,lowMean)
#         vels2.append(highMean)
        vels2.insert(0,right_vrot_incCorrected_avg)
        vels2.append(left_vrot_incCorrected_avg)

        xvalStart +=step
        xvalEnd -=step
        xvals2.insert(0,xvalStart)
        xvals2.append(xvalEnd)
        
        
#     fig = plt.figure(figsize=(12,8))
#     ax = fig.add_subplot(1,1,1)
#     plot(xvals2,vels2)
#     show()

    xData = xvals2
    yData = vels2
    fitOrder = 120
    
    # reverse it?
    if reverse:
#     xData.reverse()
        yData.reverse()
#     yData = np.array(yData)*-1

    from scipy.interpolate import interp1d
    from scipy import interpolate
    # f = interp1d(xData, yData)
    fit = interp1d(xData, yData, kind='cubic')

    xData_fit = linspace(min(xvals),max(xvals),num=100)
    yData_fit = fit(xData_fit)

    print
    print 'here is yData_fit: ', yData_fit
    print

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

    
    # calculate impact parameter and shit
    impact = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist)
    
    # RA component of impact parameter - by setting the Dec to be the same for both
    impact_RA = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_galaxy,dist)
    
    # Dec component of impact parameter - by setting the RA to be the same for both
    impact_Dec = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_galaxy,Dec_target,dist)
    
    # calculate azimuth
    az = calculateAzimuth(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist,PA)
    
    # calculate R_vir
    R_vir = calculateVirialRadius(majDiam)
    
#     inc = 20.
    
#     inc = .45
#     R_vir = 1.5*R_vir

    
#     if RA_galaxy > RA_target:
#         impact_RA = -impact_RA
    if Dec_galaxy > Dec_target:
        impact_Dec = -impact_Dec
        
    if RA_target > RA_galaxy:
        impact_RA = -impact_RA
        
#     impact_RA -= 30.
    
    
    print
    print 'impact: ',impact
    print 'impact_RA: ',impact_RA
    print 'impact_Dec: ',impact_Dec
    print 'az: ',az
    print 'R_vir: ',R_vir
    print
        
    # inclination is backwards, so flip it
    effectiveInc = 90.-inc
#     effectiveInc += 180
    
    if flipInclination:
        effectiveInc *=-1.
    
    
    # define some lists to populate later
    v_parallel_list = []
    v_parallel_inc_list = []
    vfinal_list = []
    v_projected_list = []
    
    ds_list = []
    ds_vel_list = []
    dsfinal_list = []
    dsvfinal_list = []
    
    
    # do some shit that is currently irrevelant
    zcutoff = 100
    lcutoff = 700
#     R_vir = 300
    verbose = True
    
    maxAngle = int(math.acos(impact/lcutoff) * 180./math.pi)
    print 'maxAngle with lcutoff: ',maxAngle
    maxAngle = 80
    print 'maxAngle: ',maxAngle
    

    # l is the radial component of the sightline's impact parameter
    l = impact * math.cos(az * math.pi/180.)
    print 'l: ',l

    
    # z is the 'height' above the disk component of impact
    z = impact * math.sin(az * math.pi/180.)
    print 'z: ',z
    
    
##########################################################################################
    # Define the galaxy plane
    
    majorTheta = (180. - PA) * math.pi/180
    majorPhi = (90.) * math.pi/180
    majorR = 1.
    p1_x = majorR * math.sin(majorTheta) * math.cos(majorPhi)
    p1_y = majorR * math.sin(majorTheta) * math.sin(majorPhi)
    p1_z = majorR * math.cos(majorTheta)
    p1 = np.array([p1_x, p1_y, p1_z])
    print 'p1: ',p1
    print
    
#     minorTheta = (180. - PA + 90.) * math.pi/180
    minorTheta = 90. * math.pi/180
#     minorPhi = (90. - inc) * math.pi/180
#     minorPhi = (90. - inc) * math.pi/180
    minorPhi = (effectiveInc) * math.pi/180

    minorR = 1.
    
    p2_x = minorR * math.sin(minorTheta) * math.cos(minorPhi)
    p2_y = minorR * math.sin(minorTheta) * math.sin(minorPhi)
    p2_z = minorR * math.cos(minorTheta)
    p2 = np.array([p2_x, p2_y, p2_z])
    print 'p2: ',p2
    print

    p_z = 0.
    p_y = 0.
    p_z = 0.
    p = np.array([p_z, p_y, p_z])
    
    pp1 = p1 - p
    pp2 = p2 - p
    
    print 'pp1: ',pp1
    print
    print 'pp2: ',pp2
    print
    
    # define the normal
    N = np.cross(pp2,pp1,axisa=0, axisb=0, axisc=0)
    N = N / np.linalg.norm(N)
    print 'N:' ,N
    print
    
    # Define ray -> [0,RA_dif,Dec_dif]
#     rayDirection = np.array([-1, 0, 0])
    rayDirection = np.array([1, 0, 0])
#     rayPoint = np.array([0, -95.169, -26.816])
#     rayPoint = np.array([0, 95.169, -26.816])

    # CGCG039-137
    rayPoint = np.array([0, impact_RA, impact_Dec])
    print 'rayPoint: ',rayPoint
#     print 'should be equal to : ([0, 95.169, -26.816])'
    print
    
#     rayPoint = np.array([0, 95.169, -26.816])

    # ESO343-G014
#     rayPoint = np.array([0, -464.4, 32.5])


    
##########################################################################################
##########################################################################################
    # now loop through layers of galaxy planes
    zcutoffm = 0.5
    rcutoffm = 1.5
    zcutoff = zcutoffm * R_vir
    print 'zcutoff: ',zcutoff
    print
    rcutoff = rcutoffm * R_vir
    v_proj_list = []
    intersect_list = []
    d_plot_list = []
    intersect_point_list = []
    v_list = []
    v_90_list = []
    
    # how often to sample?
    if inc <= 60:
        s = 0.1
    elif inc <=80:
        s = 0.01
    else:
        s = 0.001
    
#     for i in arange(-zcutoff,zcutoff,.1):
#     for i in arange(-99,-97.5,.0005):
    for i in arange(-zcutoff,zcutoff,s):

        # this is a point in the new, parallel but shifted plane
        planePoint = (p1-p) + (i * N)
    
        # get intersect: find_intersect(planeNormal,planePoint,rayDirection,rayPoint)
        intersect = find_intersect(N,planePoint,rayDirection,rayPoint)
        
        # this is the vector from the origin of the current plane to the intersect
        intersect_vect = intersect - (i * N)
        
        # this is the distance from the origin of the current plane to the intersect
#         p2 = intersect[0]
#         p2 = np.linalg.norm(intersect)
        p2 = np.linalg.norm(intersect_vect)

        # restrict the intersection to be within the cylinder of radius rcutoff
        if p2 <= rcutoff:
            print 'planePoint: ',planePoint
            print "intersection at", intersect
            print 'p2: ',p2

            # find the rotation velocity at this distance from the rotation curve fit center
            try:
                v_intersect = fit(p2)
                print 'v_intersect: ',v_intersect
            except Exception,e:
                # if you go beyond the fit, set velocity to 0
                v_intersect = 0
                print 'Ran out of interpolation range for {0}'.format(p2)
                print "Built in exception is {0}".format(e)
            
            #######
            #######
            #######
            # angle between sightline and vector to intersect point
        
            # unit vector towards intersect point
    #         n_p2 = intersect / np.linalg.norm(intersect)
            n_p2 = intersect_vect / np.linalg.norm(intersect_vect)
            print 'n_p2: ',n_p2
            print 'np.linalg.norm(n_p2): ',np.linalg.norm(n_p2)
        
        
            # new way of doing this:
            #
            # the result of rotating v counterclockwise by a about n is given by:
            # (cos a)v+(sin a)(n x v)
            #
            # so need to rotate by pi + pi/2 to get all the way around
            alpha = math.pi + math.pi/2
    #         alpha = math.pi/2

        
            # this is the velocity vector in the direction of intersect point, n_p2
            # edit: seems legit
            v = v_intersect * n_p2
            print 'new way: '
            print 'v: ',v
            print '||v|| : ',np.linalg.norm(v)
            v_list.append(v)
        
            # this then should be the correct rotation velocity vector, but centered at the
            # origin. So, we then need to just shift the sightline to pass through the origin
            #
            # i.e., new sightline = [1, 0 ,0 ] = rayDirection
            v_90 = math.cos(alpha) * v + math.sin(alpha) * (np.cross(N,v,axisa=0, axisb=0, axisc=0))
            print 'v_90: ',v_90
            print '||v_90|| : ',np.linalg.norm(v_90)
            print '||N||: ', np.linalg.norm(N)

            v_90_list.append(v_90)
        
            # now dot it with the sightline to get the component along
            cos_alpha = np.dot(v_90,rayDirection)
            print 'cos_alpha: ',cos_alpha
            v_proj = cos_alpha
            print 'v_proj: ',v_proj
        
    #         print 'old way:'
    
            # old way of doing this
            # cosine of angle between sightline and intersect point unit vector
    #         cos_alpha = n_p2.dot(rayDirection)
    #         alpha = math.acos(cos_alpha)
    #         print 'cos_alpha: :',cos_alpha
    #         print
    #         v_angle = math.cos(math.pi/2 - alpha)
    #         v_proj = abs(v_angle) * v_intersect
    # 
    #         print 'v_proj: ',v_proj
    #         print
    
            v_proj_list.append(v_proj)
    #         intersect_list.append(p2)
            intersect_list.append(intersect[0])
            print 'intersect[0]: ',intersect[0]
            print
            intersect_point_list.append(intersect)


            d = -planePoint.dot(N)
            d_plot_list.append(d)
        
    
#     print 'v_proj_list: ',v_proj_list
#     print 'intersect_list: ',intersect_list
    print
    
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

    # how big to make the plot boundary?
    plotExtent = 3*R_vir
    
    # how big to make the plotted cylinder?
    zHeight = zcutoff
    
    # plot velocity on the x-axis? No idea what happens if this is False...
    plotXVelocity = True
    
    # how many steps to take while plotting. each step moves the sightline forward and 
    # populates the other with that result
    steps = 50
        
#     from matplotlib import gridspec
    
    # tranpose the list of intersects for plotting
    ip_xlist,ip_ylist,ip_zlist = np.array(intersect_point_list).transpose()
    

    for i in arange(steps):
        i +=1
        print 'i: ',i
        # initial figure
        fig = plt.figure(figsize=(12,8))
        
    #     gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

        # first plot the v_proj data
        ax = fig.add_subplot(1,2,1)
        fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)


        # number of actual intercept elements to take for each step
        len_step = len(intersect_list)/steps
        
        intersect_v_list = (np.array(intersect_list)/1000)*hubbleConstant
        
        print 'len_step: ',len_step
        print 'intersect_v_list[:i*len_step] vs v_proj_list[:i*len_step]: ',intersect_v_list[:i*len_step],' , ',v_proj_list[:i*len_step]
        
        if plotXVelocity:
        
            ax.scatter(intersect_v_list[:i*len_step],v_proj_list[:i*len_step], color='black',s=10)
            ylabel(r'$\rm v_{proj} ~[km/s]$')
            xlabel(r'$\rm intersect ~[km/s]$')
        
        else:
            ax.scatter(intersect_list,v_proj_list, color='black',s=10)
            ylabel(r'$\rm v_{proj} ~(km/s)$')
            xlabel(r'$\rm intersect ~(kpc)$')
    
#         ax.scatter(intersect_v_list,v_proj_list, color='black',s=10)
        
#         v_proj_dif = max(v_proj_list) - min(v_proj_list)
        tick_num = 10
        tick_spacing = round((max(v_proj_list) - min(v_proj_list))/tick_num,2)
        print 'tick_spacing: ',tick_spacing
        print 'max(v_proj_list): ',max(v_proj_list)
        print 'min(v_proj_list): ',min(v_proj_list)
#         if tick_spacing <0.05:
#                 tick_spacing = 0.01


        xlim_pos = max(intersect_v_list) +1
        xlim_neg = min(intersect_v_list) -1
#         ylim_pos = max(v_proj_list) +2
#         ylim_neg = min(v_proj_list) -2
        ylim_pos = max(v_proj_list) + tick_spacing*2
        ylim_neg = min(v_proj_list) - tick_spacing*2
        
        print '########################################'
        print 'xlim_pos: ',xlim_pos
        print 'xlim_neg: ',xlim_neg
        print 'ylim_pos: ',ylim_pos
        print 'ylim_neg: ',ylim_neg
        print '########################################'
        print
        print 'intersect_v_list: ',intersect_v_list
        print
        print 'v_proj_list: ',v_proj_list
        print
        
        ax.set_xlim(xlim_neg, xlim_pos)
        ax.set_ylim(ylim_neg,ylim_pos)

        ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        
        # next plot the 3d fig
        ax = fig.add_subplot(1,2,2,projection='3d')

        # the galaxy plane normal
        normal = N
    
##########################################################################################
##########################################################################################

        # plot the sightline
        print 'rayPoint: ',rayPoint
        z = np.ones(1000)*rayPoint[2]
        x = np.linspace(-plotExtent, plotExtent, 1000)
        y = np.ones(1000)*rayPoint[1]
                
        len_step = len(z)/steps
        
        ax.plot(x[:i*len_step], y[:i*len_step], z[:i*len_step], color='black',lw=plotExtent/200)

        # some interesting points: 
        # v is the velocity vector in the direction of intersect point
        v = np.array(v_list[0])
    
        # v_90 is the rotation velocity vector in the direction of rotation at the intersect
        # point
        v_90 = np.array(v_90_list[0])
    
        # galaxy center
        orig = np.array([0,0,0])
    
        # intersect point
        intersect = np.array(intersect_point_list[0])
        print 'intersect: ',intersect[0],intersect[1],intersect[2]

    #     ax.plot([0,v[0]], [0,v[1]], [0,v[2]], color='green',lw=plotExtent/100)
    #     ax.plot([0,v_90[0]], [0,v_90[1]], [0,v_90[2]], color='purple',lw=plotExtent/100)
    #     ax.plot([intersect[0],v_90[0]], [intersect[1],v_90[1]], [intersect[2],v_90[2]], color='purple',lw=plotExtent/100)


        # put a star on the intersect
#         planePoint_end = [-1.18639357e-01, 6.80095455e+02, -2.46470324e+02]
#         planePoint_end2 = [1.18630006e-01, -6.79357841e+02, 2.48330210e+02]
#         ax.plot([planePoint_end[0]],[planePoint_end[1]],[planePoint_end[2]],color='red',marker='*',lw=0)
#         ax.plot([planePoint_end2[0]],[planePoint_end2[1]],[planePoint_end2[2]],color='green',marker='*',lw=0)

        ax.plot([intersect[0]],[intersect[1]],[intersect[2]],color='red',marker='*',lw=0)

        ax.set_xlim(-plotExtent, plotExtent)
        ax.set_ylim(-plotExtent, plotExtent)
        ax.set_zlim(-plotExtent, plotExtent)
        
##########################################################################################
##########################################################################################
        # plot the cylinder
    
        R = int(rcutoff)
        p0 = normal * (zHeight)
        p1 = normal * (-zHeight)
    
        tube,bottom,top = plot_cylinder(p0,p1,R)
    
        X, Y, Z = tube
        X2, Y2, Z2 = bottom
        X3, Y3, Z3 = top
    
        alphaTube = 0.2
        alphaBottom = 0.5
        alphaTop = 0.5
    
    #     ax=plt.subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, color='blue',alpha = alphaTube)
        ax.plot_surface(X2, Y2, Z2, color='blue',alpha = alphaBottom)
        ax.plot_surface(X3, Y3, Z3, color='blue',alpha = alphaTop)


        ax.set_xlabel(r'$\rm z$')
        ax.set_ylabel(r'$\rm R.A.$')
        ax.set_zlabel(r'$ Dec.$')
    
        # reverse the RA axis so negative is on the right
    #     ax = plt.gca()
        ax.invert_xaxis()

        # rotate the plot
#         ax.view_init(elev=10., azim=20)
        ax.view_init(elev=10., azim=15)

    
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
        if i == 49:
            show()
    
#         directory = '/Users/frenchd/Research/test/CGCG039-137_3/'
#         directory = '/Users/frenchd/Research/test/ESO343-G014/'
#         directory = '/Users/frenchd/Research/test/IC5325/'
#         directory = '/Users/frenchd/Research/test/MCG-03-58-009/'
#         directory = '/Users/frenchd/Research/test/NGC1566/'
#         directory = '/Users/frenchd/Research/test/NGC3513/'
#         directory = '/Users/frenchd/Research/test/NGC3633/'
#         directory = '/Users/frenchd/Research/test/NGC4536/'
#         directory = '/Users/frenchd/Research/test/NGC4939/'
        directory = '/Users/frenchd/Research/test/NGC5364/'


#         directory = '/Users/frenchd/Research/test/RFGC3781/'
        savefig("{0}{1}.jpg".format(directory,i),dpi=200,bbox_inches='tight')
        close(fig)

#         for ii in xrange(0,360,5):
#             ax.view_init(elev=10., azim=ii)
#             plt.draw()
#     
#             savefig("{0}movie{1}.jpg".format(directory,ii),dpi=120)
        
#     savefig('{0}/{1}_rotation_model_cylinder.pdf'.format(directory,galaxyName),bbox_inches='tight',format='pdf')
    
    
    
#     plt.show()

    
    # create x,y
#     xx, yy = np.meshgrid(range(-plotExtent,plotExtent), range(-plotExtent,plotExtent))
# 
# 
    # calculate corresponding z
#     total = len(d_plot_list)
#     count = 1
#     skipNum = 1
#     skipDivisor = 1
#     if total >=5:
#         skipNum = total/skipDivisor
# 
# 
#     for d in d_plot_list:
#         count +=1
#         z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]
#         print '-normal[0]: ',normal[0]
#         print 'normal[1]: ',normal[1]
#         print 'd: ',d
#         print 'normal[2]: ',normal[2]
#         print
#         print 'z:',z
#         print 'xx: ',xx
#         print
#         
        # plot the surface
#         if count % skipNum == 0:
#             ax.plot_surface(xx, yy, z)
    
    
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
    plotCircle = False
    if plotCircle:
        from mpl_toolkits.mplot3d import axes3d
        from matplotlib.patches import Circle, PathPatch
#         import matplotlib.pyplot as plt
        from matplotlib.transforms import Affine2D
        from mpl_toolkits.mplot3d import art3d
    #     import numpy as np

        def plot_vector(fig, orig, v, color='blue'):
           ax = fig.gca(projection='3d')
           orig = np.array(orig); v=np.array(v)
           ax.quiver(orig[0], orig[1], orig[2], v[0], v[1], v[2],color=color)
           ax.set_xlim(0,10);ax.set_ylim(0,10);ax.set_zlim(0,10)
           ax = fig.gca(projection='3d')  
           return fig

        def rotation_matrix(d):
            sin_angle = np.linalg.norm(d)
            if sin_angle == 0:return np.identity(3)
            d /= sin_angle
            eye = np.eye(3)
            ddt = np.outer(d, d)
            skew = np.array([[    0,  d[2],  -d[1]],
                          [-d[2],     0,  d[0]],
                          [d[1], -d[0],    0]], dtype=np.float64)

            M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
            return M

        def pathpatch_2d_to_3d(pathpatch, z, normal):
            if type(normal) is str: #Translate strings to normal vectors
                index = "xyz".index(normal)
                normal = np.roll((1.0,0,0), index)

            normal /= np.linalg.norm(normal) #Make sure the vector is normalised
            path = pathpatch.get_path() #Get the path and the associated transform
            trans = pathpatch.get_patch_transform()

            path = trans.transform_path(path) #Apply the transform

            pathpatch.__class__ = art3d.PathPatch3D #Change the class
            pathpatch._code3d = path.codes #Copy the codes
            pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color    

            verts = path.vertices #Get the vertices in 2D

            d = np.cross(normal, (0, 0, 1)) #Obtain the rotation vector    
            M = rotation_matrix(d) #Get the rotation matrix

            pathpatch._segment3d = np.array([np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts])

        def pathpatch_translate(pathpatch, delta):
            pathpatch._segment3d += delta

        def plot_plane(ax, point, normal, size=10, color='y',alpha='blue'):    
            p = Circle((0, 0), size, facecolor = color, alpha = alpha)
            ax.add_patch(p)
            pathpatch_2d_to_3d(p, z=0, normal=normal)
            pathpatch_translate(p, (point[0], point[1], point[2]))


    #     o = np.array([5,5,5])
    #     v = np.array([3,3,3])
    #     n = [0.5, 0.5, 0.5]
        o = pp1
        v = pp2
        n = N

        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')  
        plot_plane(ax, o, n, size=10,color='blue',alpha=0.7)    
        ax.set_xlim(-10,10);ax.set_ylim(-10,10);ax.set_zlim(-10,10)
        xlabel('x')
        ylabel('y')
        ax.invert_xaxis()
        
        # rotate the plot
        ax.view_init(6, 23)
    #     ax.view_init(1, 1)

        plt.draw()

        plt.show()
    

    
    
if __name__ == '__main__':
    main()
    