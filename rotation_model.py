#!/usr/bin/env python

'''
a sandbox for testing the trig based rotation model

v3: Try to include inclination as well as azimuth (12/5/17)

v5: switch to trying to use a vector method with plane/vector intersept

v6: safeguard, because v5 might actually be working

'''


import sys
import os
import csv

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


from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
import json

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

    # this one is MCG-03-58-009:
    xvals =  [\
    20.03285133378786,\
    18.390814317565052,\
    16.748777306893537,\
    15.106740301277659,\
    13.464703300221768,\
    11.822666303230218,\
    10.180629309807351,\
    8.538592319457523,\
    6.896555331685079,\
    5.254518345994366,\
    3.612481361889739,\
    1.9704443788755417,\
    0.3284073964561266,\
    -2.955666568580965,\
    -4.597703552189942,\
    -6.2397405371867425,\
    -7.881777524067014,\
    -9.523814513326412,\
    -11.165851505460584,\
    -12.807888500965182,\
    -14.449925500335862]


    vels = [\
    -136.99306357918067,\
    -161.3377309818461,\
    -182.22783345028802,\
    -155.3764195309359,\
    -143.37001914695793,\
    -140.02467556492593,\
    -143.39359012880595,\
    -164.02075738062558,\
    -141.2032703400364,\
    -134.88812007947672,\
    -100.5959044733263,\
    -63.05189832956057,\
    -26.697267487548743,\
    121.42287411965845,\
    137.10680448169478,\
    142.73067631383128,\
    142.70102986566417,\
    141.01634418274625,\
    145.24345261666713,\
    167.65392909776165,\
    162.59320274009588]

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
    filename = 'CGCG039-137-summary.json'
#     filename = 'RFGC3781-summary.json'

    
    with open(directory+filename) as data_file:
        data = json.load(data_file)    
        
        vrot_vals = data['vrot_vals']
        vrot_incCorrected_vals = data['vrot_incCorrected_vals']
        xVals = data['xVals']
        inclination = data['inclination']
        vsys_measured = data['vsys_measured']
        galaxyName = data['name']

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

    vels2 = vels
    xvals2 = xvals
    

    for i in range(100):
        vels2.insert(0,lowMean)
        vels2.append(highMean)

        xvalStart +=step
        xvalEnd -=step
        xvals2.insert(0,xvalStart)
        xvals2.append(xvalEnd)


    xData = xvals2
    yData = vels2
    fitOrder = 120
    
    # reverse it?
#     xData.reverse()
#     yData.reverse()
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

    # CGCG039-137
    impact = 98.9
    az = 71.
    inc = 72.
    PA = 157.
    
    # ESO343-G014
#     impact = 466
#     az = 74.
#     inc = 89.9
#     PA = 158.
    
    # inclination is backwards, so flip it
    effectiveInc = 90.-inc


    zcutoff = 100
    lcutoff = 700
    R_vir = 300
    verbose = True
    
    v_parallel_list = []
    v_parallel_inc_list = []
    vfinal_list = []
    v_projected_list = []
    
    ds_list = []
    ds_vel_list = []
    dsfinal_list = []
    dsvfinal_list = []
    
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
    minorPhi = (90. - inc) * math.pi/180
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
    print 'N:' ,N
    print
    
    # Define ray -> [0,RA_dif,Dec_dif]
#     rayDirection = np.array([-1, 0, 0])
    rayDirection = np.array([1, 0, 0])
#     rayPoint = np.array([0, -95.169, -26.816])
#     rayPoint = np.array([0, 95.169, -26.816])

    # CGCG039-137
    rayPoint = np.array([0, 95.169, -26.816])

    # ESO343-G014
#     rayPoint = np.array([0, -464.4, 32.5])


    
##########################################################################################
##########################################################################################
    # now loop through layers of galaxy planes
    zcutoff = 200
    v_proj_list = []
    intersect_list = []
    d_plot_list = []
    intersect_point_list = []
    for i in arange(-zcutoff,zcutoff-200,.1):
#     for i in arange(-99,-97.5,.0005):
        # this is a point in the new, parallel but shifted plane
        planePoint = (p1-p) + (i * N)
        print 'planePoint: ',planePoint
        print
    
        # get intersect: find_intersect(planeNormal,planePoint,rayDirection,rayPoint)
        intersect = find_intersect(N,planePoint,rayDirection,rayPoint)
        print "intersection at", intersect
        print
        intersect_point_list.append(intersect)
        
        # this is the vector from the origin of the current plane to the intersect
        intersect_vect = intersect - (i * N)
        
        # this is the distance from the origin of the current plane to the intersect
#         p2 = intersect[0]
#         p2 = np.linalg.norm(intersect)
        p2 = np.linalg.norm(intersect_vect)
        print 'p2: ',p2
    
        try:
            v_intersect = fit(p2)
            print 'v_intersect: ',v_intersect
            print
        except Exception,e:
            v_intersect = 0
            print 'Ran out of interpolation range for {0}'.format(p2)
            print "Built in exception is {0}".format(e)
            print
            
        #######
        #######
        #######
        # angle between sightline and vector to intersect point
        
        # normal to p2 vector, or unit vector towards intersect point
#         n_p2 = intersect / np.linalg.norm(intersect)
        n_p2 = intersect_vect / np.linalg.norm(intersect_vect)
        print 'n_p2: ',n_p2
        print
    
        # cosine of angle between sightline and intersect point unit vector
        cos_alpha = n_p2.dot(rayDirection)
        alpha = math.acos(cos_alpha)
        print 'cos_alpha: :',cos_alpha
        print
        
#         v_proj = cos_alpha * v_intersect
#         v_proj = abs(cos_alpha) * v_intersect
        v_angle = math.cos(math.pi/2 - alpha)
        v_proj = abs(v_angle) * v_intersect

        print 'v_proj: ',v_proj
        print
    
        v_proj_list.append(v_proj)
#         intersect_list.append(p2)
        intersect_list.append(intersect[0])

        d = -planePoint.dot(N)
        d_plot_list.append(d)
        
    
    print 'v_proj_list: ',v_proj_list
    print 'intersect_list: ',intersect_list
    print
    
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

    plotExtent = 450
    plotXVelocity = True
    
#     from matplotlib import gridspec
    
    # tranpose the list of intersects for plotting
    ip_xlist,ip_ylist,ip_zlist = np.array(intersect_point_list).transpose()
    
    # initial figure
    fig = plt.figure(figsize=(12,8))
#     gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 


    # first plot the v_proj data
    ax = fig.add_subplot(1,2,1)

    if plotXVelocity:
        intersect_v_list = (np.array(intersect_list)/1000)*hubbleConstant
        
        ax.scatter(intersect_v_list,v_proj_list, color='black',s=10)
#         ylabel(r'$\rm v_{proj}~(km/s)$')
        ylabel(r'$\rm v_{proj} ~[km/s]$')
        xlabel(r'$\rm intersect ~[km/s]$')
        
    else:
        ax.scatter(intersect_list,v_proj_list, color='black',s=10)
        ylabel(r'$\rm v_{proj} ~(km/s)$')
        xlabel(r'$\rm intersect ~(kpc)$')
    
    tick_spacing = round((max(v_proj_list) - min(v_proj_list))/6,1)
    if tick_spacing <=0.:
            tick_spacing = 0.01

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        
    # next plot the 3d fig
    ax = fig.add_subplot(1,2,2,projection='3d')

    # the galaxy plane normal
    normal = N

    # create x,y
    xx, yy = np.meshgrid(range(plotExtent), range(plotExtent))

    # calculate corresponding z
    total = len(d_plot_list)
    count = 1
    skipNum = 1
    skipDivisor = 1
    if total >=5:
        skipNum = total/skipDivisor
        
    for d in d_plot_list:
        count +=1
        z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]
#         print 'z:',z

        # plot the surface
        if count % skipNum == 0:
            ax.plot_surface(xx, yy, z)

    ax.set_xlabel(r'$\rm z$')
    ax.set_ylabel(r'$\rm R.A.$')
    ax.set_zlabel(r'$ Dec.$')
    
    # plot the sightline
    print 'rayPoint: ',rayPoint
    z = np.ones(1000)*rayPoint[2]
    x = np.linspace(-plotExtent, plotExtent, 1000)
    y = np.ones(1000)*rayPoint[1]
    ax.plot(x, y, z, color='black',lw=plotExtent/70)
        
    # plot the intersect ray:
#     pN = i*N
#     ax.plot([pN[0],intersect[0]],[pN[1],intersect[1]],[pN[2],intersect[2]], color='green',lw=plotExtent/70)
    
#     xlim(-15,15)
#     ylim(-15,15)
    
    
    # rotate the plot
    ax.view_init(6, 23)
#     ax.view_init(1, 1)

    plt.draw()

    # plot the fit
#     ax.yaxis.tick_right()
#     ax.yaxis.set_label_position("right")
    tight_layout()
    
    directory = '/Users/frenchd/Research/test/'
    savefig('{0}/{1}_rotation_model_2.pdf'.format(directory,galaxyName),bbox_inches='tight',format='pdf')
#     plt.show()
    
    
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
        plt.show()
    

    
    
if __name__ == '__main__':
    main()
    