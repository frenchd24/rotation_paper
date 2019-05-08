#!/Users/frenchd/anaconda2/bin/python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_nice_rotation_models.py, v1.0 07/12/18

Plot NFW and regular rotation models all in one in a publication quality plot

v1.1:
Re-plot with inward-facing tick marks

'''


import sys
import os
import csv
import time
from copy import deepcopy


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

from scipy.interpolate import interp1d
from scipy import interpolate
from scipy import interpolate, optimize
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline


from matplotlib import rc
fontScale = 18
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


'''
========================================================
'''

def cot(x):
    # cotangent
    return 1/math.tan(x)
    
    
def csc(x):
    # cosecant
    return 1/math.sin(x)
    

def inclination_error(v,dv,i,di):
    # calculates the quadrature error in the final velocity value
    # w = observed wavelength of Halpha line center
    # dw = error in wavelength of line center
    # v = systemic velocity
    # dv = error in systemic velocity

    i = i*math.pi/180.
    di = di*math.pi/180
    
    # wavelength term
    incTerm = (v * cot(i) * csc(i) * di)**2

    # v_sys term
    vsysTerm = (csc(i) * dv)**2

    totalError = math.sqrt(vsysTerm + incTerm)

    return totalError



def NFW(r,v200,c,r200):

    x = r/r200
    top = (np.log(1 + c*x) - c*x / (1 + c*x))

    bottom = x * (np.log(1 + c) - c / (1 + c))

    vr = v200 * np.sqrt(top/bottom)

    return vr



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




def main():
    # open the data files
    model_directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/SALT/NGC3633/'
    cyl_model_filename = '{0}/NGC3633-RX_J1121.2+0326_model_NFWFalse.p'.format(model_directory)
    NFW_model_filename = '{0}/NGC3633-RX_J1121.2+0326_model_NFWTrue.p'.format(model_directory)

    cyl_model_file = open(cyl_model_filename,'r')
    cyl_model = pickle.load(cyl_model_file)

    NFW_model_file = open(NFW_model_filename,'r')
    NFW_model = pickle.load(NFW_model_file)

    cyl_model_file.close()
    NFW_model_file.close()


    # which thing to plot?
    plot_velocity = False
    plot_3Dmodel = False
    plot_3Dmodel_movie = False
    
    plot_NFW_fit = False
    plot_rotation_curve = True
    
#     movie_directory = '/Users/frenchd/Research/inclination/git_inclination/thesis/DMF_thesis/NGC3633_movie/'
    movie_directory = '/Users/frenchd/Research/test/'

    save_directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/figures/'

    
    color_blue = '#436bad'     # french blue
    color_red = '#ec2d01'     # tomato red
    
    
##########################################################################################
    # get the data from the JSON files
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

#     galaxyName = 'ESO343-G014'
#     legend_loc = 'lower right'
#     x_ticks = 5.
#     y_ticks = 50.

#     galaxyName = 'CGCG039-137'
#     legend_loc = 'lower right'
#     x_ticks = 5.
#     y_ticks = 50.

#     galaxyName = 'IC5325'
#     legend_loc = 'lower right'
#     x_ticks = 2.
#     y_ticks = 50.

#     galaxyName = 'MCG-03-58-009'
#     legend_loc = 'lower left'
#     x_ticks = 5.
#     y_ticks = 50.

#     galaxyName = 'NGC1566'
#     legend_loc = 'lower right'
#     x_ticks = 2.
#     y_ticks = 100.

#     galaxyName = 'NGC3513'
#     legend_loc = 'lower left'
#     x_ticks = 2.
#     y_ticks = 50.
    
#     galaxyName = 'NGC3633'
#     legend_loc = 'lower right'
#     x_ticks = 2.
#     y_ticks = 50.
    
#     galaxyName = 'NGC4536'
#     legend_loc = 'lower right'
#     x_ticks = 4.
#     y_ticks = 100.
    
#     galaxyName = 'NGC4939'
#     legend_loc = 'lower right'
#     x_ticks = 10.
#     y_ticks = 100.
    
#     galaxyName = 'NGC5364'
#     legend_loc = 'lower left'
#     x_ticks = 4.
#     y_ticks = 50.
    
#     galaxyName = 'NGC5786'
#     legend_loc = 'lower right'
#     x_ticks = 2.
#     y_ticks = 100.
    
    galaxyName = 'UGC09760'
    legend_loc = 'lower left'
    x_ticks = 3.
    y_ticks = 50.


    agnName = 'RX_J1121.2+0326'

    
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
#     filename = 'NGC5364-summary4.json'

#     filename = '{0}-summary4.json'.format(galaxyName)
    filename = '{0}-summary6.json'.format(galaxyName)

    
    with open(directory+filename) as data_file:
        data = json.load(data_file)
        
        vrot_vals = data['vrot_vals']
        vrot_incCorrected_vals = data['vrot_incCorrected_vals']
        vrot_incCorrected_errs = data['vrot_incCorrected_errs']
        
        right_vrot_incCorrected_avg = data['right_vrot_incCorrected_avg']
        right_vrot_incCorrected_avg_err = data['right_vrot_incCorrected_avg_err']
        left_vrot_incCorrected_avg = data['left_vrot_incCorrected_avg']
        left_vrot_incCorrected_avg_err = data['left_vrot_incCorrected_avg_err']

        right_vrot_avg = data['right_vrot_avg']
        right_vrot_avg_err = data['right_vrot_avg_err']
        left_vrot_avg = data['left_vrot_avg']
        left_vrot_avg_err = data['left_vrot_avg_err']

        xVals = data['xVals']
        vsys_measured = int(round(data['vsys_measured'],0))
        vsys_measured_err = int(round(data['vsys_measured_err'],0))

        
        RA_galaxy = data['RAdeg']
        Dec_galaxy = data['DEdeg']
        dist = data['dist']
        majDiam = data['majDiam']
        inclination = data['inclination']
        di = data['di']
        PA = data['PA']
        agn = data['agn']
        
        v200 = data['v200']
        c = data['c']
        r200 = data['r200']
        
        R_vir = calculateVirialRadius(majDiam)


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
    if plot_velocity:

        # get the data
        cyl_model_intersect_v_list = cyl_model['intersect_v_list']
        cyl_model_v_proj_list = cyl_model['v_proj_list']

        NFW_model_intersect_v_list = NFW_model['intersect_v_list']
        NFW_model_v_proj_list = NFW_model['v_proj_list']
        

        # do the plotting
#         fig = plt.figure(figsize=(6.7,7.7))
#         fig = plt.figure(figsize=(8.0, 6.7))
        fig = plt.figure(figsize=(7.7,5.7))

        ax = fig.add_subplot(1,1,1)


        # first cylindrical
        ax.plot(cyl_model_intersect_v_list,
                    cyl_model_v_proj_list, 
                    color='black',
                    lw=2,
                    linestyle='dashed',
                    markersize=0,
                    label=r'$\rm Cylindrical$')
           
        # next NFW     
        ax.plot(NFW_model_intersect_v_list, 
                    NFW_model_v_proj_list, 
                    color='green',
                    lw=2,
                    markersize=0,
                    label = r'$\rm NFW$')
                
                
        ylabel(r'$\rm Projected~\emph{v}_{rot} ~[km~s^{-1}]$')
        xlabel(r'$\rm Velocity~Along~Sightline ~[km~s^{-1}]$')
            
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
    
        # y-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
    
    
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)
    #         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(-160, 0)
        xlim(-30, 30)

        savefig('{0}/{1}-{2}_model_plot3.pdf'.format(save_directory, galaxyName, agnName),format='pdf',bbox_inches='tight')


##########################################################################################
##########################################################################################

    if plot_NFW_fit:    
        # do the plotting
#         fig = plt.figure(figsize=(8.0, 6.7))
        fig = plt.figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(1,1,1)

        # lists for data
        xData = []
        yData = []
        yErrs = []
        
        # lists for absolute values of data (folded over rotation curve)
        yData_abs = []
        xData_abs = []
        
        
        for v, e, x in zip(vrot_incCorrected_vals, vrot_incCorrected_errs, xVals):
            xData.append(x)
            yData.append(v)
            yErrs.append(e)
            
            # fold it over for NFW plotting
            xData_abs.append(abs(x))
            yData_abs.append(abs(v))
            
        # turn them into numpy arrays
        xData_abs = np.array(xData_abs)
        yData_abs = np.array(yData_abs)
        yErrs = np.array(yErrs)

        # NFW fit for this rotation curve (NGC3633)
#         v200 = 111.91
#         c = 21.4
#         r200 = 60.03
        
        popt = v200, c, r200
    
        # how far to plot/extrapolate the NFW curve?
        x_lim = round(R_vir,0)+10
        x_fit = linspace(0, x_lim, num=1000)
    
#         scatter(xData, yData, color='black', s=40, lw=0, label = r'$\rm Data$')
#         scatter(xData_abs,
#                 yData_abs,
#                 color='black',
#                 s=40,
#                 lw=0,
#                 label = r'$\rm Observed~Rotation~Curve$')
                
        errorbar(xData_abs,
                yData_abs,
                yerr=yErrs,
                fmt='o',
                color='black',
                elinewidth=1,
                ms=6,
                label = r'$\rm Observed~Rotation~Curve$')

    #     plot(x_fit, NFW(x_fit, *popt), 'r-',label='fit: a={0}, rho={1}'.format(*popt))
    
        # plot the NFW fit
        plot(x_fit,
            NFW(x_fit, *popt),
            'r-',
            color='green',
            alpha=0.7,
            lw=2,
            label=r'$\rm NFW~Fit:~V200={0},~c={1},~R200={2}$'.format(round(v200,2),round(c,2),round(r200,2)))
        
        # plot the location of R_Vir
#         axvline(x=R_vir, 
#                 linewidth=3,
#                 alpha = 0.7,
#                 c=color_blue,
#                 label=r'$\rm 1 \emph{R}_{vir}$')
        
        legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)

        # x-axis
        majorLocator   = MultipleLocator(25)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(12.5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)

        # y axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        xlabel(r'$\rm R ~[kpc]$')
        ylabel(r'$\rm \emph{v}_{{rot}} ~[km s^{{-1}}]$')
    
    
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)
    #         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
#         ylim(-160, 0)
#         xlim(-30, 30)

        xlim(0, R_vir)
#         ylim(0, round(np.nanmax(yData),-1) + 15)
#         ylim(0, 250)
        ylim(0, 100)


        savefig('{0}/{1}-NFW_fit_Rvir.pdf'.format(save_directory, galaxyName),format='pdf',bbox_inches='tight')


##########################################################################################
##########################################################################################

    if plot_rotation_curve:    
        # do the plotting
#         fig = plt.figure(figsize=(8.0, 6.7))
        fig = plt.figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(1,1,1)

        # lists for data
        xData = []
        yData = []
        yErrs = []
        
        xRight = []
        xLeft = []
        
        yRight = []
        yLeft = []
        
        yErr_right = []
        yErr_left = []
        
        for v, e, x in zip(vrot_incCorrected_vals, vrot_incCorrected_errs, xVals):
            xData.append(x)
            yData.append(v)
            yErrs.append(e)
            
            if x > 0:
                xRight.append(x)
                yRight.append(v)
                yErr_right.append(e)
            if x <0:
                xLeft.append(x)
                yLeft.append(v)
                yErr_left.append(e)
            if x == 0:
                if mean(yRight) > 0:
                    if v > 0:
                        xRight.append(x)
                        yRight.append(y)
                        yErr_right.append(e)

                    else:
                        xLeft.append(x)
                        yLeft.append(y)
                        yErr_left.append(e)

                
                if mean(yRight) <0:
                    if v < 0:
                        xRight.append(x)
                        yRight.append(y)
                        yErr_right.append(e)

                    else:
                        xLeft.append(x)
                        yLeft.append(y)
                        yErr_left.append(e)

                        
                        
        # recalculate inclination error                
        inc_err_right = inclination_error(right_vrot_avg, right_vrot_avg_err, inclination, di)
        inc_err_left = inclination_error(left_vrot_avg, left_vrot_avg_err, inclination, di)

        
        left_vrot_err_final = np.sqrt(inc_err_left**2 + left_vrot_avg_err**2)
        right_vrot_err_final = np.sqrt(inc_err_right**2 + right_vrot_avg_err**2)

        print 'vsys_measured, vsys_measured_err,inclination, di: ',vsys_measured, vsys_measured_err,inclination, di
        print 'inc_err_right: ',inc_err_right
        print 'inc_err_left: ',inc_err_left
        print 'left_vrot_err_final: ',left_vrot_err_final
        print 'right_vrot_err_final: ',right_vrot_err_final
        print
        
        vrot_final = int(round(right_vrot_incCorrected_avg, 0))
        vrot_err_final = int(round(right_vrot_incCorrected_avg_err, 0))
            
            
        # turn them into numpy arrays
        xData = np.array(xData)
        yData = np.array(yData)
        yErrs = np.array(yErrs)
    
        # how far to plot/extrapolate the NFW curve?
        x_lim = int(round(max(max(xRight), max(xLeft)),0)) + 1
    
#         scatter(xData_abs,
#                 yData_abs,
#                 color='black',
#                 s=40,
#                 lw=0,
#                 label = r'$\rm Observed~Rot.~Curve$')
        
        errorbar(xData,
                yData,
                yerr=yErrs,
                fmt='o',
                color='black',
                elinewidth=1,
                ms=6)
        
        
        # plot the average velocities with errors        
        len_right = int(math.ceil(len(yRight)/2.0))
        len_left = int(math.floor(len(yLeft)/2.0))
    
        # define the outer 1/2 radius y values - observed
        outerRightX = xRight[:len_right]
        outerLeftX = xLeft[len_left:]
        
        outerRightY = yRight[:len_right]
        outerLeftY = yLeft[len_left:]
        err_outerRightY = yErr_right[:len_right]
        err_outerLeftY = yErr_left[len_left:]

        
        print 'outerRightX, outerRightX[-1]: ',outerRightX, outerRightX[-1]
        print 'outerLeftX, outerLeftX[0]: ',outerLeftX, outerLeftX[0]
        print 'right_vrot_avg: ',right_vrot_avg
        print 
        print 'right_vrot_incCorrected_avg = ', right_vrot_incCorrected_avg, right_vrot_incCorrected_avg_err
        print 'left_vrot_incCorrected_avg = ', left_vrot_incCorrected_avg, left_vrot_incCorrected_avg_err
        print
        
        # try recalculating the average rotation velocity
        new_vrot_right, new_vrot_right_err = weightedMean(outerRightY, err_outerRightY)
        new_vrot_left, new_vrot_left_err = weightedMean(outerLeftY, err_outerLeftY)
        
#         vrot_final = int(round(mean([abs(new_vrot_right), abs(new_vrot_left)]), 0))
        
        
#         vrot_err_final = int(round(mean([abs(new_vrot_right_err), abs(new_vrot_left_err)]), 0))
        
        print 'new_vrot_right = {0} pm {1} ; new_vrot_left = {2} pm {3}'.format(new_vrot_right, new_vrot_right_err, new_vrot_left, new_vrot_left_err)
        
        vrot_final = abs(vrot_final)
    
        vrot_plot_label = r'$\rm \overline{{\emph{{v}}_{{rot}}}}={0}\pm{1}~km s^{{-1}}$'.format(vrot_final, vrot_err_final)
    
    
        # old label:
#         r'$\rm Inc.~Corrected~\overline{{\emph{{v}}_{{rot}}}}={0}\pm{1}$'.format(int(round(right_vrot_incCorrected_avg, 0)), right_vrot_avg_err)
        
        # plot the average on the right
        right_vel_final = vrot_final
        left_vel_final = vrot_final

        if right_vrot_incCorrected_avg < 0:
            right_vel_final *= -1
        if left_vrot_incCorrected_avg < 0:
            left_vel_final *= -1
        
        # make a legend  for vhel
        plot(0, 0, lw=0, label=r'$\rm v_{{\rm hel}} = {0}\pm{1}~km s^{{-1}}$'.format(vsys_measured, vsys_measured_err))
        
        
        plot((outerRightX[-1], outerRightX[0]),
            (right_vel_final, right_vel_final),
            lw=2,
            c='green',
            solid_capstyle="butt",
            label=vrot_plot_label)
            
        x_err = np.array([outerRightX[-1], outerRightX[0]])
        y_err = np.array([right_vel_final, right_vel_final])
        fill_between(x_err, 
                    y_err - right_vrot_err_final,
                    y_err + right_vrot_err_final,
                    color='green',
                    alpha=0.2,
                    linewidth=0)
            
        
        # plot the std on the right
#         plot((outerRightX[-1], outerRightX[0]), (right_vrot_incCorrected_avg, right_vrot_incCorrected_avg),\
#         lw=2*right_vrot_incCorrected_avg_err,c='green',alpha=0.2,solid_capstyle="butt")
        
        # plot the average on the left
        plot((outerLeftX[0], outerLeftX[-1]), 
            (left_vel_final, left_vel_final),
            lw=2,
            c='green',
            solid_capstyle="butt")
        
        # plot the std on the left
#         plot((outerLeftX[0], outerLeftX[-1]), (left_vrot_incCorrected_avg, left_vrot_incCorrected_avg),\
#         lw=2*left_vrot_incCorrected_avg_err,c='green',alpha=0.2,solid_capstyle="butt")

        x_err = np.array([outerLeftX[0], outerLeftX[-1]])
        y_err = np.array([left_vel_final, left_vel_final])
        fill_between(x_err, 
                    y_err - left_vrot_err_final,
                    y_err + left_vrot_err_final,
                    color='green',
                    alpha=0.2,
                    linewidth=0)

#         legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)


        # x-axis
        majorLocator   = MultipleLocator(x_ticks)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(x_ticks/2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)

        # y axis
        majorLocator   = MultipleLocator(y_ticks)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(y_ticks/2)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        xlabel(r'$\rm R ~[kpc]$')
        ylabel(r'$\rm \emph{v}_{{rot}} ~[km s^{{-1}}]$')
    
    
#         leg = ax.legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)
    #         leg.get_frame().set_alpha(0.5)
    
        # make the legend
        leg = ax.legend(title=r'$\rm {0}$'.format(galaxyName),
                        scatterpoints=1,
                        prop={'size':12},
                        loc=legend_loc,
                        fancybox=True)

        leg.get_frame().set_alpha(1.0)
        leg._legend_box.align = "left"
            

        ax.grid(b=None,which='major',axis='both')
        xlim(-x_lim, x_lim)
#         ylim(-(round(np.nanmax(yData),-1) + 20), round(np.nanmax(yData),-1) + 20)
#         ylim(-200, 200)
        ylim(-(round(np.nanmax(yData),-1) + 50), round(np.nanmax(yData),-1) + 50)
        tight_layout()

        savefig('{0}/{1}-rotation_curve_nice4.pdf'.format(save_directory, galaxyName),format='pdf',bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################


    # plot the model figure cylinder thing
    if plot_3Dmodel:
        # first get the data
        z = NFW_model['z_sightline']
        y = NFW_model['y_sightline']
        x = NFW_model['x_sightline']
    
        # now the cylinder
        intersect = NFW_model['intersect']
        R = NFW_model['R']
        p0 = NFW_model['p0']
        p1 = NFW_model['p1']

        fig = plt.figure(figsize=(7.0,8.7))
        ax = fig.add_subplot(1,1,1,projection='3d')
        
    
        print 'x: ',x
        print 'y: ',y
        print 'z: ',z
        print
        
        plotExtent = 600
        
        # plot the sightline
        ax.plot(x, y, z, color='black', alpha = 0.6, lw=2)

        # intersect point
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
        
    ##########################################################################################
    ##########################################################################################
        # plot the cylinder

        tube,bottom,top = plot_cylinder(p0,p1,R)

        X, Y, Z = tube
        X2, Y2, Z2 = bottom
        X3, Y3, Z3 = top

        alphaTube = 0.15
        alphaBottom = 0.3
        alphaTop = 0.5
    
        colorTube = 'blue'
        colorBottom = 'red'
        colorTop = 'blue'

    #     ax=plt.subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, color=colorTube, alpha = alphaTube)
        ax.plot_surface(X2, Y2, Z2, color=colorBottom, alpha = alphaBottom)
        ax.plot_surface(X3, Y3, Z3, color=colorTop, alpha = alphaTop)
    
        ax.set_xlim(-plotExtent, plotExtent)
        ax.set_ylim(-plotExtent, plotExtent)
        ax.set_zlim(-plotExtent, plotExtent)

        # reverse the RA axis so negative is on the right
    #     ax = plt.gca()
        ax.invert_xaxis()

        # rotate the plot
    #         ax.view_init(elev=10., azim=5)
    #         ax.view_init(elev=15., azim=20)
        ax.view_init(elev=5., azim=5)

        yticks((-600, -400, -200, 0, 200, 400, 600), (600, 400, 200, 0, -200, -400, -600))
#         xticks((-600, -400, -200, 0, 200, 400, 600), (-600, -400, -200, 0, 200, 400, 600))
        xticks((-600, -300, 0, 300, 600), (-600, -300, 0, 300, 600))

    
        [t.set_va('center') for t in ax.get_yticklabels()]
        [t.set_ha('center') for t in ax.get_yticklabels()]
    
        [t.set_va('bottom') for t in ax.get_xticklabels()]
        [t.set_ha('left') for t in ax.get_xticklabels()]
    
        [t.set_va('center') for t in ax.get_zticklabels()]
        [t.set_ha('right') for t in ax.get_zticklabels()]
    
        ax.grid(True)
        ax.xaxis.pane.set_edgecolor('black')
        ax.yaxis.pane.set_edgecolor('black')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        ax.xaxis._axinfo['tick']['inward_factor'] = 0
        ax.xaxis._axinfo['tick']['outward_factor'] = 0.6
        ax.yaxis._axinfo['tick']['inward_factor'] = 0
        ax.yaxis._axinfo['tick']['outward_factor'] = 0.6
        ax.zaxis._axinfo['tick']['inward_factor'] = 0
        ax.zaxis._axinfo['tick']['outward_factor'] = 0.6
        ax.zaxis._axinfo['tick']['outward_factor'] = 0.6
    
        ax.set_xlabel(r'$\rm z ~ [kpc]$')
        ax.set_ylabel(r'$\rm R.A.~ [kpc]$')
        ax.set_zlabel(r'$\rm Dec.~ [kpc]$')
    
        x_label = r'$\rm z ~ [kpc]$'
        y_label = r'$\rm R.A.~ [kpc]$'
        z_label = r'$\rm Dec.~ [kpc]$'
    
        z_label = 'Dec. [kpc]'
    #         z_label = 'abcde'
    
        ax.set_zlabel(z_label, rotation=0, fontsize=14, labelpad=40)
    #         ax.set_xlabel(x_label, rotation=0, fontsize=14, labelpad=40)
    #         ax.set_ylabel(y_label, rotation=0, fontsize=14, labelpad=40)

    #         ax.xaxis.set_label_coords(10.0, -200.02)
    
        ax.xaxis.labelpad=18.0
        ax.yaxis.labelpad=1.0
        ax.zaxis.labelpad=28.0

        tight_layout()

        savefig('{0}/{1}-RX_J1121.2+0326_3Dmodel_plot4.pdf'.format(save_directory, galaxyName),format='pdf',bbox_inches='tight')



    # plot the model figure cylinder thing as a rotating movie
    if plot_3Dmodel_movie:
        steps = 100
    
        for i in arange(steps):
            i +=1
            
            color_blue = '#436bad'      # french blue
            color_red = '#ec2d01'     # tomato red

            color_purple = '#7570b3'
            color_purple2 = '#984ea3'
            color_purple3 = '#7570b3'
            color_purple4 = '#810f7c'

            color_green = '#1b9e77'
            color_orange = '#d95f02'
            color_pink = '#e7298a'
            color_lime = '#66a61e'
            color_yellow = '#e6ab02'
            color_brown = '#a6761d'
            color_coal = '#666666'
        
        
            # first get the data
            z = NFW_model['z_sightline']
            y = NFW_model['y_sightline']
            x = NFW_model['x_sightline']
        
            len_step = len(z)/steps
    
            # now the cylinder
            intersect = NFW_model['intersect']
            R = NFW_model['R']
            p0 = NFW_model['p0']
            p1 = NFW_model['p1']

            fig = plt.figure(figsize=(7.0,8.7))
            ax = fig.add_subplot(1,1,1,projection='3d')
        
    
            print 'x: ',x
            print 'y: ',y
            print 'z: ',z
            print
        
            plotExtent = 600
        
            # plot the sightline
#             ax.plot(x, y, z, color='black', alpha = 0.6, lw=2)
            ax.plot(x[:i*len_step], y[:i*len_step], z[:i*len_step], color='black', alpha = 0.7, lw=4)

            # intersect point
            print 'intersect: ',intersect[0],intersect[1],intersect[2]

        #     ax.plot([0,v[0]], [0,v[1]], [0,v[2]], color='green',lw=plotExtent/100)
        #     ax.plot([0,v_90[0]], [0,v_90[1]], [0,v_90[2]], color='purple',lw=plotExtent/100)
        #     ax.plot([intersect[0],v_90[0]], [intersect[1],v_90[1]], [intersect[2],v_90[2]], color='purple',lw=plotExtent/100)


            # put a star on the intersect
        #         planePoint_end = [-1.18639357e-01, 6.80095455e+02, -2.46470324e+02]
        #         planePoint_end2 = [1.18630006e-01, -6.79357841e+02, 2.48330210e+02]
        #         ax.plot([planePoint_end[0]],[planePoint_end[1]],[planePoint_end[2]],color='red',marker='*',lw=0)
        #         ax.plot([planePoint_end2[0]],[planePoint_end2[1]],[planePoint_end2[2]],color='green',marker='*',lw=0)

            
            ax.plot([intersect[0]],[intersect[1]],[intersect[2]], color=color_yellow, marker='*', lw=0)
        
        ##########################################################################################
        ##########################################################################################
            # plot the cylinder

            tube,bottom,top = plot_cylinder(p0,p1,R)

            X, Y, Z = tube
            X2, Y2, Z2 = bottom
            X3, Y3, Z3 = top

            alphaTube = 0.15
            alphaBottom = 0.8
            alphaTop = 0.9   

#             colorTube = 'blue'
#             colorBottom = 'red'
#             colorTop = 'blue'
    
            colorTube = color_purple3
            
            colorBottom = color_red
            colorTop = color_blue


        #     ax=plt.subplot(111, projection='3d')
            ax.plot_surface(X, Y, Z, color=colorTube, alpha = alphaTube)
            ax.plot_surface(X2, Y2, Z2, color=colorBottom, alpha = alphaBottom)
            ax.plot_surface(X3, Y3, Z3, color=colorTop, alpha = alphaTop)
    
            ax.set_xlim(-plotExtent, plotExtent)
            ax.set_ylim(-plotExtent, plotExtent)
            ax.set_zlim(-plotExtent, plotExtent)

            # reverse the RA axis so negative is on the right
        #     ax = plt.gca()
            ax.invert_xaxis()

            # rotate the plot
        #         ax.view_init(elev=10., azim=5)
        #         ax.view_init(elev=15., azim=20)
#             ax.view_init(elev=5., azim=5)

            
            ax.view_init(elev=10., azim=3.6*i)


            yticks((-600, -400, -200, 0, 200, 400, 600), (600, 400, 200, 0, -200, -400, -600))
    #         xticks((-600, -400, -200, 0, 200, 400, 600), (-600, -400, -200, 0, 200, 400, 600))
            xticks((-600, -300, 0, 300, 600), (-600, -300, 0, 300, 600))

    
            [t.set_va('center') for t in ax.get_yticklabels()]
            [t.set_ha('center') for t in ax.get_yticklabels()]
    
            [t.set_va('bottom') for t in ax.get_xticklabels()]
            [t.set_ha('left') for t in ax.get_xticklabels()]
    
            [t.set_va('center') for t in ax.get_zticklabels()]
            [t.set_ha('right') for t in ax.get_zticklabels()]
    
            ax.grid(True)
            ax.xaxis.pane.set_edgecolor('black')
            ax.yaxis.pane.set_edgecolor('black')
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False

            ax.xaxis._axinfo['tick']['inward_factor'] = 0
            ax.xaxis._axinfo['tick']['outward_factor'] = 0.6
            ax.yaxis._axinfo['tick']['inward_factor'] = 0
            ax.yaxis._axinfo['tick']['outward_factor'] = 0.6
            ax.zaxis._axinfo['tick']['inward_factor'] = 0
            ax.zaxis._axinfo['tick']['outward_factor'] = 0.6
            ax.zaxis._axinfo['tick']['outward_factor'] = 0.6
    
            ax.set_xlabel(r'$\rm z ~ [kpc]$')
            ax.set_ylabel(r'$\rm R.A.~ [kpc]$')
            ax.set_zlabel(r'$\rm Dec.~ [kpc]$')
    
            x_label = r'$\rm z ~ [kpc]$'
            y_label = r'$\rm R.A.~ [kpc]$'
            z_label = r'$\rm Dec.~ [kpc]$'
    
            z_label = 'Dec. [kpc]'
        #         z_label = 'abcde'
    
            ax.set_zlabel(z_label, rotation=0, fontsize=14, labelpad=40)
        #         ax.set_xlabel(x_label, rotation=0, fontsize=14, labelpad=40)
        #         ax.set_ylabel(y_label, rotation=0, fontsize=14, labelpad=40)

        #         ax.xaxis.set_label_coords(10.0, -200.02)
    
            ax.xaxis.labelpad=18.0
            ax.yaxis.labelpad=1.0
            ax.zaxis.labelpad=28.0

            tight_layout()

            savefig('{0}/{1}.jpg'.format(movie_directory, i), dpi=200, format='jpg', bbox_inches='tight')



    
if __name__ == '__main__':
    main()
    