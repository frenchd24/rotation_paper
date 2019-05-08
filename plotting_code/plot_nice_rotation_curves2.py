#!/Users/frenchd/anaconda2/bin/python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_nice_rotation_curve.py, v1.0 09/27/18

Plot NFW and regular rotation curves all in publication quality plots

This code deals specifically with the non-SALT galaxies - 

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

    # which thing to plot?
    plot_velocity = False
    plot_NFW_fit = True
    plot_rotation_curve = False
    
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
#     galaxyName = 'NGC2770'
#     galaxyName = 'NGC3067'
#     galaxyName = 'NGC3198'
#     galaxyName = 'NGC3351'
#     galaxyName = 'NGC3432'
#     galaxyName = 'NGC3631'
#     galaxyName = 'NGC3666'
#     galaxyName = 'NGC3726'
#     galaxyName = 'NGC4529'
#     galaxyName = 'NGC4565'
#     galaxyName = 'ESO343-G014'
#     galaxyName = 'CGCG039-137'
#     galaxyName = 'IC5325'
#     galaxyName = 'MCG-03-58-009'
#     galaxyName = 'NGC1566'
#     galaxyName = 'NGC3513'
    galaxyName = 'NGC3633'
#     galaxyName = 'NGC4536'
#     galaxyName = 'NGC4939'
#     galaxyName = 'NGC5364'
#     galaxyName = 'NGC5907'
#     galaxyName = 'NGC5951'
#     galaxyName = 'NGC6140'
#     galaxyName = 'NGC5786'
#     galaxyName = 'NGC7817'
#     galaxyName = 'UGC04238'
#     galaxyName = 'UGC06446'
#     galaxyName = 'UGC08146'
#     galaxyName = 'UGC09760'


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
        inc = data['inclination']
        inc_err = data['di']
        vsys_measured = data['vsys_measured']
        vsys_measured_err = data['vsys_measured_err']
        
        
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
        fig = plt.figure(figsize=(7.7, 5.7))

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
                
                
        ylabel(r'$\rm Projected~Rotation~Vel. ~[km~s^{-1}]$')
        xlabel(r'$\rm Vel.~Along~Sightline ~[km~s^{-1}]$')
            
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

        savefig('{0}/{1}-RX_J1121.2+0326_model_plot3.pdf'.format(save_directory, galaxyName),format='pdf',bbox_inches='tight')


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
        x_lim = 3*round(R_vir,0)+10
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
        axvline(x=R_vir, 
                linewidth=3,
                alpha = 0.7,
                c=color_blue,
                label=r'$\rm 1 \emph{R}_{vir}$')
        
#         ax.set_xscale('log')

        
        legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)

        # x-axis
#         majorLocator   = MultipleLocator(25)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(5)
#         ax.yaxis.set_major_locator(majorLocator)
#         ax.yaxis.set_major_formatter(majorFormatter)
#         ax.yaxis.set_minor_locator(minorLocator)

        # y axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
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

        xlim(0, R_vir*3)
#         ylim(0, round(np.nanmax(yData),-2) + 50)

        ylim(0,200)

        savefig('{0}/{1}-NFW_fit_Rvir_times3_2.pdf'.format(save_directory, galaxyName),format='pdf',bbox_inches='tight')


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
        yRight_err = []
        yLeft = []
        yLeft_err = []
        
        for v, e, x in zip(vrot_incCorrected_vals, vrot_incCorrected_errs, xVals):
            xData.append(x)
            yData.append(v)
            yErrs.append(e)
            
            if x > 0:
                xRight.append(x)
                yRight.append(v)
                yRight_err.append(e)
            if x <0:
                xLeft.append(x)
                yLeft.append(v)
                yLeft_err.append(e)
                
            if x == 0:
                if mean(yRight) > 0:
                    if v > 0:
                        xRight.append(x)
                        yRight.append(y)
                        yRight_err.append(e)
                    else:
                        xLeft.append(x)
                        yLeft.append(y)
                        yLeft_err.append(e)
                
                if mean(yRight) <0:
                    if v < 0:
                        xRight.append(x)
                        yRight.append(y)
                        yRight_err.append(e)
                    else:
                        xLeft.append(x)
                        yLeft.append(y)
                        yLeft_err.append(e)
            
            
            
        # turn them into numpy arrays
        xData = np.array(xData)
        yData = np.array(yData)
        yErrs = np.array(yErrs)
    
        # folded over
        xData_abs = abs(xData)
        yData_abs = abs(yData)
        yErrs_abs = abs(yErrs)
    
        # how far to plot? Rounds by 5's
#         x_lim = round(max(abs(xData)/5.),0)*5+5
        x_lim = round(max(abs(xData)/5.),0)*5+1

        # rounds to 50's
        y_lim = round(np.nanmax(yData),-2) + 50
        
        
#         errorbar(xData,
#                 yData,
#                 yerr=yErrs,
#                 fmt='o',
#                 color='black',
#                 elinewidth=1,
#                 ms=6)

        data_plot_label = r'$\rm Data$'
        data_plot = errorbar(xData_abs,
                            yData_abs,
                            yerr=yErrs,
                            fmt='o',
                            color='black',
                            elinewidth=1,
                            ms=6,
                            label=data_plot_label)
        
        
        # plot the average velocities with errors        
        len_right = int(math.ceil(len(yRight)/2.))
        len_left = int(math.floor(len(yLeft)/2.))
        
        # adjust a bit
#         len_right -=2
#         len_left +=2
#         len_right -=5
#         len_left +=5
        
        print
        print 'len_right = {0}, len_left = {1}'.format(len_right, len_left)
        print
            
        # define the outer 1/2 radius y and x values and errors - observed
        x_outer_right = xRight[:len_right]
        x_outer_left = xLeft[len_left:]
        
        y_outer_right = yRight[:len_right]
        y_outer_right_err = yRight_err[:len_right]
        y_outer_left = yLeft[len_left:]
        y_outer_left_err = yLeft_err[len_left:]
        
        # compute the weighted mean of that outside 1/2 rotation velocity        
        y_outer_right_mean = weightedMean(y_outer_right, y_outer_right_err)
        y_outer_left_mean = weightedMean(y_outer_left, y_outer_left_err)
        
        # calculate standard error in the mean for this region
        y_outer_right_sem = stats.sem(y_outer_right)
        y_outer_left_sem = stats.sem(y_outer_left)
        
        # combine inclination + vsys error with this sem error; the inc + vsys errors 
        # should be the same at every point, so just grab the first one
        right_final_err = np.sqrt(y_outer_right_sem**2 + y_outer_right_err[0]**2)
        left_final_err = np.sqrt(y_outer_left_sem**2 + y_outer_left_err[0]**2)

        # try using the weighted errors also
        right_final_werr = np.sqrt(y_outer_right_mean[1]**2 + y_outer_right_err[0]**2)
        left_final_werr = np.sqrt(y_outer_left_mean[1]**2 + y_outer_left_err[0]**2)

        # round the errors off
        right_final_err = int(round(right_final_err,0))
        left_final_err = int(round(left_final_err,0))
        
        right_final_werr = int(round(right_final_werr,0))
        left_final_werr = int(round(left_final_werr,0))
        
        
        print 'x_outer_right: ',len(x_outer_right)
        print
        print 'x_outer_left: ',len(x_outer_left)
        print
        
        
        #####
        # now the overall mean
        if galaxyName == 'NGC3067':
            region_size = 1.5
        elif galaxyName == 'NGC3666':
            region_size = 1.5
        elif galaxyName == 'NGC4529':
            region_size = 1.5
        elif galaxyName == 'NGC6140':
            region_size = 1.4
        elif galaxyName == 'NGC3726':
            region_size = 1.5
        elif galaxyName == 'UGC06446':
            region_size = 1.5
        elif galaxyName == 'UGC08146':
            region_size = 1.5
        elif galaxyName == 'UGC04238':
            region_size = 1.5
        else:
            region_size = 2.
        
        
        x_outer_half_val = (max(xData_abs) - min(xData_abs))/region_size
                
        y_outer_half = []
        x_outer_half = []
        err_outer_half = []
        for x,y,e in zip(xData_abs, yData_abs, yErrs_abs):
            if x >= x_outer_half_val:
                print 'x: ',x
                y_outer_half.append(y)
                x_outer_half.append(x)
                err_outer_half.append(e)
        
        
#         vrot, vrot_err = weightedMean(y_outer_half, x_outer_half)
        print 'y_outer_half : ',y_outer_half
        print
        print 'err_outer_half: ',err_outer_half
        print
        print
        print 'stats.sem(y_outer_half): ',stats.sem(y_outer_half)
        print 'std(y_outer_half): ',std(y_outer_half)
        print
        print 'mean(y_outer_half): ',mean(y_outer_half)
        print
        
        vrot, vrot_err = weightedMean(y_outer_half, err_outer_half)

        
        print 'vrot_err from weightedMean: ',vrot_err
        print
        vrot = int(round(vrot,0))
        vrot_err = int(round(vrot_err,0))
        
#         vrot_err_final = np.sqrt(vrot_err**2 + yErrs[0]**2)
        vrot_err_final = np.sqrt(vrot_err**2 + mean(err_outer_half)**2)

        vrot_err_final = int(round(vrot_err_final,0))
        
        
        print 'right_vrot_incCorrected_avg = ', right_vrot_incCorrected_avg, right_vrot_incCorrected_avg_err
        print 'left_vrot_incCorrected_avg = ', left_vrot_incCorrected_avg, left_vrot_incCorrected_avg_err
        print
        print 'y_outer_right_mean = {0} pm {1}'.format(y_outer_right_mean, y_outer_right_sem)
        print 'y_outer_left_mean = {0} pm {1}'.format(y_outer_left_mean, y_outer_left_sem)
        print
        print 'y_outer_right_mean = {0} pm {1} -> final err'.format(y_outer_right_mean, right_final_err)
        print 'y_outer_left_mean = {0} pm {1} -> final err'.format(y_outer_left_mean, left_final_err)
        print
        print 'y_outer_right_mean = {0} pm {1} -> final werr'.format(y_outer_right_mean, right_final_werr)
        print 'y_outer_left_mean = {0} pm {1} -> final werr'.format(y_outer_left_mean, left_final_werr)
        print
        print 'v_rot = {0} pm {1}'.format(vrot,vrot_err)
        print

        vrot_plot_label = r'$\rm \overline{{\emph{{v}}_{{rot}}}}={0}\pm{1}~km s^{{-1}}$'.format(int(round(vrot, 0)), vrot_err_final)
    
        # plot the average on the right
#         plot((x_outer_right[-1], x_outer_right[0]),
#             (y_outer_right_mean, y_outer_right_mean),
#             lw=2,
#             c='green',
#             solid_capstyle="butt",
#             label=plot_label)
#         
#         # plot the error on the right
#         plot((x_outer_right[-1],
#             x_outer_right[0]),
#             (y_outer_right_mean, y_outer_right_mean),
#             lw=2*right_final_err,
#             c='green',
#             alpha=0.2,
#             solid_capstyle="butt")
#         
#         # plot the average on the left
#         plot((x_outer_left[0],
#             x_outer_left[-1]),
#             (y_outer_left_mean, y_outer_left_mean),
#             lw=2,
#             c='green',
#             solid_capstyle="butt")
#         
#         # plot the error on the left
#         plot((x_outer_left[0],
#                 x_outer_left[-1]),
#                 (y_outer_left_mean, y_outer_left_mean),
#                 lw=2*left_final_err,
#                 c='green',
#                 alpha=0.2,
#                 solid_capstyle="butt")

        # plot the average
        vrot_plot = plot((min(x_outer_half), max(x_outer_half)),
                    (vrot, vrot),
                    lw=2,
                    c='green',
                    solid_capstyle="butt",
                    label = vrot_plot_label)
            
        x_err = np.array([min(x_outer_half), max(x_outer_half)])
        y_err = np.array([vrot, vrot])
        fill_between(x_err, 
                    y_err - vrot_err_final,
                    y_err + vrot_err_final,
                    color='green',
                    alpha=0.2,
                    linewidth=0)
            
        # plot a test for errorbar width
#         plot((max(x_outer_half)+1, max(x_outer_half)+1),
#             (vrot-vrot_err_final, vrot+vrot_err_final),
#             lw=1,
#             c=color_red,
#             alpha=1)


        # x-axis
        majorLocator   = MultipleLocator(4)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2)
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
        
        
#         leg = ax.legend([data_plot, vrot_plot],
#                         ('Data', plot_label),
#                         title=r'$\rm {0}$'.format(galaxyName),
#                         scatterpoints=1,
#                         prop={'size':14},
#                         loc='lower right',
#                         fancybox=True)

        leg = ax.legend(title=r'$\rm {0}$'.format(galaxyName),
                        scatterpoints=1,
                        prop={'size':14},
                        loc='lower right',
                        fancybox=True)

        leg.get_frame().set_alpha(1.0)
        leg._legend_box.align = "left"


#         leg = ax.legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)
        

        ax.grid(b=None,which='major',axis='both')
#         xlim(-x_lim, x_lim)
#         ylim(-(round(np.nanmax(yData),-1) + 25), round(np.nanmax(yData),-1) + 25)

        xlim(0, x_lim)
#         xlim(0, 40.)

        ylim(0, y_lim)
#         ylim(0, 200.)


        savefig('{0}/{1}-rotation_curve_nice_2.pdf'.format(save_directory, galaxyName),format='pdf',bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################

    
if __name__ == '__main__':
    main()
    