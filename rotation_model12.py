#!/Users/frenchd/anaconda2/bin/python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  rotation_model8.py, v8.0 4/27/18

v2: Calculate AGN position from coordinates. (1/30/18)

v3: Plot a cylinder instead of plane (2/07/18)

v4: Make a movie of the sightline piercing the halo (2/09/18)

v5: general updates

v6: add NFW profile fitting (03/28/18)

v7: don't remember

v8: Total rewrite. Rotates a single normal vector now to deal with orientation, works 
    much better. Also added summary files, removed some outdated stuff and added more 
    comments. (04/27/18)

'''


import sys
import os
import csv
import time
from copy import deepcopy
import datetime


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

from astropy.table import Table


from scipy.interpolate import interp1d
from scipy import interpolate
from scipy import interpolate, optimize
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

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
rc('xtick',direction='in')
rc('ytick',direction='in')


'''
========================================================
'''


def stocke_rvir_interpolate(filename):
    """
    Interpolates the Stocke et al. (2013) Figure 1 to get Rvir as a function of L/L*
    
    Parameters
    ----------
    filename    : the directory to the csv file containing an extraction of the plot data
                  Two columns - x and y, both floats
    
    Returns
    -------
    interp      : interpolation object such that interp(L) returns the corresponding Rvir
    
    
    """
    # --- read in the data  and unpack
    log_L, log_rvir = np.loadtxt(filename, delimiter=',', usecols=(0,1), unpack=True)
    
    print('log_L[0]: ',log_L[0])
    print('log_rvir[0]: ',log_rvir[0])
    print()
    
    # --- get out of log space
    L = 10**log_L
    rvir = 10**log_rvir
    
    # --- interpolate it using a cubic spline
    interp = interp1d(L, rvir, kind='cubic')
    
    return interp



def project_model_velocity(z, origin, normal, rayDirection, rayPoint, rcutoff, fit_model, spherical_halo=False, verbose=False, z_gradient=False):
    """
        Returns the projected velocity given the inputs and fit_model 
        
        Parameters
        ----------
        z           : float
                    the z height about the disk midplane
                    
        origin      : np.array [default = [0, 0, 0]]
                    origin point - the center of the galaxy in the midplane
                    
        normal      : np.array
                    galaxy plane normal vector
                    
        rayDirection: np.array [default = [1, 0 ,0]]
                    direction of the QSO vector
                    
        rayPoint    : np.array
                    a point along the QSO ray (expect [0, RA_impact, Dec_impact])
                    
        rcutoff     : float [kpc]
                    the radial extent of the model (how far out to project the disk)
                    
        fit_model   : function of 1 variable
                    the function representing the rotation model. gives rotation
                    velocity as a function of radius
        
        spherical_halo : boolean
                      True means the total distance from the origin is given to fit_model.
                      This results in a halo which is spherical instead of cylindrical
                      
                      False means the distance from the current plane origin is given.
        
        z_gradient  : boolean
                    True results in an exp{-|z|/hv} component added to produce a z-height
                    gradient assuming scale height hv
        
        
        Returns
        -------
        a dictionary containing:
        return {"intersect"         : intersect point of sightline with disk
                "intersect_dist"    : distance from center of plane to intersect
                "dist_from_origin"  : distance from ORIGIN to intersect (equal to intersect_dist for midplane)
                "v_intersect"       : velocity given by model at intersect distance
                "v_n_intersect"     : rotation velocity vector point away from center
                "v_rotation"        : rotation vector correctly oriented wrt sightline
                "v_proj"            : rotation vector projected onto sightline
        
    """
    # --- define a scale height. Only used if z_gradient = True
    h_v = 1000.

    # --- this is a point in the new, parallel but shifted plane
    planePoint = origin + (z * normal)

    # --- get intersect: find_intersect(planeNormal,planePoint,rayDirection,rayPoint)
    intersect = find_intersect(normal, planePoint, rayDirection, rayPoint)
    
    # --- this is the vector from the origin of the current plane to the intersect
    intersect_vect = intersect - (z * normal)
    
    # --- this is the distance from the origin of the current plane to the intersect
    intersect_dist = np.linalg.norm(intersect_vect)
    
    # --- this is the distance from the origin (galaxy center) to the current intersect pt
    dist_from_origin = np.linalg.norm(intersect)

    # --- restrict the intersection to be within the cylinder of radius rcutoff
    if intersect_dist <= rcutoff:
        if verbose:
            print('planePoint: ',planePoint)
            print("intersection at", intersect)
            print('intersect_dist: ',intersect_dist)

        # --- find the rotation velocity at this distance from the rotation curve fit center
        try:
            if spherical_halo:
                v_intersect = fit_model(dist_from_origin)
            else:
                v_intersect = fit_model(intersect_dist)
                if verbose:
                    print('v_intersect: ',v_intersect)
                
                if z_gradient:
#                     v_intersect = v_intersect * np.exp(-abs(z)/h_v)
                    v_intersect = v_intersect * np.exp(-abs(z)/h_v)

        except Exception,e:
            # --- if you go beyond the fit, set velocity to 0
            v_intersect = 0
            if verbose:
                print('Ran out of interpolation range for {0}'.format(intersect_dist))
                print("Built in exception is {0}".format(e))
            sys.exit()
            
        #######
        # --- angle between sightline and vector to intersect point
    
        # --- unit vector towards intersect point
        n_intersect_dist = intersect_vect / np.linalg.norm(intersect_vect)
        if verbose:
            print('n_intersect_dist: ',n_intersect_dist)
            print('np.linalg.norm(n_intersect_dist): ',np.linalg.norm(n_intersect_dist))

    
        # this is the velocity vector in the direction of intersect point, n_intersect_dist
        # edit: seems legit
        v_n_intersect = v_intersect * n_intersect_dist
        if verbose:
            print('new way: ')
            print('v_n_intersect: ',v_n_intersect)
            print('||v_n_intersect|| : ',np.linalg.norm(v_n_intersect))
#         v_list.append(v_n_intersect)


        # v_n_intersect points away from the galaxy center. Rotating this by 90 degrees
        # within the disk gives the correct vector for disk rotation, with amplitude
        # given by the model

        # the result of rotating v counterclockwise by a about n is given by:
        # (cos a)v+(sin a)(n x v)
        #
        # so need to rotate by pi + pi/2 to get all the way around
        alpha = math.pi + math.pi/2
        
        # this then should be the correct rotation velocity vector, but centered at the
        # origin. So, we then need to just shift the sightline to pass through the origin
        #
        # i.e., new sightline = [1, 0 ,0 ] = rayDirection
        v_rotation = math.cos(alpha) * v_n_intersect + math.sin(alpha) * (np.cross(normal, v_n_intersect, axisa=0, axisb=0, axisc=0))
        if verbose:
            print('v_rotation: ',v_rotation)
            print('||v_rotation|| : ',np.linalg.norm(v_rotation))
            print('||N||: ', np.linalg.norm(normal))
#         v_90_list.append(v_rotation)
    
        # now dot it with the sightline to get the component along
        v_proj = np.dot(v_rotation, rayDirection)
        if verbose:
            print('v_proj: ',v_proj)
#         v_proj_list.append(v_proj)
        
#         intersect_list.append(intersect[0])
            print('intersect: ',intersect)
            print('intersect[0]: ',intersect[0])
            print()
#         intersect_point_list.append(intersect)

#         d = -planePoint.dot(normal)
#         d_plot_list.append(d)

        return_dict = {"intersect"          :intersect,
                        "intersect_dist"    :intersect_dist, 
                        "dist_from_origin"  :dist_from_origin,
                        "v_intersect"       :v_intersect,
                        "v_n_intersect"     :v_n_intersect,
                        "v_rotation"        :v_rotation,
                        "v_proj"            :v_proj}
                        
    else:
        v_intersect = False
        v_n_intersect = False
        v_rotation = False
        v_proj = False
        return_dict = {"intersect"          :intersect,
                        "intersect_dist"    :intersect_dist, 
                        "dist_from_origin"  :dist_from_origin,
                        "v_intersect"       :v_intersect,
                        "v_n_intersect"     :v_n_intersect,
                        "v_rotation"        :v_rotation,
                        "v_proj"            :v_proj}
                        
    return return_dict



def adjust(impact,inc,az):
    a = math.sin(az*math.pi/180.) * impact
    
    ds = a * math.tan((90-inc)*math.pi/180.)
    
    return ds
    
    
    
def rotate_vector(v,k,theta):
    '''
        Rodrigues rotation: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        
        Given vector v, this rotates it about k by theta degrees
        
        returns v_rot, the rotated vector
    '''
    
    first_term = v * np.cos(theta)
    
    second_term = (np.cross(k,v,axisa=0, axisb=0, axisc=0)) * np.sin(theta)
    
    third_term = k * (np.dot(k,v)) * (1 - np.cos(theta))
    
    return first_term + second_term + third_term
    

    
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
        
    pass
    

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
        
        

# def NFW(r,a,rho):
#     
#     G = 1.
#     
#     # Hernquist:
# #     M = 4 * math.pi * rho * a**3 * ((r/a)**2) / (2 * (1 + (r/a))**2)
#     M = 4. * math.pi * rho * a**3 * (np.log(1 + r/a) - (r/a) / (1 + (r/a)))
#     
#     v = np.sqrt(G * M / r)
#     
#     return v
    
    
    
def NFW(r,v200,c,r200):

    x = r/r200
    top = (np.log(1 + c*x) - c*x / (1 + c*x))
    
    bottom = x * (np.log(1 + c) - c / (1 + c))
    
    vr = v200 * np.sqrt(top/bottom)
    
    return vr
    
    

def steidel_model(vc, inclination, hv, impact, azimuth, Dlos):
    """
        Calculate v_los based on the model of Steidel et al. (2002). Also used by
        Kacprzak et al. (2019)
        
        Parameters
        ----------
        vc     : float
                maximum rotation velocity (inclination corrected)
        
        inclination: float
                     inclination angle
                
        hv       : float
                scale height of the disk
                
        impact   : float
                impact parameter
                
        azimuth  : float
                azimuth angle
                
        Dlos     : float
                 the distance along the sightline relative to midpoint intersection
        
        Returns
        -------
        vlos    : float
                the vlos at Dlos position along the sightline
        
    """
    
    # --- convert to radians
    to_radians = np.pi/180.
    az_radians = azimuth * to_radians
    inc_radians = inclination * to_radians
    
    y0 = impact * np.sin(az_radians) / np.cos(inc_radians)
    y = (Dlos * np.sin(inc_radians)) + y0
    p = impact * np.cos(az_radians)
    
    exp_term = np.exp(-(abs(y - y0)/(hv * np.tan(inc_radians))))
    
    sqrt_term = np.sqrt(1 + (y/p)**2)
    
#     vlos = -(vc  * np.sin(inc_radians) / sqrt_term) * exp_term
#     vlos = -(vc / sqrt_term) * exp_term
    vlos = (vc  * np.sin(inc_radians) / sqrt_term) * exp_term


    return vlos
    
        
        
        
def plot_NFW(xData, yData, yErr, popt, popt_min, popt_max, x_lim, plotRvir=False):
    """
    This function makes a nice looking plot showing the NFW fit and data
    
    Parameters
    ----------
    xData   : np.array
            the x values for the observed rotation curve 
            
    yData   : np.array
            the y values for the observed rotation curve
    
    yErr    : np.array
            errors for y values
            
    popt    : np.array
            the fit values
    
    popt_min: np.array
            minimum fit values w/ errors

    popt_max: np.array
            maximum fit values w/ errors
            
    x_lim   : float
            the maximum x value to plot to
            
    plotRir : False or float
            if a float, plot a blue vertical line at Rvir. False to not plot anything
    
    Returns
    -------
    fig     : the figure object, to be shown or saved
    
    """
    color_blue = '#436bad'      # french blue
    
#     fig = plt.figure(figsize=(8,8))
    fig = plt.figure(figsize=(7.7,5.7))
    ax = fig.add_subplot(1,1,1)
    
    # --- work out the  fit  errors
    v200, c, r200 = popt
    v200_max, c_max, r200_max = popt_max
    v200_min, c_min, r200_min = popt_min
    
    v200_min_err = abs(v200 - v200_min)
    v200_max_err = abs(v200_max - v200)
    v200_err = max(v200_min_err, v200_max_err)
    
    c_min_err = abs(c - c_min)
    c_max_err = abs(c_max - c)
    c_err = max(c_min_err, c_max_err)
    
    r200_min_err = abs(r200 - r200_min)
    r200_max_err = abs(r200_max - r200)
    r200_err = max(r200_min_err,r200_max_err)
    
    v200_label = r'$\rm V200={0}~\pm~{1}$'.format(round(v200,1), round(v200_err,1))
    c_label = r'$\rm c={0}~\pm~{1}$'.format(round(c,1), round(c_err,1))
#     r200_label = r'$\rm R200={0}~\pm~{1}$'.format(round(r200,1), round(r200_err,1))
    r200_label = r'$\rm R200={0}$'.format(round(r200,1))

    x_fit = linspace(0,x_lim,num=1000)
    
#     scatter(xData, yData, color='black', s=40, lw=0, label = r'$\rm Data$')
    errorbar(xData, yData, yerr=yErr, marker='o', color='black', markersize=4, linewidth=0, elinewidth=0.5, label = r'$\rm Data$')

#     plot(x_fit, NFW(x_fit, *popt), 'r-',label='fit: a={0}, rho={1}'.format(*popt))
    plot(x_fit, NFW(x_fit, *popt), 'r-', color='green', alpha=0.7, lw=2, \
    label='NFW Fit:\n {0}\n{1}\n{2}'.format(v200_label, c_label, r200_label))
    
    if plotRvir:
        axvline(plotRvir, ymin=0, ymax=1, color=color_blue, alpha=0.85, lw=2, \
        label=r'$\rm R_{vir}$')
    
    ax.fill_between(x_fit, NFW(x_fit, *popt_min), NFW(x_fit, *popt_max), facecolor='green', alpha=0.2)
    
    legend(scatterpoints=1,prop={'size':11},loc='lower right',fancybox=True)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    
    xlim(0, x_lim)
    ylim(0, round(np.nanmax(yData),-1) + 25)
#     grid(True, alpha=0.5)

    if x_lim < 50:
        # x-axis
        majorLocator   = MultipleLocator(1)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(0.2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
    else:
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(20)
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
    ylabel(r'$\rm \emph{v}_{{rot}} ~[km~s^{{-1}}]$')
    
    return fig
    
    

def writeout_fits(xData, yData, popt, x_lim, filename):
    '''
    This function creates a text file containing the rotation curve data and
    fit parameters for both NFW and cylindrical fits
    
    xData - the x values for the observed rotation curve 
    yData - the y values for the observed rotation curve
    popt - the fit values
    x_lim - the maximum x value to plot to
    
    '''
    
    v200, c, r200 = popt
    
    t = Table()
    t.meta['comment'] = ['Fit info: ',
    'NFW Fit:','V200={0}, c={1}, R200={2}'.format(round(v200,2),round(c,2),round(r200,2))]
    
    pass



def my_round(x, prec=2, base=0.05):
    return (base * (np.array(x) / base).round()).round(prec)


def main(galaxyName):
    print("Beginning {0}".format(galaxyName))
    print()
    
    hubbleConstant = 71.0
    
    # fit an NFW profile to the rotation curve and use that? Otherwise it takes the maximal
    # rotation value and extends it forever. NFW decreases with distance
    fit_NFW = True
    
    # Makes the halo 3x4 R_vir if true. 2x3R_vir if false
    extendedHalo = False
    
    # only us the plane of the galaxy (no z-direction, or height, to the disk)
    diskOnly = False
    
    # --- Halo size:
    # --- height multiplier (* R_vir)
    zcutoffm = 3.0
    
    # --- radius multiplier (* R_vir)
    rcutoffm = 3.0
    
    # --- spherical NFW  halo? (NFW  curve applies in all directions, not just radially)
    spherical_halo = False
    
    # --- apply a exp(-|z|/hv) z gradient to  the  NFW model? Either false or a number
    # --- corresponding to h_v
    z_gradient = False
    
    # --- hv for the Steidel model
    steidel_hv = 1000.0
    
    # save the full model velocity data in a pickle file?
    save_data_pickle = True
    
    # --- save a csv file with all this in it
    save_results_csv = True
    
    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red

    # --- interpolate the Stocke et al. (2013) Lstar vs Rvir relation
    stocke_rvir_filename = '/Users/frenchd/Research/rotation_paper_data/Rvir_L_Lstar2.csv'
    stocke_rvir = stocke_rvir_interpolate(stocke_rvir_filename)
    
    # ---  Use the Stocke + 2013 Rvir calculation or the usual one?
    use_stocke_rvir = True
    
    # --- center at the midplane intersect? centers at galaxy systemic if false
    center_at_intersect = False
    
    # --- position angle error for all galaxies
    e_PA = 3.0


#     galaxyName = 'CGCG039-137'
#     galaxyName = 'ESO343-G014'
#     galaxyName = 'IC5325'
#     galaxyName = 'MCG-03-58-009'
#     galaxyName = 'NGC1566'
#     galaxyName = 'NGC3513'
#     galaxyName = 'NGC3633'
#     galaxyName = 'NGC4536'
#     galaxyName = 'NGC4939'
#     galaxyName = 'NGC5364'
#     galaxyName = 'NGC5786'
#     galaxyName = 'UGC09760'

#     galaxyName = 'NGC3198'
#     galaxyName = 'NGC4565'
#     galaxyName = 'NGC3351'
#     galaxyName = 'UGC04238'
#     galaxyName = 'NGC4529'
#     galaxyName = 'NGC6140'
#     galaxyName = 'NGC5907'
#     galaxyName = 'UGC06446'
#     galaxyName = 'NGC3631'
#     galaxyName = 'UGC06399'
#     galaxyName = 'NGC3726'
#     galaxyName = 'NGC3067'
#     galaxyName = 'NGC2770'
#     galaxyName = 'NGC3432'
#     galaxyName = 'NGC3666'
#     galaxyName = 'NGC5951'
#     galaxyName = 'NGC7817'
#     galaxyName = 'UGC08146'
#     galaxyName = 'M31'

#     galaxyName = 'z0_3528gal'


#     saveDirectory = '/Users/frenchd/Research/M31_rotation/{0}/model_v3/'.format(galaxyName)
#     saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/models_v2/{0}/'.format(galaxyName)
#     saveDirectory = '/Users/frenchd/Research/rotation_paper_data/rotation_models_v2/{0}/'.format(galaxyName)
    saveDirectory = '/Users/frenchd/Research/rotation_paper_data/rotation_models_v4/{0}/'.format(galaxyName)


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
#     directory = '/Users/frenchd/Research/M31_rotation/'
#     directory = '/Users/frenchd/Research/rotation_paper_data/summary_files_4/'
#     csv_filename = '{0}/rotation_model_summary_4.csv'.format(directory)
    directory = '/Users/frenchd/Research/rotation_paper_data/summary_files_4/'
    csv_filename = '{0}/rotation_model_summary_plots.csv'.format(directory)

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
#     filename = '{0}-summary6.json'.format(galaxyName)
#     filename = '{0}-summary7.json'.format(galaxyName)
    filename = '{0}-summary.json'.format(galaxyName)

    
    with open(directory+filename) as data_file:
        data = json.load(data_file)    
        
        vrot_vals = data['vrot_vals']
        vrot_incCorrected_vals = data['vrot_incCorrected_vals']
        right_vrot_incCorrected_avg = data['right_vrot_incCorrected_avg']
        left_vrot_incCorrected_avg = data['left_vrot_incCorrected_avg']
        
        # errors
        vrot_incCorrected_errs = data['vrot_incCorrected_errs']
        
        left_vrot_incCorrected_avg_err = data['left_vrot_incCorrected_avg_err']
        right_vrot_incCorrected_avg_err = data['right_vrot_incCorrected_avg_err']

        left_vrot_avg_err = data['left_vrot_avg_err']
        right_vrot_avg_err = data['right_vrot_avg_err']

        
        xVals = data['xVals']
        inc = data['inclination']
        e_inc = data['di']
        vsys_measured = data['vsys_measured']
        vsys_measured_err = data['vsys_measured_err']
        vsys_published = data['vsys_published']
        dvsys_published = data['dvsys_published']
        
#         inc = 40.0
        
#         galaxyName = data['name']
        RA_galaxy = data['RAdeg']
        Dec_galaxy = data['DEdeg']
        dist = data['dist']
        majDiam = data['majDiam']
        PA = data['PA']
        agn = data['agn']
        
        # --- Lstar and error calculation
        Lstar = float(data['Lstar'])
        e_Lstar = float(data['e_Lstar'])

        
        # --- calculate a few things
        if use_stocke_rvir:
            R_vir = stocke_rvir(Lstar)
#             e_R_vir_up = stocke_rvir(Lstar + e_Lstar) - R_vir
#             e_R_vir_down = R_vir - stocke_rvir(Lstar - e_Lstar)
            e_R_vir_up = R_vir + (R_vir*0.1)
            e_R_vir_up = R_vir - (R_vir*0.1)

            R_vir_old = calculateVirialRadius(majDiam)
        else:
            R_vir = calculateVirialRadius(majDiam)
            e_R_vir_up = calculateVirialRadius(majDiam + majDiam*0.1)
            e_R_vir_down = calculateVirialRadius(majDiam - majDiam*0.1)

        # --- define "vmax", the average maximum rotation velocity
        vmax_incCorrected = max(abs(right_vrot_incCorrected_avg), abs(left_vrot_incCorrected_avg))
        
        # -- remove inclination correction for the Steidel model method
        vmax = vmax_incCorrected * np.sin(inc * np.pi/180.)

        print()
        print('PA: ',PA)
        print('inc: ',inc)
        print('dist: ',dist)
        print('AGN: ',agn)
        print()
        print('Lstar: ',Lstar)
        print('stocke_rvir(Lstar) = ',stocke_rvir(Lstar))
        print()


        # which agn do you want to target? Also decide here if you want to mirror around
        # the inclination axis (i.e., so the galaxy is "facing" left vs right for PA=0)
        # 
        # reverse = True to reverse the rotation direction
        #
        # NFW_fit decides how tightly to bound the NFW profile fit. 
        # Options are: standard, tight, tighter, tightest, tightester
        
        NFW_fit = "standard"
        
        # M31
        if galaxyName == 'M31':
            flipInclination = False
            # --- reverse for NFW
            if fit_NFW:
                reverse = True
            else:
                reverse = False
            agnName = 'QSO1'
        

        # CGCG039-137
        if galaxyName == 'CGCG039-137':
            flipInclination = False
            reverse = False
#             agnName = 'RX_J1121.2+0326'
            agnName = 'SDSSJ112224.10+031802.0'


        # RFGC3781 or ESO343-G014
        if galaxyName == 'RFGC3781' or galaxyName == 'ESO343-G014':
            flipInclination = False
            reverse = False
            agnName = 'RBS1768'
        
        
        # IC5325
        if galaxyName == 'IC5325':
            flipInclination = False
            reverse = False
            agnName = 'RBS2000'
        
        
        # MCG-03-58-009
        if galaxyName == 'MCG-03-58-009':
            flipInclination = False
            
            # --- reverse for 2x3R_vir, not for NFW
            if fit_NFW:
                reverse = False
            else:
                reverse = True
            agnName = 'MRC2251-178'
        
        
        # NGC1566
        if galaxyName == 'NGC1566':
            flipInclination = True
            reverse = False
            # extendedHalo for HE0439-5254, RBS567
            # extendedHalo = True
#             agnName = 'HE0429-5343'
#             agnName = '1H0419-577'
            agnName = 'HE0435-5304'
#             agnName = 'RBS567'
#             agnName = 'HE0439-5254'


        # NGC3513
        if galaxyName  == 'NGC3513':
            flipInclination = False
            # reverse true for NFW, false for cylindrical
            if fit_NFW:
                reverse = True
            else:
                reverse = False
            agnName = 'H1101-232'

        
        # NGC3633
        if galaxyName == 'NGC3633':
            flipInclination = False
            reverse = False
            agnName = 'RX_J1121.2+0326'
#             agnName = 'SDSSJ112224.10+031802.0'
#             agnName = 'SDSSJ112005.00+041323.0'


        # NGC4536
        if galaxyName == 'NGC4536':
            # --- use NGC4536-summary7.json for NFW fits (this file has the weird data on the 
            # --- (left of NGC4536 removed, resulting in a much better NFW fit)
            # --- NGC4536-summary6.json has all the original data intact
            flipInclination = False
            reverse = False
#             agnName = '3C273.0'
            agnName = 'HE1228+0131'


        # NGC4939
        if galaxyName == 'NGC4939':
            flipInclination = False
            reverse = True
            agnName = 'PG1302-102'


        # NGC5364
        if galaxyName == 'NGC5364':
            flipInclination = False
            # --- reverse for 2x3R_vir, not reverse for NFW
            if fit_NFW:
                reverse = False
            else:
                reverse = True
            agnName = 'SDSSJ135726.27+043541.4'

        
        # NGC5786
        if galaxyName == 'NGC5786':
            flipInclination = True
            reverse = False
            agnName = 'QSO1500-4140'
        
        
        # UGC09760
        if galaxyName == 'UGC09760':
            flipInclination = True
            # --- reverse for 2x3R_vir, not reverse for NFW
            if fit_NFW:
                reverse = False
            else:
                reverse = True
            agnName = 'SDSSJ151237.15+012846.0'

##########################################################################################
        #  --- These are all the galaxies with data from the literature
##########################################################################################

        # NGC2770
        if galaxyName == 'NGC2770':
            flipInclination = True
            NFW_fit = 'tighter'
            # --- reverse for NFW, not for 2x3R_vir
            if fit_NFW:
                reverse = True
            else:
                reverse = False

#             agnName = 'FBQSJ0908+3246'
#             agnName = 'TON1009'
            agnName = 'TON1015'
#             agnName = 'SDSSJ091052.80+333008.0'
#             agnName = 'SDSSJ091127.30+325337.0'


        # NGC3067
        if galaxyName == 'NGC3067':
            flipInclination = False
            reverse = True
#             agnName = '3C232'
#             agnName = 'RX_J1002.9+3240'
            agnName = 'SDSSJ095914.80+320357.0'


        # NGC3198
        if galaxyName == 'NGC3198':
            flipInclination = True
            reverse = True
#             agnName = 'RX_J1017.5+4702'
            agnName = 'SDSSJ101622.60+470643.0'


        # NGC3351
        if galaxyName == 'NGC3351':
            flipInclination = True
            reverse = True
            agnName = 'SDSSJ104335.90+115129.0'
            
            # --- no data for any of these
#             agnName = 'SDSSJ104341.53+085558.2'
#             agnName = 'SDSSJ104709.80+130454.0'
#             agnName = 'SDSSJ104816.30+120735.0'
#             agnName = 'SDSSJ104843.50+130605.0'
#             agnName = 'SDSSJ105220.60+101751.0'


        # NGC3432
        if galaxyName == 'NGC3432':
            flipInclination = True
            NFW_fit = 'tighter'
            # --- reverse for NFW, not for 2x3R_vir
            if fit_NFW:
                reverse = True
            else:
                reverse = False

#             agnName = 'MS1047.3+3518'
#             agnName = 'CSO295'
            agnName = 'RX_J1054.2+3511'


        # NGC3631
        if galaxyName == 'NGC3631':
            flipInclination = False
            reverse = False
#             agnName = 'SDSSJ111443.70+525834.0'
#             agnName = 'RX_J1117.6+5301'
            agnName = 'SBS1116+523'
#             agnName = 'SDSSJ112448.30+531818.0'


        # NGC3666
        if galaxyName == 'NGC3666':
            flipInclination = True
            NFW_fit = 'tighter'
            # --- reverse for NFW, not for 2x3R_vir
            if fit_NFW:
                reverse = True
            else:
                reverse = False

            agnName = 'SDSSJ112439.50+113117.0'
            # --- these are probably no good
#             agnName = 'SDSSJ112632.90+120437.0'
#             agnName = 'SDSSJ112756.70+115427.0'


        # NGC3726
        if galaxyName == 'NGC3726':
            flipInclination = True
            reverse = True
            NFW_fit = 'standard'
#             agnName = 'CSO1208'
            agnName = 'RX_J1142.7+4625'


        # NGC4529
        if galaxyName == 'NGC4529':
            NFW_fit = 'tightester'
            flipInclination = False
            reverse = True
            agnName = 'MRK771'


        # NGC4565
        if galaxyName == 'NGC4565':
            flipInclination = False
            reverse = False
            agnName = 'RX_J1236.0+2641'


        # NGC5907
        if galaxyName == 'NGC5907':
            flipInclination = False
            reverse = False
#             agnName = 'SBS1503+570'
#             agnName = 'SDSSJ152053.59+571122.1'
            agnName = 'RBS1503'


        # NGC5951
        if galaxyName == 'NGC5951':
            flipInclination = False
            reverse = True
            NFW_fit = 'tightest'
            agnName = '2E1530+1511'


        # NGC6140
        if galaxyName == 'NGC6140':
            flipInclination = False
            # reverse for 2x3R_vir, not for NFW
            if fit_NFW:
                reverse = False
            else:
                reverse = True
                
            agnName = 'MRK876'
#             agnName = 'KAZ49' - don't bother
#             agnName = 'HS1626+6433' - don't bother


        # NGC7817
        if galaxyName == 'NGC7817':
            flipInclination = True
            NFW_fit = 'standard'
            # reverse for NFW, not for 2x3R_vir
            if fit_NFW:
                reverse = True
            else:
                reverse = False
                
            agnName = 'MRK335'
#             PA = 43.


        # UGC04238
        if galaxyName == 'UGC04238':
            # --- need to set tighter NFW fit bounds for this one
            NFW_fit = 'tightester2'
            flipInclination = False
            # --- reverse for 2x3R_vir, not for NFW
            if fit_NFW:
                reverse = False
            else:
                reverse = True
                
            agnName = 'PG0804+761'


        # UGC06399
        if galaxyName == 'UGC06399':
            flipInclination = False
            reverse = True
            agnName = 'SBS1116+523'
        
        
        # UGC06446
        if galaxyName == 'UGC06446':
            flipInclination = False
            reverse = False
            agnName = 'SDSSJ112448.30+531818.0'
            # --- too far -> 416.788 kpc
#             agnName = 'RX_J1117.6+5301'


        # UGC08146
        if galaxyName == 'UGC08146':
            flipInclination = False
            NFW_fit = 'tightester2'
            reverse = True
            agnName = 'PG1259+593'
            

            
        # Test with the z=0.3528 galaxy from Kacprzak 2019 (last figure in appendix)
        if galaxyName == 'z0_3528gal':
            print 'here'
            flipInclination = False
            reverse = False
            agnName = 'z0_3528'

        # --- grab the coordinates for this target
        RA_target = agn[agnName]['RAdeg']
        Dec_target = agn[agnName]['DEdeg']
        
        # --- approaching or receding side?
        side = agn[agnName]['side']

        # --- check if the sign is right for vmax
#         if side == 'receding':
        if side == 'approaching':
            vmax = float(vmax) * -1

        #  --- test which side to use error from
        
        if vmax > 0:
            if right_vrot_incCorrected_avg > 0:
                e_vmax_incCorrected = right_vrot_incCorrected_avg_err
                print('right side')
            elif left_vrot_incCorrected_avg > 0:
                e_vmax_incCorrected = left_vrot_incCorrected_avg_err
                print('left side')
            else:
                print('wtf?')

        elif vmax < 0:
            if right_vrot_incCorrected_avg < 0:
                e_vmax_incCorrected = right_vrot_incCorrected_avg_err
                print('right side')
            elif left_vrot_incCorrected_avg < 0:
                e_vmax_incCorrected = left_vrot_incCorrected_avg_err
                print('left side')
            else:
                print('wtf?')
        else:
            print('wtf?')
        
        e_vmax = e_vmax_incCorrected * np.sin(inc * np.pi/180.)
        
        print('vmax:  ',vmax)
        print('e_vmax: ',e_vmax)
        print()
        print('vmax_incCorrected: ',vmax_incCorrected)
        print('e_vmax_incCorrected: ',e_vmax_incCorrected)
        print()

        # --- append target to save directory
        saveDirectory = saveDirectory + str(agnName) + '/'
        
##########################################################################################
##########################################################################################
    # --- rename arrays something stupid
    vels = vrot_incCorrected_vals
    errs = vrot_incCorrected_errs
    xvals = deepcopy(xVals)
    vsys = vsys_measured
    
    if fit_NFW:
    
        # --- fold the data over so approaching and receding sides are both positive
        newVals = []
        newX = []
        
        outer_x = []
        outer_vels = []
        
        # --- these will be the 
        xData1 = []
        yData1 = []
        xData2 = []
        yData2 = []
        
#         print 'vels: ',vels
#         print
#         print 'xVals: ',xVals
#         print
#         print
        
        # --- stack, take absolute value and sort the x and y values
        stacked = np.dstack((xVals, vels, errs))[0]
        stacked = abs(stacked)
        stacked = stacked[stacked[:,0].argsort()]
        
        new_xs = stacked.T[0]
        new_vels = stacked.T[1]
        new_errs = stacked.T[2]

        # --- define fit initial conditions
#         a = 3.95
#         rho = 500.
        
        v200 = abs(vmax_incCorrected) / 1.1
        c = 5
        r200 = R_vir
        
        r200_lowbound = R_vir/1.2
        r200_upbound = R_vir*1.2
        v200_lowbound = abs(vmax_incCorrected)/3.0
        v200_upbound = abs(vmax_incCorrected)*1.2
        c_lowbound = 1
        c_upbound = 50
        
        r200_lowbound_max = R_vir/1.2
        r200_upbound_max = R_vir*1.2
        v200_lowbound_max = abs(vmax_incCorrected + e_vmax_incCorrected)/3.0
        v200_upbound_max = abs(vmax_incCorrected + e_vmax_incCorrected)*1.2
        c_lowbound_max = 1
        c_upbound_max = 50
        v200_max = abs(vmax_incCorrected + e_vmax_incCorrected) / 1.2

        
        r200_lowbound_min = R_vir/1.2
        r200_upbound_min = R_vir*1.2
        v200_lowbound_min = abs(vmax_incCorrected - e_vmax_incCorrected)/3.0
        v200_upbound_min = abs(vmax_incCorrected - e_vmax_incCorrected)*1.2
        c_lowbound_min = 1
        c_upbound_min = 50
        v200_min = abs(vmax_incCorrected - e_vmax_incCorrected) / 1.2

#         print 'r200_lowbound: ',r200_lowbound
#         print 'r200_upbound: ',r200_upbound
#         print 'v200_lowbound: ',v200_lowbound
#         print 'v200_upbound: ',v200_upbound
#         print 'R_vir: ',R_vir
#         print
#         print 'r200_lowbound_min: ',r200_lowbound_min
#         print 'r200_upbound_min: ',r200_upbound_min
#         print 'v200_lowbound_min: ',v200_lowbound_min
#         print 'v200_upbound_min: ',v200_upbound_min
#         print
#         print 'r200_lowbound_max: ',r200_lowbound_max
#         print 'r200_upbound_max: ',r200_upbound_max
#         print 'v200_lowbound_max: ',v200_lowbound_max
#         print 'v200_upbound_max: ',v200_upbound_max
#         print
        
        
        try:
            # --- grab just the outer fraction of the curve (don't fit the center region)
            fit_fraction = 5
            fit_xs = new_xs[int(len(new_xs)/fit_fraction):]
            fit_vels = new_vels[int(len(new_vels)/fit_fraction):]
            fit_errs = new_errs[int(len(new_errs)/fit_fraction):]
            
            print('fit_xs: ',fit_xs)
            print()
            print('new_xs: ',new_xs)
            print()
            print('fit_vels: ',fit_vels)
            print()
            print('new_vels: ',new_vels)
            print()
            
##############################################################################################
##############################################################################################
            def NFW_2(r, v200, c):
            
                r200 = R_vir

                x = r/r200
                top = (np.log(1 + c*x) - c*x / (1 + c*x))
    
                bottom = x * (np.log(1 + c) - c / (1 + c))
    
                vr = v200 * np.sqrt(top/bottom)
    
                return vr


            def ff(r,p):
                # NFW(r,v200,c,r200)
                return NFW_2(r, *p)
            
            
            def fit_bootstrap(p0, datax, datay, function, bounds, yerr_systematic=0.0):
                
                errfunc = lambda p, x, y: function(x,p) - y

                # Fit first time
                pfit, perr = optimize.leastsq(errfunc, p0, args=(datax, datay), full_output=0)
#                 pfit, perr = optimize.fmin_slsqp(errfunc, p0, args=(datax, datay), bounds=bounds)

                pfit2, perr2 = optimize.curve_fit(NFW_2,
                                                datax,
                                                datay,
                                                p0=p0,
                                                sigma=fit_errs,
                                                bounds=bounds,
                                                absolute_sigma=False)


                # Get the stdev of the residuals
#                 residuals = errfunc(pfit, datax, datay)
#                 residuals = errfunc(datax, datay, pfit)
                residuals = errfunc(pfit2, datax, datay)

                sigma_res = np.std(residuals)

                sigma_err_total = np.sqrt(sigma_res**2 + yerr_systematic**2)

                # 100 random data sets are generated and fitted
                ps = []
                for i in range(1000):

                    randomDelta = np.random.normal(0., sigma_err_total, len(datay))
                    randomdataY = datay + randomDelta

                    randomfit, randomcov = \
                            optimize.curve_fit(NFW_2,
                                            datax,
                                            randomdataY,
                                            p0=p0,
                                            sigma=fit_errs,
                                            bounds=bounds,
                                            absolute_sigma=False)

#                     randomfit, randomcov = \
#                         optimize.leastsq(errfunc, p0, args=(datax, randomdataY),\
#                                         full_output=0)

#                     randomfit, randomcov = \
#                         optimize.fmin_slsqp(errfunc, p0, args=(datax, randomdataY),\
#                                         bounds=bounds)

                    ps.append(randomfit)

                ps = np.array(ps)
                mean_pfit = np.mean(ps,0)

                # You can choose the confidence interval that you want for your
                # parameter estimates: 
                Nsigma = 1. # 1sigma gets approximately the same as methods above
                            # 1sigma corresponds to 68.3% confidence interval
                            # 2sigma corresponds to 95.44% confidence interval
                err_pfit = Nsigma * np.std(ps,0) 

                pfit_bootstrap = mean_pfit
                perr_bootstrap = err_pfit
                return pfit_bootstrap, perr_bootstrap
                
##############################################################################################
##############################################################################################
            # --- old method of fitting NFW curve. Errors don't work
            
#             popt, pcov = optimize.curve_fit(NFW,
#                                             fit_xs,
#                                             fit_vels,
#                                             p0=[v200,c,r200],
#                                             sigma=fit_errs,
#                                             absolute_sigma=False,
#                                             bounds=((v200_lowbound, c_lowbound, r200_lowbound),
#                                                     (v200_upbound, c_upbound, r200_upbound)))
# 
# 
#             popt_max, pcov_max = optimize.curve_fit(NFW,
#                                                 fit_xs,
#                                                 fit_vels,
#                                                 p0=[v200_max,c,r200],
#                                                 sigma=fit_errs,
#                                                 absolute_sigma=False,
#                                                 bounds=((v200_lowbound_max, c_lowbound_max, r200_lowbound_max),
#                                                         (v200_upbound_max, c_upbound_max, r200_upbound_max)))
# 
#             popt_min, pcov_min = optimize.curve_fit(NFW,
#                                                 fit_xs,
#                                                 fit_vels,
#                                                 p0=[v200_min,c,r200],
#                                                 sigma=fit_errs,
#                                                 absolute_sigma=False,
#                                                 bounds=((v200_lowbound_min, c_lowbound_min, r200_lowbound_min),
#                                                         (v200_upbound_min, c_upbound_min, r200_upbound_min)))

            
            pstart = [v200, c]
            bounds = ([v200_lowbound, c_lowbound], [v200_upbound, c_upbound])
#             bounds = ([v200_lowbound, v200_upbound], [c_lowbound, c_upbound])

            pfit, perr = fit_bootstrap(pstart, fit_xs, fit_vels, ff, bounds)
            
            print('pstart: ',pstart)
            print('bounds: ',bounds)
            print()
            print('alternative bootstrap method:')
            print('pfit: ',pfit)
            print('perr: ',perr)
            print()
            
            popt = [pfit[0],pfit[1],R_vir]
            popt_max = [pfit[0]+perr[0] ,pfit[1]+perr[1], R_vir]
            popt_min = [pfit[0]-perr[0] ,pfit[1]-perr[1], R_vir]
            
            
            # --- work out the  fit  errors
            v200, c, r200 = popt
            v200_max, c_max, r200_max = popt_max
            v200_min, c_min, r200_min = popt_min
    
            v200_min_err = abs(v200 - v200_min)
            v200_max_err = abs(v200_max - v200)
            v200_err = max(v200_min_err, v200_max_err)
    
            c_min_err = abs(c - c_min)
            c_max_err = abs(c_max - c)
            c_err = max(c_min_err, c_max_err)
    
            r200_min_err = abs(r200 - r200_min)
            r200_max_err = abs(r200_max - r200)
            r200_err = max(r200_min_err,r200_max_err)
            
            
        except Exception,e:
            print 'exception in curve_fit: ',e
            print
            sys.exit()
        
#         print
#         print 'now the fit. popt = {0}, pcov = {1}'.format(popt, pcov)
#         print
#         print 'Max popt = {0}, pcov = {1}'.format(popt_max, pcov_max)
#         print
#         print 'Min popt = {0}, pcov = {1}'.format(popt_min, pcov_min)
#         print
#         print 'np.sqrt(np.diag(pcov)) = ',np.sqrt(np.diag(pcov))
#         print

        # plot it
#         xData_fit = linspace(0,max(xvals),num=1000)
#         y_fit = NFW(xData_fit,*popt)
        x_lim = int(max(new_xs)) + int(max(new_xs)/2.)
        print('x_lim: ',x_lim)
        print()
        fig = plot_NFW(new_xs, new_vels, new_errs, popt, popt_min, popt_max, x_lim)
#         fig.savefig("{0}{1}_NFW_{2}.jpg".format(saveDirectory,galaxyName,x_lim),dpi=300,bbox_inches='tight')
        fig.savefig("{0}{1}_NFW_{2}.pdf".format(saveDirectory,galaxyName, int(x_lim)), bbox_inches='tight')

        x_lim = round(3*R_vir + 5,0)
#         x_lim = round(R_vir + 10,0)
        fig = plot_NFW(new_xs, new_vels, new_errs, popt, popt_min, popt_max, x_lim, plotRvir=R_vir)
#         fig.savefig("{0}{1}_NFW_{2}.jpg".format(saveDirectory,galaxyName,x_lim),dpi=300,bbox_inches='tight')
        fig.savefig("{0}{1}_NFW_{2}.pdf".format(saveDirectory,galaxyName, int(x_lim)),bbox_inches='tight')


        def fit(x):
            '''
                this function is necessary to deal with the negative side. It takes in
                the requested x value and decides whether to flip sign of the resulting
                NFW velocity
                
                e.g., if x < 0: y = y* -1
            '''

            y_val = 0
            if x >= 0:
                y_val = NFW(x, *popt)

            else:
                y = NFW(abs(x), *popt)
                
                y_val = -y
            
            if reverse:
                y_val = -y_val

            return y_val
            

        def fit_min(x):
            '''
                this function is necessary to deal with the negative side. It takes in
                the requested x value and decides whether to flip sign of the resulting
                NFW velocity
                
                e.g., if x < 0: y = y* -1
            '''

            y_val = 0
            if x >= 0:
                y_val_min = NFW(x, *popt_min)

            else:
                y_val_min = NFW(abs(x), *popt_min)
                
                y_val_min = -y_val_min
            
            if reverse:
                y_val_min = -y_val_min

            return y_val_min


        def fit_max(x):
            '''
                this function is necessary to deal with the negative side. It takes in
                the requested x value and decides whether to flip sign of the resulting
                NFW velocity
                
                e.g., if x < 0: y = y* -1
            '''

            y_val = 0
            if x >= 0:
                y_val_max = NFW(x, *popt_max)

            else:
                y_val_max = NFW(abs(x), *popt_max)
                
                y_val_min = -y_val_min
            
            if reverse:
                y_val_max = -y_val_max

            return y_val_max


    else:
    
        # this one is for 0 centered
        xvalStart = xvals[0]
        xvalEnd = xvals[-1]
 
        lowMean = mean(vels[:6])
        highMean = mean(vels[-6:])
        
        step = 5
        for i in range(250):
            vels.insert(0,right_vrot_incCorrected_avg)
            vels.append(left_vrot_incCorrected_avg)

            xvalStart +=step
            xvalEnd -=step
            xvals.insert(0,xvalStart)
            xvals.append(xvalEnd)
            
        # reverse it?
        if reverse:
    #     xData.reverse()
            vels.reverse()
            
        splineKind = 'linear'
#         splineKind = 'cubic'

        
        fit = interp1d(xvals, vels, kind=splineKind)
        popt = False
        
#         fig = plt.figure(figsize=(7,5))
        fig = plt.figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(1,1,1)
        
        x_fit = linspace(-50, 50, num=1000)

        scatter(xvals, vels, color='black', s=40, label = r'$\rm Data$')
        plot(x_fit, fit(x_fit), color='green', alpha=0.7, ms=0, lw=2, label = r'$\rm Linear-Spline ~Fit$')
        legend()
        
        # x-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        xlim(-50, 50)
        xlabel(r'$\rm R~[kpc]$')
        ylabel(r'$\rm V_{rot}$')
        
        fig.savefig("{0}{1}_splineFit_{2}.jpg".format(saveDirectory,galaxyName,splineKind),dpi=300,bbox_inches='tight')
        
##########################################################################################
##########################################################################################
##########################################################################################
    zcutoff = zcutoffm * R_vir
    rcutoff = rcutoffm * R_vir
    print('zcutoff: ',zcutoff)
    print()
    print('PA = {}, e_PA = {}'.format(PA, e_PA))
    print()

    pa_range = [(PA-e_PA), PA, (PA+e_PA)]
    results_dict = {}
    
    for pa in pa_range:
        # end of input prep #
        #####################
        # do some calculations now
        
#         RA_galaxy = 343.49872000
#         Dec_galaxy = 16.1482100
#         
#         RA_target = 343.49062000
#         Dec_target = 16.1485200

        # --- calculate impact parameter and shit
        impact = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist)
    
        # --- RA component of impact parameter - by setting the Dec to be the same for both
        impact_RA = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_galaxy,dist)
    
        # --- Dec component of impact parameter - by setting the RA to be the same for both
        impact_Dec = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_galaxy,Dec_target,dist)
    
        # --- calculate azimuth
        az = calculateAzimuth(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist,pa)
        
        if Dec_galaxy > Dec_target:
            impact_Dec = -impact_Dec

        if RA_target > RA_galaxy:
            impact_RA = -impact_RA
            
        #######################
        # --- testing....
#         impact = 203.2
#         impact_RA = 202.954
#         impact_Dec = 5.0
#         az = 88.7
        #######################


        # --- inclination is backwards, so flip it
    #     effectiveInc = 90.-inc
        effectiveInc = inc
        print('effectiveInc: ',effectiveInc)
        print()
    
        if flipInclination:
            effectiveInc *=-1.
#             impact_RA *=-1

        print()
        print('PA: ',pa)
        print('impact: ',impact)
        print('impact_RA: ',impact_RA)
        print('impact_Dec: ',impact_Dec)
        print('az: ',az)
        print('R_vir: ',R_vir)
        print('R_vir_old: ',R_vir_old)
        print()
        
        # --- define some lists to populate later
        v_parallel_list = []
        v_parallel_inc_list = []
        vfinal_list = []
        v_projected_list = []
    
        ds_list = []
        ds_vel_list = []
        dsfinal_list = []
        dsvfinal_list = []

    ##########################################################################################
    ##########################################################################################
        # --- Define the galaxy plane

        # --- start with a normal pointing along the x-axis (plane is in z-y plane)
        N = np.array([1.,0.,0.])

        # --- rotate about y for inclination
        k_inc = np.array([0.,1.,0.])
        inc_rot = effectiveInc * np.pi/180.

        # --- rotate about x for PA
        k_pa = np.array([1.,0.,0.])
        pa_rot = (90. + pa) * np.pi/180.

        N_inc = rotate_vector(N, k_inc, inc_rot)
        print('N_inc: ',N_inc)
        print()

        N_inc_pa = rotate_vector(N_inc, k_pa, pa_rot)
        print('N_inc_pa: ',N_inc_pa)
        print()
    
        normal = N_inc_pa
        origin = np.array([0.,0.,0.])

    
        # Define ray -> [0,RA_dif,Dec_dif]
        rayDirection = np.array([1, 0, 0])

        rayPoint = np.array([0, impact_RA, impact_Dec])
        print('rayPoint: ',rayPoint)
        print()
    
    ##########################################################################################
    ##########################################################################################
        
        # --- NFW model  lists
        v_proj_list = []
        v_proj_min_list = []
        v_proj_max_list = []
        intersect_list = []
        d_plot_list = []
        intersect_point_list = []
        v_list = []
        v_90_list = []
    
        # --- Steidel model
        vlos_list = []
        vlos_max_list = []
        vlos_min_list = []
        dlos_list = []
    
        # --- how often to sample?
        if inc <= 60:
            s = 0.1
        elif inc <=80:
            s = 0.01
        elif inc <= 87:
            s = 0.005
        else:
            s = 0.001
        
        if diskOnly:
            zcutoff = 0.1
            s = 0.1
        
    ##########################################################################################

        #  --- first apply the model to just the thin-disk midplane
        z = 0
    #   midplane_model_projection = project_model_velocity(z, origin, normal, rayDirection, rayPoint, rcutoff, fit, spherical_halo=False)
        midplane_model_projection = project_model_velocity(z, origin, normal, rayDirection, rayPoint, rcutoff, fit, spherical_halo=spherical_halo, z_gradient=z_gradient)


        # --- unpack the results
        midplane_intersect           = midplane_model_projection["intersect"]
        midplane_intersect_dist      = midplane_model_projection["intersect_dist"]
        midplane_dist_from_origin    = midplane_model_projection["dist_from_origin"]
        midplane_v_intersect         = midplane_model_projection["v_intersect"]
        midplane_v_n_intersect       = midplane_model_projection["v_n_intersect"]
        midplane_v_rotation          = midplane_model_projection["v_rotation"]
        midplane_v_proj              = midplane_model_projection["v_proj"]

        midplane_dlos = midplane_intersect[0]
        
        vlos = steidel_model(vmax, effectiveInc, steidel_hv, impact, az, midplane_dlos)
        print('vlos: ',vlos)
        print('vmax: ',vmax)
    
        print('midplane_intersect: ',midplane_intersect)
        print('midplane_intersect_dist: ',midplane_intersect_dist)
        print('midplane_dist_from_origin: ',midplane_dist_from_origin)
        print('midplane_v_intersect: ',midplane_v_intersect)
        print('midplane_v_n_intersect: ',midplane_v_n_intersect)
        print('midplane_v_rotation: ',midplane_v_rotation)
        print('midplane_v_proj: ',midplane_v_proj)
        print('origin: ',origin)
        print('normal: ',normal)
        print('rayDirection: ',rayDirection)
        print('rayPoint: ',rayPoint)
        print()
        
    ##########################################################################################
        # --- now loop through layers of galaxy planes

        # --- the NFW model
        for i in arange(-zcutoff, zcutoff, s):
            model_projection    = project_model_velocity(i, origin, normal, rayDirection, rayPoint, rcutoff, fit, spherical_halo=spherical_halo, z_gradient=z_gradient)
            intersect           = model_projection["intersect"]
            intersect_dist      = model_projection["intersect_dist"]
            dist_from_origin    = model_projection["dist_from_origin"]
            v_intersect         = model_projection["v_intersect"]
            v_n_intersect       = model_projection["v_n_intersect"]
            v_rotation          = model_projection["v_rotation"]
            v_proj              = model_projection["v_proj"]
            
            # ---  deal with velocity errors
            model_projection_min = project_model_velocity(i, origin, normal, rayDirection, rayPoint, rcutoff, fit_min, spherical_halo=spherical_halo, z_gradient=z_gradient)
            model_projection_max = project_model_velocity(i, origin, normal, rayDirection, rayPoint, rcutoff, fit_max, spherical_halo=spherical_halo, z_gradient=z_gradient)
            v_proj_min           = model_projection_min["v_proj"]
            v_proj_max           = model_projection_max["v_proj"]
        
            # --- Where is the zero point? Center either at the midplane intersect or the
            # --- galaxy systemic velocity
            if center_at_intersect:
#                 intersect_xaxis = intersect[0]
                intersect_xaxis = intersect[0] - midplane_dlos
            else:
#                 intersect_xaxis = intersect[0] - midplane_dlos
                intersect_xaxis = intersect[0]

            if v_proj:
                v_proj_list.append(v_proj)
                v_proj_min_list.append(v_proj_min)
                v_proj_max_list.append(v_proj_max)
                intersect_list.append(intersect_xaxis)
                intersect_point_list.append(intersect)
                v_list.append(v_n_intersect)
                v_90_list.append(v_rotation)

    ##########################################################################################
    
        # --- Steidel model:
        if inc <= 30:
            steidel_step = 0.001
        elif inc <=50:
            steidel_step = 0.01
        else:
            steidel_step = 0.1
    
#         print 'min(intersect_list): ',min(intersect_list)
#         print 'max(intersect_list): ',max(intersect_list)
#         print 'intersect[0]: ',intersect[0]
#         print 'vmax: ',vmax
#         print
    
        # --- define the model range
        if center_at_intersect:
            steidel_min = -600
            steidel_max = 600
        else:
            steidel_min = -midplane_dlos - 600
            steidel_max = -midplane_dlos + 600
            
        
        # --- check if there has been a change in the sign of vmax (i.e., flipping from
        # --- approaching to receding due to PA change
        med_v_proj = np.median(v_proj_list)
        eff_vmax = vmax
        if med_v_proj > 0:
            if eff_vmax < 0:
                eff_vmax *= -1
                
        elif med_v_proj < 0:
            if eff_vmax > 0:
                eff_vmax *= -1
        else:
            print('med_v_proj = {}, eff_vmax = {}. Exiting...'.format(med_v_proj, eff_vmax))
            sys.exit()

        print('med_v_proj = {}, eff_vmax = {}'.format(med_v_proj, eff_vmax))
    
        for i in arange(steidel_min, steidel_max+steidel_step, steidel_step):
#             effective_r = np.sqrt(i**2 + intersect_dist**2)
#             effective_vc = fit(effective_r)
    #         print "effective_r = {0}, effective_vc = {1}".format(effective_r, effective_vc)
        
            # --- the model
            # -- defined: steidel_model(vc, inclination, hv, impact, azimuth, Dlos):
            if eff_vmax >0:
                eff_vmax_max = eff_vmax + e_vmax
                eff_vmax_min = eff_vmax - e_vmax

            else:
                eff_vmax_max = eff_vmax - e_vmax
                eff_vmax_min = eff_vmax + e_vmax
                
                
            vlos     = steidel_model(eff_vmax,     abs(effectiveInc), steidel_hv, impact, az, i)
            vlos_max = steidel_model(eff_vmax_max, abs(effectiveInc), steidel_hv, impact, az, i)
            vlos_min = steidel_model(eff_vmax_min, abs(effectiveInc), steidel_hv, impact, az, i)

            # --- center at midplane intersect or at galaxy systemic
            if center_at_intersect:
#                 intersect_xaxis = i + midplane_dlos
                intersect_xaxis = i

            else:
                intersect_xaxis = i + midplane_dlos
                
            if vlos > 200:
                print('vlos > 200: eff_vmax = {}, effectiveInc = {}, impact = {}, az = {}, i = {}'.format(eff_vmax, effectiveInc, impact, az, i))
                print()

            vlos_list.append(vlos)
            vlos_max_list.append(vlos_max)
            vlos_min_list.append(vlos_min)
            dlos_list.append(intersect_xaxis)
        
        print('eff_vmax = {}'.format(eff_vmax))
        print('median, max of vlos_list = {}, {}'.format(np.median(vlos_list), np.max(vlos_list)))
        print('np.median(vlos_max_list) = {}'.format(np.median(vlos_max_list)))
        print('np.median(vlos_min_list) = {}'.format(np.median(vlos_min_list)))
        print()
        # --- add to the results dictionary
        results_dict[pa] = {'steidel_v':vlos_list, 
                            'steidel_max_v':vlos_max_list,
                            'steidel_min_v':vlos_min_list,
                            'steidel_d':dlos_list, 
                            'v_proj_list':v_proj_list,
                            'v_proj_max_list':v_proj_max_list,
                            'v_proj_min_list':v_proj_min_list,
                            'intersect_list':intersect_list,
                            'intersect_point_list':intersect_point_list,
                            'midplane_dlos':midplane_dlos,
                            'normal':normal,
                            'rayPoint':rayPoint}


##########################################################################################
##########################################################################################
##########################################################################################
    # --- DO THE PLOTTING

##########################################################################################
##########################################################################################
    # --- define a bunch of colors
    color_blue = '#436bad' # french blue
    color_red = '#ec2d01'  # tomato red
    
    color_purple = '#7570b3'
    color_green = '#1b9e77'
    color_orange = '#d95f02'
    color_yellow = '#e6ab02'
    color_pink   = '#e7298a'
    color_blue_purple = '#7570b3'
    
    color_red2   = '#e41a1c'
    color_blue2 = '#377eb8'
    color_green2 = '#4daf4a'
    color_purple2 = '#984ea3'
    color_orange2 = '#ff7f00'
    
    color_green3  = '#1b9e77'
    color_orange3 = '#d95f02'
    color_purple3 = '#7570b3'
    color_pink3   = '#e7298a'
    color_lime3   = '#66a61e'
    color_gold3   = '#e6ab02'
    color_yellow3 = '#a6761d'
    color_grey3   = '#666666'
    
    color_light_green = '#e5f5f9'
    color_med_green = '#99d8c9'
    color_dark_green = '#2ca25f'
    
    color_light_orange = '#fff7bc'
    color_med_orange = '#fec44f'
    color_dark_orange = '#d95f0e'
    
    # --- alpha for error region shading
    alpha_error = 0.3

    # --- center PA
    center_results          = results_dict[pa_range[1]]
    center_steidel_vs       = center_results['steidel_v']            # previously vlos_list
    center_steidel_max_vs   = center_results['steidel_max_v']        # using vmax+e_vmax
    center_steidel_min_vs   = center_results['steidel_min_v']        # using vmax-e_vmax
    center_steidel_ds       = center_results['steidel_d']            # previously dlos_list
    center_NFW_vs           = center_results['v_proj_list']          # previously v_proj_list
    center_NFW_max_vs       = center_results['v_proj_max_list']      # using fit_max NFW profile
    center_NFW_min_vs       = center_results['v_proj_min_list']      # using fit_min NFW profile
    center_intersects       = center_results['intersect_list']       # previously intersect_list
    center_intersect_ps     = center_results['intersect_point_list'] # previously intersect_point_list
    center_midplane_ds      = center_results['midplane_dlos']        # previously midplane_dlos
    center_normals          = center_results['normal']               # previously normal
    center_rayPoints        = center_results['rayPoint']             # previously rayPoint

    # --- min PA
    low_results         = results_dict[pa_range[0]]
    low_steidel_vs      = low_results['steidel_v']            # previously vlos_list
    low_steidel_max_vs  = low_results['steidel_max_v']        # using vmax+e_vmax
    low_steidel_min_vs  = low_results['steidel_min_v']        # using vmax-e_vmax
    low_steidel_ds      = low_results['steidel_d']            # previously dlos_list
    low_NFW_vs          = low_results['v_proj_list']          # previously v_proj_list
    low_NFW_max_vs      = low_results['v_proj_max_list']      # using fit_max NFW profile
    low_NFW_min_vs      = low_results['v_proj_min_list']      # using fit_min NFW profile
    low_intersects      = low_results['intersect_list']       # previously intersect_list
    low_intersect_ps    = low_results['intersect_point_list'] # previously intersect_point_list
    low_midplane_ds     = low_results['midplane_dlos']        # previously midplane_dlos
    low_normals         = low_results['normal']               # previously normal
    low_rayPoints       = low_results['rayPoint']             # previously rayPoint

    #  --- max PA
    high_results        = results_dict[pa_range[2]]
    high_steidel_vs     = high_results['steidel_v']            # previously vlos_list
    high_steidel_max_vs = high_results['steidel_max_v']       # using vmax+e_vmax
    high_steidel_min_vs = high_results['steidel_min_v']       # using vmax-e_vmax
    high_steidel_ds     = high_results['steidel_d']            # previously dlos_list
    high_NFW_vs         = high_results['v_proj_list']          # previously v_proj_list
    high_NFW_max_vs     = high_results['v_proj_max_list']      # using fit_max NFW profile
    high_NFW_min_vs     = high_results['v_proj_min_list']      # using fit_min NFW profile
    high_intersects     = high_results['intersect_list']       # previously intersect_list
    high_intersect_ps   = high_results['intersect_point_list'] # previously intersect_point_list
    high_midplane_ds    = high_results['midplane_dlos']        # previously midplane_dlos
    high_normals        = high_results['normal']               # previously normal
    high_rayPoints      = high_results['rayPoint']             # previously rayPoint


    # --- Define the extent of the plot boundary
    plotExtent = round(1.5*rcutoff,-1)
    
    if plotExtent <= 400:
        plotExtent = 400
    elif plotExtent > 400 and plotExtent <=800:
        plotExtent = 800
    else:
        plotExtent = 1000
        
    print('Beginning Plotting')
    print('rcutoff: ',rcutoff)
    print()
    print()
    
    # --- how big to make the plotted cylinder?
    zHeight = zcutoff
    
    # --- plot velocity on the x-axis? Don't do this.
    plot_x_velocity = False
    
    # --- how many steps to take while plotting. Each step moves the sightline forward and 
    # --- populates the v_proj plot with that result
    steps = 1
    
    # tranpose the list of intersect points for plotting
    ip_xlist, ip_ylist, ip_zlist = np.array(center_intersect_ps).transpose()

    for i in arange(steps):
        i +=1
        print('i: ',i)

#         fig = plt.figure(figsize=(9,5))
#         fig = plt.figure(figsize=(9,5))
        fig = plt.figure(figsize=(7.7,5.7))

        # first plot the v_proj data
        ax = fig.add_subplot(1,1,1)
        
#         agnName = agnName.replace('_','\_')
#         fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)

        # --- number of actual intercept elements to take for each step
        len_step = len(center_intersects)/steps
        print('len_step: ',len_step)

        # --- x axis in kpc or km/s
        if plot_x_velocity:
            NFW_x_intersects   = (np.array(center_intersects)/1000.)*hubbleConstant
            low_x_intersects   = (np.array(low_intersects)/1000.)*hubbleConstant
            high_x_intersects  = (np.array(high_intersects)/1000.)*hubbleConstant
        else:
            NFW_x_intersects   = np.array(center_intersects)
            low_x_intersects   = np.array(low_intersects)
            high_x_intersects  = np.array(high_intersects)
            
        # --- the errors don't have exactly the same range, find the overlap to plot the
        # --- shading between; returns (the intersect  values, the indices of overlap from 
        # --- the first array, the indices of overlap from the second array)
#         x_intersects       = np.array(center_intersects)
#         low_x_intersects   = np.array(low_intersects)
#         high_x_intersects  = np.array(high_intersects)

        max_x = np.min([np.max(low_x_intersects), np.max(high_x_intersects)])
        min_x = np.max([np.min(low_x_intersects), np.min(high_x_intersects)])
        x_range = np.arange(min_x, max_x)
        
        steidel_max_x = np.min([np.max(low_steidel_ds), np.max(high_steidel_ds)])
        steidel_min_x = np.max([np.min(low_steidel_ds), np.min(high_steidel_ds)])
        steidel_x_range = np.arange(steidel_min_x, steidel_max_x)
        
        steidel_lower_err_vs = []
        steidel_higher_err_vs = []
        NFW_lower_err_vs = []
        NFW_higher_err_vs = []
        
        print('len of steidel min, center, max vs: ',len(low_steidel_min_vs), len(center_steidel_vs), len(high_steidel_min_vs))
        print('len of NFW min, center, max vs: ',len(low_NFW_min_vs), len(center_NFW_vs), len(high_NFW_min_vs))

        # --- check if the PA - e_PA (the lower PA) results in higher or lower absolute
        # --- velocities for Steidel model
        if np.max(np.absolute(low_steidel_vs)) < np.max(np.absolute(center_steidel_vs)):
            print('np.max(np.absolute(low_steidel_vs)) < np.max(np.absolute(center_steidel_vs)) :',np.max(np.absolute(low_steidel_vs)), np.max(np.absolute(center_steidel_vs)))

            # --- 'low' is below. Thus 'high' is above. The lowest value then
            # --- corresponds to the min vmax on the low PA set.
            steidel_lower_err_vs = low_steidel_min_vs
            steidel_higher_err_vs = high_steidel_max_vs
            
            steidel_lower_ds = low_steidel_ds
            steidel_higher_ds = high_steidel_ds
            
        else:
            # --- 'low' is above. Thus 'high' is below
            steidel_lower_err_vs = high_steidel_min_vs
            steidel_higher_err_vs = low_steidel_max_vs
            
            steidel_lower_ds = high_steidel_ds
            steidel_higher_ds = low_steidel_ds
        
        # --- check if the PA - e_PA (the lower PA) results in higher or lower absolute
        # --- velocities for NFW model
        if np.max(np.absolute(low_NFW_vs)) < np.max(np.absolute(center_NFW_vs)):
            print('np.max(np.absolute(low_NFW_vs)) < np.max(np.absolute(center_NFW_vs)) :',np.max(np.absolute(low_NFW_vs)), np.max(np.absolute(center_NFW_vs)))
            # --- 'low' is below. Thus 'high' is above. The lowest value then
            # --- corresponds to the min vmax on the low PA set.
            NFW_lower_err_vs = low_NFW_min_vs
            NFW_higher_err_vs = high_NFW_max_vs
            
            NFW_lower_x_intercepts = low_x_intersects
            NFW_higher_x_intercepts = high_x_intersects

        else:
            print('np.max(np.absolute(low_NFW_vs)) > np.max(np.absolute(center_NFW_vs)) :',np.max(np.absolute(low_NFW_vs)), np.max(np.absolute(center_NFW_vs)))
            # --- 'low' is above. Thus 'high' is below
            print('low is above')
            NFW_lower_err_vs = high_NFW_min_vs
            NFW_higher_err_vs = low_NFW_max_vs
        
            NFW_lower_x_intercepts = high_x_intersects
            NFW_higher_x_intercepts = low_x_intersects
            
        
        #  --- interpolate the min and max velocity curves
        low_steidel_error_interp  = interp1d(steidel_lower_ds,  steidel_lower_err_vs,  kind='cubic')
        high_steidel_error_interp = interp1d(steidel_higher_ds, steidel_higher_err_vs, kind='cubic')

        low_NFW_error_interp  = interp1d(NFW_lower_x_intercepts,  NFW_lower_err_vs,  kind='cubic')
        high_NFW_error_interp = interp1d(NFW_higher_x_intercepts, NFW_higher_err_vs, kind='cubic')

#         low_steidel_error_interp  = interp1d(low_steidel_ds,  low_steidel_vs,  kind='cubic')
#         high_steidel_error_interp = interp1d(high_steidel_ds, high_steidel_vs, kind='cubic')
# 
#         low_NFW_error_interp  = interp1d(low_x_intersects,  low_NFW_vs,  kind='cubic')
#         high_NFW_error_interp = interp1d(high_x_intersects, high_NFW_vs, kind='cubic')

        print('NFW_lower_x_intercepts: ', type(NFW_lower_x_intercepts))
        print('NFW_higher_x_intercepts: ',type(NFW_higher_x_intercepts))
        print()
        

        ##################################################################################
        # --- plot it
#         ax.scatter(x_intersects[:i*len_step], center_NFW_vs[:i*len_step], color='black', s=2)
        ax.plot(NFW_x_intersects[:i*len_step], center_NFW_vs[:i*len_step], color=color_blue, lw=1.5,
                label =  r'$\rm NFW ~ Model$')

        # --- plot the error regions
#         ax.plot(low_x_intersects[:i*len_step],  low_NFW_vs[:i*len_step],  color='black', ls='dashed', lw=0.4)
#         ax.plot(high_x_intersects[:i*len_step], high_NFW_vs[:i*len_step], color='black', ls='dashed', lw=0.4)

#         ax.plot(x_range, low_NFW_error_interp(x_range),  color='black', ls='dashed', lw=1.4)
#         ax.plot(x_range, high_NFW_error_interp(x_range), color='black', ls='dashed', lw=1.4)

        # --- shade the error region
        ax.fill_between(x_range, low_NFW_error_interp(x_range), high_NFW_error_interp(x_range), facecolor=color_blue, alpha=alpha_error)

        ylabel(r'$\rm Velocity ~[km~s^{-1}]$')
        
        ##################################################################################
        # --- plot the Steidel model results
#         ax.scatter(center_steidel_ds, center_steidel_vs, color=color_blue, s=2, label = 'Steidel Model')
        ax.plot(center_steidel_ds, center_steidel_vs, color=color_red, lw=1.5, ls='dotted',
                    label = r'$\rm Steidel~et~al.~(2002) ~ Model$')

        # --- plot the Steidel model error regions
#         ax.plot(low_steidel_ds,  low_steidel_vs,  color=color_blue, ls='dotted', label = 'Steidel Model', lw=0.4)
#         ax.plot(high_steidel_ds, high_steidel_vs, color=color_blue, ls='dotted', label = 'Steidel Model', lw=0.4)

#         ax.plot(steidel_x_range, low_steidel_error_interp(steidel_x_range),  color=color_blue, ls='dotted', label = 'Steidel Model', lw=1.4)
#         ax.plot(steidel_x_range, high_steidel_error_interp(steidel_x_range), color=color_blue, ls='dotted', label = 'Steidel Model', lw=1.4)
        
        # --- shade the error region
        ax.fill_between(steidel_x_range, low_steidel_error_interp(steidel_x_range), high_steidel_error_interp(steidel_x_range), facecolor=color_red, alpha=alpha_error)
        ##################################################################################

        # --- decide on the xlabel
        if plot_x_velocity:
            xlabel(r'$\rm Intersect ~[km~s^{-1}]$')
        else:
            xlabel(r'$\rm D_{los} ~[kpc]$')
    
        # --- Define the minimum and maximum extent for the x and y axes
        max_velocity = max(np.concatenate((center_NFW_vs, low_NFW_vs, high_NFW_vs, center_steidel_vs, low_steidel_vs, high_steidel_vs)))
        min_velocity = min(np.concatenate((center_NFW_vs, low_NFW_vs, high_NFW_vs, center_steidel_vs, low_steidel_vs, high_steidel_vs)))
        
        max_x = max(np.concatenate((NFW_x_intersects, center_steidel_ds)))
        min_x = min(np.concatenate((NFW_x_intersects, center_steidel_ds)))
        
        
        tick_num = 4.
        x_tick_num = 6.
        tick_spacing = round((max_velocity - min_velocity)/tick_num,-1)
        x_tick_spacing = round((max_x - min_x)/x_tick_num,-1)
        
        y_extent = max_velocity - min_velocity
        x_extent = max_x - min_x
        tick_spacing = round(y_extent + y_extent/4., -1)/tick_num
        x_tick_spacing = round(x_extent + x_extent/4., -1)/x_tick_num

        print('tick_spacing: ',tick_spacing)
        print('max(center_NFW_vs): ',max(center_NFW_vs))
        print('min(center_NFW_vs): ',min(center_NFW_vs))

#         xlim_pos = max_x + x_tick_spacing/2
#         xlim_neg = min_x - x_tick_spacing/2
        xlim_pos = steidel_max
        xlim_neg = steidel_min
        ylim_pos = max_velocity + tick_spacing/2
        ylim_neg = min_velocity - tick_spacing/2
        
        print('########################################')
        print('xlim_pos: ',xlim_pos)
        print('xlim_neg: ',xlim_neg)
        print('ylim_pos: ',ylim_pos)
        print('ylim_neg: ',ylim_neg)
        print('########################################')
        print
#         print 'x_intersects: ',x_intersects
#         print
#         print 'center_NFW_vs: ',center_NFW_vs
#         print
        
#         ax.set_xlim(xlim_neg, xlim_pos)
#         ax.set_ylim(ylim_neg,ylim_pos)

#         xlim_neg2 = math.floor(xlim_neg - tick_spacing)
#         xlim_pos2 = math.floor(xlim_pos + tick_spacing)
#         ylim_neg2 = math.floor(ylim_neg - tick_spacing)
#         ylim_pos2 = math.floor(ylim_pos + tick_spacing)

        
#         ax.set_xlim(math.floor(xlim_neg), math.floor(xlim_pos))
#         ax.set_ylim(math.floor(ylim_neg), math.floor(ylim_pos))

#         ax.set_xlim(int(xlim_neg), int(xlim_pos))
#         ax.set_ylim(int(ylim_neg), int(ylim_pos))

#         ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        x_tick_spacing = (int(xlim_pos) - int(xlim_neg))/x_tick_num
        tick_spacing =  (int(ylim_pos) - int(ylim_neg))/tick_num
        
        x_tick_spacing = 200.0
        
#         print 'tick_spacing = ',tick_spacing
#         print 'after tick_spacing: ',tick_spacing

        if tick_spacing < 25:
            tick_spacing = 25.0
        
        elif tick_spacing == 0.0:
#             tick_spacing = round((max_velocity - min_velocity)/tick_num,5)
            tick_spacing = my_round((max_velocity - min_velocity)/tick_num, base=0.05)
        else:
            tick_spacing = my_round(tick_spacing, base=50)


        # x-axis
        majorLocator   = MultipleLocator(x_tick_spacing)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(x_tick_spacing/4)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(tick_spacing)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(tick_spacing/5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
#         ax.set_xbound(lower=int(xlim_neg), upper=int(xlim_pos))
#         xlim(int(xlim_neg), int(xlim_pos))
        ax.set_xbound(lower=math.floor(xlim_neg), upper=math.floor(xlim_pos))
        
        y_lower = int(ylim_neg)
        y_upper = int(ylim_pos)
        
        # --- try to make sure the y limits are at least +/- 25 in the right direction
        # --- +1 in either direction to make the ticks look better
        if y_lower < 0:
            if abs(y_lower) < 25:
                y_lower = -26.0
        
        if y_upper > 0:
            if y_upper < 25:
                y_upper = 26.0
                

        ylim(y_lower, y_upper)
        xlim(-600.0, 600.0)

        if tick_spacing <1:
            ylim(ylim_neg-ylim_neg/10, ylim_pos + ylim_pos/10)
            
#         ax.legend(loc='lower right', borderpad=0.4, fontsize=12, fancybox=True)
        ax.legend(borderpad=0.4, fontsize=12, fancybox=True)
#         ax.grid(True, alpha=0.5)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')

#         savefig("{0}{1}a.jpg".format(saveDirectory, i), dpi=300, bbox_inches='tight')
        savefig("{0}{1}a.pdf".format(saveDirectory, i), bbox_inches='tight')
        close(fig)
    
##########################################################################################
##########################################################################################
        # next plot the 3d fig
        fig = plt.figure(figsize=(8,10))
        
        # first plot the v_proj data
        ax = fig.add_subplot(1,1,1,projection='3d')

        # plot the sightline
        z = np.ones(1000)*rayPoint[2]
        x = np.linspace(-plotExtent, plotExtent, 1000)
        y = np.ones(1000)*rayPoint[1]
        
        z_sightline = z
        y_sightline = y
        x_sightline = x
                
        len_step = len(z)/steps
        
        ax.plot(x[:i*len_step], y[:i*len_step], z[:i*len_step], color='black', alpha = 0.6, lw=plotExtent/190)

        # some interesting points: 
        # v is the velocity vector in the direction of intersect point
#         v = np.array(v_list[0])
    
        # v_90 is the rotation velocity vector in the direction of rotation at the intersect
        # point
#         v_90 = np.array(v_90_list[0])
    
        # galaxy center
#         orig = np.array([0,0,0])
    
        # intersect point
        intersect = np.array(intersect_point_list[0])
        print('intersect: ',intersect[0],intersect[1],intersect[2])

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
        # plot the cylinder
    
        R = int(rcutoff)
        p0 = normal * (zHeight)
        p1 = normal * (-zHeight)
    
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
        ax.plot_surface(X,  Y,  Z,  color=colorTube,   alpha = alphaTube)
        ax.plot_surface(X2, Y2, Z2, color=colorBottom, alpha = alphaBottom)
        ax.plot_surface(X3, Y3, Z3, color=colorTop,    alpha = alphaTop)
        
        ax.set_xlim(-plotExtent, plotExtent)
        ax.set_ylim(-plotExtent, plotExtent)
        ax.set_zlim(-plotExtent, plotExtent)
    
        # reverse the RA axis so negative is on the right
    #     ax = plt.gca()
        ax.invert_xaxis()

        # rotate the plot
#         ax.view_init(elev=-65., azim=-39)
        ax.view_init(elev=5, azim=45)

        print('R: ',R)
        print('p0: ',p0)
        print('p1: ',p1)
#         print('tube: ',tube)
#         print('bottom: ',bottom)
#         print("top: ",top)
        print()

        # reverse the X-axis ticks without actually reversing the x-axis
        if plotExtent <= 400:
            yticks((-400, -200, 0, 200, 400), (400, 200, 0, -200, -400))
            xticks((-400, -200, 0, 200, 400), (-400, -200, 0, 200, 400))
        elif plotExtent >400 and plotExtent <=800:
            yticks((-800, -400, 0, 400, 800), (800, 400, 0, -400, -800))
            xticks((-800, -400, 0, 400, 800), (-800, -400, 0, 400, 800))
        else:
            yticks((-1000, -500, 0, 500, 1000), (1000, 500, 0, -500, -1000))
            xticks((-1000, -500, 0, 500, 1000), (-1000, -500, 0, 500, 1000))

        [t.set_va('center') for t in ax.get_yticklabels()]
        [t.set_ha('center') for t in ax.get_yticklabels()]
        
        [t.set_va('bottom') for t in ax.get_xticklabels()]
        [t.set_ha('left') for t in ax.get_xticklabels()]
        
        [t.set_va('center') for t in ax.get_zticklabels()]
        [t.set_ha('right') for t in ax.get_zticklabels()]
        
#         ax.grid(True)
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
        
        ax.set_zlabel(z_label, rotation=0, fontsize=14, labelpad=40)
#         ax.set_xlabel(x_label, rotation=0, fontsize=14, labelpad=40)
#         ax.set_ylabel(y_label, rotation=0, fontsize=14, labelpad=40)

#         ax.xaxis.set_label_coords(10.0, -200.02)
        
        ax.xaxis.labelpad=18.0
        ax.yaxis.labelpad=1.0
        ax.zaxis.labelpad=28.0

        tight_layout()


        # disable auto rotation
#         ax.xaxis.set_rotate_label(False)
#         ax.yaxis.set_rotate_label(False)

#         yticks((-1000, -500, 0, 500, 1000), (1000, 500, 0, -500, -1000))

#         xticks((-400, -200, 0, 200, 400), (-400, -200, 0, 200, 400))
#         ax.set_zticks((-400, -200, 0, 200, 400), (-400, -200, 0, 200, 400))

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#         if i == 2:
#             show()
    
#         savefig("{0}{1}b.jpg".format(saveDirectory, i), dpi=300, bbox_inches='tight')
        savefig("{0}{1}b.pdf".format(saveDirectory, i), bbox_inches='tight', pad_inches=0.2)
        close(fig)
#         show()
    
##########################################################################################
##########################################################################################
    # write out a summary txt file
    
    summary_file = open(saveDirectory + 'summary.txt','wt')
    
    summary_file.write('Galaxy name: {0}\n'.format(galaxyName))
    summary_file.write('Target name: {0}\n'.format(agnName))
    summary_file.write('Impact: {0}\n'.format(impact))
    summary_file.write('Azimuth: {0}\n'.format(az))
    summary_file.write('R_vir: {0}\n'.format(R_vir))
    summary_file.write('Inclination: {0} pm {1}\n'.format(inc, e_inc))
    summary_file.write('PA: {0} +/- {1}\n'.format(PA, e_PA))
    summary_file.write('\n')
    summary_file.write('Max V_rot / sin(i): {0} \pm {1}\n'.format(vmax_incCorrected, e_vmax_incCorrected))
    summary_file.write('\n')
    summary_file.write('flipInclination: {0}\n'.format(flipInclination))
    summary_file.write('reverse: {0}\n'.format(reverse))
    summary_file.write('fit_NFW: {0}\n'.format(fit_NFW))
    summary_file.write('Fit popt: [V200, c, R200] = {0}\n'.format(popt,))
    summary_file.write('Fit e_popt: [e_V200, e_c, e_R200] = [{0}, {1}, {2}]\n'.format(v200_err, c_err, r200_err))
    summary_file.write('Fit popt_min: [V200, c, R200] = {0}\n'.format(popt_min))
    summary_file.write('Fit popt_max: [V200, c, R200] = {0}\n'.format(popt_max))
    summary_file.write('\n')
    summary_file.write('\n')
    summary_file.write('------ NFW Model -------\n')
    summary_file.write('\n')
    summary_file.write('NFW V_proj range: [{0}, {1}]\n'.format(min(center_NFW_vs), max(center_NFW_vs)))
    summary_file.write('\n')
    summary_file.write('NFW D_los range: [{0}, {1}]\n'.format(min(NFW_x_intersects), max(NFW_x_intersects)))
    summary_file.write('\n')
    summary_file.write('\n')
    summary_file.write('NFW V_proj range lower error: [{0}, {1}]\n'.format(min(NFW_lower_err_vs), max(NFW_lower_err_vs)))
    summary_file.write('\n')
    summary_file.write('NFW D_los range lower error: [{0}, {1}]\n'.format(min(NFW_lower_x_intercepts), max(NFW_lower_x_intercepts)))
    summary_file.write('\n')
    summary_file.write('NFW V_proj range upper error: [{0}, {1}]\n'.format(min(NFW_higher_err_vs), max(NFW_higher_err_vs)))
    summary_file.write('\n')
    summary_file.write('NFW D_los range upper error: [{0}, {1}]\n'.format(min(NFW_higher_x_intercepts), max(NFW_higher_x_intercepts)))
    summary_file.write('\n')
    summary_file.write('-------- Steidel Model ---------\n')
    summary_file.write('\n')
    summary_file.write('Steidel Model V_proj range: [{0}, {1}]\n'.format(min(center_steidel_vs), max(center_steidel_vs)))
    summary_file.write('\n')
    summary_file.write('Steidel Model D_los range: [{0}, {1}]\n'.format(min(center_steidel_ds), max(center_steidel_ds)))
    summary_file.write('\n')
    summary_file.write('\n')
    summary_file.write('Steidel Model V_proj range lower error: [{0}, {1}]\n'.format(min(steidel_lower_err_vs), max(steidel_lower_err_vs)))
    summary_file.write('\n')
    summary_file.write('Steidel Model D_los range lower error: [{0}, {1}]\n'.format(min(steidel_lower_ds), max(steidel_lower_ds)))
    summary_file.write('\n')
    summary_file.write('Steidel Model V_proj range upper error: [{0}, {1}]\n'.format(min(steidel_higher_err_vs), max(steidel_higher_err_vs)))
    summary_file.write('\n')
    summary_file.write('Steidel Model D_los range upper error: [{0}, {1}]\n'.format(min(steidel_higher_ds), max(steidel_higher_ds)))
    summary_file.write('\n')

    # this part sums the x and y velocities. For example, if you're at -10 km/s along the 
    # sightline and the measured velocity is -10, the velocity at that point is -20 km/s.
    if plot_x_velocity:
        totals = []
        for x, y in zip(NFW_x_intersects, center_NFW_vs):
            totals.append(x+y)
    
        totals.sort()
        total_nozeros = []
        for i in totals:
            if i != 0:
                total_nozeros.append(i)
    
        total_min = min(total_nozeros)
    
        summary_file.write('\n')
        summary_file.write('Combined velocity range: [{0}, {1}]'.format(total_min, max(totals)))
        summary_file.write('\n')

    summary_file.close()
    
    ####################
    # --- Save results to a csv file
    
    if save_results_csv:
        # --- check if the csv file exists already
        csv_file_exists = os.path.exists(csv_filename)
        print('csv_file_exists: ',csv_file_exists)
        
        # --- open the csv file in append mode ('a')
        csv_file = open(csv_filename, 'a')

        # --- find the bounds of each error array (i.e., the min and max values of error regions)
        e_vs = [min(NFW_lower_err_vs), max(NFW_lower_err_vs), min(NFW_higher_err_vs), max(NFW_higher_err_vs)]
        e_ds = [min(NFW_lower_x_intercepts), max(NFW_lower_x_intercepts), min(NFW_higher_x_intercepts), max(NFW_higher_x_intercepts)]
        
        NFW_final_vs = [min(min(e_vs), min(center_NFW_vs)), max(max(e_vs), max(center_NFW_vs))]
        NFW_final_ds = [min(min(e_ds), min(NFW_x_intersects)), max(max(e_ds), max(NFW_x_intersects))]
        
        e_steidel_vs = [min(steidel_lower_err_vs), max(steidel_lower_err_vs), min(steidel_higher_err_vs), max(steidel_higher_err_vs)]
        e_steidel_ds = [min(steidel_lower_ds), max(steidel_lower_ds) , min(steidel_higher_ds), max(steidel_higher_ds)]
        
        steidel_final_vs = [min(min(e_steidel_vs), min(center_steidel_vs)), max(max(e_steidel_vs), max(center_steidel_vs))]
        steidel_final_ds = [min(min(e_steidel_ds), min(center_steidel_ds)), max(max(e_steidel_ds), max(center_steidel_ds))]

        
        fieldnames = ['Galaxy',
                      'Targetname',
                      'impact',
                      'azimuth',
                      'side',
                      'MajDiam',
                      'Rvir',
                      'Rvir_stocke',
                      'Lstar',
                      'e_Lstar',
                      'inc',
                      'e_inc',
                      'PA',
                      'e_PA',
                      'Vhel_measured',
                      'e_Vhel_measured',
                      'Vhel_published',
                      'e_Vhel_published',
                      'right_vrot_incCorrected_avg',
                      'left_vrot_incCorrected_avg',
                      'vmax',
                      'e_vmax',
                      'vmax_incCorrected',
                      'e_vmax_incCorrected',
                      'V200',
                      'c',
                      'R200',
                      'V200_min',
                      'c_min',
                      'R200_min',
                      'V200_max',
                      'c_max',
                      'R200_max',
                      'NFW_final_vs',
                      'NFW_final_ds',
                      'steidel_final_vs',
                      'steidel_final_ds',
                      'min_NFW_vs',
                      'max_NFW_vs',
                      'min_NFW_ds',
                      'max_NFW_ds',
                      'e_min_NFW_vs',
                      'e_max_NFW_vs',
                      'e_min_NFW_ds',
                      'e_max_NFW_ds',
                      'min_steidel_vs',
                      'max_steidel_vs',
                      'min_steidel_ds',
                      'max_steidel_ds',
                      'e_min_steidel_vs',
                      'e_max_steidel_vs',
                      'e_min_steidel_ds',
                      'e_max_steidel_ds',
                      'zcutoffm',
                      'rcutoffm',
                      'spherical_halo',
                      'z_gradient',
                      'steidel_hv',
                      'date']
                      
        write_list = [galaxyName,
                    agnName,
                    np.around(impact,2),
                    az,
                    side,
                    majDiam,
                    np.around(R_vir_old,2),
                    np.around(R_vir,2),
                    Lstar,
                    e_Lstar,
                    effectiveInc,
                    e_inc,
                    PA,
                    e_PA,
                    vsys_measured,
                    vsys_measured_err,
                    vsys_published,
                    dvsys_published,
                    right_vrot_incCorrected_avg,
                    left_vrot_incCorrected_avg,
                    np.around(vmax,2),
                    np.around(e_vmax, 2),
                    np.around(vmax_incCorrected,2),
                    np.around(e_vmax_incCorrected, 2),
                    np.around(popt[0],2),
                    np.around(popt[1],2),
                    np.around(popt[2],2),
                    np.around(popt_min[0],2),
                    np.around(popt_min[1],2),
                    np.around(popt_min[2],2),
                    np.around(popt_max[0],2),
                    np.around(popt_max[1],2),
                    np.around(popt_max[2],2),
                    np.around(NFW_final_vs,2),
                    np.around(NFW_final_ds,2),
                    np.around(steidel_final_vs,2),
                    np.around(steidel_final_ds,2),
                    np.around(min(center_NFW_vs),2),
                    np.around(max(center_NFW_vs),2),
                    np.around(min(NFW_x_intersects),2),
                    np.around(max(NFW_x_intersects),2),
                    np.around(min(e_vs),2),
                    np.around(max(e_vs),2),
                    np.around(min(e_ds),2),
                    np.around(max(e_ds),2),
                    np.around(min(center_steidel_vs),2),
                    np.around(max(center_steidel_vs),2),
                    np.around(min(center_steidel_ds),2),
                    np.around(max(center_steidel_ds),2),
                    np.around(min(e_steidel_vs),2),
                    np.around(max(e_steidel_vs),2),
                    np.around(min(e_steidel_ds),2),
                    np.around(max(e_steidel_ds),2),
                    zcutoffm,
                    rcutoffm,
                    spherical_halo,
                    z_gradient,
                    steidel_hv,
                    [str(datetime.datetime.now())]]
        
        
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        
        # --- write headers if the file is new
        if not csv_file_exists:
            headers = dict((n,n) for n in fieldnames)
            writer.writerow(headers)
       
        row = dict((f,o) for f,o in zip(fieldnames,write_list))
        writer.writerow(row)
        csv_file.close()
        
    #######################
    # --- save pickle data?
    if save_data_pickle:
        pickle_filename = '{0}/{1}-{2}_model_NFW{3}.p'.format(saveDirectory, galaxyName, agnName, fit_NFW)
    
        pickle_file = open(pickle_filename,'wt')
        
        d = {}
        
        # --- distances/velocity along the sightline
        d['NFW_x_intersects'] = NFW_x_intersects
        d['NFW_lower_x_intercepts'] = NFW_lower_x_intercepts
        d['NFW_higher_x_intercepts'] = NFW_higher_x_intercepts
        
        d['center_steidel_ds'] = center_steidel_ds
        d['steidel_lower_ds'] = steidel_lower_ds
        d['steidel_higher_ds'] = steidel_higher_ds
        
        # --- projected rotation velocity
        d['center_NFW_vs'] = center_NFW_vs
        d['NFW_higher_err_vs'] = NFW_higher_err_vs
        d['NFW_lower_err_vs'] = NFW_lower_err_vs
        
        d['steidel_center_vs'] = center_steidel_vs
        d['steidel_lower_err_vs'] = steidel_lower_err_vs
        d['steidel_higher_err_vs'] = steidel_higher_err_vs

        # --- measured max velocities with errors
        d['vmax'] = vmax
        d['e_vmax'] = e_vmax
        d['vmax_incCorrected'] = vmax_incCorrected
        d['e_vmax_incCorrected'] = e_vmax_incCorrected
        
        # --- now the sightline
        d['z_sightline'] = z_sightline
        d['y_sightline'] = y_sightline
        d['x_sightline'] = x_sightline
        
        # --- now the cylinder
        d['intersect'] = intersect
        d['R'] = R
        d['p0'] = p0
        d['p1'] = p1
        
        d['popt'] = popt

        pickle.dump(d, pickle_file)
        pickle_file.close()

    print('Finished {0}'.format(galaxyName))
    print()
    
if __name__ == '__main__':
#     galaxy_names = ['CGCG039-137',
#                     'ESO343-G014',
#                     'IC5325',
#                     'MCG-03-58-009',
#                     'NGC1566',
#                     'NGC3513',
#                     'NGC3633',
#                     'NGC4536',
#                     'NGC4939',
#                     'NGC5364',
#                     'NGC5786',
#                     'UGC09760',
#                     'NGC3198',
#                     'NGC4565',
#                     'NGC3351',
#                     'UGC04238',
#                     'NGC4529',
#                     'NGC6140',
#                     'NGC5907',
#                     'UGC06446',
#                     'NGC3631',
#                     'NGC3726',
#                     'NGC3067',
#                     'NGC2770',
#                     'NGC3432',
#                     'NGC3666',
#                     'NGC5951',
#                     'NGC7817',
#                     'UGC08146']

    # --- have at least 2 sightlines
#     galaxy_names = ['CGCG039-137',
#                     'NGC1566',
#                     'NGC3633',
#                     'NGC4536',
#                     'NGC2770',
#                     'NGC3067',
#                     'NGC3198',
#                     'NGC3432',
#                     'NGC3631',
#                     'NGC3666',
#                     'NGC3726',
#                     'NGC5907']

    # --- have at least 3 sightlines
#     galaxy_names = ['NGC1566',
#                     'NGC3633',
#                     'NGC2770',
#                     'NGC3067',
#                     'NGC3432',
#                     'NGC3631',
#                     'NGC3666',
#                     'NGC5907']

    # --- have at least 4 sightline
#     galaxy_names = ['NGC1566',
#                     'NGC2770',
#                     'NGC3631']

    # --- have 5 sightlines
#     galaxy_names = ['NGC1566',
#                     'NGC2770']

    galaxy_names = ['NGC3633']
    
    
    for name in galaxy_names:
        main(name)
    