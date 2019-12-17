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
    
    print 'log_L[0]: ',log_L[0]
    print 'log_rvir[0]: ',log_rvir[0]
    print
    
    # --- get out of log space
    L = 10**log_L
    rvir = 10**log_rvir
    
    # --- interpolate it using a cubic spline
    interp = interp1d(L, rvir, kind='cubic')
    
    return interp



def project_model_velocity(z, origin, normal, rayDirection, rayPoint, rcutoff, fit_model, sphericalHalo=True, verbose=False):
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
        
        sphericalhalo : boolean
                      True means the total distance from the origin is given to fit_model.
                      This results in a halo which is spherical instead of cylindrical
                      
                      False means the distance from the current plane origin is given.
        
        
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
    

    # this is a point in the new, parallel but shifted plane
    planePoint = origin + (z * normal)

    # get intersect: find_intersect(planeNormal,planePoint,rayDirection,rayPoint)
    intersect = find_intersect(normal, planePoint, rayDirection, rayPoint)
    
    # this is the vector from the origin of the current plane to the intersect
    intersect_vect = intersect - (z * normal)
    
    # this is the distance from the origin of the current plane to the intersect
    intersect_dist = np.linalg.norm(intersect_vect)
    
    # this is the distance from the origin (galaxy center) to the current intersect pt
    dist_from_origin = np.linalg.norm(intersect)

    # restrict the intersection to be within the cylinder of radius rcutoff
    if intersect_dist <= rcutoff:
        if verbose:
            print 'planePoint: ',planePoint
            print "intersection at", intersect
            print 'intersect_dist: ',intersect_dist

        # find the rotation velocity at this distance from the rotation curve fit center
        try:
            if sphericalHalo:
                v_intersect = fit_model(dist_from_origin)
            else:
                v_intersect = fit_model(intersect_dist)
                if verbose:
                    print 'v_intersect: ',v_intersect
        except Exception,e:
            # if you go beyond the fit, set velocity to 0
            v_intersect = 0
            if verbose:
                print 'Ran out of interpolation range for {0}'.format(intersect_dist)
                print "Built in exception is {0}".format(e)
            sys.exit()
            
        #######
        #######
        #######
        # angle between sightline and vector to intersect point
    
        # unit vector towards intersect point
        n_intersect_dist = intersect_vect / np.linalg.norm(intersect_vect)
        if verbose:
            print 'n_intersect_dist: ',n_intersect_dist
            print 'np.linalg.norm(n_intersect_dist): ',np.linalg.norm(n_intersect_dist)

    
        # this is the velocity vector in the direction of intersect point, n_intersect_dist
        # edit: seems legit
        v_n_intersect = v_intersect * n_intersect_dist
        if verbose:
            print 'new way: '
            print 'v_n_intersect: ',v_n_intersect
            print '||v_n_intersect|| : ',np.linalg.norm(v_n_intersect)
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
            print 'v_rotation: ',v_rotation
            print '||v_rotation|| : ',np.linalg.norm(v_rotation)
            print '||N||: ', np.linalg.norm(normal)
#         v_90_list.append(v_rotation)
    
        # now dot it with the sightline to get the component along
        v_proj = np.dot(v_rotation, rayDirection)
        if verbose:
            print 'v_proj: ',v_proj
#         v_proj_list.append(v_proj)
        
#         intersect_list.append(intersect[0])
            print 'intersect: ',intersect
            print 'intersect[0]: ',intersect[0]
            print
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
        vmax     : float
                maximum rotation velocity (inclination corrected)
                
        hv       : float
                scale height of the disk
                
        impact   : float
                impact parameter
                
        azimuth  : float
                azimuth angle
                
        y        : float
                the y component
        
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
    
    vlos = -(vc  * np.sin(inc_radians) / sqrt_term) * exp_term
#     vlos = -(vc / sqrt_term) * exp_term

    return vlos
    
        
        
        
def plot_NFW(xData, yData, popt, x_lim):
    """
    This function makes a nice looking plot showing the NFW fit and data
    
    Parameters
    ----------
    xData   : np.array
            the x values for the observed rotation curve 
            
    yData   : np.array
            the y values for the observed rotation curve
            
    popt    : np.array
            the fit values
            
    x_lim   : float
            the maximum x value to plot to
    
    Returns
    -------
    fig     : the figure object, to be shown or saved
    
    """
    
#     fig = plt.figure(figsize=(8,8))
    fig = plt.figure(figsize=(7.7,5.7))
    ax = fig.add_subplot(1,1,1)
    
    v200, c, r200 = popt
    
    x_fit = linspace(0,x_lim,num=1000)
    
    scatter(xData, yData, color='black', s=40, lw=0, label = r'$\rm Data$')
#     plot(x_fit, NFW(x_fit, *popt), 'r-',label='fit: a={0}, rho={1}'.format(*popt))
    plot(x_fit, NFW(x_fit, *popt), 'r-', color='green', alpha=0.7, lw=2, \
    label='NFW Fit: V200={0}, c={1}, R200={2}'.format(round(v200,2),round(c,2),round(r200,2)))
    
    legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)
    
    xlim(0, x_lim)
    ylim(0, round(np.nanmax(yData),-1) + 15)

    # x-axis
    majorLocator   = MultipleLocator(25)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_minor_locator(minorLocator)

    # y axis
    majorLocator   = MultipleLocator(50)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(10)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_minor_locator(minorLocator)

    xlabel(r'$\rm R ~[kpc]$')
    ylabel(r'$\rm \emph{v}_{{rot}} ~[km s^{{-1}}]$')
    
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




def main():
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
    zcutoffm = 2.0
    
    # --- radius multiplier (* R_vir)
    rcutoffm = 3.0
    
    # --- hv for the Steidel model
    steidel_hv = 1000.0
    
    # save the full model velocity data in a pickle file?
    save_data_pickle = True
    
    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red
    
    # --- interpolate the Stocke et al. (2013) Lstar vs Rvir relation
    stocke_rvir_filename = '/Users/frenchd/Research/rotation_paper_data/Rvir_L_Lstar2.csv'
    stocke_rvir = stocke_rvir_interpolate(stocke_rvir_filename)
    
    # ---  Use the Stocke + 2013 Rvir calculation or the usual one?
    use_stocke_rvir = True
    
    # --- center at the midplane intersect? centers at galaxy systemic if false
    center_at_intersect = False


    galaxyName = 'CGCG039-137'
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

#     saveDirectory = '/Users/frenchd/Research/M31_rotation/{0}/model_v3/'.format(galaxyName)
#     saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/models_v2/{0}/'.format(galaxyName)
    saveDirectory = '/Users/frenchd/Research/rotation_paper_data/rotation_models_v2/{0}/'.format(galaxyName)


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
    directory = '/Users/frenchd/Research/rotation_paper_data/summary_files/'


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
        
        xVals = data['xVals']
        inc = data['inclination']
        vsys_measured = data['vsys_measured']
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
            e_R_vir_up = stocke_rvir(Lstar + e_Lstar) - R_vir
            e_R_vir_down = R_vir - stocke_rvir(Lstar - e_Lstar)
        else:
            R_vir = calculateVirialRadius(majDiam)
            e_R_vir_up = calculateVirialRadius(majDiam + majDiam*0.1)
            e_R_vir_down = calculateVirialRadius(majDiam - majDiam*0.1)

        # --- define "vmax", the average maximum rotation velocity
        vmax_incCorrected = max(abs(right_vrot_incCorrected_avg), abs(left_vrot_incCorrected_avg))
        
        # -- remove inclination correction for the Steidel model method
        vmax = vmax_incCorrected * np.sin(inc * np.pi/180.)
        
        print
        print 'PA: ',PA
        print 'inc: ',inc
        print 'dist: ',dist
        print 'AGN: ',agn
        print
        print 'vmax_incCorrected: ',vmax_incCorrected
        print 'vmax: ',vmax
        print
        print 'Lstar: ',Lstar
        print 'stocke_rvir(Lstar) = ',stocke_rvir(Lstar)
        print
        

        # which agn do you want to target? Also decide here if you want to mirror around
        # the inclination axis (i.e., so the galaxy is "facing" left vs right for PA=0)
        # 
        # reverse = True to reverse the rotation direction
        #
        # NFW_fit decides how tightly to bound the NFW profile fit. 
        # Options are: standard, tight, tighter, tightest, tightester
        
        NFW_fit = "standard"

        # CGCG039-137
        if galaxyName == 'CGCG039-137':
            flipInclination = False
            reverse = False
            agnName = 'RX_J1121.2+0326'
#             agnName = 'SDSSJ112224.10+031802.0'

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
        
        # M31
        if galaxyName == 'M31':
            flipInclination = False
            # --- reverse for NFW
            if fit_NFW:
                reverse = True
            else:
                reverse = False
            agnName = 'QSO1'
        
        
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
            agnName = 'HE0439-5254'
#             agnName = 'HE0435-5304'
#             agnName = 'RBS567'
#             agnName = 'HE0429-5343'
#             agnName = '1H0419-577'


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
            agnName = 'SDSSJ112005.00+041323.0'
#             agnName = 'RX_J1121.2+0326'
#             agnName = 'SDSSJ112224.10+031802.0'


        # NGC4536
        if galaxyName == 'NGC4536':
            # --- use NGC4536-summary7.json for NFW fits (this file has the weird data on the 
            # --- (left of NGC4536 removed, resulting in a much better NFW fit)
            # --- NGC4536-summary6.json has all the original data intact
            flipInclination = False
            reverse = False
            agnName = '3C273.0'
#             agnName = 'HE1228+0131'


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

        # NGC3198
        if galaxyName == 'NGC3198':
            flipInclination = True
            reverse = True
            agnName = 'RX_J1017.5+4702'
#             agnName = 'SDSSJ101622.60+470643.0'


        # NGC4565
        if galaxyName == 'NGC4565':
            flipInclination = False
            reverse = False
            agnName = 'RX_J1236.0+2641'


        # NGC3351
        if galaxyName == 'NGC3351':
            flipInclination = True
            reverse = True
            agnName = 'SDSSJ104335.90+115129.0'
#             agnName = 'SDSSJ104341.53+085558.2'
#             agnName = 'SDSSJ104709.80+130454.0'
#             agnName = 'SDSSJ104816.30+120735.0'
#             agnName = 'SDSSJ104843.50+130605.0'
#             agnName = 'SDSSJ105220.60+101751.0'


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


        # NGC4529
        if galaxyName == 'NGC4529':
            NFW_fit = 'tightester'
            flipInclination = False
            reverse = True
            agnName = 'MRK771'


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


        # NGC5907
        if galaxyName == 'NGC5907':
            flipInclination = False
            reverse = False
            agnName = 'SBS1503+570'
#             agnName = 'SDSSJ152053.59+571122.1'
#             agnName = 'RBS1503'


        # UGC06446
        if galaxyName == 'UGC06446':
            flipInclination = False
            reverse = False
            agnName = 'SDSSJ112448.30+531818.0'
#             agnName = 'RX_J1117.6+5301'


        # NGC3631
        if galaxyName == 'NGC3631':
            flipInclination = False
            reverse = False
            agnName = 'SDSSJ111443.70+525834.0'
#             agnName = 'RX_J1117.6+5301'
#             agnName = 'SBS1116+523'
#             agnName = 'SDSSJ112448.30+531818.0'


        # UGC06399
        if galaxyName == 'UGC06399':
            flipInclination = False
            reverse = True
            agnName = 'SBS1116+523'


        # NGC3726
        if galaxyName == 'NGC3726':
            flipInclination = True
            reverse = True
            NFW_fit = 'standard'
            agnName = 'CSO1208'
#             agnName = 'RX_J1142.7+4625'


        # NGC3067
        if galaxyName == 'NGC3067':
            flipInclination = False
            reverse = True
            agnName = '3C232'
#             agnName = 'RX_J1002.9+3240'
#             agnName = 'SDSSJ095914.80+320357.0'


        # NGC2770
        if galaxyName == 'NGC2770':
            flipInclination = True
            NFW_fit = 'tighter'
            # --- reverse for NFW, not for 2x3R_vir
            if fit_NFW:
                reverse = True
            else:
                reverse = False

            agnName = 'FBQSJ0908+3246'
#             agnName = 'TON1009'
#             agnName = 'TON1015'
#             agnName = 'SDSSJ091052.80+333008.0'
#             agnName = 'SDSSJ091127.30+325337.0'


        # NGC3432
        if galaxyName == 'NGC3432':
            flipInclination = True
            NFW_fit = 'tighter'
            # --- reverse for NFW, not for 2x3R_vir
            if fit_NFW:
                reverse = True
            else:
                reverse = False

            agnName = 'MS1047.3+3518'
#             agnName = 'CSO295'
#             agnName = 'RX_J1054.2+3511'


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
#             agnName = 'SDSSJ112632.90+120437.0'
#             agnName = 'SDSSJ112756.70+115427.0'


        # NGC5951
        if galaxyName == 'NGC5951':
            flipInclination = False
            reverse = True
            NFW_fit = 'tightest'
            agnName = '2E1530+1511'


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
            PA = 43.
        
        # UGC08146
        if galaxyName == 'UGC08146':
            flipInclination = False
            NFW_fit = 'tightester2'
            reverse = True
            agnName = 'PG1259+593'

        # --- grab the coordinates for this target
        RA_target = agn[agnName]['RAdeg']
        Dec_target = agn[agnName]['DEdeg']
        
        # --- approaching or receding side?
        side = agn[agnName]['side']


##########################################################################################
##########################################################################################
    # --- rename arrays something stupid
    vels = vrot_incCorrected_vals
    xvals = deepcopy(xVals)
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
    
    if fit_NFW:
#         reverse it?
#         if reverse:
#         xData.reverse()
#             vels.reverse()
#         yData = np.array(yData)*-1
    
    
        # --- fold the data over so approaching and receding sides are both positive
        newVals = []
        newX = []
        
        # --- these will be the 
        xData1 = []
        yData1 = []
        xData2 = []
        yData2 = []
        
        print 'vels: ',vels
        print
        print 'xVals: ',xVals
        print
        print
        
        for v,x in zip(vels, xVals):
            xData1.append(x)
            yData1.append(v)

            newX.append(abs(x))
            newVals.append(abs(v))
    
        newX = np.array(newX)[len(newX)/3:]
        newVals = np.array(newVals)[len(newVals)/3:]
        
        print 'newX: ',newX
        print
        print 'newVals: ',newVals
        print
        
        a = 3.95
        rho = 500.
        
        v200 = abs(right_vrot_incCorrected_avg) / 1.2
        c = 5
        r200 = R_vir
        
        r200_lowerbound = R_vir/1.2
        r200_upperbound = R_vir*1.2
        v200_lowerbound = abs(right_vrot_incCorrected_avg)/2.0
        v200_upperbound = abs(right_vrot_incCorrected_avg)
        c_lowerbound = 1
        c_upperbound = 25
        
        print 'r200_lowerbound: ',r200_lowerbound
        print 'r200_upperbound: ',r200_upperbound
        print 'v200_lowerbound: ',v200_lowerbound
        print 'v200_upperbound: ',v200_upperbound
        print 'R_vir: ',R_vir
        print
        
        
        # --- slightly tighter (e.g., for NGC2770)
#         if NFW_fit == 'tight':
#             r200_lowerbound = 10
#             r200_upperbound = 250
#             v200_lowerbound = 10
#             v200_upperbound = 130
#             c_lowerbound = 1
#             c_upperbound = 35
#         
#             v200 = 50
#             c = 10
#             r200 = R_vir
        

        # --- tighter bounds (e.g., for UGC09760)
#         if NFW_fit == 'tighter':
#             r200_lowerbound = 10
#             r200_upperbound = 250
#             v200_lowerbound = 10
#             v200_upperbound = 115
#             c_lowerbound = 1
#             c_upperbound = 35
#             
#             v200 = 50
#             c = 10
#             r200 = R_vir


        # --- tightest bounds (e.g., for UGC04238, UGC09760)
#         if NFW_fit == 'tightest':
#             r200_lowerbound = 10
#             r200_upperbound = 250
#             v200_lowerbound = 10
#             v200_upperbound = 90
#             c_lowerbound = 1
#             c_upperbound = 35
#             
#             v200 = 50
#             c = 10
#             r200 = R_vir
            
        # --- more than tightest bounds (e.g., for UGC08146)
#         if NFW_fit == 'tightester':
#             r200_lowerbound = 10
#             r200_upperbound = 250
#             v200_lowerbound = 10
#             v200_upperbound = 80
#             c_lowerbound = 1
#             c_upperbound = 35
#             
#             v200 = 50
#             c = 10
#             r200 = R_vir
            
        # for UGC04238
#         if NFW_fit == 'tightester2':
#             r200_lowerbound = 10
#             r200_upperbound = 250
#             v200_lowerbound = 10
#             v200_upperbound = 75
#             c_lowerbound = 1
#             c_upperbound = 35
#             
#             v200 = 50
#             c = 10
#             r200 = R_vir
            
        # for NGC4536
#         if NFW_fit == 'NGC4536':
#             r200_lowerbound = 20
#             r200_upperbound = 200
#             v200_lowerbound = 150
#             v200_upperbound = 350
#             c_lowerbound = 1
#             c_upperbound = 50
#             
#             v200 = 150
#             c = 30
#             r200 = R_vir
        
        try:
            print
            print 'NFW : ', NFW
            print 'newX : ',np.shape(newX)
            print 'newVals : ',np.shape(newVals)
            print
            popt, pcov = optimize.curve_fit(NFW, newX, newVals, p0=[v200,c,r200], \
            bounds=((v200_lowerbound, c_lowerbound, r200_lowerbound), (v200_upperbound, c_upperbound, r200_upperbound)))
            
        except Exception,e:
            print 'exception in curve_fit: ',e
            print
            sys.exit()
        
        print
        print 'popt: ',popt
        print
        print
        
        print 'now the fit. popt = {0}, pcov = {0}'.format(popt,pcov)
        print
        print 'np.sqrt(np.diag(pcov)) = ',np.sqrt(np.diag(pcov))
        print
        
        # plot it
#         xData_fit = linspace(0,max(xvals),num=1000)
#         y_fit = NFW(xData_fit,*popt)
        x_lim = int(max(newX)) + int(max(newX)/4.)
        fig = plot_NFW(newX, newVals, popt, x_lim)
        fig.savefig("{0}{1}_NFW_{2}.jpg".format(saveDirectory,galaxyName,x_lim),dpi=300,bbox_inches='tight')
        
#         x_lim = 500
        x_lim = round(R_vir + 10,0)
        print 'x_lim: ',x_lim
        fig = plot_NFW(newX, newVals, popt, x_lim)
        fig.savefig("{0}{1}_NFW_{2}.jpg".format(saveDirectory,galaxyName,x_lim),dpi=300,bbox_inches='tight')
        
        
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

    else:
    
        # this one is for 0 centered
        xvalStart = xvals[0]
        xvalEnd = xvals[-1]
 
        lowMean = mean(vels[:6])
        highMean = mean(vels[-6:])
    
        print 'lowMean: ',lowMean
        print 'highMean: ',highMean
        print
        print 'right_vrot_incCorrected_avg: ',right_vrot_incCorrected_avg
        print 'left_vrot_incCorrected_avg: ',left_vrot_incCorrected_avg
        print

        print 'vels[0]: ',vels[0]
        print 'vels[-1]: ',vels[-1]
        print
        
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
    
    # --- calculate impact parameter and shit
    impact = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist)
    
    # --- RA component of impact parameter - by setting the Dec to be the same for both
    impact_RA = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_target,Dec_galaxy,dist)
    
    # --- Dec component of impact parameter - by setting the RA to be the same for both
    impact_Dec = calculateImpactParameter(RA_galaxy,Dec_galaxy,RA_galaxy,Dec_target,dist)
    
    # --- calculate azimuth
    az = calculateAzimuth(RA_galaxy,Dec_galaxy,RA_target,Dec_target,dist,PA)

        
#     inc = .45
#     R_vir = 1.5*R_vir

#     if RA_galaxy > RA_target:
#         impact_RA = -impact_RA
    if Dec_galaxy > Dec_target:
        impact_Dec = -impact_Dec

    if RA_target > RA_galaxy:
        impact_RA = -impact_RA
            
    
    print
    print 'impact: ',impact
    print 'impact_RA: ',impact_RA
    print 'impact_Dec: ',impact_Dec
    print 'az: ',az
    print 'R_vir: ',R_vir
    print
        
    # --- inclination is backwards, so flip it
#     effectiveInc = 90.-inc
    effectiveInc = inc
    print 'effectiveInc: ',effectiveInc
    print
        
    if flipInclination:
        effectiveInc *=-1.
    
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
    pa_rot = (90. + PA) * np.pi/180.

    N_inc = rotate_vector(N, k_inc, inc_rot)
    print 'N_inc: ',N_inc
    print

    N_inc_pa = rotate_vector(N_inc, k_pa, pa_rot)
    print 'N_inc_pa: ',N_inc_pa
    print
    
    normal = N_inc_pa
    origin = np.array([0.,0.,0.])

    
    # Define ray -> [0,RA_dif,Dec_dif]
    rayDirection = np.array([1, 0, 0])

    rayPoint = np.array([0, impact_RA, impact_Dec])
    print 'rayPoint: ',rayPoint
    print

    
##########################################################################################
##########################################################################################
    # --- now loop through layers of galaxy planes
        
    zcutoff = zcutoffm * R_vir
    rcutoff = rcutoffm * R_vir
    print 'zcutoff: ',zcutoff
    print
    
    v_proj_list = []
    intersect_list = []
    d_plot_list = []
    intersect_point_list = []
    v_list = []
    v_90_list = []
    
    vlos_list = []
    dlos_list = []
    
    # --- how often to sample?
    if inc <= 60:
        s = 0.1
    elif inc <=80:
        s = 0.01
    elif inc <= 87:
        s = 0.005
    else:
        s = 0.0001
        
    if diskOnly:
        zcutoff = 0.1
        s = 0.1
        
        
##########################################################################################

    #  --- first apply the model to just the thin-disk midplane
    z = 0
#   midplane_model_projection = project_model_velocity(z, origin, normal, rayDirection, rayPoint, rcutoff, fit, sphericalHalo=False)
    midplane_model_projection = project_model_velocity(z, origin, normal, rayDirection, rayPoint, rcutoff, fit, sphericalHalo=False)

    # --- unpack the results
    midplane_intersect           = midplane_model_projection["intersect"]
    midplane_intersect_dist      = midplane_model_projection["intersect_dist"]
    midplane_dist_from_origin    = midplane_model_projection["dist_from_origin"]
    midplane_v_intersect         = midplane_model_projection["v_intersect"]
    midplane_v_n_intersect       = midplane_model_projection["v_n_intersect"]
    midplane_v_rotation          = midplane_model_projection["v_rotation"]
    midplane_v_proj              = midplane_model_projection["v_proj"]
    
    
    midplane_dlos = midplane_intersect[0]
    
    print 'midplane_intersect: ',midplane_intersect
    print 'midplane_intersect_dist: ',midplane_intersect_dist
    print 'midplane_dist_from_origin: ',midplane_dist_from_origin
    print 'midplane_v_intersect: ',midplane_v_intersect
    print 'midplane_v_n_intersect: ',midplane_v_n_intersect
    print 'midplane_v_rotation: ',midplane_v_rotation
    print 'midplane_v_proj: ',midplane_v_proj
    print 'origin: ',origin
    print 'normal: ',normal
    print 'rayDirection: ',rayDirection
    print 'rayPoint: ',rayPoint
    print
        
##########################################################################################

    # the NFW model
    for i in arange(-zcutoff, zcutoff, s):
        model_projection    = project_model_velocity(i, origin, normal, rayDirection, rayPoint, rcutoff, fit, sphericalHalo=False)
        intersect           = model_projection["intersect"]
        intersect_dist      = model_projection["intersect_dist"]
        dist_from_origin    = model_projection["dist_from_origin"]
        v_intersect         = model_projection["v_intersect"]
        v_n_intersect       = model_projection["v_n_intersect"]
        v_rotation          = model_projection["v_rotation"]
        v_proj              = model_projection["v_proj"]
        
        
        # --- Where is the zero point? Center either at the midplane intersect or the
        # --- galaxy systemic velocity
        if center_at_intersect:
            intersect_xaxis = intersect[0]
        else:
            intersect_xaxis = intersect[0] - midplane_dlos

        if v_proj:
            v_proj_list.append(v_proj)
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
    
#     for i in arange(-rcutoff, rcutoff, steidel_step):
    print 'min(intersect_list): ',min(intersect_list)
    print 'max(intersect_list): ',max(intersect_list)
    print 'intersect[0]: ',intersect[0]
    print
    
    # --- define the model range
    steidel_min = -600
    steidel_max = 600
#     steidel_min = min(intersect_list)
#     steidel_max = max(intersect_list)

    # --- check if the sign is right for vmax
    if side == 'receding':
        vmax = float(vmax) * -1
    
    
    for i in arange(steidel_min, steidel_max, steidel_step):
        # --- steidel_model(vmax, inclination, hv, impact, azimuth, Dlos)
        
        effective_r = np.sqrt(i**2 + intersect_dist**2)
        effective_vc = fit(effective_r)
#         print "effective_r = {0}, effective_vc = {1}".format(effective_r, effective_vc)
        
        # --- the model
#         vlos = steidel_model(effective_vc, effectiveInc, hv, impact, az, i)

        vlos = steidel_model(vmax, effectiveInc, steidel_hv, impact, az, i)
        print 'dlos: ',i
        print 'vlos: ',vlos
        print 'vmax: ',vmax
        print 'inc: ',inc
        print 'az: ',az
        print 'hv: ',steidel_hv
        print

        # --- center at midplane intersect or at galaxy systemic
        if center_at_intersect:
            intersect_xaxis = i + midplane_dlos
        else:
            intersect_xaxis = i

        vlos_list.append(vlos)
        dlos_list.append(intersect_xaxis)

# 
#     for i in arange(steidel_min, steidel_max, steidel_step):
#         # --- steidel_model(vmax, inclination, hv, impact, azimuth, Dlos)
#         
#         effective_r = np.sqrt(i**2 + intersect_dist**2)
#         effective_vc = fit(effective_r)
# #         print "effective_r = {0}, effective_vc = {1}".format(effective_r, effective_vc)
#         
#         # --- the model
# #         vlos = steidel_model(effective_vc, effectiveInc, hv, impact, az, i)
# 
#         vlos = steidel_model(vmax, effectiveInc+5, steidel_hv, impact, az, i)
#         print 'dlos: ',i
#         print 'vlos: ',vlos
#         print 'vmax: ',vmax
#         print 'inc: ',inc
#         print 'az: ',az
#         print 'hv: ',steidel_hv
#         print
# 
#         # --- center at midplane intersect or at galaxy systemic
#         if center_at_intersect:
#             intersect_xaxis = i + midplane_dlos
#         else:
#             intersect_xaxis = i
# 
#         vlos_list.append(vlos)
#         dlos_list.append(intersect_xaxis)
            
            
            
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

    # --- Define the extent of the plot boundary
    plotExtent = round(1.5*rcutoff,-1)
    
    if plotExtent <400:
        plotExtent = 400
    if plotExtent >400:
        plotExtent = 800

    print 'beginning plotExtent: ',plotExtent
    print 'rcutoff: ',rcutoff
    print
    print
    print
    print
    print
    print
    
    # how big to make the plotted cylinder?
    zHeight = zcutoff
    
    # plot velocity on the x-axis? No idea what happens if this is False...
    plot_x_velocity = False
    
    # how many steps to take while plotting. each step moves the sightline forward and 
    # populates the other with that result
    steps = 1
    
    # tranpose the list of intersects for plotting
    ip_xlist, ip_ylist, ip_zlist = np.array(intersect_point_list).transpose()
    
    # --- some colors
    color_blue = '#436bad'  # french blue
    color_red = '#ec2d01'   # tomato red

    for i in arange(steps):
        i +=1
        print 'i: ',i
        # initial figure
        fig = plt.figure(figsize=(12,8))
        
        # first plot the v_proj data
        ax = fig.add_subplot(1,2,1)
        
#         agnName = agnName.replace('_','\_')
#         fig.suptitle(r'$\rm {0} - {1}:~ {2} x {3}R_{{vir}}$'.format(galaxyName,agnName,zcutoffm,rcutoffm), fontsize=16)


        # number of actual intercept elements to take for each step
        len_step = len(intersect_list)/steps
        
        # --- x axis in kpc or km/s
        if plot_x_velocity:
            intersect_x_list = (np.array(intersect_list)/1000.)*hubbleConstant
        else:
            intersect_x_list = np.array(intersect_list)
        
        print 'len_step: ',len_step
#         print 'intersect_x_list[:i*len_step] vs v_proj_list[:i*len_step]: ',intersect_x_list[:i*len_step],' , ',v_proj_list[:i*len_step]
        
        # --- plot it
        ax.scatter(intersect_x_list[:i*len_step], v_proj_list[:i*len_step], color='black', s=10)
        ylabel(r'$\rm V_{proj} ~[km~s^{-1}]$')
        
        ##################################################################################
        # --- plot the Steidel model results
#         print "dlos_list[:i*len_step]: ",dlos_list[:i*len_step]
#         print
#         print "vlos_list[:i*len_step]: ",vlos_list[:i*len_step]
#         print
        ax.scatter(dlos_list, vlos_list, color=color_blue, s=15, label = 'Steidel Model')


        # --- decide on the xlabel
        if plot_x_velocity:
            xlabel(r'$\rm Intersect ~[km~s^{-1}]$')
        else:
            xlabel(r'$\rm Intersect ~[kpc]$')
    
        # --- Define the minimum and maximum extent for the x and y axes
        max_velocity = max(np.concatenate((v_proj_list, vlos_list)))
        min_velocity = min(np.concatenate((v_proj_list, vlos_list)))
        
        max_x = max(np.concatenate((intersect_x_list, dlos_list)))
        min_x = min(np.concatenate((intersect_x_list, dlos_list)))
        
        
        tick_num = 9.
        x_tick_num = 6.
        tick_spacing = round((max_velocity - min_velocity)/tick_num,-1)
        x_tick_spacing = round((max_x - min_x)/x_tick_num,-1)
        
        y_extent = max_velocity - min_velocity
        x_extent = max_x - min_x
        tick_spacing = round(y_extent + y_extent/4., -1)/tick_num
        x_tick_spacing = round(x_extent + x_extent/4., -1)/x_tick_num

        print 'tick_spacing: ',tick_spacing
        print 'max(v_proj_list): ',max(v_proj_list)
        print 'min(v_proj_list): ',min(v_proj_list)

#         xlim_pos = max_x + x_tick_spacing/2
#         xlim_neg = min_x - x_tick_spacing/2
        xlim_pos = steidel_max
        xlim_neg = steidel_min
        ylim_pos = max_velocity + tick_spacing/2
        ylim_neg = min_velocity - tick_spacing/2
        
        print '########################################'
        print 'xlim_pos: ',xlim_pos
        print 'xlim_neg: ',xlim_neg
        print 'ylim_pos: ',ylim_pos
        print 'ylim_neg: ',ylim_neg
        print '########################################'
        print
#         print 'intersect_x_list: ',intersect_x_list
#         print
#         print 'v_proj_list: ',v_proj_list
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

        if tick_spacing == 0.0:
            tick_spacing = round((max_velocity - min_velocity)/tick_num,5)
            print 'inside'
        
        print
        print
        print 'int(xlim_neg): ', int(xlim_neg)
        print 'int(xlim_pos): ', int(xlim_pos)
        print 'int(ylim_neg): ',int(ylim_neg)
        print 'int(ylim_pos): ',int(ylim_pos)
        print
        print 'x_tick_spacing: ',x_tick_spacing
        print 'tick_spacing: ',tick_spacing
        print 'xlim_neg: ',xlim_neg
        print 'xlim_pos: ',xlim_pos
        print
        print
        print
        print

        # x-axis
        majorLocator   = MultipleLocator(x_tick_spacing)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(x_tick_spacing/2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(tick_spacing)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(tick_spacing/2)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
#         ax.set_xbound(lower=int(xlim_neg), upper=int(xlim_pos))
#         xlim(int(xlim_neg), int(xlim_pos))
        ax.set_xbound(lower=math.floor(xlim_neg), upper=math.floor(xlim_pos))

        ylim(int(ylim_neg), int(ylim_pos))
        xlim(int(xlim_neg), int(xlim_pos))
        
        if tick_spacing <1:
            ylim(ylim_neg-ylim_neg/10, ylim_pos + ylim_pos/10)
        
        # next plot the 3d fig
        ax = fig.add_subplot(1,2,2,projection='3d')

    
##########################################################################################
##########################################################################################

        # plot the sightline
        print 'rayPoint: ',rayPoint
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
        ax.view_init(elev=10., azim=15)
#         ax.view_init(elev=10., azim=10)


        # reverse the X-axis ticks without actually reversing the x-axis
        if plotExtent <= 400:
            yticks((-400, -200, 0, 200, 400), (400, 200, 0, -200, -400))
            xticks((-400, -200, 0, 200, 400), (-400, -200, 0, 200, 400))
        else:
            yticks((-800, -400, 0, 400, 800), (800, 400, 0, -400, -800))
            xticks((-800, -400, 0, 400, 800), (-800, -400, 0, 400, 800))
        
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
    
        savefig("{0}{1}.jpg".format(saveDirectory, i), dpi=200, bbox_inches='tight')
        close(fig)
    
##########################################################################################
##########################################################################################
    # write out a summary txt file
    
    summary_file = open(saveDirectory + 'summary.txt','wt')
    
    summary_file.write('Galaxy name: {0}\n'.format(galaxyName))
    summary_file.write('Target name: {0}\n'.format(agnName))
    summary_file.write('R_vir: {0}\n'.format(R_vir))
    summary_file.write('Inclination: {0}\n'.format(inc))
    summary_file.write('Max V_rot / sin(i): {0}\n'.format(right_vrot_incCorrected_avg))
    summary_file.write('\n')
    summary_file.write('flipInclination: {0}\n'.format(flipInclination))
    summary_file.write('reverse: {0}\n'.format(reverse))
    summary_file.write('fit_NFW: {0}\n'.format(fit_NFW))
    summary_file.write('Fit popt: [V200, c, R200] = {0}\n'.format(popt))
    summary_file.write('\n')
    summary_file.write('Allowed Y-velocity range: [{0}, {1}]\n'.format(min(v_proj_list), max(v_proj_list)))
    summary_file.write('\n')
    summary_file.write('Allowed X-velocity range: [{0}, {1}]\n'.format(min(intersect_x_list), max(intersect_x_list)))
    summary_file.write('\n')
    summary_file.write('-----------------\n')
    summary_file.write('\n')
    summary_file.write('Steidel Model range: [{0}, {1}]\n'.format(min(vlos_list), max(vlos_list)))
    summary_file.write('\n')
    summary_file.write('Steidel Model x-range: [{0}, {1}]\n'.format(min(dlos_list), max(dlos_list)))

    # this part sums the x and y velocities. For example, if you're at -10 km/s along the 
    # sightline and the measured velocity is -10, the velocity at that point is -20 km/s.
    if plot_x_velocity:
        totals = []
        for x, y in zip(intersect_x_list, v_proj_list):
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

    
    # save pickle data?
    if save_data_pickle:
        pickle_filename = '{0}/{1}-{2}_model_NFW{3}.p'.format(saveDirectory, galaxyName, agnName, fit_NFW)
    
        pickle_file = open(pickle_filename,'wt')
        
        d = {}
        
        # velocity along the sightline
        d['intersect_x_list'] = intersect_x_list
        
        # projected rotation velocity
        d['v_proj_list'] = v_proj_list
        
        # now the sightline
        d['z_sightline'] = z_sightline
        d['y_sightline'] = y_sightline
        d['x_sightline'] = x_sightline
        
        # now the cylinder
        d['intersect'] = intersect
        d['R'] = R
        d['p0'] = p0
        d['p1'] = p1
        
        d['popt'] = popt

        pickle.dump(d, pickle_file)
        pickle_file.close()

    
    
if __name__ == '__main__':
    main()
    