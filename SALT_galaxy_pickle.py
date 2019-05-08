#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  SALT_galaxy_pickle.py, v1.0 02/01/18

Pickle all the data for each SALT galaxy


NOT FINISHED

'''


import csv

# from pylab import *
# from math import *
from utilities import *
# from scipy import stats
import getpass
# import math
import pickle
import json
# import io
import numpy as np



def main():
    
    # make a dictionary of all the info
    SALTpickle = {}
    
    # CGCG039-137
    target = 'RX_J1121.2+0326'
    RA_target = agn[target]['RAdeg']
    Dec_target = agn[target]['DEdeg']
    
    RA_galaxy = 170.36231
    Dec_galaxy = 3.44491
    dist = 101.21
    majDiam = 26.35
    
    impact = 98.9
    R_vir = 166.09
    az = 71.
    inc = 63.
    PA = 157.
    
    SALTpickle[
    
    # ESO343-G014
    target = 'RBS1768'
    RA_target = 324.7079167
    Dec_target = -38.47777778
    
    RA_galaxy = 324.43825
    Dec_galaxy = -38.49256
    dist = 126.07
    majDiam = 45.23
    
    impact = 465.6
    az = 74.
    inc = 89.9
    PA = 158.
    
    # IC5325
    target = 'RBS2000'
    RA_target = 351.18625
    Dec_target = -40.68027778
    
    RA_galaxy = 352.18096
    Dec_galaxy = -41.33347
    dist = 18.1
    majDiam = 20.45
    
    impact = 314.335827
    az = 64.1
    inc = 25.
    PA = 15.
    
    # MCG-03-58-009
    target = 'MRC2251-178'
    RA_target = 343.5245833
    Dec_target = -17.58194444
    
    RA_galaxy = 343.42021
    Dec_galaxy = -17.47889
    dist = 142.
    majDiam = 75.31
    
    impact = 355.0699641
    az = 71.
    inc = 49.
    PA = 27.
    
    # NGC1566
    target = '1H0419-577'
    RA_target = 66.50291667
    Dec_target = -57.20055556
    
    RA_galaxy = 65.00175
    Dec_galaxy = -54.93781
    dist = 7.19
    majDiam = 15.34
    
    impact = 302.77587
    az = 9.8
    inc = 48.
    PA = 170.

    # NGC1566
    target = 'HE0429-5343'
    RA_target = 67.66666667
    Dec_target = -53.61555556
    
    RA_galaxy = 65.00175
    Dec_galaxy = -54.93781
    dist = 7.19
    majDiam = 15.34
    
    impact = 256.2063291
    az = 60.1
    inc = 48.
    PA = 170.

    # NGC1566
    target = 'RBS567'
    RA_target = 69.91125
    Dec_target = -53.19194444
    
    RA_galaxy = 65.00175
    Dec_galaxy = -54.93781
    dist = 7.19
    majDiam = 15.34
    
    impact = 422.6192722
    az = 69.3
    inc = 48.
    PA = 170.