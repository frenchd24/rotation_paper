#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: line_fit_prepare.py, v1 05/01/18

Takes all the spec-TARGET-HI-1215_VELOCITY files and creates:

spec-TARGET-HI-1215_VELOCTY_BINNING.txt, which only has wavelength, flux, and error
    columns, binning according to BINNING (probably 1, 2, 5)
    
Also creates a TARGET-VELOCTITY_BINNING.pars file for VoigtFit to use

    
'''

import sys
import os
import csv
# import string
import warnings
import numpy as np
# import atpy
import getpass
from utilities import *
import math
from astropy.io import ascii
from astropy import table

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################


    
def main():
    # check which computer we're on, and grab the appropriate file
    
    directory = '/Users/frenchd/Research/inclination/line_fits/'
    outDirectory = '/Users/frenchd/Research/inclination/line_fits_prepared2/'
#     outDirectory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/line_fits2/'

    files = os.listdir(directory)
    files2 = []
    
    # ignore hidden files
    for f in files:
        if f[0] != '.':
            files2.append(f)
            
    
    for f in files2:
        print 'f: ',f
        
        content = np.loadtxt(directory+f, unpack=True)
        
        f = f.replace('.txt','')
        
        print 'f: ',f
        print 'content: ',content[0]
        
        new_f1 = f+'_1.txt'
#         new_f2 = directory+f+'_2.txt'
#         new_f5 = directory+f+'_5.txt'
        
        np.savetxt(outDirectory+new_f1, np.transpose([content[0],content[1],content[2]]))
#         np.savetxt(new_f2, np.transpose([content[0],content[5],content[2]]))
#         np.savetxt(new_f5, np.transpose([content[0],content[6],content[10]]))
        
                
        # now the .pars files
        
        # first get the target name
        f_name = f.replace('spec-','').replace('-HI-1215','')
        
        # now the redshift
        velocity = f[-4:].replace('_','')
        
        c = 2.99792*10**5
        
        line_z = float(velocity) / c
        
        # write out the .pars file
        pars_file = open(outDirectory + f_name + '.pars', 'wt')
        
        pars_file.write('name : {0} \n'.format(f_name))
        pars_file.write('z_sys:  {0} \n'.format(line_z))
        pars_file.write('norm_method:  linear   # or spline \n')
        pars_file.write('save : \n')
        pars_file.write('\n')
        pars_file.write('\n')

        pars_file.write('# Load spectra: \n')
        pars_file.write('#     filename     spectral resolution in km/s\n')
        pars_file.write('data  {0}   16.9 \n'.format(new_f1))
        pars_file.write('\n')
        pars_file.write('\n')
        
        pars_file.write("# Include optional commands to the fit, e.g., rebin=2, method='nelder' ...\n")
        pars_file.write('fit-options rebin=2 \n')
        pars_file.write('\n')
    

        pars_file.write('# Continuum Fitting using Chebyshev Polynomials:\n')
        pars_file.write('C_order = 3 \n')
        pars_file.write('\n')
        pars_file.write('\n')
        
        
        pars_file.write('# Define the lines that should be included in the fit:\n')
        pars_file.write('# The default velocity span is 500 km/s but can specified\n')
        pars_file.write('# for each individual lines\n')
        pars_file.write('lines HI_1  velspan=500 \n')
        pars_file.write('\n')
        pars_file.write('# Define components in redshift space:\n')
        pars_file.write('#          ion   z    b   logN\n')
        pars_file.write('#component ___  ___  ___  ____\n')
        pars_file.write('component  HI  0.     20.  13.3   velocity var_z=True \n')
        pars_file.write('\n')
        pars_file.write('\n')
        
        pars_file.write('# Define components using the interactive mode:\n')
        pars_file.write('# interactive HI_1 \n')
        pars_file.write('\n')
        pars_file.write('\n')

        pars_file.write('# --- Output Commands: \n')
        pars_file.write('\n')
        pars_file.write('output velocity \n')
        pars_file.write('\n')
        pars_file.write('# Print total abundances for each ion: \n')
        pars_file.write('total \n')
        pars_file.write('\n')
        
        pars_file.close()
        
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
