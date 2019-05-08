#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: pickleSALT.py v 1.0 12/31/2017

Comes from: buildDataLists3.py, v 5.0 12/02/2015

Makes a pickle file with all the info from:
salt_sightlines_all_results_include.csv



'''

import sys
import os
import csv

from pylab import *
# import atpy
from math import *
from utilities import *
import getpass
import pickle

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from matplotlib import rc


###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    if getpass.getuser() == 'frenchd':
        filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/salt_sightlines_all_results_include.csv'
        pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/pickleSALT.p'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    theFile = open(filename,'rU')
    pickleFile = open(pickleFilename,'wt')
    
    reader = csv.DictReader(theFile)

    
    # overall structure: fullDict is a dictionary with all the lines and their data in it
    # separated into 'associated' and 'ambiguous' as the two keys. Associated contains
    # all the lists of data for lines associated with a galaxy. Ambiguous contains all
    # the lists of data for lines not unambiguously associated (could be many galaxies
    # or none)
    
    fullDict = {}
    
    # all the lists to be used for associated lines
    targetNameL = []
    centerL = []
    galaxyNameL = []
    environmentL = []
    RA_agnL = []
    Dec_agnL = []
    RA_galL = []
    Dec_galL = []
    likelihoodL = []
    likelihood_cusL = []
    virialRadiusL = []
    cusL = []
    impactParameterL = []
    redshiftDistancesL = []
    vcorrGalaxyL = []
    radialVelocityL = []
    vel_diffL = []
    distGalaxyL = []
    AGNL = []
    majorAxisL = []
    minorAxisL = []
    inclinationL = []
    positionAngleL = []
    azimuthL = []
    RC3flagL = []
    RC3typeL = []
    RC3incL = []
    RC3paL = []
    morphologyL = []
    final_morphologyL = []
    galaxyRedshiftL = []
    AGNredshiftL = []
    spectrumStatusL = []
    includeL = []
    include_virL = []
    include_customL = []
    Lya_vL = []
    vlimitsL = []
    Lya_WL = []
    NaL = []
    bL = []
    identifiedL = []
    commentL = []
    
    
    for l in reader:
        #grab all the values
        targetName = l['targetName']
        center = eval(l['center'])
        galaxyName = l['galaxyName']
        environment = eval(l['environment'])
        degreesJ2000RA_DecAGN = eval(l['degreesJ2000RA_DecAGN'])
        degreesJ2000RA_DecGalaxy = eval(l['degreesJ2000RA_DecGalaxy'])
        likelihood = l['likelihood']
        likelihood_cus = l['likelihood_1.5']
        virialRadius = l['virialRadius']
        cus = l['d^1.5']
        impactParameter = l['impactParameter (kpc)']
        redshiftDistances = l['redshiftDistances']
        vcorrGalaxy = l['vcorrGalaxy (km/s)']
        radialVelocity = l['radialVelocity (km/s)']
        vel_diff = l['vel_diff']
        distGalaxy = l['distGalaxy (Mpc)']
        AGN = l['AGN S/N']
        majorAxis = l['majorAxis (kpc)']
        minorAxis = l['minorAxis (kpc)']
        inclination = l['inclination (deg)']
        positionAngle = l['positionAngle (deg)']
        azimuth = l['azimuth (deg)']
        RC3flag = l['RC3flag']
        RC3type = l['RC3type']
        RC3inc = l['RC3inc (deg)']
        RC3pa = l['RC3pa (deg)']
        morphology = l['morphology']
        final_morphology = l['final_morphology']
        galaxyRedshift = l['galaxyRedshift']
        AGNredshift = l['AGNredshift']
        spectrumStatus = l['spectrumStatus']
        include = l['include']
        include_vir = l['include_vir']
        include_custom = l['include_custom']
        Lya_v = eval(l['Lya_v'])
        vlimits = eval(l['vlimits'])
        Lya_W = eval(l['Lya_W'])
        Na = eval(l['Na'])
        b = eval(l['b'])
        identified = l['identified']
        comment = l['comment']
        
        RA_agn, Dec_agn = degreesJ2000RA_DecAGN
        RA_gal, Dec_gal = degreesJ2000RA_DecGalaxy

        targetNameL.append(targetName)
        centerL.append(center)
        galaxyNameL.append(galaxyName)
        environmentL.append(environment)
        RA_agnL.append(RA_agn)
        Dec_agnL.append(Dec_agn)
        RA_galL.append(RA_gal)
        Dec_galL.append(Dec_gal)
        likelihoodL.append(likelihood)
        likelihood_cusL.append(likelihood_cus)
        virialRadiusL.append(virialRadius)
        cusL.append(cus)
        impactParameterL.append(impactParameter)
        redshiftDistancesL.append(redshiftDistances)
        vcorrGalaxyL.append(vcorrGalaxy)
        radialVelocityL.append(radialVelocity)
        vel_diffL.append(vel_diff)
        distGalaxyL.append(distGalaxy)
        AGNL.append(AGN)
        majorAxisL.append(majorAxis)
        minorAxisL.append(minorAxis)
        inclinationL.append(inclination)
        positionAngleL.append(positionAngle)
        azimuthL.append(azimuth)
        RC3flagL.append(RC3flag)
        RC3typeL.append(RC3type)
        RC3incL.append(RC3inc)
        RC3paL.append(RC3pa)
        morphologyL.append(morphology)
        final_morphologyL.append(final_morphology)
        galaxyRedshiftL.append(galaxyRedshift)
        AGNredshiftL.append(AGNredshift)
        spectrumStatusL.append(spectrumStatus)
        includeL.append(include)
        include_virL.append(include_vir)
        include_customL.append(include_custom)
        Lya_vL.append(Lya_v)
        vlimitsL.append(vlimits)
        Lya_WL.append(Lya_W)
        NaL.append(Na)
        bL.append(b)
        identifiedL.append(identified)
        commentL.append(comment)


    # populate the dictionary

        
    fullDict['targetName'] = targetNameL
    fullDict['center'] = centerL
    fullDict['galaxyName'] = galaxyNameL
    fullDict['environment'] = environmentL
    fullDict['RA_agn'] = RA_agnL
    fullDict['Dec_agn'] = Dec_agnL
    fullDict['RA_gal'] = RA_galL
    fullDict['Dec_gal'] = Dec_galL
    fullDict['likelihood'] = likelihoodL
    fullDict['likelihood_cus'] = likelihood_cusL
    fullDict['virialRadius'] = virialRadiusL
    fullDict['cus'] = cusL
    fullDict['impact'] = impactParameterL
    fullDict['redshiftDistances'] = redshiftDistancesL
    fullDict['vcorr'] = vcorrGalaxyL
    fullDict['radialVelocity'] = radialVelocityL
    fullDict['vel_diff'] = vel_diffL
    fullDict['distGalaxy'] = distGalaxyL
    fullDict['AGNsn'] = AGNL
    fullDict['majorAxis'] = majorAxisL
    fullDict['minorAxis'] = minorAxisL
    fullDict['inclination'] = inclinationL
    fullDict['PA'] = positionAngleL
    fullDict['azimuth'] = azimuthL
    fullDict['RC3flag'] = RC3flagL
    fullDict['RC3type'] = RC3typeL
    fullDict['RC3inc'] = RC3incL
    fullDict['RC3pa'] = RC3paL
    fullDict['morphology'] = morphologyL
    fullDict['final_morphology'] = final_morphologyL
    fullDict['galaxyRedshift'] = galaxyRedshiftL
    fullDict['AGNredshift'] = AGNredshiftL
    fullDict['spectrumStatus'] = spectrumStatusL
    fullDict['include'] = includeL
    fullDict['include_vir'] = include_virL
    fullDict['include_custom'] = include_customL
    fullDict['Lya_v'] = Lya_vL
    fullDict['vlimits'] = vlimitsL
    fullDict['Lya_W'] = Lya_WL
    fullDict['Na'] = NaL
    fullDict['b'] = bL
    fullDict['identified'] = identifiedL
    fullDict['comment'] = commentL
    

##########################################################################################
##########################################################################################
##########################################################################################

    
    pickle.dump(fullDict,pickleFile)
    pickleFile.close()
    theFile.close()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    