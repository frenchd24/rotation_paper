#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotAzimuthMap_rotation.py, v 1.0 1/01/18



Take from: plotAzimuthMap2.py, v 2.1 1/20/17

Plot absorption properties on a map showing their distribution around a central galaxy

v2: Separate blue and red into separate plots, so now it's a plane of 6. Size points by 
    EW (10/06/2016)
    
    
v2.1: Use open symbols for rho <= R_vir to be consistent with the other plots. -> /plots6/
    (1/20/17)

'''

import sys
import os
import csv
from scipy import stats

from pylab import *
# import atpy
from math import *
from utilities import *
import getpass
import pickle
from matplotlib.patches import Ellipse


# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})

fontScale = 18
rc('text', usetex=True)
rc('font', size=18, family='serif', weight='normal')
rc('xtick.major',size=8,width=0.6)
rc('xtick.minor',size=5,width=0.6)
rc('ytick.major',size=8,width=0.6)
rc('ytick.minor',size=5,width=0.6)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
rc('axes',labelsize = fontScale)
rc('xtick', labelsize = fontScale)
rc('ytick',labelsize = fontScale)
# rc('font', weight = 450)
# rc('axes',labelweight = 'bold')
rc('axes',linewidth = 1)

    

###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    if getpass.getuser() == 'frenchd':
#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/Users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/pilot_paper_code/plots6/'
#         WS09data = '/Users/frenchd/Research/inclination/git_inclination/WS2009_lya_data.tsv'
        
#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/pickleSALT.p'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/figures/'

#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALTcut.p'
        pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'
        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/'


    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # use the old pickle file to get the full galaxy dataset info
    pickleFile = open(pickleFilename,'rU')
    fullDict = pickle.load(pickleFile)
    
    pickleFile.close()
    
    
    # save each plot?
    save = False
    
#     results = open(resultsFilename,'rU')
#     reader = csv.DictReader(results)
    
#     WS = open(WS09data,'rU')
#     WSreader = csv.DictReader(WS,delimiter=';')
    
    virInclude = False
    cusInclude = False
    finalInclude = 1
    
    maxEnv = 3000
    minL = 0.001
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # all the lists to be used for associated lines
    raList = []
    decList = []
    lyaVList = []
    lyaWList = []
    lyaErrList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    incList = []
    fancyIncList = []
    cosIncList = []
    cosFancyIncList = []
    paList = []
    vcorrList = []
    majList = []
    difList = []
    envList = []
    morphList = []
    m15List = []
    virList = []
    likeList = []
    likem15List = []
    
    
    # WS lists
#     WSvcorr = []
#     WSdiam = []
#     WSimpact =[]
#     WSew = []
#     WSvel = []
#     WSlya = []
#     WSvel_dif = []
#     WSvir = []
#     WSlike = []
#     
#     l_min = 0.001
# 
#     for w in WSreader:
#         vcorr = w['HV']
#         diam = w['Diam']
#         rho = w['rho']
#         ew = w['EWLya']
#         vel = w['LyaVel']
#         lya = w['Lya']
#         
#         if lya == 'Lya  ' and isNumber(diam) and isNumber(ew) and isNumber(rho):
#             if float(rho) <=500.0:
#                 # this is a single galaxy association
#                 vir = calculateVirialRadius(float(diam))
#                 
#                 vel_dif = float(vcorr) - float(vel)
#     
#                 # try this "sphere of influence" value instead
#                 m15 = float(diam)**1.5
# 
#                 # first for the virial radius
#                 likelihood = math.exp(-(float(rho)/vir)**2) * math.exp(-(vel_dif/200.)**2)
#                 
#                 if vir>= float(rho):
#                     likelihood = likelihood*2
#                     
#                 # then for the second 'virial like' m15 radius
#                 likelihoodm15 = math.exp(-(float(rho)/m15)**2) * math.exp(-(vel_dif/200.)**2)
#                 
#                 if m15>= float(rho):
#                     likelihoodm15 = likelihoodm15*2
#                     
#                 if likelihood <= likelihoodm15:
#                     likelihood = likelihoodm15
#                     
#                 WSlike.append(likelihood)
#                 
# #                 l_min=0
#                 
#                 if likelihood >= l_min:
#                 
#                     WSvcorr.append(float(vcorr))
#                     WSdiam.append(float(diam))
#                     WSvir.append(vir)
#                     WSimpact.append(float(rho))
#                     WSew.append(float(ew))
#                     WSvel.append(float(vel))
#                     WSlya.append(lya)
#                     WSvel_dif.append(vel_dif)
    
    
    
    targetNameL= fullDict['targetName']
    galaxyNameL = fullDict['galaxyName']
    environmentL = fullDict['environment']
    RA_agnL = fullDict['RA_agn']
    Dec_agnL = fullDict['Dec_agn']
    RA_galL = fullDict['RA_gal']
    Dec_galL = fullDict['Dec_gal']
    likelihoodL = fullDict['likelihood']
    likelihood_cusL = fullDict['likelihood_cus']
    virialRadiusL = fullDict['virialRadius']
    cusL = fullDict['cus']
    impactParameterL = fullDict['impact']
    vcorrL = fullDict['vcorr']
    radialVelocityL = fullDict['radialVelocity']
    vel_diffL = fullDict['vel_diff']
    distGalaxyL = fullDict['distGalaxy']
    majorAxisL = fullDict['majorAxis']
    minorAxisL = fullDict['minorAxis']
    inclinationL = fullDict['inclination']
    positionAngleL = fullDict['PA']
    azimuthL = fullDict['azimuth']
    RC3flagL = fullDict['RC3flag']
    RC3typeL = fullDict['RC3type']
    RC3incL = fullDict['RC3inc']
    RC3paL = fullDict['RC3pa']
    final_morphologyL = fullDict['final_morphology']
    includeL = fullDict['include']
    include_virL = fullDict['include_vir']
    include_customL = fullDict['include_custom']
    Lya_vL = fullDict['Lya_v']
    vlimitsL = fullDict['vlimits']
    Lya_WL = fullDict['Lya_W']
    NaL = fullDict['Na']
    bL = fullDict['b']
    identifiedL = fullDict['identified']
    
    print 'initial len(Lya_vL): ',len(Lya_vL)
    print

    i = -1
    for inc,inc_vir,inc_cus in zip(includeL,include_virL,include_customL):
        i+=1
        go = False
        if match:
            if virInclude == inc_vir and cusInclude == inc_cus:
                go = True
            else:
                go = False
                
        else:
            if virInclude and inc_vir:
                go = True
                
            elif cusInclude and inc_cus:
                go = True
                
            elif finalInclude and inc:
                go = True
            
            else:
                go = False
        
        if go:
            RA_agn = RA_agnL[i]
            Dec_agn = Dec_agnL[i]
            RA_gal = RA_galL[i]
            Dec_gal = Dec_galL[i]
            lyaV = Lya_vL[i]
            lyaW = Lya_WL[i]
            lyaW_err = lyaW*0.1
            env = environmentL[i]
            galaxyName = galaxyNameL[i]
            impact = impactParameterL[i]
            galaxyDist = distGalaxyL[i]
            pa = positionAngleL[i]
            RC3pa = RC3paL[i]
            morph = final_morphologyL[i]
            vcorr = vcorrL[i]
            maj = majorAxisL[i]
            minor = minorAxisL[i]
            inc = inclinationL[i]
            az = azimuthL[i]
            b = bL[i]
            b_err = b*0.1
            na = NaL[i]
            na_err = na*0.1
            likelihood = likelihoodL[i]
            likelihoodm15 = likelihood_cusL[i]
            virialRadius = virialRadiusL[i]
            m15 = cusL[i]
            vel_diff = vel_diffL[i]
            
            if isNumber(RC3pa) and not isNumber(pa):
                pa = RC3pa
            
            if isNumber(inc):
                cosInc = cos(float(inc) * pi/180.)
                
                if isNumber(maj) and isNumber(minor):
                    q0 = 0.2
                    fancyInc = calculateFancyInclination(maj,minor,q0)
                    cosFancyInc = cos(fancyInc * pi/180)
                else:
                    fancyInc = -99
                    cosFancyInc = -99
            else:
                cosInc = -99
                inc = -99
                fancyInc = -99
                cosFancyInc = -99
            
            # all the lists to be used for associated lines
            if float(env) <= maxEnv and float(likelihood) >= minL:
                raList.append(RA_gal)
                decList.append(Dec_gal)
                lyaVList.append(float(lyaV))
                lyaWList.append(float(lyaW))
                lyaErrList.append(float(lyaW_err))
                naList.append(na)
                bList.append(float(b))
                impactList.append(float(impact))
                print az
                azList.append(float(az))
                incList.append(float(inc))
                fancyIncList.append(fancyInc)
                cosIncList.append(cosInc)
                cosFancyIncList.append(cosFancyInc)
                paList.append(pa)
                vcorrList.append(vcorr)
                majList.append(maj)
                difList.append(float(vel_diff))
                envList.append(float(env))
                morphList.append(morph)
                m15List.append(m15)
                virList.append(virialRadius)
                likeList.append(likelihood)
                likem15List.append(likelihoodm15)


#     results.close()
#     WS.close()
            
    # lists for the full galaxy dataset

    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    
    print 'len(lyaVList): ',len(lyaVList)

########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_azimuthMap = True
    save = False
    
    if plot_azimuthMap:
        
        # azimuths for each inclination interval
        azBlue_1 = []
        azRed_1 = []
        azBlue_2 = []
        azRed_2 = []
        azBlue_3 = []
        azRed_3 = []
        
        # scaled sizes for EW-scaled markers
        size_1b = []
        size_1r = []
        size_2b = []
        size_2r = []
        size_3b = []
        size_3r = []
        
        # new coordinates of absorption for plotting (w/ galaxy at (0,0))
        yList_1r = []
        xList_1r = []
        yList_2r = []
        xList_2r = []
        yList_3r = []
        xList_3r = []
        
        yList_1b = []
        xList_1b = []
        yList_2b = []
        xList_2b = []
        yList_3b = []
        xList_3b = []
        
        # save rho/R_vir also
        rR_1r = []
        rR_2r = []
        rR_3r = []
        
        rR_1b = []
        rR_2b = []
        rR_3b = []
        
        largestEW = max(lyaWList)
        smallestEW = min(lyaWList)
        maxSize = 400
        minSize = 20
                
        # calculate the position on the sky for each absorption feature wrt to the galaxy
        for r,d,i,a,fInc,vir,dif,ew in zip(raList,decList,impactList,azList,fancyIncList,virList,difList,lyaWList):
            if float(a) >= 0 and float(fInc)>=0 and float(fInc)<=90:

                # y coordinate
                y = (float(i)/float(vir)) * sin((a*pi)/180.)
            
                # x coordinate
                x = (float(i)/float(vir)) * cos(a*pi/180.)
                
                # rho / Rvir
                rhoRvir = float(i)/float(vir)
                
                # new size for the marker point
                newSize = ((float(ew) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
                    
                if fInc <=40:
                    if dif >0:
                        #blue absorber
                        azBlue_1.append(a)
                        xList_1b.append(x)
                        yList_1b.append(y)
                        size_1b.append(newSize)
                        rR_1b.append(rhoRvir)

                    else:
                        azRed_1.append(a)
                        xList_1r.append(x)
                        yList_1r.append(y)
                        size_1r.append(newSize)
                        rR_1r.append(rhoRvir)

                if fInc > 40 and fInc <=65:
                    if dif >0:
                        #blue absorber
                        azBlue_2.append(a)
                        yList_2b.append(y)
                        xList_2b.append(x)
                        size_2b.append(newSize)
                        rR_2b.append(rhoRvir)

                    else:
                        azRed_2.append(a)
                        yList_2r.append(y)
                        xList_2r.append(x)
                        size_2r.append(newSize)
                        rR_2r.append(rhoRvir)

                if fInc > 65:
                    if dif >0:
                        #blue absorber
                        azBlue_3.append(a)
                        yList_3b.append(y)
                        xList_3b.append(x)
                        size_3b.append(newSize)
                        rR_3b.append(rhoRvir)

                    else:
                        azRed_3.append(a)
                        yList_3r.append(y)
                        xList_3r.append(x)
                        size_3r.append(newSize)
                        rR_3r.append(rhoRvir)

            else:
                print 'float(a) <0: ',r,d,i,a,fInc

        # calculate the average red vs blue azimuth line
        blueAvg1 = mean(azBlue_1)
        print 'blueAvg1: ',blueAvg1
        redAvg1 = mean(azRed_1)
        print 'redAvg1: ',redAvg1
        print
        blueAvg2 = mean(azBlue_2)
        print 'blueAvg2: ',blueAvg2
        redAvg2 = mean(azRed_2)
        print 'redAvg2: ',redAvg2
        print
        blueAvg3 = mean(azBlue_3)
        print 'blueAvg3: ',blueAvg3
        redAvg3 = mean(azRed_3)
        print 'redAvg3: ',redAvg3
        print
        
        xyBlueAvg1 = (500.* cos(blueAvg1 * pi/180.), 500.* sin(blueAvg1 * pi/180.))
        xyRedAvg1 = (500.* cos(redAvg1 * pi/180.), 500.* sin(redAvg1 * pi/180.))
        print 'xyBlueAvg1: ',xyBlueAvg1
        print 'xyRedAvg1: ',xyRedAvg1
        xyBlueAvg2 = (500.* cos(blueAvg2 * pi/180.), 500.* sin(blueAvg2 * pi/180.))
        xyRedAvg2 = (500.* cos(redAvg2 * pi/180.), 500.* sin(redAvg2 * pi/180.))
        print 'xyBlueAvg2: ',xyBlueAvg2
        print 'xyRedAvg2: ',xyRedAvg2
        xyBlueAvg3 = (500.* cos(blueAvg3 * pi/180.), 500.* sin(blueAvg3 * pi/180.))
        xyRedAvg3 = (500.* cos(redAvg3 * pi/180.), 500.* sin(redAvg3 * pi/180.))
        print 'xyBlueAvg3: ',xyBlueAvg3
        print 'xyRedAvg3: ',xyRedAvg3

##########################################################################################

        # plot the distributions
        fig = figure(figsize=(14.5,11))
        alpha = 0.6
        bSymbol = 'D'
        rSymbol = 'o'
    
        ax1 = fig.add_subplot(2,3,1)
        subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.001,hspace=0.001)

        # Low inclinations = circular mostly
        e = Ellipse(xy=(0,0), width=1.0, height=0.9, angle=0)
#         legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
                    
        # no transparency
        e.set_alpha(0.3)
        ax1.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label1 = r'$\rm Inc \leq 40$'
        
#         y coordinate
#         y = (float(i)/float(vir)) * sin((a*pi)/180.)
#     
#         x coordinate
#         x = (float(i)/float(vir)) * cos(a*pi/180.)
        
        for x, y, a, s,r in zip(xList_1b,yList_1b,azBlue_1,size_1b,rR_1b):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if r <=1:
                fc = 'none'
                ec = 'blue'
                ax1.scatter(x,y,c='blue',marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = 'blue'
                ec = 'black'
                ax1.scatter(x,y,c='blue',marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)

#         ax1.scatter(xList_1r,yList_1r,c='red',marker=rSymbol,alpha=alpha,s=50)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax1.text(0.715, 0.93, label1, transform=ax1.transAxes, fontsize=15,verticalalignment='top', bbox=props)

        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax1.get_xticklabels(), visible=False)
    
        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)

        xlim(0,2.0)
        ylim(0,2.5)
        ylabel(r'$\rm \rho / R_{vir}$')

#         ax1.get_yaxis().set_tick_params(which='both', direction='out')

##########################################################################################
        ax2 = fig.add_subplot(2,3,2,sharey=ax1)
        
        # medium inclination = ellipse
        e = Ellipse(xy=(0,0), width=1.0, height=0.6, angle=0)
                    
        # no transparency
        e.set_alpha(0.3)
        ax2.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label2 = r'$\rm 40 < Inc \leq 65$'
        
        for x, y, a, s, r in zip(xList_2b,yList_2b,azBlue_2,size_2b,rR_2b):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if r <=1:
                fc = 'none'
                ec = 'blue'
                ax2.scatter(x,y,c='blue',marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = 'blue'
                ec = 'black'
                ax2.scatter(x,y,c='blue',marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
        
#         ax2.scatter(xList_2b,yList_2b,c='blue',marker=bSymbol,alpha=alpha,s=size_2b)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax2.text(0.58, 0.93, label2, transform=ax2.transAxes, fontsize=15,verticalalignment='top', bbox=props)
    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax2.yaxis.set_major_locator(majorLocator)
        ax2.yaxis.set_major_formatter(majorFormatter)
        ax2.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax2.get_xticklabels(), visible=False)

        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax2.yaxis.set_major_locator(majorLocator)
        ax2.yaxis.set_major_formatter(majorFormatter)
        ax2.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax2.get_xticklabels(), visible=False)

        xlim(0,2.0)
        ylim(0,2.5)
        plt.setp(ax2.get_yticklabels(), visible=False)
        xlabel(r'$\rm \rho / R_{vir}$')

#         ax2.yaxis.tick_right()
#         ax2.get_yaxis().set_tick_params(which='both', direction='out')

##########################################################################################
        ax3 = fig.add_subplot(2,3,3,sharey=ax1)
        
        # plot the flat galaxy line
        e = Ellipse(xy=(0,0), width=1.0, height=0.3, angle=0)
                    
        # no transparency
        e.set_alpha(0.3)
        ax3.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label3 = r'$\rm Inc > 65$'
        
        for x, y, a, s, r in zip(xList_3b,yList_3b,azBlue_3,size_3b,rR_3b):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if r <=1:
                fc = 'none'
                ec = 'blue'
                ax3.scatter(x,y,c='blue',marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = 'blue'
                ec = 'black'
                ax3.scatter(x,y,c='blue',marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
        
#         ax3.scatter(xList_3b,yList_3b,c='blue',marker=bSymbol,alpha=alpha,s=size_3b)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax3.text(0.715, 0.93, label3, transform=ax3.transAxes, fontsize=15,verticalalignment='top', bbox=props)
    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax3.yaxis.set_major_locator(majorLocator)
        ax3.yaxis.set_major_formatter(majorFormatter)
        ax3.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax3.get_xticklabels(), visible=False)

        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax3.yaxis.set_major_locator(majorLocator)
        ax3.yaxis.set_major_formatter(majorFormatter)
        ax3.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax3.get_yticklabels(), visible=False)

        xlim(0,2.0)
        ylim(0,2.5)
#             xlabel(r'$\rm \rho / R_{vir}$')


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

        ax4 = fig.add_subplot(2,3,4,sharex=ax1)

        # Low inclinations = circular mostly
        e = Ellipse(xy=(0,0), width=1.0, height=0.9, angle=0)
#         legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
                    
        # no transparency
        e.set_alpha(0.3)
        ax4.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label1 = r'$\rm Inc \leq 40$'

        for x, y, a, s, r in zip(xList_1r,yList_1r,azRed_1,size_1r,rR_1r):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if r <=1:
                fc = 'none'
                ec = 'red'
                ax4.scatter(x,y,c='red',marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = 'red'
                ec = 'black'
                ax4.scatter(x,y,c='red',marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)

#         ax4.scatter(xList_1r,yList_1r,c='red',marker=rSymbol,alpha=alpha,s=size_1r)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax4.text(0.715, 0.93, label1, transform=ax4.transAxes, fontsize=15,verticalalignment='top', bbox=props)

        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax4.yaxis.set_major_locator(majorLocator)
        ax4.yaxis.set_major_formatter(majorFormatter)
        ax4.yaxis.set_minor_locator(minorLocator)
        ax4.set_xticks([0.0,0.5,1.0,1.5])
    
        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax4.yaxis.set_major_locator(majorLocator)
        ax4.yaxis.set_major_formatter(majorFormatter)
        ax4.yaxis.set_minor_locator(minorLocator)
        ax4.set_yticks([0.0,0.5,1.0,1.5,2.0])

        xlim(0,2.0)
        ylim(0,2.5)
        ylabel(r'$\rm \rho / R_{vir}$')

#         ax1.get_yaxis().set_tick_params(which='both', direction='out')

##########################################################################################
        ax5 = fig.add_subplot(2,3,5,sharey=ax4,sharex=ax2)
        
        # medium inclination = ellipse
        e = Ellipse(xy=(0,0), width=1.0, height=0.6, angle=0)
                    
        # no transparency
        e.set_alpha(0.3)
        ax5.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label5 = r'$\rm 40 < Inc \leq 65$'
#         ax5.scatter(xList_2b,yList_2b,c='blue',marker=bSymbol,alpha=alpha,s=50)

        for x, y, a, s, r in zip(xList_2r,yList_2r,azRed_2,size_2r,rR_2r):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if rhoRvir <=1:
                fc = 'none'
                ec = 'red'
                ax5.scatter(x,y,c='red',marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = 'red'
                ec = 'black'
                ax5.scatter(x,y,c='red',marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)

#         ax5.scatter(xList_2r,yList_2r,c='red',marker=rSymbol,alpha=alpha,s=size_2r)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax5.text(0.58, 0.93, label5, transform=ax5.transAxes, fontsize=15,verticalalignment='top', bbox=props)
    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax5.yaxis.set_major_locator(majorLocator)
        ax5.yaxis.set_major_formatter(majorFormatter)
        ax5.yaxis.set_minor_locator(minorLocator)
        ax5.set_xticks([0.0,0.5,1.0,1.5])

        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax5.yaxis.set_major_locator(majorLocator)
        ax5.yaxis.set_major_formatter(majorFormatter)
        ax5.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax5.get_yticklabels(), visible=False)
        ax5.set_yticks([0.0,0.5,1.0,1.5,2.0])

        xlim(0,2.0)
        ylim(0,2.5)
        xlabel(r'$\rm \rho / R_{vir}$')

#         ax2.yaxis.tick_right()
#         ax2.get_yaxis().set_tick_params(which='both', direction='out')

##########################################################################################
        ax6 = fig.add_subplot(2,3,6,sharex=ax3,sharey=ax5)
        # plot the flat galaxy line
        e = Ellipse(xy=(0,0), width=1.0, height=0.3, angle=0)
                    
        # no transparency
        e.set_alpha(0.3)
        ax6.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label6 = r'$\rm Inc > 65$'
#         ax6.scatter(xList_3b,yList_3b,c='blue',marker=bSymbol,alpha=alpha,s=50)

        for x, y, a, s, r in zip(xList_3r,yList_3r,azRed_3,size_3r,rR_3r):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)

            if r <=1:
                fc = 'none'
                ec = 'red'
                ax6.scatter(x,y,c='red',marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = 'red'
                ec = 'black'
                ax6.scatter(x,y,c='red',marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)

#         ax6.scatter(xList_3r,yList_3r,c='red',marker=rSymbol,alpha=alpha,s=size_3r)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax6.text(0.715, 0.93, label6, transform=ax6.transAxes, fontsize=15,verticalalignment='top', bbox=props)
    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax6.yaxis.set_major_locator(majorLocator)
        ax6.yaxis.set_major_formatter(majorFormatter)
        ax6.yaxis.set_minor_locator(minorLocator)
        ax6.set_xticks([0.0,0.5,1.0,1.5,2.0])
    
        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax6.yaxis.set_major_locator(majorLocator)
        ax6.yaxis.set_major_formatter(majorFormatter)
        ax6.yaxis.set_minor_locator(minorLocator)
        ax6.set_yticks([0.0,0.5,1.0,1.5,2.0])
        plt.setp(ax6.get_yticklabels(), visible=False)

        xlim(0,2.0)
        ylim(0,2.5)
#             xlabel(r'$\rm \rho / R_{vir}$')

        if save:
            savefig('{0}/azimuthMap_separate2.pdf'.format(saveDirectory),\
            format='pdf',bbox_inches='tight')
        else:
            show()


###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    