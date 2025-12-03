from __future__ import print_function
from __future__ import division
######### EPACRIS:autoT.py ######################################################
#                                                                               #
# Routine to plot and compare temperature profiles                              #
#                                                                               #
# v1: written by M.Scheucher 03/2022                                            #
#                                                                               #
#################################################################################

from builtins import str
from builtins import range
from past.utils import old_div
import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator,LogLocator,AutoMinorLocator)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re
from datetime import datetime
from scipy.stats import binned_statistic

#--- Overwriting Hatch-Linewidth in backend python.lib -------------------
import matplotlib.backends.backend_pdf
#from matplotlib import six
from matplotlib.backends.backend_pdf import Name, Op
from matplotlib.transforms import Affine2D
from matplotlib.path import Path



######### USER INTERFACE  #################################################
#-------- Define Output file name -----------------------------------------
#outfile = 'potTemp-testing' #Name of outfile
#outfile = 'TOI712b' #Name of outfile
#outfile = 'TOI-comp' #Name of outfile
#outfile = 'Jupiter-test' #Name of outfile
#outfile = 'GJ436-guillot' #Name of outfile
#outfile = 'GJ1214-RC' #Name of outfile
#outfile = 'K218b-RC' #Name of outfile
#outfile = 'GJ1214-climate' #Name of outfile
#outfile = 'K218b-climate' #Name of outfile
outfile = 'RC_debug-climate' #Name of outfile
outfile = 'RC_debug-GJ1214' #Name of outfile
outfile = 'RC_valid-Jupiter' #Name of outfile
#-------- Define Range of simulation folders ------------------------------
simu = [
##        'INPUTS/gen_e8',
#        'INPUTS/gen_e9',
#        'INPUTS/Jupiter',
#        'INPUTS/Earth1976',
#        'INPUTS/Earth1986',

#        'COMPARISON/Helios-GJ436b',
#        'COMPARISON/Atmo-GJ436b',
#        'COMPARISON/Helios-GJ1214-pure',
#        'COMPARISON/Helios-GJ1214-iso',
#        'COMPARISON/Helios-GJ1214-fw',
#        'COMPARISON/Helios-GJ1214-bw',
##        'COMPARISON/Tsai-2021-K218b',
#        'COMPARISON/Scheucher-k218-H2He1x',
#        'COMPARISON/Scheucher-k218-H2He50x',
##        'COMPARISON/Scheucher-k218-H2HeOcean',
#        'COMPARISON/Jupiter-orig',

##        'gj1214-RC-time',
##        'gj1214-RC-jacob',
#        'rteq_GJ1214-Teff775-z0.13-rjacob50',
#        'gj1214-time',
#        'gj1214-jacob',
#        'gj1214-time-tinit1000',
#        'gj1214-heng',
#        'gj1214-toon',
#        'gj1214-jacob-Toon',

#        'k218b-RC-time',
#        'k218b-RC-jacob',
#        'k218b-time',
#        'k218b-time-50x',
#        'k218b-jacob',
#        'RC_debug',
#        'RC_debug-20',
#        'RC_debug-500',

#        'Jupiter-Renyu',
#        'jupiter-time',
#        'jupiter-jacob',
#        'jupiter-jacob-Toon',

#        'jupiter-guillot',

        'jupiter-valid1-Guillot',
        'jupiter-valid2-Texp8',
        'jupiter-valid3-lowTiso',
        'jupiter-valid4-highTiso',
        'jupiter-valid5-TPJupiter',
        'jupiter-valid6-Tiso2k',

        'jupiter-valid1time-Guillot',
        'jupiter-valid2time-Texp8',
        'jupiter-valid3time-lowTiso',
        'jupiter-valid4time-highTiso',
        'jupiter-valid5time-TPJupiter',
        'jupiter-valid6time-Tiso2k',
        
        'jupiter-guillot',
        'INPUTS/gen_e8',
        'INPUTS/Tiso110',
        'INPUTS/Tiso1000',
        'INPUTS/Jupiter',
        'INPUTS/Tiso2000',

        ]
simulabels = [
#        r'HELIOS',
#        r'ATMO',
#        r'T$_{eff}$ 712 K',
#        r'HELIOS: pure absorption',
#        r'HELIOS: isotropic scattering',
#        r'HELIOS: forward scattering',
#        r'HELIOS: backward scattering',
#        r'Rayleigh scattering',
#        r'1D-TERRA: metal 1x',
#        r'1D-TERRA: metal 50x',
#        r'HELIOS: metal 100x T$_{int}$=30k',
#        r'metal 1x',
#        r'metal 1x T$_{int}$=30K',
#        r'metal 100x T$_{int}$=30K',
#        r'Solar',
#        r'Solar; T$_{int}$=30K',
        ]

colorcycle = 6 #max index of colors vector to be used; all if -1
uselabels = 0 #Use simulabels as plotlabels, 0:automatic labels generation from simu[]
selection =[-1] #Selection of simulations to compare *** ALL=[-1] ***
#selection =[0,1] #Selection of simulations to compare *** ALL=[-1] ***

#-------- Reference profiles to plot --------------------------------------
#JPL Gattaca file structure:
printJUPITER = 0    #EPACRIS Tstart
#DLR WARG file structure:
printUSstd1976 = 0 #Overplot Temperature profile from US Standard atmosphere 1976
printMIPAS = 0 #Overplot Temperature profile from US Standard atmosphere 1976
printACEFTS = 0 #Overplot Temperature profile from US Standard atmosphere 1976
printMARSref = 0    #Overplot Mars Reference atmosphere from Haberle 2017
marsrefname = 'Data/Mars/Reference_atm_clear_Haberle17.dat'
printMARS = 0    #Overplot Mars data from Mariner 9
printVIRA = 0    #Overplot Venus International Reference Atmosphere
viraname = 'Data/Venus/vira/T_VIRA1.txt'
printmendonca = 0 #Overplot cloud-free Venus from Mendonca+ 2015
mendname = 'Data/Venus/vira/mendonca-data.csv'
printhaus = 0 #Overplot Heating/Cooling for cloud-free Venus from Haus+ 2015 provided by R. Haus
printhausold = 0 #Overplot Heating/Cooling for cloud-free Venus taken directly from Haus+ 2015

#-------- Selection of Parameters to plot ---------------------------------
#--- LATEX STYLE SUPPORT --------------------------------------------------
#Define Range of parameters (NZ)
species = [ 'Z [km]','P [bar]','T [K]']
#mselection =[11,7,8,9] #Specify what to plot. *** ALL=[-1] ***
mselection =[-1] #Specify what to plot. *** ALL=[-1] ***
yselection = 1 #Specify which species to plot on the y-axis
ysecondary = -1 #Specify a secondary y-axis on the right. *** NONE= -1 ***
#-------- Plots, Files and Path Definitions -------------------------------
plotcols_profiles = 1 #number of columns of plots
overplot = 0 #Overplot all in one plot?
width = 5 #Plotwidth in inches
height = 5 #Plotheight for each figure
spaceh = 0.3 #Vertical fractional space between figures
spacew = 0.3 #0.8 #Horizontal fractional space between figures
home = '/scratch_lg_edge/renyu-group/mscheucher/Epacris/Results/'  #home directory
path = home  #
temppath = 'NewTemperature.dat' #same for temp 
plotpath = path + '../plots/'  #where to plot outfile
######### END USER INTERFACE  #############################################



print('\n**********************************')
print('Start: easyT.py')
print('**********************************')

plt.close("all")
plt.ioff()


######## Additional Definitions ###################################
#rc('text',usetex=True)
rc('font',family='serif')
rc('text.latex', preamble=r'\usepackage{color}')
#plt.rcParams['text.usetex'] = True

avogadro = 6.022e+23

icolor=colorcycle
if(icolor<0.0):icolor=20
colors=['C'+str(x) for x in range(icolor)]
#colors=['dimgrey','black','darkgrey','silver','C2','C1','C2']
#colors=['C2','C0','C1']
lines = ['-','--',':','-.']
markers = ['', '-', '+', 'xxx', '\\', '*', 'o', 'O', '.', '/', '//']


######## Processing User Inputs ###################################
if (selection[0] == -1):  compared = simu
else: compared = [simu[l] for l in selection]
if (mselection[0] == -1): mselection = [l for l in range(2,len(species))]
if (len(mselection) <3): plotcols_profiles = len(mselection)    
params = [species[l] for l in mselection]

######## Tempering with directories, if needed ####################
def assure_path_exists(path):
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)
        print('\n**********************************')
        print('<<<created folder>>>: %s \n' %dir)
        print('**********************************\n')

def delete_folder(path, count):
    dir = os.path.dirname(path)
    print(dir)
    print(count)
    if (count == 2):
        print('deleting: %s' %dir)
        shutil.rmtree(dir)

assure_path_exists(plotpath)

############### Read all Data ############################################
def read_alldata():
    global lendata
    global data
    lendata = []
    data = []
    iall=0
    tp1 = datetime.now()
    for runs in compared[:]: 
        lendata.append(1)
        data.append(1)
        filename = path+runs+'/'+temppath
        print(('reading: ',runs))
        temp = pd.read_csv(filename,delim_whitespace=True, header=None, skiprows=0)
        lendata[iall] = len(temp.values[:,0])
        data[iall] = temp.values[:,0:3]
        print(("N_layers: ",lendata[iall]))
        iall=iall+1
    tp2 = datetime.now()
    print('h/km, P/bar, T/K:')
    print(data[0][0:2][:])
    print('read_alldata() took %s [h:m:s]\n' %(tp2-tp1))
# END #### Read all Data ##########################################


######## Plot all Profiles ########################################
def plot_profiles():
        t0 = datetime.now()
        fig = plt.figure() 
        ax = fig.add_subplot(111)
        fig.subplots_adjust(hspace=spaceh, wspace=spacew)
        if (ysecondary != -1):ax2 = ax.twinx()
        im=0
        for plots in params[:]:
            legs= []
            labels= []
            legs2= []
            labels2= []
            isel=im
            im=im+1
            i=0
            if (im > 1 and overplot == 0):
                if (ysecondary != -1):
                    n=len(fig.axes)
                    for j in range(0,n,2):
                        k=int(old_div((n),(2*plotcols_profiles)))
                        fig.axes[j].change_geometry(k+1, plotcols_profiles, old_div(j,2)+1)
                        fig.axes[j+1].change_geometry(k+1, plotcols_profiles, old_div(j,2)+1)
                    ax = fig.add_subplot(k+1, plotcols_profiles, old_div(n,2)+1)
                    ax2 = ax.twinx() 
                else:    
                    n=len(fig.axes)
                    for j in range(n):
                        k=int(old_div(n,plotcols_profiles))
                        fig.axes[j].change_geometry(k+1, plotcols_profiles, j+1)
                    ax = fig.add_subplot(k+1, plotcols_profiles, n+1)

            xmin = ymin = 1e20 #SM2021
            xmax = ymax = -1e20 #SM2021
            for runs in compared[:]: 
                print(runs,': ',plots)
                tc1 = datetime.now()
                xdata = []
                ydata = []
                y2data = []

#                print('data[:,%s]: ' %(plots), data[:,mselection[isel]])
                for ii in range(0,lendata[i]):
                    xdata.append(1)
                    ydata.append(1)
                    y2data.append(1)
                    xdata[ii]=data[i][ii][mselection[isel]] 
                    if(yselection == 0): ydata[ii]=data[i][ii][yselection]
                    else: ydata[ii]=10**(data[i][ii][yselection])/10**5 
                    if (ysecondary != -1):
                        y2data[ii]=data[i][ii][ysecondary] 
                print('xdata[:]: ', xdata[:])
                print('ydata[:]: ', ydata[:])

                ax.set_title('Title',alpha=0.0) #for plot(bbox_inches('tight')) if legend is above figure
	        
                if (yselection == 1): 
                    ax.set_yscale('log')
                    ax.yaxis.set_major_locator(LogLocator(base=10,numticks=15))
                    ax.yaxis.set_minor_locator(LogLocator(base=10))
                ax.set_xscale('linear') 
                if (np.min(xdata)<xmin): xmin = np.min(xdata) #SM2021
                if (np.max(xdata)>xmax): xmax = np.max(xdata) #SM2021 
                if (np.min(ydata)<ymin): ymin = np.min(ydata) #SM2021
                if (np.max(ydata)>ymax): ymax = np.max(ydata) #SM2021
                ax.set_xlabel(r'%s' %plots)
                ax.set_ylabel(r'%s' %(species[yselection]))
	            
                tc2 = datetime.now()
                print('setup took %s [h:m:s]' %(tc2-tc1))
                mylinestyle = lines[int(math.floor(old_div(i,len(colors))))]
                mycolscheme = colors[i%len(colors)]
                if(printmendonca): mycolscheme = colors[(i+1)%len(colors)]
                if(printhaus): mycolscheme = colors[(i+1)%len(colors)]
                if (uselabels == 1): mylabel = ['%s'%(simulabels[i]),'']
                else: mylabel = ['%s'%(runs),'']
                if (overplot == 1): 
                    mycolscheme = colors[isel%len(colors)]
                    mylabel = ['%s'%(plots)]
	
                ax.plot(xdata,ydata,label=mylabel[0],linestyle=mylinestyle,color=mycolscheme,linewidth=1)
                if (ysecondary != -1):
                    ax2.plot(y2data,y2data,alpha=0)
                    ax2.set_ylabel(r'%s' %(species[ysecondary]))
                    if (ysecondary == 0): ax2.set_yscale('log')
                    ax2.set_ylim([y2data[-1],y2data[0]])
#	        if (mselection[isel]>7 and mselection[isel]<10): ax.set_xlim(0,50)
	        #ax.legend(loc='upper right',bbox_to_anchor=(1.70,1.02),ncol=1, fancybox=True, frameon=False, fontsize=6)
                ax.legend(loc='best', fancybox=True, frameon=False, fontsize=9)
	        #print('\n')
	        #break    
                i=i+1

            if (ax.set_xscale == "linear"): ax.set_xlim([xmin-10,xmax+10]) #SM2021 
            else: ax.set_xlim([0.9*xmin,1.1*xmax]) #SM2021
            ax.set_ylim([0.9*ymin,1.1*ymax]) #SM2021
            ax.set_ylim([1.1*ymax,0.9*ymin]) #SM2021

            if(printVIRA):
                if(ysecondary==10):ax2.set_ylim([0,124])
                if(yselection==10):ax.set_ylim([0,124])

            if(printhaus):
                if (mselection[isel]>7 and mselection[isel]<10):
#               if(ysecondary==10):ax2.set_ylim([0,124])
#               if(yselection==10):ax.set_ylim([20,90])
                    ax.set_xlim([0,15])
                    ax.set_ylim([10**1,5*10**(-4)])

	#sys.exit()
        t1 = datetime.now()
        print('Creation took %s [h:m:s]' %(t1-t0))
	
        tpl0 = datetime.now()
	
        if (overplot == 0):
            k=int(old_div((im-1.),plotcols_profiles) +1.)
            fig.set_size_inches(width*plotcols_profiles + spacew*width*(plotcols_profiles-1), height*k + spaceh*height*(k-1)) #whatever looks good and keeps the ratio for different figure numbers
        else:
            fig.set_size_inches(width, height)
        extension = '_plot.pdf' 
        plotfile = plotpath + outfile + extension
        plt.grid(True, linestyle=':')
        plt.savefig(plotfile, bbox_inches='tight')
        plt.close(fig)
        tpl1 = datetime.now()
        print('Plotting took %s [h:m:s]' %(tpl1-tpl0))
        print('\n**********************************')
        print('output: %s' %plotfile)
        print('**********************************')
# END ##### Plot all Profiles #####################################




###################################################################
########                ###########################################
########  MAIN PROGRAM  ###########################################
########                ###########################################
###################################################################

read_alldata()

if (printUSstd1976): USstd76()
if (printMIPAS): MIPAS()
if (printACEFTS): ACEFTS()
if (printMARSref): MARSref()
if (printMARS): MARS()
if (printVIRA): VIRA()
if (printmendonca): Mendonca()
if (printhaus): Haus()
if (printhausold): Hausold()

plot_profiles()

###################################################################
########                ###########################################
########  END PROGRAM   ###########################################
########                ###########################################
###################################################################
