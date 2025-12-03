######### 1DTERRA-plots.py ######################################################
#                                                                               #
# Routine to compare various data vectors from different files in separate plots#
#                                                                               #
# v1.0 written by M.Scheucher 11/2018                                           #
#                                                                               #
#        -->  sim: [sun,adleo, ..] -> Output/Diagnostics -> file.dat            #
# HOME  |                                                                       #
#        -->  plots ->  outfile.pdf                                             #
#                                                                               #
# - Simulation folders are stored in 'simu' in USER INTERFACE                   #
# - Available plot parameters are stored in 'species' in USER INTERFACE         #
# - 'mselection' specifies which parameters from 'species' to plot              #
# - SET column-plot Y-RANGE in L:~390                                           #
#                                                                               #
# To change the structure, carefully review and modify paths in USER INTERFACE. #
#                                                                               #
#################################################################################

import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import math
import matplotlib.ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re
from datetime import datetime
from scipy.stats import binned_statistic

#--- Overwriting Hatch-Linewidth in backend python.lib -------------------
import matplotlib.backends.backend_pdf
from matplotlib import six
from matplotlib.backends.backend_pdf import Name, Op
from matplotlib.transforms import Affine2D
from matplotlib.path import Path



######### USER INTERFACE  #################################################
#Define Range of simulation folders
simu = [
        'base-5-none',
        'rebin-quads-5-none',

        ]
simulabels = [
        r'SOLAR-MRAC',
        r'REDFOX',
        ]
uselabels = 0 #Use simulabels as plotlabels, 0:automatic labels generation from simu[]
selection =[-1] #Selection of simulations to compare *** ALL=[-1] ***
outfile = 'validate-FUP' #Name of outfile
#outfile = 'REDFOX-MRAC32-FUPDN-comparison' #Name of outfile
#outfile = 'Rayleigh-FUPDN-comparison' #Name of outfile
#-------- Selection of Parameters to plot ---------------------------------
#--- LATEX STYLE SUPPORT --------------------------------------------------
#Define Range of parameters (NZ)
mrac = 0 #old MRAC routine?
solar = 0 #old SOLAR routine?
solir = [38,25] #38 SOLAR and 25 MRAC bands
if(solar): solir[0]=38
if(mrac): solir[1]=25
rtselect = 1 #plot Solar(0) or IR(1) bands. *** BOTH=-1 ***
fupdnselect = [0] #plot Up(0) or Down(1) fluxes. *** BOTH=-1 ***
bandselect = [-1] #Specify which bands to plot. *** ALL=[-1] ***
#-------- Plots, Files and Path Definitions -------------------------------
plotcols_profiles = 5 #number of columns of plots
overplot = 0 #Overplot everything in one plot?
width = 5 #Plotwidth in inches
height = 5 #Plotheight for each figure
spaceh = 0.3 #Vertical fractional space between figures
spacew = 0.3 #0.8 #Horizontal fractional space between figures
home = '/raid/user/m.scheucher/c4/' #Home Directory
sim = 'models/1D-TERRA/Experiments/'  #Data path within home directory
path = home + sim #+ '/' #CHECK: Full path to Molecule dirs
solpath = 'Output/Diagnostic/SW-FUPDN.dat'  #path within each sim directory
if(solar): solpath = 'Output/Diagnostic/SOLAR-FUPDN.dat'  #path within each sim directory
irpath = 'Output/Diagnostic/LW-FUPDN.dat' #same for temp 
if(mrac): irpath = 'Output/Diagnostic/MRAC-FUPDN.dat' #same for temp 
plotpath = path + 'plots/'  #where to plot outfile
######### END USER INTERFACE  #############################################



print '\n**********************************'
print 'Start: C4-Chem.py'
print '**********************************'

plt.close("all")
plt.ioff()


######## Additional Definitions ###################################
rc('text',usetex=True)
rc('font',family='serif')
rc('text.latex', preamble='\usepackage{color}')

avogadro = 6.022e+23
#colors = [
#'lightgrey','dimgrey','silver','grey','darkgrey','black',
#'olive','yellowgreen','lawngreen','sage','darkgreen','lime',
#'darkorange','orange','darkgoldenrod','goldenrod','gold','yellow',
#'lightcoral','indianred','firebrick','darkred','red','tomato',
#'aqua','dodgerblue','deepskyblue','cornflowerblue','darkblue','blue',
#'indigo','darkviolet','mediumorchid','purple','magenta','deeppink',
#] #'white'

colors = [
'black','orange','darkgreen','darkgoldenrod','darkred',
'indigo','aqua','orange','purple','olive','darkblue',
] #'white'
lines = ['-','--',':','-.']
markers = ['', '-', '+', 'xxx', '\\', '*', 'o', 'O', '.', '/', '//']


######## Processing User Inputs ###################################
if (selection[0] == -1):  compared = simu
else: compared = [simu[l] for l in selection]
if (rtselect == -1): rtselect = sum([solir[l] for l in range(0,len(solir))])
else: rtselect = solir[rtselect]
if (fupdnselect[0] == -1): fupdnselect = [0,1]
if (bandselect[0] == -1): bandselect = [l for l in range(1,rtselect+1)]

params = len(fupdnselect) * (bandselect[-1])

print 'compared=',compared
print 'rtselect=',rtselect
print 'fupdnselect=',fupdnselect
print 'bandselect=',bandselect
print 'Total number of plots=',params


######## Tempering with directories, if needed ####################
def assure_path_exists(path):
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)
        print '\n**********************************'
        print ('<<<created folder>>>: %s \n' %dir)
        print '**********************************\n'

def delete_folder(path, count):
    dir = os.path.dirname(path)
    print (dir)
    print (count)
    if (count == 2):
        print ('deleting: %s' %dir)
        shutil.rmtree(dir)

assure_path_exists(plotpath)
#for filename in glob.glob(impath):
#    count = count +1
#    delete_folder(filename, count)
# END #### Various Functions ######################################

######## Read all Data ############################################
def read_alldata():
    global lendata
    global updata
    global dndata
    lendata = []
    updata = []
    dndata = []
    iall=0
    tp1 = datetime.now()
    for runs in compared[:]: 
        lendata.append(1)
        updata.append(1)
        dndata.append(1)
        solname = path+runs+'/'+solpath
        irname = path+runs+'/'+irpath
        print 'reading: ',runs
        solfluxes = pd.read_csv(solname,delim_whitespace=True, header=None, skiprows=1)
        lendata[iall] = len(solfluxes.values[::2,0]) 
        irfluxes = pd.read_csv(irname,delim_whitespace=True, header=None, skiprows=1)
        if (rtselect==solir[0]+solir[1]):
            updata[iall] = np.hstack((solfluxes.values[-2::-2,2::],irfluxes.values[0::2,2::]))
            dndata[iall] = np.hstack((solfluxes.values[-1::-2,2::],irfluxes.values[1::2,2::]))
#           if (mrac == 1):
#               updata[iall] = np.hstack((solfluxes.values[-2::-2,2::],irfluxes.values[0::2,-1:1:-1]))
#               dndata[iall] = np.hstack((solfluxes.values[-1::-2,2::],irfluxes.values[1::2,-1:1:-1]))
        if (rtselect==solir[0]):
            updata[iall] = solfluxes.values[-2::-2,2::]
            dndata[iall] = solfluxes.values[-1::-2,2::]
        if (rtselect==solir[1]):
            updata[iall] = irfluxes.values[0::2,2::]
            dndata[iall] = irfluxes.values[1::2,2::]
#           if (mrac == 1):
#               updata[iall] = irfluxes.values[0::2,-1:1:-1]
#               dndata[iall] = irfluxes.values[1::2,-1:1:-1]
        print "N_layers: ",lendata[iall]
        iall=iall+1
    tp2 = datetime.now()
#    print 'P/bar, FSDN, FSUP, FIRDN, FIRUP, FTOT, Alt/km, T/K:'
    print 'FUP (%s): %s'%(lendata[0],updata[0][0][:])
    print 'FDN (%s): %s'%(lendata[0],dndata[0][0][:])
    print 'read_alldata() took %s [h:m:s]\n' %(tp2-tp1)
# END #### Read all Data ##########################################


######## Plot all Profiles ########################################
def plot_profiles():
	t0 = datetime.now()
	fig = plt.figure() 
	ax = fig.add_subplot(111)
	fig.subplots_adjust(hspace=spaceh, wspace=spacew)
        
	im=0
	for plots in range(params):
            isel=im
	    im=im+1
            band=im
            if(len(fupdnselect)==2): band=int(isel/2)+1
	    i=0
	    if (im > 1 and overplot == 0):
	        n=len(fig.axes)
	        for j in range(n):
                    k=int(n/plotcols_profiles)
	            fig.axes[j].change_geometry(k+1, plotcols_profiles, j+1)
	        ax = fig.add_subplot(k+1, plotcols_profiles, n+1)
	
	    for runs in compared[:]: 
	        print runs,': ',plots
	        tc1 = datetime.now()
                xdata = []
                ydata=np.arange(lendata[i])
	        ax.set_yscale('linear')
	        ax.set_ylim([ydata[0],ydata[-1]])


                for ii in range(0,lendata[i]):
                    xdata.append(1)
                    if (len(fupdnselect)==1):
                        if (fupdnselect[0]==0): 
                            xdata[ii]=updata[i][ii][isel] 
	                    ax.set_xlabel(r'Fup [W/m$^2$]')
                        else: 
                            xdata[ii]=dndata[i][ii][isel]
	                    ax.set_xlabel(r'Fdn [W/m$^2$]')
                    else:
                        if (isel%2==0): 
                            xdata[ii]=updata[i][ii][int(isel/2)] 
	                    ax.set_xlabel(r'Fup [W/m$^2$]')
                        else: 
                            xdata[ii]=dndata[i][ii][int(isel/2)]
	                    ax.set_xlabel(r'Fdn [W/m$^2$]')

                if (rtselect==solir[0]): ax.set_title('SW Band %s'%(band)) 
                elif (rtselect==solir[1]): ax.set_title('LW Band %s'%(band))
                else: 
                    if (band<=solir[0]): ax.set_title('SW Band %s'%(band))
                    else: ax.set_title('LW Band %s'%(band-38))
                if (sum(xdata[:])==0.0): ax.set_xscale('linear')
                else: ax.set_xscale('log') 
                ax.set_ylabel(r'Layer' )
	            
	        tc2 = datetime.now()
	        print 'setup took %s [h:m:s]' %(tc2-tc1)
	    
                mylinestyle = lines[i%len(lines)]
	        mycolscheme = colors[i%len(colors)]
                if (uselabels == 1): mylabel = ['%s'%(simulabels[i]),'']
                else: mylabel = ['%s'%(runs),'']

	        if (overplot == 1): 
	            mycolscheme = colors[isel%len(colors)]
	            mylabel = ['%s'%(plots)]
	
                ax.plot(xdata,ydata,label=mylabel[0],linestyle=mylinestyle,color=mycolscheme,linewidth=1)
	        #ax.legend(loc='upper right',bbox_to_anchor=(1.70,1.02),ncol=1, fancybox=True, frameon=False, fontsize=6)
                ax.legend(loc='best', fancybox=True, frameon=False, fontsize=6)
	        #print '\n'
	        #break    
	        i=i+1

	#sys.exit()
	t1 = datetime.now()
	print 'Creation took %s [h:m:s]' %(t1-t0)
	
	tpl0 = datetime.now()
	
	if (overplot == 0):
            k=int((im-1.)/plotcols_profiles +1.)
	    fig.set_size_inches(width*plotcols_profiles + spacew*width*(plotcols_profiles-1), height*k + spaceh*height*(k-1)) #whatever looks good and keeps the ratio for different figure numbers
	else:
	    fig.set_size_inches(width, height)
	
        extension = '_plot.pdf' 
	plotfile = plotpath + outfile + extension
	plt.savefig(plotfile, bbox_inches='tight')
	plt.close(fig)
	
	tpl1 = datetime.now()
	print 'Plotting took %s [h:m:s]' %(tpl1-tpl0)
	
	print '\n**********************************'
	print 'output: %s' %plotfile
	print '**********************************'
# END ##### Plot all Profiles #####################################




###################################################################
########                ###########################################
########  MAIN PROGRAM  ###########################################
########                ###########################################
###################################################################

read_alldata()

plot_profiles()

###################################################################
########                ###########################################
########  END PROGRAM   ###########################################
########                ###########################################
###################################################################
