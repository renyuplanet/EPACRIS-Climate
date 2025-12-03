######### 1DTERRA:climata-diag-comp.py ##########################################
#                                                                               #
# Routine to compare and analyze climate relevant data:                         #
#       * Heating / Cooling rates against literature                            #
#       * Temperature Profiles against literature and standard-atmospheres      #
#       * Radiative up/dn fluxes                                                #
#                                                                               #
# vr2.0 written by M.Scheucher 05/2020                                          #
#                                                                               #
#        -->  Exp: [sun,adleo, ..] -> Output/Diagnostics -> file.dat            #
# HOME  |                                                                       #
#        -->  Exp./plots ->  outfile.pdf                                        #
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
#       'mars-quiescent',
#       'mars-quiescent-updated',
        'Mars-Ctl-2021',
#       'mars-iup',
#       'mars-iup-ic',
#       'mars-iup-frel-ic-t300',
#       'mars-iup-frel2-ic-t300',
#       'mars-iup-frel2-t300',
       'mars-t300-evo0-01',
       'mars-t300-evo0-02',
       'mars-t300-evo0-11',
       'mars-t300-evo0-12',
       'mars-t300-evo0-21',
       'mars-t300-evo0-22',
       'mars-t300-evo0-31',
       'mars-t300-evo0-32',
       'mars-t300-evo1-01',
       'mars-t300-evo1-02',
       'mars-t300-evo1-11',
       'mars-t300-evo1-12',
       'mars-t300-evo1-21',
       'mars-t300-evo1-22',
       'mars-t300-evo1-31',
       'mars-t300-evo1-32',
        ]
simulabels = [
        r'1D-TERRA Mars control',
        r'CO$_2$ = 0.2 bars, init',
#        r'GJ176: CO$_2$ = 0.1 bars, init',
        r'2$^{nd}$ iteration',
        r'3$^{rd}$ iteration',
        r'4$^{th}$ iteration',
        r'5$^{th}$ iteration',
#        r'CO$_2$ = 1e-2, TRAPPIST',
#        r'CO$_2$ = 1e-1, TRAPPIST',
#        r'CO$_2$ = 2e-1, TRAPPIST',
#        r'CO$_2$ = 0.0, cyc 1',
#        r'CO$_2$ = 1e-2, cyc 1',
#        r'CO$_2$ = 1e-1, cyc 1',
#        r'CO$_2$ = 2e-1, cyc 1',
#        r'CO$_2$ = 1e-2, GJ176',
#        r'CO$_2$ = 1e-1, GJ176',
#        r'CO$_2$ = 2e-1, GJ176',
        ]
colorcycle = 7 #max index of colors vector to be used; all if -1
uselabels = 0 #Use simulabels as plotlabels, 0:automatic labels generation from simu[]
selection =[-1] #Selection of simulations to compare *** ALL=[-1] ***
#selection =[55,20] #Selection of simulations to compare *** ALL=[-1] ***
printUSstd1976 = 0 #Overplot Temperature profile from US Standard atmosphere 1976
printMIPAS = 0 #Overplot Temperature profile from US Standard atmosphere 1976
printACEFTS = 0 #Overplot Temperature profile from US Standard atmosphere 1976
printMARSref = 1    #Overplot Mars Reference atmosphere from Haberle 2017
marsrefname = 'Data/Mars/Reference_atm_clear_Haberle17.dat'
printMARS = 1    #Overplot Mars data from Mariner 9
printMCS_north = 1    #Overplot Mars data from Mars Climate Sounder 29.3 North
printMCS_south = 1    #Overplot Mars data from Mars Climate Sounder 59.7 South
mcsrefname = 'Data/Mars/MCS-data.txt'
printVIRA = 0    #Overplot Venus International Reference Atmosphere
viraname = 'Data/Venus/vira/T_VIRA1.txt'
printmendonca = 0 #Overplot cloud-free Venus from Mendonca+ 2015
mendname = 'Data/Venus/vira/mendonca-data.csv'
printhaus = 0 #Overplot Heating/Cooling for cloud-free Venus from Haus+ 2015 provided by R. Haus
printhausold = 0 #Overplot Heating/Cooling for cloud-free Venus taken directly from Haus+ 2015
outfile = 'Mars-Ctl-test' #Name of outfile
outfile = 'mars-test' #Name of outfile


#-------- Selection of Parameters to plot ---------------------------------
#--- LATEX STYLE SUPPORT --------------------------------------------------
#Define Range of parameters (NZ)
species = [ 'P/hPa','Fs-dn','Fs-up','Fth-dn','Fth-up','Fs-net','Fth-net','Ftot','Rheat','Rcool','Z [km]','T [K]']
#mselection =[11,7,8,9] #Specify what to plot. *** ALL=[-1] ***
mselection =[-1] #Specify what to plot. *** ALL=[-1] ***
mselection =[7,11] #Specify what to plot. *** ALL=[-1] ***
#mselection =[11] #Specify what to plot. *** ALL=[-1] ***
yselection = 0 #Specify what from species to plot on the y-axis
ysecondary = -1 #Specify a secondary y-axis on the right. *** NONE= -1 ***
#species = [ 'P/bar','Fs-dn','Fs-up','Fth-dn','Fth-up','Ftot','Fs-net','Fth-net','Alt [km]','T [K]' ]    
#mselection =[9,1,2,6,5,3,4,7] #Specify what to plot. *** ALL=[-1] ***
#yselection = 0 #Specify what from species to plot on the y-axis
#ysecondary = 8 #Specify a secondary y-axis on the right. *** NONE= -1 ***
#-------- Plots, Files and Path Definitions -------------------------------
plotcols_profiles = 3 #number of columns of plots
overplot = 0 #Overplot all in one plot?
width = 5 #Plotwidth in inches
height = 5 #Plotheight for each figure
spaceh = 0.3 #Vertical fractional space between figures
spacew = 0.3 #0.8 #Horizontal fractional space between figures
home = '/raid/user/m.scheucher/1D-TERRA/Experiments/'  #home directory
path = home  #
fluxpath = 'Output/Diagnostic/fluxes.dat'  #path within each sim directory
temppath = 'Output/Clima/Tprofile.dat' #same for temp 
plotpath = path + 'plots/'  #where to plot outfile
######### END USER INTERFACE  #############################################



print '\n**********************************'
print 'Start: climate-diag-comp.py'
print '**********************************'

plt.close("all")
plt.ioff()


######## Additional Definitions ###################################
rc('text',usetex=True)
rc('font',family='serif')
rc('text.latex', preamble='\usepackage{color}')

avogadro = 6.022e+23

icolor=colorcycle
if(icolor<0.0):icolor=20
colors=['C'+str(x) for x in range(icolor)]
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

######## US Standard Atmosphere 1976 ##############################
def USstd76():
    #Temperature parameterization after Report NASA-TM-X-74335,
    #U.S. STANDARD ATMOSPHERE 1976
    print 'Creating US Standard atmosphere 1976'
    global Tstd
    global Hstd
    npts = 110
    Tstd = np.zeros(npts)
    Hdata = 110 #data[2][0][6] 
    Hstd = np.linspace(0,Hdata,num=npts,endpoint=True)
    #Parameterization < 86km:T = Tb + Lb * ( H - Hb ); b from [0,7]
    r0 = 6356.776 #Earth radius for H->Z
    Hb = [0,11,20,32,47,51,71,84.852]
    Hb = [(r0*Hb[i])/(r0-Hb[i]) for i in range(0,len(Hb))] #H->Z
    Lb = [-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0]
    Tb = np.zeros(len(Hb))
    Tb[0] = 288.15
    for i in range(1,len(Hb)):
        Tb[i] = Tb[i-1]+Lb[i-1]*(Hb[i]-Hb[i-1]) 
    #Temp 86-91km:
    Tt = 186.8673
    #Parameterization 91-110km
    Tc = 263.1905
    AA = -76.3232
    a  = -19.9429

    #NOW let's calculate the Temperatures:
    for i in range(0,len(Hstd)):
        for j in range (1,len(Hb)):
            if (Hstd[i]<Hb[j]):
                Tstd[i] = Tb[j-1]+Lb[j-1]*(Hstd[i]-Hb[j-1])
                break
        if (Hstd[i]>=Hb[-1] and Hstd[i]<91.0):
            Tstd[i] = Tt
        if (Hstd[i]>=91.0 and Hstd[i]<=110.0):
            Tstd[i] = Tc + AA * math.sqrt(1 - (((Hstd[i]-91)/a)**2))

#   print 'Height:',Hdata
#   print 'Hb:',Hb[-1]
#   print 'Tb:',Tb
#   print 'Temp - US Standard: ',Tstd
#   print 'H - US Standard: ',Hstd
# END #### US Standard Atmosphere 1976 ############################

######## MIPAS Data ###############################################
def MIPAS():
    print 'Creating Earth MIPAS Data'
    global Tmipas
    global T1mipas
    global T2mipas
    global Hmipas
    global Pmipas

    mipasname='Data/Earth/Mipas.dat'

    mipas = pd.read_csv(mipasname,delim_whitespace=True,header=None,skiprows=1,names=['h','p','t','low','high','noise'])
    Tmipas=mipas.values[:,2]
    T1mipas=mipas.values[:,3]
    T2mipas=mipas.values[:,4]
    Hmipas=mipas.values[:,0]
    Pmipas=mipas.values[:,1]*1e-3

######## ACE-FTS Data #############################################
def ACEFTS():
    print 'Creating Earth ACE-FTS Data'
    global Tacefts
    global T1acefts
    global T2acefts
    global Hacefts
    global Pacefts

    aceftsname='Data/Earth/ACE-FTS.dat'

    acefts = pd.read_csv(aceftsname,delim_whitespace=True,header=None,skiprows=1,names=['h','p','t','low','high','noise'])
    Tacefts=acefts.values[:,2]
    T1acefts=acefts.values[:,3]
    T2acefts=acefts.values[:,4]
    Hacefts=acefts.values[:,0]
    Pacefts=acefts.values[:,1]*1e-0


# END #### US Standard Atmosphere 1976 ############################

######## Mars Reference Atmosphere ################################
def MARSref():
    #Haberle 2017 from MCS data
    print 'Creating Mars Reference Atmosphere'
    global Tmarsref
    global Hmarsref
    global Pmarsref

    haberle = pd.read_csv(marsrefname,delim_whitespace=True,header=None,skiprows=7,names=['z','T','p','rho'])
    Tmarsref=haberle.values[:,1]
    Hmarsref=haberle.values[:,0]
    Pmarsref=haberle.values[:,2]*1e-2 #Pa to mbar
    Rmarsref=haberle.values[:,3]*1e-2 #Pa to mbar

#   print 'Pvira ',Pvira
#   print 'Hvira ',Hvira
#   print 'Tvira ',Tvira[:,1]

########## Mars Data ##############################################

def MARS():
    #Mariner 9 IRIS data
    global Tmarsmax
    global Tmarsmin
    global Tmarsmean
    global Pmars

    Tmarsmin = np.array([150.0,157.0,168.0,160.0,160.0])
    Tmarsmax = np.array([170.0,174.0,180.0,193.0,208.0])
    Tmarsmean = (Tmarsmax+Tmarsmin)*0.5
    Pmars = [0.3e-0,0.5e-0,1.e-0,2e-0,5e-0]

######## Mars Climate Sounder Data ################################
def MCS_North():
    #Haberle 2017 from MCS data
    print 'Creating MCS North Reference'
    global Tmcs_n
    global Tlmcs_n
    global Thmcs_n
    global Pmcs_n
    global Smcs_n

    mcs_n = pd.read_csv(mcsrefname,delim_whitespace=False,header=None,skiprows=40,nrows=105,usecols=[1,2,3],dtype=float)
    Tmcs_n=mcs_n.values[:,1]
    Pmcs_n=mcs_n.values[:,0]*1e-2 #Pa to mbar
    Smcs_n=mcs_n.values[:,2] #sigma range
    Tlmcs_n=mcs_n.values[:,1]-abs(mcs_n.values[:,2])
    Thmcs_n=mcs_n.values[:,1]+abs(mcs_n.values[:,2])

#   print 'P MCS_North ',Pmcs_n
#   print 'T MCS_North ',Tmcs_n
#   print 'Sigma MCS_North ',Smcs_n
#   print 'T low MCS_North ',Tlmcs_n[10:-10]
#   print 'T high MCS_North ',Thmcs_n[10:-10]

def MCS_South():
    #Haberle 2017 from MCS data
    print 'Creating MCS North Reference'
    global Tmcs_s
    global Tlmcs_s
    global Thmcs_s
    global Pmcs_s
    global Smcs_s

    mcs_s = pd.read_csv(mcsrefname,delim_whitespace=False,header=None,skiprows=25798,nrows=105,usecols=[1,2,3],dtype=float)
    Tmcs_s=mcs_s.values[:,1]
    Pmcs_s=mcs_s.values[:,0]*1e-2 #Pa to bar
    Smcs_s=mcs_s.values[:,2] #sigma range
    Tlmcs_s=mcs_s.values[:,1]-abs(mcs_s.values[:,2])
    Thmcs_s=mcs_s.values[:,1]+abs(mcs_s.values[:,2])

    print 'P MCS_South ',Pmcs_s
    print 'T MCS_South ',Tmcs_s
    print 'Sigma MCS_South ',Smcs_s
    print 'T low MCS_South ',Tlmcs_s[0:-10]
    print 'T high MCS_South ',Thmcs_s[0:-10]

# END #### Mars ###################################################

######## Venus International Reference Atmosphere (VIRA) ##########
def VIRA():
    #Temperature parameterization after Report NASA-TM-X-74335,
    print 'Creating Venus Reference Atmosphere (VIRA)'
    global Tvira
    global Hvira
    global Pvira

    vira = pd.read_csv(viraname,delim_whitespace=True,header=None,skiprows=2,names=['p','z','T30','T45','T60','T75','T85'])
    Tvira=vira.values[:,2:6]
    Hvira=vira.values[:,1]
    Pvira=vira.values[:,0]*1e2

#   print 'Pvira ',Pvira
#   print 'Hvira ',Hvira
#   print 'Tvira ',Tvira[:,1]

# END #### Venus International Reference Atmosphere (VIRA) ########

######### Cloud Free Venus reference from Mendonca+ 2015 ###########
def Mendonca():
    print 'Creating Mendonca+ 2015 Temp profile'
    global Tmend
    global Pmend

    mend = pd.read_csv(mendname,delim_whitespace=False,header=None,skiprows=1,names=['T','p'])
    Tmend=mend.values[:,0]
    Pmend=mend.values[:,1]*1e2

#   print 'Pmend ',Pmend
#   print 'Tmend ',Tmend

# END ####  

######### Gas-only Venus Heating/Cooling from Haus+ 2015 ###########
def Haus(): 
    #Data provided by Rainer Haus in 2019
    print 'Creating Haus+ 2015 Heating/Cooling rates'
    global h_Haus 
    global p_Haus 
    global HEAT_Haus
    global COOL_Haus

    pzname = 'Data/Venus/vira/Fig24_Pressure-Altitude.dat' #New, provided by Rainer Haus
    heatname = 'Data/Venus/vira/Fig24_Heating.dat' #New, provided by Rainer Haus
    coolname = 'Data/Venus/vira/Fig24_Cooling.dat' #New, provided by Rainer Haus

    pz = pd.read_csv(pzname,delim_whitespace=True,header=None,skiprows=1,names=['p','h'])
    p_Haus=pz.values[:,0]*1e-0
    h_Haus=pz.values[:,1]
    heat = pd.read_csv(heatname,delim_whitespace=True,header=None,skiprows=0,names=['heat','h'])
    HEAT_Haus=heat.values[:,0]
    cool = pd.read_csv(coolname,delim_whitespace=True,header=None,skiprows=0,names=['cool','h'])
    COOL_Haus=cool.values[:,0]*(-1.)

# END ####  
def Hausold(): 
    #Data extracted from Paper Haus+ 2015
    print 'Creating old Haus+ 2015 Heating/Cooling rates'
    global hHhaus 
    global hChaus 
    global pHhaus 
    global pChaus 
    global HEAThaus
    global COOLhaus

    heatname = 'Data/Venus/vira/Haus15-heating.csv' #Old taken from Paper
    coolname = 'Data/Venus/vira/Haus15-cooling.csv' #Old taken from Paper

    heat = pd.read_csv(heatname,delim_whitespace=False,header=None,skiprows=1,names=['heat','h','p'])
    HEAThaus=heat.values[:,0]
    hHhaus=heat.values[:,1]
    pHhaus=heat.values[:,2]*1e2
    cool = pd.read_csv(coolname,delim_whitespace=False,header=None,skiprows=1,names=['cool','h','p'])
    COOLhaus=cool.values[:,0]
    hChaus=cool.values[:,1]
    pChaus=cool.values[:,2]*1e2

# END ####  


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
        filename = path+runs+'/'+fluxpath
        tempfname = path+runs+'/'+temppath
        print 'reading: ',runs
        fluxes = pd.read_csv(filename,delim_whitespace=True, header=None, skiprows=3,names=['P','FSdn','FSup','FIRdn','FIRup','FSnet','FIRnet','Ftot','Rheat','Rcool'])
        lendata[iall] = len(fluxes.values[:,0])
        #print 'Fluxes shape:',fluxes.shape
        fluxes.values[:,1:8] *=1e-3
        fluxes.values[:,9] = abs(fluxes.values[:,9])
        temp = pd.read_csv(tempfname,delim_whitespace=True, header=None, skiprows=1)
#       print 'Temp shape:',temp.shape
        data[iall] = np.hstack((fluxes.values[:,:],temp.values[:,1:3]))
        print "N_layers: ",lendata[iall]
        iall=iall+1
    tp2 = datetime.now()
#    print 'P/bar, FSDN, FSUP, FIRDN, FIRUP, FTOT, Alt/km, T/K:'
#    print data[0][0][:]
    print 'read_alldata() took %s [h:m:s]\n' %(tp2-tp1)
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
                        k=int((n)/(2*plotcols_profiles))
	                fig.axes[j].change_geometry(k+1, plotcols_profiles, j/2+1)
	                fig.axes[j+1].change_geometry(k+1, plotcols_profiles, j/2+1)
	            ax = fig.add_subplot(k+1, plotcols_profiles, n/2+1)
	            ax2 = ax.twinx() 
	        else:    
	            n=len(fig.axes)
	            for j in range(n):
                        k=int(n/plotcols_profiles)
	                fig.axes[j].change_geometry(k+1, plotcols_profiles, j+1)
	            ax = fig.add_subplot(k+1, plotcols_profiles, n+1)
            #Plot US Standard atmosphere 1976:
            if (printUSstd1976):            
                if (plots==species[11]):
	            ax.set_xlim([np.min(Tstd)-30,np.max(Tstd)+50])
                    print 'Plotting US Standard atmosphere 1976'
                    if (yselection==10):
                        ax.plot(Tstd,Hstd,label='US Standard 1976',linestyle=lines[1],color='blue',linewidth=1)
                    if (ysecondary==10):
                        ax2.plot(Tstd,Hstd,label='US Standard 1976',linestyle=lines[1],color='blue',linewidth=1)

            #Plot Mipas Data:
            if (printMIPAS):            
                if (plots==species[11]):
                    print 'Plotting MIPAS Data'
                    if (yselection==10):
                        ax.plot(Tmipas,Hmipas,color='black',linestyle=lines[1],lw=1)
                        ax.fill_betweenx(Hmipas,T1mipas,T2mipas,label='MIPAS range',hatch='\ ',color='black',linewidth=0,alpha=0.3)
                    if (yselection==0):
                        ax.plot(Tmipas,Pmipas,color='black',linestyle=lines[1],lw=1)
                        ax.fill_betweenx(Pmipas,T1mipas,T2mipas,label='MIPAS range',hatch='\ ',color='black',linewidth=0,alpha=0.3)

            #Plot ACEFTSData:
            if (printACEFTS):            
                if (plots==species[11]):
                    print 'Plotting ACEFTS Data'
                    if (yselection==10):
                        ax.plot(Tacefts,Hacefts,color='grey',linestyle=lines[2],lw=1)
                        ax.fill_betweenx(Hacefts,T1acefts,T2acefts,label='ACE-FTS range',hatch='/ ',color='grey',linewidth=0,alpha=0.3)
                    if (yselection==0):
                        ax.plot(Tacefts,Pacefts,color='grey',linestyle=lines[2],lw=1)
                        ax.fill_betweenx(Pacefts,T1acefts,T2acefts,label='ACE-FTS range',hatch='/ ',color='grey',linewidth=0,alpha=0.3)

            #Plot Mars Mariner 9 data:
            if (printMARS):            
                if (plots==species[11]):
                    print 'Plotting MARS Mariner 9 Data'
                    if (ysecondary==0):
                        #ax2.plot(Tmarsmean,Pmars,label='Mariner 9',color='green',linewidth=1,linestyle=lines[1])
                        ax2.fill_betweenx(Pmars,Tmarsmin,Tmarsmax,label='Mariner 9',color='green',linewidth=0,alpha=0.3)
                    if (yselection==0):
                        #ax.plot(Tmarsmean,Pmars,label='Mariner 9',color='green',linewidth=1,linestyle=lines[1])
                        #for nn in range(len(Pmars)):
                            #ax.hlines(Pmars[nn],Tmarsmin[nn],Tmarsmax[nn],color='green',linewidth=1)
                        ax.fill_betweenx(Pmars,Tmarsmin,Tmarsmax,label='Mariner 9',color='goldenrod',linewidth=0,alpha=0.5)

            #MCS Data:
            if (printMCS_north):
                if (plots==species[11]):
                    print 'Plotting MCS North Data'
                    if (ysecondary==0):
                        ax2.fill_betweenx(Pmcs_n[13:-12],Tlmcs_n[13:-12],Thmcs_n[13:-12],label='MCS 31 Oct 2011,Lat -59.7',color='goldenrod',linewidth=0,alpha=0.3)
                    if (yselection==0):
                        #ax.errorbar(Tmcs_n,Pmcs_n,xerr=Smcs_n,fmt=lines[1],linewidth=0.5,elinewidth=0.3,capsize=2,markeredgewidth=0.5)
                        ax.fill_betweenx(Pmcs_n[13:-10],Tlmcs_n[13:-10],Thmcs_n[13:-10],label='MCS 31 Oct 2011,Lat -59.7',color='dimgrey',linewidth=0,alpha=0.5)
            if (printMCS_south):
                if (plots==species[11]):
                    print 'Plotting MCS North Data'
                    if (ysecondary==0):
                        ax2.fill_betweenx(Pmcs_s[13:-14],Tlmcs_s[13:-14],Thmcs_s[13:-14],label='MCS 31 Oct 2011,Lat 29.3',color='black',linewidth=0,alpha=0.3)
                    if (yselection==0):
                        #ax.errorbar(Tmcs_n,Pmcs_n,xerr=Smcs_n,fmt=lines[1],linewidth=0.5,elinewidth=0.3,capsize=2,markeredgewidth=0.5)
                        ax.fill_betweenx(Pmcs_s[3:-14],Tlmcs_s[3:-14],Thmcs_s[3:-14],label='MCS 31 Oct 2011,Lat 29.3',color='black',linewidth=0,alpha=0.5)
            #Plot Mars Reference (Haberle2017):
            if (printMARSref):            
                if (plots==species[11]):
                    print 'Plotting Mars Reference Atmosphere (Haberle2017)'
                    if (yselection==10):
                        ax.plot(Tmarsref,Hmarsref,label='Haberle(2017) global',linestyle=lines[1],color='black',linewidth=1)
                    if (yselection==0):
                        #ax.plot(Tmarsref,Pmarsref,label='Haberle 2017 global',linestyle=lines[1],color='black',linewidth=1)
                        ax.fill_betweenx(Pmarsref,Tmarsref-1.0,Tmarsref+1.0,label='Haberle 2017 global',linewidth=0,color='black',alpha=0.7)


            #Plot Venus VIRA:
            if (printVIRA):            
                if (plots==species[11]):
                    print 'Plotting VENUS Reference Atmosphere (VIRA)'
                    if (yselection==10):
                        ax.plot(Tvira[:,0],Hvira,label='VIRA 30$^\circ$',linestyle=lines[1],color='blue',linewidth=1)
                    if (yselection==0):
                        ax.plot(Tvira[:,0],Pvira,label='VIRA 30$^\circ$',linestyle=lines[1],color='blue',linewidth=1)

            #Plot Mendonca+ 2015:
            if (printmendonca):            
                if (plots==species[11]):
                    print 'Plotting Mendonca+ 2015 Cloud-free Venus'
                    if (ysecondary==0):
                        ax2.plot(Tmend,Pmend,label='Mendonca+ 2015, cloud-free',linestyle=lines[2],color='black',linewidth=1)
                    if (yselection==0):
                        ax.plot(Tmend,Pmend,label='Mendonca+ 2015, cloud-free',linestyle=lines[2],color='black',linewidth=1)

            #Plot Haus+ 2015:
            if (printhaus):            
                if (plots==species[8]):
                    print 'Plotting Haus+ 2015 Heating rates'
                    if (yselection==0):
                        ax.plot(HEAT_Haus,p_Haus,label='Haus+ 2015, Gas-only',linestyle=lines[1],color='blue',linewidth=1)
                    if (yselection==10):
                        ax.plot(HEAT_Haus,h_Haus,label='Haus+ 2015, Gas-only',linestyle=lines[1],color='blue',linewidth=1)
                if (plots==species[9]):
                    print 'Plotting Haus+ 2015 Cooling rates'
                    if (yselection==0):
                        ax.plot(COOL_Haus,p_Haus,label='Haus+ 2015, Gas-only',linestyle=lines[1],color='blue',linewidth=1)
                    if (yselection==10):
                        ax.plot(COOL_Haus,h_Haus,label='Haus+ 2015, Gas-only',linestyle=lines[1],color='blue',linewidth=1)

            #Plot Haus+ 2015:
            if (printhausold):            
                if (plots==species[8]):
                    print 'Plotting old Haus+ 2015 Heating rates'
                    if (ysecondary==0):
                        ax2.plot(HEAThaus,pHhaus,label='old Haus+ 2015, Gas-only',linestyle=lines[2],color=colors[-1],linewidth=1)
                    if (yselection==0):
                        ax.plot(HEAThaus,pHhaus,label='old Haus+ 2015, Gas-only',linestyle=lines[2],color=colors[-1],linewidth=1)
                if (plots==species[9]):
                    print 'Plotting old Haus+ 2015 Cooling rates'
                    if (ysecondary==0):
                        ax2.plot(COOLhaus,pChaus,label='old Haus+ 2015, Gas-only',linestyle=lines[2],color=colors[-1],linewidth=1)
                    if (yselection==0):
                        ax.plot(COOLhaus,pChaus,label='old Haus+ 2015, Gas-only',linestyle=lines[2],color=colors[-1],linewidth=1)



            xmin = ymin = 1e20 #SM2021
            xmax = ymax = -1e20 #SM2021
	    for runs in compared[:]: 
	        print runs,': ',plots
	        tc1 = datetime.now()
                xdata = []
                ydata = []
                y2data = []

#                print 'data[:,%s]: ' %(plots), data[:,mselection[isel]]
                for ii in range(0,lendata[i]):
                    xdata.append(1)
                    ydata.append(1)
                    y2data.append(1)
                    xdata[ii]=data[i][ii][mselection[isel]] 
                    ydata[ii]=data[i][ii][yselection]*1e3 
	            if (ysecondary != -1):
                        y2data[ii]=data[i][ii][ysecondary] 

                ax.set_title('Title',alpha=0.0) #for plot(bbox_inches('tight')) if legend is above figure
	        
	        if (yselection == 0): ax.set_yscale('log')
                if (mselection[isel]==3 or mselection[isel]==8 or mselection[isel]==9): ax.set_xscale('log')
#               if (mselection[isel]==3): ax.set_xscale('log')
                else: ax.set_xscale('linear') 
#SM2021	        #ax.set_ylim([ydata[-1],ydata[0]])
                if (np.min(xdata)<xmin): xmin = np.min(xdata) #SM2021
                if (np.max(xdata[:])>xmax): xmax = np.max(xdata) #SM2021
                if (np.min(ydata)<ymin): ymin = np.min(ydata) #SM2021
                if (np.max(ydata)>ymax): ymax = np.max(ydata) #SM2021
#SM2021         #if (ax.set_xscale == "linear"): ax.set_xlim([np.min(xdata)-10,np.max(xdata)+10])
#SM2021         #else: ax.set_xlim([0.9*np.min(xdata),1.1*np.max(xdata)])
	        if (mselection[isel]>0 and mselection[isel]<8): 
                    ax.set_xlabel(r'%s [W/m$^2$]' %plots)
                elif (mselection[isel]>7 and mselection[isel]<10): 
                    ax.set_xlabel(r'%s [K/day]' %plots)
                else:
                    ax.set_xlabel(r'%s' %plots)
                ax.set_ylabel(r'%s' %(species[yselection]))
	            
	        tc2 = datetime.now()
	        print 'setup took %s [h:m:s]' %(tc2-tc1)
	    
                mylinestyle = lines[int(math.floor(i/len(colors)))]
	        mycolscheme = colors[i%len(colors)]
                if(printmendonca): mycolscheme = colors[(i+1)%len(colors)]
                if(printhaus): mycolscheme = colors[(i+1)%len(colors)]
                if (uselabels == 1): mylabel = ['%s'%(simulabels[i]),'']
                else: mylabel = ['%s'%(runs),'']

	        if (overplot == 1): 
	            mycolscheme = colors[isel%len(colors)]
	            mylabel = ['%s'%(plots)]
	
                ax.plot(xdata,ydata,label=mylabel[0],linestyle=mylinestyle,color=mycolscheme,linewidth=1)
                #ax.plot(xdata,ydata,label=mylabel[0],linestyle=lines[1],color='red',linewidth=2)
	        if (ysecondary != -1):
	            ax2.plot(y2data,y2data,alpha=0)
	            ax2.set_ylabel(r'%s' %(species[ysecondary]))
                    if (ysecondary == 0): ax2.set_yscale('log')
	            ax2.set_ylim([y2data[-1],y2data[0]])
	                            
#	        if (mselection[isel]>7 and mselection[isel]<10): ax.set_xlim(0,50)
	        #ax.legend(loc='upper right',bbox_to_anchor=(1.70,1.02),ncol=1, fancybox=True, frameon=False, fontsize=6)
                ax.legend(loc='best', fancybox=True, frameon=False, fontsize=9)
	        #print '\n'
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

if (printUSstd1976): USstd76()
if (printMIPAS): MIPAS()
if (printACEFTS): ACEFTS()
if (printMARSref): MARSref()
if (printMARS): MARS()
if (printMCS_north): MCS_North()
if (printMCS_south): MCS_South()
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
