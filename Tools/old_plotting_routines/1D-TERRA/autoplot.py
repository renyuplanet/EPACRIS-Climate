######### autoplot.py ###########################################################
#                                                                               #
# Routine to plot 1D-TERRA atm profiles from 'composition_vmr/molcc/mmr.dat'    #
#                                                                               #
# v6.0 written by M.Scheucher 06/2020                                           #
# - Updated from v5.0 to work with BLACKWOLF                                    #
# - This is essentially a flexible version of 1DTERRA/Plotting/py4chem.py       #
# <email:scheucher@tu-berlin.de>  <markus.scheucher@icloud.com>                 #
#                                                                               #
# Assumed Directory structure follows the new 1D-TERRA -Code structure:         #
#                                                                               #
#        -->  Exp: [sun,adleo, ..] -> Output/Photo ->  composition_xxx.dat      #
# HOME  |                                                                       #
#        -->  output/plots ->  outfile.pdf, columnfile.pdf                      #
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
import matplotlib.ticker as ticker
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

######### General Purpose #################################################
# Is this the autoplot routine? (auto=True), Or used to compare runs? (auto=False)
auto=True

######### Path Definitions ################################################
path = '/raid/user/m.scheucher/1D-TERRA/Experiments/'  #Molecules Data path within home directory
plotpath = path + 'plots/'  #where to plot outfile
datapath = 'Output/Photo/composition_vmr.dat'  #path within each sim directory
molccpath = 'Output/Photo/composition_molcc.dat' #same for mol/ccm
mfracpath = 'Output/Photo/composition_mmr.dat' #same for mass fractions

if(auto): #current directory for autoplot:
    expdir = re.split('/',os.getcwd())
    expdir = expdir[-1]
    path = '../' 
    plotpath = 'Output/plots/'  #where to plot outfile
######### USER INTERFACE  #################################################
#if auto=False: Define Range of simulation folders
simu = [
        'earth',
        ]
simulabels = [

        'Earth-Ctl',
        ]
colorcycle = 8 #max index of colors vector to be used; all if -1
uselabels = 1 #Use simulabels as plotlabels, 0:automatic labels generation from simu[]
selection =[-1] #Selection of simulations to compare *** ALL=[-1] ***
outfile = 'compared'
if(auto): 
    outfile = expdir #Name of outfile
    uselabels = 0
    colorcycle = 1
#overplotting literature data (here Meadows 2018 Temp):
printMeadows = 0 
meadowsname = 'Data/CenB/Meadows-2018-Temp.csv'

#--------------------------------------------------------------------------
#-------- Selection of Parameters to plot ---------------------------------
#--------------------------------------------------------------------------
# Example adapted/resorted from a composition_vmr header.
# Refer to species in this spelling (List not complete - check network.f!!):
#--------------------------------------------------------------------------
#   Height  Layer P Weight Density ... atm. attributes
#   T ... Temperature plot
#   UVA UVB UVC FUV ... Radiation dosages
#--- Molecules ------------------------------------------------------------
#   C C2 CO CO2 CH3 CH4 CH 3CH2 C2H2O CS C2H C2H3 C2H2 C2H2N C2H4 C2H5 CH2CCH2 CH3C2H C3H5 C3H6 C3H7 C3H2 1CH2 C4H2 C3H8 CH3ONO CH3ONO2 CH3Cl CH2Cl CS2 C5H4 CH3CO CH3OOH CH3S CH4S C2H5CHO C2H4NH C3H3 CH3NH2 COCl2 COCl C4H C2H6 CH3O2 C2HO CH3OH CN CH3CHO C2H3O CNO C2H5O 
#   Cl ClO Cl2O ClONO2 Cl2O2 ClOO Cl2 ClONO ClS2 ClCO3 ClSO2 Cl2S2       
#   H H2 HO2 HSO2 HSO3 H2O2 HNO HNO2 HNO3 HCl HOCl HCS HCOO HCOOH HO2NO2 H2SO4 H3CO HCN HCNO H2CO HCO HSO H2O H2S HS 
#   N N2 NO NO2 NO3 N2O NH NH2 NH3 NOCl N2H3 N2H4 N2H2 N2O5            
#   O O1D O2 O3 OH OCS OClO OCS2 OSCl 
#   S S2 S3 S4 S5 S6 S7 S8 SO SO2 SO3 SCl SCl2 S2O S2O2 1SO2 3SO2 
#   SO4AER
#   Ar
#--- Groups of Molecules --------------------------------------------------
#   NOx HOx ClOx 
#--- Molecule Ratios ------------------------------------------------------
#   OH/HO2 NO/NO2 Cl/ClO Clo/ClO2 O/O1D
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

#Define Range of parameters (NZ)
#mselection =['T','H2O'] #Specify what to plot. *** ALL=[-1] ***
mselection =[-1] #Specify what to plot. *** ALL=[-1] ***
yselection = 'P' #Specify y-axis (Height, P, or Density)!ONLY 1 choice possible
ysecondary = '' #Specify a secondary y-axis on the right (see choices above). *** NONE= '' ***
#-------- Plots, Files and Path Definitions -------------------------------
vmrmolcc = 0  #0:plot vmr profiles, 1:plot molec/ccm profiles, 2:Mass Fractions, -1: plot all
plotcolumns = 1 #1:plot column densities. (update 2020: AUTO prints it to stdout, but doesn't plot) 
plotprofiles = 1 #1:Plot profiles
plotcols_profiles = 4 #number of columns of plots when plptting profiles
showmarkers = 0 #1:Show markers in profile plots
overplot = 0 #Overplot species concentrations on same plot?
width = 3.5 #Plotwidth in inches
width_colplot = 10*width #width for columns plot
height = 3.5 #Plotheight for each figure
height_colplot = 2*height #width for columns plot
spaceh = 0.30 #Vertical fractional space between figures
spacew = 0.30 #0.8 #Horizontal fractional space between figures
columnfile = outfile + '-columns' #Column Densities plot
######### END USER INTERFACE  #############################################



print '\n**********************************'
print 'Start: autoplot.py'
print '**********************************'

plt.close("all")
plt.ioff()

######## Additional Definitions ###################################
rc('text',usetex=True)
rc('font',family='serif')
rc('text.latex', preamble='\usepackage{color}')

avogadro = 6.022e+23 #needed for some comversions

icolor=colorcycle
if(icolor<0.0):icolor=20
colors=['C'+str(x) for x in range(icolor)]
lines = ['-','--',':','-.']
markers = ['', '-', '+', 'xxx', '\\', '*', 'o', 'O', '.', '/', '//']

######## Functions for Column Plot Settings #######################
# Not used at the moment (06/2020)!!!

#--- Defining ax.xaxis formatting style ---------------------------
class CustomTicker(ticker.LogFormatterSciNotation):
    def __call__(self, x, pos=None):
        if (x>0.01 and x<1000):
            return "{x:g}".format(x=x)
        else:
            return ticker.LogFormatterSciNotation.__call__(self,x,pos=None)

class CustomLinTicker(ticker.ScalarFormatter):
    def __call__(self, x, pos=None):
        if (x>0.01 and x<1000):
            return "{x:g}".format(x=x)
        else:
            return "{x:%.1e}".format(x=x)

#--- Overwriting Hatch-Linewidth in backend python.lib ------------
def setCustomHatchWidth(customHatchWidth):

    def _writeHatches(self):
        hatchDict = dict()
        sidelen = 72.0
        for hatch_style, name in six.iteritems(self.hatchPatterns):
            ob = self.reserveObject('hatch pattern')
            hatchDict[name] = ob
            res = {'Procsets':
                   [Name(x) for x in "PDF Text ImageB ImageC ImageI".split()]}
            self.beginStream(
                ob.id, None,
                {'Type': Name('Pattern'),
                 'PatternType': 1, 'PaintType': 1, 'TilingType': 1,
                 'BBox': [0, 0, sidelen, sidelen],
                 'XStep': sidelen, 'YStep': sidelen,
                 'Resources': res})

            # lst is a tuple of stroke color, fill color,
            # number of - lines, number of / lines,
            # number of | lines, number of \ lines
            rgb = hatch_style[0]
            self.output(rgb[0], rgb[1], rgb[2], Op.setrgb_stroke)
            if hatch_style[1] is not None:
                rgb = hatch_style[1]
                self.output(rgb[0], rgb[1], rgb[2], Op.setrgb_nonstroke,
                            0, 0, sidelen, sidelen, Op.rectangle,
                            Op.fill)

            self.output(customHatchWidth, Op.setlinewidth)

            # TODO: We could make this dpi-dependent, but that would be
            # an API change
            self.output(*self.pathOperations(
                Path.hatch(hatch_style[2]),
                Affine2D().scale(sidelen),
                simplify=False))
            self.output(Op.stroke)

            self.endStream()
        self.writeObject(self.hatchObject, hatchDict)

    matplotlib.backends.backend_pdf.PdfFile.writeHatches = _writeHatches
# END #### Functions for Column Plot settings #####################

####  Tempering with directories, if needed #######################
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
# END #### Directory edits ########################################

######## Read Meadows 2018 Temp profile ###########################
def Meadows():
    print 'Creating Meadows+ 2018 Temp profile'
    global Tmeadows
    global Pmeadows

    mead = pd.read_csv(meadowsname,delim_whitespace=False,header=None,skiprows=1,names=['T','p'])
    Tmeadows = mead.values[:,0]
    Pmeadows = mead.values[:,1]/100.
# END #### Meadows 2018 Temp profile ##############################

######## Read all Data ############################################
def read_alldata():
    global lendata
    global alldata
    global molccdata
    global mfracdata
    global titles
    global units
    lendata = []
    alldata = []
    molccdata = []
    mfracdata = []
    iall=0
    tp1 = datetime.now()
    for runs in compared[:]: 
        lendata.append(1)
        alldata.append(1)
        molccdata.append(1)
        mfracdata.append(1)
        filename = path+runs+'/'+datapath
        molccfname = path+runs+'/'+molccpath
        mfracfname = path+runs+'/'+mfracpath
        print 'reading: ',runs
        pandahead = pd.read_csv(filename,delim_whitespace=True,nrows=1)
        titles = pandahead.columns.values[:]
        pandaunits = pd.read_csv(filename,delim_whitespace=True,skiprows=1,nrows=1)
        units = pandaunits.columns.values[:]
        panda = pd.read_csv(filename,delim_whitespace=True, header=None, skiprows=2, names=titles)
        lendata[iall] = len(panda.values[:,0]) 
        alldata[iall] = panda.values[:,:] 
        molcc = pd.read_csv(molccfname,delim_whitespace=True, header=None, skiprows=2, names=titles)
        molccdata[iall] = molcc.values[:,:] 
        for imod in range(0,lendata[iall]):
            molccdata[iall][imod][list(titles).index('Density')] *= 1e+3 / avogadro * molccdata[iall][imod][list(titles).index('Weight')]
        mfrac = pd.read_csv(mfracfname,delim_whitespace=True, header=None, skiprows=2, names=titles)
        mfracdata[iall] = mfrac.values[:,:] 
        for imod in range(0,lendata[iall]):
            mfracdata[iall][imod][list(titles).index('Density')] *= 1e+3
        print "N_layers: ",lendata[iall]
        iall=iall+1
    tp2 = datetime.now()
    print 'read_alldata() took %s [h:m:s]\n' %(tp2-tp1)
# END #### Read all Data ##########################################


######## Plot Columns #############################################
# Needs some edits to work again with 1D-TERRA (06/2020)
def plot_columns():
    tpl0 = datetime.now()
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    n=0
    for l in range(len(molecules)): 
        if (list(titles).index(molecules[l]) > list(titles).index('T') and list(titles).index(molecules[l]) < list(titles).index('Weight')):n=n+1
    ind = np.arange(n)
    barwidth = 1./(len(compared)+1)
    legs= []
    legs2= []
    labels= []
    labels2= []
    i=0
    for runs in compared[:]: 
        data=alldata[i]
        column = []
        xlabels = []
        #For conformity for calling special_profiles()
        mylinestyle = lines[0]
        mycolscheme = colors[i%len(colors)]
        dz = np.zeros(lendata[i])
        mymarker = ''
        if (uselabels == 1): mylabel = ['%s'%(simulabels[i]),'']
        else: mylabel = ['%s'%(runs),'']
        #END
        isel=0
        icol=0
	print '\n**********************************'
        print runs
	for plots in molecules[:]:
            todo=0

            dz[0] = data[0,list(titles).index('Height')]
            for n in range(1,lendata[i]):
                dz[n] = data[n,list(titles).index('Height')] - data[n-1,list(titles).index('Height')] - dz[n-1]
            dz[:] *= 2.0

            conversion = data[:,list(titles).index('Density')]* avogadro/data[:,list(titles).index('Weight')]*1.e-3 * (dz[:]) * 10**5 / 2.7 * 10**(-16)
            if (list(titles).index(plots) > list(titles).index('T') and list(titles).index(plots) < list(titles).index('Weight')):
                column.append(sum(data[:,list(titles).index(plots)] * conversion))
                todo=1
            if (plots == 'NOx'):
                column.append(sum((data[:,list(titles).index('NO')]+data[:,list(titles).index('NO2')]+data[:,list(titles).index('N')]+data[:,list(titles).index('NO3')]) * conversion))
                todo=1
            if (plots == 'HOx'):
                column.append(sum((data[:,list(titles).index('OH')]+data[:,list(titles).index('HO2')]+data[:,list(titles).index('H')]) * conversion))
                todo=1
            if (plots == 'ClOx'):
                column.append(sum((data[:,list(titles).index('Cl')]+data[:,list(titles).index('ClO')]+data[:,list(titles).index('OClO')]) * conversion))
                todo=1
            if (plots == 'UVA' or plots == 'UVB' or plots == 'UVC' or plots == 'FUV'): todo=2    
            if (todo==1):
                print plots,' Column:   ',column[icol],'[DU]'
                xlabels.append(plots)
                icol=icol+1
            if (todo == 2): 
                surf_dosage=(data[0,list(titles).index(plots)]-0.5*(data[1,list(titles).index(plots)]-data[0,list(titles).index(plots)]))
                if (surf_dosage < 0.0): surf_dosage = 0.0
                print plots,'dosage:   ',surf_dosage,'[W/m$^2$]'    
            isel=isel+1

        ax.bar(ind+i*barwidth,column,barwidth,label=mylabel[0],color=mycolscheme)
        ax.legend(loc='best', fancybox=True, frameon=True, fontsize=9)
        ax.set_yscale('log')

#--- SET COLUMNS Y-LIMITS HERE ---#
#       ax.set_ylim(10**(6),10**(-4))
#---------------------------------#
        i=i+1
    
    ax.set_xticks(ind+i*0.5*barwidth)
    ax.set_xticklabels(xlabels, rotation=70, horizontalalignment='center')
    ax.set_ylabel('Column Amount [DU]')
    fig.set_size_inches(width_colplot, height_colplot)
    extension = '__plot.pdf' 
    plotfile = plotpath + columnfile + extension

    setCustomHatchWidth(0.5)
    if(not auto): plt.savefig(plotfile, bbox_inches='tight')
    plt.close(fig)
	
    tpl1 = datetime.now()
    if(not auto):
        print '\nPlotting took %s [h:m:s]' %(tpl1-tpl0)
	
        print '\n**********************************'
        print 'output: %s' %plotfile
    print '**********************************'
# END ##### Plot Columns ##########################################


######## Plot all Profiles ########################################
def plot_profiles():
	t0 = datetime.now()
	fig = plt.figure() 
	ax = fig.add_subplot(111)
	fig.subplots_adjust(hspace=spaceh, wspace=spacew)
	if (ysecondary != ''):ax2 = ax.twinx()
	
	im=0
	for plots in molecules[:]:
            legs= []
            labels= []
            legs2= []
            labels2= []
            isel=im
            irank=isel+3
	    im=im+1
	    i=0
	    if (im > 1 and overplot == 0):
	        if (ysecondary != ''):
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

            #Plot Meadows+ 2018 results
            if(printMeadows):
                if(plots=='T'):
                    if(yselection=='P'):ax.plot(Tmeadows,Pmeadows,label='Meadows+ 2018, Earth-like',linestyle=lines[2],color='black',linewidth=1)
                    if(ysecondary=='P'):ax2.plot(Tmeadows,Pmeadows,label='Meadows+ 2018, Earth-like',linestyle=lines[2],color='black',linewidth=1)

            xmin = 1e100
            xmax = -1     
	    for runs in compared[:]: 
                if(not auto): print runs,': ',plots
	        tc1 = datetime.now()
                if (vmrmolcc == 0): data = alldata[i]
                elif (vmrmolcc == 1): data = molccdata[i] 
                else: data = mfracdata[i]

                ratio=False

	        if (plots == 'NOx'):
	            xdata=data[:,list(titles).index('N')]+data[:,list(titles).index('NO')]+data[:,list(titles).index('NO2')]+data[:,list(titles).index('NO3')]
	        elif (plots == 'HOx'):
	            xdata=data[:,list(titles).index('H')]+data[:,list(titles).index('OH')]+data[:,list(titles).index('HO2')]
	        elif (plots == 'ClOx'):
	            xdata=data[:,list(titles).index('Cl')]+data[:,list(titles).index('ClO')]+data[:,list(titles).index('OClO')]
                elif (plots == 'OH/HO2'):
                    xdata=data[:,list(titles).index('OH')]/data[:,list(titles).index('HO2')]
                    ratio=True
                elif (plots == 'NO/NO2'):
                    xdata=data[:,list(titles).index('NO')]/data[:,list(titles).index('NO2')]
                    ratio=True
                elif (plots == 'Cl/ClO'):
                    xdata=data[:,list(titles).index('Cl')]/data[:,list(titles).index('ClO')]
                    ratio=True
                elif (plots == 'ClO/ClO2'):
                    xdata=data[:,list(titles).index('ClO')]/data[:,list(titles).index('OClO')]
                    ratio=True
                elif (plots == 'O/O1D'):
                    xdata=data[:,list(titles).index('O')]/data[:,list(titles).index('O1D')]
                    ratio=True
	        else:
	            xdata=data[:,list(titles).index(plots)]

#                print 'SUM(data[:,%s]): ' %(plots), sum(xdata[:])
       
                ax.set_title('Title',alpha=0.0) #for plot(bbox_inches('tight')) if legend is above figure
	        
	        if (yselection == 'P' or yselection == 'Density'): 
                    ax.set_yscale('log')
                    ax.yaxis.set_major_formatter(CustomTicker())
	        if ( plots != 'T' and plots != 'Weight'): 
                    if (sum(xdata[:]) == 0.0): ax.set_xscale('linear')
                    else: 
                        ax.set_xscale('log') 
                        ax.xaxis.set_major_formatter(CustomTicker())
	        ax.set_ylim([data[0,list(titles).index(yselection)],data[-1,list(titles).index(yselection)]])
	        if (not ratio and plots!= 'Weight' and plots != 'Density'): 
                    if (vmrmolcc == 0): ax.set_xlabel(r'%s [vmr]' %plots)
                    elif (vmrmolcc == 1): ax.set_xlabel(r'%s [Molec. / cm$^3$]' %plots)
                    else: ax.set_xlabel(r'%s [mmr]' %plots)
                elif (plots == 'Weight' or plots == 'Density'): ax.set_xlabel(r'%s %s' %(titles[list(titles).index(plots)],units[list(titles).index(plots)]))  
	        else: ax.set_xlabel(r'%s' %plots)
                if (plots=='UVA' or plots=='UVB' or plots=='UVC' or plots == 'FUV'):
                    ax.set_xscale('linear')
                    ax.set_xlabel(r'%s [W/m$^2$]' %plots)
                    surf_dosage=(data[0,list(titles).index(plots)]-0.5*(data[1,list(titles).index(plots)]-data[0,list(titles).index(plots)]))
                    if (surf_dosage < 0.0): surf_dosage = 0.0
                    ax.set_title(r'%s dosage: %s [W/m$^2$]' %(plots,surf_dosage),alpha=1.0)
                ax.set_ylabel(r'%s %s' %(yselection,units[list(titles).index(yselection)]))
            
                if(np.min(xdata[:])<xmin):xmin = np.min(xdata[:])
                if(np.max(xdata[:])>xmax):xmax = np.max(xdata[:])

                if(xmin > 1e-100):
                    if ((xmax/xmin) < 20 and xmin>0.01 and xmax<1e3):
                        ax.set_xscale('linear')
#                       ax.xaxis.set_major_locator(ticker.AutoLocator())
#                       ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

	        tc2 = datetime.now()
                if(not auto): print 'setup took %s [h:m:s]' %(tc2-tc1)
	    
                mylinestyle = lines[int(math.floor(i/len(colors)))]
	        mycolscheme = colors[i%len(colors)]
                mymarker = ''
                if (uselabels == 1): mylabel = ['%s'%(simulabels[i]),'']
                else: mylabel = ['%s'%(runs),'']

	        if (overplot == 1): 
	            mycolscheme = colors[isel%len(colors)]
	            mylabel = ['%s'%(plots)]
	
                if (showmarkers==0):mymarker=''
                ax.plot(xdata,data[:,list(titles).index(yselection)],label=mylabel[0],marker=mymarker,linestyle=mylinestyle,color=mycolscheme,linewidth=1)
	        if (ysecondary != ''):
	            ax2.plot(xdata,data[:,list(titles).index(ysecondary)],alpha=0)
	            ax2.set_ylabel(r'%s %s' %(ysecondary,units[list(titles).index(ysecondary)]))
	            if (ysecondary == 'P' or ysecondary == 'Density'): ax2.set_yscale('log')
	            ax2.set_ylim([data[0,list(titles).index(ysecondary)],data[-1,list(titles).index(ysecondary)]])
	                            
                #ax.legend(loc='upper right',bbox_to_anchor=(1.70,1.02),ncol=1, fancybox=True, frameon=False, fontsize=6)
                ax.legend(loc='best', fancybox=True, frameon=True, fontsize=6)
	        #print '\n'
	        #break    
	        i=i+1

#               ax.xaxis.get_major_formatter().set_scientific(True)
#                ax.xaxis.set_major_formatter(CustomLinTicker())
#               ax.set_xlim(5*10**(-2),6*10**(-2))
	#sys.exit()
#--- SET COLUMNS Y-LIMITS HERE ---#
#       ax.set_ylim(10**(6),10**(-4))
#       ax.set_xlim(150,700)
#---------------------------------#
        ax.grid(which='major',linestyle=':')
	t1 = datetime.now()
	print 'Creation took %s [h:m:s]' %(t1-t0)
	
	tpl0 = datetime.now()
	
	if (overplot == 0):
            k=int((im-1.)/plotcols_profiles +1.)
	    fig.set_size_inches(width*plotcols_profiles + spacew*width*(plotcols_profiles-1), height*k + spaceh*height*(k-1)) #whatever looks good and keeps the ratio for different figure numbers
	else:
	    fig.set_size_inches(width, height)
	
        if (vmrmolcc == 0): extension = '-vmr__plot.pdf' 
        elif (vmrmolcc == 1): extension = '-molcc__plot.pdf'
        else: extension = '-mmr__plot.pdf'
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

if (selection[0] == -1):  compared = simu
else: compared = [simu[l] for l in selection]
if(auto): compared = [expdir for i in range(0,1)]

read_alldata() #... reading data files


#------------------------------------------------------------------
if (len(mselection) == 1 and mselection[0] != -1):plotcols_profiles = 1
if (mselection[0] == -1): molecules = [titles[l] for l in range(3,len(titles))]
else: molecules = mselection    
if (len(molecules) <(plotcols_profiles+1)): plotcols_profiles = len(mselection)    
#------------------------------------------------------------------


if (printMeadows): Meadows() #... for overplotting literature data


#--- Now Plotting -------------------------------------------------
if (vmrmolcc == -1):
    for ivmr in range(0,3):
        vmrmolcc = ivmr
        if (plotprofiles == 1): plot_profiles()
else: 
    if (plotprofiles == 1): plot_profiles()

if (plotcolumns ==1): plot_columns()

###################################################################
########                ###########################################
########  END PROGRAM   ###########################################
########                ###########################################
###################################################################
