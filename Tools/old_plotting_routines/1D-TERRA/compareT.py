#!/usr/bin/python

'''
Routine to plot temperature of one or multiple model simulations
  The data which can be displayed are: 
  * Temperature profile, T, in K
  * Scatter plot and/or line plot of surface temperature, Tsurf, in K  
  The data are taken from ../Experiments/EXP/Output/Clima/Tprofile.dat

General notes
  * always put "-" at the beginning of an optional parameter
  * for more than one experiment use " or ' to give them as one argument
  * have a look on a few examples below

Usage: 
  python compareT.py "EXP1 EXP2"
  * this will display a standard profiles of T
  python compareT.py "EXP1 EXP2" -Tsurf
  * this will display a plot only with surface temperatures

Optional parameter for data:
  -Tsurf        Show only surface temperature (y-axis will be surface temperature)
  -height       y-axis will be in height(km) not pressure(hPa)
  -expdir       Experiment directory. Default is: "../Experiments"
  -subdir       This directory will be added before experiment name
                EXAMPLE: python compareT.py "EXP1" -subdir="Earth"
                        -> data taken from ../Experiments/Earth/EXP1/Output

Optional parameter for Tsurf plot:
  -xlab         Set x-label (default is Number)
  -xtick        Set x-ticklabels (default is 0,1...n)
  -xlog         x-axis in log scale (default is linear)  
  -marker       Style of marker (default is )
  -line         Line plot of surface temperatures
  -noscatter    Omit scatter plot

Optional parameter for plot style:
  -figsize      Set figure size in inches -figsize="xsize ysize"
  -title        Set title for entire plot: use e.g. \$ \lamda \$ for math 
  -label        Set label for each experiment: use ~ for space and \$ for math
  -color        Set color for each experiment (will be repeated if ncolor<nexp) 
                Default: "C0 C1 .. C9" (tableau palette) 
                (see here https://matplotlib.org/3.1.0/gallery/color/named_colors.html)
  -ls           Set line style for each experiment (will be repeated if nls<nexp)
                Default: "- -- -. :"
  -lw           Set line width for each experiment (will be repeated if nlw<nexp)
  -ylim         Set limits of y-axis: for height in km e.g. -ylim="0 100" 
                                      for pressure in hPa e.g. -ylim="1000 0.1"
  -xlim         Set limits of x-axis. This will be the same for each subplot
  -sqaxis       x-axis will only show data range for displayed heights
  -legendloc    Location of legend inside subplot (default=0: best)
  -nolegend     No legend will be displayed

Output:
  When nothing is set the plot will be displayed (graphical environment required!)
  otherwise use
  -out="path/filename.extension" (The extension can be e.g. pdf, png or svg)

Examples:
python compareT.py "S0.8 S1.0 S1.2" -subdir="Earth" -out="Plots/Earth_T_Sx.pdf"
python compareT.py "S0.8 S1.0 S1.2" -Tsurf -xtick="0.8 1.0 1.2" -xlab="Instellation"

F. Wunderlich, May 2020
'''

import sys
import matplotlib
import numpy as np
import matplotlib.colors as mplc

cols=['darkmagenta','midnightblue','darkblue','C0', 'C2', 'orange','C3','darkred','k']
cmap = mplc.LinearSegmentedColormap.from_list('Tmap', cols)
norm = mplc.Normalize(vmin=170,vmax=400)


# experiments folder
rtdir='../Experiments/'

def getProf(ifile,getH=False):
  '''
  Get profile and p or z of a given species or species combination
  Input
    ifile:  input file (path+allchemmolcc.dat or path+allchemvmr.dat)
    specie: single specie or combination (NOx, HOx)
  '''
  
  data=np.loadtxt(ifile,skiprows=1)
  if getH:
    y=np.squeeze(data[:,1])
  else:
    y=np.squeeze(data[:,0])*1000.
  
  spec= np.squeeze(data[:,2])
  
  return spec,y

# plot #######################################################

exps   = sys.argv[1] 
exps = np.array(exps.split())
nexp=exps.size

# Optional parameters
subdir=''
label=np.copy(exps)
color=['C'+str(x) for x in range(10)]
savefig=False
tmpdat=False
getH=False
ls=['-','--','-.',':']
surfT=False
xlab = 'Experiment Number'
xtick=range(nexp)
xlog=False
yset=False
sqaxis=False
yset=False
xset=False
setlab=False
ls=['-','--','-.',':']
lw=[2,2]
loc=0
ofn=''
title=''
figsize=(4,4)
setlw=False
setcol=False
legend=True
line=False
marker='o'
setcol=False
scatter=True

while len(sys.argv)>2:
  # Optional Tsurf
  if '-Tsurf' == sys.argv[2]:
    surfT=True
  elif '-xlog' == sys.argv[2]:
    xlog=True
  elif '-xlab' in sys.argv[2]:
    xlab = sys.argv[2][6:]
  elif '-xtick' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][7:]
    xtick=np.array(sys.argv[2].split()).astype('float')
  elif '-line' == sys.argv[2]:
    line=True
  elif '-noscatter' == sys.argv[2]:
    scatter=False
  elif '-marker' in sys.argv[2]:
    marker = sys.argv[2][8:]
  #Optional all others
  elif '-subdir' in sys.argv[2]:
    subdir = sys.argv[2][8:]
    subdir+='/'
  elif '-expdir' in sys.argv[2]:
    rtdir = sys.argv[2][8:]
    rtdir+='/'
  elif '-label' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][7:]
    label=np.array(sys.argv[2].split())
    label=[lab.replace('~' , ' ') for lab in label]
    label=[lab.replace('\$' , '$') for lab in label]
    setlab=True
  elif '-color' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][7:]
    color=np.array(sys.argv[2].split())
    setcol=True
  elif '-ls' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][4:]
    ls=np.array(sys.argv[2].split())
  elif '-lw' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][4:]
    lw=np.array(sys.argv[2].split()).astype('float')
    setlw=True
  elif '-title' in sys.argv[2]:
    title = sys.argv[2][7:]
    title=title.replace('\$' , '$')
  elif '-ylim' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][6:]
    ylim=np.array(sys.argv[2].split()).astype('float')
    yset=True
  elif '-xlim' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][6:]
    xlim=np.array(sys.argv[2].split()).astype('float')
    xset=True
  elif '-figsize' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][9:]
    figsize=np.array(sys.argv[2].split()).astype('float')
  elif '-nolegend' == sys.argv[2]:
    legend=False
  elif '-legendloc' in sys.argv[2]:
    loc = int(sys.argv[2][11:])
  elif '-out' in sys.argv[2]:
    savefig=True
    ofn = sys.argv[2][5:]
    matplotlib.use('Agg')
  elif '-sqaxis' == sys.argv[2]:
    sqaxis=True
  elif '-height' == sys.argv[2]:
    getH=True
  else:
    print('Parameter: '+sys.argv[2]+' not found')
    exit()
  
  del sys.argv[2]

import matplotlib.pyplot as plt
# ------------------------------------------------------------

if (not surfT):
  xlab="Temperature [K]"

bbox_props = dict(fc="white", ec="gray", lw=1)

if getH:
  if not yset:
    ylim=[0,80]
  ylab = 'Height [km]'
elif surfT:
  ylab = 'Surface Temperature[K]'
else:
  ylab = 'Pressure [hPa]'

fnexp=''
for i in range(nexp):
  fnexp+='_'+exps[i]

ts=[]

fig, axs = plt.subplots(figsize=figsize,nrows=1,ncols=1, sharey=True)

for j in range(nexp):

  ifile=rtdir+subdir+exps[j]+'/Output/Clima/Tprofile.dat'
  spec, p = getProf(ifile,getH=getH)

  # surface temperature, only data
  if surfT:
    ts=ts+[float(spec[-1])]
  # T profile plotting
  else:
    if len(label[j])>40:
      lab = label[j][0:18]+'...'+label[j][-18:]
    else:
      lab = label[j]

    axs.semilogy(spec,p,label=lab,lw=lw[j%len(lw)],
                 ls=ls[j%len(ls)],color=color[j%len(color)])

    # get y limits
    if not yset:
      if getH:
        ylim=[p.min(),p[:-3].max()]
        if ylim[0]<2:
          ylim[0]=0
      else:
        ylim=[p.max(),p[:-3].min()]
        y1=np.log10(ylim[0]) - np.round(np.log10(ylim[0]))
        if (y1<0) and (y1>-0.1):
          ylim[0]=10**np.round(np.log10(ylim[0]))
    
    ytk=np.array([10**x for x in range(5,-10,-1)])
    ytk=ytk[ytk<=ylim[0]]
    ytk=ytk[ytk>=ylim[1]]
    
    print(exps[j]+', surface temperature = '+str(spec[-1])+' K')
  if getH:
    axs.set_yscale('linear')


# T surf plotting
if surfT:
  if setcol:
    if line:
      axs.plot(xtick,ts,color=color[0],lw=lw[0],ls=ls[0])
    if scatter:
      axs.scatter(xtick,ts,c=color[0],marker=marker)
  else:
    if line:
      axs.plot(xtick,ts,color='black',lw=lw[0],ls=ls[0])
    if scatter:
      axs.scatter(xtick,ts,c=ts,cmap=cmap,norm=norm,marker=marker)
  if xlog:
    axs.set_xscale('log')

axs.set_ylabel(ylab)
axs.set_xlabel(xlab)

axs.grid(True,linestyle=':')
if not surfT:
  axs.set_ylim(ylim)

if not getH and not surfT:
  plt.yticks(ytk,ytk)
  
if surfT:
  plt.axhline(y=273.15,ls='--',color='k')
  plt.axhline(y=373.15,ls='--',color='k')

if (not surfT) and (legend):
  plt.legend(loc=loc)

plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=0.4)

#title
if len(title)>0:
  plt.suptitle(title)
  plt.subplots_adjust(top=0.93)

# save or show
if savefig:
  plt.savefig(ofn)
  print(ofn)
else:
  plt.show()
plt.close()
