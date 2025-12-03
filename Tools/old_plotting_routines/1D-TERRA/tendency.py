#!/usr/bin/python

'''
Routine to plot tendency of surface temperature or mixing ratios of any specie
for one or several experiments
  * Surface temperature for all coupled steps will be shown
  * For species there are two modes:
    * -steps  All coupled steps for 3 layers but all species (default)
    * -dt     Time dependent last coupled steps for 10 layers and selected species 
  * The data are taken from ../Experiments/EXP/Output/Photo/tendency.dat
    and ../Experiments/EXP/Output/Clima/clima_allout.tab

General notes
  * always put "-" at the beginning of an optional parameter
  * for more than one experiment use " or ' to give them as one argument
  * have a look on a few examples below

Usage: 
  python tendency.py "EXP1 EXP2"
  * this will display mixing ratio tendency of selected species
  python tendency.py "EXP1 EXP2" "T"
  * this will display tendebcy of surface temperature

Select other species:
  python compareExp.py "EXP1" "H2O CH4 CO2 O3"
   or
  python compareExp.py "EXP1" -specie="H2O CH4 CO2 O3"

Optional parameter for data:
  -steps        For species: all coupled steps for 3 layers but all species (default)
  -time         For species: only last coupled step for 10 layers and selected species
  -layer        Layer number 
                * For -step mode select betwen 1,2,3 for 101 layers -> 1,51,101
                * For -dt mode select betwen 1-10 for 101 layers -> 1,11,21,..101
  -expdir       Experiment directory. Default is: "../Experiments"
  -subdir       This directory will be added before experiment name
                EXAMPLE: python compareFluxes.py "EXP1" -subdir="Earth"
                        -> data taken from ../Experiments/Earth/EXP1/Output

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
  -nxny         Number of subplots in x and y direction (-nxny="2 5")
  -legendloc    Location of legend inside subplot (default=0: best)
  -legsubpl     On which subplot the legend will be displayed
  -nolegend     No legend will be displayed 

Output:
  When nothing is set the plot will be displayed (graphical environment required!)
  otherwise use
  -out="path/filename.extension" (The extension can be e.g. pdf, png or svg)

Examples:
python compareFluxes.py "S0.8 S1.0 S1.2" -subdir="Earth" -out="Plots/Earth_Flux_Sx.pdf"
python compareFluxes.py "S0.8 S1.0 S1.2" -height -ylim="0 80"
python compareFluxes.py "S0.8 S1.0 S1.2" -color="blue green red" -ls="-" -lw=2
                     -label="S\$_{sun}\$~=~0.8 S\$_{sun}\$~=~1.0 S\$_{sun}\$~=~1.2"

F. Wunderlich, May 2020
'''

import sys
import matplotlib
import numpy as np


# experiments folder
rtdir='../Experiments/'

exp   = sys.argv[1] 
exp = np.array(exp.split())
nexp=exp.size

# Species (Optional)
rspecs = "H2O CH4 O3 CO2 CO O2"
if (len(sys.argv)>2):
  sys2=sys.argv[2]
  if (sys2[0:1] != "-"):
    rspecs = sys.argv[2]
    del sys.argv[2]

# Optional parameter
layer=np.arange(1,10,2).astype('int')
ncoup=1
subdir=''
label=np.copy(exp)
color=['C'+str(x) for x in range(10)]
savefig=False
getH=False
sqaxis=False
ls=['-','--','-.',':']
lw=[2,2]
xset=False
yset=False
title=''
setlegsp=False
legend=True
nset=False
loc=0
setfigsize=False
nasx=True
unit='s'

while len(sys.argv)>2:
  # Optional subdir
  if '-subdir' in sys.argv[2]:
    subdir = sys.argv[2][8:]
    subdir+='/'
  elif '-label' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][7:]
    label=np.array(sys.argv[2].split())
    label=[lab.replace('~' , ' ') for lab in label]
    label=[lab.replace('\$' , '$') for lab in label]
  elif '-color' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][7:]
    color=np.array(sys.argv[2].split())
  elif '-ls' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][4:]
    ls=np.array(sys.argv[2].split())
  elif '-lw' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][4:]
    lw=np.array(sys.argv[2].split()).astype('float')
    setlw=True
  elif '-layer' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][7:]
    layer=np.array(sys.argv[2].split()).astype('int')
  elif '-time' in sys.argv[2]:
    nasx=False
  elif '-steps' in sys.argv[2]:
    nasx=True
  elif '-unit' in sys.argv[2]:
    unit = sys.argv[2][6:]
  elif '-ylim' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][6:]
    ylim=np.array(sys.argv[2].split()).astype('float')
    yset=True
  elif '-xlim' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][6:]
    xlim=np.array(sys.argv[2].split()).astype('float')
    xset=True
  elif '-title' in sys.argv[2]:
    title = sys.argv[2][7:]
    title=title.replace('\$' , '$')
  elif '-out' in sys.argv[2]:
    savefig=True
    ofn = sys.argv[2][5:]
    matplotlib.use('Agg')
  elif '-sqaxis' == sys.argv[2]:
    sqaxis=True
  elif '-height' == sys.argv[2]:
    getH=True
  elif '-legsubpl' in sys.argv[2]:
    legsp = int(sys.argv[2][10:])
    setlegsp=True
  elif '-legendloc' in sys.argv[2]:
    loc = int(sys.argv[2][11:])
  elif '-nolegend' == sys.argv[2]:
    legend=False
  elif '-nxny' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][6:]
    nxny=np.array(sys.argv[2].split()).astype('int')
    nset=True
  elif '-specie' in sys.argv[2]:
    rspecs = sys.argv[2][8:]
  elif '-figsize' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][9:]
    figsize=np.array(sys.argv[2].split()).astype('float')
    setfigsize=True
  else:
    print('Parameter: ',sys.argv[2],' not found')
    exit()
  
  del sys.argv[2]

import matplotlib.pyplot as plt

# define some variables
rspecs = np.array(rspecs.split())
nsp=rspecs.size

if (rspecs[0]=="T"):
  Tprof=True
  nasx=True
else:
  Tprof=False

nl=layer.size
rtdir=rtdir+subdir

if nsp == 15:
  nx=3
  ny=5
elif nsp==3:
  nx=1
  ny=3
else:
  nx=np.round(np.sqrt(nsp)).astype('int')
  ny=np.ceil(1.*nsp/nx).astype('int')

ylab='vmr'

if Tprof:
  label=exp
  ylab='Surface Temperature [K]'
  nl=1
  
elif nasx:
  ifile = rtdir+exp[0]+'/Output/Tendency/steps1.dat'
  layer=np.array([1,2,3])
  nl=layer.size
  label=['Surface','Middle Atmosphere','Top of Atmosphere']
else:
  ifile = rtdir+exp[0]+'/Output/Tendency/steps_nc001_L01.dat'

# read species names
if not Tprof:
  f = open(ifile, 'r') 
  specnames = np.array(f.readline().split()[1:])
  f.close()

# plotting
bbox_props = dict(fc="white", ec="gray", lw=1)
fig, axs = plt.subplots(figsize=(5*ny+1,3*nx+1),nrows=nx,ncols=ny, sharex=True)
if nsp==1:
  axs = [axs]
else:
  axs=axs.flatten()

# loop over layer and species
for e in range(nexp):
  for i in range(nsp):
  
    for j in range(nl):
      
      if not Tprof:
        ispec = np.where(rspecs[i]==specnames)[0][0]
      
      if Tprof:
        try:
          ifile = rtdir+exp[e]+'/Output/Clima/clima_allout'
          fo=open(ifile, "r")
          fo.close()
        except:
          ifile = rtdir+exp[e]+'/Output/Clima/clima_allout.tab'
      elif nasx:
        ifile = rtdir+exp[e]+'/Output/Tendency/steps'+str(layer[j])+'.dat'
      else:
        fn='steps_nc'+str(ncoup).zfill(3)+'_L'+str(layer[j]).zfill(2)+'.dat'
        ifile = rtdir+exp[e]+'/Output/Tendency/'+fn
      
      if Tprof:
        vmr=[]
        n=0
        fo = open(ifile, "r")
        lines = fo.readlines()
        try:
          lines=np.asarray(lines)
        except:
          print('Last part of lines displayed due to memory error')
          nl=len(lines)
          lines=np.asarray(lines[-nl/3:])
        strfnd='JCONV='
        iTln=np.flatnonzero(np.core.defchararray.find(lines,strfnd)!=-1)
        vmr=(np.char.partition(lines[iTln],' T(ND)=')[:,2]).astype('float')
        
        #with open(ifile) as openfile:
        #  for line in openfile:
        #    iT0 = line.find(' T(ND)=')
        #    if iT0>=0:
        #      vmr=vmr+[float(line[iT0+7:])]
        #      n+=1

        nn=range(n)  
        
      elif nasx:
        try:
          nn,vmr = np.loadtxt(ifile,unpack=True,usecols=(1,ispec))
        except:
          nn,vmrt = np.loadtxt(ifile,unpack=True,usecols=(1,ispec),dtype='str')
          nn=nn.astype('float')
          vmr=np.empty(vmrt.size)
          for v in range(vmr.size):
            try:
              vmr[v] = float(vmrt[i])
            except:
              vmr[v] = 1e-99
      else:
        p,t,vmr = np.loadtxt(ifile,unpack=True,usecols=(1,2,ispec))
        p=p[0]/1000
        t=np.cumsum(t)

      if Tprof:
        axs[i].plot(vmr,label=label[e],lw=lw[e%len(lw)],
                      ls='-',color=color[e%len(color)])
        axs[i].set_xlabel('Steps')
      elif nasx:
        if nexp>1:
          lab=exp[e]+' ('+label[j]+')'
        else:
          lab=label[j]
        axs[i].semilogy(vmr,label=lab,lw=lw[e%len(lw)],ls=ls[e%len(ls)],
                                      color=color[j%len(color)])
        axs[i].set_xlabel('Steps')
      else:
        if nexp>1:
          lab=exp[e]+' ('+str(p)+'hPa)'
        else:
          lab=str(p)+'hPa'

        if unit=="y":
          yins=60.*60.*24.*365.25
          t=t/yins
        axs[i].loglog(t,vmr,label=lab,lw=lw[e%len(lw)],
                      ls=ls[e%len(ls)],color=color[j%len(color)])
        axs[i].set_xlabel('Time['+unit+']')
    
    axs[i].set_ylabel(ylab)
    axs[i].text(0.15, 0.9, rspecs[i],transform=axs[i].transAxes,
                      bbox=bbox_props,ha='right')
    if not nasx:
      if unit=="y":
        tmin = 1/yins
      else:
        tmin=1e-3
      axs[i].set_xlim([tmin,t[-1]])

plt.legend(loc=0)
plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=0.4)
if savefig:
  if len(ofn)>0: plt.savefig('../Experiments/plots/'+str(ofn))
  else: plt.savefig('pltTend'+str(nsp)+'_'+exp[0]+'.png')
else:
  plt.show()
plt.close()
