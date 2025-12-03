#!/usr/bin/python

'''
Routine to plot emax and dt of one or several model simulations
  The data are taken from ../Experiments/EXP/Output/Tendency/emax.dat
  emax: Maximum relative change of all species and layer in %
        ABS(USOL_i(jSpec,jLayer)/USOL_i+1(jSpec,jLayer))*100
  dt:   Time step in s or y

General notes
  * always put "-" at the beginning of an optional parameter
  * for more than one experiment use " or ' to give them as one argument
  * have a look on a few examples below
  * horizontal line for emax is at 0.01 % (convergence if dt is big enough)
  * horizontal line for dt is at 1 year (convergence dt defined in input.dat, normally ~1Gyrs)

Usage: 
  python emax.py "EXP1 EXP2"
  * this will display emax and dt evolution

Optional parameter for data:
  -last         Will only show emax of last 
  -expdir       Experiment directory. Default is: "../Experiments"
  -subdir       This directory will be added before experiment name
                EXAMPLE: python compareFluxes.py "EXP1" -subdir="Earth"
                        -> data taken from ../Experiments/Earth/EXP1/Output

Optional parameter for plot style:
  -unit         Set -unit='y' for years
  -figsize      Set figure size in inches -figsize="xsize ysize"
  -title        Set title for entire plot: use e.g. \$ \lamda \$ for math 
  -label        Set label for each experiment: use ~ for space and \$ for math
  -color        Set color for each experiment (will be repeated if ncolor<nexp) 
                Default: "C0 C1 .. C9" (tableau palette) 
                (see here https://matplotlib.org/3.1.0/gallery/color/named_colors.html)
  -ls           Set line style for each experiment (will be repeated if nls<nexp)
                Default: "-"
  -lw           Set line width for each experiment (will be repeated if nlw<nexp)
  -xlim         Set limits of x-axis. This will be the same for each subplot
  -ylim         Set limits of y-axis. This will be the same for each subplot

Output:
  When nothing is set the plot will be displayed (graphical environment required!)
  otherwise use
  -out="path/filename.extension" (The extension can be e.g. pdf, png or svg)

Examples:
python emax.py "S0.8 S1.0 S1.2" -subdir="Earth" -out="Plots/Earth_emax.pdf"
python emax.py "S0.8 S1.0 S1.2" -unit="y" -last

F. Wunderlich, May 2020
'''

import sys
import numpy as np
import matplotlib

# experiments folder
rtdir='../Experiments/'

exp   = sys.argv[1] 
exp = np.array(exp.split())
nexp=exp.size

# Optional parameters
subdir=''
label=np.copy(exp)
color=['C'+str(x) for x in range(10)]
savefig=False
ls=['-']
lw=[2,2]
title=''
legend=True
nset=False
loc=0
setfigsize=False
last=False
unit='s'
yset=False
xset=False

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
  elif '-unit' in sys.argv[2]:
    unit = sys.argv[2][6:]
  elif '-out' in sys.argv[2]:
    savefig=True
    ofn = sys.argv[2][5:]
    matplotlib.use('Agg')
  elif '-last' in sys.argv[2]:
    last=True
  elif '-legendloc' in sys.argv[2]:
    loc = int(sys.argv[2][11:])
  elif '-figsize' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][9:]
    figsize=np.array(sys.argv[2].split()).astype('float')
    setfigsize=True
  else:
    print('Parameter: ',sys.argv[2],' not found')
    exit()
  
  del sys.argv[2]

import matplotlib.pyplot as plt

if not setfigsize:
  figsize=[7,5]

fig, ax = plt.subplots(figsize=figsize,ncols=1,nrows=2)
for i in range(nexp):
  try:
    nn,dt,emax = np.loadtxt(rtdir+subdir+exp[i]+'/Output/Tendency/emax.dat',
                            unpack=True,skiprows=1,usecols=(0,2,3))
  except:
    nn,dt,emax = np.loadtxt(rtdir+subdir+exp[i]+'/Output/Tendency/emax.dat',
                            unpack=True,skiprows=1,usecols=(0,2,3),dtype='str')
    for j in range(nn.size):
      try:
        nn[j]=float(nn[j])
      except:
        nn[j] = 0.0
      try:
        dt[j]=float(dt[j])
      except:
        dt[j] = 0.0
      try:
        emax[j]=float(emax[j])
      except:
        emax[j] = 0.0
    nn=np.array(nn).astype('float')
    dt=np.array(dt).astype('float')
    emax=np.array(emax).astype('float') 
              
  
  dts=np.sum(dt)
  yins=31557600.
  print('%s total time running: %.1e s (%.1f a)' %(exp[i],dts,dts/yins))
  
  if (last):
    ii   = np.where((nn[1:]-nn[0:-1])<0)[0]
    if ii.size:
      imin = np.max(ii)+1
    else:
      imin=0
    nn = nn[imin:]
    dt = dt[imin:]
    emax = emax[imin:]
  else:
    idxcp=np.where((nn[:-1]-nn[1:])>0)[0]+1
    idxcp=np.append(idxcp,nn.size)
    for j in range(idxcp.size-1):
      nn[idxcp[j]:idxcp[j+1]] = nn[idxcp[j]-1]+nn[idxcp[j]:idxcp[j+1]]

  if len(exp[i])>40:
    lab = exp[i][0:18]+'...'+exp[i][-18:]
  else:
    lab = exp[i]

  if unit=='y':
    dt=dt/yins
  p = ax[0].semilogy(nn, dt,lw=lw[i%len(lw)],ls=ls[i%len(ls)],
                                   color=color[i%len(color)])
  ax[1].semilogy(nn, emax*100,label=lab,lw=lw[i%len(lw)],ls=ls[i%len(ls)],
                                   color=color[i%len(color)])
  
  
  if (not last):
    for xc in idxcp[:-1]:
      ax[0].axvline(x=nn[xc], color=p[0].get_color(), linestyle=':')
      ax[1].axvline(x=nn[xc], color=p[0].get_color(), linestyle=':')

if unit=='y':
  yins=1

ax[0].axhline(yins,ls=":",color="black")
ax[0].axhline(1e7*yins,ls=":",color="black")
ax[1].axhline(1e-2,ls=":",color="black")

if xset:
  ax[0].set_xlim(xlim)
  ax[1].set_xlim(xlim)
if yset:
  ax[0].set_ylim(ylim)
  ax[1].set_ylim(ylim)

#title
if len(title)>0:
  plt.suptitle(title)
  plt.subplots_adjust(top=0.93)

ax[0].set_xlabel('Step')
ax[0].set_ylabel('dt ('+unit+')')
plt.legend()
ax[1].set_ylabel('Emax (%)')
ax[1].set_xlabel('Step')
fig.tight_layout()

# save or show
if savefig:
  plt.savefig(ofn)
  print(ofn)
else:
  plt.show()
plt.close()
