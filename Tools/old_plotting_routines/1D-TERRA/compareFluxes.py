#!/usr/bin/python

'''
Routine to plot fluxes of one or multiple model simulations
  * UVA (317.5 - 405 nm)
  * UVB (281.7 - 317.5 nm)
  * UVC (175 - 281.7 nm)
  * FUV (100 - 175 nm)
  The data are taken from ../Experiments/EXP/Output/Photo/composition_vmr.dat

General notes
  * always put "-" at the beginning of an optional parameter
  * for more than one experiment use " or ' to give them as one argument
  * have a look on a few examples below

Usage: 
  python compareFluxes.py "EXP1 EXP2"
  * this will display profiles of all UV fluxes

Select fluxes:
  python compareExp.py "EXP1 EXP2" "UVA UVB"

Optional parameter for data:
  -height       y-axis will be in height(km) not pressure(hPa)
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
  -sqaxis       x-axis will only show data range for displayed heights
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

def getProf(ifile,specname,getH=False):
  '''
  Get profile and p or z of a given species or species combination
  Input
    ifile:  input file (path+allchemmolcc.dat or path+allchemvmr.dat)
    specie: single specie or combination (NOx, HOx)
  '''
  
  f = open(ifile, 'r')
  specnames = np.array(f.readline().split()[:])
  f.close()
  data=np.loadtxt(ifile,skiprows=2)
  if getH:
    y=np.squeeze(data[:,0])
  else:
    y=np.squeeze(data[:,2])
  
  k=np.where(specname==specnames)[0]
  if k:
    spec= np.squeeze(data[:,k])
  else:
    spec=np.zeros(y.size)
    print(specname, 'not present in data!')
  
  return spec,y

# plot #######################################################

exps   = sys.argv[1]

rspecs = "UVA UVB UVC FUV"
rspecs_cp = np.array(["UVA","UVB","UVC","FUV"])

uvtit=np.array(["UVA (317.5 - 405 nm)","UVB (281.7 - 317.5 nm)","UVC (175 - 281.7 nm)","FUV (100 - 175 nm)"])

if (len(sys.argv)>2):
  sys2=sys.argv[2]
  if (sys2[0:1] != "-"):
    rspecs = sys.argv[2]
    del sys.argv[2]

exps = np.array(exps.split())
nexp=exps.size

rspecs = np.array(rspecs.split())
nsp=rspecs.size

idxtit=[True]*4
for i in range(rspecs_cp.size):
  idxtit[i] = rspecs_cp[i] in rspecs

uvtit = uvtit[idxtit]

# Optional parameters
subdir=''
label=np.copy(exps)
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
  elif '-figsize' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][9:]
    figsize=np.array(sys.argv[2].split()).astype('float')
    setfigsize=True
  else:
    print('Parameter: ',sys.argv[2],' not found')
    exit()
  
  del sys.argv[2]

import matplotlib.pyplot as plt

if nset:
  nx=nxny[0]
  ny=nxny[1]
else:
  nx=1
  ny=nsp

bbox_props = dict(fc="white", ec="gray", lw=1)

fnexp=''
for i in range(nexp):
  fnexp+='_'+exps[i]

if not setfigsize:
  figsize=[3*ny+1,3*nx+1]

fig, axs = plt.subplots(figsize=figsize,nrows=nx,ncols=ny, sharey=True)
if nsp==1:
  axs = [axs]
else:
  axs=axs.flatten()

for i in range(nsp):
 
  for j in range(nexp):
    try:
      # load data
      ifile=rtdir+subdir+exps[j]+'/Output/Photo/composition_vmr.dat'
      spec, p = getProf(ifile,rspecs[i],getH=getH)
    except:
      ifile=rtdir+subdir+exps[j]+'/Output/Photo/composition_vmr_tmp.dat'
      spec, p = getProf(ifile,rspecs[i],getH=getH)
    

    if len(label[j])>40:
      lab = label[j][0:18]+'...'+label[j][-18:]
    else:
      lab = label[j]
    
    axs[i].loglog(spec,p,label=lab,lw=lw[j%len(lw)],ls=ls[j%len(ls)],
                                   color=color[j%len(color)])
    xlab = 'Flux [W/m$^2$]'

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
  
  if (i%ny == 0):
    if getH:
      axs[i].set_ylabel('Height [km]')
    else:
      axs[i].set_ylabel('Pressure [hPa]')
  
  #if (i>=(nx-1)*ny):
  axs[i].set_xlabel(xlab)
  
  axs[i].grid(True,linestyle=':')
  axs[i].set_ylim(ylim)
  if xset:
    axs[k].set_xlim(xlim)
  
  if (not xset) and (np.min(spec)<1e-30):
    axs[i].set_xlim([1e-30,1])

  axs[i].text(0.95, 0.9, uvtit[i],transform=axs[i].transAxes,
                    bbox=bbox_props,ha='right')

  if getH:
      axs[i].set_yscale('linear')

if not getH:
  plt.yticks(ytk,ytk)

if not setlegsp:
  legsp=i

if legend:
  axs[legsp].legend(loc=loc)

plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=0.4)

#title
if len(title)>0:
  plt.suptitle(title)
  plt.subplots_adjust(top=0.93)

if savefig:
  plt.savefig(ofn)
  print(ofn)
else:
  plt.show()
plt.close()
