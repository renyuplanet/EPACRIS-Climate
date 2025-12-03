#!/usr/bin/python

'''
Routine to plot output of one or multiple model simulations
  The data which can be displayed are: 
  * Profiles of any available specie in VMR or molec/cm3
  * Temperature profile, T, in K
  * Eddy diffusion coefficients profile, Edd, in cm2/s
  The data are taken from ../Experiments/EXP/Output

General notes
  * always put "-" at the beginning of an optional parameter
  * for more than one experiment use " or ' to give them as one argument
  * have a look on a few examples below

Usage: 
  python compareExp.py "EXP1 EXP2"
  * this will display a standard profiles of T and some species

Select other species:
  python compareExp.py "EXP1 EXP2" "T Edd H2O CH4 CO2 O3"
   or
  python compareExp.py "EXP1 EXP2" -specie="T Edd H2O CH4 CO2 O3"
   or
  python compareExp.py "EXP1 EXP2" "all"  
  (all species will be plotted)

Optional parameter for data:
  -molcc        Profile will be in molec/cm3
  -height       y-axis will be in height(km) not pressure(hPa)
  -expdir       Experiment directory. Default is: "../Experiments"
  -subdir       This directory will be added before experiment name
                EXAMPLE: python compareExp.py "EXP1" -subdir="Earth"
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
  -spaghetti    Several species in one single plot w/o subpots
  -defaultxaxis x-axis not extend if less than one order of magnitude is displayed  
  -sqaxis       x-axis will only show data range for displayed heights
  -legendloc    Location of legend inside subplot (default=0: best)
  -legsubpl     On which subplot the legend will be displayed
  -customleg    Use a custom legend using "color1 lw1 ls1 color2 lw2 ls2 ..."
  -nolegend     No legend will be displayed

Compare to data or observations
  -data        Give path and beginning of filename to data directory (e.g. observations or a paper)
               -data="Earth/MIPAS_obs" will use Earth/MIPAS_obs_vmr.dat
  -shade       Shade for observation range   

Output:
  When nothing is set the plot will be displayed (graphical environment required!)
  otherwise use
  -out="path/filename.extension" (The extension can be e.g. pdf, png or svg)

Examples:
python compareExp.py "S0.8 S1.0 S1.2" "T H2O CH4 O3" -subdir="Earth" -out="Plots/Earth_Specs_Sx.pdf"
python compareExp.py "S0.8 S1.0 S1.2" "T H2O CH4 O3" -molcc -height -ylim="0 80"
python compareExp.py "S0.8 S1.0 S1.2" "T H2O CH4 O3" -color="blue green red" -ls="-" -lw=2
                     -label="S\$_{sun}\$~=~0.8 S\$_{sun}\$~=~1.0 S\$_{sun}\$~=~1.2"

F. Wunderlich, May 2020
'''

import os
import sys
import matplotlib
import numpy as np
from matplotlib.ticker import ScalarFormatter

# experiments folder
rtdir='../Experiments/'

def getProf(ifile,specname,getH=False,skiprows=2):
  '''
  Get profile and p or z of a given species or species combination
  Input
    ifile:  input file (path+allchemmolcc.dat or path+allchemvmr.dat)
    specie: single specie or combination (NOx, HOx)
  '''
  
  sumspecnames = {'NOx':['N','NO2','NO','NO3'],
                  'HOx':['H','OH','HO2'],
                  'HxOx':['H','OH','HO2','H2O2'],
                  'HNOx':['HNO','HNO2','HNO3','HO2NO2'],
                  'ClOx':['Cl','ClO','ClO2','Cl2O2'],
                  'CxHx':['C4H61','C4H63','C4H62','C3H5',
                          'C3H2','C3H3','C3H8','C2H4','C2H2',
                          'C6H2','C6H6','C4H10','C2H5','C2H3',
                          'CH','C2H','C6H3','2C3H4','C8H2',
                          'C4H9','C4H8','C4H3','C4H2','C4H5',
                          'C4H4','C3H7','C2H6','1C3H4','C3H6']}
  
  # get species name and find type of input data
  #try:
  f = open(ifile, 'r')
  specnames = np.array(f.readline().split()[:])
  unitsline = f.readline()
  units     = np.array(unitsline.split()[:])
  f.close()
  #except:
  #  specnames=np.array([''])
  #  units=''
  #  print ifile, 'not found, use -tmp to show temporary result'
    

  getrnx=False
  getrmnx=False
  if ('Min' in unitsline and 'Max' in unitsline):
    getrnx=True
    if ('Mean' in unitsline):
      getrmnx=True
    else:
      getrmnx=False
    
  
  # get data
  #try:
  data=np.loadtxt(ifile,skiprows=skiprows,dtype='str')
  #except:
  #  data=np.array([[np.nan]])
  
  if getH:
    ihei=np.where(specnames=='Height')[0]
    y=np.squeeze(data[:,ihei]).astype('float')
  else:
    try:
      ipress=np.where(specnames=='P')[0][0]
      y=np.squeeze(data[:,ipress]).astype('float')
    except:
      ipress=np.where(specnames=='Pressure')[0][0]
      y=np.squeeze(data[:,ipress]).astype('float')*1e-3
  
  if specname in sumspecnames:
    i=0
    for subspecname in sumspecnames[specname]:
      k=np.where(subspecname==specnames)[0]
      if (i==0):
        if k:
          spec= np.squeeze(data[:,k]).astype('float')
        else:
          print(subspecname+ ' not present in data!')
      else:
        if k:
          spec+= np.squeeze(data[:,k]).astype('float')
        else:
          print(subspecname+ ' not present in data!')  
      i+=1
    
    if not 'spec' in locals():
      spec=np.zeros(y.size) 
    
    spec/=i 
  else:
    k=np.where(specname==specnames)[0]

    spec=[]
    if getrmnx and k:
      spec += [np.squeeze(data[:,(3*k-4)]).astype('float')]
      spec += [np.squeeze(data[:,(3*k-3)]).astype('float')]
      spec += [np.squeeze(data[:,(3*k-2)]).astype('float')]
    elif getrnx and k:
      spec += [np.squeeze(data[:,(2*k-2)]).astype('float')]
      spec += [np.squeeze(data[:,(2*k-1)]).astype('float')]
    elif k:
      spec= np.squeeze(data[:,k]).astype('float')
    else:
      spec=np.zeros(y.size)
      print(specname + ' not present in data!')

  return np.abs(np.array(spec)),y

def getErr(ifn,spec):
  '''
  Get x,y and corresponding error of obs_errobar file
  '''
  
  dat = np.loadtxt(ifn,dtype='str')
  dat = dat[(dat[:,0]==spec),:] 
  
  x=dat[:,1].astype('float')
  y=dat[:,2].astype('float')
  xerr=np.transpose(dat[:,3:5].astype('float'))
  yerr=np.transpose(dat[:,5:].astype('float'))
  
  return x,y,xerr,yerr

# plot #######################################################

# Read in variables -----------------------------------------

# Experiments
expin   = sys.argv[1] 
exps = np.array(expin.split())
nexp=exps.size

for e in exps:
  if '*' in e:
    edir=e.replace('*','')
    dirs = os.listdir('../Experiments/'+edir)
    addexp=''
    for d in dirs:
      addexp=addexp+edir+d+' '
    
    expin = expin.replace(e,addexp)

exps = np.array(expin.split())
nexp=exps.size

# Species (Optional)
rspecsdf = "T H2O CO2 N2 CO O2 O3 CH4 OH NO NO2 N2O HNO3 H2 H O SO2 OCS HCl CH3Cl"
rspecs=''

if (len(sys.argv)>2):
  sys2=sys.argv[2]
  if (sys2[0:1] != "-"):
    rspecs = sys.argv[2]
    del sys.argv[2]

# Optional parameters
subdir=''
label=np.copy(exps)
colors=['C'+str(x) for x in range(10)]
savefig=False
datfn=[]
molcc=False
utype='vmr'
getH=False
sqaxis=False
spaghetti=False
getdat=False
yset=False
xset=False
setlab=False
addh=0
ls=['-','--','-.',':']
lw=[2,2]
loc=0
ofn=''
title=''
figsize=(0,0)
extxaxis=True
shade=False
setlw=False
setcol=False
customleg=False
setlegsp=False
nset=False
legend=True


while len(sys.argv)>2:
  # Optional subdir
  if '-subdir' in sys.argv[2]:
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
    colors=np.array(sys.argv[2].split())
    setcol=True
  elif '-specie' in sys.argv[2]:
    rspecs = sys.argv[2][8:]
  elif '-ls' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][4:]
    ls=np.array(sys.argv[2].split())
  elif '-lw' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][4:]
    lw=np.array(sys.argv[2].split()).astype('float')
    setlw=True
  elif '-addh' in sys.argv[2]:
    addh = float(sys.argv[2][6:])
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
  elif '-data' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][6:]
    datfn=np.array(sys.argv[2].split())
    getdat=True
  elif '-legsubpl' in sys.argv[2]:
    legsp = int(sys.argv[2][10:])
    setlegsp=True
  elif '-customleg' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][11:]
    clegs=np.array(sys.argv[2].split())
    customleg=True
  elif '-nolegend' == sys.argv[2]:
    legend=False
  elif '-legendloc' in sys.argv[2]:
    loc = int(sys.argv[2][11:])
  elif '-obs' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][5:]
    datfn=np.array(sys.argv[2].split())
    getdat=True
  elif '-out' in sys.argv[2]:
    savefig=True
    ofn = sys.argv[2][5:]
    matplotlib.use('Agg')
  elif '-shade' == sys.argv[2]:
    shade=True
  elif '-defaultxaxis' == sys.argv[2]:
    extxaxis=False
  elif '-molcc' == sys.argv[2]:
    utype='molcc'
    molcc=True
  elif '-sqaxis' == sys.argv[2]:
    sqaxis=True
  elif '-height' == sys.argv[2]:
    getH=True
  elif '-spaghetti' == sys.argv[2]:
    spaghetti=True
  elif '-nxny' in sys.argv[2]:
    sys.argv[2] = sys.argv[2][6:]
    nxny=np.array(sys.argv[2].split()).astype('int')
    nset=True
  else:
    print('Parameter: '+sys.argv[2]+' not found')
    exit()
  
  del sys.argv[2]

import matplotlib.pyplot as plt
# ------------------------------------------------------------

if (rspecs=='all'):
  rspecs=[]
  for i in range(nexp):
    ifile=rtdir+subdir+exps[i]+'/Output/Photo/composition_vmr.dat'
    f = open(ifile, 'r')
    rspecs += f.readline().split()[3:-7]
    f.close()
  
  rspecs=' '.join(set(rspecs))
  

if rspecs=='':
  if spaghetti:
    rspecs=rspecsdfspag
  else:
    rspecs=rspecsdf
rspecs = np.array(rspecs.split())
nsp=rspecs.size

if getdat:
  ndat=datfn.size
else:
  ndat=1
  
if not setlab and getdat:
  label=np.append(label,datfn)

if getH:
  ylab = 'Height [km]'
else:
  ylab = 'Pressure [hPa]'

if nsp == 15:
  nx=3
  ny=5
elif nsp==3:
  nx=1
  ny=3
elif nsp==8:
  nx=2
  ny=4  
elif nsp==10:
  nx=2
  ny=5 
elif nsp==18:
  nx=3
  ny=6 
else:
  nx=np.round(np.sqrt(nsp)).astype('int')
  ny=np.ceil(1.*nsp/nx).astype('int')

if nset:
  nx=nxny[0]
  ny=nxny[1]

bbox_props = dict(fc="white", ec="gray", lw=1,alpha=0.75)

if figsize[0]==0:
  if spaghetti:
    figsize=(6,6)
  else:
    figsize=(3*ny+1,3*nx+1)

fnexp=''
for i in range(nexp):
  fnexp+='_'+exps[i]
  
if spaghetti:
  fig, axs = plt.subplots(figsize=figsize)
  axs = [axs]
else:
  fig, axs = plt.subplots(figsize=figsize,nrows=nx,ncols=ny, 
                          sharey=True)
  if nsp==1:
    axs = [axs]
  else:
    axs=axs.flatten()

# Loop over species ====================================================
for i in range(nsp):
  
  if spaghetti:
    k = 0
  else:
    k = i
  
  # Loop over experiments ----------------------------------------------
  for j in range(nexp):

    if (rspecs[i]=='Edd') or (rspecs[i]=='Eddy_diff.'):
      ifile=rtdir+subdir+exps[j]+'/Output/Photo/atm_composition.dat'
      spec, p = getProf(ifile,'Eddy_diff.',getH=getH,skiprows=1)
    else:
      try:
        ifile=rtdir+subdir+exps[j]+'/Output/Photo/composition_'+utype+'.dat'
        spec, p = getProf(ifile,rspecs[i],getH=getH)
      except:
        ifile=rtdir+subdir+exps[j]+'/Output/Photo/composition_'+utype+'_tmp.dat' 
        print(ifile + 'used, run not converged!')
        spec, p = getProf(ifile,rspecs[i],getH=getH)
    p=p+addh

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
      
    # change species number to subscript
    spname=rspecs[i]
    for nn in range(1,10):
      spname=spname.replace(str(nn),'$_'+str(nn)+'$')
    
    if spname=="T":
      spname="Temperature"

    if spname=="Edd":
      spname="Eddy diffusion"
    
    # squeeze axis if desired
    if sqaxis:
      if getH:
        take=p<ylim[1] 
      else:
        take=p>ylim[1]
      take[np.max(np.where(take))+1] = True
      spec=spec[take]
      p=p[take]
    
    # x label
    if molcc:
      xlab = 'Molecules cm$^{-3}$'
    else:
      xlab = 'Volume Mixing Ratio'
      
    # add legend if spaghetti
    if spaghetti:
      color=colors[i]
      lab = spname
      if j>0:
        lab=''
      if j==0 and i==0 and nexp>1:
        legexps=[]
        for l in range(nexp):
          legexp, = axs[0].plot(np.nan,np.nan,lw=lw[l%len(lw)],
                                ls=ls[l%len(ls)],color='black')
          legexps=legexps+[legexp]
        leg1=axs[0].legend(legexps,label, loc='upper right')
        axs[0].add_artist(leg1)
          
    else:
    # color and label for non-spaghetti plot
      color=colors[j%len(colors)]
      try:
        if len(label[j])>40:
          lab = label[j][0:18]+'...'+label[j][-18:]
        else:
          lab = label[j]
      except:
        lab=''  
    
    # MAIN PLOT #####################################
    axs[k].loglog(spec,p,label=lab,lw=lw[j%len(lw)],
                  ls=ls[j%len(ls)],color=color)
    ##################################################
    
    # scaling of axis
    if getH:
      axs[k].set_yscale('linear')
    if rspecs[i]=='T':
      xlab = 'Temperature [K]'
      axs[k].set_xscale('linear')
    if rspecs[i]=='Edd':
      xlab = 'Eddy diffusion [cm$^2$s$^{-1}$]'
    
    # ylabel
    if (i%ny == 0):
      axs[k].set_ylabel(ylab)
   
  # Loop over data -----------------------------------------------------
  if getdat:
    for d in range(ndat):
      if molcc:
        ifile = 'Data/'+datfn[d]+'_molcc.dat'
      else:
        ifile = 'Data/'+datfn[d]+'_vmr.dat'

      # get obs range
      getrnx=False
      getrmnx=False
      getprof=False

      if 'errorbar' in ifile:
        x,y,xerr,yerr=getErr(ifile,rspecs[i])
        geterrxy=True
      else:
        dat, p = getProf(ifile,rspecs[i],getH=getH)
        geterrxy=False
      
        dat[dat<=0.] = np.nan
      
        if len(dat)==3:
          datm,datn,datx=dat[0,:],dat[1,:],dat[2,:]
          getrmnx=True
        elif len(dat)==2:
          datm=np.array(np.nan)
          datn,datx=dat[0,:],dat[1,:]
          getrnx=True
        else:
          getprof=True
  
      j=nexp+d
       
      if (not setcol) and (ndat==1):
        color='k'
      else:
        if spaghetti:
          color=colors[i]
        else:
          color=colors[(j)%len(colors)]
        
      
      if setlw:
        lwi=lw[j%len(lw)]
      else:
        lwi=1
      
      try:
        if len(label[j])>40:
          lab = label[j][0:18]+'...'+label[j][-18:]
        else:
          lab = label[j]
      except:
        lab=''
      
      # errorbar plot
      if getrmnx or getrnx:
        if not shade:
          for e in range(datn.size):
            if e==0:
              axs[k].plot([datn[e],datx[e]],[p[e],p[e]],label=lab,
                          lw=lwi,color=color)
            else:
              axs[k].plot([datn[e],datx[e]],[p[e],p[e]],
                          lw=lwi,color=color)
        else:
          axs[i].fill_betweenx(p,datn, datx, lw=0,
                              color='gray', alpha=0.20,label=lab)
        if getrmnx:
          axs[k].loglog(datm,p,label=lab,lw=1,
                     ls='-',color='black')
      elif (geterrxy):
        axs[k].errorbar(x,y,xerr=xerr,yerr=yerr,color=color,marker='o',
                        ms=5,lw=lwi,ls='None')
      else:
        axs[k].loglog(dat,p,label=lab,lw=lw[j%len(lw)],
                  ls=ls[j%len(ls)],color=color)
      
    if getH:
      axs[k].set_yscale('linear')
    
  
  # For all experiments and data ---------------------------------------
  axs[k].set_xlabel(xlab)
  axs[k].grid(True,linestyle=':')
  axs[k].set_ylim(ylim)
  if xset:
    axs[k].set_xlim(xlim)
  
  if (extxaxis) and (not 'Temperature' in axs[k].get_xlabel()):
    if not xset:
      xlim = np.log10(axs[k].get_xlim())
      if xlim[1]-xlim[0] < 1.1:
        xmean=(xlim[1]+xlim[0])/2
        xlim[0] = 10**(xmean-0.55)
        xlim[1] = 10**(xmean+0.55)
        
        axs[k].set_xlim(xlim)
  
  if not spaghetti:
    axs[k].text(0.95, 0.9, spname,transform=axs[k].transAxes,
                      bbox=bbox_props,ha='right')

if not getH:
  plt.yticks(ytk,ytk)
 
plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=0.4)

# legend   
if not setlegsp:
  legsp=k

if customleg:
  from matplotlib.lines import Line2D
  custom_lines = []
  for i in range(len(clegs)):
    ccol,clw,cls = clegs[i].split(',')
    custom_lines += [Line2D([0], [0], color=ccol, lw=clw,ls=cls)]
  axs[legsp].legend(custom_lines, label,loc=loc)
else:
  if legend:
    axs[legsp].legend(loc=loc)

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
