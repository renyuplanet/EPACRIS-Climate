#!/usr/bin/python

'''
Routine to plot loss and production rates of a specie of one or several experiments 
  The data are taken from ../Experiments/EXP/Output/PAP/rate_coefficients.dat

General notes
  * always put "-" at the beginning of an optional parameter
  * for more than one experiment use " or ' to give them as one argument
  * have a look on a few examples below

Usage: 
  python plotProdLossRates.py "EXP1" "O3"
  * this will display a production and loss rates for ozone

Optional parameter for data:
  -coeff        Plot rate coeffcients (without number density)
  -rmax         Number of maximum reactions shown (default=8)
  -lref         Reference layer: largest reactions from this layer
                Default=-1: largest reactions from any layer from surface to ToA
  -height       y-axis will be in height(km) not pressure(hPa)
  -subdir       This directory will be added before experiment name
                EXAMPLE: python compareT.py "EXP1" -subdir="Earth"
                        -> data taken from ../Experiments/Earth/EXP1/Output

Optional parameter for plot style:
  -figsize      Set figure size in inches -figsize="xsize ysize"
  -title        Set title for entire plot: use e.g. \$ \lamda \$ for math 
  -ls           Set line style for each experiment (will be repeated if nls<nexp)
                Default: "- -- -. :"
  -lw           Set line width for each experiment
  -ylim         Set limits of y-axis: for height in km e.g. -ylim="0 100" 
                                      for pressure in hPa e.g. -ylim="1000 0.1"
  -xliml        Set limits of loss x-axis.
  -xlimp        Set limits of production x-axis.
  -legendloc    Location of legend inside subplot (default=0: best)

Output:
  When nothing is set the plot will be displayed (graphical environment required!)
  otherwise use
  -out="path/filename.extension" (The extension can be e.g. pdf, png or svg)

Examples:
python plotProdLossRates.py "S1.0" "O3" -subdir="Earth" -out="Plots/Earth_LP_O3.pdf"
python plotProdLossRates.py "S0.8 S1.0" "O3" -rmax=5 -coeff

F. Wunderlich, May 2020
'''

import numpy as np
import matplotlib
import sys
from matplotlib import colors as mcolors


def getData(indir,layernum,spec,getH=False):
  '''
  Return pressure in hPa for given layernumber
  '''
  
  # total number density
  ifile2='../Experiments/'+indir+'/Output/Photo/atm_composition.dat'
  nden = np.loadtxt(ifile2,skiprows=1,usecols=(4),unpack=True)
  
  # pressure and number density of specie
  ifile='../Experiments/'+indir+'/Output/Photo/composition_molcc.dat'
  
  f = open(ifile, 'r')
  specnames = np.array(f.readline().split()[:])
  units     = np.array(f.readline().split()[:])
  f.close()
  
  k=np.where(spec==specnames)[0][0]
  
  if getH:
    iz=0
  else:
    iz=2
  
  z,spec = np.loadtxt(ifile,skiprows=2,usecols=(iz,k),unpack=True)
  
  
  return z[layernum-1],nden[layernum-1],spec[layernum-1]
  

def getProdLossRates(spec,indir,layernum,getRates=True,
                     layersum=False,idx_sorting=0,getResidual=False):
  
  '''
  Get loss and production rates or coefficients from model run (PAP files) 
  for given specie
  
  Input
    spec:        specie name (string)
    indir:       path of model experiment (single string or array of strings)
    layernum:    number of model layer (integer scalar or array) 
    getRates:    True: rates will be returned; False: coefficients will be returned
    layersum:    True: sum over all layer; False: all layer will be returned
    idx_sorting: if indir is array: after which idex would be sorted
  
  Output
    scalar, array or matrix (nr,nf,nl) of production rates or coefficients
    scalar, array or matrix (nr,nf,nl) of loss rates or coefficients
    string scalar or array of production reaction (nr,nl)
    string scalar or array of loss reaction (nr,nl)
  '''
  
  def skipuntil(f,fnd):
    '''
    Skip lines until string found in line
    '''
    line=''
    while not fnd in line:
      line=f.readline()

    return line
  
  layernum = np.array(layernum,ndmin=1)
  indir    = np.array(indir,ndmin=1)
  
  nlayer = layernum.size
  nfiles = indir.size
  
  for l in range(nlayer):
  
    for i in range(nfiles):
      f = open('../Experiments/'+indir[i]+'/Output/PAP/PAP_'+str(layernum[l]).zfill(3)+'.dat', 'r')
      
      if not getRates:
        coeff = np.loadtxt('../Experiments/'+indir[i]+'/Output/PAP/rate_coefficients.dat')
      
      lline = skipuntil(f,'reactions')
      nr = int(lline[0:5])
      
      losst = []
      prodt = []
      lossr   = []
      prodr   = []
      
      for j in range(nr):
        line = f.readline()
        rn = line[28:33].strip()
        
        rr = line[34:-1].strip() # whole reaction
        rl = rr[:rr.find('>')-1].split() # prodcution
        rp = rr[rr.find('>')+1:].split() # loss
        
        if spec in rp:
          prodt += [rr]
          
          if getRates:
            prodr   += [float(line[0:25])]
          else:
            prodr   += [coeff[layernum[l],int(rn[1:])]]
        
        if spec in rl:
          losst += [rr]
          if getRates:
            lossr   += [float(line[0:25])]
          else:
            lossr   += [coeff[layernum[l],int(rn[1:])]]
      
      f.close()
      
      if ((i == 0) and (l == 0)):
        lossrates = np.empty((len(lossr),nfiles,nlayer))
        prodrates = np.empty((len(prodr),nfiles,nlayer))
        losstxt   = np.empty((len(lossr),nlayer)).astype("str")
        prodtxt   = np.empty((len(prodr),nlayer)).astype("str")
      
      lossrates[0:len(lossr),i,l] = lossr
      prodrates[0:len(prodr),i,l] = prodr
    
    prodtxt[0:len(prodt),l] = prodt
    losstxt[0:len(losst),l] = losst
    
    
    if getResidual:
      residual = np.sum(prodrates,axis=0)-np.sum(lossrates,axis=0)
      if np.sum(residual) > 0:
        if l==0:
          lossrates = np.append(lossrates,residual[np.newaxis,:,:],axis=0)
          losstxt   = np.append(losstxt,[["residual"]*nlayer],axis=0)
        else:
          lossrates[-1,:,l] = residual[:,l]
      else:
        residual = np.abs(residual)
        if l==0:
          prodrates = np.append(prodrates,residual[np.newaxis,:,:],axis=0)
          prodtxt   = np.append(prodtxt,[["residual"]*nlayer],axis=0)
        else:
          prodrates[-1,:,l] = residual[:,l]
    
  # sorting

  if lref==-1:
    nloss=lossrates[:,0,0].size
    nprod=prodrates[:,0,0].size
    sidx_l = []
    sidx_p = []
    for n in range(nloss):
      for ll in range(nlayer):
        sidx = np.argsort(lossrates[:,0,ll])[::-1]
        if sidx[n] not in sidx_l:
          sidx_l += [sidx[n]]

    for n in range(nprod):
      for ll in range(nlayer):
        sidx = np.argsort(prodrates[:,0,ll])[::-1]
        if sidx[n] not in sidx_p:
          sidx_p += [sidx[n]]
    
    lossrates[:,:,:] = lossrates[sidx_l,:,:]
    losstxt[:,:] = losstxt[sidx_l,:]
    prodrates[:,:,:] = prodrates[sidx_p,:,:]
    prodtxt[:,:] = prodtxt[sidx_p,:]
  else:
    sidx = np.argsort(lossrates[:,0,lref])[::-1]
    
    lossrates[:,:,:] = lossrates[sidx,:,:]
    losstxt[:,:] = losstxt[sidx,:]
    
    sidx = np.argsort(prodrates[:,0,lref])[::-1]
    prodrates[:,:,:] = prodrates[sidx,:,:]
    prodtxt[:,:] = prodtxt[sidx,:]
    
  if layersum:
    prodrates = np.sum(prodrates,axis=2)[:,:,np.newaxis]
    lossrates = np.sum(lossrates,axis=2)[:,:,np.newaxis]
    prodtxt   = prodtxt[:,0:1]
    losstxt   = losstxt[:,0:1]
  
  return prodrates,lossrates,prodtxt,losstxt


# preparation---------------------------------------------------------------------

exps=sys.argv[1]
spec=sys.argv[2]

exps = np.array(exps.split())
nexp=exps.size

if nexp == 1:
  title = exps[0]
else:
  title = ''


colors=['C'+str(x) for x in range(10)]
cols = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Sort colors by hue, saturation, value and name.
by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                for name, color in cols.items())
sorted_names = [name for hsv, name in by_hsv]
colors=colors+sorted_names[17:]


# default values
nlayer=100
layernum=np.arange(1,nlayer)
rmax=8
getRates=True
subdir=''
xunit='s'
label=np.copy(exps)
getH=False
ylab='Pressure [hPa]'
yset=False
xsetl=False
xsetp=False
xlab=r'Molec/cm$^3$/s'
savefig=False
lref=-1
ls=['-','--',':',';']
lw=1
figsize=(0,0)
loc=0

while len(sys.argv)>3:
  # Optional subdir
  if '-subdir' in sys.argv[3]:
    subdir = sys.argv[3][8:]
    subdir+='/'
  elif '-coeff' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][3:]
    getRates=False
    xlab='Rate coefficient'
  elif '-nlayer' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][8:]
    nlayer=int(sys.argv[3])
  elif '-rmax' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][6:]
    rmax=int(sys.argv[3])
  elif '-lref' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][6:]
    lref=int(sys.argv[3])
  elif '-height' == sys.argv[3]:
    getH=True
    ylab='Height [km]'
  elif '-label' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][7:]
    label=np.array(sys.argv[3].split())
    label=[lab.replace('~' , ' ') for lab in label]
    label=[lab.replace('\$' , '$') for lab in label]
    setlab=True
  elif '-ls' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][4:]
    ls=np.array(sys.argv[3].split())
  elif '-lw' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][4:]
    lw=float(sys.argv[3])
  elif '-out' in sys.argv[3]:
    savefig=True
    ofn=sys.argv[3][5:]
    matplotlib.use('Agg')
  elif '-figsize' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][9:]
    figsize=np.array(sys.argv[3].split()).astype('float')
  elif '-title' in sys.argv[3]:
    title = sys.argv[3][7:]
    title=title.replace('\$' , '$')
  elif '-legendloc' in sys.argv[3]:
    loc = int(sys.argv[3][11:])
  elif '-ylim' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][6:]
    ylim=np.array(sys.argv[3].split()).astype('float')
    yset=True
  elif '-xliml' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][7:]
    xliml=np.array(sys.argv[3].split()).astype('float')
    xsetl=True
  elif '-xlimp' in sys.argv[3]:
    sys.argv[3] = sys.argv[3][7:]
    xlimp=np.array(sys.argv[3].split()).astype('float')
    xsetp=True
  else:
    print('Parameter: ',sys.argv[3],' not found')
    exit()
  
  del sys.argv[3]

import matplotlib.pyplot as plt

nl = layernum.size
layernum=np.arange(1,nlayer)

nx = 1
ny = 2

#if nexp > 1:
#  rmax=0
#else:
#  label[:] = 'Sum'

prodrates,lossrates,prodtxt,losstxt = getProdLossRates(spec,exps,layernum,getRates=getRates)
if figsize[0]==0:
  figsize=(5*ny+1,4*nx+1)

fig, axs = plt.subplots(figsize=figsize,nrows=nx,ncols=ny)#, sharey=True)
axs=axs.flatten()

for j in range(nexp):
  
  p,nden,specnden=getData(subdir+exps[j],layernum,spec,getH=getH)
  
  if getRates:
    prodrates[:,j,:] = prodrates[:,j,:]*nden*1.e-9
    lossrates[:,j,:] = lossrates[:,j,:]*nden*1.e-9
  
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
  
  # Production
  nprod = prodrates[:,0,0].size
  if rmax!=-1: nprod=min(rmax*1,nprod)
  
  for i in range(nprod):
    if j==0:
     lab=prodtxt[i,j]
    else:
     lab=''
    
    axs[0].loglog(prodrates[i,j,:],p,label=lab,ls=ls[j],color=colors[i],lw=lw)

  if nexp==1:
    slab='Sum'
  else:
    slab='Sum ('+label[j]+')'
  
  psum = np.sum(prodrates[:,j,:],axis=0)
  #px   = 10**(np.log10(np.max(psum))+0.7)
  #pn   = px * 1e-9
  axs[0].set_title('Sources of '+spec)
  axs[0].loglog(psum,p,'k',ls=ls[j],label=slab,lw=lw)
  #axs[j].set_xlim([pn,px])
  
  # Sinks
  nloss = lossrates[:,0,0].size
  if rmax!=-1: nloss=min(rmax*1,nloss)
  for i in range(nloss):
    if j==0:
     lab=losstxt[i,j]
    else:
     lab=''

    axs[1].loglog(lossrates[i,j,:],p,ls=ls[j],label=lab,color=colors[i],lw=lw)

  ssum = np.sum(lossrates[:,j,:],axis=0)
  #sx   = 10**(np.log10(np.max(ssum))+0.7)
  #sn   = sx * 1e-9
  axs[1].set_title('Sinks of '+spec)
  axs[1].loglog(ssum,p,'k',ls=ls[j],label=slab,lw=lw)
  #axs[j+1].set_xlim([sn,sx])
  
  for k in range(2):
    axs[k].legend(prop={'size': 8},loc=loc) 
    axs[k].grid(True,linestyle=':')
    axs[k].set_ylim(ylim)
    #axs[k].set_yticks(ytk,ytk)
    axs[k].set_ylabel(ylab)
    axs[k].set_xlabel(xlab)
    if getH:
      axs[k].set_yscale('linear')
    
  
  if xsetp:
    axs[0].set_xlim(xlimp)
  if xsetl:
    axs[1].set_xlim(xliml)

  

#title
if len(title)>0:
  plt.suptitle(title)
  #plt.subplots_adjust(top=0.93)

#plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=0.4)

if savefig:
  plt.savefig(ofn)
  plt.close()
else:
  plt.show()
