#!/usr/bin/python

'''
Routine returns if experiment is running in queue and prints last line 
of standard output in command window. 
 * Input can be the experiment name or the queue ID.
 * Giving queue ID may take a while depending on the number of experiments
   in your directory. Use the -subdir option if possible.

Usage:
  python status.py "EXP1 EXP2"
   or
  python status.py "queueID1 queueID2"

Optional:
  -expdir     Experiment directory. Default is: "../Experiments"
  -subdir     This directory will be added before experiment name
              EXAMPLE: python compareExp.py "EXP1" -subdir="Earth"
              -> data taken from ../Experiments/Earth/EXP1/Output
  -log

F. Wunderlich, May 2020
'''

import os
import sys
import numpy as np
import subprocess
import glob
import fnmatch


def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results



# experiments folder
rtdir='../Experiments/'

# Experiment
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
    

# Optional parameters
subdir=''


while len(sys.argv)>2:
  # Optional subdir
  if '-subdir' in sys.argv[2]:
    subdir = sys.argv[2][8:]
    subdir+='/'
  elif '-expdir' in sys.argv[2]:
    rtdir = sys.argv[2][8:]
    rtdir+='/'
  else:
    print('Parameter: ',sys.argv[2],' not found')
    exit()
  
  del sys.argv[2]

# ------------------------------------------------------------

for i in range(nexp):
  try:
    qid=int(exps[i])
    stdoutfile=recursive_glob(rtdir+subdir,'*stdout*'+str(qid)+'*')[0]
    exps=exps.astype('S999')
    exps[i]='/'.join(stdoutfile.split('/')[2:-1])
    
    exps=exps.astype('str')
  except:
    files = os.listdir(rtdir+subdir+'/'+exps[i])
    qid=[]
    for f in files:
      if 'stdout' in f:
      
        fsplit = f.split('_')
        qid+=[int(fsplit[2])]
   
    qid = max(qid)
  
  qstat = subprocess.check_output("qstat")
  qstat = str(qstat)
  qstat = qstat.split('\n')
  
  qlst = []
  for j in range(2,len(qstat)-1):
    qlst += [int(qstat[j][0:7])]
  
  if qid in qlst:
    ii = np.where(qid==np.array(qlst))[0][0]
    print(exps[i] + ' -> ' +qstat[ii+2])
  else:
    print(exps[i] +' not running!')
  
  stoutfile = glob.glob(rtdir+subdir+'/'+exps[i]+'/stdout*'+str(qid)+'*')[0]
  for line in open(stoutfile):
     pass
  print(line)
  
  f=open(rtdir+subdir+'/'+exps[i]+'/Output/err.log')
  errlines=f.readlines()
  for e in errlines:
    if 'Error' in ''.join(errlines):
      print(e[0:-1])
  
  print('-'*73)
  
