#!/bin/python

import numpy, pylab, glob
import optparse

systems=['bi', 'car', 'apo']
for sys in systems:
    print sys
    progress=dict()
    file='%s/%s_paths.txt' % (sys, sys)
    fhandle=open(file)
    fluxes=numpy.loadtxt('%s/%s_fluxes.txt' % (sys,sys))
    for (n, path) in enumerate(fhandle.readlines()):
        flux=fluxes[n]/fluxes[0]
        states=path.split()
        ligand='eval-types'
        aucs=[]
        lowci=[]
        hici=[]
        percent=[]
        for (k, state) in enumerate(states):
            percent=round(float(k)/len(states), 2)
            state=int(state)
            if percent not in progress.keys():
                progress[percent]=dict()
                progress[percent]['aucs']=[]
                progress[percent]['lowci']=[]
                progress[percent]['hici']=[]
            mfile='./%s/%s/new-matlab-%s-%s-%s-aucs-95ci.dat' % (sys, ligand, sys, state, ligand)
            value=numpy.loadtxt(mfile, usecols=(0,))
            aucs.append(value)
            progress[percent]['aucs'].append(value)
            data=numpy.loadtxt(mfile, usecols=(1,))
            lowci.append(value-data)
            progress[percent]['lowci'].append(value-data)
            data=numpy.loadtxt(mfile, usecols=(2,))
            hici.append(data-value)
            progress[percent]['hici'].append(data-value)

    pylab.figure()
    low_values=dict()
    hi_values=dict()
    for x in sorted(progress.keys()):
        low_values[x]=numpy.sqrt(numpy.sum([i**2 for i in progress[x]['lowci']]))
        hi_values[x]=numpy.sqrt(numpy.sum([i**2 for i in progress[x]['hici']]))
    lower=[low_values[j] for j in sorted(low_values.keys())]
    upper=[hi_values[j] for j in sorted(hi_values.keys())]
    pylab.errorbar(sorted(progress.keys()), [numpy.mean(progress[i]['aucs']) for i in sorted(progress.keys())], yerr=[lower, upper])
    pylab.show()
    #pylab.plot(sorted(progress.keys()), [low_values[i] for i in sorted(progress.keys())], 'bo') 
    #pylab.plot(sorted(progress.keys()), [hi_values[i] for i in sorted(progress.keys())], 'ro') 
