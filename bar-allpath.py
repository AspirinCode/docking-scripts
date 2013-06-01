
from fitting import *
import random
import optparse
import pylab
import matplotlib.gridspec as gridspec
import numpy
import os
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from numpy import linalg

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                        'size'   : 16}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 14,
          'legend.linewidth': 2}
pylab.rcParams.update(params)

systems=['bi', 'car', 'apo']
auclabels=dict()
auclabels['types']='          Agonist vs. Antag.'
auclabels['antagonist']='       Antag. vs. Decoy'
auclabels['agonist']='       Agonist vs. Decoy'

aucnames=['agonist', 'antagonist', 'types'] 
colors=['red', 'blue', 'purple' ]
allvalues=dict()
ranges=dict()
index=0
fig=pylab.figure(figsize=(8,6))
gs=gridspec.GridSpec(1, 2, width_ratios=[2,1])
ax1=pylab.subplot(gs[0] )
ax2=pylab.subplot(gs[1])
for (n, aucname) in enumerate(aucnames):
    print "on %s" % aucname
    allvalues[aucname]=dict()
    start=index
    fhandle=open('%s_allvalues.dat' % aucname)
    for line in fhandle.readlines():
        key=line.split()[0]
        value=float(line.split()[1])
        try:
            test=int(key)
            if test in allvalues[aucname].keys():
                allvalues[aucname][test].append(value)
            else:
                allvalues[aucname][test]=[]
                allvalues[aucname][test].append(value)
        except ValueError:
            if key in allvalues[aucname].keys():
                allvalues[aucname][key].append(value)
            else:
                allvalues[aucname][key]=[]
                allvalues[aucname][key].append(value)
    print "loaded file data"
    avgs=[]
    errs=[]
    for x in sorted(allvalues[aucname].keys()):
        avgs.append(numpy.mean(allvalues[aucname][x]))
        errs.append(numpy.std(allvalues[aucname][x]))
        index+=1
    width=0.8
    if aucname=='types':
        start=0
        index=len(allvalues[aucname].keys())
        ax2.bar(left=range(start, index), height=avgs,
                alpha=0.9,  width=width, yerr=errs, color=colors[n], ecolor='k', label=auclabels[aucname])
    else:
        ax1.bar(left=range(start, index), height=avgs,
                alpha=0.9, width=width, yerr=errs,
            color=colors[n], ecolor='k', label=auclabels[aucname])
    ranges[aucname]=range(start,index)
    pylab.hold(True)

#ax2.plot(range(0, 4), [0.5]*4,'k--', label='Random')
#ax1.plot(range(0, 7), [0.5]*7, 'k--')
#lg=ax1.legend(loc=9)
#lg.draw_frame(False)
#lg=ax2.legend(loc=9)
#lg.draw_frame(False)
#axis_label=[]
#ax1.set_xlim(0,6)
#ax2.set_xlim(0,3)
ax1.set_ylim(0.4, 1.0)
ax2.set_ylim(0.4, 1.0)
#ax1.set_ylabel('AUCs')
#ax1.xaxis.set_ticks(numpy.arange(0,6))
#ax1.xaxis.set_ticklabels([' ', auclabels[aucnames[0]],' ', ' ', auclabels[aucnames[1]],' '])

#ax2.yaxis.set_ticks(numpy.arange(0.4, 1.0))
#ax2.yaxis.set_ticklabels(['   '*6])
#ax2.xaxis.set_ticks(range(0, 3))
#ax2.xaxis.set_ticklabels(['   ', auclabels[aucnames[2]], ' '])
#pylab.savefig('ligand_state_aucs.png', dpi=300)
pylab.show()

