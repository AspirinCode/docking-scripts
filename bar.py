
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
labels=['Agonist-bound', 'Inv. Agonist-bound', 'Apo']
auclabels=dict()
auclabels['types']='          Agonist vs. Antag.'
auclabels['antagonist']='       Antag. vs. Decoy'
auclabels['agonist']='       Agonist vs. Decoy'


aucnames=['agonist', 'antagonist', 'types'] 
systems=['bi', 'car', 'apo']
colors=['red', 'blue', 'grey' ]
allvalues=dict()
ranges=dict()
index=0
#f, (ax1, ax2)=pylab.subplots(1, 2, sharey=True )
fig=pylab.figure(figsize=(8,6))
gs=gridspec.GridSpec(1, 2, width_ratios=[2,1])
#ax1=fig.add_subplot(gs[0] )
#ax2=fig.add_subplot(gs[1])
ax1=pylab.subplot(gs[0] )
ax2=pylab.subplot(gs[1])
for aucname in aucnames:
    start=index
    allvalues[aucname]=dict()
    for (n, sys) in enumerate(systems):
        if aucname=='agonist':
            ax1.bar(left=range(-3,-2), height=[0,], color=colors[n], label=labels[n])
        type='path'
        dir='./%s/structural-data/' % sys
        pops=numpy.loadtxt('%s/Populations.dat' % dir)
        states=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, aucname), usecols=(0,), dtype=int)
        aucs=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, aucname), usecols=(1,))
        pops=pops[states]
        if sys=='bi':
            cutoff=min(pops)
        else:
            frames=numpy.where(pops>= cutoff)[0]
            pops=pops[frames]
            states=states[frames]
            aucs=aucs[frames]
        allvalues[aucname][sys]=aucs
        index+=1
    avgs=[]
    errs=[]
    for x in ['bi', 'car', 'apo']:
        avgs.append(numpy.mean(allvalues[aucname][x]))
        errs.append(numpy.std(allvalues[aucname][x]))
    width=0.8
    if aucname=='types':
        ax2.bar(left=range(0,3), height=avgs,  width=width, yerr=errs,
            color=['r', 'b', 'grey'], ecolor='k')
    else:
        ax1.bar(left=range(start,index), height=avgs,  width=width, yerr=errs,
            color=['r', 'b', 'grey'], ecolor='k')
    ranges[aucname]=range(start,index)
    pylab.hold(True)

ax2.plot(range(0, 4), [0.5]*4,'k--', label='Random')
ax1.plot(range(0, 7), [0.5]*7, 'k--')
lg=ax1.legend(loc=9)
lg.draw_frame(False)
lg=ax2.legend(loc=9)
lg.draw_frame(False)
axis_label=[]
ax1.set_xlim(0,6)
ax2.set_xlim(0,3)
ax1.set_ylim(0.4, 1.0)
ax2.set_ylim(0.4, 1.0)
ax1.set_ylabel('AUCs')
ax1.xaxis.set_ticks(numpy.arange(0,6))
ax1.xaxis.set_ticklabels([' ', auclabels[aucnames[0]],' ', ' ', auclabels[aucnames[1]],' '])

ax2.yaxis.set_ticks(numpy.arange(0.4, 1.0))
ax2.yaxis.set_ticklabels(['   '*6])
ax2.xaxis.set_ticks(range(0, 3))
ax2.xaxis.set_ticklabels(['   ', auclabels[aucnames[2]], ' '])
pylab.savefig('ligand_state_aucs.png', dpi=300)
pylab.show()

