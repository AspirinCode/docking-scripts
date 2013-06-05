import random
import optparse
import pylab
import numpy
import os
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from numpy import linalg


def bar_progress(data):
    auclabels=dict()
    auclabels['types']='      Discrimination'
    auclabels['antagonist']='Antagonist'
    auclabels['agonist']='Agonist'

    systems=['bi', 'car', 'apo']
    labels=['Agonist-bound', 'Inv. Agnonist bound', 'Apo']
    colors=['red', 'blue', 'grey' ]
    names=['inactive', 'inter', 'active']
    index=0
    fig=pylab.figure(figsize=(8,6))
    #gs=gridspec.GridSpec(1, 2, width_ratios=[2,1])
    ax1=fig.add_subplot(321, aspect=7)
    ax2=fig.add_subplot(322, aspect=7)
    ax3=fig.add_subplot(323, aspect=7)
    ax4=fig.add_subplot(324, aspect=7)
    ax5=fig.add_subplot(325, aspect=7)
    ax6=fig.add_subplot(326, aspect=7)
    ax=[ax1, ax3, ax5, ax2, ax4, ax6]
    width=0.8
    count=0
    for aucname in ['agonist', 'antagonist']:
        limits=(0.6,1.0)
        for (n, sys) in enumerate(systems):
            avgs=[]
            errs=[]
            for key in names:
                avgs.append(numpy.mean(data[aucname][sys][key])) 
                errs.append(numpy.std(data[aucname][sys][key]))
            if aucname=='agonist':
                ax[count].bar(left=range(-3,-2), height=[0,], color=colors[n],
                        label=labels[n])
                lg=ax[count].legend(loc=9)
                lg.draw_frame(False)

            ax[count].bar(left=range(0, len(names)), height=avgs,  
                    width=width, yerr=errs,  color=colors[n], ecolor='k')
            ax[count].set_ylim(limits[0], limits[1])
            ax[count].yaxis.set_ticks(numpy.arange(limits[0],limits[1]+0.05,0.1))
            ax[count].xaxis.set_ticks(numpy.arange(0,4))
            ax[count].xaxis.set_ticklabels([' ', 'Inactive',  ' ', 'Active'])
            count+=1
    pylab.setp( ax1.get_xticklabels(), visible=False)
    pylab.setp( ax2.get_xticklabels(), visible=False)
    pylab.setp( ax3.get_xticklabels(), visible=False)
    pylab.setp( ax4.get_xticklabels(), visible=False)
    pylab.text(0.5, 1.15, 'Agonist vs. Decoys',
                     horizontalalignment='center',
                              fontsize=16, transform = ax1.transAxes)
    pylab.text(0.5, 1.15, 'Antagonist vs. Decoys',
                     horizontalalignment='center',
                              fontsize=16, transform = ax2.transAxes)
    ax3.set_ylabel('AUCs')
    ax5.set_xlabel('Path Progress')
    ax6.set_xlabel('Path Progress')
    fig.subplots_adjust(wspace=0.001)
    pylab.savefig('ligand_path_aucs.png', dpi=300)
    ###########
    aucname='types'
    limits=(0.4,1.0)
    fig=pylab.figure()
    ax1=fig.add_subplot(311, aspect=5)
    ax2=fig.add_subplot(312, aspect=5)
    ax3=fig.add_subplot(313, aspect=5)
    count=0
    ax=[ax1, ax2, ax3]
    for (n, sys) in enumerate(systems):
        avgs=[]
        errs=[]
        for key in names:
            avgs.append(numpy.mean(data[aucname][sys][key])) 
            errs.append(numpy.std(data[aucname][sys][key]))
        if sys=='bi':
            ax[count].plot(range(0,5), [0.5]*5, 'k--', label='Random')
            lg=ax1.legend(loc=9)
            lg.draw_frame(False)
        else:
            ax[count].plot(range(0,5), [0.5]*5, 'k--', label='Random')
        ax[count].bar(left=range(0, len(names)), height=avgs,  
                width=width, yerr=errs, color=colors[n], ecolor='k')
        ax[count].set_ylim(limits[0], limits[1])
        ax[count].yaxis.set_ticks(numpy.arange(limits[0],limits[1]+0.05,0.1))
        ax[count].xaxis.set_ticks(numpy.arange(0,3))
        ax[count].xaxis.set_ticklabels([' ', 'Inactive',  ' ', 'Active'])
        count+=1
    pylab.text(0.5, 1.15, 'Agonist vs. Inv. Ag.',
                     horizontalalignment='center',
                              fontsize=16, transform = ax1.transAxes)
    pylab.setp( ax1.get_xticklabels(), visible=False)
    pylab.setp( ax2.get_xticklabels(), visible=False)
    ax3.set_xlabel('Path Progress')
    ax2.set_ylabel('AUCs')
    #ax1.set_ylabel('AUCs')
    pylab.savefig('types_path_aucs.png', dpi=300)
    pylab.show()

