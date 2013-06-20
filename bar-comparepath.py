
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
                        'size'   : 14}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 12,
          'legend.linewidth': 2}
pylab.rcParams.update(params)

def main(refname, ligandname, keynumber):
    colors=['blue', 'grey', 'purple', 'red']
    labelnames=['Inactive Xtal', 'MSM Ensemble', 'Shaw Random', 'Active Xtal']
    auclabels=dict()
    auclabels['types']='Agonist vs. Antag.'
    auclabels['antagonist']='Antag. vs. Decoy'
    auclabels['agonist']='Agonist vs. Decoy'
    
    aucnames=['agonist', 'antagonist', 'types'] 
    allvalues=dict()
    ranges=dict()
    fig=pylab.figure(figsize=(9,7))
    #gs=gridspec.GridSpec(1, 3,  width_ratios=[2,2, 1])
    #ax1=pylab.subplot(gs[0] )
    #ax2=pylab.subplot(gs[1])
    #ax2=pylab.subplot(gs[2])
    ax1=pylab.subplot(131)
    ax2=pylab.subplot(132)
    ax3=pylab.subplot(133)
    axes=[ax1, ax2, ax3]
    for (n, aucname) in enumerate(aucnames):
        index=0
        print "on %s" % aucname
        allvalues[aucname]=dict()
        start=index
        file='%s_%s_%s_all%svalues.dat' % (aucname, refname, ligandname, keynumber)
        os.system('sed "s/\[//g"  < %s |  sed "s/\]//g" | sed "s/,//g" > mod-%s'
                % (file, file))
        fhandle=open('mod-%s' % file)
        for line in fhandle.readlines():
            key=float(line.split()[0])
            for value in line.split()[1:]:
                value=float(value)
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
        labels=[]
        inactive=numpy.loadtxt('xtal-docking/%s/new-matlab-2rh1-%s-aucs-95ci.dat' % (aucname, aucname), usecols=(0,))
        active=numpy.loadtxt('xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat'
                % (aucname, aucname), usecols=(0,))
        hi=numpy.loadtxt('xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat'
                % (aucname, aucname), usecols=(2,))
        activeerr=hi-active
        hi=numpy.loadtxt('xtal-docking/%s/new-matlab-2rh1-%s-aucs-95ci.dat'
                % (aucname,aucname), usecols=(2,))
        inactiveerr=hi-inactive
        labels.append(0)
        avgs.append(inactive)
        errs.append(inactiveerr)
        index+=1
        for x in sorted(allvalues[aucname].keys()):
            avgs.append(numpy.mean(allvalues[aucname][x]))
            errs.append(numpy.std(allvalues[aucname][x]))
            labels.append(1)
            index+=1
        shawauc=numpy.mean(numpy.loadtxt('bi/bi_shaw_%s_auc_ci.txt' % aucname,
            usecols=(1,)))
        shawerr=numpy.std(shawauc)
        labels.append(2)
        avgs.append(shawauc)
        errs.append(shawerr)
        index+=1
        labels.append(3)
        avgs.append(active)
        errs.append(activeerr)
        index+=1
        width=0.8
        labelcolors=[colors[i] for i in labels]
        axes[n].bar(left=range(start, index), height=avgs,
                    alpha=0.9, width=width, yerr=errs,
                color=labelcolors, ecolor='k') #, label=auclabels[aucname])
        ranges[aucname]=range(start,index)
        pylab.hold(True)

        axes[n].xaxis.set_ticks(numpy.arange(0,index))
        labels=[' ', ' ', '   Inactive', '', '',  '', ' ', '            Active']
        axes[n].xaxis.set_ticklabels(labels)
        axes[n].set_xlabel('States Along Path')
    ax1.bar(left=range(0,1), height=[0,], width=width, yerr=[0,],
            color=colors[0],label=labelnames[0])
    ax1.bar(left=range(0,1), height=[0,], width=width, yerr=[0,],
            color=colors[1],label=labelnames[1])
    ax1.bar(left=range(0,1), height=[0,], width=width, yerr=[0,],
            color=colors[2],label=labelnames[2])
    ax1.bar(left=range(0,1), height=[0,], width=width, yerr=[0,],
            color=colors[3],label=labelnames[3])
    ax1.set_ylabel('AUCs')
    fig.text(0.13, 0.93, auclabels['agonist'], fontsize=16)
    fig.text(0.43, 0.93, auclabels['antagonist'], fontsize=16)
    fig.text(0.68, 0.93, auclabels['types'], fontsize=16)
    ax1.yaxis.set_ticks(numpy.arange(0.6, 1.05, 0.1))
    ax2.yaxis.set_ticks(numpy.arange(0.6, 1.05, 0.1))
    lg=ax1.legend(loc=9)
    lg.draw_frame(False)
    ax1.set_ylim(0.6, 1.0)
    ax2.set_ylim(0.6, 1.0)
    ax1.set_xlim(0, index)
    ax2.set_xlim(0, index)
    ax3.set_xlim(0, index)
    fig.subplots_adjust(hspace = 0.3, wspace = 0.3) 
    # types graph
    length=len(allvalues['types'].keys())
    ax3.set_ylim(0.4, 1.0)
    pylab.savefig('compare_%s_progress_%s_aucs_op%s.png' % (ligandname, refname, keynumber), dpi=300)
    pylab.show()


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-r', '--refname', dest='refname',
                      help='refname allvalues AUC set')
    parser.add_option('-l', '--ligandname', dest='ligandname',
                      help='ligand or all')
    parser.add_option('-k', '--keynumber', dest='keynumber',
                      help='# OPs (1- bulge, 2 - conn, 3 - npxxy, 4 - h36')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(refname=options.refname, ligandname=options.ligandname, keynumber=options.keynumber)


