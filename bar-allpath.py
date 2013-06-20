
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

def main(refname, keynumber):
    systems=['bi', 'car', 'apo']
    auclabels=dict()
    auclabels['types']='Agonist vs.\n Antag.'
    auclabels['antagonist']='Antag. vs. Decoy'
    auclabels['agonist']='Agonist vs. Decoy'
    
    aucnames=['agonist', 'antagonist', 'types'] 
    colors=['red', 'blue', 'purple' ]
    allvalues=dict()
    ranges=dict()
    index=0
    fig=pylab.figure(figsize=(9,7))
    gs=gridspec.GridSpec(1, 2, width_ratios=[2,1])
    ax1=pylab.subplot(gs[0] )
    ax2=pylab.subplot(gs[1])
    for (n, aucname) in enumerate(aucnames):
        print "on %s" % aucname
        allvalues[aucname]=dict()
        start=index
        file='%s_%s_all%svalues.dat' % (aucname, refname, keynumber)
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
    
    length=len(allvalues['agonist'].keys())+len(allvalues['agonist'].keys())
    ax1.set_ylabel('AUCs')
    ax1.set_xlabel('States Along Path')
    ax1.xaxis.set_ticks(numpy.arange(0,length))
    labels=[' ', ' ', 'Inactive', '', '',  '', ' ', '   Active']
    labels2=[' ', ' ', '    ', '     Inactive', '',  '', ' ', '                Active']
    ax1.xaxis.set_ticklabels(labels+labels2)
    ax1.yaxis.set_ticks(numpy.arange(0.6, 1.05, 0.1))
    lg=ax1.legend(loc=9)
    lg.draw_frame(False)
    ax1.set_ylim(0.6, 1.0)
    # types graph
    length=len(allvalues['types'].keys())
    ax2.set_ylim(0.4, 1.0)
    ax2.plot(range(0, (length+1)), [0.5]*(length+1),'k--', label='Random')
    ax2.set_xlim(0, length)
    ax2.xaxis.set_ticks(numpy.arange(0,length))
    ax2.xaxis.set_ticklabels(labels)
    ax2.set_xlabel('States Along Path')
    lg=ax2.legend(loc=9)
    lg.draw_frame(False)
    #ax1.plot(range(0, 7), [0.5]*7, 'k--')
    #axis_label=[]
    #ax1.set_xlim(0,6)
    #ax2.set_xlim(0,3)
    
    #ax2.yaxis.set_ticks(numpy.arange(0.4, 1.0))
    #ax2.yaxis.set_ticklabels(['   '*6])
    #ax2.xaxis.set_ticks(range(0, 3))
    #ax2.xaxis.set_ticklabels(['   ', auclabels[aucnames[2]], ' '])
    #pylab.savefig('ligand_state_aucs.png', dpi=300)
    pylab.savefig('all_progress_%s_aucs_op%s.png' % (refname, keynumber), dpi=300)
    pylab.show()


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-r', '--refname', dest='refname',
                      help='refname allvalues AUC set')
    parser.add_option('-k', '--keynumber', dest='keynumber',
                      help='# OPs (1- bulge, 2 - conn, 3 - npxxy, 4 - h36')
    #parser.add_option('-x', '--testname', dest='testname',
    #                  help='test AUC set')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(refname=options.refname, keynumber=options.keynumber)


