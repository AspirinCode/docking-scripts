#!/bin/python

import numpy, pylab, glob
import optparse
import operator

def main(sys):
    ligand='eval-types'
    if sys=='apo':
        target=103
    if sys=='bi':
        target=21
    if sys=='car':
        target=44
    print "%s random states" % target
    aucs=dict()
    best_aucs=dict()
    ci_best_aucs=dict()
    print "on %s" % ligand
    high_aucs=dict()
    low_aucs=dict()
    ci_high_aucs=dict()
    ci_low_aucs=dict()

    x_aucs=dict()
    x_hici=dict()
    x_lowci=dict()
    alllowci=dict()
    allhici=dict()
    allaucs=dict()
    xtal=['2rh1', '3p0g']
    for state in xtal:
        mfile='../xtal-docking/%s/new-matlab-%s-%s-aucs-95ci.dat' % (ligand, state, ligand)
        value=numpy.loadtxt(mfile, usecols=(0,))
        x_aucs[state]=value
        data=numpy.loadtxt(mfile, usecols=(1,))
        x_lowci[state]=data
        data=numpy.loadtxt(mfile, usecols=(2,))
        x_hici[state]=data
        print "%s %s AUC %s" % (state, ligand, x_hici[state])
    aucs=dict()
    lowci=dict()
    hici=dict()
    states=range(1,target)
    for state in states:
        state=int(state)
        mfile='./%s/new-matlab-%s-%s-rand-%s-aucs-95ci.dat' % (ligand, sys, state, ligand)
        value=numpy.loadtxt(mfile, usecols=(0,))
        aucs[state]=value
        allaucs[state]=value
        data=numpy.loadtxt(mfile, usecols=(1,))
        lowci[state]=data
        alllowci[state]=data
        data=numpy.loadtxt(mfile, usecols=(2,))
        hici[state]=data
        allhici[state]=data

    # have states in aucs, xtal in x_aucs
    for x in allaucs.keys():
        if allaucs[x] > x_aucs['3p0g']:
            if x in high_aucs.keys():
                if allaucs[x]>high_aucs[x]:
                    high_aucs[x]=allaucs[x]
            else:
                high_aucs[x]=allaucs[x]
        if alllowci[x] > x_hici['3p0g']:
            if x in ci_high_aucs.keys():
                if allaucs[x]>ci_high_aucs[x]:
                    ci_high_aucs[x]=alllowci[x]
            else:
                ci_high_aucs[x]=alllowci[x]

        if allaucs[x] < x_aucs['2rh1']:
            if x in low_aucs.keys():
                if allaucs[x]<low_aucs[x]:
                    low_aucs[x]=allaucs[x]
            else:
                low_aucs[x]=allaucs[x]
        if allhici[x] < x_lowci['2rh1']:
            if x in ci_low_aucs.keys():
                if allaucs[x]<ci_low_aucs[x]:
                    ci_low_aucs[x]=allhici[x]
            else:
                ci_low_aucs[x]=allhici[x]
 
    print "CI best %s AUCs" % len(ci_best_aucs.keys())
    print "not-CI best %s AUCs" % len(best_aucs.keys())
    ohandle=open('%s_rando_types_auc_ci.txt' % (sys), 'w')
    for key in sorted(allaucs.keys()):
        ohandle.write('%s\t%s\t%s\t%s\n' % (key, allaucs[key], alllowci[key],
            allhici[key]))
    print "CI high %s AUCs" % ligand
    print len(ci_high_aucs.keys()), set(ci_high_aucs.keys()), [ci_high_aucs[j]
            for j in ci_high_aucs.keys()]
    print "CI low %s AUCs" % ligand
    print len(ci_low_aucs.keys()), set(ci_low_aucs.keys()), [ci_low_aucs[j]
            for j in ci_low_aucs.keys()]
    print "high %s AUCs" % ligand
    print len(high_aucs.keys()), set(high_aucs.keys()), [high_aucs[j] for j in
            high_aucs.keys()]
    print "low %s AUCs" % ligand
    print len(low_aucs.keys()), set(low_aucs.keys()), [low_aucs[j] for j in low_aucs.keys()]


    
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys)
 
